#!/usr/bin/env python3
r"""
Charge-sector-blocked fast reimplementation of the rung-15 moment pipeline.
===========================================================================

moments_pipeline.instance_moments diagonalizes the FULL 2^Nc-dimensional
Hamiltonian and drags dense 2^Nc x 2^Nc matrices through the derivative
tower -- measured ~37 min/instance at Nc = 12 on this machine.  But every
object in that pipeline is graded by the conserved U(1) charge, so the whole
computation decomposes over popcount sectors (largest block at Nc = 12 is
C(12,6) = 924; eigh(924) is trivial vs eigh(4096)).  This module is a
sector-blocked reimplementation with IDENTICAL semantics: same seed -> same
draw_couplings call -> same disorder instance -> same mu_k numbers (asserted
against the dense path in test_moments_fast.py at Nc = 6, 8, 10).

Derivation of the block structure (each step asserted in the tests):

 1. [H, Q] = 0 (dirac.py), so in the occupation basis H = direct-sum of
    blocks H_q over the popcount-q sectors S_q, dim(S_q) = C(Nc, q).
 2. c_i lowers Q by one: its only nonzero block maps S_q -> S_{q-1} (one
    signed entry per column that has mode i occupied; columns without mode i
    are annihilated).  The derivative chain c^(n) = i[H, c^(n-1)] preserves
    that corner: [H, X] with block-diagonal H and X in the (q-1, q) corner is
    H_{q-1} X_q - X_q H_q, still in the (q-1, q) corner.
 3. M_n = sum_i c_i^dag c_i^(n) is therefore block-diagonal:
    (M_n)_q = sum_i (C_{i,q})^dag (C^(n)_{i,q}) with C_{i,q} the S_q -> S_{q-1}
    block.  Since C_{i,q} is a partial signed permutation (injective on the
    occupied-i columns), the chain step n = 1 and every accumulation
    C^dag_{i,q} X are pure indexing operations -- only chain steps n >= 2 need
    real matrix products.
 4. The unitary particle-hole map C (dirac.particle_hole) is a signed
    COMPLEMENT permutation: C|s> = sign(s) |~s|, so it maps S_q -> S_{Nc-q}.
    Conjugating a block-diagonal X therefore stays block-diagonal, but the
    sector-q block of C X C^dag comes from the MIRROR block q' = Nc - q:
        (C X C^dag)_q = B_{q'} X_{q'} B_{q'}^dag ,   q' = Nc - q,
    with B_{q'} the (signed permutation) S_{q'} -> S_q block of C.  The CP
    projections X_pm = (X +- C X C^dag)/2 thus couple mirror sector PAIRS
    (q, Nc - q) -- there are no off-diagonal blocks anywhere, the mixing is
    entirely between the diagonal blocks of the pair.  Implemented as signed
    index gymnastics (no matmul); asserted against dense cp_conjugate.
 5. Traceless subtraction is a scalar shift of each diagonal block; the
    Hermitian/anti-Hermitian split acts blockwise; V = direct-sum V_q, so
    Wt = V^dag op V is block-diagonal and
        mu_k = sum_q sum_{ab in S_q} W_q[a,b] (E_qa - E_qb)^k .
    The h/a split commutes with the eigenbasis rotation, so we compute
    Wt = V_q^dag op_q V_q ONCE and split Wt afterwards (algebraically
    identical, half the matmuls).

Exact-degeneracy caveat (why "identical" still means ~1e-12, not bitwise):
on a C-invariant instance the mirror sectors are exactly isospectral
(E_q = E_{Nc-q} pairwise), and the dense eigh is free to rotate eigenvectors
inside those cross-sector degenerate pairs.  The moments are invariant under
any rotation within a degenerate eigenspace -- the weight P_E op P_{E'} has a
basis-independent Frobenius norm and omega depends only on the eigenvalues --
so blocked and dense mu_k agree to roundoff (gated at 1e-9 relative in the
tests, observed ~1e-12).  Individual matrix elements of Wt are NOT
comparable between the two paths; only the mu_k are.

What is exact vs measured here: the block decomposition (points 1-5) is an
exact operator statement, asserted numerically; all physics output (mu_k,
w2, kurt and their disorder error bars) is measured ED data identical in
meaning to moments_pipeline.py.

Cost: the dense path is O(n_couplings * 4^Nc) assembly + O(8^Nc) linear
algebra; the blocked path replaces 8^Nc = (sum_q C(Nc,q))^3-ish work by
sum_q C(Nc,q)^3-ish work (ratio ~200x at Nc = 12) and needs only
sum_q C(Nc,q)^2 = C(2Nc, Nc) storage per operator (43 MB at Nc = 12,
642 MB at Nc = 14) instead of 4^Nc.
"""

from __future__ import annotations

import argparse
import json
import time
from math import comb
from pathlib import Path

import numpy as np

from dirac import draw_couplings, apply_monomial, _mode_bit, _below_mask

# Optional path to a concurrently-computed dense-ED partial (per-seed
# instance moments) for the live cross-check; the committed
# moments_nc12.json already records the completed 4-seed check at 5e-15.
DENSE_PARTIAL = None


# --------------------------------------------------------------------------
# sector bookkeeping
# --------------------------------------------------------------------------
def sector_indices(Nc: int):
    """(idx, pos, sec): idx[q] = ascending full-basis indices of the popcount-q
    sector; pos[s] = local position of state s inside its sector; sec[s] = its
    popcount.  Ascending order matters: it makes the blocked H_q literally the
    slice H[np.ix_(idx[q], idx[q])] of the dense builder's H."""
    dim = 1 << Nc
    states = np.arange(dim, dtype=np.int64)
    sec = np.bitwise_count(states).astype(np.int64)
    idx = [np.nonzero(sec == q)[0] for q in range(Nc + 1)]
    pos = np.empty(dim, dtype=np.int64)
    for q in range(Nc + 1):
        pos[idx[q]] = np.arange(len(idx[q]))
    return idx, pos, sec


def blocked_hamiltonian(Nc: int, p: int, couplings: dict, idx, pos, sec):
    """[H_0, ..., H_Nc]: the charge-sector blocks of dirac_hamiltonian for the
    same coupling table, assembled monomial-by-monomial with the SAME
    apply_monomial as the dense builder (asserted equal to dense slices)."""
    dim = 1 << Nc
    states = np.arange(dim, dtype=np.int64)
    blocks = [np.zeros((len(idx[q]), len(idx[q])), dtype=complex)
              for q in range(Nc + 1)]
    rows_l, cols_l, vals_l = [], [], []
    for (A, B), J in couplings.items():
        valid, rows, signs = apply_monomial(Nc, A, B, states)
        r, c = rows[valid], states[valid]
        rows_l.append(r)
        cols_l.append(c)
        vals_l.append(J * signs[valid])
    r = np.concatenate(rows_l) if rows_l else np.empty(0, dtype=np.int64)
    c = np.concatenate(cols_l) if cols_l else np.empty(0, dtype=np.int64)
    v = np.concatenate(vals_l) if vals_l else np.empty(0, dtype=complex)
    # charge conservation of every monomial (|A| = |B|): row sector == col sector
    assert np.array_equal(sec[r], sec[c]), "monomial broke charge conservation"
    for q in range(Nc + 1):
        m = sec[c] == q
        np.add.at(blocks[q], (pos[r[m]], pos[c[m]]), v[m])
    return blocks


def blocked_annihilators(Nc: int, idx, pos):
    """cmaps[i][q] = (cl, rl, sg): the S_q -> S_{q-1} block of c_i as a partial
    signed permutation -- column cl[j] (local, sector q) maps to row rl[j]
    (local, sector q-1) with sign sg[j].  Built from the same apply_monomial
    as the dense annihilators."""
    cmaps = []
    for i in range(Nc):
        per_q = {}
        for q in range(1, Nc + 1):
            valid, rows, signs = apply_monomial(Nc, (), (i,), idx[q])
            per_q[q] = (pos[idx[q][valid]], pos[rows[valid]], signs[valid])
        cmaps.append(per_q)
    return cmaps


def blocked_particle_hole(Nc: int, idx, pos):
    """ph[q] = (r, sg): the S_q -> S_{Nc-q} block of the unitary C as a signed
    permutation -- local column c of sector q maps to local row r[c] of sector
    Nc-q with sign sg[c].  Same right-to-left factor loop as
    dirac.particle_hole (C|s> = sign(s)|~s>, asserted in the tests)."""
    dim = 1 << Nc
    states = np.arange(dim, dtype=np.int64)
    perm = states.copy()
    sign = np.ones(dim)
    for i in reversed(range(Nc)):
        bit = _mode_bit(Nc, i)
        below = _below_mask(Nc, i)
        sign = sign * (1.0 - 2.0 * (np.bitwise_count(perm & below) & 1))
        perm = perm ^ bit
    assert np.array_equal(perm, states ^ (dim - 1)), "C is not the complement map"
    return {q: (pos[perm[idx[q]]], sign[idx[q]]) for q in range(Nc + 1)}


def cp_conjugate_blocked(Nc: int, ph, X_blocks):
    """Blocks of C X C^dag for block-diagonal X: the sector-q output block is
    the signed-permuted MIRROR block q' = Nc - q of X (derivation pt. 4)."""
    out = []
    for q in range(Nc + 1):
        qp = Nc - q
        r, sg = ph[qp]
        M = (sg[:, None] * sg[None, :]) * X_blocks[qp]
        Y = np.empty_like(M)
        Y[np.ix_(r, r)] = M
        out.append(Y)
    return out


# --------------------------------------------------------------------------
# blocked derivative tower
# --------------------------------------------------------------------------
def blocked_mn_tower(Nc: int, H_blocks, cmaps, n_max: int, keep=None):
    """{n: [blocks q = 0..Nc]} of M_n = sum_i c_i^dag c_i^(n) for n in `keep`
    (default: all 0..n_max).  Chain step n = 1 and every accumulation
    C^dag X are indexing operations (C_{i,q} is a partial signed permutation);
    only chain steps n >= 2 pay for matmuls."""
    if keep is None:
        keep = set(range(n_max + 1))
    keep = set(keep)
    dims = [b.shape[0] for b in H_blocks]
    towers = {n: [np.zeros((d, d), dtype=complex) for d in dims] for n in keep}
    for q in range(1, Nc + 1):
        Hq, Hqm = H_blocks[q], H_blocks[q - 1]
        for i in range(Nc):
            cl, rl, sg = cmaps[i][q]
            if 0 in keep:       # C^dag C = occupation projector of mode i
                towers[0][q][cl, cl] += 1.0
            if n_max == 0:
                continue
            # n = 1: C1 = i(H_{q-1} C0 - C0 H_q), via column/row selection
            C = np.zeros((dims[q - 1], dims[q]), dtype=complex)
            C[:, cl] = Hqm[:, rl] * sg[None, :]
            C[rl, :] -= sg[:, None] * Hq[cl, :]
            C *= 1j
            if 1 in keep:
                towers[1][q][cl, :] += sg[:, None] * C[rl, :]
            for n in range(2, n_max + 1):
                C = 1j * (Hqm @ C - C @ Hq)
                if n in keep:
                    towers[n][q][cl, :] += sg[:, None] * C[rl, :]
    return towers


# --------------------------------------------------------------------------
# the pipeline, blocked
# --------------------------------------------------------------------------
def instance_moments_blocked(Nc: int, p: int, seed: int, n_list=(2, 3),
                             n_max=None):
    """Drop-in replacement for moments_pipeline.instance_moments: same seed ->
    same couplings -> same {(n, 'right'/'wrong'): {mu0, mu1, mu2, mu4,
    escale}} numbers (to roundoff; see the degeneracy caveat above)."""
    if n_max is None:
        n_max = max(n_list)
    rng = np.random.default_rng(seed)
    couplings = draw_couplings(Nc, p, rng, "c_symmetric")
    dim = 1 << Nc
    idx, pos, sec = sector_indices(Nc)
    Hb = blocked_hamiltonian(Nc, p, couplings, idx, pos, sec)
    Es, Vs = [], []
    sum_E2 = 0.0
    for q in range(Nc + 1):
        E, V = np.linalg.eigh(Hb[q])
        Es.append(E)
        Vs.append(V)
        sum_E2 += float(np.sum(E ** 2))
    escale = float(np.sqrt(sum_E2 / dim))
    cmaps = blocked_annihilators(Nc, idx, pos)
    ph = blocked_particle_hole(Nc, idx, pos)
    towers = blocked_mn_tower(Nc, Hb, cmaps, n_max, keep=set(n_list))
    del Hb

    out = {}
    for n in n_list:
        tr = complex(sum(np.trace(b) for b in towers[n]))
        Mb = [towers[n][q] - (tr / dim) * np.eye(towers[n][q].shape[0])
              for q in range(Nc + 1)]
        Mc = cp_conjugate_blocked(Nc, ph, Mb)
        Mp = [0.5 * (Mb[q] + Mc[q]) for q in range(Nc + 1)]
        Mm = [0.5 * (Mb[q] - Mc[q]) for q in range(Nc + 1)]
        cp_target = (-1) ** (n + 1)
        chans = {"right": Mp if cp_target == +1 else Mm,
                 "wrong": Mm if cp_target == +1 else Mp}
        for tag, opb in chans.items():
            mu = {k: 0.0 for k in (0, 1, 2, 4)}
            for q in range(Nc + 1):
                Wt = Vs[q].conj().T @ opb[q] @ Vs[q]
                Wt_h = 0.5 * (Wt + Wt.conj().T)
                Wt_a = 0.5 * (Wt - Wt.conj().T)
                W = (np.abs(Wt_h) ** 2 + np.abs(Wt_a) ** 2) / dim
                om = Es[q][:, None] - Es[q][None, :]
                for k in (0, 1, 2, 4):
                    mu[k] += float(np.sum(W * om ** k))
            out[(n, tag)] = dict(mu0=mu[0], mu1=mu[1], mu2=mu[2], mu4=mu[4],
                                 escale=escale)
    return out


def rows_from_instances(Nc, p, inst_list, n_list=(2, 3)):
    """ensemble_moments-format rows (same keys, same ddof=1 SEM math) from a
    list of instance_moments dicts."""
    rows = {key: [] for n in n_list for key in [(n, "right"), (n, "wrong")]}
    for inst in inst_list:
        for key, d in inst.items():
            w2 = d["mu2"] / (d["mu0"] * d["escale"] ** 2)
            kurt = d["mu4"] * d["mu0"] / d["mu2"] ** 2
            odd = abs(d["mu1"]) / np.sqrt(d["mu0"] * d["mu2"])
            rows[key].append((d["mu0"], w2, kurt, odd))
    result = []
    for (n, tag), vals in rows.items():
        arr = np.array(vals)
        mean = arr.mean(axis=0)
        sem = arr.std(axis=0, ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 \
            else 0 * mean
        result.append(dict(
            Nc=Nc, p=p, n=n, channel=tag, n_inst=len(inst_list),
            mu0=mean[0], mu0_err=float(sem[0]),
            w2=mean[1], w2_err=float(sem[1]),
            kurt=mean[2], kurt_err=float(sem[2]),
            odd_moment_check=float(mean[3])))
    return result


def ensemble_moments_blocked(Nc: int, p: int, n_inst: int, seed0: int = 100,
                             n_list=(2, 3)):
    """Blocked twin of moments_pipeline.ensemble_moments (identical row
    format and statistics)."""
    insts = [instance_moments_blocked(Nc, p, seed0 + j, n_list)
             for j in range(n_inst)]
    return rows_from_instances(Nc, p, insts, n_list)


# --------------------------------------------------------------------------
# cost model (for the Nc = 14 go/no-go and honest timing extrapolation)
# --------------------------------------------------------------------------
def _flop_units(Nc: int, n_list=(2, 3)):
    """Relative cost in complex-MAC units of one blocked instance: derivative
    chains (2 matmuls per chain step n >= 2, per mode, per sector) + the
    weight stage (2 matmuls of d^3 per sector per (n, channel)) + eigh
    (~9 d^3 per sector, crude).  Used only for RATIOS between Nc values."""
    dims = [comb(Nc, q) for q in range(Nc + 1)]
    n_max = max(n_list)
    chain = sum(dims[q - 1] ** 2 * dims[q] + dims[q - 1] * dims[q] ** 2
                for q in range(1, Nc + 1))
    tower = Nc * (n_max - 1) * chain
    d3 = sum(d ** 3 for d in dims)
    weights = 2 * len(n_list) * 2 * d3
    eigh = 9 * d3
    return tower + weights + eigh


def peak_memory_gb(Nc: int, n_list=(2, 3)):
    """Peak complex128 working set of one blocked instance, counted:
    H blocks + eigenvector blocks + len(n_list) tower block-sets (each
    sum_q C(Nc,q)^2 = C(2Nc, Nc) entries) + 3 transient chain/CP matrices of
    the largest block size."""
    per_op = comb(2 * Nc, Nc)
    dmax = comb(Nc, Nc // 2)
    entries = (3 + len(n_list)) * per_op + 3 * dmax ** 2
    return 16 * entries / 1e9


# --------------------------------------------------------------------------
# production main: Nc = 12 ensemble, dense cross-check, fixed-p trend
# --------------------------------------------------------------------------
def _load_dense_partial(path=DENSE_PARTIAL):
    if path is None:
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def crosscheck_against_dense(inst_by_seed, dense_rows):
    """Compare blocked instance moments to concurrently computed dense ones
    (moments_nc12_partial.json rows).  Returns (n_checked, max_rel)."""
    max_rel = 0.0
    n_checked = 0
    for rec in dense_rows or []:
        seed = rec["seed"]
        if seed not in inst_by_seed:
            continue
        mine = inst_by_seed[seed]
        for (n, tag), d in mine.items():
            ref = rec[f"{n}:{tag}"]
            for key in ("mu0", "mu2", "mu4", "escale"):
                rel = abs(d[key] - ref[key]) / max(abs(ref[key]), 1e-300)
                max_rel = max(max_rel, rel)
            # mu1 is a zero: compare on the pipeline's normalized odd scale
            scale = np.sqrt(d["mu0"] * d["mu2"])
            max_rel = max(max_rel, abs(d["mu1"] - ref["mu1"]) / scale)
        n_checked += 1
    return n_checked, max_rel


def trend_table(new_rows, n_list=(2, 3)):
    """Fixed-p trend of the shape invariants vs 1/Nc: committed moments.json
    rows (Nc = 6, 8, 10) + the new rows.  Adds the per-Nc CP-channel w2
    splitting (wrong - right) with quadrature errors -- a CAVEAT applies: the
    two channels share disorder instances, so quadrature ignores their
    correlation; for the new rows the paired per-instance SEM is also given
    and is the number to trust there."""
    committed = Path(__file__).resolve().parent / "moments.json"
    with open(committed) as fh:
        old_rows = json.load(fh)
    rows = old_rows + list(new_rows)
    by = {}
    for r in rows:
        by[(r["Nc"], r["n"], r["channel"])] = r
    ncs = sorted({r["Nc"] for r in rows})
    trend = []
    for n in n_list:
        for Nc in ncs:
            right = by[(Nc, n, "right")]
            wrong = by[(Nc, n, "wrong")]
            split = wrong["w2"] - right["w2"]
            split_err = float(np.hypot(wrong["w2_err"], right["w2_err"]))
            trend.append(dict(
                n=n, Nc=Nc, lam=16.0 / Nc, inv_Nc=1.0 / Nc,
                n_inst=right["n_inst"],
                w2_right=right["w2"], w2_right_err=right["w2_err"],
                w2_wrong=wrong["w2"], w2_wrong_err=wrong["w2_err"],
                kurt_right=right["kurt"], kurt_right_err=right["kurt_err"],
                kurt_wrong=wrong["kurt"], kurt_wrong_err=wrong["kurt_err"],
                w2_split=split, w2_split_err=split_err))
    return trend


def paired_split_sem(inst_list, n):
    """Mean +- SEM (ddof=1) of the PER-INSTANCE w2 difference (wrong - right)
    -- the correlation-honest error bar for the channel splitting."""
    diffs = []
    for inst in inst_list:
        w2 = {tag: inst[(n, tag)]["mu2"] /
              (inst[(n, tag)]["mu0"] * inst[(n, tag)]["escale"] ** 2)
              for tag in ("right", "wrong")}
        diffs.append(w2["wrong"] - w2["right"])
    arr = np.array(diffs)
    sem = arr.std(ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0.0
    return float(arr.mean()), float(sem)



def trend_plot(payload, out_path):
    """moments_trend.png: w2 and kurt per CP channel vs 1/Nc at fixed p = 4
    (left edge = the lambda -> 0 matching corner).  Reads the payload written
    by main(); an Nc = 14 single-realization point, if present, is drawn as an
    open marker (no error bar -- one instance)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    trend = payload["trend"]
    nc14 = payload.get("nc14_single_instance")
    colors = {"right": "#2166ac", "wrong": "#b2182b"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    for col, n in enumerate(sorted({t["n"] for t in trend})):
        rows = sorted([t for t in trend if t["n"] == n], key=lambda t: t["Nc"])
        x = [t["inv_Nc"] for t in rows]
        for row, obs in enumerate(("w2", "kurt")):
            ax = axes[row][col]
            for tag in ("right", "wrong"):
                ax.errorbar(x, [t[f"{obs}_{tag}"] for t in rows],
                            yerr=[t[f"{obs}_{tag}_err"] for t in rows],
                            fmt="o-", ms=4, lw=1.0, capsize=3,
                            color=colors[tag], label=f"{tag} CP channel")
                if nc14 is not None and f"{n}:{tag}" in nc14:
                    ax.plot(1.0 / 14.0, nc14[f"{n}:{tag}"][obs], "o",
                            ms=6, mfc="none", color=colors[tag])
            ax.grid(True, alpha=0.2)
            if row == 0:
                ax.set_title(f"$n = {n}$   (CP$_{{target}} = {(-1)**(n+1):+d}$)")
            else:
                ax.set_xlabel(r"$1/N_c$   ($\lambda = 16/N_c$; "
                              r"left edge = matching corner)")
        axes[0][col].set_ylabel(r"$w_2 = \mu_2/(\mu_0\,\langle H^2\rangle)$")
        axes[1][col].set_ylabel(r"kurt $= \mu_4\mu_0/\mu_2^2$")
    axes[0][0].legend(frameon=False, fontsize=9)
    if nc14 is not None:
        axes[0][1].text(0.03, 0.95, "open markers: single $N_c=14$ instance\n"
                        "(no error bar)", transform=axes[0][1].transAxes,
                        fontsize=8, color="0.4", va="top")
    fig.suptitle("CP-resolved $M_n$ shape invariants, fixed $p=4$ "
                 "(charge-sector-blocked ED)", y=0.995)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--nc14", action="store_true",
                    help="attempt ONE Nc=14 instance (seed 100) if the "
                         "measured Nc=12 timing extrapolates below ~40 min")
    ap.add_argument("--n-inst", type=int, default=None,
                    help="override the timing-based Nc=12 instance count")
    ap.add_argument("--dense-partial", default=None,
                    help="path to a concurrent dense-ED partial JSON for the live cross-check")
    ap.add_argument("--plot", action="store_true",
                    help="regenerate moments_trend.png from the existing "
                         "moments_nc12.json and exit (no recompute)")
    args = ap.parse_args()

    if args.plot:
        here = Path(__file__).resolve().parent
        with open(here / "moments_nc12.json") as fh:
            payload = json.load(fh)
        trend_plot(payload, here / "moments_trend.png")
        return

    Nc, p, n_list, seed0 = 12, 4, (2, 3), 100
    here = Path(__file__).resolve().parent
    out_path = here / "moments_nc12.json"

    print(f"blocked moment pipeline: Nc = {Nc}, p = {p}, "
          f"lambda = {p * p / Nc:.3f}, largest block C(12,6) = {comb(12, 6)}")
    print(f"predicted peak memory: {peak_memory_gb(Nc):.2f} GB\n")

    # -- production ensemble, timing-adaptive instance count ---------------
    t0 = time.perf_counter()
    insts = [instance_moments_blocked(Nc, p, seed0, n_list)]
    t_first = time.perf_counter() - t0
    if args.n_inst is not None:
        n_inst = args.n_inst
    elif t_first > 90.0:
        n_inst = 12
    elif t_first > 45.0:
        n_inst = 16
    elif t_first > 15.0:
        n_inst = 20
    else:
        n_inst = 40          # match the committed Nc = 6, 8 statistics
    print(f"instance seed={seed0}: {t_first:.1f} s  "
          f"(dense path measured ~37 min) -> n_inst = {n_inst}")
    times = [t_first]
    for j in range(1, n_inst):
        t0 = time.perf_counter()
        insts.append(instance_moments_blocked(Nc, p, seed0 + j, n_list))
        times.append(time.perf_counter() - t0)
        print(f"instance seed={seed0 + j}: {times[-1]:.1f} s", flush=True)
    t_mean = float(np.mean(times))

    rows = rows_from_instances(Nc, p, insts, n_list)
    print(f"\n{'Nc':>3} {'n':>2} {'chan':>6} {'mu0':>12} "
          f"{'w2':>18} {'kurt':>18}")
    for r in rows:
        print(f"{r['Nc']:>3} {r['n']:>2} {r['channel']:>6} "
              f"{r['mu0']:>12.4e} {r['w2']:>10.4f} +- {r['w2_err']:<6.4f} "
              f"{r['kurt']:>10.3f} +- {r['kurt_err']:<6.3f}")

    # -- cross-check against the concurrent dense job ----------------------
    inst_by_seed = {seed0 + j: insts[j] for j in range(n_inst)}
    dense_rows = _load_dense_partial(args.dense_partial)
    if dense_rows:
        n_chk, max_rel = crosscheck_against_dense(inst_by_seed, dense_rows)
        print(f"\ndense cross-check: {n_chk} seed(s) "
              f"{sorted(r['seed'] for r in dense_rows)}, "
              f"max relative deviation = {max_rel:.2e}")
        assert max_rel < 1e-8, "blocked path disagrees with dense ED"
        crosscheck = dict(n_seeds=n_chk,
                          seeds=sorted(r["seed"] for r in dense_rows),
                          max_rel=max_rel)
    else:
        print("\ndense cross-check: partial file absent/empty -- skipped "
              "(gated instead by the Nc = 6/8/10 dense-agreement tests)")
        crosscheck = dict(n_seeds=0, seeds=[], max_rel=None)

    # -- fixed-p trend toward the lambda -> 0 corner -----------------------
    trend = trend_table(rows, n_list)
    print(f"\nfixed-p trend, p = {p} (lambda = 16/Nc):")
    print(f"{'n':>2} {'Nc':>3} {'lam':>5} {'w2 right':>18} {'w2 wrong':>18} "
          f"{'split':>16} {'kurt right':>16} {'kurt wrong':>16}")
    for t in trend:
        print(f"{t['n']:>2} {t['Nc']:>3} {t['lam']:>5.2f} "
              f"{t['w2_right']:>10.4f} +- {t['w2_right_err']:<6.4f}"
              f"{t['w2_wrong']:>10.4f} +- {t['w2_wrong_err']:<6.4f}"
              f"{t['w2_split']:>9.4f} +- {t['w2_split_err']:<5.4f}"
              f"{t['kurt_right']:>9.3f} +- {t['kurt_right_err']:<5.3f}"
              f"{t['kurt_wrong']:>9.3f} +- {t['kurt_wrong_err']:<5.3f}")
    paired = {}
    for n in n_list:
        m, s = paired_split_sem(insts, n)
        paired[str(n)] = dict(mean=m, sem=s)
        print(f"  Nc=12 paired per-instance w2 split, n={n}: "
              f"{m:.4f} +- {s:.4f}  (correlation-honest SEM)")

    # -- optional single Nc = 14 instance ----------------------------------
    nc14_result = None
    ratio = _flop_units(14, n_list) / _flop_units(12, n_list)
    t14_pred = t_mean * ratio
    print(f"\nNc=14 extrapolation: {t_mean:.0f} s/instance x flop ratio "
          f"{ratio:.1f} = {t14_pred / 60:.0f} min predicted, "
          f"peak memory {peak_memory_gb(14):.2f} GB")
    if args.nc14:
        if t14_pred < 40 * 60:
            t0 = time.perf_counter()
            inst14 = instance_moments_blocked(14, p, seed0, n_list)
            t14 = time.perf_counter() - t0
            nc14_result = {"seed": seed0, "t_seconds": t14}
            print(f"Nc=14 seed={seed0} done in {t14 / 60:.1f} min "
                  f"(SINGLE realization -- no error bars):")
            for (n, tag), d in inst14.items():
                w2 = d["mu2"] / (d["mu0"] * d["escale"] ** 2)
                kurt = d["mu4"] * d["mu0"] / d["mu2"] ** 2
                nc14_result[f"{n}:{tag}"] = dict(
                    d, w2=w2, kurt=kurt)
                print(f"  n={n} {tag:>6}: mu0={d['mu0']:.4e} "
                      f"w2={w2:.4f} kurt={kurt:.3f}")
        else:
            print("Nc=14 skipped: predicted time exceeds the 40 min budget")
    else:
        print("Nc=14 not attempted (pass --nc14)")

    payload = dict(
        meta=dict(Nc=Nc, p=p, n_list=list(n_list), ensemble="c_symmetric",
                  seeds=list(range(seed0, seed0 + n_inst)),
                  mean_seconds_per_instance=t_mean,
                  method="charge-sector-blocked ED (moments_fast.py), "
                         "semantics identical to moments_pipeline"),
        rows=rows,
        paired_w2_split_nc12=paired,
        dense_crosscheck=crosscheck,
        trend=trend,
        nc14_single_instance=nc14_result)
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
