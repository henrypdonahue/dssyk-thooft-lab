#!/usr/bin/env python3
r"""
Rung 27: freeness onset -- the early/late free-product transition, measured.
============================================================================

Cui-Kolchmeyer 2607.13665 (Sec. 4) argue that past the scrambling time the
algebras of early and late observer operators combine as a FREE PRODUCT: all
mixed moments are governed by free independence with respect to the trace
state.  At T_B = infinity the SYK state IS the normalized trace tau = Tr/dim,
so the claim becomes a concrete, instance-level statement about mixed moments
of a(0) and b(t) -- directly measurable by ED.  Nobody has exhibited (or
refuted) this transition in a microscopic dS-candidate model; this module is
the first data (road_map rung 27).

The design dodges the graded-freeness subtlety of odd (fermionic) variables:
the HEADLINE pair is bosonic,

    E = i psi_0 psi_1,      F(t) = i psi_2(t) psi_3(t),

both Hermitian, traceless, and squaring to 1 (spectrum +-1, half-half).

  * t = 0: E and F live on disjoint mode pairs, so they COMMUTE and the
    trace factorizes -- classical independence.  The spectrum of E + F is
    then EXACTLY the binomial law {-2: 1/4, 0: 1/2, +2: 1/4}.
  * Free independence: for two symmetries with spectrum +-1 the pair
    distribution is completely determined by the traces of alternating
    words, and E, F free  <=>  spec(E + F) is the ARCSINE law on (-2, 2),
    density 1/(pi sqrt(4 - x^2))  (the free convolution
    Bernoulli [+] Bernoulli; the classical convolution is the binomial).

So the Wasserstein-1 distances of spec(E + F(t)) to the binomial and to the
arcsine law are a complete, parameter-free order parameter for the
classical -> free crossover of the pair (E, F(t)).  (W1 = integrated
|CDF difference|: exact and stable for atomic spectra, where the KS sup is
ill-conditioned at shared atoms.)  Alongside, the alternating-word battery
(free prediction 0 for every entry):

    m_k(t)  = tau[(E F(t))^k],  k = 1, 2, 3
              (m_1 = the 2-pt function -- decays by ordinary thermalization;
               m_2 = the even OTOC, classical value 1 at t = 0;
               m_3 starts at 0 and must STAY 0 under freeness),
    w_4(t)  = tau[a b(t) a b(t)],  a = psi_0, b = psi_2
              (the standard Majorana OTOC, = -1 at t = 0; the odd-variable
               counterpart, reported as auxiliary data).

Everything is evaluated in the eigenbasis of one disorder realization
(F(t) is a phase-Hadamard product -- O(dim^2) per time), disorder-averaged
with SEM error bars; per-instance KS distances (freeness is claimed per
instance; singlet self-averaging is the measured 1/C(N,p)).

What a verdict looks like
-------------------------
  * Onset: W1_bin rises from 0 while W1_arc falls from its t = 0 value
    W1(binomial, arcsine) = 1 - 4(sqrt(2)-1)/pi ~= 0.4727 (closed form,
    asserted against quadrature in the tests); the crossover time t_x
    (W1_arc = W1_bin) against 1/lambda_L ~ 1/(2 scriptJ) [rung 25:
    lambda_L -> 2 scriptJ at T_B = infinity].
  * Depth: the late-time plateau of W1_arc and of |m_k|.  Exact freeness at
    finite dim is impossible (random-matrix fluctuations of the spectrum
    give W1 ~ dim^{-1/2}-ish floors), so the honest statement is the
    N-trend of the plateau: falling with dim supports "free product in the
    large-N limit"; an N-independent floor refutes it.

Time is reported in units of 1/scriptJ, scriptJ = sqrt(2) p (the repo's
ED-locked conventions bridge, duality/LARGEQ.md).

Run:  python3 freeness.py [--quick]     (production writes freeness.json)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from pauli_strings import majorana_hamiltonian_fast
from syk import majorana_operators

TJ_GRID = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0,
           2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 14.0, 20.0)


# ---------------------------------------------------------------------------
# reference laws for spec(E + F)
# ---------------------------------------------------------------------------
def arcsine_cdf(x: np.ndarray) -> np.ndarray:
    """CDF of the arcsine law on (-2, 2): the free convolution of two
    symmetric +-1 Bernoullis."""
    return 0.5 + np.arcsin(np.clip(x, -2.0, 2.0) / 2.0) / np.pi


def binomial_cdf(x: np.ndarray) -> np.ndarray:
    """CDF of the classical convolution: atoms {-2: 1/4, 0: 1/2, +2: 1/4}."""
    return np.where(x < -2.0, 0.0,
                    np.where(x < 0.0, 0.25,
                             np.where(x < 2.0, 0.75, 1.0)))


_W1_GRID = np.linspace(-2.0005, 2.0005, 8001)   # midpoint-free of the atoms
_W1_DX = _W1_GRID[1] - _W1_GRID[0]


def w1_distance(eigs: np.ndarray, cdf) -> float:
    """Wasserstein-1 distance of a spectrum to a reference law on [-2, 2]:
    integral of |empirical CDF - cdf| on a fine grid.  Exact up to O(grid
    spacing) and, unlike the KS sup, well-conditioned when both measures
    carry atoms at the same points."""
    ecdf = np.searchsorted(np.sort(eigs), _W1_GRID, side="right") / len(eigs)
    return float(np.sum(np.abs(ecdf - cdf(_W1_GRID))) * _W1_DX)


# W1(binomial, arcsine) = 1 - 4(sqrt(2)-1)/pi: the t = 0 value of W1_arc
# and the natural scale of both distances (asserted vs quadrature in tests).
W1_BINOMIAL_VS_ARCSINE = 1.0 - 4.0 * (np.sqrt(2.0) - 1.0) / np.pi


# ---------------------------------------------------------------------------
# one disorder instance
# ---------------------------------------------------------------------------
def _dense(op) -> np.ndarray:
    return op.toarray() if hasattr(op, "toarray") else np.asarray(op)


def instance_operators(N: int, p: int, rng):
    """Eigen-decompose one realization; return (lam, ops-in-eigenbasis) for
    E = i psi0 psi1, F = i psi2 psi3, a = psi0, b = psi2."""
    H = majorana_hamiltonian_fast(N, p, rng)
    lam, V = np.linalg.eigh(H)
    psis = majorana_operators(N)
    E = 1j * (psis[0] @ psis[1]).toarray()
    F = 1j * (psis[2] @ psis[3]).toarray()
    a = psis[0].toarray()
    b = psis[2].toarray()
    Vh = V.conj().T
    return lam, (Vh @ E @ V, Vh @ F @ V, Vh @ a @ V, Vh @ b @ V)


def heisenberg(op_tilde: np.ndarray, lam: np.ndarray, t: float) -> np.ndarray:
    """op(t) in the eigenbasis: elementwise phase product, O(dim^2)."""
    phase = np.exp(1j * t * (lam[:, None] - lam[None, :]))
    return phase * op_tilde


def measure_instance(lam, ops, ts, snap_idx=()) -> dict:
    """All observables for one realization on a time grid (real time t).
    snap_idx: grid indices at which the full spectrum of E + F(t) is also
    returned (for the JSON-driven histogram figure)."""
    Et, Ft0, at, bt0 = ops
    dim = len(lam)
    out = dict(m1=[], m2=[], m3=[], w4=[], W1_arc=[], W1_bin=[])
    snaps = {}
    for i, t in enumerate(ts):
        Ft = heisenberg(Ft0, lam, t)
        M = Et @ Ft                      # alternating-word workhorse
        M2 = M @ M
        out["m1"].append(np.trace(M).real / dim)
        out["m2"].append(np.trace(M2).real / dim)
        out["m3"].append(np.sum(M2 * M.T).real / dim)      # Tr(M^3)
        bt = heisenberg(bt0, lam, t)
        Mo = at @ bt
        out["w4"].append(np.sum(Mo * Mo.T).real / dim)     # Tr(Mo^2)
        eigs = np.linalg.eigvalsh(Et + Ft)
        out["W1_arc"].append(w1_distance(eigs, arcsine_cdf))
        out["W1_bin"].append(w1_distance(eigs, binomial_cdf))
        if i in snap_idx:
            snaps[i] = eigs
    out = {k: np.array(v) for k, v in out.items()}
    out["snaps"] = snaps
    return out


# ---------------------------------------------------------------------------
# disorder-averaged point
# ---------------------------------------------------------------------------
_HIST_EDGES = np.linspace(-2.02, 2.02, 82)


def run_point(N: int, p: int, n_inst: int, seed: int,
              tJ_grid=TJ_GRID, snap_idx=None) -> dict:
    scriptJ = np.sqrt(2.0) * p
    ts = np.array(tJ_grid) / scriptJ
    if snap_idx is None:      # spectra kept at t = 0, t ~ t_x (1.25/J), t_max
        mid = int(np.argmin(np.abs(np.array(tJ_grid) - 1.25)))
        snap_idx = (0, mid, len(tJ_grid) - 1)
    rng = np.random.default_rng(seed)
    acc = {}
    hists = {i: np.zeros(len(_HIST_EDGES) - 1) for i in snap_idx}
    t0 = time.time()
    for _ in range(n_inst):
        lam, ops = instance_operators(N, p, rng)
        r = measure_instance(lam, ops, ts, snap_idx)
        for i, eigs in r.pop("snaps").items():
            hists[i] += np.histogram(eigs, bins=_HIST_EDGES)[0]
        for k, v in r.items():
            acc.setdefault(k, []).append(v)
    sem = 1.0 / np.sqrt(max(n_inst - 1, 1))
    out = dict(N=N, p=p, n_inst=n_inst, seed=seed, scriptJ=scriptJ,
               dim=1 << (N // 2), tJ=list(tJ_grid),
               runtime_s=round(time.time() - t0, 1),
               hist_edges=_HIST_EDGES.tolist(),
               spec_hists={str(tJ_grid[i]): h.tolist()
                           for i, h in hists.items()})
    for k, rows in acc.items():
        A = np.vstack(rows)
        out[k] = A.mean(axis=0).tolist()
        out[k + "_sem"] = (A.std(axis=0, ddof=1) * sem).tolist()

    # verdict scalars ------------------------------------------------------
    # late-time plateaus: mean over the last 3 grid points
    for k in ("W1_arc", "W1_bin", "m2", "w4"):
        out[f"plateau_{k}"] = float(np.mean(np.abs(out[k][-3:])))
    # crossover time t_x (in tJ units): first sign change of W1_arc - W1_bin,
    # linearly interpolated
    d = np.array(out["W1_arc"]) - np.array(out["W1_bin"])
    tJ = np.array(out["tJ"])
    out["t_x_J"] = None
    for i in range(1, len(d)):
        if d[i - 1] > 0.0 >= d[i]:
            frac = d[i - 1] / (d[i - 1] - d[i])
            out["t_x_J"] = float(tJ[i - 1] + frac * (tJ[i] - tJ[i - 1]))
            break
    return out


def main(quick: bool = False) -> dict:
    if quick:
        points = [dict(N=12, p=4, n_inst=8, seed=270)]
    else:
        points = [dict(N=12, p=4, n_inst=96, seed=270),
                  dict(N=16, p=4, n_inst=48, seed=271),
                  dict(N=20, p=4, n_inst=24, seed=272),
                  dict(N=24, p=4, n_inst=6, seed=273)]
    results = []
    for kw in points:
        print(f"running {kw} ...", flush=True)
        r = run_point(**kw)
        results.append(r)
        print(f"  done in {r['runtime_s']} s;  plateau W1_arc = "
              f"{r['plateau_W1_arc']:.4f}, W1_bin = {r['plateau_W1_bin']:.4f}, "
              f"|m2| = {r['plateau_m2']:.4f}, |w4| = {r['plateau_w4']:.4f}, "
              f"t_x J = {r['t_x_J']}", flush=True)
    payload = dict(
        provenance="freeness.py -- rung 27 (CK 2607.13665 Sec. 4): "
                   "classical->free crossover of (E, F(t)) at T_B = infinity; "
                   "W1 distances of spec(E+F(t)) to binomial (classical) and "
                   "arcsine (free) laws; alternating-word battery; "
                   "time in 1/scriptJ units, scriptJ = sqrt(2) p",
        w1_binomial_vs_arcsine=W1_BINOMIAL_VS_ARCSINE,
        points=results)
    if not quick:
        path = Path(__file__).resolve().parent / "freeness.json"
        path.write_text(json.dumps(payload, indent=1))
        print(f"wrote {path}")
    return payload


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
