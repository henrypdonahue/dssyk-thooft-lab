#!/usr/bin/env python3
r"""
The large-q ladder kernel, DERIVED: the skeleton behind the M^e/M^o sets.
=========================================================================

largeq_anchor.py TRANSCRIBES the quantization sets M^e/M^o from
Choi-Mezei-Sarosi (CMS) Eq. (3.3).  This module DERIVES them in-repo,
from the Maldacena-Stanford (MS, 1604.07818) finite-temperature large-q
ladder-kernel eigenproblem, with no root-finding on the CMS conditions
anywhere in the pipeline.  The derived and transcribed sets are then
cross-checked against each other (largeq_kernel.json, test_largeq_kernel).

WHY.  road_map.md Part III names the O(1/q) residues on the
exact-in-temperature skeleton as the open task of the lambda -> 0 corner.
Any perturbation theory on that skeleton needs the skeleton itself as a
derived object -- an eigenproblem one can perturb -- not as a transcribed
root condition.  This module is Phase 0 of that program: the LEADING-order
skeleton, derived.  The O(1/q) residues on it remain open; what is added
here is the coded substrate plus the derived-vs-transcribed cross-check.

THE DERIVATION CHAIN (every equation fetched and read from the sources)

* MS (2.18)-(2.19): the finite-T large-q classical solution
  e^{g(tau)} = [cos(pi v/2) / cos(pi v (1/2 - |tau|/beta))]^2, with
  beta scriptJ = pi v / cos(pi v/2).  All temperature dependence is in v.

* MS (3.112): the large-q symmetric-kernel eigenvalue equation
  K-tilde Psi = k Psi, side rails sgn/2, rung e^g.  Applying
  d_theta3 d_theta4 e^{-g/2} turns it into a differential equation.
  With the eigenvalue parameterized as k = 2/[h(h-1)] and a Fourier
  ansatz in the center-of-mass angle y, MS (3.113)-(3.114) give the
  EXACT eigenproblem at any coupling (beta = 2pi units):

      Psi(theta_1, theta_2) = e^{-i n y} psi_n(x) / sin(x-tilde/2),
      x-tilde = v x + (1-v) pi,
      [ n^2 + 4 d_x^2 - h(h-1) v^2 / sin^2(x-tilde/2) ] psi_n(x) = 0,

  with Dirichlet boundary condition psi(x = 0) = 0 and the parity
  constraint from MS (3.115)-(3.116): psi symmetric about x = pi for
  even n, antisymmetric for odd n.

* The mode condition is k(h) = 1 (ladder poles at 1/(1-k); MS Sec. 3.5:
  "the dimensions are determined by k_c(h_m) = 1").  Setting k = 1 fixes
  h(h-1) = 2 in the potential.  The REMAINING spectral parameter is the
  center-of-mass frequency: in the antipodal variable z = x-tilde - pi,
  the eigenproblem becomes a Schroedinger problem on z in (-v pi, v pi),

      [ - d^2/dz^2 + 1/(2 cos^2(z/2)) ] psi = (m / 2v)^2 psi,
      psi(+-v pi) = 0,             m = v n-tilde,  n-tilde = n / v,

  and the Dirichlet spectrum, split by parity about z = 0 (even -> M^e,
  odd -> M^o), IS the would-be tower.  This is the same object CMS reach:
  their (3.1)-(3.2) eigenfunctions are the MS (3.117)-(3.118)
  hypergeometrics in elementary form, and their (3.3) root conditions are
  the Dirichlet condition evaluated on the closed forms.  (Fixing integer
  n and quantizing h by the boundary condition is the same spectrum read
  along the other axis: k(h(m)) = 1 at m = n.)

* Equivalent PDE route: CMS (2.9)-(2.10) [= MS Sec. 3.4] invert the
  kernel to -2/(v^2 P-tilde^2) d_theta1 d_theta2 with
  P-tilde = 1/sin(x-tilde/2); separation of variables on
  [-d_y^2/4 + v^2 d_xtilde^2 - v^2 P-tilde^2/2] gives the identical
  x-problem.  One skeleton, two source presentations.

DISCRETIZATION (the derived object)

Sine-Galerkin on the Dirichlet interval: basis
sin(k pi (z + L) / 2L) / sqrt(L), L = v pi, k = 1..K.  The kinetic term
is diagonal, the potential matrix is assembled by Gauss-Legendre
quadrature (the potential is analytic on the closed interval for every
v < 1), and the operator block-diagonalizes exactly by parity: odd k are
even about z = 0 (-> M^e), even k are odd (-> M^o).  scipy.linalg.eigh
per block; m = 2 v sqrt(eigenvalue).  Convergence is checked by grid
doubling K -> 2K (largeq_kernel.json records the shifts).

THE SPECIAL POLES, DERIVED

CMS (3.13): after the Sommerfeld-Watson deformation the 4-point function
has Regge poles ONLY at m = 0 (Hamiltonian exchange, the entire TOC) and
m = +-v (the scramblon, OTOC only) -- both OFF the quantization sets.
The kernel derivation locates both as exact closed-form solutions of the
same ODE that fail the Dirichlet condition:

    psi(z) = tan(z/2):  eigenvalue 0    => m = 0    (odd slot),
    psi(z) = sec(z/2):  eigenvalue 1/4  => m = +-v  (even slot),

while the closed-form eigenfunctions DEGENERATE at those slots
(psi^e_{m=0} == 0 and psi^o_{m=v} == 0 identically -- CMS's "triviality"
remark below (3.3)).  In particular m = v solves the M^o root condition
FORMALLY but is excluded because its eigenfunction vanishes; the
transcription in largeq_anchor returns roots m > v only, which this
module verifies is the correct exclusion, not an artifact.

CONVENTIONS.  Everything here is in the source papers' own variables
(beta = 2pi, temperature parameter v); nothing depends on the repo bridge
scriptJ = sqrt(2) p or on the chord lambda convention.  To enter from a
physical beta*scriptJ, use largeq_anchor.v_of_betaJ; the bridge is not
re-derived here (road_map.md "Convention landmine").

HONEST SCOPE.  Leading order in 1/q only.  The O(1/q) residues organized
on this skeleton -- the actual Part III task -- are NOT computed here.
The deliverable is the skeleton as a derived, perturbable eigenproblem,
gated against the transcription across v in {0.05, ..., 0.95}.

Output: largeq_kernel.json.  Runtime: seconds (unloaded machine).
"""

from __future__ import annotations

import json
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent

# k = 2/[h(h-1)] = 1  <=>  h(h-1) = 2: the k(h) = 1 condition, fixed once
HH1_AT_K_EQUALS_1 = 2.0

V_SCAN = tuple(round(0.05 * i, 2) for i in range(1, 20))   # 0.05 .. 0.95
M_MAX = 12.0
K_DEFAULT = 400          # sine modes; v = 0.95 needs the headroom (see json)


# ---------------------------------------------------------------------------
# the derived eigenproblem [MS (3.113)-(3.114) at k = 1, antipodal variable]
# ---------------------------------------------------------------------------
def potential(z, hh1: float = HH1_AT_K_EQUALS_1):
    """W(z) = h(h-1)/(4 cos^2(z/2)): MS (3.114)'s h(h-1) v^2/sin^2(xt/2)
    term in the variable z = x-tilde - pi, after dividing by 4 v^2.  At
    k = 1, h(h-1) = 2 and W = 1/(2 cos^2(z/2))."""
    return hh1 / (4.0 * np.cos(np.asarray(z, float) / 2.0) ** 2)


@lru_cache(maxsize=8)
def _gauss(nq: int):
    t, w = leggauss(nq)
    return t, w


def _parity_blocks(v: float, K: int):
    """Sine-Galerkin blocks of -d^2/dz^2 + W(z) on z in (-L, L), L = v pi,
    Dirichlet ends.  Basis e_k = sin(k pi (z+L)/2L)/sqrt(L); odd k are even
    about z = 0, even k are odd (W is even, so the blocks are exact).
    Returns dict parity -> (H, ks)."""
    L = v * np.pi
    nq = 2 * K + 96
    t, w = _gauss(nq)                      # z = L t, all L factors cancel
    alpha = np.pi * (t + 1.0) / 2.0
    W = potential(L * t)
    r = np.arange(2 * K + 1)
    C = np.cos(np.outer(r, alpha)) @ (w * W)   # C_r = int cos(r alpha) W dt
    out = {}
    for name, rem in (("e", 1), ("o", 0)):
        ks = np.arange(1, K + 1)
        ks = ks[ks % 2 == rem]
        H = np.diag((ks * np.pi / (2.0 * L)) ** 2).astype(float)
        H += 0.5 * (C[np.abs(ks[:, None] - ks[None, :])]
                    - C[ks[:, None] + ks[None, :]])
        out[name] = (H, ks)
    return out


def derived_quantization_sets(v: float, m_max: float = M_MAX,
                              K: int = K_DEFAULT):
    """The DERIVED M^e/M^o: Dirichlet spectrum of the k = 1 large-q kernel
    eigenproblem, m = 2 v sqrt(eigenvalue), split by parity.  Same
    signature contract as largeq_anchor.quantization_sets, which this must
    (and does) reproduce -- but derived, not transcribed."""
    blocks = _parity_blocks(v, K)
    sets = {}
    for name, (H, _) in blocks.items():
        lam = eigh(H, eigvals_only=True)
        m = 2.0 * v * np.sqrt(np.maximum(lam, 0.0))
        sets[name] = [float(x) for x in m[m <= m_max]]
    return sets["e"], sets["o"]


def derived_modes(v: float, m_max: float = M_MAX, K: int = K_DEFAULT):
    """Roots AND eigenfunctions: returns (z, wz, modes) with z the
    Gauss-Legendre nodes on (-v pi, v pi), wz the matching weights for
    int dz, and modes a list of dicts (m, parity, psi) sorted by m.  The
    psi samples are EXACTLY unit L2 norm (orthonormal sine basis, so the
    eigh vector norm is the continuum norm); sign fixed to rise from the
    left end."""
    L = v * np.pi
    nq = 2 * K + 96
    t, w = _gauss(nq)
    z, wz = L * t, L * w
    blocks = _parity_blocks(v, K)
    modes = []
    for name, (H, ks) in blocks.items():
        lam, vec = eigh(H)
        m = 2.0 * v * np.sqrt(np.maximum(lam, 0.0))
        basis = np.sin(np.outer(z + L, ks * np.pi / (2.0 * L))) / np.sqrt(L)
        for i in np.nonzero(m <= m_max)[0]:
            psi = basis @ vec[:, i]
            i0 = int(np.argmax(np.abs(psi) > 1e-3 * np.max(np.abs(psi))))
            if psi[i0] < 0:
                psi = -psi
            modes.append(dict(m=float(m[i]), parity=name, psi=psi))
    modes.sort(key=lambda d: d["m"])
    return z, wz, modes


# ---------------------------------------------------------------------------
# CMS closed forms [their (3.2), (3.4)] -- cross-check targets, never inputs
# ---------------------------------------------------------------------------
def cms_psi_e(m: float, v: float, x):
    """CMS Eq. (3.2), even branch:
    psi^e_m(x) = (m/v) cos(m(pi-x)/2) + sin(m(pi-x)/2) tan(v(pi-x)/2)."""
    a = 0.5 * (np.pi - np.asarray(x, float))
    return m * np.cos(m * a) / v + np.sin(m * a) * np.tan(v * a)


def cms_psi_o(m: float, v: float, x):
    """CMS Eq. (3.2), odd branch:
    psi^o_m(x) = (m/v) sin(m(pi-x)/2) - cos(m(pi-x)/2) tan(v(pi-x)/2)."""
    a = 0.5 * (np.pi - np.asarray(x, float))
    return m * np.sin(m * a) / v - np.cos(m * a) * np.tan(v * a)


def cms_norm(m: float, v: float, parity: str) -> float:
    """CMS Eq. (3.4): N^{e,o}_m = (m-v)(m+v)(m pi -+ sin(m pi))/(m v^2),
    the closed-form value of int_0^{2pi} psi(x)^2 dx at a quantized m."""
    s = -1.0 if parity == "e" else 1.0
    return (m - v) * (m + v) * (m * np.pi + s * np.sin(m * np.pi)) / (m * v ** 2)


# ---------------------------------------------------------------------------
# the special poles m = 0, +-v [CMS (3.13)]: exact non-Dirichlet solutions
# ---------------------------------------------------------------------------
def zero_mode(z):
    """psi(z) = tan(z/2): EXACT solution of -psi'' + W psi = lam psi at
    lam = 0, i.e. m = 0 (Hamiltonian exchange slot).  Odd parity; violates
    the Dirichlet condition, so it is NOT in M^o.  It is CMS's psi^o_{m=0}
    (their psi^e_{m=0} vanishes identically)."""
    return np.tan(np.asarray(z, float) / 2.0)


def scramblon_mode(z):
    """psi(z) = sec(z/2): EXACT solution of -psi'' + W psi = lam psi at
    lam = 1/4, i.e. m = v (the scramblon slot; n-tilde = +-1 are the
    would-be SL(2) modes).  Even parity; violates the Dirichlet condition,
    so it is NOT in M^e.  It is CMS's psi^e_{m=v} (their psi^o_{m=v}
    vanishes identically -- which is why m = v solves the M^o root
    condition formally yet carries no mode)."""
    return 1.0 / np.cos(np.asarray(z, float) / 2.0)


def ode_residual(psi_fn, lam: float, v: float, n_pts: int = 41) -> float:
    """max |(-psi'' + W psi - lam psi)| / max |lam psi + psi| on an interior
    grid, psi'' by 4th-order central differences (independent of any
    hand-derived derivative, so the closed-form identities are tested
    against the ODE itself, not against themselves)."""
    L = v * np.pi
    z = np.linspace(-0.8 * L, 0.8 * L, n_pts)
    h = 1e-3
    d2 = (-psi_fn(z + 2 * h) + 16 * psi_fn(z + h) - 30 * psi_fn(z)
          + 16 * psi_fn(z - h) - psi_fn(z - 2 * h)) / (12.0 * h ** 2)
    res = -d2 + potential(z) * psi_fn(z) - lam * psi_fn(z)
    scale = np.max(np.abs(psi_fn(z))) * max(abs(lam), 1.0)
    return float(np.max(np.abs(res)) / scale)


# ---------------------------------------------------------------------------
# derived vs transcribed
# ---------------------------------------------------------------------------
def _gate_row(v: float, derived_pair, m_max: float = M_MAX):
    """One row of the gate: derived (this module) against transcribed
    (largeq_anchor.quantization_sets) at temperature parameter v.  Roots
    within 1e-6 of m_max are dropped on both sides (truncation guard).
    Returns the row dict; raises on a count mismatch, which would be a
    FINDING, not a tolerance issue."""
    from largeq_anchor import quantization_sets
    d_e, d_o = derived_pair
    t_e, t_o = quantization_sets(v, m_max=m_max)
    cut = m_max - 1e-6
    row = dict(v=v)
    for name, dd, tt in (("e", d_e, t_e), ("o", d_o, t_o)):
        dd = [m for m in dd if m < cut]
        tt = [m for m in tt if m < cut]
        if len(dd) != len(tt):
            raise AssertionError(
                f"count mismatch at v={v}, branch {name}: "
                f"derived {len(dd)} vs transcribed {len(tt)}")
        rel = np.abs(np.array(dd) - np.array(tt)) / np.array(tt)
        row[f"derived_{name}"] = dd
        row[f"transcribed_{name}"] = tt
        row[f"max_rel_dev_{name}"] = float(np.max(rel))
    return row


def compare_at_v(v: float, K: int = K_DEFAULT, m_max: float = M_MAX):
    """Convenience wrapper: derive at K, then gate against the
    transcription."""
    return _gate_row(v, derived_quantization_sets(v, m_max=m_max, K=K),
                     m_max=m_max)


def _shift_between(pair1, pair2, m_max: float = M_MAX) -> float:
    """max |m(K) - m(2K)| over both branches: the spectral convergence
    measure the JSON records per v."""
    cut = m_max - 1e-6
    out = 0.0
    for a, b in zip(pair1, pair2):
        a = [m for m in a if m < cut]
        b = [m for m in b if m < cut]
        n = min(len(a), len(b))
        out = max(out, float(np.max(np.abs(np.array(a[:n]) - np.array(b[:n])))))
    return out


def doubling_shift(v: float, K: int = K_DEFAULT, m_max: float = M_MAX) -> float:
    """Convenience wrapper: derive at K and 2K, return the max root
    shift."""
    return _shift_between(derived_quantization_sets(v, m_max=m_max, K=K),
                          derived_quantization_sets(v, m_max=m_max, K=2 * K),
                          m_max=m_max)


def eigenfunction_check(v: float, K: int = K_DEFAULT, n_roots: int = 2):
    """Optional cross-check (CMS Sec. 3.2): numerical eigenvectors against
    the CMS (3.2) closed forms at the first n_roots of each branch, plus
    the Gauss-Legendre norm int_0^{2pi} psi^2 dx against the CMS (3.4)
    closed form.  Returns a list of dicts (v, m, parity, max_abs_dev,
    norm_rel_dev)."""
    z, wz, modes = derived_modes(v, K=K)
    x = np.pi + z / v                       # z = v (x - pi)
    rows, seen = [], {"e": 0, "o": 0}
    for md in modes:
        p = md["parity"]
        if seen[p] >= n_roots:
            continue
        seen[p] += 1
        fn = cms_psi_e if p == "e" else cms_psi_o
        ref = fn(md["m"], v, x)
        nrm2 = float(wz @ ref ** 2)         # int_{-L}^{L} psi(z)^2 dz
        ref_n = ref / np.sqrt(nrm2)
        if float(md["psi"] @ (wz * ref_n)) < 0:
            ref_n = -ref_n
        dev = float(np.max(np.abs(md["psi"] - ref_n)))
        # int_0^{2pi} psi^2 dx = (1/v) int psi(z)^2 dz vs CMS (3.4)
        norm_x = nrm2 / v
        rows.append(dict(v=v, m=md["m"], parity=p, max_abs_dev=dev,
                         norm_rel_dev=float(abs(norm_x / cms_norm(md["m"], v, p)
                                                - 1.0))))
    return rows


# ---------------------------------------------------------------------------
# producer
# ---------------------------------------------------------------------------
def run(out_path=None, K: int = K_DEFAULT) -> dict:
    t0 = time.time()
    per_v, shifts = [], []
    for v in V_SCAN:
        d1 = derived_quantization_sets(v, K=K)
        d2 = derived_quantization_sets(v, K=2 * K)
        row = _gate_row(v, d1)
        row["doubling_shift"] = _shift_between(d1, d2)
        shifts.append(row["doubling_shift"])
        per_v.append(row)
    max_dev = max(max(r["max_rel_dev_e"], r["max_rel_dev_o"]) for r in per_v)

    eig_rows = []
    for v in (0.3, 0.6, 0.9):
        eig_rows += eigenfunction_check(v, K=K)

    v_sp = 0.6
    special = dict(
        residual_zero_mode=ode_residual(zero_mode, 0.0, v_sp),
        residual_scramblon=ode_residual(scramblon_mode, 0.25, v_sp),
        psi_e_trivial_at_m0=float(np.max(np.abs(
            cms_psi_e(0.0, v_sp, np.linspace(0.3, 5.9, 41))))),
        psi_o_trivial_at_mv=float(np.max(np.abs(
            cms_psi_o(v_sp, v_sp, np.linspace(0.3, 5.9, 41))))),
        min_root_over_v=min(min(r["derived_e"][0], r["derived_o"][0]) / r["v"]
                            for r in per_v),
        note="m = 0 and m = +-v are exact non-Dirichlet solutions (the CMS "
             "(3.13) Regge-pole slots); the quantized tower stays at m > v "
             "in every row, so the transcription's m > v exclusion is "
             "derived, not assumed")

    # CMS's one printed example, v = 3/5 (text below their Eq. (3.3))
    ex_e, ex_o = derived_quantization_sets(0.6, m_max=10.5, K=K)
    printed = dict(v=0.6,
                   Me_2dp=[round(m, 2) for m in ex_e[:5]],
                   Mo_2dp=[round(m, 2) for m in ex_o[:5]],
                   cms_printed_Me=[1.35, 3.16, 5.10, 7.07, 9.06],
                   cms_printed_Mo=[2.23, 4.13, 6.09, 8.06, 10.05])

    payload = dict(
        provenance="largeq_kernel.py -- the leading finite-T large-q ladder "
                   "kernel eigenproblem of MS 1604.07818 (3.112)-(3.118) at "
                   "k = 2/[h(h-1)] = 1, discretized by parity-split "
                   "sine-Galerkin + eigh, derived quantization sets gated "
                   "against the CMS 1912.00004 (3.3) transcription in "
                   "largeq_anchor.py; scope: leading order ONLY, the O(1/q) "
                   "residues on this skeleton remain the open Part III task",
        method=dict(basis="sine (Dirichlet) on z in (-v pi, v pi)",
                    K=K, K_doubled=2 * K, quadrature=f"Gauss-Legendre {2*K+96}",
                    eigensolver="scipy.linalg.eigh per parity block",
                    m_from_eigenvalue="m = 2 v sqrt(lambda)"),
        v_scan=list(V_SCAN), m_max=M_MAX,
        per_v=per_v,
        achieved=dict(max_rel_dev_overall=max_dev,
                      max_doubling_shift=max(shifts),
                      target_rel=1e-8),
        eigenfunction_check=eig_rows,
        special_poles=special,
        cms_printed_example=printed,
        runtime_s=round(time.time() - t0, 1))
    out = Path(out_path) if out_path else HERE / "largeq_kernel.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"derived vs transcribed over v in {{{V_SCAN[0]} .. {V_SCAN[-1]}}}: "
          f"max rel dev {max_dev:.2e} (target 1e-8)")
    print(f"grid doubling K = {K} -> {2*K}: max root shift "
          f"{max(shifts):.2e}")
    ed = max(r["max_abs_dev"] for r in eig_rows)
    en = max(r["norm_rel_dev"] for r in eig_rows)
    print(f"eigenfunctions vs CMS (3.2) closed forms: max dev {ed:.2e}; "
          f"norms vs CMS (3.4): max rel dev {en:.2e}")
    print(f"special poles: ODE residuals zero-mode {special['residual_zero_mode']:.1e}, "
          f"scramblon {special['residual_scramblon']:.1e}; "
          f"min root / v = {special['min_root_over_v']:.3f} (> 1)")
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


if __name__ == "__main__":
    run()
