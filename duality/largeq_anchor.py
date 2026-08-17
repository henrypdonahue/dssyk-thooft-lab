#!/usr/bin/env python3
r"""
The large-q (lambda -> 0) anchor: verified closed forms and their structure.
============================================================================

Rung 17 of road_map.md.  Every formula here was transcribed from the source
papers with equation-level verification (see the provenance block at the
bottom) and is checked against this repo's ED bench
(syk-self-averaging/largeq_bench.json) in test_largeq.py.

Conventions bridge (exact, no approximation):
    repo (psi^2 = 1, Var J = p!/N^{p-1})  <->  Maldacena-Stanford 1604.07818
    (psi^2 = 1/2, <j^2> = J^2 (q-1)!/N^{q-1}, script-J^2 = q J^2 2^{1-q}):
        rescaling psi = sqrt(2) psi_MS forces  J_MS^2 = q 2^q,  hence
        script-J^2 = 2 q^2  ->  script-J = sqrt(2) q  in repo units.
    (Streicher 1911.10171 and Choi-Mezei-Sarosi 1912.00004 use the same
    script-J; all three papers' v solve the same equation.)

The structural result this module encodes (the rung-17 finding):

  * The large-q singlet 4-point function is known in closed form at ALL
    temperatures (Streicher 3.3/3.10; CMS 3.16/3.20), with every
    temperature-dependence absorbed into v.
  * Its analytic structure in the Sommerfeld-Watson plane is MEROMORPHIC
    with contributions ONLY from m = 0 (Hamiltonian exchange; the entire TOC)
    and m = +-v (the scramblon, OTOC only).  The would-be bilinear tower --
    the quantization sets M^e/M^o below, the all-temperature analog of
    k(h) = 1 -- CANCELS at leading order, and its couplings are
    O(1/q^2) (Gross-Rosenhaus 1702.08016 Eq. 2.15).
  * Therefore: THE lambda -> 0 ANCHOR OF THE 'T HOOFT MATCH CANNOT COME FROM
    THE LEADING LARGE-q FOUR-POINT FUNCTION.  The meson tower first appears
    in the O(1/q) corrections, organized on the exact-in-betaJ skeleton
    (the M^e/M^o sets) -- whose residues at finite temperature are
    unpublished.  Computing them IS rung 17 proper.
  * Encouragingly for the correspondence: M^e/M^o alternate (even/odd), the
    exact CP-alternation slot of the 't Hooft tower; and at T = infinity the
    sets degenerate to integers (in 2pi/beta units) -- consistent with
    BOOST.md: physical masses are not boost-frequency peak positions.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq


# ---------------------------------------------------------------------------
# the interpolation parameter v(beta J)  [MS (2.19); Streicher (2.12)]
# ---------------------------------------------------------------------------
def v_of_betaJ(betaJ: float) -> float:
    """Solve beta*scriptJ = pi v / cos(pi v/2), v in (0, 1).
    v -> 0 at weak coupling / high T; v -> 1 in the conformal limit."""
    if betaJ <= 0:
        return 0.0
    f = lambda v: np.pi * v / np.cos(np.pi * v / 2.0) - betaJ
    return brentq(f, 0.0, 1.0 - 1e-12)


# ---------------------------------------------------------------------------
# two-point function [MS (2.18); T = inf real-time limit]
# ---------------------------------------------------------------------------
def g_thermal_euclidean(tau: float, beta: float, scriptJ: float) -> float:
    """e^{g(tau)} of MS (2.18), 0 <= tau <= beta."""
    v = v_of_betaJ(beta * scriptJ)
    num = np.cos(np.pi * v / 2.0)
    den = np.cos(np.pi * v * (0.5 - abs(tau) / beta))
    return (num / den) ** 2


def G_infT_repo(t, p: int):
    """T = infinity real-time two-point function in REPO conventions
    (psi^2 = 1, Var J = p!/N^{p-1}):

        G(t) = [sech(scriptJ t)]^{2/p},   scriptJ = sqrt(2) p.

    (The beta J -> 0 limit of MS (2.18)-(2.19) continued to real time; the
    exponentiated large-q resummation G = e^{g/q} at leading order.)
    Corrections are O(1/p) and O(lambda = p^2/N)."""
    scriptJ = np.sqrt(2.0) * p
    return (1.0 / np.cosh(scriptJ * np.asarray(t, float))) ** (2.0 / p)


# ---------------------------------------------------------------------------
# the would-be tower: quantization sets [CMS (3.3)] and their limits
# ---------------------------------------------------------------------------
def quantization_sets(v: float, m_max: float = 12.0):
    """The all-temperature analog of k(h) = 1 [CMS Eq. (3.3)]:

        M^e:  m cot(m pi/2) + v tan(v pi/2) = 0
        M^o:  m tan(m pi/2) - v tan(v pi/2) = 0

    Returns (M_even, M_odd) with roots m > v up to m_max.  As v -> 1
    (conformal) the union approaches the integers > 1; as v -> 0 (T = inf)
    M^e -> odd integers and M^o -> even integers."""
    c = v * np.tan(np.pi * v / 2.0)
    fe = lambda m: m / np.tan(np.pi * m / 2.0) + c
    fo = lambda m: m * np.tan(np.pi * m / 2.0) - c
    roots_e, roots_o = [], []
    grid = np.linspace(max(v, 1e-9) + 1e-9, m_max, 24001)
    for f, out in [(fe, roots_e), (fo, roots_o)]:
        vals = np.array([f(m) for m in grid])
        for i in range(len(grid) - 1):
            a, b = vals[i], vals[i + 1]
            # skip the poles of tan/cot (sign flips with huge magnitude)
            if not (np.isfinite(a) and np.isfinite(b)):
                continue
            if abs(a) > 1e3 or abs(b) > 1e3:
                continue
            if a == 0.0:
                out.append(float(grid[i]))
            elif a * b < 0:
                out.append(float(brentq(f, grid[i], grid[i + 1])))
    return roots_e, roots_o


# ---------------------------------------------------------------------------
# the conformal-side tower the O(1/q) corrections must connect to
# [Gross-Rosenhaus 1702.08016 (2.14)-(2.15); MS (3.94), Delta -> 0]
# ---------------------------------------------------------------------------
def gr_tower(q: float, n_max: int = 8):
    """h_n = 2n + 1 + 2 eps_n with eps_n = (1/q)(2n^2+n+1)/(2n^2+n-1), and
    the O(1/q^2) squared OPE couplings c_n^2 [GR Eq. (2.15)] -- the
    conformal-limit data whose finite-temperature continuation is the
    unpublished content of rung 17."""
    from math import gamma, sqrt, pi
    rows = []
    for n in range(1, n_max + 1):
        eps = (2 * n * n + n + 1) / (2 * n * n + n - 1) / q
        h = 2 * n + 1 + 2 * eps
        c2 = (eps ** 2 * n * (1 + 2 * n)
              / ((n * (1 + 2 * n) + 1) * (n * (1 + 2 * n) - 1))
              * sqrt(pi) * gamma(2 * n + 1)
              / (gamma(2 * n + 0.5) * 2 ** (4 * n - 2)))
        rows.append(dict(n=n, h=h, eps=eps, c2=c2,
                         m2_bulk=2 * n * (2 * n + 1)))
    return rows


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------
PROVENANCE = """
Equation-verified sources (arXiv, fetched 2026-07-21):
  1604.07818 (Maldacena-Stanford): (2.2)-(2.3) model+variance; (2.16)
    Liouville + scriptJ def; (2.18)-(2.19) thermal e^g and v; (3.113)-(3.120)
    exact large-q kernel eigenproblem at any coupling; (3.76) k_c = 2/[h(h-1)]
    at q = inf; (3.145) exact q=inf CoM propagator.
  1911.10171 (Streicher): (2.1) conventions (psi^2 = 1 THERE TOO, but
    different Var normalization); (3.3) TOC and (3.10) OTOC closed forms,
    uniformly valid in betaJ; (3.4) 'Hamiltonian + identity exhaust the TOC'.
  1912.00004 (Choi-Mezei-Sarosi): (3.2)-(3.3) eigenfunctions + quantization;
    (3.13) Sommerfeld-Watson deformation, poles only at m = 0, +-v, no cuts
    (their Fig. 1); (3.16)/(3.20) final TOC/OTOC.
  1702.08016 (Gross-Rosenhaus): (2.14)-(2.15) the large-q conformal tower
    h_n and c_n^2 ~ 1/q^2.
The repo <-> MS convention bridge scriptJ = sqrt(2) p is exact and is
validated against ED (<H^2>/dim and G(t)) in test_largeq.py.
"""

if __name__ == "__main__":
    print("v(beta scriptJ):", [(bj, round(v_of_betaJ(bj), 4))
                               for bj in (0.1, 1.0, 10.0, 100.0)])
    print("\nQuantization sets (the would-be tower), showing the even/odd")
    print("alternation (the CP slot) and the T=inf integer degeneration:")
    for v in (0.9, 0.5, 0.1):
        Me, Mo = quantization_sets(v, m_max=8.0)
        print(f"  v={v}:  M^e = {[round(m, 3) for m in Me[:4]]}   "
              f"M^o = {[round(m, 3) for m in Mo[:4]]}")
    print("\nGR conformal tower at q = 8 (the data the O(1/q) continuation")
    print("must connect to; note c_n^2 ~ 1/q^2 -- the tower is turned OFF at")
    print("leading order):")
    for r in gr_tower(8.0, 4):
        print(f"  n={r['n']}: h = {r['h']:.4f}  c_n^2 = {r['c2']:.2e}")
    print(PROVENANCE)
