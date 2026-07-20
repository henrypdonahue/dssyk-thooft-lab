#!/usr/bin/env python3
r"""
Exact disorder statistics of the adjoint bilinear: Eq. (3.28) beyond ED.
========================================================================

The symmetry-violation observable of self_averaging.py is

    B_jk = Tr(psi_j H psi_k H)/dim,      H = sum_c J_c Gamma_c,
    Gamma_c = i^{p/2} psi_{c1}...psi_{cp}   (ascending p-subsets c),
    J_c iid N(0, sigma^2),  sigma^2 = p!/N^{p-1}.

B_jk is QUADRATIC in the Gaussian couplings, so its full disorder statistics
reduce, by Wick's theorem, to Majorana trace combinatorics -- computable in
closed form at ANY (N, p), including the double-scaling regime p ~ sqrt(N),
N ~ 10^3 that instance-level ED (dim 2^{N/2}) can never reach.  This module
carries the module's Eq. (3.28) test into that regime.

Derivation (asserted step-by-step against string enumeration in the tests)
--------------------------------------------------------------------------
The normalized Majorana trace of a monomial vanishes unless every index
appears an even number of times.  With T_cd = Tr(psi_j Gamma_c psi_k
Gamma_d)/dim and j != k, nonvanishing requires the symmetric difference
c XOR d = {j, k}, i.e.

    d = (c \ {j}) u {k}   (j in c, k not in c),      or j <-> k.

Writing c = e u {j}, d = e u {k} (|e| = p-1) and pulling psi_j, psi_k out of
the ordered products gives

    T_cd = i^p sigma_c sigma_d (-1)^{(p-1)(p-2)/2},     T_dc = T_cd,

with sigma's = +-1 reordering signs.  Hence T_cd^2 = T_cd T_dc = +1 for every
admissible pair.  E[B_jk] = sigma^2 sum_c T_cc = 0 (no admissible c), and the
two surviving Wick pairings give, EXACTLY,

    Var(B_jk) = sigma^4 [ sum_{cd} T_cd^2 + sum_{cd} T_cd T_dc ]
              = 4 sigma^4 C(N-2, p-1),

so the typical single-instance symmetry violation is

    rms(B_jk) = 2 (p!/N^{p-1}) sqrt(C(N-2, p-1)).

The same algebra on the diagonal (T_cc = -1 for j in c, +1 otherwise -- i.e.
psi_j Gamma_c psi_j Gamma_c = -+1) gives

    E[B_jj] = sigma^2 [ C(N,p) - 2 C(N-1,p-1) ] = sigma^2 C(N,p) (1 - 2p/N):

the ensemble-mean diagonal VANISHES at p = N/2 -- deriving, rather than just
flagging, the "ratio artifact" self_averaging.py warns about.

What this settles
-----------------
Along the double-scaling line p = sqrt(lambda N):

    ln rms(B_jk) = -(1/4) sqrt(N) ln N (1 + o(1))        (at lambda = 1),

super-exponential in sqrt(N): smaller than e^{-a sqrt N} for EVERY a, exactly
as the paper's Eq. (3.28) bound demands -- now verified by exact computation
in the regime the bound is about, not extrapolated from N <= 18 fits.  At
fixed p (the regime where the 't Hooft comparison lives) the same formula is
the power law rms ~ 2 p!/sqrt((p-1)!) N^{-(p-1)/2}: symmetry emergence is
polynomially slow exactly where the spectrum match is to be made -- the
quantitative content of correction #2's footnote.
"""

from __future__ import annotations

from itertools import combinations
from math import comb, log, sqrt

from syk import coupling_variance
from pauli_strings import string_product, majorana_string

_I_POW = [1, 1j, -1, -1j]


# ---------------------------------------------------------------------------
# closed forms
# ---------------------------------------------------------------------------
def var_offdiag_exact(N: int, p: int) -> float:
    """Var(B_jk) over disorder for j != k:  4 sigma^4 C(N-2, p-1), exact."""
    sigma2 = coupling_variance(N, p)
    return 4.0 * sigma2 ** 2 * comb(N - 2, p - 1)


def rms_offdiag_exact(N: int, p: int) -> float:
    """Typical symmetry violation sqrt(Var(B_jk)) = 2 sigma^2 sqrt(C(N-2,p-1))."""
    return sqrt(var_offdiag_exact(N, p))


def ln_rms_offdiag(N: int, p: int) -> float:
    """ln rms(B_jk) evaluated stably at large (N, p) via lgamma."""
    from math import lgamma
    ln_sigma2 = lgamma(p + 1) - (p - 1) * log(N)
    ln_binom = lgamma(N - 1) - lgamma(p) - lgamma(N - p)   # C(N-2, p-1)
    return log(2.0) + ln_sigma2 + 0.5 * ln_binom


def mean_diag_exact(N: int, p: int) -> float:
    """E[B_jj] = sigma^2 C(N,p) (1 - 2p/N); vanishes at p = N/2 (the derived
    origin of the ratio artifact)."""
    sigma2 = coupling_variance(N, p)
    return sigma2 * comb(N, p) * (1.0 - 2.0 * p / N)


# ---------------------------------------------------------------------------
# brute-force string enumeration (verifies every step of the derivation)
# ---------------------------------------------------------------------------
def _trace_of_product(N: int, index_lists) -> complex:
    """Normalized trace of a product of Majorana monomials, via Pauli strings:
    nonzero only if the total string is the identity, in which case it is the
    accumulated phase."""
    s = (0, 0, 0)
    for idx in index_lists:
        for a in idx:
            s = string_product(s, majorana_string(N, a))
    phase, x, z = s
    return _I_POW[phase] if (x == 0 and z == 0) else 0.0


def var_offdiag_enumerated(N: int, p: int, j: int = 0, k: int = 1) -> float:
    """Var(B_jk) summed literally over ALL pairs (c, d) of p-subsets:
    sigma^4 [sum T_cd^2 + sum T_cd T_dc], each T computed as a Majorana trace.
    Exponentially expensive -- for tests at small (N, p) only."""
    sigma2 = coupling_variance(N, p)
    prefac = (p // 2) % 4     # i^{p/2} per Gamma
    subsets = list(combinations(range(N), p))
    T = {}
    for c in subsets:
        for d in subsets:
            t = _trace_of_product(N, [(j,), c, (k,), d])
            t = (_I_POW[(2 * prefac) % 4]) * t     # the two i^{p/2} factors
            if abs(t) > 1e-12:
                assert abs(t.imag) < 1e-12, "T_cd must be real"
                T[(c, d)] = t.real
    s_sq = sum(t * t for t in T.values())
    s_cross = sum(t * T.get((d, c), 0.0) for (c, d), t in T.items())
    return sigma2 ** 2 * (s_sq + s_cross)


def mean_diag_enumerated(N: int, p: int, j: int = 0) -> float:
    """E[B_jj] = sigma^2 sum_c T_cc, enumerated."""
    sigma2 = coupling_variance(N, p)
    prefac = (p // 2) % 4
    total = 0.0
    for c in combinations(range(N), p):
        t = _trace_of_product(N, [(j,), c, (j,), c])
        t = (_I_POW[(2 * prefac) % 4]) * t
        total += t.real
    return sigma2 * total


# ---------------------------------------------------------------------------
# the double-scaling scan: Eq. (3.28) where ED cannot go
# ---------------------------------------------------------------------------
def double_scaling_scan(lam: float, ns) -> list:
    """Exact ln rms(B_jk) along p = nearest even sqrt(lambda N)."""
    from self_averaging import even_p_near
    rows = []
    for N in ns:
        p = even_p_near((lam * N) ** 0.5)
        if p >= N - 1 or p < 2:
            continue
        rows.append(dict(N=N, p=p, ln_rms=ln_rms_offdiag(N, p)))
    return rows


if __name__ == "__main__":
    print("Exact vs ED (values from results.json where available):")
    import json
    import pathlib
    rj = pathlib.Path(__file__).parent / "results.json"
    measured = {}
    if rj.exists():
        data = json.loads(rj.read_text())
        for block in data.get("C_symmetry_violation", {}).values():
            for r in block["rows"]:
                measured[(r["N"], r["p"])] = (r["rms_off"], r["rms_off_err"])
    print(f"  {'N':>4} {'p':>3} {'exact rms':>12} {'ED rms_off':>12} {'ratio':>7}")
    for (N, p) in [(10, 4), (12, 4), (14, 4), (16, 4), (18, 4),
                   (10, 6), (14, 6), (18, 6)]:
        ex = rms_offdiag_exact(N, p)
        if (N, p) in measured:
            m, err = measured[(N, p)]
            print(f"  {N:>4} {p:>3} {ex:>12.4e} {m:>12.4e} {m/ex:>7.4f}")
        else:
            print(f"  {N:>4} {p:>3} {ex:>12.4e} {'--':>12}")

    print("\nDerived: E[B_jj] = sigma^2 C(N,p)(1 - 2p/N) -> the p = N/2 ratio")
    print("artifact is the exact vanishing of the mean diagonal:")
    for (N, p) in [(12, 4), (12, 6), (16, 8)]:
        print(f"  N={N:>3} p={p}:  E[B_jj] = {mean_diag_exact(N, p):+.4e}"
              + ("   <- p = N/2" if p * 2 == N else ""))

    print("\nDouble-scaling line p = sqrt(N) (lambda = 1), exact -- the regime")
    print("ED cannot reach.  Eq. (3.28) demands ln rms / sqrt(N) -> -infinity:")
    rows = double_scaling_scan(1.0, [16, 36, 64, 100, 196, 400, 784, 1600])
    print(f"  {'N':>5} {'p':>4} {'ln rms(B_jk)':>14} {'ln rms/sqrtN':>13} "
          f"{'-(1/4)sqrtN lnN':>16}")
    for r in rows:
        N = r["N"]
        print(f"  {N:>5} {r['p']:>4} {r['ln_rms']:>14.2f} "
              f"{r['ln_rms']/sqrt(N):>13.3f} {-0.25*sqrt(N)*log(N):>16.2f}")
    print("\n  ln rms / sqrt(N) decreases without bound -> the violation is")
    print("  smaller than e^{-a sqrt N} for EVERY a: Eq. (3.28) holds, with the")
    print("  sharp form  ln rms = -(1/4) sqrt N ln N (1 + o(1))  at lambda = 1.")
