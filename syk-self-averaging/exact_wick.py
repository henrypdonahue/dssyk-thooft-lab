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

Extension: exact disorder statistics of the singlet fourth moment m4
--------------------------------------------------------------------
m4 = Tr(H^4)/dim is QUARTIC in the couplings, so E[m4] is a 3-pairing Wick
sum and Var(m4) an 8th-order one (105 pair partitions) -- still pure Majorana
trace combinatorics, exact at any (N, p).  Two facts drive everything (both
asserted against string enumeration in test_exact_wick_m4.py):

  * exchange sign:  Gamma_c Gamma_d = (-1)^{|c n d|} Gamma_d Gamma_c  for
    p even, so  Tr(Gamma_c Gamma_d Gamma_c Gamma_d)/dim = (-1)^{|c n d|};
  * support: a product of Gamma's has nonzero normalized trace iff the
    symmetric difference of the index sets is empty, and then equals +-1.

With C = C(N,p), sigma^2 = p!/N^{p-1}, and the crossing sum

    K(N,p) = sum_{j=0}^p (-1)^j C(p,j) C(N-p, p-j),

Isserlis' theorem (exact at coincident indices; the diagonal c=d=e=f is the
threefold pairing, i.e. E[J^4] = 3 sigma^4, automatically) gives

    E[m4] = sigma^4 [ 2 C^2 + C K ],       E[m4]/E[m2]^2 = 2 + K/C.

K/C is the chord q-parameter:  K/C -> exp(-lambda_chord) = exp(-2 p^2/N)
(measured: |K/C - q| = O(N^-2) at fixed p; -> e^{-2 lambda} with O(1/p)
error along p = sqrt(lambda N)), so E[m4]/E[m2]^2 -> 2 + q is the chord
count for the 4th moment (2 non-crossing diagrams + 1 crossing, weight q).

Var(m4): of the 105 pairings of the 8 couplings, 9 reproduce E[m4]^2; the
72 with a single cross-pair force both crossing subsets equal (support
argument) and sum to 8 sigma^8 C (2C+K)^2; the 24 four-cross pairings
reduce, by the dihedral symmetry of the trace, to 8 sigma^8 (Z + 2W), where
(|x1 XOR x2| = 2m fixes |x1 n x2| = p - m, whence the sign in W)

    Z = sum_m C(N,2m) [C(2m,m) C(N-2m,p-m)]^2   (# 4-tuples, empty symdiff),
    W = sum_m (-1)^{p-m} C(N,2m) [C(2m,m) C(N-2m,p-m)]^2.

So, EXACTLY and in closed binomial form (no residual pattern sum):

    Var(m4)    = 8 sigma^8 [ C (2C+K)^2 + Z + 2W ],
    relVar(m4) = 8/C  +  8 (Z + 2W) / (C^2 (2C+K)^2).

The first term is exactly 8/C at every (N, p), so relVar(m4) =
4 relVar(m2) x (1 + O(4^p/C)).  Fixed p: relVar(m4) -> 8 p!/N^p (measured
ratio to 8/C: 1.070 at (12,4), 1.00013 at (60,4)).  Along p = sqrt(N) the
correction dies super-exponentially (< 1e-6 by N = 100), so

    ln relVar(m4) = -(1/2) sqrt(N) ln N (1 + o(1)),

the SAME Eq. (3.27) super-exponential law and leading coefficient as
relVar(m2) = 2/C(N,p), with relVar(m4)/relVar(m2) -> 4: the 4th moment
extends the m2 self-averaging anchor at a fixed factor.  Verified: brute
enumeration (shared code: none) at (6,2), (8,2), (6,4), (8,4) to < 1e-9;
ED disorder ensembles at N = 12, 14 (250/200 instances) within 1 sigma.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, product
from math import comb, log, sqrt

from syk import coupling_variance
from pauli_strings import (_I_POW, string_product, majorana_string,
                           monomial_string)


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


def ln_mean_diag(N: int, p: int) -> float:
    """ln E[B_jj] evaluated stably at large (N, p) via lgamma (needs
    p < N/2 so the log argument is positive)."""
    from math import lgamma, log1p
    ln_sigma2 = lgamma(p + 1) - (p - 1) * log(N)
    ln_binom = lgamma(N + 1) - lgamma(p + 1) - lgamma(N - p + 1)
    return ln_sigma2 + ln_binom + log1p(-2.0 * p / N)


# ---------------------------------------------------------------------------
# closed forms: the singlet fourth moment m4 = Tr(H^4)/dim
# ---------------------------------------------------------------------------
def crossing_sum(N: int, p: int) -> int:
    """K(N,p) = sum_j (-1)^j C(p,j) C(N-p,p-j), exact integer: the signed count
    of p-subset pairs by intersection size, with the exchange sign (-1)^|c n d|."""
    return sum((-1) ** j * comb(p, j) * comb(N - p, p - j) for j in range(p + 1))


def crossing_q(N: int, p: int) -> float:
    """Normalized crossing sum K/C(N,p): the finite-(N,p) chord q-parameter.
    -> exp(-lambda_chord) = exp(-2 p^2/N) in the double-scaled limit."""
    return crossing_sum(N, p) / comb(N, p)


def mean_m4_exact(N: int, p: int) -> float:
    """E[m4] = sigma^4 [2 C(N,p)^2 + C(N,p) K(N,p)], exact (Isserlis; the
    degenerate c=d=e=f diagonal with E[J^4] = 3 sigma^4 is included)."""
    sigma2 = coupling_variance(N, p)
    C = comb(N, p)
    return sigma2 ** 2 * C * (2 * C + crossing_sum(N, p))


def quad_null_count(N: int, p: int) -> int:
    """Z(N,p) = # ordered 4-tuples of p-subsets with empty symmetric difference
    = sum_m C(N,2m) [C(2m,m) C(N-2m,p-m)]^2  (m = half the pair symdiff size)."""
    return sum(comb(N, 2 * m) * (comb(2 * m, m) * comb(N - 2 * m, p - m)) ** 2
               for m in range(min(p, N // 2) + 1))


def quad_null_alternating(N: int, p: int) -> int:
    """W(N,p) = sum over the same 4-tuples of (-1)^{|x1 n x2|}; since
    |x1 XOR x2| = 2m fixes |x1 n x2| = p-m, this is the Z sum with (-1)^{p-m}."""
    return sum((-1) ** (p - m) * comb(N, 2 * m)
               * (comb(2 * m, m) * comb(N - 2 * m, p - m)) ** 2
               for m in range(min(p, N // 2) + 1))


def var_m4_exact(N: int, p: int) -> float:
    """Var(m4) = 8 sigma^8 [C (2C+K)^2 + Z + 2W], exact at polynomial cost:
    72 single-cross Wick pairings give the first term, 24 four-cross ones Z+2W."""
    sigma2 = coupling_variance(N, p)
    C = comb(N, p)
    K = crossing_sum(N, p)
    return 8.0 * sigma2 ** 4 * (C * (2 * C + K) ** 2
                                + quad_null_count(N, p)
                                + 2 * quad_null_alternating(N, p))


def _relvar_m4_ratio(N: int, p: int):
    """Exact integer (numerator, denominator) of relVar(m4); sigma^8 cancels."""
    C = comb(N, p)
    K = crossing_sum(N, p)
    num = 8 * (C * (2 * C + K) ** 2
               + quad_null_count(N, p) + 2 * quad_null_alternating(N, p))
    den = C * C * (2 * C + K) ** 2
    return num, den


def relvar_m4_exact(N: int, p: int) -> float:
    """relVar(m4) = 8/C + 8(Z+2W)/(C^2 (2C+K)^2), exact.  Leading term is
    exactly 8/C(N,p) = 4 relVar(m2) at every (N, p); big-int division keeps
    this correct far beyond float range of the intermediates."""
    num, den = _relvar_m4_ratio(N, p)
    return num / den


def ln_relvar_m4(N: int, p: int) -> float:
    """ln relVar(m4) evaluated exactly via big-integer logs (any N, p)."""
    num, den = _relvar_m4_ratio(N, p)
    return log(num) - log(den)


def relvar_m4_double_scaling_scan(lam: float, ns) -> list:
    """Exact ln relVar(m4) along p = nearest even sqrt(lambda N), with the
    m2 comparison 2/C(N,p) -- the Eq. (3.27) law at the 4th moment."""
    from self_averaging import even_p_near
    rows = []
    for N in ns:
        p = even_p_near((lam * N) ** 0.5)
        if p >= N - 1 or p < 2:
            continue
        rows.append(dict(N=N, p=p, ln_relvar_m4=ln_relvar_m4(N, p),
                         ln_relvar_m2=log(2) - log(comb(N, p))))
    return rows


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
# brute-force m4 statistics (shares NO code with the closed forms above:
# literal Pauli-string traces + the Gaussian moment formula E[J^k]=(k-1)!!s^k)
# ---------------------------------------------------------------------------
def _m4_monomials(N: int, p: int, filter_xor: bool = False) -> dict:
    """m4 = Tr(H^4)/dim as a COLLECTED polynomial in the couplings: key = the
    sorted 4-tuple of subset indices (with multiplicity), value = the integer
    sum of literal normalized traces (each +-1, computed by string products).

    filter_xor=True skips tuples whose total Majorana support (XOR of the four
    subset bitmasks) is nonzero -- those traces vanish identically (asserted on
    random samples in the tests); needed to make (8,4) affordable."""
    subsets = list(combinations(range(N), p))
    strings = [monomial_string(N, c) for c in subsets]
    n = len(subsets)
    if filter_xor:
        masks = [sum(1 << a for a in c) for c in subsets]
        by_xor = {}
        for i in range(n):
            for j in range(n):
                by_xor.setdefault(masks[i] ^ masks[j], []).append((i, j))
        tuples = ((i, j, k, l) for pairs in by_xor.values()
                  for (i, j) in pairs for (k, l) in pairs)
    else:
        tuples = product(range(n), repeat=4)
    prefac = _I_POW[(4 * (p // 2)) % 4]      # the four i^{p/2} factors (= +1)
    coeffs = {}
    for (i, j, k, l) in tuples:
        phase, x, z = string_product(string_product(strings[i], strings[j]),
                                     string_product(strings[k], strings[l]))
        if x == 0 and z == 0:
            t = complex(prefac * _I_POW[phase])
            assert t == 1 or t == -1, "empty-support Gamma product must be +-1"
            key = tuple(sorted((i, j, k, l)))
            coeffs[key] = coeffs.get(key, 0) + int(t.real)
    return coeffs


def _gaussian_moment(k: int, sigma2: float) -> float:
    """E[J^k] for J ~ N(0, sigma2): (k-1)!! sigma^k for even k, else 0."""
    if k % 2:
        return 0.0
    df = 1
    for i in range(k - 1, 0, -2):
        df *= i
    return df * sigma2 ** (k // 2)


def _mean_var_m4_enumerated(N: int, p: int, filter_xor: bool = False):
    """(E[m4], Var(m4)) by generic polynomial algebra: square the collected
    polynomial and apply the moment formula variable by variable.  Only pairs
    of monomials with the SAME odd-exponent variable set survive (odd Gaussian
    moments vanish), so the pair loop groups by that signature.  Exponential
    cost in (N, p) -- for verification at small sizes only."""
    sigma2 = coupling_variance(N, p)
    monos = []
    for key, c in _m4_monomials(N, p, filter_xor).items():
        cnt = Counter(key)
        odd = frozenset(v for v, mult in cnt.items() if mult % 2)
        monos.append((cnt, odd, c))
    mean = 0.0
    for cnt, odd, c in monos:
        if not odd:
            mom = 1.0
            for mult in cnt.values():
                mom *= _gaussian_moment(mult, sigma2)
            mean += c * mom
    groups = {}
    for mono in monos:
        groups.setdefault(mono[1], []).append(mono)
    e2 = 0.0
    for group in groups.values():
        for cnt1, _, c1 in group:
            for cnt2, _, c2 in group:
                mom = 1.0
                for mult in (cnt1 + cnt2).values():
                    mom *= _gaussian_moment(mult, sigma2)
                e2 += c1 * c2 * mom
    return mean, e2 - mean * mean


def mean_m4_enumerated(N: int, p: int, filter_xor: bool = False) -> float:
    """E[m4] summed literally over subset 4-tuples (traces via Pauli strings,
    expectations via the Gaussian moment formula -- no Wick pairing taxonomy)."""
    return _mean_var_m4_enumerated(N, p, filter_xor)[0]


def var_m4_enumerated(N: int, p: int, filter_xor: bool = False) -> float:
    """Var(m4) by squaring the literal coupling polynomial of m4."""
    return _mean_var_m4_enumerated(N, p, filter_xor)[1]


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


# ---------------------------------------------------------------------------
# the fixed-p scan: symmetry emergence where the 't Hooft comparison lives
# ---------------------------------------------------------------------------
def fixed_p_scan(p: int, ns) -> list:
    """Exact symmetry violation at FIXED p as lambda = p^2/N -> 0 (rung 16):
    the regime of the spectrum match, NOT the double-scaling regime.  Per N:
    lambda, ln rms(B_jk), the local slope d ln rms / d ln N (asymptote
    -(p-1)/2), and the normalized per-element violation
    ln[rms(B_jk)/E[B_jj]] (asymptote slope -(p+1)/2)."""
    rows = []
    prev = None
    for N in ns:
        lr = ln_rms_offdiag(N, p)
        row = dict(N=N, p=p, lam=p * p / N, ln_rms=lr,
                   ln_ratio=lr - ln_mean_diag(N, p), slope=None)
        if prev is not None:
            row["slope"] = (lr - prev["ln_rms"]) / (log(N) - log(prev["N"]))
        rows.append(row)
        prev = row
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

    from math import exp, factorial
    print("\nFixed p, lambda = p^2/N -> 0 (rung 16) -- the regime where the")
    print("'t Hooft comparison lives.  Emergence is only POLYNOMIAL here:")
    print("rms ~ 2 p!/sqrt((p-1)!) N^{-(p-1)/2}, and the per-element violation")
    print("rms(B_jk)/E[B_jj] ~ 2 p!/sqrt((p-1)!) N^{-(p+1)/2}:")
    for p in (4, 6):
        print(f"  p = {p}  (slope asymptotes: rms -(p-1)/2 = {-(p-1)/2:.1f}, "
              f"ratio -(p+1)/2 = {-(p+1)/2:.1f})")
        print(f"  {'N':>6} {'lambda':>8} {'rms(B_jk)':>11} "
              f"{'d ln rms/d ln N':>16} {'rms/E[B_jj]':>12}")
        for r in fixed_p_scan(p, [20, 40, 80, 160, 320, 640, 1280]):
            slope = "--" if r["slope"] is None else f"{r['slope']:.3f}"
            print(f"  {r['N']:>6} {r['lam']:>8.3f} {exp(r['ln_rms']):>11.3e} "
                  f"{slope:>16} {exp(r['ln_ratio']):>12.3e}")
    pref = 2.0 * factorial(4) / sqrt(factorial(3))
    n_at = {eps: (pref / eps) ** (2.0 / 5.0) for eps in (1e-3, 1e-6, 1e-9)}
    print("\n  The consequence nobody had drawn: along double scaling the")
    print("  violation is super-exponential, but at fixed p = 4 a per-element")
    print("  violation of 1e-3 / 1e-6 / 1e-9 needs N ~ "
          + " / ".join(f"{n_at[e]:.0f}" for e in (1e-3, 1e-6, 1e-9))
          + ",")
    print("  i.e. the emergent U(N) symmetry is quantitatively WEAKEST exactly")
    print("  where the spectrum match is to be made (lambda -> 0 at fixed p).")

    print("\n" + "=" * 70)
    print("m4 extension: E[m4]/E[m2]^2 = 2 + K/C and the chord q-parameter")
    print("=" * 70)
    print(f"  {'N':>6} {'p':>3} {'K/C':>11} {'exp(-2p^2/N)':>13} {'diff':>10}")
    for (N, p) in [(100, 4), (1000, 4), (10000, 4),
                   (100, 10), (400, 20), (1600, 40), (6400, 80)]:
        s = crossing_q(N, p)
        q = exp(-2.0 * p * p / N)
        print(f"  {N:>6} {p:>3} {s:>11.6f} {q:>13.6f} {s - q:>+10.2e}")
    print("  -> K/C -> q = exp(-lambda_chord), lambda_chord = 2p^2/N (NOT p^2/N):")
    print("     O(N^-2) convergence at fixed p, O(1/p) along p = sqrt(lambda N).")

    print("\nrelVar(m4) exact vs its leading term 8/C(N,p) = 4 relVar(m2):")
    print(f"  {'N':>4} {'p':>3} {'relVar(m4)':>13} {'8/C':>13} {'ratio':>8}")
    for (N, p) in [(12, 4), (14, 4), (20, 4), (40, 4), (60, 4), (30, 6)]:
        rv = relvar_m4_exact(N, p)
        lead = 8.0 / comb(N, p)
        print(f"  {N:>4} {p:>3} {rv:>13.4e} {lead:>13.4e} {rv / lead:>8.5f}")

    print("\nDouble-scaling line p = sqrt(N): relVar(m4) vs relVar(m2) = 2/C:")
    rows = relvar_m4_double_scaling_scan(1.0, [16, 64, 100, 400, 1600, 6400])
    print(f"  {'N':>5} {'p':>3} {'ln relVar(m4)':>14} {'/sqrtN':>8} "
          f"{'rv4/rv2':>8}")
    for r in rows:
        N = r["N"]
        ratio = exp(r["ln_relvar_m4"] - r["ln_relvar_m2"])
        print(f"  {N:>5} {r['p']:>3} {r['ln_relvar_m4']:>14.2f} "
              f"{r['ln_relvar_m4'] / sqrt(N):>8.3f} {ratio:>8.4f}")
    print("\n  relVar(m4) = 4 relVar(m2) (1 + o(1)): the 4th moment self-averages")
    print("  with the SAME -(1/2) sqrt N ln N super-exponential law as m2 --")
    print("  Eq. (3.27) extends to m4 at a fixed factor of 4.")
