#!/usr/bin/env python3
r"""
The singlet tower as lambda -> 0: a conformal OPE plus a ladder rung.
=====================================================================

singlet_tower.py writes every Majorana singlet channel
A_n = sum_i psi_i (ad_H)^n psi_i (n <= 8) at N = infinity as a conserved
polynomial in H plus matter chords O_m of dimension m = 2..n, each dressed
by a polynomial F_{n,m}(E, E') acting from the left (E) and right (E').
It measured, at small lambda, that the even chords fade and the odd ones
survive (2, 3, 5, 7), the pattern of the large-q SYK bilinear spectrum,
and road_map.md Part III item 3 asks whether the chords behave as the
conformal operators of that spectrum.  This module answers in the
semiclassical limit.

THE LIMIT
---------
Write sigma = (E + E')/2 and omega = E - E' (the frequency: ad_H acts on a
dressed chord as omega).  At small lambda a chord of dimension m is the
ladder sech^(2m)(a t), a = sqrt(lambda) (horizon_tower.py), so omega ~ a
and sigma ~ 1.  Scale omega = a w.  Then, through every level computed,

    F_{n,m}(E, E')  =  a^(n-m) P_{n,m}(w) (1 + O(a))

with P_{n,m} a polynomial of the parity of n - m and NO sigma dependence
(the sigma terms vanish as powers of a).  The coefficients of P_{n,m} are
rationals (chord_limits(): exact rational arithmetic at
lambda = 2^-9 .. 2^-14 and q within 1e-60 of e^-lambda, polynomial
extrapolation to lambda = 0, rational reconstruction; the residuals,
below 1e-13, are recorded).  The same holds for the norms: norm_m = nu_m lambda^(2 floor(m/2))
(1 + O(lambda)) with nu_2..nu_6 = 1, 2, 52/5, 88/7, 4120/7.

RESULT 1: WITHOUT THE REMAINDER, A_n IS AN SL(2) OPE
----------------------------------------------------
The chord rule of singlet_tower.py has one non-pairwise ingredient, the
"remainder": when the first and last slots of Omega_k contract with each
other, the enclosing chord sees the grade of what it encloses, and the
O(1) part is 2 lambda times an X-weight on each enclosed pair (the
r >= 1 entries of the word tables).  Switch it off (NoRemainder).  Then
every P_{n,m} is EXACTLY the prediction of a one-dimensional conformal
OPE

    sum_i psi_i(0) psi_i(t)  =  sum_m C_m (du)^m 1F1(m; 2m; du d_u) O_m(u),
    u = tanh(a t),  du = u(t) - u(0),  O_m^(t) = (du/dt)^m O_m^(u),

with each O_m a primary of weight h = m and Delta_psi -> 0
(frame_dressings()).  The top coefficient is the plain block coefficient,
P_{n,m} = 4 C(n, m) (m)_(n-m) / (2m)_(n-m) w^(n-m) + ..., and every lower
coefficient is the frame correction; e.g. A_4 > (36/5 w^2 + 16/5) O_2,
A_6 > (60/7 w^4 + 96/7 w^2 + 64/7) O_2.  The frame u = tanh(a t) is the
one in which the chord two-point function sech^(2m)(a t) is a conformal
two-point function; u = e^(2at) is Moebius-equivalent and gives the same
numbers, and u = tan(a t) flips the sign of every frame correction (it
fails).  The match is exact on all 28 coefficients below the diagonal
through n = 7 (the tests re-derive the 17 through n = 6;
conformal_tower.json has n = 7).

RESULT 2: THE REMAINDER IS THE FIRST LADDER RUNG
------------------------------------------------
With the remainder on, the chord answer differs from the OPE only in the
lower coefficients, and at leading order the difference is one ladder
rung acting on the full bilocal,

    B(0, t)  =  B_OPE(0, t) - 2 a^2 Int_{0<t1<t2<t} W(t2 - t1) B(t1, t2),

where B(t1, t2) = sum_j psi_j(t1) psi_j(t2) (so the total translation
e^(i t1 ad_H) acts on it).  This is the shape of the linearized large-q
Liouville equation d1 d2 dg = -2 J^2 e^g dg, whose rung is e^g = sech^2.
At W = 1 it reproduces every coefficient through n = 7 except the three
that see the tau^2 term of W.  Those need a sector-dependent tau^2 weight,
-5/6 for O_2 (A_6 and A_7 agree) and -9/10 for O_3, where Liouville's
e^g = sech^2 would give -1 (rung_tau2_fit()).  So the rung's leading form
is identified and its subleading form is not.  The remainder is what
makes the full channels non-conformal at lambda -> 0; its weight is O(1)
there (it moves the O_5 share of A_5 by about 10%).

RESULT 3: EXACT SEMICLASSICAL SHARES, AND WHY THE ODD CHORDS SURVIVE
-------------------------------------------------------------------
With P_{n,m}, the norms and the ladder moments <w^(2k)>_m =
(-d^2/dt^2)^k sech^(2m)(t) at 0, the lambda -> 0 shares of the dynamical
weight are exact rationals (limiting_shares()):

    A_5 -> 2047/2058 O_3 + 11/2058 O_5                       (0.5345% O_5)
    A_6 -> 20443/20608 O_3 + 165/20608 O_5                   (0.8007% O_5)
    A_7 -> 540648/548261 O_3 + 7578/548261 O_5 + 5/78323 O_7 (0.0064% O_7)

singlet_tower.py's values at lambda = 1/1280 approach these.  The weight
of dressed O_m in A_n scales as lambda^(n - m + 2 floor(m/2)): lambda^(n-1)
for odd m, lambda^n for even m.  So the even chords fade BECAUSE of the
norm law norm_m ~ lambda^(2 floor(m/2)); given the law, the survival
pattern is derived, not observed.

RESULT 4: THE T = infinity OPE COEFFICIENTS ARE THE IR ONES, TIMES (1-k_c)^2
---------------------------------------------------------------------------
In block normalization C_m^2 = 16 nu_m / m!^2 lambda^(2 floor(m/2) - m)
(ope_coefficients()).  The even ones are O(1): C_2^2 = 4, C_4^2 = 13/45,
C_6^2 = 103/5670.  The odd ones are O(1/lambda): lambda C_3^2 = 8/9,
lambda C_5^2 = 22/1575, lambda C_7^2 = 1/2205.  Against Gross-Rosenhaus
(1702.08016 (2.15)), the q -> infinity couplings c_h^2 of the IR
conformal tower, exactly for h = 3, 5, 7,

    C_h^2 = N (1 - k_c(h))^2 c_h^2,   k_c(h) = 2/(h(h-1)),   N = 2 p^2/lambda,

with k_c the q = infinity conformal ladder-kernel eigenvalue
(Maldacena-Stanford (3.76)).  The form was read off from the h = 5 and 7
ratios (27/20 and 10/7 are (s-1)/s over its h = 3 value, s = h(h-1)/2);
the absolute normalization at h = 3, exactly N, is an independent check.
It is not derived.  It predicts norm_9 -> (21482496/143) lambda^8, recorded in the
artifact for a future n = 9 computation.

SCOPE.  Majorana A_n, strict N = infinity, the leading order in lambda.
"Exact" for the limits means exact rational arithmetic plus polynomial
extrapolation and rational reconstruction, with residuals recorded; it
is not a symbolic proof.  Results 1 and 2 are verified identities on
the levels computed, not derivations.

Output: conformal_tower.json.  `python conformal_tower.py` prints the
n <= 5 comparison in seconds; `--run [KMAX]` writes the artifact
(KMAX = 4, n <= 6, about a minute; KMAX = 5, n = 7, builds the K = 16
word tables first, about 4 minutes on one core).
"""

from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import singlet_tower as st

HERE = Path(__file__).resolve().parent

LAMS = tuple(Fraction(1, 2 ** j) for j in range(9, 15))


# ---------------------------------------------------------------------------
# exact power series in t with coefficients polynomial in w
# (a series is a list over t-degree of {w-power: Fraction}; complex units
#  are tracked by working with F_n directly, see _to_dressings)
# ---------------------------------------------------------------------------
def _padd(p, r, c=1):
    out = dict(p)
    for k, v in r.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def _pmul(p, r):
    out: dict = {}
    for i, a in p.items():
        for j, b in r.items():
            out[i + j] = out.get(i + j, 0) + a * b
    return {k: v for k, v in out.items() if v != 0}


def _smul(s1, s2, nmax):
    out = [dict() for _ in range(nmax + 1)]
    for i, p in enumerate(s1):
        if not p:
            continue
        for j, r in enumerate(s2):
            if i + j > nmax or not r:
                continue
            out[i + j] = _padd(out[i + j], _pmul(p, r))
    return out


def _spow(s, k, nmax):
    out = [dict() for _ in range(nmax + 1)]
    out[0] = {0: Fraction(1)}
    for _ in range(k):
        out = _smul(out, s, nmax)
    return out


def _tanh_series(nmax):
    """tanh t = sum T_k t^k, exact rationals (from tanh' = 1 - tanh^2)."""
    T = [Fraction(0)] * (nmax + 1)
    if nmax >= 1:
        T[1] = Fraction(1)
    for n in range(1, nmax):
        sq = sum(T[i] * T[n - i] for i in range(n + 1))
        T[n + 1] = ((1 if n == 0 else 0) - sq) / (n + 1)
    return T


def _binom_series(alpha_c, alpha_w, sign, nmax):
    """(1 + sign*u)^(alpha_c + alpha_w * I*w) as a u-series whose
    coefficients are polynomials in (I w); the imaginary unit is carried
    as a formal variable y = I*w, so a coefficient {k: c} means sum c y^k."""
    out = [dict() for _ in range(nmax + 1)]
    out[0] = {0: Fraction(1)}
    cur = {0: Fraction(1)}
    for k in range(1, nmax + 1):
        # multiply by (alpha - (k-1)) / k * sign
        fac = {0: Fraction(alpha_c) - (k - 1), 1: Fraction(alpha_w)}
        cur = {e: v * sign / k for e, v in _pmul(cur, fac).items()}
        out[k] = cur
    return out


def frame_dressings(nmax: int = 8) -> dict:
    """The lambda -> 0 dressings predicted by the SL(2) OPE in the frame
    u = tanh(a t), a = 1: {(n, m): {w-power: Fraction}} with the
    normalization F_{m,m} = 4 (the coefficient of the new chord).

    O_m^(u)(u) e^(i w t(u)) = (1 - u)^(-m - y/2) (1 + u)^(-m + y/2), y = i w,
    because dt/du = 1/(1 - u^2) and t(u) = artanh(u).  The OPE is applied
    around u = 0 and du = tanh t is substituted."""
    T = _tanh_series(nmax)
    tser = [dict() for _ in range(nmax + 1)]
    for k, v in enumerate(T):
        if v:
            tser[k] = {0: v}
    out = {}
    for m in range(2, nmax + 1):
        A = _binom_series(-m, Fraction(-1, 2), -1, nmax)   # (1-u)^(-m-y/2)
        Bs = _binom_series(-m, Fraction(1, 2), 1, nmax)    # (1+u)^(-m+y/2)
        f = _smul(A, Bs, nmax)                             # u-series, coeffs in y
        # 1F1(m; 2m; du d_u) f at u = 0:  sum_k (m)_k/(k!(2m)_k) du^k f^(k)(0)
        #   = sum_k (m)_k/(2m)_k f_k du^(m+k)   (f^(k)(0)/k! = f_k)
        b = [dict() for _ in range(nmax + 1)]
        for k in range(0, nmax - m + 1):
            ck = Fraction(_rising(m, k), _rising(2 * m, k))
            term = _spow(tser, m + k, nmax)
            for n, p in enumerate(term):
                if p:
                    b[n] = _padd(b[n], _pmul(p, {e: ck * v for e, v in f[k].items()}))
        # b[n] is the t^n coefficient as a polynomial in y = i w.  The
        # dressing is F_n = n! (t^n coeff) / i^n, normalized to F_m = 4
        # (so F_n = 4 (n!/m!) (b_n / b_m) i^(m - n)).
        norm = b[m].get(0)
        for n in range(m, nmax + 1):
            F = {}
            for e, v in b[n].items():
                # y^e i^(m - n) = i^(e + m - n) w^e ; e + m - n is even
                if (e + m - n) % 2:
                    raise AssertionError("parity")
                sign = _sgn((e + m - n) // 2)
                F[e] = 4 * factorial(n) * v * sign / (norm * factorial(m))
            out[(n, m)] = {e: v for e, v in F.items() if v != 0}
    return out


def _sgn(k: int) -> int:
    """(-1)^k as an int for any integer k (a negative power would be a float)."""
    return -1 if k % 2 else 1


def _rising(x: int, k: int) -> int:
    out = 1
    for j in range(k):
        out *= x + j
    return out


def block_coefficient(n: int, m: int) -> Fraction:
    """The top coefficient of P_{n,m}: 4 C(n, m) (m)_(n-m) / (2m)_(n-m)."""
    k = n - m
    return Fraction(4 * comb(n, m) * _rising(m, k), _rising(2 * m, k))


# ---------------------------------------------------------------------------
# the ladder rung
# ---------------------------------------------------------------------------
def ladder_dressings(base: dict, nmax: int, rung=(1,), kappa=-2,
                     iterations: int = None) -> dict:
    """Solve B = B_base + K B sector by sector, with
    (K B)(0, t) = kappa Int_{0<t1<t2<t} W(t2 - t1) e^(i w t1) b(t2 - t1),
    a = 1, W(tau) = sum rung[k] tau^k.  base: {(n, m): {w-power: c}}.
    In dressing language, with c_n the t^n coefficients of W(tau) b(tau),
    F_N = F_N^base + kappa sum_n c_n n! i^(-n-2) w^(N-2-n).  kappa = -2 is
    the sign of the linearized Liouville rung, -2 J^2 e^g."""
    iterations = iterations or nmax
    out = {}
    for m in range(2, nmax + 1):
        F0 = {n: dict(base.get((n, m), {})) for n in range(nmax + 1)}
        F = {n: dict(p) for n, p in F0.items()}
        for _ in range(iterations):
            # beta_n = F_n i^n / n!  -> keep as (real poly, power of i)
            new = {}
            for N in range(nmax + 1):
                acc = dict(F0[N])
                for n in range(0, N - 1):
                    # c_n = sum_k rung[k] beta_(n-k); beta_j = F_j i^j / j!
                    for k, wk in enumerate(rung):
                        j = n - k
                        if j < 0 or not wk or not F.get(j):
                            continue
                        # term: kappa * wk * F_j i^j / j! * n! * i^(-n-2) * w^(N-2-n)
                        ph = j - n - 2          # power of i; even by parity
                        if ph % 2:
                            continue
                        coef = wk * (kappa * _sgn(ph // 2)
                                     * Fraction(factorial(n), factorial(j)))
                        acc = _padd(acc, {e + N - 2 - n: coef * v
                                          for e, v in F[j].items()})
                new[N] = acc
            F = new
        for n in range(m, nmax + 1):
            out[(n, m)] = {e: v for e, v in F[n].items() if v != 0}
    return out


# ---------------------------------------------------------------------------
# the chord engine's lambda -> 0 limit
# ---------------------------------------------------------------------------
class NoRemainder(st.TowerChord):
    """The chord rule with the grade remainder switched off: only the
    r = 0 entries of the word tables (no powers of 2 lambda)."""

    def omega_one(self, k, A):
        def fn():
            q, mu = self.q, self.mu
            return sum(c * q ** nx * (-mu if crossed else mu)
                       for (crossed, nx, r), c in st.omega_one_table(k, A)
                       if r == 0)
        return self._memo(("o1", k, A), fn)

    def omega_word(self, k1, k2, A, B):
        def fn():
            q, mu = self.q, self.mu
            wt = (mu + mu * mu, mu * mu - mu, mu * mu, -mu * mu)
            return sum(c * q ** nx * wt[kind]
                       for (kind, nx, r), c in st.omega_word_table(k1, k2, A, B)
                       if r == 0)
        return self._memo(("o2", k1, k2, A, B), fn)


def _q_of(lam: Fraction) -> Fraction:
    """A rational within 1e-60 of exp(-lam)."""
    import mpmath as mp
    with mp.workdps(90):
        v = mp.exp(-mp.mpf(lam.numerator) / lam.denominator)
        return Fraction(mp.nstr(v, 80)).limit_denominator(10 ** 60)


def _scaled_dressings(T, n: int) -> dict:
    """{(m, w-power): c / lambda^((n-m-power)/2)} for the sigma-free part of
    F_{n,m}, plus {(m, 'sigma', a, b): ...} for the sigma terms (which must
    vanish in the limit)."""
    out = {}
    lam = T.lam
    for m, f in T.a_kernels(n).items():
        # E = s + w/2, E' = s - w/2:  E^i E'^j -> sum over binomials
        poly: dict = {}
        for (i, j), c in f.items():
            for x in range(i + 1):
                for y in range(j + 1):
                    a_pow = x + y                       # sigma power
                    b_pow = (i - x) + (j - y)           # w power
                    v = c * comb(i, x) * comb(j, y) * Fraction(1, 2 ** b_pow) \
                        * (-1) ** (j - y)
                    poly[(a_pow, b_pow)] = poly.get((a_pow, b_pow), 0) + v
        k = n - m
        for (a_pow, b_pow), c in poly.items():
            if c == 0:
                continue
            if a_pow == 0 and (k - b_pow) % 2 == 0:
                out[(m, b_pow)] = c / lam ** ((k - b_pow) // 2)
            else:
                # scaled size of a sigma term: c * lam^((b-k)/2); report its
                # square to stay rational
                out[(m, "sigma", a_pow, b_pow)] = c * c * lam ** (b_pow - k)
    return out


def _extrapolate(xs, ys) -> Fraction:
    """Lagrange interpolation through (xs, ys), evaluated at 0, exactly."""
    tot = Fraction(0)
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        L = Fraction(1)
        for j, xj in enumerate(xs):
            if j != i:
                L *= Fraction(-xj) / (xi - xj)
        tot += yi * L
    return tot


def _reconstruct(v: Fraction, max_den: int = 10 ** 5) -> tuple:
    r = v.limit_denominator(max_den)
    return r, abs(float(v - r))


def chord_limits(kmax: int = 3, remainder: bool = True, lams=LAMS) -> dict:
    """The lambda -> 0 limits of the scaled dressings P_{n,m} (n <= kmax+2)
    and of the norms nu_m, from exact rational chord sums.

    Returns {"dressings": {(n, m): {w-power: Fraction}}, "norms": {m: Fraction},
    "max_residual": worst |extrapolated - reconstructed rational|,
    "sigma_max": the largest extrapolated sigma term (squared, scaled),
    "mu0": {(n, m): limit of sector weight / lambda^(n-1)}}."""
    cls = st.TowerChord if remainder else NoRemainder
    rows = []
    for lam in lams:
        T = cls(lam, _q_of(lam))
        T.decompose(kmax, check_len=0)
        d = {}
        for n in range(3, kmax + 3):
            for key, v in _scaled_dressings(T, n).items():
                d[(n,) + key] = v
            for m in T.a_kernels(n):
                d[("mu0", n, m)] = T.sector_moment(n, 0, m) / lam ** (n - 1)
        for m, v in T.norm.items():
            d[("norm", m)] = v / lam ** (2 * (m // 2))
        rows.append(d)
    keys = set(rows[0])
    for r in rows[1:]:
        keys |= set(r)
    xs = list(lams)
    dress: dict = {}
    norms: dict = {}
    mu0: dict = {}
    worst = 0.0
    sigma_max = 0.0
    for key in keys:
        ys = [r.get(key, Fraction(0)) for r in rows]
        L = _extrapolate(xs, ys)
        if key[0] == "norm":
            r, res = _reconstruct(L)
            norms[key[1]] = r
            worst = max(worst, res)
        elif key[0] == "mu0":
            mu0[(key[1], key[2])] = L
        elif len(key) > 3 and key[2] == "sigma":
            sigma_max = max(sigma_max, abs(float(L)))
        else:
            n, m, b = key
            r, res = _reconstruct(L)
            worst = max(worst, res)
            if r != 0:
                dress.setdefault((n, m), {})[b] = r
    # the top chord of each level: constant dressing 4
    for n in range(2, kmax + 3):
        dress[(n, n)] = {0: Fraction(4)}
    return dict(dressings=dress, norms=norms, mu0=mu0, max_residual=worst,
                sigma_max=sigma_max)


# ---------------------------------------------------------------------------
# shares and OPE coefficients
# ---------------------------------------------------------------------------
def ladder_moment(m: int, k: int) -> Fraction:
    """<w^k> of the Delta = m ladder sech^(2m)(t): (-1)^(k/2) d^k/dt^k at 0."""
    if k % 2:
        return Fraction(0)
    # sech^(2m)(t) = exp(2m log sech t); log sech t = -sum ... use series
    nmax = k
    # log cosh t = sum_{j>=1} c_j t^(2j); derive from tanh series: d/dt log cosh = tanh
    T = _tanh_series(nmax + 1)
    logc = [Fraction(0)] * (nmax + 1)
    for j in range(1, nmax + 1):
        logc[j] = T[j - 1] / j
    # g = -2m log cosh t ; exp(g) series
    g = [-2 * m * c for c in logc]
    e = [Fraction(0)] * (nmax + 1)
    e[0] = Fraction(1)
    for n in range(1, nmax + 1):           # e' = g' e
        e[n] = sum(j * g[j] * e[n - j] for j in range(1, n + 1)) / n
    return e[k] * factorial(k) * (-1) ** (k // 2)


def sector_weight(P: dict, m: int) -> Fraction:
    """<|P(w)|^2> over the Delta = m ladder (P of one parity, real)."""
    tot = Fraction(0)
    for i, a in P.items():
        for j, b in P.items():
            tot += a * b * ladder_moment(m, i + j)
    return tot if tot >= 0 else -tot


def limiting_shares(dressings: dict, norms: dict, nmax: int) -> dict:
    """{n: {"lead_power": p, "shares": {m: Fraction}, "fading": [m...]}}:
    the lambda -> 0 shares of A_n's dynamical weight.  Sector m weighs
    norm_m <|F|^2> ~ lambda^(2 floor(m/2) + n - m)."""
    out = {}
    for n in range(3, nmax + 1):
        ws = {}
        for m in range(2, n + 1):
            P = dressings.get((n, m))
            if not P or m not in norms:
                continue
            ws[m] = (2 * (m // 2) + n - m, norms[m] * sector_weight(P, m))
        lead = min(p for p, _ in ws.values())
        tot = sum(v for p, v in ws.values() if p == lead)
        out[n] = dict(lead_power=lead,
                      shares={m: v / tot for m, (p, v) in ws.items() if p == lead},
                      fading=sorted(m for m, (p, _) in ws.items() if p > lead))
    return out


def ope_coefficients(norms: dict) -> dict:
    """C_m^2 lambda^(m - 2 floor(m/2)) = 16 nu_m / m!^2 (GR block normalization,
    O_m unit-normalized in the u-frame)."""
    return {m: Fraction(16) * v / factorial(m) ** 2 for m, v in norms.items()}


def gr_ratio(n: int) -> Fraction:
    """Gross-Rosenhaus c_n^2 / c_1^2 at q -> infinity (1702.08016 (2.15);
    largeq_anchor.gr_tower), exactly: the Gamma ratio is rational."""
    def part(k):
        eps = Fraction(2 * k * k + k + 1, 2 * k * k + k - 1)
        s = k * (1 + 2 * k)
        # sqrt(pi) Gamma(2k+1) / (Gamma(2k+1/2) 2^(4k-2)), Gamma(2k+1/2) =
        # (4k)! sqrt(pi) / (4^(2k) (2k)!)
        g = Fraction(factorial(2 * k) * 4 ** (2 * k) * factorial(2 * k),
                     factorial(4 * k) * 2 ** (4 * k - 2))
        return eps * eps * Fraction(s, (s + 1) * (s - 1)) * g
    return part(n) / part(1)


def gr_c2_q2(n: int) -> Fraction:
    """q^2 c_n^2 at q -> infinity, the absolute GR coupling of the h = 2n + 1
    primary (1702.08016 (2.15)), exactly."""
    s = n * (2 * n + 1)
    eps_q = Fraction(s + 1, s - 1)
    g = Fraction(4 * factorial(2 * n) ** 2, factorial(4 * n))
    return eps_q * eps_q * Fraction(s, (s + 1) * (s - 1)) * g


def kc(h: int) -> Fraction:
    """The q = infinity conformal ladder-kernel eigenvalue 2/[h(h-1)]
    (Maldacena-Stanford (3.76))."""
    return Fraction(2, h * (h - 1))


def odd_ope_from_gr(h: int) -> Fraction:
    """lambda C_h^2 predicted by C_h^2 = N (1 - k_c(h))^2 c_h^2(GR), with
    N = 2 p^2 / lambda (q = p): 2 (1 - k_c)^2 q^2 c^2."""
    return 2 * (1 - kc(h)) ** 2 * gr_c2_q2((h - 1) // 2)


def rung_tau2_fit(frame: dict, chord_on: dict, nmax: int) -> list:
    """For each coefficient that the W = 1 ladder misses, the tau^2 weight
    w2 of W(tau) = 1 + w2 tau^2 that would fix it (the miss is linear in
    w2).  Returns [(n, m, w-power, w2)]."""
    l0 = ladder_dressings(frame, nmax, rung=(1, 0, 0))
    l1 = ladder_dressings(frame, nmax, rung=(1, 0, 1))
    out = []
    for key, P in sorted(chord_on.items()):
        if key[0] <= key[1]:
            continue
        for e in sorted(set(P) | set(l0[key])):
            v0, v1, d = l0[key].get(e, 0), l1[key].get(e, 0), P.get(e, 0)
            if v0 != d:
                out.append((key[0], key[1], e,
                            (d - v0) / (v1 - v0) if v1 != v0 else None))
    return out


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def _fmt(d: dict) -> dict:
    return {str(k): str(v) for k, v in sorted(d.items())}


def compare(kmax: int = 3, verbose: bool = True) -> dict:
    nmax = kmax + 2
    frame = frame_dressings(nmax)
    off = chord_limits(kmax, remainder=False)
    on = chord_limits(kmax, remainder=True)
    frame_match = {str(k): off["dressings"].get(k) == frame[k]
                   for k in sorted(frame) if k[0] > k[1] and k[0] <= nmax}
    ladder = ladder_dressings(frame, nmax, rung=(1,))
    ladder_match = {str(k): on["dressings"].get(k) == ladder[k]
                    for k in sorted(ladder) if k[0] > k[1] and k[0] <= nmax}
    top = {str((n, m)): frame[(n, m)].get(n - m) == block_coefficient(n, m)
           for (n, m) in frame if n >= m}
    shares = limiting_shares(on["dressings"], on["norms"], nmax)
    shares_engine = {}
    for n in range(3, nmax + 1):
        lead = {m: v for (nn, m), v in on["mu0"].items() if nn == n}
        tot = sum(lead.values())
        shares_engine[n] = {m: float(v / tot) for m, v in lead.items() if v / tot > 1e-12}
    out = dict(
        kmax=kmax, nmax=nmax,
        frame_dressings={str(k): _fmt(v) for k, v in sorted(frame.items())
                         if k[0] <= nmax},
        chord_off={str(k): _fmt(v) for k, v in sorted(off["dressings"].items())},
        chord_on={str(k): _fmt(v) for k, v in sorted(on["dressings"].items())},
        norms={str(m): str(v) for m, v in sorted(on["norms"].items())},
        norms_off_equal=on["norms"] == off["norms"],
        max_residual=max(on["max_residual"], off["max_residual"]),
        sigma_max=max(on["sigma_max"], off["sigma_max"]),
        top_is_block_coefficient=top,
        frame_matches_no_remainder=frame_match,
        ladder_w1_matches_full=ladder_match,
        shares={str(n): dict(lead_power=r["lead_power"],
                             shares={str(m): str(v) for m, v in r["shares"].items()},
                             shares_float={str(m): float(v) for m, v in r["shares"].items()},
                             fading=r["fading"])
                for n, r in shares.items()},
        shares_engine_mu0={str(n): {str(m): v for m, v in r.items()}
                           for n, r in shares_engine.items()},
        ope_c2_scaled={str(m): str(v) for m, v in ope_coefficients(on["norms"]).items()},
    )
    c2 = ope_coefficients(on["norms"])
    out["rung_tau2_fit"] = [dict(n=n, m=m, w_power=e,
                                 w2=None if w2 is None else str(w2))
                            for n, m, e, w2 in rung_tau2_fit(frame, on["dressings"], nmax)]
    out["odd_ope_vs_gr"] = {
        str(h): dict(chord=str(c2[h]), n_one_minus_kc_sq_gr=str(odd_ope_from_gr(h)),
                     equal=c2[h] == odd_ope_from_gr(h))
        for h in sorted(c2) if h % 2}
    # the relation's prediction for the first odd chord not computed
    h_next = nmax + 1 if (nmax + 1) % 2 else nmax + 2
    nu_next = odd_ope_from_gr(h_next) * factorial(h_next) ** 2 / 16
    out["prediction"] = dict(
        h=h_next, lam_C2=str(odd_ope_from_gr(h_next)), nu=str(nu_next),
        nu_float=float(nu_next),
        note=f"norm_{h_next} -> nu lambda^{2 * (h_next // 2)} if the relation holds")
    if 5 in c2:
        out["odd_ratio_5_3"] = dict(chord=str(c2[5] / c2[3]), gr=str(gr_ratio(2)),
                                    chord_over_gr=str(c2[5] / c2[3] / gr_ratio(2)))
    if 7 in c2:
        out["odd_ratio_7_3"] = dict(chord=str(c2[7] / c2[3]), gr=str(gr_ratio(3)),
                                    chord_over_gr=str(c2[7] / c2[3] / gr_ratio(3)))
    if verbose:
        print(f"n <= {nmax}: residual of rational reconstruction "
              f"{out['max_residual']:.1e}, largest sigma term {out['sigma_max']:.1e}")
        print("  top coefficient = SL(2) block coefficient:", all(top.values()))
        print("  no-remainder chords = tanh-frame OPE:", all(frame_match.values()))
        bad = [k for k, v in ladder_match.items() if not v]
        print("  full chords = OPE + one rung (W = 1):", "all" if not bad else f"all but {bad}")
        for n, r in out["shares"].items():
            print(f"  A_{n}: " + ", ".join(f"O_{m} {s}" for m, s in r["shares"].items())
                  + (f"  (fading {r['fading']})" if r["fading"] else ""))
        print("  norms nu_m:", out["norms"])
        if "odd_ratio_5_3" in out:
            print("  C_5^2/C_3^2 =", out["odd_ratio_5_3"])
        print("  C_h^2 = N (1 - k_c(h))^2 c_h^2(GR):",
              {h: r["equal"] for h, r in out["odd_ope_vs_gr"].items()})
        if out["rung_tau2_fit"]:
            print("  tau^2 rung weight each missed coefficient needs:",
                  [(r["n"], r["m"], r["w_power"], r["w2"]) for r in out["rung_tau2_fit"]])
    return out


def run(kmax: int = 4, out_path=None, workers: int = 1) -> dict:
    t0 = time.time()
    if kmax >= 5:
        # the K >= 16 word tables: one enumeration per size (singlet_tower)
        st.build_word_tables(st.needed_word_keys(kmax, 0), workers)
    payload = compare(kmax)
    payload["runtime_s"] = round(time.time() - t0, 1)
    payload["provenance"] = (
        "conformal_tower.py --run -- lambda -> 0 limits of the singlet-tower "
        "dressings from exact rational chord sums at lambda = 2^-9..2^-14, "
        "against the tanh-frame SL(2) OPE and the one-rung ladder")
    out = Path(out_path) if out_path else HERE / "conformal_tower.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
        fh.write("\n")
    print(f"wrote {out}")
    return payload


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", nargs="?", const=4, type=int, metavar="KMAX")
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()
    if args.run is not None:
        run(args.run, workers=args.workers)
    else:
        compare(3)
