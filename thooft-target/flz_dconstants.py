#!/usr/bin/env python3
r"""
The FLZ spectral-determinant constants d+- at alpha = 0, from this
repo's own spectrum.
====================================================================

This module consumes the last unconsumed external anchor in
external_anchors.json: flz_d_constants_alpha0, the twelve-digit
constants d+ and d- of Fateev-Lukyanov-Zamolodchikov (FLZ,
arXiv:0905.2280), Eq. (5.1):

    d+ = 0.963178456398 ,    d- = 0.433582639833      (alpha1 = alpha2 = 0).

WHAT d+- ACTUALLY ARE (pinned from the FLZ full text, 2026-08-27)
-----------------------------------------------------------------
NOT wavefunction integrals.  The anchor's original status note guessed
"endpoint constants"; the paper says otherwise.  FLZ define spectral
determinants over the two CP towers (their Eq. 2.18),

    D+(lambda) = prod_{m>=0} (1 - lambda/lambda_{2m})   e^{lambda/lambda_{2m}} ,
    D-(lambda) = e^{2 lambda}
                 prod_{m>=0} (1 - lambda/lambda_{2m+1}) e^{lambda/lambda_{2m+1}} ,

and d+- are the overall constants in their lambda -> -infinity
asymptotics (their Eq. 4.16),

    D+-(lambda) ~ d+- (8 pi e^{-2+gamma_E})^lambda (-lambda)^{lambda - 1/8 +- 1/4}
                  exp[ sum_{k>=1} F^{(k)}_{+-}(L) lambda^{-k} ] .

FLZ state they have NO analytic expression for d+ or d- individually;
only the ratio is exact (their Eq. 4.20):

    d-/d+ = sqrt(2)/pi .

What makes the constants computable from a spectrum is their Eq. (4.21),
a pair of convergent products over the eigenvalues alone:

    d+ = Gamma(3/8)/sqrt(2 pi) * prod_{m>=0} (m + 3/8) / lambda_{2m} ,
    d- = Gamma(7/8)/sqrt(2 pi) * prod_{m>=0} (m + 7/8) / lambda_{2m+1} ,

(lambda_n = two_lambda_n / 2 in this repo's reporting convention
2 lambda_n = M_n^2/(pi g^2); even FLZ index n = 2m is the symmetric
sector, odd n = 2m+1 the antisymmetric one).  FLZ obtained their
Eq. (5.1) numbers by evaluating exactly these products on their
numerical spectrum.  We do the same on ours.

So consuming this anchor certifies the SPECTRUM through a new
functional -- an infinite product with Gamma-function grading, sensitive
to every computed level that enters the head (n <= n_cut = 80) or the
residual fit window (n <= 90), and in particular to the m + 3/8 vs
m + 7/8 offset between the towers -- fully independent of the power-sum
rules G^(s) already tested.  (Levels above n = 90 enter only through
the asymptotic tail; the G^(s) sum rules cover them.)  It does NOT
certify the eigenfunctions; the eigenfunction-level checks remain the
internal ones in eigenfunctions.py (Parseval, weak-form Rayleigh,
cross-solver overlap).

THE TAIL (transcribed FLZ Eq. 1.6)
----------------------------------
With nu = n + 3/4, FLZ's large-n expansion at alpha = 0 is

    2 lambda_n = nu - 2/(3 pi^6 nu^3)
                 + [2 (-1)^{n+1} / (pi^4 nu^2)]
                   * { 1 - 4 log(pi e^{gamma_E - 1/2} nu) / (pi^2 nu) }
                 + O(log^2 n / n^4) .

Each product factor is (m+a)/lambda = 1/(1 + delta_n/(2(m+a))) with
a = 3/8 or 7/8 and delta_n = 2 lambda_n - nu, so truncating the product
at n_cut discards ~ sum delta_nu/nu ~ 1/(pi^4 nu_cut^2): 7e-7 at
n_cut = 80.  The tail is therefore summed explicitly with delta from
Eq. (1.6) (a one-signed, absolutely convergent sum; the alternating
(-1)^{n+1} has fixed sign within one parity tower).

MEASURED-RESIDUAL ACCELERATION (optional, on by default)
--------------------------------------------------------
Eq. (1.6) is exact only through O(log n / n^3); its own error,
O(log^2 n / n^4), leaves ~2e-10 in the tail at n_cut = 80.  The last
two digits come from fitting the measured residual

    r(n) = (2 lambda_n^{computed} - nu) - delta^{(1.6)}(n)
         ~ (u + v L + w L^2) / nu^4 ,     L = log nu ,

per sector on 30 <= n <= 90 (default), and adding the fitted r into the
tail.  This is the same drift-fit pattern as jacobi_solver.sum_rule_G.
The fit sees only this solver's own spectrum -- the anchor values enter
nowhere in the computation, only in the asserting tests.

WHAT IS TRANSCRIBED vs DERIVED
------------------------------
Transcribed from FLZ (via the arXiv full text): Eqs. (2.18), (4.16),
(4.20), (4.21), (5.1), (1.6), quoted above.  The anchor numbers
themselves live in external_anchors.json, the one provenance-checked
copy; this module deliberately does not duplicate them.  Derived here:
everything computed, from the repo's eigenvalues only.

MEASURED AGREEMENT (deterministic; asserted with margin in
test_thooft.py, section "FLZ spectral-determinant constants")
-------------------------------------------------------------
Chebyshev solver, n_cut = 80, tail + residual fit on 30..90, deviation
from the FLZ 12-digit values (d+ / d-):

    K = 100:  -1.5e-09 / -6.6e-10
    K = 150:  -1.3e-10 / -6.0e-11
    K = 200:  -2.4e-11 / -1.1e-11
    K = 300:  -1.5e-12 / -7.6e-13
    K = 400:  +3.1e-13 / +2.4e-14
    K = 600:  +6.8e-13 / +1.9e-13
    K = 800:  +6.7e-13 / +2.0e-13

The ladder CONVERGES.  thooft_spectrum.solve_sector solves the
reciprocal problem at alpha = 0, so the low levels no longer sit on a
roundoff floor that grows with K, and the K = 300 row is basis
truncation of the levels near n = 90, not a floor.  (The direct solve
did have a floor there that grew with K; ERRATA E16.)

What limits the default configuration is the tail model: with the head
cut at n = 80 the exact ratio d-/d+ = sqrt(2)/pi is met to 1.8e-13.
Deepening the head and the residual-fit window with the basis removes
that too.  At K = 800:

    n_cut = 120, fit 40..150:   ratio to 3e-14
    n_cut = 160, fit 60..200:   ratio to 2e-15
    n_cut = 200, fit 80..260:   ratio to 6e-16

and the two constants settle at

    d+ = 0.9631784563984(1),      d- = 0.43358263983325(1),

both inside the +-5e-13 rounding of FLZ's printed 12 digits (which
truncate these), with one to two digits added.  The anchor is consumed
at the full precision FLZ published; the scatter over the three
configurations above (1.2e-13 in d+, 2e-14 in d-) is the quoted
uncertainty.

Ladder of the tail treatment at K = 300, n_cut = 80, d+ deviation:

    bare truncation      -7.3e-07
    + Eq. (1.6) tail     -1.5e-10
    + residual fit       -1.5e-12

Jacobi cross-solver at alpha = 0 (same trial space, independent
Gagliardo-quadrature assembly): K = 60 lands at -1.5e-08, K = 100 at
-1.5e-09, tracking its documented quadrature floor, with the ratio at
2.1e-10 and 9.4e-12.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln

from thooft_spectrum import solve_sector

PI = np.pi
GAMMA_E = float(np.euler_gamma)

#: numerator offsets in FLZ Eq. (4.21): m + 3/8 (even tower), m + 7/8 (odd)
_OFFSET = {0: 0.375, 1: 0.875}

#: the Eq. (1.6) tail is summed term-by-term out to this many product
#: factors; the remainder beyond it is ~1/(8 pi^4 m_max^2) ~ 3e-14
_M_TAIL = 2_000_000


def flz_delta_asymptotic(n: np.ndarray) -> np.ndarray:
    """Transcribed FLZ Eq. (1.6): delta_n = 2 lambda_n - (n + 3/4), valid
    through O(log n / n^3), error O(log^2 n / n^4)."""
    n = np.asarray(n, dtype=float)
    nu = n + 0.75
    t3 = -2.0 / (3.0 * PI ** 6 * nu ** 3)
    sign = np.where(np.asarray(n) % 2 == 0, -1.0, 1.0)   # (-1)^(n+1)
    t2 = (2.0 * sign / (PI ** 4 * nu ** 2)
          * (1.0 - 4.0 * np.log(PI * np.exp(GAMMA_E - 0.5) * nu)
             / (PI ** 2 * nu)))
    return t3 + t2


def _tail_log_sum(parity: int, m_start: int,
                  resid_coef: np.ndarray | None = None) -> float:
    """sum_{m >= m_start} log[(m+a)/lambda_n], n = 2m + parity, with
    lambda from Eq. (1.6) (+ the fitted residual model if given).

    Each factor is -log(1 + delta_n/(2(m+a))) since
    lambda_n = (nu + delta_n)/2 and 2(m+a) = nu.  Chunked so the peak
    footprint stays small."""
    total = 0.0
    chunk = 250_000
    for lo in range(m_start, _M_TAIL, chunk):
        m = np.arange(lo, min(lo + chunk, _M_TAIL), dtype=float)
        n = 2.0 * m + parity
        nu = n + 0.75
        d = flz_delta_asymptotic(n)
        if resid_coef is not None:
            L = np.log(nu)
            d = d + (resid_coef[0] + resid_coef[1] * L
                     + resid_coef[2] * L * L) / nu ** 4
        total += float(np.sum(-np.log1p(d / nu)))
    return total


def d_constant_sector(two_lambda: np.ndarray, parity: int, n_cut: int = 80,
                      tail: bool = True, resid_fit: bool = True,
                      fit_window: tuple[int, int] = (30, 90)) -> float:
    """One FLZ Eq. (4.21) product from one sector's computed spectrum.

    two_lambda : ascending 2*lambda_n of one CP sector (m-th entry is FLZ
                 level n = 2m + parity), as returned by
                 thooft_spectrum.solve_sector or
                 jacobi_solver.solve_sector_jacobi at alpha = 0.
    parity     : 0 (symmetric tower -> d+), 1 (antisymmetric -> d-).
    n_cut      : head of the product taken from the computed levels
                 (FLZ index n <= n_cut).
    tail       : append the Eq. (1.6) asymptotic tail (else bare
                 truncation, ~7e-7 low at n_cut = 80).
    resid_fit  : accelerate the tail with the measured residual of
                 Eq. (1.6), fitted as (u + v L + w L^2)/nu^4 on
                 fit_window (needs computed levels covering the window).

    Returns d+ or d- as a float.  Deterministic.
    """
    a = _OFFSET[parity]
    tl = np.asarray(two_lambda, dtype=float)
    lam = tl / 2.0
    ms = np.arange(len(lam), dtype=float)
    ns = 2.0 * ms + parity
    head = ns <= n_cut
    if not np.any(~head):
        raise ValueError(f"n_cut = {n_cut} consumes all {len(lam)} computed "
                         "levels; no reliable head/tail split -- raise the "
                         "basis size or lower n_cut")
    log_prod = float(np.sum(np.log((ms[head] + a) / lam[head])))
    m_start = int(np.sum(head))

    coef = None
    if resid_fit:
        lo, hi = fit_window
        sel = (ns >= lo) & (ns <= hi)
        if int(np.sum(sel)) < 6:
            raise ValueError(f"fit_window {fit_window} covers only "
                             f"{int(np.sum(sel))} levels of parity {parity}; "
                             "need >= 6")
        nu = ns[sel] + 0.75
        r = (tl[sel] - nu) - flz_delta_asymptotic(ns[sel])
        L = np.log(nu)
        X = np.vstack([np.ones_like(nu), L, L * L]).T / nu[:, None] ** 4
        coef = np.linalg.lstsq(X, r, rcond=None)[0]

    if tail:
        log_prod += _tail_log_sum(parity, m_start, resid_coef=coef)
    log_pref = gammaln(a) - 0.5 * np.log(2.0 * PI)
    return float(np.exp(log_pref + log_prod))


def d_constants(num_basis: int = 300, n_cut: int = 80, tail: bool = True,
                resid_fit: bool = True,
                fit_window: tuple[int, int] = (30, 90),
                solver: str = "chebyshev") -> tuple[float, float]:
    """(d+, d-) from a fresh alpha = 0 solve of both CP sectors.

    solver = "chebyshev" (thooft_spectrum, exact matrix elements) or
    "jacobi" (jacobi_solver at alpha = 0: same trial space, independent
    quadrature assembly -- the cross-check path)."""
    if solver == "chebyshev":
        def solve(parity):
            ev, _, _ = solve_sector(num_basis, parity=parity,
                                    eigvals_only=True)
            return ev
    elif solver == "jacobi":
        from jacobi_solver import solve_sector_jacobi

        def solve(parity):
            ev, _, _ = solve_sector_jacobi(num_basis, parity, 0.0,
                                           eigvals_only=True)
            return ev
    else:
        raise ValueError(f"unknown solver {solver!r}")
    dp = d_constant_sector(solve(0), 0, n_cut, tail, resid_fit, fit_window)
    dm = d_constant_sector(solve(1), 1, n_cut, tail, resid_fit, fit_window)
    return dp, dm


if __name__ == "__main__":
    import json
    from pathlib import Path

    anchors = json.loads((Path(__file__).resolve().parent
                          / "external_anchors.json").read_text())
    ref = anchors["flz_d_constants_alpha0"]
    dp_flz, dm_flz = float(ref["d+"]), float(ref["d-"])
    ratio_exact = np.sqrt(2.0) / PI

    print("FLZ spectral-determinant constants d+- at alpha = 0")
    print(f"  anchor (FLZ Eq. 5.1): d+ = {dp_flz:.12f}   d- = {dm_flz:.12f}")
    print(f"  exact ratio (Eq. 4.20): sqrt(2)/pi = {ratio_exact:.14f}")
    print(f"  anchor self-consistency: d-/d+ - sqrt(2)/pi = "
          f"{dm_flz / dp_flz - ratio_exact:+.1e}\n")

    print("Chebyshev solver, n_cut = 80, Eq. (1.6) tail + residual fit:")
    print("  (converges in K; beyond K ~ 400 the 80-level head and its tail"
          " model set the last digit)")
    ev300 = {}
    for K in (100, 150, 200, 300, 600):
        evs = {par: solve_sector(K, parity=par, eigvals_only=True)[0]
               for par in (0, 1)}
        if K == 300:
            ev300 = evs        # reused by the ladder below: one solve
        dp = d_constant_sector(evs[0], 0)
        dm = d_constant_sector(evs[1], 1)
        print(f"  K = {K:3d}:  d+ dev = {dp - dp_flz:+.2e}   "
              f"d- dev = {dm - dm_flz:+.2e}   "
              f"ratio - sqrt(2)/pi = {dm / dp - ratio_exact:+.2e}")

    print("\nTail-treatment ladder at K = 300 (d+ deviation):")
    for tail, fit, label in ((False, False, "bare truncation  "),
                             (True, False, "+ Eq. (1.6) tail "),
                             (True, True, "+ residual fit   ")):
        dp = d_constant_sector(ev300[0], 0, tail=tail, resid_fit=fit)
        print(f"  {label}: {dp - dp_flz:+.2e}")

    print("\nDeep head and fit window at K = 800 (the exact ratio is the gauge):")
    ev800 = {par: solve_sector(800, parity=par, eigvals_only=True)[0]
             for par in (0, 1)}
    for n_cut, win in ((120, (40, 150)), (160, (60, 200)), (200, (80, 260))):
        dp = d_constant_sector(ev800[0], 0, n_cut=n_cut, fit_window=win)
        dm = d_constant_sector(ev800[1], 1, n_cut=n_cut, fit_window=win)
        print(f"  n_cut = {n_cut}, fit {win[0]}..{win[1]}:  d+ = {dp:.15f}   "
              f"d- = {dm:.15f}   ratio - sqrt(2)/pi = {dm / dp - ratio_exact:+.1e}")

    print("\nJacobi cross-solver (alpha = 0, quadrature-independent):")
    dp, dm = d_constants(num_basis=60, n_cut=50, fit_window=(30, 60),
                         solver="jacobi")
    print(f"  K =  60:  d+ dev = {dp - dp_flz:+.2e}   "
          f"d- dev = {dm - dm_flz:+.2e}   "
          f"ratio - sqrt(2)/pi = {dm / dp - ratio_exact:+.2e}")
