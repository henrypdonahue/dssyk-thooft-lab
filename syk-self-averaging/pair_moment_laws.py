#!/usr/bin/env python3
r"""
The exact 1/N laws of the pair-channel moments, by interpolation.
=================================================================

road_map.md lists the "analytic 1/N chord correction" as the open theory
half of the tower extraction.  That is the law behind the measured
fixed-p tail slopes 2.96, 3.90, 4.82, and behind the double-scaling
c/(pN).

This module closes the fixed-p half, without a chord derivation, by exact
computation.

THE OBSERVATION.  Let

    P(N) = 2 mu_k(N) C(N,p)^m        (m = (k+2n)/2 subsets per word)

Then P is a POLYNOMIAL in N of degree at most m p.

Proof.  Each Wick/signature term of the exact assembly
(exact_pair_moments.py) is an integer signed count of the form
sum_patterns prod_cells C(n_generic, s), with at most m p elements
placed.  So it is a polynomial in n_generic = N - n_dist of degree at
most m p.  And the denominator collected over a word is exactly
C(N,p)^m.

So m p + 1 exact samples determine mu_k(N) at fixed p IN CLOSED FORM.
Every 1/N coefficient is then an exact rational number, not a fit.  In
particular the leading law

    mu_k(N) = A / N^d + B / N^(d+1) + ...,   d = m p - deg P

FINDINGS.  All exact.  The artifact is pair_moment_laws.json, pinned by
test_pair_moment_laws.py.

* k = 2 deficits are d = n + 2 at every computed (n, p): n = 1..3 across
  p = 2..12, and n = 4 at p = 2 and 4.  So the float campaign's fitted
  tail slopes 2.96 / 3.90 / 4.82 are the shadows of the exact powers
  N^-3 / N^-4 / N^-5, and the n = 4 law is N^-6.

* Closed forms in p, certified out of sample at p = 14:

      A(n=1, p) = -16 p^2 (p-1)
      A(n=2, p) = -64 p^3 (p-1)(4p-5)

  A(n=3, p) has degree above 5 in p.  The minimum-degree interpolant
  fails the p = 14 check, and that is reported rather than hidden.  Its
  exact values through p = 12 are in the artifact.

* k = 4 deficits are n + 3, for n = 1 and 2.  So mu_4 dies one power
  faster than mu_2 in the pair channel.

  And the exact leading coefficient of mu_4(n=1) EQUALS that of
  mu_2(n=2), at every shared p.  That is consistent with ad-shifting
  between channels, since ad(X_1) = X_2 plus a cross term that enters
  beyond leading order.  But it is a measured identity of the exact
  coefficients, not a derived one.

* Double scaling at lambda_chord = 2, so N = p^2, exact.
  c_2(p) = |mu_2(n=2)| p N runs 28.918 -> 30.529 over p = 4..16, and it is
  still rising.  That continues the committed float scan's 28.9 -> 29.8.
  It is not settled by p = 16.  c_3(p) runs 203.5 -> 206.6 over
  p = 4..12, and it is not monotone at small p.

COST.  The sample grid needs m p + 2 exact points.  The expensive regime
is large m, where the matchings go as (K-1)!!, crossed with large p,
where the DP tables go as (p+1)^m.  The default campaign is the
affordable grid.  `--deep` adds the heavier laws: n = 3 at p = 10 and 12,
and n = 4 at p = 4.

Output: pair_moment_laws.json.
"""

from __future__ import annotations

import argparse
import json
import time
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

from exact_pair_moments import (_DP_CACHE, _EXP_CACHE, channel_pair_moments)

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# exact polynomial helpers (coefficient lists, index = power, Fractions)
# ---------------------------------------------------------------------------
def _poly_mul(a: list, b: list) -> list:
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return out


def _poly_pow(a: list, k: int) -> list:
    out = [Fraction(1)]
    for _ in range(k):
        out = _poly_mul(out, a)
    return out


def _poly_eval(a: list, x: Fraction) -> Fraction:
    acc = Fraction(0)
    for c in reversed(a):
        acc = acc * x + c
    return acc


def _binom_poly(p: int) -> list:
    """C(N, p) as a polynomial in N: prod_{i<p} (N - i) / p!."""
    out = [Fraction(1)]
    for i in range(p):
        out = _poly_mul(out, [Fraction(-i), Fraction(1)])
    return [c / factorial(p) for c in out]


def _newton_to_monomial(vals: list, n0: int) -> list:
    """Exact monomial coefficients of the unique polynomial through
    (n0 + j, vals[j]): Newton forward differences, then expansion of the
    falling-factorial basis C(N - n0, j)."""
    diffs, row = [], [Fraction(v) for v in vals]
    while row:
        diffs.append(row[0])
        row = [row[i + 1] - row[i] for i in range(len(row) - 1)]
        if row and all(x == 0 for x in row):
            diffs.extend([Fraction(0)] * 0)
            break
    poly = [Fraction(0)]
    basis = [Fraction(1)]                     # C(N - n0, 0)
    for j, dj in enumerate(diffs):
        if dj:
            term = [c * dj for c in basis]
            poly = [a + b for a, b in
                    zip(poly + [Fraction(0)] * (len(term) - len(poly)),
                        term + [Fraction(0)] * (len(poly) - len(term)))]
        # basis_{j+1} = basis_j * (N - n0 - j) / (j + 1)
        basis = [c / (j + 1) for c in
                 _poly_mul(basis, [Fraction(-(n0 + j)), Fraction(1)])]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def _series_at_infinity(P: list, D: list, orders: int = 2):
    """(d, coeffs): P/D = sum_j coeffs[j] / N^(d+j), exact, with
    d = deg D - deg P.  Power-series division on reversed coefficients."""
    dP, dD = len(P) - 1, len(D) - 1
    a = list(reversed(P))                     # a[0] = leading
    b = list(reversed(D))
    out = []
    for _ in range(orders):
        c = a[0] / b[0]
        out.append(c)
        a = [ai - c * bi for ai, bi in
             zip(a + [Fraction(0)] * len(b), b + [Fraction(0)] * len(a))]
        a.pop(0)
    return dD - dP, out


# ---------------------------------------------------------------------------
# the law extraction
# ---------------------------------------------------------------------------
def mu_exact(N: int, p: int, n: int, k: int) -> Fraction:
    """One exact pair-channel moment; caches cleared so long scans keep
    flat memory (the caches key on N and are useless across N)."""
    out = channel_pair_moments(N, p, n, ks=(k,), exact=True)[f"mu{k}"]
    _DP_CACHE.clear()
    _EXP_CACHE.clear()
    return out


def law(p: int, n: int, k: int = 2, orders: int = 2,
        progress: bool = False) -> dict:
    """The exact closed form of mu_k(N) at fixed p: interpolate
    P(N) = 2 mu_k C(N,p)^m through m p + 2 integer points (one more than
    the degree bound, so the top Newton difference vanishing certifies
    the bound), and expand at N = infinity.

    Returns dict with the integer polynomial P (as strings for JSON
    safety), the deficit d, and the exact leading coefficients:
    mu_k = A/N^d + B/N^(d+1) + ...  A check evaluation at an
    out-of-sample N is included."""
    m = (k + 2 * n) // 2
    deg_bound = m * p
    n0 = 2 * p + 2
    Ns = list(range(n0, n0 + deg_bound + 2))
    t0 = time.time()
    vals = []
    for N in Ns:
        mu = mu_exact(N, p, n, k)
        P_val = 2 * mu * Fraction(comb(N, p)) ** m
        if P_val.denominator != 1:
            raise RuntimeError(f"P({N}) not an integer at (p,n,k)="
                               f"({p},{n},{k}): {P_val}")
        vals.append(P_val.numerator)
        if progress:
            print(f"    sampled N={N} ({time.time() - t0:.1f} s)",
                  flush=True)
    P = _newton_to_monomial(vals, n0)
    if len(P) - 1 > deg_bound:
        raise RuntimeError(f"degree {len(P) - 1} exceeds the bound "
                           f"{deg_bound} at (p,n,k)=({p},{n},{k})")
    D = [2 * c for c in _poly_pow(_binom_poly(p), m)]
    d, coeffs = _series_at_infinity(P, D, orders)
    # out-of-sample certification: the closed form must reproduce a direct
    # exact evaluation at an N the interpolation never saw
    N_chk = n0 + deg_bound + 5
    direct = mu_exact(N_chk, p, n, k)
    interp = _poly_eval(P, Fraction(N_chk)) / _poly_eval(D, Fraction(N_chk))
    if direct != interp:
        raise RuntimeError(f"closed form fails out-of-sample at N={N_chk}, "
                           f"(p,n,k)=({p},{n},{k})")
    return dict(p=p, n=n, k=k, m=m, deg_P=len(P) - 1, deficit=d,
                A=str(coeffs[0]), B=str(coeffs[1]) if orders > 1 else None,
                A_float=float(coeffs[0]),
                P_coeffs=[str(c) for c in P],
                out_of_sample_N=N_chk, runtime_s=round(time.time() - t0, 1))


def double_scaling_exact(n: int, ps) -> list:
    """Exact mu_2 on the line lambda_chord = 2 (N = p^2):
    c(p) = |mu_2| p N, the constant the committed float scan fits."""
    rows = []
    for p in ps:
        N = p * p
        mu = mu_exact(N, p, n, 2)
        rows.append(dict(p=p, N=N, mu2=float(mu),
                         c=float(abs(mu) * p * N)))
        print(f"  lambda_chord=2, n={n}: p={p:2d} N={N:3d}  "
              f"mu2={float(mu):+.6e}  c={rows[-1]['c']:.3f}", flush=True)
    return rows


def _lagrange_in_p(pts: list) -> list:
    """Exact monomial coefficients of the minimum-degree polynomial
    through the (p, value) points."""
    poly = [Fraction(0)]
    for i, (xi, yi) in enumerate(pts):
        basis = [Fraction(1)]
        denom = Fraction(1)
        for j, (xj, _) in enumerate(pts):
            if j == i:
                continue
            basis = _poly_mul(basis, [Fraction(-xj), Fraction(1)])
            denom *= xi - xj
        term = [c * yi / denom for c in basis]
        poly = [a + b for a, b in
                zip(poly + [Fraction(0)] * (len(term) - len(poly)),
                    term + [Fraction(0)] * (len(poly) - len(term)))]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def a_closed_form(n: int, laws_rows: list, p_check: int = None,
                  k: int = 2) -> dict:
    """Candidate closed form of the leading coefficient A(p) at fixed n:
    exact minimum-degree polynomial interpolation over the sampled p
    grid, then an OUT-OF-SAMPLE certification at p_check -- a fresh law()
    extraction at a p the interpolation never saw.  A match makes the
    closed form exact at every sampled point plus one blind point; a
    mismatch is reported as not-identified (degree above the sample
    depth), never silently dropped."""
    pts = sorted(((r["p"], Fraction(r["A"])) for r in laws_rows
                  if r["n"] == n and r["k"] == k), key=lambda t: t[0])
    if len(pts) < 4:
        return dict(n=n, identified=False, reason="fewer than 4 p samples")
    poly = _lagrange_in_p(list(pts))
    out = dict(n=n, k=k, p_samples=[int(x) for x, _ in pts],
               A_poly_coeffs=[str(c) for c in poly],
               degree=len(poly) - 1)
    if p_check is not None:
        predicted = _poly_eval(poly, Fraction(p_check))
        fresh = Fraction(law(p_check, n, k)["A"])
        out["p_check"] = p_check
        out["identified"] = bool(predicted == fresh)
        if not out["identified"]:
            out["reason"] = (f"degree > {len(poly) - 1}: predicted "
                             f"{predicted} vs computed {fresh} at "
                             f"p = {p_check}")
    else:
        out["identified"] = False
        out["reason"] = "no out-of-sample check run"
    return out


def run(deep: bool = False, out_path=None) -> dict:
    t0 = time.time()
    laws = []
    plan = [(n, p, 2) for n in (1, 2) for p in (2, 4, 6, 8, 10, 12)]
    plan += [(3, p, 2) for p in (2, 4, 6, 8)]
    plan += [(1, p, 4) for p in (2, 4, 6)] + [(2, p, 4) for p in (2, 4)]
    plan += [(4, 2, 2)]
    if deep:
        plan += [(3, 10, 2), (3, 12, 2), (4, 4, 2)]
    for (n, p, k) in plan:
        print(f"law n={n} p={p} k={k} ...", flush=True)
        row = law(p, n, k, progress=(n >= 4 or (n == 3 and p >= 8)))
        laws.append(row)
        print(f"  deficit={row['deficit']}  A={row['A']}  "
              f"({row['runtime_s']} s)", flush=True)
    ds = dict(n2=double_scaling_exact(2, (4, 6, 8, 10, 12, 14, 16)),
              n3=double_scaling_exact(3, (4, 6, 8, 10, 12)))
    # closed forms of A(p): out-of-sample certification at p = 14
    # (n = 1, 2 always; n = 3 only in deep, where the p grid is deep
    # enough for a candidate and the p = 14 law is affordable)
    A_of_p = [a_closed_form(1, laws, p_check=14),
              a_closed_form(2, laws, p_check=14)]
    if deep:
        print("closed form n=3: out-of-sample law at p=14 ...", flush=True)
        A_of_p.append(a_closed_form(3, laws, p_check=14))
    payload = dict(
        provenance="pair_moment_laws.py -- exact closed-form 1/N laws of "
                   "the pair-channel moments at fixed p by rational "
                   "interpolation (P(N) = 2 mu_k C(N,p)^m is a polynomial, "
                   "degree bound m p, certified out-of-sample); the "
                   "fixed-p half of road_map's 'analytic 1/N chord "
                   "correction', now exact",
        runtime_s=round(time.time() - t0, 1),
        laws=laws, double_scaling=ds, A_of_p=A_of_p)
    out = Path(out_path) if out_path else HERE / "pair_moment_laws.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--deep", action="store_true",
                    help="add the heavier laws: n = 3 at p = 10, 12 and "
                         "n = 4 at p = 4")
    args = ap.parse_args()
    run(deep=args.deep)
