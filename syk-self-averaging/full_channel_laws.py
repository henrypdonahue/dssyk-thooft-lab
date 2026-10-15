#!/usr/bin/env python3
r"""
The FULL A_n channel at fixed p, exactly: the corner where the match lives.
==========================================================================

road_map.md names fixed p, lambda -> 0 as the regime where the meson
spectrum match has to be made, and lists the order-of-limits corner as
open.  The ED bench (mn_majorana_bench.py) reaches N <= 20 there and
finds the raw A_2 width w2 = mu_2/mu_0 RISING along fixed p = 4 as
lambda_chord = 2p^2/N falls, 0.30 -> 1.20 over N = 10 -> 20.  That was
read as "the tower channel turns on toward the corner".

The exact engine (exact_pair_moments.py) evaluates BOTH halves of the
channel at any N, with no ED:

    <A_n ad^k A_n> = N mu_k^diag + N(N-1) mu_k^pair ,

so the fixed-p scan continues past the bench's window in seconds.  The
bench sits on the rising side of a PEAK.  At p = 4 the full width tops
out at N = 24 (lambda_chord = 4/3) and then dies; at p = 6 the peak is at
N ~ 48-72 (lambda_chord ~ 1-1.5).  Past the peak the exact law is

    w2_full(N; p) = 8 p (p - 2) / N + O(1/N^2)   (chord units, <H^2> = 1),

with the leading coefficient a certified closed form: the interpolated
laws at p = 2, 4, 6 give 0, 64, 192, and the form 8p(p-2) is certified
out of sample at p = 8, 10, 12 (384, 640, 960).  Both mu_0 and mu_2 of
the full channel are exact polynomials in N over C(N,p)^m, so this is
interpolation, not a fit, like pair_moment_laws.py.  The two leading
coefficients are closed forms on their own,

    mu_0^full = 64 p^2 (p-1)^2 / N^2 + ...,
    mu_2^full = 512 p^3 (p-1)^2 (p-2) / N^3 + ...,

certified on the same p grid.

READ IT IN PHYSICAL UNITS.  Chord units measure frequencies against the
extensive bandwidth: <H^2>/dim = C(N,p) p!/N^(p-1) ~ N at fixed p in the
paper's normalization Var(J) = p!/N^(p-1).  So the physical mean-square
frequency of the A_2 channel tends to the N-INDEPENDENT constant

    <omega^2>_{A_2} -> 8 p (p - 2)  (units J = 1),

i.e. omega_rms -> 8 at p = 4, 13.9 at p = 6, ~2 sqrt(2) p at large p --
a finite O(J) width at fixed p, and in the paper's dictionary
(alpha-bar = g-bar^2 N = p^2) a width^2 of 8 alpha-bar (1 - 2/p).  That
number is the spec: a DSSYK <-> 't Hooft dictionary at fixed p has to
reproduce a bilinear channel whose rms frequency is 2 sqrt(2 alpha-bar)
in the flat-space corner, and whose width in bandwidth units peaks at
lambda_chord ~ 4/3 before the corner.  Nothing here says what the
spectral SHAPE is; the boost calibration (duality/BOOST.md) says the
nodes, not the moments, carry the mass data.

The exact per-instance identity behind the channel is A_2 = 4 sum_i H_i^2
(conserved_points.py), H_i the part of H containing flavor i.

READ THE PEAK WITH ITS DENOMINATOR.  mu_0 here is identity-subtracted,
and A_2 overlaps the conserved H^2.  Removing that component exactly
(duality/singlet_tower.py: finite_n_a2_split, fixed_p_dynamic) leaves a
width with no maximum beyond the smallest sizes: it sits near 2 at the
small-N end, on the N = infinity law 2(1 - q^2), and falls monotonically
(from N = 12 at p = 4, from N = 16 at p = 6) to the same 8p(p-2)/N.  The peak at N = 24 is the conserved share of mu_0
draining away (85 % at N = 10, 19 % at N = 24, 0.01 % at N = 1024 for
p = 4).  The asymptotic law and the physical width 8p(p-2) J^2 are
untouched; the peak is a property of the normalization.

What this corrects.  "The tower channel turns on toward the corner"
(syk README, duality/CHORD.md, the bench docstring and its test) is the
rising side of the peak seen through a window that ends at N = 20.  The
bench rows are unchanged and still asserted; the statement about the
corner is replaced by the law above.

Output: full_channel_laws.json (scans + laws), pinned by
test_full_channel_laws.py.  Runtime: a few minutes.
"""

from __future__ import annotations

import argparse
import json
import time
from fractions import Fraction
from math import comb
from pathlib import Path

import exact_pair_moments as epm
from pair_moment_laws import (_binom_poly, _newton_to_monomial, _poly_eval,
                              _poly_pow, _series_at_infinity)

HERE = Path(__file__).resolve().parent


def full_moments(N: int, p: int, n: int, ks=(0, 2)) -> dict:
    """Exact (Fractions) full-channel moments N mu^diag + N(N-1) mu^pair,
    chord units, both orderings symmetrized, mu_0 tracelessed."""
    pr = epm.channel_pair_moments(N, p, n, marks=('i', 'i', 'j', 'j'),
                                  ks=ks, exact=True)
    dg = epm.channel_pair_moments(N, p, n, marks=('i', 'i', 'i', 'i'),
                                  ks=ks, exact=True)
    epm._DP_CACHE.clear()
    epm._EXP_CACHE.clear()
    return {k: N * dg[f"mu{k}"] + N * (N - 1) * pr[f"mu{k}"] for k in ks}


def sigma_h2(N: int, p: int) -> Fraction:
    """<H^2>/dim in the paper's normalization: C(N,p) p!/N^(p-1)."""
    from math import factorial
    return Fraction(comb(N, p) * factorial(p), N ** (p - 1))


def scan(p: int, n: int, Ns, ks=(0, 2), progress: bool = True) -> list:
    rows = []
    for N in Ns:
        t0 = time.time()
        f = full_moments(N, p, n, ks)
        w2 = f[2] / f[0]
        rows.append(dict(N=N, lam_chord=2 * p * p / N, w2=float(w2),
                         w2_physical=float(w2 * sigma_h2(N, p)),
                         N_w2=float(w2 * N), runtime_s=round(time.time() - t0, 2)))
        if progress:
            print(f"  p={p} n={n} N={N:5d} lam={rows[-1]['lam_chord']:.3f} "
                  f"w2={rows[-1]['w2']:.5f}  N*w2={rows[-1]['N_w2']:.3f}",
                  flush=True)
    return rows


def full_law(p: int, n: int, k: int, orders: int = 2) -> dict:
    """Exact 1/N expansion of the full-channel mu_k at fixed p by
    interpolating P(N) = 2 mu_k^full C(N,p)^m (m = (k+2n)/2) through
    m p + 2 integer points, certified at an out-of-sample N (the same
    construction as pair_moment_laws.law)."""
    m = (k + 2 * n) // 2
    deg_bound = m * p
    n0 = 2 * p + 2
    Ns = list(range(n0, n0 + deg_bound + 2))
    N_chk = n0 + deg_bound + 5
    t0 = time.time()
    vals = []
    for N in Ns:
        mu = full_moments(N, p, n, (k,))[k]
        P_val = 2 * mu * Fraction(comb(N, p)) ** m
        if P_val.denominator != 1:
            raise RuntimeError(f"P({N}) not an integer at (p,n,k)="
                               f"({p},{n},{k})")
        vals.append(P_val.numerator)
    P = _newton_to_monomial(vals, n0)
    if len(P) - 1 > deg_bound:
        raise RuntimeError(f"degree {len(P) - 1} exceeds bound {deg_bound}")
    D = [2 * c for c in _poly_pow(_binom_poly(p), m)]
    d, coeffs = _series_at_infinity(P, D, orders)
    direct = full_moments(N_chk, p, n, (k,))[k]
    interp = _poly_eval(P, Fraction(N_chk)) / _poly_eval(D, Fraction(N_chk))
    if direct != interp:
        raise RuntimeError(f"closed form fails out of sample at N={N_chk}")
    return dict(p=p, n=n, k=k, m=m, deg_P=len(P) - 1, deficit=d,
                A=str(coeffs[0]), B=str(coeffs[1]) if orders > 1 else None,
                P_coeffs=[str(c) for c in P], out_of_sample_N=N_chk,
                runtime_s=round(time.time() - t0, 1))


def width_law(p: int, n: int = 2) -> dict:
    """N w2_full -> L(p): the ratio of the leading coefficients of mu_2
    and mu_0 when their deficits differ by exactly one."""
    l2 = full_law(p, n, 2)
    l0 = full_law(p, n, 0)
    gap = l2["deficit"] - l0["deficit"]
    if gap < 1:
        raise RuntimeError(f"w2 does not vanish at p={p}: deficits "
                           f"{l2['deficit']}, {l0['deficit']}")
    # gap > 1 means the 1/N coefficient of w2 is exactly zero (p = 2, where
    # 8p(p-2) = 0 and w2 falls faster than 1/N)
    L = Fraction(l2["A"]) / Fraction(l0["A"]) if gap == 1 else Fraction(0)
    return dict(p=p, n=n, L=str(L), L_float=float(L),
                predicted_8p_pm2=8 * p * (p - 2), mu2_law=l2, mu0_law=l0)


def closed_form_L(p: int) -> int:
    return 8 * p * (p - 2)


def closed_form_A0(p: int) -> int:
    """Leading coefficient of the full-channel mu_0 (n = 2): 64 p^2 (p-1)^2 / N^2."""
    return 64 * p * p * (p - 1) ** 2


def closed_form_A2(p: int) -> int:
    """Leading coefficient of the full-channel mu_2 (n = 2):
    512 p^3 (p-1)^2 (p-2) / N^3 = 8p(p-2) x the mu_0 coefficient."""
    return 512 * p ** 3 * (p - 1) ** 2 * (p - 2)


def run(out_path=None) -> dict:
    t0 = time.time()
    scans = {}
    grids = {4: [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 40, 48, 64, 96,
                 128, 192, 256, 384, 512, 768, 1024],
             6: [14, 16, 18, 20, 22, 24, 28, 32, 36, 40, 48, 56, 64, 72, 84,
                 96, 128, 144, 192, 288, 384, 576]}
    for p, Ns in grids.items():
        for n in (2, 3):
            print(f"scan p={p} n={n}", flush=True)
            rows = scan(p, n, Ns)
            peak = max(rows, key=lambda r: r["w2"])
            scans[f"p{p}_n{n}"] = dict(p=p, n=n, rows=rows,
                                      peak_N=peak["N"],
                                      peak_lam_chord=peak["lam_chord"],
                                      peak_w2=peak["w2"])
    print("laws", flush=True)
    laws = []
    for p in (2, 4, 6, 8, 10, 12):
        w = width_law(p, 2)
        laws.append(w)
        print(f"  p={p}: N w2 -> {w['L']}  (8p(p-2) = {w['predicted_8p_pm2']})",
              flush=True)
    identified = all(Fraction(w["L"]) == closed_form_L(w["p"])
                     for w in laws)
    coeffs_identified = all(
        Fraction(w["mu0_law"]["A"]) == closed_form_A0(w["p"])
        and Fraction(w["mu2_law"]["A"]) == closed_form_A2(w["p"])
        for w in laws)
    payload = dict(
        provenance="full_channel_laws.py -- exact fixed-p scans of the "
                   "FULL A_n channel width (N mu^diag + N(N-1) mu^pair, "
                   "exact Fractions, exact_pair_moments engine) and the "
                   "exact leading law N w2 -> 8p(p-2), certified out of "
                   "sample at p = 8, 10, 12; peak positions on the scan grid",
        runtime_s=round(time.time() - t0, 1), scans=scans, width_laws=laws,
        L_closed_form="8 p (p - 2)", L_identified=identified,
        A0_closed_form="64 p^2 (p-1)^2", A2_closed_form="512 p^3 (p-1)^2 (p-2)",
        coefficients_identified=coeffs_identified,
        L_samples=[2, 4, 6], L_checks=[8, 10, 12])
    out = Path(out_path) if out_path else HERE / "full_channel_laws.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()
    if args.run:
        run()
    else:
        for N in (16, 24, 32, 64, 256):
            f = full_moments(N, 4, 2)
            print(f"p=4 N={N}: w2_full = {float(f[2] / f[0]):.5f}, "
                  f"N w2 = {float(N * f[2] / f[0]):.3f} (law: 64)")
