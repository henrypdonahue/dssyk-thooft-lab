#!/usr/bin/env python3
r"""
The FULL A_2 singlet channel at N = infinity: an exact chord law.
=================================================================

The conserved-channel theorem (chord_charges.py) freezes the i != j pair
amplitude of A_n = sum_i psi_i (ad_H)^n psi_i at N = infinity.  The
diagonal channel (diagonal_chord.py) does not freeze.  Each of those is
a per-flavor object; the physical channel is their flavor sum,

    <A_n ad^k A_n> = N mu_k^diag + N(N-1) mu_k^pair ,

and until now its N = infinity value was a Richardson extrapolation of
the exact finite-N series, w2 ~ 1.31(5) at lambda_chord = 2, with "the
pair channel's 1/N coefficient" named as the missing chord-side piece.

It is not missing.  Sum the flavor FIRST.  The exact operator identity
(syk-self-averaging/conserved_points.py)

    A_2 = 4 sum_i H_i^2 = 4 sum_{c,c'} |c cap c'| J_c J_c' psi_c psi_c'

removes the fermion marks altogether: A_2 is H^2 with each pair of
H-chords weighted by the number of Majorana indices they SHARE.  In the
double-scaled limit the pairwise intersections |c cap c'| of distinct
random p-subsets are independent Poisson(lambda/2) variables (that is
where q = E[(-1)^X] = e^{-lambda} comes from), so every full-channel word
<Tr[H^a A_2 H^b A_2]> is a finite sum over chord matchings of the a + b + 4
H-nodes with the standard crossing weight q per crossing pair, and for the
two weighted pairs (the two chords leaving each A_2 insertion) the
intersection-moment factors

    E[X (-1)^{sX}]   = (lambda/2) (-q)^s ,
    E[X^2 (-1)^{sX}] = mu + mu^2  (s = 0),   q (mu^2 - mu)  (s = 1),
    mu = lambda/2, s = 1 iff the two chords cross,

the second line for the matchings that contract both slots of one
insertion into both slots of the other.  A self-contraction inside one
insertion is the c-number 4 sum_c p J_c^2, i.e. <A_2>, which the traceless
moments discard.  The mu_k assembly is the same binomial ad-expansion as
chord_moments.py.  Everything is a polynomial in lambda and q:

    mu_0 = 4 lambda [lambda (1+q) + 2 (1-q)]
    mu_2 = 16 lambda (1-q) [1 - q^2 + 2 lambda q]
    mu_4 = lambda^2 (96q^5 - 128q^3 - 128q^2 + 160q)
           + lambda (-48q^6 + 112q^4 + 64q^3 - 144q^2 - 64q + 80)

At lambda_chord = 2: mu_0 = 32 exactly, mu_2 = 42.141005,

    w2_full(lambda = 2) = 1.316906 ,   kurt = 3.852 .

GATES (test_full_channel_chord.py).
  * The exact finite-N full channel on lambda_chord = 2 (exact Fractions,
    syk-self-averaging engine, p = 4..40) extrapolates in 1/p to
    mu_0 = 32.0002, mu_2 = 42.1410, w2 = 1.31690; on lambda_chord = 1
    (p <= 20) to 10.5282, 16.1853, 1.53737 against 10.5285, 16.1866,
    1.53741.  The p <= 10 Richardson estimate 1.31(5) is superseded.
  * Semiclassics: w2 = 4 lambda - 4 lambda^2 + ..., kurt -> 7/2 as
    lambda -> 0, i.e. exactly w2 = 2 lambda Delta, kurt = 3 + 1/Delta at
    Delta = 2: the full A_2 channel is a dimension-2 matter chord in the
    large-q limit, as psi d^2 psi should be.
  * The operator identity itself, dense at (8,4), (10,4), (12,4).

THE CONSERVED ADMIXTURE, RESOLVED (singlet_tower.py).  Mean subtraction
removes the identity from A_2, not its overlap with H^2.  That overlap
is exactly

    <A_2 (H^2 - 1)> = 2 lambda (1 - q),      <(H^2 - 1)^2> = 1 + q,

and what is left after removing it is, word for word, a matter chord of
dimension Delta = 2 (a theorem; the proof is a sign-reversing involution
on the matchings and sits in singlet_tower.py's docstring):

    mu_0 = 4 lambda^2 (1-q)^2 / (1+q)          [conserved, H^2 - 1]
         + 8 lambda (1 - q^2 + 2 lambda q) / (1+q)   [dynamical],
    w2_dyn = mu_2 / mu_0^dyn = 2 (1 - q^2)     exactly, every lambda.

So the w2 above is the Delta = 2 matter-chord width 2(1 - q^Delta)
diluted by conserved weight: 1.963369 against 1.316906 at
lambda_chord = 2, where 32.9 % of mu_0 is the H^2 component.  The
semiclassical gate below (w2 -> 4 lambda, kurt -> 7/2) is the small-
lambda face of an identity that holds at every lambda.
dynamic_channel() carries the split.

WHAT IT MEANS FOR THE FLAVOR SPLIT.  Per flavor, N mu^diag and
N(N-1) mu^pair are each O(sqrt N) on a double-scaling line (they are
O(1/p) per flavor times N ~ p^2/lambda) and cancel identically at that
order -- the pair channel's double-scaling constant IS the diagonal
channel's Delta-slope (syk-self-averaging/double_scaling_law.py).  The
O(1) remainder is the full channel above, and it belongs to neither half.
"All of the tower's N = infinity weight is diagonal" was that bookkeeping
read as physics; the physical statement is this closed form.

SCOPE.  n = 2, strict N = infinity at fixed lambda; corrections are
O(1/p) (hypergeometric against Poisson), so at fixed p the exact Wick
engine remains the tool (syk-self-averaging/full_channel_laws.py: the
fixed-p width dies as 8p(p-2)/N).  For n >= 3 the same construction
works once the conserved admixtures of H, H^2, ... are projected out;
that is singlet_tower.py, which finds every A_n computed (n <= 6) to be
a polynomial in H plus matter chords of dimension 2..n.

Output: full_channel_chord.json (the closed forms, a lambda scan, the
two gate series with their 1/p limits), ~1 min for the finite-N gates.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# the matching sum
# ---------------------------------------------------------------------------
def _matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for k in range(1, len(items)):
        b = items[k]
        rest = items[1:k] + items[k + 1:]
        for rem in _matchings(rest):
            yield [(a, b)] + rem


def _crosses(c1, c2) -> bool:
    (a, b), (c, d) = sorted(c1), sorted(c2)
    return (a < c < b) != (a < d < b)


def word_amplitude(a: int, b: int, lam, q):
    """<Tr[H^a A_2 H^b A_2]> in chord units at N = infinity (traceless,
    i.e. without the in-insertion self-contractions), for numeric or
    sympy lam, q.  The two weighted pairs carry the Poisson intersection
    moments; every other crossing pair carries q."""
    mu = lam / 2
    nodes = ['H'] * a + ['c1', 'c2'] + ['H'] * b + ['d1', 'd2']
    idx = list(range(len(nodes)))
    pos = {nodes[i]: i for i in idx if nodes[i] != 'H'}
    tot = 0
    for m in _matchings(idx):
        chords = [tuple(sorted(pr)) for pr in m]
        of = {}
        for ci, (x, y) in enumerate(chords):
            of[x] = ci
            of[y] = ci
        e = tuple(sorted((of[pos['c1']], of[pos['c2']])))
        f = tuple(sorted((of[pos['d1']], of[pos['d2']])))
        if e[0] == e[1] or f[0] == f[1]:
            continue                      # self-contraction: <A_2> part
        crossing = set()
        for i in range(len(chords)):
            for j in range(i + 1, len(chords)):
                if _crosses(chords[i], chords[j]):
                    crossing.add((i, j))
        w = 16
        for g in crossing:
            if g not in (e, f):
                w = w * q
        if e == f:
            w = w * ((q * (mu * mu - mu)) if e in crossing else (mu + mu * mu))
        else:
            for g in (e, f):
                w = w * ((-mu * q) if g in crossing else mu)
        tot = tot + w
    return tot


def mu_full(k: int, lam, q):
    """mu_k of the traceless full A_2 channel (both orderings coincide by
    cyclicity for the flavor-summed object)."""
    return sum((-1) ** s * comb(k, s) * word_amplitude(k - s, s, lam, q)
               for s in range(k + 1))


# ---------------------------------------------------------------------------
# closed forms
# ---------------------------------------------------------------------------
def mu0_closed(lam, q):
    return 4 * lam * (lam * (1 + q) + 2 * (1 - q))


def mu2_closed(lam, q):
    return 16 * lam * (1 - q) * (1 - q * q + 2 * lam * q)


def mu4_closed(lam, q):
    return (lam ** 2 * (96 * q ** 5 - 128 * q ** 3 - 128 * q ** 2 + 160 * q)
            + lam * (-48 * q ** 6 + 112 * q ** 4 + 64 * q ** 3 - 144 * q ** 2
                     - 64 * q + 80))


def channel(lam: float) -> dict:
    q = math.exp(-lam)
    m0, m2, m4 = mu0_closed(lam, q), mu2_closed(lam, q), mu4_closed(lam, q)
    return dict(lam=lam, q=q, mu0=m0, mu2=m2, mu4=m4, w2=m2 / m0,
                kurt=m4 * m0 / m2 ** 2)


def dynamic_channel(lam: float) -> dict:
    """The full A_2 channel with its conserved H^2 component removed
    (the theorem of singlet_tower.py): what remains is a Delta = 2 matter
    chord, so its width is 2(1 - q^2) and its kurtosis is the O_Delta
    closed form at Delta = 2."""
    q = math.exp(-lam)
    c = channel(lam)
    stat = 4 * lam * lam * (1 - q) ** 2 / (1 + q)
    dyn = 8 * lam * (1 - q * q + 2 * lam * q) / (1 + q)
    return dict(lam=lam, q=q, mu0=c["mu0"], mu0_stationary=stat,
                mu0_dynamic=dyn, stationary_fraction=stat / c["mu0"],
                h2_coefficient=2 * lam * (1 - q) / (1 + q),
                w2_raw=c["w2"], w2_dynamic=c["mu2"] / dyn,
                kurt_dynamic=c["mu4"] * dyn / c["mu2"] ** 2)


def symbolic():
    """The closed forms re-derived from the matching sum in sympy."""
    import sympy as sp
    lam, q = sp.symbols("lam q", positive=True)
    out = {}
    for k in (0, 2, 4):
        out[k] = sp.expand(mu_full(k, lam, q))
    checks = dict(
        mu0=sp.simplify(out[0] - mu0_closed(lam, q)) == 0,
        mu2=sp.simplify(out[2] - mu2_closed(lam, q)) == 0,
        mu4=sp.simplify(out[4] - mu4_closed(lam, q)) == 0)
    w2 = sp.series((out[2] / out[0]).subs(q, sp.exp(-lam)), lam, 0, 3)
    kurt = sp.series((out[4] * out[0] / out[2] ** 2).subs(q, sp.exp(-lam)),
                     lam, 0, 2)
    return dict(mu={k: str(v) for k, v in out.items()}, closed_form_checks=checks,
                w2_small_lambda=str(w2), kurt_small_lambda=str(kurt))


# ---------------------------------------------------------------------------
# gates: the exact finite-N full channel on a double-scaling line
# ---------------------------------------------------------------------------
def finite_n_full(N: int, p: int, ks=(0, 2)) -> dict:
    sys.path.insert(0, str(HERE.parent / "syk-self-averaging"))
    import exact_pair_moments as epm
    pr = epm.channel_pair_moments(N, p, 2, marks=('i', 'i', 'j', 'j'),
                                  ks=ks, exact=True)
    dg = epm.channel_pair_moments(N, p, 2, marks=('i', 'i', 'i', 'i'),
                                  ks=ks, exact=True)
    epm._DP_CACHE.clear()
    epm._EXP_CACHE.clear()
    return {k: N * dg[f"mu{k}"] + N * (N - 1) * pr[f"mu{k}"] for k in ks}


def line_series(lam: float, ps, progress: bool = True) -> list:
    rows = []
    for p in ps:
        N = 2 * p * p / lam
        if abs(N - round(N)) > 1e-9:
            raise ValueError(f"2p^2/lambda not integer at p={p}")
        N = int(round(N))
        t0 = time.time()
        f = finite_n_full(N, p)
        rows.append(dict(p=p, N=N, mu0=float(f[0]), mu2=float(f[2]),
                         w2=float(f[2] / f[0]), mu0_exact=str(f[0]),
                         mu2_exact=str(f[2]), runtime_s=round(time.time() - t0, 1)))
        if progress:
            print(f"  lam={lam:g} p={p:2d} N={N:4d}: mu0={rows[-1]['mu0']:.5f} "
                  f"mu2={rows[-1]['mu2']:.5f} w2={rows[-1]['w2']:.6f}",
                  flush=True)
    return rows


def inverse_p_limit(rows: list, key: str, orders=(3, 4, 5), p_min: int = 8) -> dict:
    """Least-squares polynomial in 1/p; central = highest order with the
    smallest residual, spread over orders.  An estimate."""
    import numpy as np
    pts = [(r["p"], r[key]) for r in rows if r["p"] >= p_min]
    x = np.array([1.0 / p for p, _ in pts])
    y = np.array([v for _, v in pts])
    fits = []
    for o in orders:
        if len(pts) <= o + 1:
            continue
        A = np.vander(x, o + 1, increasing=True)
        c = np.linalg.lstsq(A, y, rcond=None)[0]
        fits.append(dict(order=o, limit=float(c[0]),
                         max_residual=float(np.max(np.abs(A @ c - y)))))
    best = min(fits, key=lambda f: f["max_residual"])
    return dict(fits=fits, central=best["limit"],
                spread=max(f["limit"] for f in fits) - min(f["limit"] for f in fits))


def run(out_path=None, ps2=(4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32),
        ps1=(4, 6, 8, 10, 12, 14, 16, 18)) -> dict:
    t0 = time.time()
    print("symbolic closed forms", flush=True)
    sym = symbolic()
    lams = [4.0, 3.2, 2.6667, 2.2857, 2.0, 1.7778, 1.6, 1.0, 0.5, 0.25, 0.1,
            0.05, 0.02]
    scan = [channel(l) for l in lams]
    gates = {}
    for lam, ps in ((2.0, ps2), (1.0, ps1)):
        print(f"finite-N gate on lambda_chord = {lam}", flush=True)
        rows = line_series(lam, ps)
        cf = channel(lam)
        gates[f"lam{lam:g}"] = dict(
            lam=lam, rows=rows,
            limits={key: inverse_p_limit(rows, key) for key in ("mu0", "mu2", "w2")},
            closed_form=dict(mu0=cf["mu0"], mu2=cf["mu2"], w2=cf["w2"]))
        for key in ("mu0", "mu2", "w2"):
            lim = gates[f"lam{lam:g}"]["limits"][key]
            print(f"  {key}: 1/p limit {lim['central']:.6f} (spread "
                  f"{lim['spread']:.1e})  closed form {cf[key]:.6f}", flush=True)
    payload = dict(
        provenance="full_channel_chord.py -- the FULL A_2 singlet channel at "
                   "N = infinity by the intersection-weighted chord rule "
                   "(A_2 = 4 sum |c cap c'| J J psi psi, Poisson(lambda/2) "
                   "intersections), closed forms in (lambda, q), sympy-checked; "
                   "gated by the exact finite-N full channel on two "
                   "double-scaling lines and by the semiclassical limit",
        runtime_s=round(time.time() - t0, 1), symbolic=sym, scan=scan,
        gates=gates,
        closed_forms=dict(mu0="4 lam [lam (1+q) + 2 (1-q)]",
                          mu2="16 lam (1-q) [1 - q^2 + 2 lam q]",
                          mu4="lam^2 (96q^5 - 128q^3 - 128q^2 + 160q) + lam "
                              "(-48q^6 + 112q^4 + 64q^3 - 144q^2 - 64q + 80)"),
        lambda2=channel(2.0))
    out = Path(out_path) if out_path else HERE / "full_channel_chord.json"
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
        for lam in (4.0, 2.0, 1.0, 0.5, 0.1):
            c = channel(lam)
            d = dynamic_channel(lam)
            print(f"lambda_chord={lam:<4}: mu0={c['mu0']:.6f} mu2={c['mu2']:.6f} "
                  f"w2={c['w2']:.6f} kurt={c['kurt']:.5f}   conserved H^2 share "
                  f"{d['stationary_fraction']:.4f}, dynamical w2={d['w2_dynamic']:.6f} "
                  f"kurt={d['kurt_dynamic']:.5f}")
