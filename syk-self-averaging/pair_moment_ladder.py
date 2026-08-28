#!/usr/bin/env python3
r"""
The mu_6 extension and the word-length ladder of the exact 1/N laws.
====================================================================

road_map.md Part III left one explicitly unlocked computation on the
tower-extraction item: extending the exact pair-moment laws to mu_6,
which needs K = 12 H-node words -- (K-1)!! = 10395 matchings per word,
m = 6 chord canonicalization.  This module runs it, plus the two gaps
the committed campaign reported rather than hid: the unidentified
A(n=3) closed form, and the still-rising double-scaling constant
c_2(p).

Everything here is the same exact machinery as pair_moment_laws.py
(rational interpolation of P(N) = 2 mu_k C(N,p)^m, certified out of
sample); this module only adds the campaign, a multiprocessing mapper,
and the ladder bookkeeping.  Three engine changes made K = 12
affordable: the per-word canonical-signature multiset is now hoisted
out of the (N, p) evaluation, the chord canonicalization runs on the
cached automorphism coset of the crossing set instead of all m!
permutations, and the DP-target reconstruction aggregates membership
patterns first.

THE LADDER.  Write K = k + 2n for the H-word length of the mu_k moment
of the A_n channel.  The committed laws showed deficits n + 2 (k = 2)
and n + 3 (k = 4), and one coefficient coincidence,
A(n=1, k=4) = A(n=2, k=2).  Both facts are faces of a K-organized
structure this campaign measures out to K = 12:

* the DEFICIT RULE candidate, deficit = K/2 + 1 -- each extra ad_H,
  and each extra pair of spectral-moment powers, would cost exactly
  one power of 1/N at leading order;

* the LADDER STEPS A(n, k) vs A(n+1, k-2) at fixed K.  The committed
  data already shows the coincidence is not rung-wide
  (A(2,4) != A(3,2) at every shared p), so the open question the new
  rungs answer is whether the n = 1 -> 2 step stays exact at k = 6,
  8, 10 while the deeper steps stay broken.

FINDINGS (all exact unless labeled as an estimate; the artifact is
pair_moment_ladder.json, pinned by test_pair_moment_ladder.py)

* THE DEFICIT RULE HOLDS EVERYWHERE.  Every law computed here and in
  the committed campaign -- k = 2..10, K = 4..12, 45 laws -- has
  deficit exactly K/2 + 1.  The mu_6 laws (n = 1 at p = 2..12(14),
  n = 2 at p = 2, 4, n = 3 at p = 2) have deficits 5, 6, 7; the whole
  K = 12 frontier at p = 2 -- (n,k) = (5,2), (4,4), (3,6), (2,8),
  (1,10) -- has deficit 7.  So mu_2 ~ N^-(n+2) is certified through
  n = 5, mu_4 ~ N^-(n+3) through n = 4, mu_6 ~ N^-(n+4) through n = 3.

* THE n = 1 -> 2 STEP IS EXACT; EVERY DEEPER STEP IS BROKEN.
  A(n=1, k) = A(n=2, k-2) at every rung and every shared p:
  K = 8 (p = 2..14), K = 10 (p = 2, 4), K = 12 (p = 2; both equal
  -155713536).  Every step with n >= 2 on the left differs, at every
  rung and every p.  At K = 10, p = 2 the coefficients over n = 1..4
  are 21 : 21 : 14 : 9 times -2^17; at K = 12, p = 2 they are
  -155713536, -155713536, -97255424, -54001664, -32505856 over
  n = 1..5.  So ad(X_1) = X_2 holds at leading 1/N order exactly,
  while the cross terms in ad(X_n) = X_{n+1} + ... enter at leading
  order for every n >= 2.  (At p = 2 the subleading B is shared across
  the n = 1 -> 2 step as well; at p = 4 it is not.  Recorded, not
  generalized.)

* A(n=3, k=2) IS IDENTIFIED.  With the p grid deepened to 14 and 16
  the minimum-degree interpolant is degree 7, and it certifies out of
  sample at p = 18:

      A(n=3, k=2; p) = -256 p^4 (p-1) (34 p^2 - 108 p + 91).

  With A(1) = -16 p^2 (p-1) and A(2) = -64 p^3 (p-1)(4p-5) the k = 2
  family reads A(n) = -2^(2n+2) p^(n+1) (p-1) Q_n(p), deg Q_n = n-1,
  and the rung degrees run 3, 5, 7 for K = 4, 6, 8.  The other K = 8
  family exceeded degree 5 too (recorded as such at p = 14); --refit
  over eight samples identifies it, certified at p = 18:

      A(n=1, k=6; p) = A(n=2, k=4; p) = -256 p^4 (p-1) (34 p^2 - 93 p + 64),

  the same -2^8 p^4 (p-1) prefactor and the same leading 34 p^2 as the
  A(3,2) quadratic, differing only below.  Every rung family through
  K = 8 is now a certified closed form.

* DOUBLE SCALING, lambda_chord = 2.  Exact c_2(p) = |mu_2| p N runs
  30.326, 30.529 (the committed p = 14, 16, reproduced) then 30.696,
  30.835, 30.953, 31.054, 31.141, 31.218 over p = 18..28, then (via
  --ds-extend) 31.344, 31.445, 31.528, 31.596, 31.654 over p = 32..48,
  monotone with shrinking increments.  The sliding-window 1/p
  extrapolation is recorded as an ESTIMATE with its window spread:
  c_2(inf) = 32.334, with the order-2, -3 and -4 windows agreeing to
  1.4e-3 and the two order-4 windows to 1e-6.  It is not 97/3 (8e-4
  below the stable value); no closed form is identified.  c_3(p)
  continues 208.33, 209.76 at p = 14, 16, still rising.

* Every K = 12 evaluation is gated end to end by the exact A_1
  flavor-sum identity at k = 10,
  (N-1) mu_10^pair(n=1) + mu_10^diag(n=1) = 0, exactly, in Fractions
  (pair = -9656064/117649, diag = 9656064/16807 at N = 8, p = 2).

Cost: the full campaign is a few hours on 8 cores (--workers).  The
default run covers everything above; --deep adds the p = 4 K = 12 laws.

Output: pair_moment_ladder.json.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from fractions import Fraction
from math import comb
from pathlib import Path

from exact_pair_moments import channel_pair_moments
from pair_moment_laws import (_lagrange_in_p, _mu_sample, _poly_eval, law,
                              mu_exact)

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# campaign plan
# ---------------------------------------------------------------------------
# (n, k, p): new laws by rung K = k + 2n.  Cheapest rung first, so a
# regression in the new k >= 6 assembly path surfaces in minutes.
PLAN_K8 = [(1, 6, p) for p in (2, 4, 6, 8, 10, 12)] \
        + [(2, 4, p) for p in (6, 8, 10, 12)] \
        + [(3, 2, p) for p in (14, 16)]
PLAN_K10 = [(1, 8, 2), (1, 8, 4), (2, 6, 2), (2, 6, 4), (3, 4, 2), (3, 4, 4)]
PLAN_K12 = [(1, 10, 2), (2, 8, 2), (3, 6, 2), (4, 4, 2), (5, 2, 2)]
PLAN_K12_DEEP = [(5, 2, 4), (3, 6, 4)]

# double-scaling extension on lambda_chord = 2, i.e. N = p^2.  p = 14, 16
# for n = 2 overlap the committed exact scan (consistency gate); the rest
# is new territory.
DS_N2_PS = (14, 16, 18, 20, 22, 24, 26, 28)
DS_N3_PS = (14, 16)


def _committed_laws() -> list:
    """The committed pair_moment_laws.json rows, for rung bookkeeping and
    closed-form fits that pool old and new samples."""
    with open(HERE / "pair_moment_laws.json") as fh:
        return json.load(fh)["laws"]


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------
def a1_identity_k10(N: int = 8, p: int = 2) -> dict:
    """The exact end-to-end gate of the K = 12 machinery: A_1 = -2pH per
    instance forces (N-1) mu_k^pair(n=1) + mu_k^diag(n=1) = 0 at every
    N -- here at k = 10, whose words are the full K = 12 scope.  The
    residual must be EXACTLY zero in Fractions."""
    pair = channel_pair_moments(N, p, 1, ks=(10,), exact=True)["mu10"]
    diag = channel_pair_moments(N, p, 1, marks=('i', 'i', 'i', 'i'),
                                ks=(10,), exact=True)["mu10"]
    resid = (N - 1) * pair + diag
    return dict(N=N, p=p, k=10, exact_zero=bool(resid == 0),
                residual=str(resid), pair=str(pair), diag=str(diag))


# ---------------------------------------------------------------------------
# ladder bookkeeping
# ---------------------------------------------------------------------------
def rung_table(all_laws: list) -> dict:
    """Group the exact leading data by rung K = k + 2n: deficits and the
    adjacent-step coefficient comparisons A(n, k) vs A(n+1, k-2)."""
    by_key = {(r["n"], r["k"], r["p"]): r for r in all_laws}
    rungs: dict = {}
    for (n, k, p), r in by_key.items():
        K = k + 2 * n
        rungs.setdefault(K, []).append(
            dict(n=n, k=k, p=p, deficit=r["deficit"], A=r["A"]))
    checks = []
    for (n, k, p), r in sorted(by_key.items()):
        other = by_key.get((n + 1, k - 2, p))
        if other is not None and k >= 4:
            checks.append(dict(
                K=k + 2 * n, p=p,
                lhs=f"A(n={n},k={k})", rhs=f"A(n={n + 1},k={k - 2})",
                lhs_A=r["A"], rhs_A=other["A"],
                equal=bool(Fraction(r["A"]) == Fraction(other["A"]))))
    deficit_rule = all(row["deficit"] == K // 2 + 1
                       for K, rows in rungs.items() for row in rows)
    return dict(deficit_rule_K_half_plus_1=bool(deficit_rule),
                rungs={str(K): sorted(rows, key=lambda r: (r["n"], r["p"]))
                       for K, rows in sorted(rungs.items())},
                adjacent_step_checks=checks)


def a_closed_form_pooled(n: int, k: int, rows: list, p_check: int,
                         mapper=None):
    """(form, fresh_law): closed form of A(p) at fixed (n, k) from the
    pooled sample rows (committed + new), certified by a FRESH law()
    extraction at an out-of-sample p.  Mirrors
    pair_moment_laws.a_closed_form, with the (n, k) channel explicit;
    the fresh check law is returned too, since it is real data."""
    pts = sorted(((r["p"], Fraction(r["A"])) for r in rows
                  if r["n"] == n and r["k"] == k), key=lambda t: t[0])
    if len(pts) < 4:
        return dict(n=n, k=k, identified=False,
                    reason="fewer than 4 p samples"), None
    poly = _lagrange_in_p(list(pts))
    out = dict(n=n, k=k, p_samples=[int(x) for x, _ in pts],
               A_poly_coeffs=[str(c) for c in poly],
               degree=len(poly) - 1, p_check=p_check)
    fresh = law(p_check, n, k, mapper=mapper)
    predicted = _poly_eval(poly, Fraction(p_check))
    out["identified"] = bool(predicted == Fraction(fresh["A"]))
    out["A_at_p_check"] = fresh["A"]
    if not out["identified"]:
        out["reason"] = (f"degree > {len(poly) - 1}: predicted {predicted} "
                         f"vs computed {fresh['A']} at p = {p_check}")
    return out, fresh


# ---------------------------------------------------------------------------
# double scaling
# ---------------------------------------------------------------------------
def ds_extension(n: int, ps, mapper=None) -> list:
    """Exact mu_2 on lambda_chord = 2 (N = p^2), c(p) = |mu_2| p N."""
    args = [(p * p, p, n, 2) for p in ps]
    if mapper is None:
        mus = [mu_exact(*a) for a in args]
    else:
        mus = list(mapper(_mu_sample, args))
    rows = []
    for p, mu in zip(ps, mus):
        rows.append(dict(p=p, N=p * p, mu2=float(mu),
                         c=float(abs(mu) * p * p * p)))
        print(f"  lambda_chord=2, n={n}: p={p:2d} N={p * p:4d}  "
              f"mu2={float(mu):+.6e}  c={rows[-1]['c']:.3f}", flush=True)
    return rows


def _window_fit(window: list) -> float:
    """Exact Vandermonde solve of c(p) = x0 + x1/p + ... through the
    (p, c) window; returns float(x0).  Window length sets the order."""
    F = Fraction
    npts = len(window)
    rows_ = [[F(1, p ** j) for j in range(npts)] + [F(c)]
             for (p, c) in window]
    for i in range(npts):
        piv = rows_[i][i]
        rows_[i] = [v / piv for v in rows_[i]]
        for j in range(npts):
            if j != i and rows_[j][i]:
                f = rows_[j][i]
                rows_[j] = [a - f * b for a, b in zip(rows_[j], rows_[i])]
    return float(rows_[0][npts])


def c_limit_estimate(rows: list) -> dict:
    """Sliding-window 1/p extrapolation of c(p) = c_inf + a/p + ...:
    orders 2..4 (window lengths 3..5) on the deepest windows, with the
    spread across windows and orders as the honest error bar.  An
    ESTIMATE, not a law; a ln(p)/p term in the true asymptotics would
    bias it, and the artifact records the raw exact points so anyone can
    refit."""
    pts = sorted((r["p"], r["c"]) for r in rows)
    if len(pts) < 4:
        return dict(estimated=False, reason="too few points")
    fits = []
    for width in (3, 4, 5):
        if len(pts) < width + 1:
            continue
        fits.append(dict(order=width - 1,
                         window=[p for p, _ in pts[-width:]],
                         c_inf=_window_fit(pts[-width:])))
        fits.append(dict(order=width - 1,
                         window=[p for p, _ in pts[-width - 1:-1]],
                         c_inf=_window_fit(pts[-width - 1:-1])))
    vals = [f["c_inf"] for f in fits]
    return dict(estimated=True, fits=fits,
                c_inf_central=fits[-2]["c_inf"] if fits else None,
                spread=max(vals) - min(vals),
                note="sliding-window 1/p fits, orders 2-4; estimate, "
                     "not exact; central = deepest widest window")


# ---------------------------------------------------------------------------
# campaign
# ---------------------------------------------------------------------------
def run(deep: bool = False, workers: int = None, out_path=None) -> dict:
    t0 = time.time()
    workers = workers or max(1, mp.cpu_count() - 1)
    ctx = mp.get_context("spawn")
    pool = ctx.Pool(processes=workers)
    mapper = pool.map
    out = Path(out_path) if out_path else HERE / "pair_moment_ladder.json"
    payload = dict(provenance=(
        "pair_moment_ladder.py -- the mu_6 / K = 12 extension of the "
        "exact 1/N pair-moment laws, the word-length-ladder structure "
        "(deficit = K/2 + 1; the n = 1 -> 2 coefficient step exact, the "
        "others not), the A(n=3) and A(n=1,k=6) closed forms, and the "
        "double-scaling c_2(p) extension; same exact machinery as "
        "pair_moment_laws.py, certified out of sample per law"))

    def checkpoint():
        payload["runtime_s"] = round(time.time() - t0, 1)
        with open(out, "w") as fh:
            json.dump(payload, fh, indent=1)

    try:
        print(f"[gate] A_1 flavor-sum identity at k=10 (K=12 words) ...",
              flush=True)
        payload["a1_gate_k10"] = a1_identity_k10()
        if not payload["a1_gate_k10"]["exact_zero"]:
            raise RuntimeError("K=12 gate failed: "
                               + payload["a1_gate_k10"]["residual"])
        print("  exact zero: True", flush=True)
        checkpoint()

        laws = []
        phases = [("K=8", PLAN_K8), ("K=10", PLAN_K10), ("K=12", PLAN_K12)]
        if deep:
            phases.append(("K=12 deep", PLAN_K12_DEEP))
        for phase, plan in phases:
            for (n, k, p) in plan:
                print(f"[{phase}] law n={n} k={k} p={p} ...", flush=True)
                row = law(p, n, k, mapper=mapper)
                laws.append(row)
                payload["laws"] = laws
                print(f"  deficit={row['deficit']}  A={row['A']}  "
                      f"({row['runtime_s']} s)", flush=True)
                checkpoint()

        committed = _committed_laws()
        payload["ladder"] = rung_table(committed + laws)
        checkpoint()

        # closed forms: A(1,6) over its new p grid, checked at p = 14;
        # A(2,4) over committed + new, checked at p = 14; A(3,2) over
        # committed + new (p up to 16), checked at p = 18.
        forms = []
        for (n, k, p_chk) in ((1, 6, 14), (2, 4, 14), (3, 2, 18)):
            print(f"[closed form] A(n={n},k={k}), out-of-sample p={p_chk} "
                  f"...", flush=True)
            form, fresh = a_closed_form_pooled(n, k, committed + laws,
                                               p_chk, mapper=mapper)
            if fresh is not None:       # the fresh check law is real data
                laws.append(fresh)
            payload["laws"] = laws
            forms.append(form)
            payload["A_closed_forms"] = forms
            print(f"  identified={form['identified']} "
                  f"degree={form.get('degree')}", flush=True)
            checkpoint()

        # refresh the rung table with the check laws included
        payload["ladder"] = rung_table(committed + laws)

        print("[double scaling] n=2 extension ...", flush=True)
        ds2 = ds_extension(2, DS_N2_PS, mapper=mapper)
        payload["double_scaling"] = dict(n2=ds2)
        payload["double_scaling"]["c2_limit"] = c_limit_estimate(ds2)
        checkpoint()
        print("[double scaling] n=3 extension ...", flush=True)
        ds3 = ds_extension(3, DS_N3_PS, mapper=mapper)
        payload["double_scaling"]["n3"] = ds3
        checkpoint()
    finally:
        pool.close()
        pool.join()

    checkpoint()
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


def ds_extend(ps, workers: int = None, out_path=None) -> dict:
    """Deepen the n = 2 double-scaling scan in an EXISTING artifact:
    compute the new exact points, merge, refit the limit estimate, and
    rewrite.  Used to push c_2(p) beyond the default campaign grid
    without re-running the laws."""
    out = Path(out_path) if out_path else HERE / "pair_moment_ladder.json"
    with open(out) as fh:
        payload = json.load(fh)
    have = {r["p"] for r in payload["double_scaling"]["n2"]}
    todo = tuple(p for p in ps if p not in have)
    workers = workers or max(1, mp.cpu_count() - 1)
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=min(workers, max(1, len(todo)))) as pool:
        rows = ds_extension(2, todo, mapper=pool.map)
    payload["double_scaling"]["n2"].extend(rows)
    payload["double_scaling"]["n2"].sort(key=lambda r: r["p"])
    payload["double_scaling"]["c2_limit"] = c_limit_estimate(
        payload["double_scaling"]["n2"])
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"extended n2 double scaling to p={max(have | set(ps))}; "
          f"c_inf estimate: {payload['double_scaling']['c2_limit']}")
    return payload


def refit_closed_form(n: int, k: int, sample_ps, p_check: int,
                      workers: int = None, out_path=None) -> dict:
    """Deepen the p grid of one A(p) family in an EXISTING artifact and
    re-certify: compute any missing sample laws, fit the minimum-degree
    interpolant through all of them, certify by a fresh law at p_check,
    replace the family's A_closed_forms entry, and rewrite.  Used when a
    family's degree exceeds the campaign grid (the K = 8 rung needs eight
    samples)."""
    out = Path(out_path) if out_path else HERE / "pair_moment_ladder.json"
    with open(out) as fh:
        payload = json.load(fh)
    laws = payload["laws"]
    have = {(r["n"], r["k"], r["p"]) for r in laws} | {
        (r["n"], r["k"], r["p"]) for r in _committed_laws()}
    missing = [p for p in sample_ps if (n, k, p) not in have]
    workers = workers or max(1, mp.cpu_count() - 1)
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=workers) as pool:
        for p in missing:
            print(f"[refit] sample law n={n} k={k} p={p} ...", flush=True)
            row = law(p, n, k, mapper=pool.map)
            laws.append(row)
            print(f"  deficit={row['deficit']}  A={row['A']}  "
                  f"({row['runtime_s']} s)", flush=True)
        rows = [r for r in _committed_laws() + laws
                if r["n"] == n and r["k"] == k and r["p"] in set(sample_ps)]
        print(f"[refit] closed form A(n={n},k={k}) over p={sorted(sample_ps)}"
              f", out-of-sample p={p_check} ...", flush=True)
        form, fresh = a_closed_form_pooled(n, k, rows, p_check,
                                           mapper=pool.map)
    if fresh is not None:
        laws.append(fresh)
    forms = [f for f in payload.get("A_closed_forms", [])
             if not (f["n"] == n and f["k"] == k)]
    forms.append(form)
    payload["A_closed_forms"] = forms
    payload["laws"] = laws
    payload["ladder"] = rung_table(_committed_laws() + laws)
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"  identified={form['identified']} degree={form.get('degree')}",
          flush=True)
    return form


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--deep", action="store_true",
                    help="add the p = 4 K = 12 laws (hours)")
    ap.add_argument("--workers", type=int, default=None,
                    help="pool size (default: cpu_count - 1)")
    ap.add_argument("--ds-extend", type=str, default=None,
                    help="comma-separated extra p values for the n = 2 "
                         "double-scaling line, merged into the existing "
                         "artifact (e.g. 32,36,40,44,48)")
    ap.add_argument("--refit", type=str, default=None,
                    help="n,k:p1,p2,...:p_check -- deepen one A(p) family's "
                         "sample grid and re-certify (e.g. "
                         "1,6:2,4,6,8,10,12,14,16:18)")
    args = ap.parse_args()
    if args.ds_extend:
        ds_extend(tuple(int(x) for x in args.ds_extend.split(",")),
                  workers=args.workers)
    elif args.refit:
        nk, ps, pc = args.refit.split(":")
        n_, k_ = (int(x) for x in nk.split(","))
        refit_closed_form(n_, k_, tuple(int(x) for x in ps.split(",")),
                          int(pc), workers=args.workers)
    else:
        run(deep=args.deep, workers=args.workers)
