#!/usr/bin/env python3
r"""
The singlet channels in the 't Hooft limit: exact finite-p shape.
================================================================

Miyashita-Sekino-Susskind compare with the 't Hooft model at fixed
p^2 = alpha-bar, N -> infinity, "since that's where the calculations are
done" (2607.05678 Sec. 4.1; 2506.18054v3 Sec. 5).  Everything at
N = infinity in duality/ is at fixed lambda instead.  This module gives
the exact shape of the Majorana singlet channels in the fixed-p limit and
compares it with the lambda -> 0 ladders of duality/HORIZON.md.

THE FIXED-p SIDE (exact).  full_channel_laws.full_law interpolates the
exact finite-N moments of the full A_n channel in N and certifies the
closed form out of sample.  Its leading coefficients A_k (mu_k ~ A_k
N^-(d_0 + k/2)) give the N = infinity, fixed-p moments in physical units
(coupling Var(J) = p!/N^(p-1), J = 1):

    <omega^k> / <omega^0> = A_k / A_0       (<H^2>/dim -> N).

For A_2 the identity-subtracted mu_0 carries no conserved weight in this
limit (the H^2 share falls like N^-2, full_channel_laws.py), so these are
the moments of the propagating channel.  For A_3 mu_0 keeps a conserved
polynomial in H, so the A_3 entry uses mu_4/mu_2, which no conserved
piece reaches.

THE LADDER SIDE.  As lambda -> 0 the A_2 channel is a dimension-2 matter
chord and A_3 tends to the dimension-3 chord O_3 (duality/singlet_tower.py),
each the ladder sech^(2 Delta)(scriptJ t) with scriptJ^2 = 2 p^2 J^2 in
this Majorana normalization.  Its cumulants are 2 Delta, 4 Delta, 32 Delta
(times scriptJ^2, ^4, ^6), so <omega^2> = 2 Delta scriptJ^2,
<omega^4> = (4 Delta + 12 Delta^2) scriptJ^4, <omega^6> = (32 Delta +
120 Delta^2 + 120 Delta^3) scriptJ^6, and <omega^4>/<omega^2> =
4 p^2 (1 + 3 Delta): 28 p^2 at Delta = 2, 40 p^2 at Delta = 3.

RESULTS (exact rationals; see LAWS for what each form was fit on and
checked on):

    A_2:  <omega^2> = 8 p (p - 2)                  (ladder 8 p^2)
          <omega^4> = 32 p^2 (p - 2)(7 p - 15)     (ladder 224 p^4)
          kurtosis  = 7/2 - 1/(2 (p - 2))          (ladder 7/2)
          mu6/mu2^3 = 47/2 - 37/(4(p-2)) + 5/(4(p-2)^2)   (the limit 47/2 is
                      imposed in the fit, not found)
    A_3:  <omega^4>/<omega^2> = 4 p (10 p - 21)    (a quadratic through the
                      computed p, checked at p = 10; at those p the value is
                      still between the Delta = 2 and Delta = 3 ladders)

HOW TO READ IT.  Both limits sit at p -> infinity with p^2/N -> 0, so
agreement at leading order in 1/p is the expected large-p <-> small-lambda
correspondence, and this is a consistency check of it, not news about the
paper's finite alpha-bar.  The new content is the exact size of the
finite-p corrections, and they are large: <omega^2> is (1 - 2/p) of its
ladder value (half at p = 4), <omega^4> is (1 - 2/p)(1 - 15/(7p)) of it
(0.68 at p = 12).  At the p where a comparison with the 't Hooft model
would be made the channels are the ladders only roughly.  Moments do not
fix analytic structure, so none of this decides whether the finite-p
spectrum is gapped.

Output: thooft_limit.json, which stores every exact leading coefficient.
`python thooft_limit.py` prints the tables from it; `--run` recomputes all
laws from scratch (several hours: the A_3 fourth moments and the A_2
sixth moments at p = 8 and 10 dominate); `--from-cache FILE` assembles the
artifact from coefficients full_law already produced.
"""

from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path

from full_channel_laws import full_law

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# the ladder side
# ---------------------------------------------------------------------------
def ladder_moments(delta, p) -> dict:
    """<omega^2>, <omega^4>, <omega^6> of sech^(2 Delta)(scriptJ t) with
    scriptJ^2 = 2 p^2 (physical units, J = 1)."""
    J2 = Fraction(2 * p * p)
    d = Fraction(delta)
    return {2: 2 * d * J2,
            4: (4 * d + 12 * d * d) * J2 ** 2,
            6: (32 * d + 120 * d * d + 120 * d ** 3) * J2 ** 3}


# ---------------------------------------------------------------------------
# the fixed-p side
# ---------------------------------------------------------------------------
def leading(p: int, n: int, k: int) -> dict:
    law = full_law(p, n, k)
    return dict(deficit=law["deficit"], A=Fraction(law["A"]),
                out_of_sample_N=law["out_of_sample_N"])


def a2_shape(p: int, with_mu6: bool = False) -> dict:
    """Exact fixed-p, N -> infinity moments of the full A_2 channel."""
    L = {k: leading(p, 2, k) for k in ((0, 2, 4, 6) if with_mu6 else (0, 2, 4))}
    d0 = L[0]["deficit"]
    for k, v in L.items():
        if v["deficit"] != d0 + k // 2:
            raise RuntimeError(f"deficit pattern broken at p={p}, k={k}")
    om = {k: L[k]["A"] / L[0]["A"] for k in L if k}
    out = dict(p=p, omega2=om[2], omega4=om[4], kurtosis=om[4] / om[2] ** 2)
    if 6 in om:
        out["omega6"] = om[6]
        out["hyper"] = om[6] / om[2] ** 3
    return out


def a3_ratio(p: int) -> dict:
    """<omega^4>/<omega^2> of the full A_3 channel at fixed p, N -> infinity:
    the ratio of the leading mu_4 and mu_2 coefficients (no mu_0, so no
    conserved piece enters)."""
    l2, l4 = leading(p, 3, 2), leading(p, 3, 4)
    if l4["deficit"] != l2["deficit"] + 1:
        raise RuntimeError(f"A_3 deficits {l2['deficit']}, {l4['deficit']} at p={p}")
    return dict(p=p, ratio=l4["A"] / l2["A"])


# ---------------------------------------------------------------------------
# the closed forms (fit on the first p values, certified on the rest)
# ---------------------------------------------------------------------------
def omega2_closed(p):
    return 8 * p * (p - 2)


def omega4_closed(p):
    return 32 * p * p * (p - 2) * (7 * p - 15)


def kurtosis_closed(p):
    return Fraction(7, 2) - Fraction(1, 2 * (p - 2))


def hyper_closed(p):
    x = p - 2
    return Fraction(47, 2) - Fraction(37, 4 * x) + Fraction(5, 4 * x * x)


def a3_ratio_closed(p):
    return 4 * p * (10 * p - 21)


LAWS = dict(
    omega2=dict(form="8 p (p-2)", fit_on=[2, 4, 6], checked_on=[8, 10, 12],
                source="full_channel_laws.py"),
    kurtosis=dict(form="7/2 - 1/(2 (p-2))", fit_on=[4, 6], checked_on=[8, 10, 12]),
    hyper=dict(form="47/2 - 37/(4 (p-2)) + 5/(4 (p-2)^2)",
               fit_on=[4, 6], checked_on=[8, 10],
               note="two parameters fit with the ladder limit 47/2 imposed"),
    a3_ratio=dict(form="4 p (10 p - 21)", fit_on=[4, 6, 8], checked_on=[10],
                  note="a quadratic through three points; p = 10 is the check"))

GRID = {(2, 0): (4, 6, 8, 10, 12), (2, 2): (4, 6, 8, 10, 12), (2, 4): (4, 6, 8, 10, 12),
        (2, 6): (4, 6, 8, 10), (3, 2): (4, 6, 8, 10), (3, 4): (4, 6, 8, 10)}


def compute_coefficients(cache=None) -> dict:
    """{"p,n,k": {"deficit", "A"}} over GRID, from full_law or a cache of
    its earlier output."""
    out = {}
    for (n, k), ps in GRID.items():
        for p in ps:
            key = f"{p},{n},{k}"
            if cache and key in cache:
                out[key] = dict(deficit=cache[key]["deficit"], A=str(cache[key]["A"]))
            else:
                t0 = time.time()
                law = full_law(p, n, k)
                out[key] = dict(deficit=law["deficit"], A=law["A"],
                                runtime_s=round(time.time() - t0, 1))
                print(f"  p={p} n={n} k={k}: deficit {law['deficit']}  A = {law['A']}", flush=True)
    return out


def tables(coeffs: dict) -> dict:
    """The physical moments and shape ratios from the leading coefficients."""
    A = {key: Fraction(v["A"]) for key, v in coeffs.items()}
    D = {key: v["deficit"] for key, v in coeffs.items()}
    a2, a3 = [], []
    for p in GRID[(2, 0)]:
        k0, k2, k4 = (f"{p},2,{k}" for k in (0, 2, 4))
        assert D[k2] == D[k0] + 1 and D[k4] == D[k0] + 2, p
        w2, w4 = A[k2] / A[k0], A[k4] / A[k0]
        row = dict(p=p, omega2=str(w2), omega4=str(w4), kurtosis=str(w4 / w2 ** 2),
                   kurtosis_float=float(w4 / w2 ** 2),
                   omega2_closed_ok=w2 == omega2_closed(p),
                   omega4_closed_ok=w4 == omega4_closed(p),
                   kurtosis_closed_ok=w4 / w2 ** 2 == kurtosis_closed(p),
                   omega2_over_ladder=float(w2 / ladder_moments(2, p)[2]),
                   omega4_over_ladder=float(w4 / ladder_moments(2, p)[4]))
        k6 = f"{p},2,6"
        if k6 in A:
            assert D[k6] == D[k0] + 3, p
            h = A[k6] * A[k0] ** 2 / A[k2] ** 3
            row.update(hyper=str(h), hyper_float=float(h), hyper_closed_ok=h == hyper_closed(p))
        a2.append(row)
    for p in GRID[(3, 2)]:
        k2, k4 = f"{p},3,2", f"{p},3,4"
        if k2 not in A or k4 not in A:
            continue
        assert D[k4] == D[k2] + 1, p
        r = A[k4] / A[k2]
        a3.append(dict(p=p, ratio=str(r), ratio_closed_ok=r == a3_ratio_closed(p),
                       ratio_over_ladder3=float(r / (40 * p * p)),
                       ratio_over_ladder2=float(r / (28 * p * p))))
    return dict(a2=a2, a3=a3)


def run(out_path=None, cache=None) -> dict:
    t0 = time.time()
    coeffs = compute_coefficients(cache)
    payload = dict(
        provenance="thooft_limit.py -- exact fixed-p, N -> infinity moments "
                   "of the full singlet channels (leading 1/N coefficients of "
                   "full_channel_laws.full_law, each certified out of sample "
                   "in N) against the lambda -> 0 ladders of the N = infinity "
                   "decomposition",
        coefficients=coeffs, laws=LAWS, **tables(coeffs),
        runtime_s=round(time.time() - t0, 1))
    out = Path(out_path) if out_path else HERE / "thooft_limit.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    return payload


def report(payload: dict) -> None:
    print("A_2 at fixed p, N -> infinity (physical units, J = 1):")
    for r in payload["a2"]:
        h = f"  mu6/mu2^3 = {r['hyper']}" if "hyper" in r else ""
        print(f"  p={r['p']:2d}: <w^2> = {r['omega2']}  <w^4> = {r['omega4']}  "
              f"kurtosis = {r['kurtosis']}{h}")
    print("  ladder (Delta = 2): kurtosis 7/2, mu6/mu2^3 47/2")
    print("A_3: <w^4>/<w^2>, ladder 40 p^2 (Delta = 3) against 28 p^2 (Delta = 2):")
    for r in payload["a3"]:
        print(f"  p={r['p']:2d}: {r['ratio']}  ({r['ratio_over_ladder3']:.3f} x 40p^2)")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--from-cache", type=str, default=None)
    args = ap.parse_args()
    if args.run or args.from_cache:
        cache = None
        if args.from_cache:
            cache = {}
            for path in args.from_cache.split(","):
                with open(path) as fh:
                    cache.update(json.load(fh))
        report(run(cache=cache))
    else:
        with open(HERE / "thooft_limit.json") as fh:
            report(json.load(fh))
