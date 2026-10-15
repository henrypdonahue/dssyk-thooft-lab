#!/usr/bin/env python3
r"""
The double-scaling constants of the pair channel have closed forms.
===================================================================

road_map.md left one number on the tower-extraction item as "an
estimate, not a law": the double-scaling constant of the pair channel,
c_2(p) = |mu_2^pair(N = p^2, p; n = 2)| p N along lambda_chord = 2, exact
to p = 48 (31.654) and still rising, with a windowed extrapolation
c_2(inf) = 32.334(1) and no closed form identified.  Its n = 3 partner
c_3(p) was known to p = 16 (209.76), still rising, with no target.

THE LAW.  On the line lambda_chord = lambda (N = 2p^2/lambda,
q = e^{-lambda}) the pair channel's constant equals the Delta_psi -> 0
slope of the N = infinity DIAGONAL channel (duality/diagonal_chord.py):

    c_{n,k}(lambda) := lim_{p->inf} |mu_k^pair| p N = lambda P_K(q),
    mu_k^diag(q, Delta) = lambda P_K(q) Delta + O(Delta^2),
    K = k + 2n,

with P_K depending on the word length K only (the ad-shift structure of
the ladder, exact at N = infinity), and, from the symbolic walker
(sympy on diag_word_amplitude, exact):

    P_2  = 2
    P_4  = 4 (1-q)
    P_6  = 4 (1-q)^2 (5 + 3q)
    P_8  = 4 (1-q)^3 (35 + 56q + 35q^2 + 10q^3)
    P_10 = 4 (1-q)^4 (294 + 840q + 1134q^2 + 951q^3 + 525q^4 + 189q^5 + 35q^6)

So c_2(lambda) = 4 lambda (1-q)^2 (5+3q).  At lambda = 2 that is
32.3341891, against the artifact's 32.334(1), and against a global
polynomial-in-1/p fit of the 18 exact points (orders 4-6: 32.33419,
stable to 1e-5).  It is NOT 97/3.

Equivalently, on every double-scaling line

    (N-1) mu_k^pair + mu_k^diag = O(Delta^2):

the O(Delta) parts of the pair and diagonal channels cancel, which is
the n >= 2 shadow of the exact A_1 identity (N-1) mu^pair + mu^diag = 0.
That cancellation is the mechanism; it is measured here at finite N
(mu_2^diag p -> P_6 from above, (N-1) mu_2^pair p -> -P_6 from below, the
sum x p^2 converging to a finite O(Delta^2) coefficient), and it is what
a chord-side derivation of the pair channel's 1/N term has to produce.
The full-channel width is the ratio of the two O(Delta^2) remainders, so
it is not fixed by the P_K alone: on lambda_chord = 2 the exact series
to p = 40 gives w2_full(inf) = 1.31691(1) (diagonal_chord.py).

EVIDENCE, all exact-arithmetic sequences with polynomial-in-1/p limits:

  * K = 6: lambda = 2 (p <= 48, committed): 32.33419 vs 32.33419;
    lambda = 1 (p <= 32): 9.75548 vs 9.75548;  lambda = 4 (p <= 40):
    77.9436 vs 77.9436.  Three lines, one function.
  * K = 8: lambda = 2 (p <= 32): the c_3(p) sequence heads to
    lambda P_8 = 223.66; lambda = 1, 4 likewise (see the artifact for
    the limit estimates and their window spreads).
  * The small-lambda expansion of the law fixes the large-p growth of
    every exact rational 1/N coefficient of the fixed-p laws
    (pair_moment_laws.py): mu_2 = -(lambda^2/ (2 p^3)) ... more
    precisely mu_2^pair(N = 2p^2/lambda) = -2 lambda^2 (1-q)^2 (5+3q)/p^3
    = -(16 lambda^4 - 22 lambda^5 + 55/3 lambda^6 - ...)/p^3, so the
    N^-(4+m) coefficient C_m(p) grows as -gamma_m 2^(4+m) p^(5+2m):
    A -> -256 p^5 (the certified A(2,2) = -64 p^3 (p-1)(4p-5) does), and
    B/p^7 -> +704.  The committed B(2,2; p) at p = 4..12 fits

        B(2,2; p) = 64 p^3 (p-1) (11 p^3 - 36 p^2 + 29 p - 3),

    with integer coefficients from a five-point fit, leading 704 p^7 as
    predicted, and it reproduces the committed p = 2 value (-512) out of
    sample.  A new certified closed form, and the first subleading one.

Scope, stated plainly.  The identification c_{n,k} = lambda P_K is a
measured law: exact sequences, three lines, six digits at K = 6, and
the finite-N cancellation seen directly.  It is not a derivation.  The
derivation is the open theory item that road_map.md keeps: the 1/N
correction to the pair chord rules, which must reproduce P_K.

Output: double_scaling_law.json (the P_K, the sequences on three lines
for K = 6 and 8, their limits, the cancellation table, the B closed
form), pinned by test_double_scaling_law.py.  The producer is
checkpointed and resumable (--resume); the K = 8 points at large p cost
minutes each (7 workers: about an hour for the default grid).
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

from pair_moment_laws import _mu_sample

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "duality"))


# ---------------------------------------------------------------------------
# the closed forms
# ---------------------------------------------------------------------------
# P_K(q) = 4 (1-q)^(K/2-1) R_K(q) for K >= 4; coefficient lists of R_K,
# index = power of q (from the symbolic walker; recomputed for K <= 8 by
# test_double_scaling_law.py, checked numerically at K = 10)
_R = {
    4: [1],
    6: [5, 3],
    8: [35, 56, 35, 10],
    10: [294, 840, 1134, 951, 525, 189, 35],
}


def P_K(K: int, q: float):
    """The Delta-slope polynomial of the N = infinity diagonal channel:
    mu_k^diag -> lambda P_K(q) Delta, K = k + 2n."""
    if K == 2:
        return 2 * (q ** 0)
    R = sum(c * q ** j for j, c in enumerate(_R[K]))
    return 4 * (1 - q) ** (K // 2 - 1) * R


def c_closed_form(lam: float, K: int) -> float:
    """c_{n,k}(lambda) = lambda P_K(e^-lambda)."""
    return lam * P_K(K, math.exp(-lam))


def delta_slope_symbolic(n: int, k: int):
    """The exact Delta -> 0 slope of mu_k^diag / lambda as a sympy
    polynomial in q, from the parity-law walker run symbolically.  Seconds
    for K <= 8; minutes at K = 10."""
    import sympy as sp
    import diagonal_chord as dc
    q, d = sp.symbols("q d", positive=True)
    ks = tuple(sorted(set((0, 2, k))))
    mu = sp.expand(dc.channel_moments_diag(q, d, n, ks=ks)[f"mu{k}"])
    s = sp.diff(mu, d).subs(d, 0).subs(sp.log(q), -sp.Symbol("lam"))
    return sp.nsimplify(sp.expand(s / sp.Symbol("lam")), rational=True)


def B22_closed_form(p: int) -> int:
    """The N^-5 coefficient of mu_2^pair(n = 2) at fixed p."""
    return 64 * p ** 3 * (p - 1) * (11 * p ** 3 - 36 * p ** 2 + 29 * p - 3)


def lambda_expansion(order: int = 8):
    """Taylor coefficients gamma_m of 2 lambda^2 (1-q)^2 (5+3q) = sum
    gamma_j lambda^j, and the predicted large-p growth of the N^-(4+m)
    coefficients: C_m(p)/p^(5+2m) -> -gamma_{4+m} 2^(4+m)."""
    import sympy as sp
    lam = sp.symbols("lam")
    F = 2 * lam ** 2 * (1 - sp.exp(-lam)) ** 2 * (5 + 3 * sp.exp(-lam))
    ser = sp.series(F, lam, 0, order + 1).removeO()
    out = []
    for m in range(0, order - 3):
        g = ser.coeff(lam, 4 + m)
        out.append(dict(m=m, gamma=str(g),
                        C_m_over_p_power=str(-g * 2 ** (4 + m)),
                        power=5 + 2 * m))
    return out


# ---------------------------------------------------------------------------
# the exact sequences on the double-scaling lines
# ---------------------------------------------------------------------------
def line_N(lam: float, p: int) -> int:
    N = 2 * p * p / lam
    if abs(N - round(N)) > 1e-9:
        raise ValueError(f"2p^2/lambda not an integer at p={p}, lambda={lam}")
    return int(round(N))


def limit_estimate(rows: list, min_points: int = 5) -> dict:
    """Polynomial-in-1/p least-squares limits of c(p), orders 2..5 on the
    deepest window, central = highest order that is residual-stable;
    spread over orders and over the two deepest windows.  Estimate."""
    import numpy as np
    pts = sorted((r["p"], r["c"]) for r in rows)
    if len(pts) < min_points:
        return dict(estimated=False, n_points=len(pts))
    fits = []
    for lo in (0, 1, 2):
        sub = pts[lo:]
        if len(sub) < min_points:
            continue
        x = np.array([1.0 / p for p, _ in sub])
        y = np.array([c for _, c in sub])
        for order in (2, 3, 4, 5):
            if len(sub) <= order + 1:
                continue
            A = np.vander(x, order + 1, increasing=True)
            coef = np.linalg.lstsq(A, y, rcond=None)[0]
            res = float(np.max(np.abs(A @ coef - y)))
            fits.append(dict(order=order, p_min=sub[0][0], c_inf=float(coef[0]),
                             max_residual=res))
    best = min(fits, key=lambda f: f["max_residual"])
    highs = [f for f in fits if f["order"] >= 3]
    spread = max(f["c_inf"] for f in highs) - min(f["c_inf"] for f in highs)
    return dict(estimated=True, n_points=len(pts), fits=fits,
                c_inf_central=best["c_inf"], spread=float(spread),
                note="least-squares polynomial in 1/p, orders 2-5, three "
                     "windows; central = smallest-residual fit; estimate")


DEFAULT_GRID = {
    # (lambda, n): p values.  lambda = 2, n = 2 is imported from the ladder.
    (1.0, 2): [6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32],
    (4.0, 2): [6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40],
    (1.0, 3): [6, 8, 10, 12, 14, 16, 18, 20],
    (2.0, 3): [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32],
    (4.0, 3): [6, 8, 10, 12, 14, 16, 18, 20, 24],
}


def _write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)


def _ladder_n2_rows() -> list:
    """The committed lambda_chord = 2, n = 2 sequence (pair_moment_ladder.json
    + pair_moment_laws.json), re-keyed to this artifact's row format."""
    rows = {}
    for name in ("pair_moment_laws.json", "pair_moment_ladder.json"):
        path = HERE / name
        if not path.exists():
            continue
        with open(path) as fh:
            d = json.load(fh)
        for r in d.get("double_scaling", {}).get("n2", []):
            rows[r["p"]] = dict(lam=2.0, n=2, p=r["p"], N=r["N"],
                                mu2=r["mu2"], c=r["c"], source=name)
    return [rows[p] for p in sorted(rows)]


def cancellation_table(ps=(4, 6, 8, 10, 12, 14, 16, 20, 24)) -> list:
    """(N-1) mu_2^pair + mu_2^diag and the mu_0 analog on lambda_chord = 2,
    exact at finite N, next to the closed forms they approach."""
    import exact_pair_moments as epm
    out = []
    for p in ps:
        N = p * p
        pr = epm.channel_pair_moments(N, p, 2, marks=('i', 'i', 'j', 'j'),
                                      ks=(0, 2), exact=True)
        dg = epm.channel_pair_moments(N, p, 2, marks=('i', 'i', 'i', 'i'),
                                      ks=(0, 2), exact=True)
        epm._DP_CACHE.clear()
        epm._EXP_CACHE.clear()
        row = dict(p=p, N=N)
        for k in (0, 2):
            s = (N - 1) * pr[f"mu{k}"] + dg[f"mu{k}"]
            row[f"diag_mu{k}_times_p"] = float(dg[f"mu{k}"] * p)
            row[f"pair_mu{k}_Nm1_times_p"] = float((N - 1) * pr[f"mu{k}"] * p)
            row[f"sum_mu{k}_times_p2"] = float(s * p * p)
        out.append(row)
    return out


def run(workers: int = None, resume: bool = False, out_path=None,
        grid: dict = None) -> dict:
    t0 = time.time()
    out = Path(out_path) if out_path else HERE / "double_scaling_law.json"
    grid = grid or DEFAULT_GRID
    payload = None
    if resume and out.exists():
        with open(out) as fh:
            payload = json.load(fh)
    if payload is None:
        payload = dict(sequences={}, runtime_s=0.0)
    seqs = payload.setdefault("sequences", {})
    have = {(float(r["lam"]), int(r["n"]), int(r["p"]))
            for key, rows in seqs.items() for r in rows}
    tasks = []
    for (lam, n), ps in grid.items():
        for p in ps:
            if (lam, n, p) not in have:
                tasks.append((line_N(lam, p), p, n, 2, lam))
    workers = workers or max(1, mp.cpu_count() - 1)
    print(f"{len(tasks)} exact points to compute on {workers} workers",
          flush=True)
    if tasks:
        ctx = mp.get_context("spawn")
        # heaviest first so the pool tail is short
        order = sorted(tasks, key=lambda t: -(t[1] ** 2 * t[2]))
        with ctx.Pool(processes=min(workers, len(order))) as pool:
            for (N, p, n, k), mu in pool.imap_unordered(
                    _mu_sample, [t[:4] for t in order]):
                lam = 2 * p * p / N
                key = f"lam{lam:g}_n{n}"
                seqs.setdefault(key, []).append(dict(
                    lam=lam, n=n, p=p, N=N, mu2=float(mu), mu2_exact=str(mu),
                    c=float(abs(mu) * p * N)))
                seqs[key].sort(key=lambda r: r["p"])
                print(f"  lam={lam:g} n={n} p={p:2d} N={N:5d} "
                      f"c={float(abs(mu) * p * N):.6f}", flush=True)
                payload["runtime_s"] = round(
                    payload.get("runtime_s", 0.0) + 0, 1)
                _write_json(out, payload)
    seqs["lam2_n2"] = _ladder_n2_rows()
    # limits and closed forms
    limits = {}
    for key, rows in seqs.items():
        lam = rows[0]["lam"]
        n = rows[0]["n"]
        K = 2 + 2 * n
        est = limit_estimate(rows)
        est.update(lam=lam, n=n, K=K, closed_form=c_closed_form(lam, K),
                   closed_form_expr=f"lambda P_{K}(e^-lambda)")
        if est.get("estimated"):
            est["rel_dev"] = abs(est["c_inf_central"] - est["closed_form"]) \
                / est["closed_form"]
        limits[key] = est
    payload["limits"] = limits
    payload["P_K"] = {str(K): dict(R_coeffs=_R[K], prefactor=f"4 (1-q)^{K // 2 - 1}")
                      for K in _R}
    payload["P_2"] = 2
    payload["lambda_expansion"] = lambda_expansion()
    payload["B22_closed_form"] = "64 p^3 (p-1)(11 p^3 - 36 p^2 + 29 p - 3)"
    payload["B22_check"] = _b22_check()
    print("cancellation table", flush=True)
    payload["cancellation_lambda2"] = cancellation_table()
    payload["provenance"] = (
        "double_scaling_law.py -- the pair channel's double-scaling "
        "constants c_{n,k}(lambda) = lim |mu_k^pair| p N on lambda_chord = "
        "lambda, exact sequences (pair_moment_laws.mu_exact) on three lines, "
        "against the closed forms lambda P_K(e^-lambda) from the Delta-slope "
        "of the N = infinity diagonal chord channel; the finite-N "
        "diagonal-pair O(Delta) cancellation; the B(2,2;p) closed form")
    payload["runtime_s"] = round(payload.get("runtime_s", 0.0)
                                 + time.time() - t0, 1)
    _write_json(out, payload)
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


def _b22_check() -> dict:
    """The committed B(2,2; p) values against the closed form."""
    path = HERE / "pair_moment_laws.json"
    with open(path) as fh:
        d = json.load(fh)
    rows = [(l["p"], Fraction(l["B"])) for l in d["laws"]
            if l["n"] == 2 and l["k"] == 2]
    return dict(fit_ps=[p for p, _ in rows if p >= 4],
                check_p=[p for p, _ in rows if p < 4],
                matches={str(p): (B == B22_closed_form(p)) for p, B in rows},
                all_match=all(B == B22_closed_form(p) for p, B in rows))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()
    if args.run or args.resume:
        run(workers=args.workers, resume=args.resume)
    else:
        for lam in (1.0, 2.0, 4.0):
            print(f"lambda={lam}: c_2 -> {c_closed_form(lam, 6):.6f}   "
                  f"c_3 -> {c_closed_form(lam, 8):.6f}")
        print("B(2,2;p) closed form matches committed laws:", _b22_check())
