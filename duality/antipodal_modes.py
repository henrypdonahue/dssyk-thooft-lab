#!/usr/bin/env python3
r"""
Mode separation of the antipodal growth: a mechanism, tested.
=============================================================

The deep campaign left two open items.  The flipped mode's rate, and a
mechanism for the transient window.

The natural candidate had a sharp signature.  In the semiclassical
skeleton the growing tower e^{n v thL} rotates by e^{-i n pi} under the
half-fake-period shift.  Odd harmonics flip, even ones do not.  So the
window would end when the unflippable n = 2 harmonic, at rate exactly
2 kappa_1, overtakes.

Two instruments, one verdict.

1. SHIFT-CIRCLE HARMONIC ANALYSIS (primary, and it needs no exponential
   fitting in thL).  Sample the curve at shifts Delta along the PHYSICAL
   half of the fake circle, Delta in [0, pi/v] -- the range this module's
   companion (signflip_probe.g2_along_shift) certifies, short of the
   2-pt pole at Delta = pi/2 + pi/v -- and fit c_0 + sum_n A_n e^{-i n v
   Delta} by least squares.  If the tower were harmonic, |A_n(thL)|
   would grow at rates kappa_n = n kappa_1 with |A_2| << |A_1|.

   MEASURED (two- and three-harmonic models both recorded, because the
   answer depends on the model at the 10% level and the residuals are
   5-9% of the range -- the response is not a clean sum of a few
   v-harmonics): the second harmonic is subdominant, |A_2/A_1| ~
   0.3-0.6, and grows faster than the first, kappa_2/kappa_1 ~ 1.2-1.5
   at lambda = 1 and 0.6.  That sits between the harmonic tower's 2 and
   the ratio ~1 an earlier version of this instrument reported.

   That earlier instrument sampled the FULL fake circle and took a DFT.
   Its samples past the 2-pt pole (Delta > pi/2 + pi/v) are two orders
   of magnitude larger than the physical ones, so every DFT bin was the
   transform of that spike: the "|c_2/c_1| ~ 1, kappa_2/kappa_1 ~ 1"
   refutation was an artifact and is withdrawn.  The corrected verdict
   is weaker in both directions: the even-harmonic mechanism is neither
   refuted nor established.  A subdominant, faster-growing second
   harmonic is there; it is not at rate 2 kappa_1 and it does not carry
   the weight a clean tower would need.

2. TWO-MODE VARIABLE-PROJECTION FIT (secondary).  It establishes the
   SIGNS robustly: A1 > 0, the flipped mode, with a smaller A2 < 0,
   unflipped and growing, wherever the flip operates.

   But the rate surface is degenerate.  The free fit wanders, and
   lambda = 0.6 collapses outright, which is kept here as the recorded
   failure.  Rate extraction from these curves is unreliable, and the
   band's continuous content is why.

NET.  The anti-scrambling epoch ends because an unflipped growing
component overtakes.  But that component is NOT the semiclassical second
harmonic.  The flipped mode's dS-rate question stays open, and it needs a
method that resolves the band's continuous content -- working directly
with Lehmann data, for example.

Output: antipodal_modes.json.  Asserted in test_antipodal_modes.py.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

from chord import beta_chord_from_betaJ
from chord_fold import CONFIG_CROSSED, ChordCorrelators
from largeq_anchor import v_of_betaJ

HERE = Path(__file__).resolve().parent


def _curve(cc, Nflav, delta, thLs, config=CONFIG_CROSSED):
    t1, t2, t3, t4 = config
    return np.array([cc.connected_ratio(
        (t1 + delta + 1j * t, t2 + delta + 1j * t, t3, t4),
        'c', Nflav).real for t in thLs])


def two_mode_fit(thLs: np.ndarray, r: np.ndarray,
                 kappa_grid=np.linspace(0.15, 1.2, 106),
                 constrained: bool = True):
    """Variable projection: for each kappa, solve linearly for
    (c0, A1, A2) in r = c0 + A1 e^{k thL} + A2 e^{k2 thL}; return the
    best kappa and amplitudes.  constrained: k2 = 2k; free: scan k2 too.

    kappa_at_grid_max flags fits pinned at the scan boundary (censored,
    not optima -- several committed free fits and the lambda = 0.8
    constrained fit saturate the 1.2 ceiling; kappa_over_v/t_star
    derived from a flagged fit inherit the censoring)."""
    best = None
    k2_grid = np.linspace(0.3, 2.4, 71)
    for k in kappa_grid:
        for k2 in ([2.0 * k] if constrained else k2_grid):
            if k2 <= k + 0.05:
                continue
            A = np.vstack([np.ones_like(thLs), np.exp(k * thLs),
                           np.exp(k2 * thLs)]).T
            coef, res, *_ = np.linalg.lstsq(A, r, rcond=None)
            resid = float(np.sqrt(np.mean((A @ coef - r) ** 2)))
            if best is None or resid < best["resid"]:
                best = dict(kappa=float(k), kappa2=float(k2),
                            c0=float(coef[0]), A1=float(coef[1]),
                            A2=float(coef[2]), resid=resid)
    best["kappa_at_grid_max"] = bool(
        best["kappa"] >= float(kappa_grid[-1]) - 1e-12
        or (not constrained
            and best["kappa2"] >= float(k2_grid[-1]) - 1e-12))
    return best


def separate(lam: float, betaJ: float = 2.0, p: int = 4,
             nmax: int = 90) -> dict:
    """Two-mode separation at one lambda: antipodal and unshifted."""
    v = v_of_betaJ(betaJ)
    Nflav = 2.0 * p * p / lam
    cc = ChordCorrelators(lam, p, beta_chord_from_betaJ(betaJ, lam), nmax)
    thLs = np.linspace(0.2, 2.4, 12)
    out = dict(lam=lam, betaJ=betaJ, v=v)
    for name, delta in [("antipodal", math.pi / v), ("unshifted", 0.0)]:
        r = _curve(cc, Nflav, delta, thLs)
        con = two_mode_fit(thLs, r, constrained=True)
        free = two_mode_fit(thLs, r, constrained=False)
        rng = float(np.max(np.abs(r - np.mean(r))))
        row = dict(curve=r.tolist(), thLs=thLs.tolist(),
                   constrained=con, free=free,
                   resid_over_range_con=con["resid"] / rng,
                   resid_over_range_free=free["resid"] / rng)
        if name == "antipodal" and con["A1"] > 0 > con["A2"]:
            row["t_star"] = math.log(con["A1"] / abs(con["A2"])) / con["kappa"]
        out[name] = row
    out["kappa_over_v"] = out["antipodal"]["constrained"]["kappa"] / v
    out["kappa_over_veff"] = (out["antipodal"]["constrained"]["kappa"]
                              / (v / (1.0 + v)))
    return out


def harmonics(lam: float, betaJ: float = 2.0, p: int = 4, nmax: int = 90,
              n_delta: int = 9, thLs=(0.8, 1.2, 1.6, 2.0),
              n_harm=(2, 3)) -> dict:
    """The clean experiment: sample the curve at n_delta shifts along the
    physical half-circle Delta in [0, pi/v] and least-squares fit
    c_0 + sum_{n=1..n_harm} A_n e^{-i n v Delta} at each thL.  Records, for
    each model size, |A_n(thL)|, the per-harmonic growth rates kappa_n
    (mean log-slope of |A_n| in thL), the ratio kappa_2/kappa_1, the
    weight ratios |A_2/A_1|, and the fit residual over the curve's range.
    Samples beyond the 2-pt pole at Delta = pi/2 + pi/v are excluded by
    construction (they dominated the earlier full-circle DFT)."""
    v = v_of_betaJ(betaJ)
    Nflav = 2.0 * p * p / lam
    cc = ChordCorrelators(lam, p, beta_chord_from_betaJ(betaJ, lam), nmax)
    thLs = np.array(thLs)
    deltas = np.linspace(0.0, math.pi / v, n_delta)
    R = np.empty((len(thLs), n_delta))
    for j, d in enumerate(deltas):
        R[:, j] = _curve(cc, Nflav, float(d), thLs)
    out = dict(lam=lam, betaJ=betaJ, v=v, n_delta=n_delta,
               deltas=deltas.tolist(), thLs=thLs.tolist(),
               pole_delta=math.pi / 2.0 + math.pi / v, R=R.tolist(),
               models={})
    for nh in n_harm:
        A = np.column_stack([np.exp(-1j * n * v * deltas) for n in range(nh + 1)])
        amps = np.empty((len(thLs), nh + 1))
        resid = []
        for i in range(len(thLs)):
            c = np.linalg.lstsq(A, R[i].astype(complex), rcond=None)[0]
            amps[i] = np.abs(c)
            pred = (A @ c).real
            rng = R[i].max() - R[i].min()
            resid.append(float(np.sqrt(np.mean((pred - R[i]) ** 2)) / rng))
        m = dict(resid_over_range=resid)
        for n in range(nh + 1):
            m[f"abs_A{n}"] = amps[:, n].tolist()
        for n in range(1, nh + 1):
            rates = np.diff(np.log(np.maximum(amps[:, n], 1e-300))) / np.diff(thLs)
            m[f"kappa_{n}"] = float(np.mean(rates))
        m["A2_over_A1"] = (amps[:, 2] / np.maximum(amps[:, 1], 1e-300)).tolist()
        m["kappa2_over_kappa1"] = (m["kappa_2"] / m["kappa_1"]
                                   if m["kappa_1"] > 0 else float("nan"))
        out["models"][str(nh)] = m
    # headline numbers from the two-harmonic model, with the three-harmonic
    # model's values as the model-dependence band
    m2, m3 = out["models"]["2"], out["models"]["3"]
    out["A2_over_A1_range"] = [min(min(m2["A2_over_A1"]), min(m3["A2_over_A1"])),
                               max(max(m2["A2_over_A1"]), max(m3["A2_over_A1"]))]
    out["kappa2_over_kappa1_range"] = sorted([m2["kappa2_over_kappa1"],
                                              m3["kappa2_over_kappa1"]])
    return out


def main():
    t0 = time.time()
    payload = dict(
        provenance="antipodal_modes.py -- mode separation of the "
                   "antipodal growth: half-circle harmonic analysis on "
                   "Delta in [0, pi/v] (primary; two- and three-harmonic "
                   "least-squares models, no exponential fitting in thL) "
                   "plus a two-mode variable-projection fit (secondary); "
                   "odd harmonics flip, even do not (DISCRIMINATOR.md)",
        harmonics=[harmonics(lam) for lam in (1.5, 1.0, 0.6)],
        points=[separate(lam) for lam in (1.5, 1.2, 1.0, 0.8, 0.6)])
    payload["runtime_s"] = round(time.time() - t0, 1)
    out = HERE / "antipodal_modes.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    for h in payload["harmonics"]:
        for nh, m in h["models"].items():
            print(f"harmonics lam={h['lam']} ({nh} harmonics): "
                  f"kappa_1={m['kappa_1']:.3f} kappa_2={m['kappa_2']:.3f} "
                  f"(ratio {m['kappa2_over_kappa1']:.2f})  |A2/A1|: "
                  f"{m['A2_over_A1'][0]:.3f} -> {m['A2_over_A1'][-1]:.3f}  "
                  f"resid/range {max(m['resid_over_range']):.3f}  v={h['v']:.3f}")
    for pt in payload["points"]:
        a = pt["antipodal"]["constrained"]
        f = pt["antipodal"]["free"]
        print(f"lam={pt['lam']}: kappa={a['kappa']:.3f} "
              f"(free k2/k={f['kappa2']/f['kappa']:.2f})  "
              f"A1={a['A1']:+.3f} A2={a['A2']:+.3f}  "
              f"resid con/free={a['resid']/max(f['resid'],1e-12):.2f}  "
              f"kappa/v={pt['kappa_over_v']:.2f}  "
              f"t*={pt['antipodal'].get('t_star', float('nan')):.2f}")
    print(f"wrote {out} ({payload['runtime_s']} s)")


if __name__ == "__main__":
    main()
