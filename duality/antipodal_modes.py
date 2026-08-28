#!/usr/bin/env python3
r"""
Mode separation of the antipodal growth: a mechanism, tested and refuted.
=========================================================================

The deep campaign left two open items.  The flipped mode's rate, and a
mechanism for the transient window.

The natural candidate had a sharp signature.  In the semiclassical
skeleton the growing tower e^{n v thL} rotates by e^{-i n pi} under the
half-fake-period shift.  Odd harmonics flip, even ones do not.  So the
window would end when the unflippable n = 2 harmonic, at rate exactly
2 kappa_1, overtakes.

Two instruments, one verdict.

1. SHIFT-CIRCLE FOURIER ANALYSIS (primary, and it needs no exponential
   fitting).  Sample the curve at M shifts around the full fake circle
   and take the DFT in Delta.  If the tower were harmonic, |c_n(thL)|
   would grow at rates kappa_n = n kappa_1.

   MEASURED: kappa_2/kappa_1 is about 0.9-1.2 on the rows with a
   measurable growth window (kappa_1 > 0.05), with comparable weights
   |c_2/c_1| about 0.7-1.0.  The lam = 1.5 row has both rates negative,
   so the modes decay and the ratio is undefined as a growth comparison;
   that gate excludes it in the tests.

   So the harmonic mechanism is REFUTED at finite lambda, where it is
   testable.  The reason is structural.  The exact shift-response has
   CONTINUOUS frequency content, set by the energy band, because each
   Lehmann line contributes its own e^{-dtau dE}.  It does not have
   v-harmonics.  The harmonic tower is only the strict lambda -> 0
   skeleton.  This matches the earlier finding that the phase period (v)
   and the magnitude growth (kappa > v) decouple at finite lambda.

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
              M: int = 8, thLs=(0.8, 1.2, 1.6, 2.0)) -> dict:
    """The clean experiment: sample the curve at M shifts around the FULL
    fake circle, Delta_j = j (2 pi/v)/M, and DFT in Delta.  Each shift
    harmonic n contributes A_n e^{...thL} e^{-i n v Delta}: the DFT
    separates them with no exponential fitting.  Tests: |c_1| flips sign
    under the half-period (built in), |c_2| is shift-invariant; the
    growth rate of each |c_n| in thL, kappa_n, should satisfy
    kappa_n ~ n * kappa_1 if the tower is harmonic."""
    v = v_of_betaJ(betaJ)
    Nflav = 2.0 * p * p / lam
    cc = ChordCorrelators(lam, p, beta_chord_from_betaJ(betaJ, lam), nmax)
    thLs = np.array(thLs)
    deltas = np.arange(M) * (2.0 * math.pi / v) / M
    R = np.empty((len(thLs), M))
    for j, d in enumerate(deltas):
        R[:, j] = _curve(cc, Nflav, float(d), thLs)
    # Under the model R_j = sum_n A_n e^{-i n v Delta_j}, fft bin n holds
    # A_{-n mod M}, not A_n; _curve takes .real, so |A_{-n}| = |A_n| and
    # every magnitude diagnostic below is unaffected.  Dropping .real (a
    # complex-ratio extension) would need np.fft.ifft * M instead.
    C = np.fft.fft(R, axis=1) / M          # c_n(thL), n = 0..M-1
    out = dict(lam=lam, betaJ=betaJ, v=v, M=M, thLs=thLs.tolist())
    for n in (0, 1, 2, 3):
        cn = C[:, n]
        out[f"abs_c{n}"] = np.abs(cn).tolist()
        if n >= 1:
            # per-harmonic growth rate from consecutive |c_n| ratios
            rates = np.diff(np.log(np.maximum(np.abs(cn), 1e-300))) \
                / np.diff(thLs)
            out[f"kappa_{n}"] = float(np.mean(rates))
    # phase test: c_1 at the half-period equals -c_1 at 0 by construction;
    # the physical content is the RELATIVE harmonic weights vs thL
    out["c2_over_c1_last"] = float(np.abs(C[-1, 2]) / max(np.abs(C[-1, 1]),
                                                          1e-300))
    out["c2_over_c1_first"] = float(np.abs(C[0, 2]) / max(np.abs(C[0, 1]),
                                                          1e-300))
    return out


def main():
    t0 = time.time()
    payload = dict(
        provenance="antipodal_modes.py -- mode separation of the "
                   "antipodal growth: shift-circle Fourier analysis "
                   "(primary; no exponential fitting) plus a two-mode "
                   "variable-projection fit (secondary); odd harmonics "
                   "flip, even do not (DISCRIMINATOR.md)",
        harmonics=[harmonics(lam) for lam in (1.5, 1.0, 0.6)],
        points=[separate(lam) for lam in (1.5, 1.2, 1.0, 0.8, 0.6)])
    payload["runtime_s"] = round(time.time() - t0, 1)
    out = HERE / "antipodal_modes.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    for h in payload["harmonics"]:
        print(f"harmonics lam={h['lam']}: kappa_1={h['kappa_1']:.3f} "
              f"kappa_2={h['kappa_2']:.3f} (ratio "
              f"{h['kappa_2']/max(h['kappa_1'],1e-9):.2f})  "
              f"|c2/c1|: {h['c2_over_c1_first']:.3f} -> "
              f"{h['c2_over_c1_last']:.3f}  v={h['v']:.3f}")
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
