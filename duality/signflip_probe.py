#!/usr/bin/env python3
r"""
Can a contour shift flip the scrambling sign?  An exploratory probe.
====================================================================

The dS-dictionary problem, measured in this repo: under the flat time map
DSSYK scrambles (Re a = -1 at every temperature), a dS observer needs the
opposite sign, and the Euclidean fold does not deliver it.  This module
tests the next candidate: displace the probe pair by HALF THE FAKE PERIOD
around the thermal circle (Delta = pi/v in Euclidean angle; the fake
period is 2 pi/v).  Geometrically this is an antipodal shift on the fake
disk.  It is also what the Narovlansky-Verlinde dressing does: their
operators sit at +-i beta_dS/4, i.e. half a period apart, and they call
the two insertions antipodal.

Three levels, one JSON:

1. Closed form (large-q, `antiscrambling.otoc_ratio`).  Shift the
   (theta1, theta2) pair by Delta and extract the growth coefficient
   a(Delta) by a complex fit in the late-time window.  Prediction from
   the pole structure: a(Delta) = a(0) exp(-i v Delta), so at
   Delta = pi/v the coefficient flips exactly: Re a goes from -1 to +1.
   The probe verifies the rotation law and the flip.

2. Two-point physicality.  The shifted W-V separation is
   theta = pi/2 + Delta.  The continued large-q 2-pt stays real and
   positive up to its pole at theta = pi + pi/v, so the antipodal point
   sits INSIDE the physical window, at margin pi/2 from the pole.
   The probe records the values and the margin.

3. Exact chord theory at finite lambda (`chord_fold.ChordCorrelators`).
   Same shifted contours, exact in lambda, truncation-checked.  The
   observable: does Re F/F_d move DOWN with Lorentzian time (scrambling)
   at Delta = 0 and UP (anti-scrambling side) at Delta = pi/v?  Plus the
   fitted coefficient where the window allows.

Honesty: a sign flip under a contour shift is NECESSARY for this
candidate map, not sufficient.  It does not by itself build the dS
dictionary: the full correlator suite must stay consistent, and the
algebra story (why the observer's operators live there) is NV's
construction to complete.  This probe measures whether the mechanism
survives exactness in lambda -- the question the closed form cannot
answer.

Output: signflip_probe.json; asserted in test_signflip.py; narrated in
DISCRIMINATOR.md.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

from antiscrambling import G2_thermal, otoc_ratio
from chord_fold import CONFIG_CROSSED, ChordCorrelators, thetas_to_z
from largeq_anchor import v_of_betaJ

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# level 1: closed form
# ---------------------------------------------------------------------------
def growth_coefficient_shifted(v: float, delta: float) -> complex:
    """Coefficient of exp(v thL) in the closed-form crossed ratio, with the
    (theta1, theta2) pair shifted by `delta` in Euclidean angle.  Complex
    least-squares fit on the basis {1, e^{v thL}, e^{-v thL}} over a late
    window (the same technique antiscrambling.py validates at delta = 0)."""
    t1, t2, t3, t4 = CONFIG_CROSSED
    thLs = np.array([4.0, 5.0, 6.0, 7.0]) / v
    vals = np.array([complex(otoc_ratio(t1 + delta + 1j * t, t2 + delta + 1j * t,
                                        t3, t4, v)) for t in thLs])
    A = np.vstack([np.ones_like(thLs), np.exp(v * thLs),
                   np.exp(-v * thLs)]).T
    coef, *_ = np.linalg.lstsq(A.astype(complex), vals, rcond=None)
    return complex(coef[1])


def closed_form_scan(betaJ: float, n_delta: int = 9) -> dict:
    """a(Delta) over Delta in [0, pi/v]: the rotation law and the flip."""
    v = v_of_betaJ(betaJ)
    a0 = growth_coefficient_shifted(v, 0.0)
    deltas = np.linspace(0.0, math.pi / v, n_delta)
    rows = []
    for d in deltas:
        a = growth_coefficient_shifted(v, float(d))
        pred = a0 * np.exp(-1j * v * d)
        rows.append(dict(delta=float(d), re_a=a.real, im_a=a.imag,
                         rotation_law_error=abs(a - pred) / abs(a0)))
    return dict(betaJ=betaJ, v=v, re_a_unshifted=a0.real,
                re_a_antipodal=rows[-1]["re_a"], rows=rows)


# ---------------------------------------------------------------------------
# level 2: 2-pt physicality along the shift
# ---------------------------------------------------------------------------
def g2_along_shift(betaJ: float, q_syk: int = 6, n_pts: int = 12) -> dict:
    """The continued large-q 2-pt at the shifted W-V separation
    pi/2 + Delta, up to the antipodal point.  Pole sits at pi + pi/v:
    margin pi/2 at the endpoint."""
    v = v_of_betaJ(betaJ)
    deltas = np.linspace(0.0, math.pi / v, n_pts)
    vals = [float(np.real(G2_thermal(0.5 * math.pi + d, betaJ, q_syk)))
            for d in deltas]
    return dict(betaJ=betaJ, v=v, deltas=deltas.tolist(), g2=vals,
                pole_at=math.pi + math.pi / v,
                endpoint_margin=math.pi / 2.0,
                all_finite_positive=bool(all(np.isfinite(vals))
                                         and all(x > 0 for x in vals)))


# ---------------------------------------------------------------------------
# level 3: exact chord theory
# ---------------------------------------------------------------------------
def chord_scan(lam: float, betaJ: float, p: int = 4, nmax: int = 70,
               thL_grid=(0.0, 0.4, 0.8, 1.2)) -> dict:
    """Connected crossed ratio vs Lorentzian theta-time at Delta = 0,
    pi/(2v), pi/v -- exact in lambda.  Reports the trend direction and a
    truncation check at the antipodal endpoint."""
    v = v_of_betaJ(betaJ)
    Nflav = 2.0 * p * p / lam
    from chord import beta_chord_from_betaJ
    beta_c = beta_chord_from_betaJ(betaJ, lam)
    cc = ChordCorrelators(lam, p, beta_c, nmax)
    t1, t2, t3, t4 = CONFIG_CROSSED
    out_rows = []
    for d in (0.0, 0.5 * math.pi / v, math.pi / v):
        curve = []
        for thL in thL_grid:
            cfg = (t1 + d + 1j * thL, t2 + d + 1j * thL, t3, t4)
            val = cc.connected_ratio(cfg, 'c', Nflav)
            curve.append(dict(thL=thL, re=val.real, im=val.imag))
        res = [c["re"] for c in curve]
        out_rows.append(dict(delta=float(d), curve=curve,
                             re_move=res[-1] - res[0],
                             direction="up" if res[-1] > res[0] else "down"))
    # truncation check at the antipodal endpoint, largest thL
    d = math.pi / v
    cfg = (t1 + d + 1.2j, t2 + d + 1.2j, t3, t4)
    v1 = cc.connected_ratio(cfg, 'c', Nflav)
    cc2 = ChordCorrelators(lam, p, beta_c, nmax + 20)
    v2 = cc2.connected_ratio(cfg, 'c', Nflav)
    drift = abs(v1 - v2) / max(abs(v2), 1e-12)
    return dict(lam=lam, betaJ=betaJ, p=p, nmax=nmax, v=v, Nflav=Nflav,
                rows=out_rows, truncation_drift=float(drift))


def main():
    t0 = time.time()
    betaJ = 2.0
    payload = dict(
        provenance="signflip_probe.py -- exploratory probe of the "
                   "antipodal half-fake-period shift as a candidate "
                   "OTOC-sign-flipping map (DISCRIMINATOR.md); necessary, "
                   "not sufficient",
        closed_form=[closed_form_scan(bJ) for bJ in (0.5, 2.0, 5.0)],
        g2=g2_along_shift(betaJ),
        chord=[chord_scan(lam, betaJ, nmax=80,
                          thL_grid=(0.0, 0.5, 1.0, 1.5, 2.0))
               for lam in (2.0, 1.6, 1.2, 1.0, 0.8, 0.6)])
    payload["runtime_s"] = round(time.time() - t0, 1)
    out = HERE / "signflip_probe.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    for cf in payload["closed_form"]:
        print(f"betaJ={cf['betaJ']}: Re a unshifted = "
              f"{cf['re_a_unshifted']:+.6f} -> antipodal = "
              f"{cf['re_a_antipodal']:+.6f}")
    print(f"g2 along shift finite and positive: "
          f"{payload['g2']['all_finite_positive']}")
    for ch in payload["chord"]:
        dirs = {r['delta']: r['direction'] for r in ch['rows']}
        print(f"chord lam={ch['lam']}: directions {dirs} "
              f"drift={ch['truncation_drift']:.1e}")
    print(f"wrote {out} ({payload['runtime_s']} s)")


if __name__ == "__main__":
    main()
