#!/usr/bin/env python3
r"""
Can a contour shift flip the scrambling sign?  A probe.
=======================================================

Here is the dS-dictionary problem, as this repo measures it.  Under the
flat time map DSSYK scrambles, with a(v) < 0 at every temperature.  A dS
observer needs the opposite sign.  The Euclidean fold does not deliver
it.

This module tests the next candidate.  Displace the probe pair by HALF
THE FAKE PERIOD around the thermal circle, Delta = pi/v in Euclidean
angle.  The fake period is 2 pi/v.

Geometrically that is an antipodal shift on the fake disk.  It is also
what the Narovlansky-Verlinde dressing does: their operators sit at
+-i beta_dS/4, which is half a period apart, and they call the two
insertions antipodal.

THREE LEVELS, ONE JSON
----------------------
1. Closed form (large-q, `antiscrambling.otoc_ratio`).  Shift the
   (theta1, theta2) pair by Delta.  Extract the growth coefficient
   a(Delta) by a complex fit in the late-time window.

   The pole structure predicts a(Delta) = a(0) exp(-i v Delta).  So at
   Delta = pi/v the coefficient flips exactly: Re a goes from a(0) to
   -a(0).

   Read that for what it is.  Delta enters the growing mode's phase, so
   half a period negates it, and the algebra is short.  The probe
   verifies the rotation law and the flip.  The content is in levels 2
   and 3.

2. Two-point physicality.  The shifted W-V separation is
   theta = pi/2 + Delta.  The continued large-q 2-pt function stays real
   and positive up to its pole at theta = pi + pi/v.  So the antipodal
   point sits INSIDE the physical window, with a margin of pi/2 from the
   pole.  The probe records the values and the margin.

3. Exact chord theory at finite lambda (`chord_fold.ChordCorrelators`).
   The same shifted contours, exact in lambda, truncation-checked.

   The observable: does Re F/F_d move DOWN with Lorentzian time
   (scrambling) at Delta = 0, and UP (the anti-scrambling side) at
   Delta = pi/v?  Plus the fitted coefficient wherever the window allows.

HONESTY
-------
A sign flip under a contour shift is NECESSARY for this candidate map.
It is not sufficient.  It does not by itself build the dS dictionary.

The full correlator suite must stay consistent, and the algebra story --
why the observer's operators live there -- is NV's construction to
complete.

What this probe measures is whether the mechanism survives exactness in
lambda.  That is the question the closed form cannot answer.

Output: signflip_probe.json.  Asserted in test_signflip.py.  Narrated in
DISCRIMINATOR.md.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

from antiscrambling import G2_thermal, otoc_ratio
from chord import Contour, Sector0, beta_chord_from_betaJ, q_of_lambda
from chord_fold import CONFIG_CROSSED, CONFIG_TOC, ChordCorrelators
from largeq_anchor import v_of_betaJ

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# level 1: closed form
# ---------------------------------------------------------------------------
def growth_coefficient_shifted(v: float, delta: float,
                               config=CONFIG_CROSSED) -> complex:
    """Coefficient of exp(v thL) in the closed-form crossed ratio, with the
    (theta1, theta2) pair shifted by `delta` in Euclidean angle.  The fit
    is a complex least-squares fit on the basis
    {1, e^{v thL}, e^{-v thL}} over a late window.  antiscrambling.py
    validates the same technique at delta = 0."""
    t1, t2, t3, t4 = config
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


# ---------------------------------------------------------------------------
# level 4: the antipodal frame is thermal at the tomperature (exact)
# ---------------------------------------------------------------------------
def detailed_balance(lam: float, betaJ: float, p: int = 4,
                     nmax: int = 120) -> dict:
    """The antipodal frame is thermal at beta_eff = beta + 2 dtau.

    The statement has two parts with different logical status, and the
    returned dict separates them:

    * The per-line LAW.  In the Lehmann representation, shifting one
      probe operator by Euclidean time dtau multiplies each line (j, k)
      by exp(-dtau (Ej - Ek)).  With the symmetric vertex M this forces
      detailed balance at beta_eff = beta + 2 dtau for any dtau.  That
      is arithmetic, not a measurement.  `max_line_deviation` only
      confirms the construction (float-level).
    * The independent CHECK (`engine_check_dev`).  The same shifted
      correlator is evaluated on the contour engine, the stack walker
      with matter-sector `expm_multiply`.  That code path contains no
      Lehmann weights.  The result is compared against the
      per-line-shifted sum, at two Lorentzian times and both shift
      signs.  This is what makes the Lehmann representation (and so the
      law) a verified property of the chord theory rather than of this
      function.

    Physics: the sign-flipping shift is dtau = beta_fake/2.  So at the
    T = infinity hologram point (beta -> 0) the antipodal frame obeys
    detailed balance at the tomperature, beta_eff = beta_fake.  (This is
    what NV's +-i beta/4 splitting produces, with their beta_dS
    identified as the fake temperature.)"""
    v = v_of_betaJ(betaJ)
    q = q_of_lambda(lam)
    beta_c = beta_chord_from_betaJ(betaJ, lam)
    dtau_c = 0.5 * beta_c / v                 # half the fake period
    s = Sector0(nmax, q)
    v0 = s.V[0, :]
    M = s.V.T @ ((q ** ((1.0 / p) * np.arange(nmax + 1)))[:, None] * s.V)
    # Wightman ordering (the object KMS is about):
    # G_W(t) = sum_jk W_jk e^{i(Ej-Ek)t},  W_jk = v0_j v0_k M_jk e^{-b Ej}
    W = np.outer(v0 * np.exp(-beta_c * s.E), v0) * M
    dE = np.subtract.outer(s.E, s.E)
    Wsh = W * np.exp(-dtau_c * dE)
    beta_eff = beta_c + 2.0 * dtau_c
    mask = np.abs(W) > 1e-10 * np.abs(W).max()
    ratio = np.log(np.abs(Wsh[mask] / Wsh.T[mask]))
    ideal = -beta_eff * dE[mask]
    dev = float(np.max(np.abs(ratio - ideal)))
    # independent check: contour engine vs the per-line-shifted Lehmann
    # sum, R(u) = sum W e^{u dE} = walker[('A', u), ('A', beta - u)]
    nmax_chk = min(nmax, 60)
    s2 = Sector0(nmax_chk, q)
    v02 = s2.V[0, :]
    M2 = s2.V.T @ ((q ** ((1.0 / p) * np.arange(nmax_chk + 1)))[:, None]
                   * s2.V)
    W2 = np.outer(v02 * np.exp(-beta_c * s2.E), v02) * M2
    dE2 = np.subtract.outer(s2.E, s2.E)
    eng = Contour(q, nmax_chk, {'A': 1.0 / p})
    dev_eng = 0.0
    for t in (0.4, 1.1):
        for sgn in (+1.0, -1.0):
            u = sgn * dtau_c + 1j * t
            lehmann = complex(np.sum(W2 * np.exp(u * dE2)))
            walker = eng.correlator([('A', u), ('A', beta_c - u)])
            dev_eng = max(dev_eng, abs(walker - lehmann)
                          / max(abs(lehmann), 1e-300))
    return dict(lam=lam, betaJ=betaJ, v=v, beta_c=beta_c,
                dtau_c=dtau_c, beta_eff=beta_eff,
                beta_fake_c=beta_c / v,
                max_line_deviation=dev,
                engine_check_dev=float(dev_eng))


# ---------------------------------------------------------------------------
# level 5: self-calibrated shift at strong coupling (honest negative)
# ---------------------------------------------------------------------------
def calibrated_scan(lam: float, betaJ: float, p: int = 4,
                    nmax: int = 80) -> dict:
    """Measure an effective growth exponent kappa from the unshifted curve
    (three-point log-difference on thL = 0.5..2.0).  Then try the flip
    with Delta = pi/kappa as well as the large-q Delta = pi/v.

    CALIBRATION.  The same estimator applied to the large-q closed form
    itself (antiscrambling.otoc_ratio) on the same window gives
    kappa_closed = 0.94 at betaJ = 2, not v = 0.47: the cosh(v thL)
    growing/decaying pair contaminates early differences, and the closed
    form only reaches v for thL >~ 6.  So kappa is a window-dependent
    effective rate, recorded next to its closed-form calibration
    (kappa_closed_same_window) and their ratio.  The exact-chord values
    (0.83-0.88) sit BELOW the closed form's own early-window value; they
    are not evidence that the magnitude rate exceeds v.  An earlier
    version of this docstring read kappa > v as "phase and magnitude
    rates decouple at finite lambda"; that reading is withdrawn.

    Two findings survive.  (1) At strong coupling both shifts fail: the
    flip is genuinely semiclassical, not a rate-calibration artifact.
    (2) Where the flip works (lam <= 1.7), it works with pi/v and NOT with
    pi/kappa_eff -- the phase period of the growing mode is the large-q
    fake period, and pi/kappa_eff (~3.8) lands near pi/(2v) = 3.3 where
    the rotation law puts Re a = 0, so its failure is expected."""
    v = v_of_betaJ(betaJ)
    Nflav = 2.0 * p * p / lam
    cc = ChordCorrelators(lam, p, beta_chord_from_betaJ(betaJ, lam), nmax)
    t1, t2, t3, t4 = CONFIG_CROSSED
    h, thLs = 0.5, (0.5, 1.0, 1.5, 2.0)

    def three_point(rs):
        d1, d2, d3 = rs[1] - rs[0], rs[2] - rs[1], rs[3] - rs[2]
        return 0.5 * (math.log(abs(d2 / d1)) + math.log(abs(d3 / d2))) / h

    r = [cc.connected_ratio((t1 + 1j * t, t2 + 1j * t, t3, t4),
                            'c', Nflav).real for t in thLs]
    kappa = three_point(r)
    from antiscrambling import otoc_ratio
    r_closed = [otoc_ratio(t1 + 1j * t, t2 + 1j * t, t3, t4, v).real
                for t in thLs]
    kappa_closed = three_point(r_closed)
    out = dict(lam=lam, betaJ=betaJ, v_largeq=v, kappa=float(kappa),
               kappa_closed_same_window=float(kappa_closed),
               kappa_over_closed=float(kappa / kappa_closed),
               estimator="three-point log-difference, thL = 0.5..2.0 "
                         "(window-dependent; see docstring)")
    for name, d in [("largeq", math.pi / v), ("calibrated", math.pi / kappa)]:
        rs = [cc.connected_ratio((t1 + d + 1j * t, t2 + d + 1j * t, t3, t4),
                                 'c', Nflav).real for t in (0.0, 1.0, 2.0)]
        out[name] = dict(delta=float(d), re_move=rs[-1] - rs[0],
                         direction="up" if rs[-1] > rs[0] else "down")
    return out


# ---------------------------------------------------------------------------
# level 6: TOC control and configuration robustness
# ---------------------------------------------------------------------------
# A GENUINELY inequivalent crossed configuration.  The previous choice
# (1.6pi, 0.6pi, 1.1pi, 0.1pi) was CONFIG_CROSSED rigidly translated by
# 0.1 pi, so it shares every translation invariant and -- now that
# otoc_ratio is correctly translation invariant (antiscrambling.py, the
# 2026-08-24 gauge fix) -- carried provably zero information.  This one
# differs in all three invariants:
#   CONFIG_CROSSED: th12- = 1.00 pi, th34- = 1.00 pi, th12+ - th34+ = 1.00 pi
#   CONFIG_CROSSED_B: 1.20 pi, 1.10 pi, 0.90 pi
# and it keeps the crossed ordering th1 > th3 > th2 > th4 inside the
# principal window.  Being asymmetric, its a(0) has a genuine imaginary
# part -- which is exactly why it is a real test of the rotation law.
CONFIG_CROSSED_B = (1.7 * math.pi, 0.5 * math.pi, 1.2 * math.pi,
                    0.1 * math.pi)


def toc_control(lam: float, betaJ: float, p: int = 4,
                nmax: int = 80) -> dict:
    """The time-ordered 4-pt under the same shift: it must stay bounded
    (no exponential growth) while the crossed one flips.  Reports the max
    magnitude ratio shifted/unshifted over the Lorentzian grid."""
    v = v_of_betaJ(betaJ)
    Nflav = 2.0 * p * p / lam
    cc = ChordCorrelators(lam, p, beta_chord_from_betaJ(betaJ, lam), nmax)
    t1, t2, t3, t4 = CONFIG_TOC
    d = math.pi / v
    m0, m1 = [], []
    for thL in (0.0, 1.0, 2.0):
        u = cc.connected_ratio((t1 + 1j * thL, t2 + 1j * thL, t3, t4),
                               'u', Nflav)
        s_ = cc.connected_ratio((t1 + d + 1j * thL, t2 + d + 1j * thL,
                                 t3, t4), 'u', Nflav)
        m0.append(abs(u))
        m1.append(abs(s_))
    return dict(lam=lam, max_unshifted=max(m0), max_shifted=max(m1),
                bounded=bool(max(m1) < 20.0 * max(max(m0), 1.0)))


def config_b_scan(betaJ: float) -> dict:
    """Second, GENUINELY inequivalent crossed configuration (see
    CONFIG_CROSSED_B): the closed-form rotation law
    a(Delta) = a(0) exp(-i v Delta) is configuration-independent, so
    Re a flips sign at Delta = pi/v wherever it starts -- here from a
    complex a(0), which the symmetric primary configuration cannot
    exhibit.  The exact chord flip also repeats at lambda = 1."""
    v = v_of_betaJ(betaJ)
    t1, t2, t3, t4 = CONFIG_CROSSED_B
    a0 = growth_coefficient_shifted(v, 0.0, CONFIG_CROSSED_B)
    a1 = growth_coefficient_shifted(v, math.pi / v, CONFIG_CROSSED_B)
    lam, p = 1.0, 4
    Nflav = 2.0 * p * p / lam
    cc = ChordCorrelators(lam, p, beta_chord_from_betaJ(betaJ, lam), 80)
    d = math.pi / v
    rs0 = [cc.connected_ratio((t1 + 1j * t, t2 + 1j * t, t3, t4),
                              'c', Nflav).real for t in (0.0, 1.0, 2.0)]
    rs1 = [cc.connected_ratio((t1 + d + 1j * t, t2 + d + 1j * t, t3, t4),
                              'c', Nflav).real for t in (0.0, 1.0, 2.0)]
    return dict(betaJ=betaJ, re_a0=a0.real, re_a_antipodal=a1.real,
                rotation_error=abs(a1 + a0) / abs(a0),
                chord_lam1_unshifted="down" if rs0[-1] < rs0[0] else "up",
                chord_lam1_antipodal="up" if rs1[-1] > rs1[0] else "down")


# ---------------------------------------------------------------------------
# the DEEP campaign (--deep): the hologram point, the betaJ = pi shape
# prediction, the MSS-evasion locator, and the observer rate
# ---------------------------------------------------------------------------
def _ratio_z(cc, zs_path, Nflav):
    """Connected crossed ratio from explicit path-ordered slot times
    (works at beta_c = 0, where the theta parametrization degenerates)."""
    F = -cc.f4(['A', 'B', 'A', 'B'], zs_path)
    z1, z3, z2, z4 = zs_path
    Fd = cc.g2(z1, z2) * cc.g2(z3, z4)
    return Nflav * (F - Fd) / Fd


def hologram_point(lam: float, p: int = 4, nmax: int = 100,
                   ts=tuple(0.25 * k for k in range(15))) -> dict:
    """The flip at the ACTUAL hologram point, beta_c = 0.  The tomperature
    is beta_f = pi/Jc there (equivalently lambda_L = 2 Jc = 2 pi/beta_f).
    W-pair regulator: Euclidean split 2*eps with eps = 0.15 beta_f
    (stated convention).  Finding: the antipodal frame anti-scrambles in
    a TRANSIENT window whose extent grows as lambda falls -- scrambling
    reasserts at later times, consistent with a semiclassical mechanism
    fought by lambda-corrections.

    The window end is the maximum of the antipodal curve, located by a
    parabola through the three grid points around the largest sample
    (t_window), with the grid-quantized value kept as t_window_grid and a
    `censored` flag if the curve is still rising at the last grid point
    (the grid runs to t = 3.5, past every window through lambda = 0.3)."""
    Jc = math.sqrt(lam) * math.exp(lam / 8.0)
    beta_f = math.pi / Jc
    tau_s = 0.5 * beta_f
    eps = 0.15 * beta_f
    Nflav = 2.0 * p * p / lam
    cc = ChordCorrelators(lam, p, 0.0, nmax)

    def curve(shift):
        out = []
        for t in ts:
            zs = [-1j * (shift + eps) + t, -1j * eps,
                  -1j * (shift - eps) + t, +1j * eps]
            out.append(float(_ratio_z(cc, zs, Nflav).real))
        return out

    un, an = curve(0.0), curve(tau_s)
    # the anti-scrambling window: the maximum of the antipodal curve
    imax = 0
    for i in range(1, len(ts)):
        if an[i] > an[i - 1]:
            imax = i
        else:
            break
    t_window_grid = ts[imax]
    censored = imax == len(ts) - 1
    if 0 < imax < len(ts) - 1:
        # parabola through (t_{i-1}, t_i, t_{i+1}): vertex position
        (x0, y0), (x1, y1), (x2, y2) = [(ts[j], an[j]) for j in (imax - 1, imax, imax + 1)]
        denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
        a2 = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
        b1 = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom
        t_window = float(-b1 / (2.0 * a2)) if a2 < 0 else float(t_window_grid)
    else:
        t_window = float(t_window_grid)
    cc2 = ChordCorrelators(lam, p, 0.0, nmax + 20)
    zs = [-1j * (tau_s + eps) + 1.0, -1j * eps,
          -1j * (tau_s - eps) + 1.0, +1j * eps]
    r2 = _ratio_z(cc2, zs, Nflav)
    drift = abs(_ratio_z(cc, zs, Nflav) - r2) / max(abs(r2), 1e-12)
    return dict(lam=lam, beta_f=beta_f, tau_s=tau_s, eps=eps,
                ts=list(ts), unshifted=un, antipodal=an,
                t_window=t_window, t_window_grid=float(t_window_grid),
                censored=bool(censored), Jc=Jc,
                window_in_Jc_units=t_window * Jc,
                window_grid_in_Jc_units=float(t_window_grid) * Jc,
                window_quantization_Jc=0.125 * Jc,
                truncation_drift=float(drift))


def thermal_shape(lam: float, p: int = 4, nmax: int = 140) -> dict:
    """The SHAPE of the antipodal observer's 2-pt at the hologram point.

    The observable is <0| e^{zT} O e^{-zT} O |0> with z = tau_s + i t: the
    T = infinity 2-pt continued to the complex Euclidean time z, Euclidean
    shift tau_s = beta_fake/2 (the antipodal point) plus real time t.  Its
    large-q form is G(z) = [1/cos(Jc z)]^{2 Delta}: in real time that is
    the familiar sech(Jc t), in Euclidean time the fake-period growth
    1/cos(Jc tau) with its singularity at tau = beta_fake/2.  So at the
    antipodal shift the target is |1/cos(Jc(tau_s + i t))|^{2 Delta}
    = |sech(Jc(t - i tau_s))|^{2 Delta}, which at t = 0 is the sec(pi/2)
    singularity itself; the exact chord value there is finite because the
    chord spectrum is bounded (an edge cutoff, not a fold pathology).

    Measured verdict (signflip_deep.json): for t >= 0.25/Jc the antipodal
    frame's 2-pt IS the continued large-q form, with a deviation that
    SHRINKS toward semiclassics (0.46 -> 0.18 over lambda = 2 -> 0.25);
    the Euclidean 2-pt likewise tracks the fake-period growth toward the
    singularity (deviation 0.49 -> 0.26 at tau <= 0.45 beta_fake).  The
    naive Gibbs guess (the observer thermal at beta J = pi in the ordinary
    finite-temperature large-q form, MS Eq. 2.18 continued to real time)
    fails, worse toward semiclassics: the tomperature thermality is the
    detailed-balance identity of the shifted sector, not a Gibbs state.

    An earlier version of this function continued the real-time form with
    the roles of t and tau_s swapped, cosh(Jc(tau_s + i t)), and tested
    'periodicity' in real time; both instruments were miscalibrated, and
    the 'shapes fail worse toward semiclassics' negative they produced is
    withdrawn (the historical numbers stay in the git history)."""
    Jc = math.sqrt(lam) * math.exp(lam / 8.0)
    beta_f = math.pi / Jc
    tau_s = 0.5 * beta_f
    q = q_of_lambda(lam)
    dpsi = 1.0 / p
    s = Sector0(nmax, q)
    v0 = s.V[0, :]
    M = s.V.T @ ((q ** (dpsi * np.arange(nmax + 1)))[:, None] * s.V)

    def g_obs_at(z):
        return complex((v0 * np.exp(z * s.E)) @ M @ (v0 * np.exp(-z * s.E)))

    # real time along the antipodal shift; t = 0 is the sec(pi/2) point and
    # is excluded from the shape comparison (kept in the curve)
    ts = np.linspace(0.0, 2.0 / Jc, 25)
    g_obs = np.array([g_obs_at(tau_s + 1j * t) for t in ts])
    z = Jc * (tau_s + 1j * ts)
    g_true = np.abs(1.0 / np.cos(z)) ** (2 * dpsi)
    g_true[0] = float("nan")                   # the singular point
    win = ts >= 0.25 / Jc                      # past the edge cutoff
    dev_sech = float(np.max(np.abs(np.abs(g_obs[win]) - g_true[win])
                            / g_true[win]))
    # the naive Gibbs guess (negative control): the full MS finite-T form
    # at beta J = pi, continued to real time
    v_f = v_of_betaJ(math.pi)
    g_th = np.abs(np.cos(np.pi * v_f / 2.0)
                  / np.cos(np.pi * v_f / 2.0 - 1j * np.pi * v_f * ts / beta_f)
                  ) ** (2 * dpsi)
    dev_gibbs = float(np.max(np.abs(np.abs(g_obs[win]) - g_th[win]) / g_th[win]))
    # the Euclidean fake-period growth toward the singularity
    taus = np.linspace(0.0, 0.45 * beta_f, 10)
    g_eu = np.array([abs(g_obs_at(tau)) for tau in taus])
    g_eu_true = (1.0 / np.cos(Jc * taus)) ** (2 * dpsi)
    dev_eu = float(np.max(np.abs(g_eu - g_eu_true) / g_eu_true))
    return dict(lam=lam, v_f=v_f, beta_f=beta_f,
                ts=[float(x) for x in ts],
                g_obs_abs=[float(abs(x)) for x in g_obs],
                g_sech_pred=[float(x) for x in g_true],
                g_gibbs_pred=[float(x) for x in g_th],
                shape_window="t >= 0.25/Jc",
                max_rel_dev_sech=dev_sech,
                max_rel_dev_gibbs=dev_gibbs,
                taus=[float(x) for x in taus],
                g_euclid_abs=[float(x) for x in g_eu],
                g_euclid_pred=[float(x) for x in g_eu_true],
                max_rel_dev_euclid=dev_eu,
                target="|1/cos(Jc (tau_s + i t))|^(2 Delta), t >= 0.25/Jc")


def mss_strip(betaJ: float = 2.0, n_y: int = 33) -> dict:
    """Where does the observer frame break the chaos-bound hypotheses?
    Scan the closed-form growth coefficient around the observer's full
    fake circle: Re a(Delta + y) = Re[a(0) e^{-i v (Delta+y)}] oscillates,
    so there is NO Euclidean placement of the four operators at which the
    observer's regularized OTOC has the definite-sign structure the MSS
    proof requires globally -- the correlator is unbounded on the strip
    (|a| is constant while the phase winds).  The located evasion: the
    observer's effective thermality is a detailed-balance statement, not
    the strip-analytic KMS structure the bound assumes."""
    v = v_of_betaJ(betaJ)
    ys = np.linspace(0.0, 2.0 * math.pi / v, n_y)
    re_a = []
    for y in ys:
        a = growth_coefficient_shifted(v, math.pi / v + float(y))
        re_a.append(float(a.real))
    re_a = np.array(re_a)
    return dict(betaJ=betaJ, v=v, ys=ys.tolist(), re_a=re_a.tolist(),
                oscillates=bool(re_a.max() > 0.5 and re_a.min() < -0.5),
                abs_a_constant=bool(np.allclose(
                    [abs(complex(growth_coefficient_shifted(v, d)))
                     for d in (0.0, math.pi / v, 1.7 * math.pi / v)],
                    abs(complex(growth_coefficient_shifted(v, 0.0))),
                    rtol=1e-6)))


def observer_rate(lam: float, betaJ: float = 2.0, p: int = 4,
                  nmax: int = 80) -> dict:
    """Growth exponent of the FLIPPED mode vs the dS expectation
    2 pi/beta_eff.  Fit kappa_obs from the antipodal curve where the flip
    is robust; compare against 2 pi v/beta_c (= 2 pi/beta_fake) and
    2 pi/beta_eff."""
    v = v_of_betaJ(betaJ)
    Nflav = 2.0 * p * p / lam
    beta_c = beta_chord_from_betaJ(betaJ, lam)
    cc = ChordCorrelators(lam, p, beta_c, nmax)
    t1, t2, t3, t4 = CONFIG_CROSSED
    d = math.pi / v
    h, thLs = 0.5, (0.5, 1.0, 1.5, 2.0)
    r = [cc.connected_ratio((t1 + d + 1j * t, t2 + d + 1j * t, t3, t4),
                            'c', Nflav).real for t in thLs]
    d1, d2, d3 = r[1] - r[0], r[2] - r[1], r[3] - r[2]
    kap = 0.5 * (math.log(abs(d2 / d1)) + math.log(abs(d3 / d2))) / h
    # rates are per theta_L; convert to chord time: t_c = beta_c thL/(2 pi)
    to_time = 2.0 * math.pi / beta_c
    return dict(lam=lam, betaJ=betaJ,
                kappa_obs_time=kap * to_time,
                rate_fake=2.0 * math.pi / (beta_c / v),
                rate_eff=2.0 * math.pi / (beta_c + beta_c / v))


def main_deep():
    t0 = time.time()
    payload = dict(
        provenance="signflip_probe.py --deep -- does the antipodal map "
                   "pan out: the hologram point (transient window), the "
                   "betaJ = pi shape prediction, the MSS-evasion "
                   "locator, and the observer rate",
        hologram=[hologram_point(lam) for lam in (1.0, 0.6, 0.4, 0.3)],
        shape=[thermal_shape(lam) for lam in (2.0, 1.0, 0.5, 0.25)],
        mss=mss_strip(),
        rate=[observer_rate(lam) for lam in (1.0, 0.8, 0.6)])
    payload["runtime_s"] = round(time.time() - t0, 1)
    out = HERE / "signflip_deep.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    for h in payload["hologram"]:
        print(f"hologram lam={h['lam']}: window t*Jc = "
              f"{h['window_in_Jc_units']:.2f}  drift={h['truncation_drift']:.1e}")
    for s_ in payload["shape"]:
        print(f"shape lam={s_['lam']}: dev vs continued large-q form = "
              f"{s_['max_rel_dev_sech']:.3f}  vs Gibbs (neg. control) = "
              f"{s_['max_rel_dev_gibbs']:.3f}  Euclidean fake-period growth "
              f"dev = {s_['max_rel_dev_euclid']:.3f}")
    print(f"mss: coefficient phase winds = {payload['mss']['oscillates']}, "
          f"|a| constant = {payload['mss']['abs_a_constant']}")
    for r in payload["rate"]:
        print(f"rate lam={r['lam']}: kappa_obs={r['kappa_obs_time']:.3f} "
              f"vs 2pi/beta_fake={r['rate_fake']:.3f} "
              f"2pi/beta_eff={r['rate_eff']:.3f}")
    print(f"wrote {out} ({payload['runtime_s']} s)")


def main():
    t0 = time.time()
    betaJ = 2.0
    payload = dict(
        provenance="signflip_probe.py -- the antipodal half-fake-period "
                   "shift as a candidate OTOC-sign-flipping map "
                   "(DISCRIMINATOR.md); necessary, not sufficient",
        closed_form=[closed_form_scan(bJ) for bJ in (0.5, 2.0, 5.0)],
        g2=g2_along_shift(betaJ),
        chord=[chord_scan(lam, betaJ, nmax=80,
                          thL_grid=(0.0, 0.5, 1.0, 1.5, 2.0))
               for lam in (2.0, 1.8, 1.7, 1.6, 1.5, 1.2, 1.0, 0.8, 0.6)],
        detailed_balance=[detailed_balance(lam, bJ)
                          for lam in (2.0, 1.0)
                          for bJ in (0.02, 2.0)],
        calibrated=[calibrated_scan(lam, betaJ)
                    for lam in (4.0, 3.0, 2.0, 1.0, 0.6)],
        toc=[toc_control(lam, betaJ) for lam in (2.0, 1.0)],
        config_b=config_b_scan(betaJ))
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
    import sys
    if "--deep" in sys.argv:
        main_deep()
    else:
        main()
