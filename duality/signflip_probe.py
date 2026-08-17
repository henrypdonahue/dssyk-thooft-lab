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


# ---------------------------------------------------------------------------
# level 4: the antipodal frame is thermal at the tomperature (exact)
# ---------------------------------------------------------------------------
def detailed_balance(lam: float, betaJ: float, p: int = 4,
                     nmax: int = 120) -> dict:
    """Exact spectral statement: shifting one probe operator by Euclidean
    time dtau multiplies each Lehmann line (j, k) by exp(-dtau (Ej - Ek)),
    so the shifted-frame weights obey detailed balance at
    beta_eff = beta + 2 dtau -- an exact identity at ANY lambda, checked
    per line here.  Physics: the sign-flipping shift is dtau =
    beta_fake/2, so at the T = infinity hologram point (beta -> 0) the
    antipodal frame satisfies EXACT detailed balance at the tomperature:
    beta_eff = beta_fake.  (This is also what NV's +-i beta/4 splitting
    produces, with their beta_dS identified as the fake temperature.)"""
    from chord import Sector0, beta_chord_from_betaJ, q_of_lambda
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
    # shift one operator by Euclidean dtau (increasing the separation):
    # each line picks up e^{-dtau (Ej - Ek)}
    Wsh = W * np.exp(-dtau_c * dE)
    beta_eff = beta_c + 2.0 * dtau_c
    mask = np.abs(W) > 1e-10 * np.abs(W).max()
    ratio = np.log(np.abs(Wsh[mask] / Wsh.T[mask]))
    ideal = -beta_eff * dE[mask]
    dev = float(np.max(np.abs(ratio - ideal)))
    return dict(lam=lam, betaJ=betaJ, v=v, beta_c=beta_c,
                dtau_c=dtau_c, beta_eff=beta_eff,
                beta_fake_c=beta_c / v,
                max_line_deviation=dev)


# ---------------------------------------------------------------------------
# level 5: self-calibrated shift at strong coupling (honest negative)
# ---------------------------------------------------------------------------
def calibrated_scan(lam: float, betaJ: float, p: int = 4,
                    nmax: int = 80) -> dict:
    """Measure the exact growth exponent kappa from the unshifted curve
    (three-point log-difference), then try the flip with Delta = pi/kappa
    as well as the large-q Delta = pi/v.  Two findings.  (1) At strong
    coupling both shifts fail: the flip is genuinely semiclassical, not a
    rate-calibration artifact.  (2) Where the flip works (lam <= 1.7), it
    works with pi/v and NOT with pi/kappa, even though kappa > v: the
    phase period of the growing mode is set by the large-q v while its
    magnitude grows faster -- phase and magnitude rates decouple at
    finite lambda.  kappa itself is a new measured number."""
    from chord import beta_chord_from_betaJ
    v = v_of_betaJ(betaJ)
    Nflav = 2.0 * p * p / lam
    cc = ChordCorrelators(lam, p, beta_chord_from_betaJ(betaJ, lam), nmax)
    t1, t2, t3, t4 = CONFIG_CROSSED
    h, thLs = 0.5, (0.5, 1.0, 1.5, 2.0)
    r = [cc.connected_ratio((t1 + 1j * t, t2 + 1j * t, t3, t4),
                            'c', Nflav).real for t in thLs]
    d1, d2, d3 = r[1] - r[0], r[2] - r[1], r[3] - r[2]
    kappa = 0.5 * (math.log(abs(d2 / d1)) + math.log(abs(d3 / d2))) / h
    out = dict(lam=lam, betaJ=betaJ, v_largeq=v, kappa=float(kappa))
    for name, d in [("largeq", math.pi / v), ("calibrated", math.pi / kappa)]:
        rs = [cc.connected_ratio((t1 + d + 1j * t, t2 + d + 1j * t, t3, t4),
                                 'c', Nflav).real for t in (0.0, 1.0, 2.0)]
        out[name] = dict(delta=float(d), re_move=rs[-1] - rs[0],
                         direction="up" if rs[-1] > rs[0] else "down")
    return out


# ---------------------------------------------------------------------------
# level 6: TOC control and configuration robustness
# ---------------------------------------------------------------------------
CONFIG_CROSSED_B = (1.6 * math.pi, 0.6 * math.pi, 1.1 * math.pi,
                    0.1 * math.pi)


def toc_control(lam: float, betaJ: float, p: int = 4,
                nmax: int = 80) -> dict:
    """The time-ordered 4-pt under the same shift: it must stay bounded
    (no exponential growth) while the crossed one flips.  Reports the max
    magnitude ratio shifted/unshifted over the Lorentzian grid."""
    from chord import beta_chord_from_betaJ
    v = v_of_betaJ(betaJ)
    Nflav = 2.0 * p * p / lam
    cc = ChordCorrelators(lam, p, beta_chord_from_betaJ(betaJ, lam), nmax)
    from chord_fold import CONFIG_TOC
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
    """Second crossed configuration: the closed-form rotation law
    a(Delta) = a(0) exp(-i v Delta) is configuration-independent, so
    Re a flips sign at Delta = pi/v wherever it starts; and the exact
    chord flip repeats at lambda = 1."""
    v = v_of_betaJ(betaJ)
    t1, t2, t3, t4 = CONFIG_CROSSED_B

    def coeff(delta):
        thLs = np.array([4.0, 5.0, 6.0, 7.0]) / v
        vals = np.array([complex(otoc_ratio(t1 + delta + 1j * t,
                                            t2 + delta + 1j * t,
                                            t3, t4, v)) for t in thLs])
        A = np.vstack([np.ones_like(thLs), np.exp(v * thLs),
                       np.exp(-v * thLs)]).T
        c, *_ = np.linalg.lstsq(A.astype(complex), vals, rcond=None)
        return complex(c[1])

    a0 = coeff(0.0)
    a1 = coeff(math.pi / v)
    from chord import beta_chord_from_betaJ
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
# the DEEP campaign (--deep): the hologram point, the beta.J = pi shape
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
                   ts=(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)) -> dict:
    """The flip at the ACTUAL hologram point, beta_c = 0.  The tomperature
    is beta_f = pi/Jc there (equivalently lambda_L = 2 Jc = 2 pi/beta_f).
    W-pair regulator: Euclidean split 2*eps with eps = 0.15 beta_f
    (stated convention).  Finding: the antipodal frame anti-scrambles in
    a TRANSIENT window whose extent grows as lambda falls -- scrambling
    reasserts at later times, consistent with a semiclassical mechanism
    fought by lambda-corrections."""
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
    # the anti-scrambling window: largest t up to which the antipodal
    # curve is still rising
    t_window = 0.0
    for i in range(1, len(ts)):
        if an[i] > an[i - 1]:
            t_window = ts[i]
        else:
            break
    cc2 = ChordCorrelators(lam, p, 0.0, nmax + 20)
    zs = [-1j * (tau_s + eps) + 1.0, -1j * eps,
          -1j * (tau_s - eps) + 1.0, +1j * eps]
    drift = abs(_ratio_z(cc, zs, Nflav) - _ratio_z(cc2, zs, Nflav)) / max(
        abs(_ratio_z(cc2, zs, Nflav)), 1e-12)
    return dict(lam=lam, beta_f=beta_f, tau_s=tau_s, eps=eps,
                ts=list(ts), unshifted=un, antipodal=an,
                t_window=t_window, Jc=Jc,
                window_in_Jc_units=t_window * Jc,
                truncation_drift=float(drift))


def thermal_shape(lam: float, p: int = 4, nmax: int = 140) -> dict:
    """The SHAPE of the antipodal observer's 2-pt at the hologram point.

    Measured verdict (both candidate shapes FAIL, worse as lambda -> 0):
    the shifted 2-pt exceeds unity and grows along the shift --
    |G_obs(0)| ~ 2.3 at lambda = 0.25 -- because the antipodal excursion
    is a FOLD-ADJACENT contour: its e^{+tau(Ej-Ek)} weights are dominated
    by the spectral edge, exactly like the Harlow-Zhao fold
    (ANTISCRAMBLING.md decider #2, fold_edge_probe).  Neither the
    continued sech nor the Gibbs betaJ = pi shape describes it, and the
    fake-period recurrence is absent.  What survives edge domination is
    the detailed-balance RATIO identity (exact at any lambda) and the
    trend-level OTOC flip: shapes do not.  Both failed targets and the
    periodicity deviation are recorded as the measured negatives."""
    from chord import Sector0, q_of_lambda
    Jc = math.sqrt(lam) * math.exp(lam / 8.0)
    beta_f = math.pi / Jc
    tau_s = 0.5 * beta_f
    q = q_of_lambda(lam)
    dpsi = 1.0 / p
    s = Sector0(nmax, q)
    v0 = s.V[0, :]
    M = s.V.T @ ((q ** (dpsi * np.arange(nmax + 1)))[:, None] * s.V)

    def g_obs_at(t):
        return complex((v0 * np.exp((tau_s + 1j * t) * s.E)) @ M
                       @ (v0 * np.exp(-(tau_s + 1j * t) * s.E)))

    ts = np.linspace(0.0, 2.0 / Jc, 25)
    g_obs = np.array([g_obs_at(t) for t in ts])
    # (a) the continued-sech target
    z = Jc * (tau_s + 1j * ts)
    g_sech = np.abs(1.0 / np.cosh(z)) ** (2 * dpsi)
    dev_sech = float(np.max(np.abs(np.abs(g_obs) - g_sech)
                            / np.maximum(g_sech, 1e-12)))
    # (b) the naive Gibbs guess (negative control)
    v_f = v_of_betaJ(math.pi)
    g_th = (np.cos(np.pi * v_f / 2.0)
            / np.cosh(np.pi * v_f * ts / beta_f)) ** (2 * dpsi)
    dev_gibbs = float(np.max(np.abs(np.abs(g_obs) - g_th) / g_th))
    # periodicity with the fake period
    tp = np.linspace(0.0, 0.9 * beta_f, 10)
    per = [abs(g_obs_at(t + beta_f)) / max(abs(g_obs_at(t)), 1e-300)
           for t in tp]
    per_dev = float(np.max(np.abs(np.array(per) - 1.0)))
    return dict(lam=lam, v_f=v_f, beta_f=beta_f,
                ts=[float(x) for x in ts],
                g_obs_abs=[float(abs(x)) for x in g_obs],
                g_sech_pred=[float(x) for x in g_sech],
                g_gibbs_pred=[float(x) for x in g_th],
                max_rel_dev_sech=dev_sech,
                max_rel_dev_gibbs=dev_gibbs,
                fake_periodicity_dev=per_dev)


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
    from chord import beta_chord_from_betaJ
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
                   "beta.J = pi shape prediction, the MSS-evasion "
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
        print(f"shape lam={s_['lam']}: dev vs continued-sech = "
              f"{s_['max_rel_dev_sech']:.3f}  vs Gibbs (neg. control) = "
              f"{s_['max_rel_dev_gibbs']:.3f}  fake-periodicity dev = "
              f"{s_['fake_periodicity_dev']:.3f}")
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
