"""Tests for antiscrambling.py.

Four things: the closed-form identities, the scrambling-sign verdict, the
Harlow-Zhao fold map, and the ED ground truth in
syk-self-averaging/fold_bench.json.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from antiscrambling import (G2_thermal, growth_coefficient,
                            growth_coefficient_numeric, fold_map,
                            fold_map_negT, f_toc, g_otoc, lyapunov_rate,
                            otoc_ratio, prefactor_C, squeeze, toc_ratio)
from largeq_anchor import g_thermal_euclidean, v_of_betaJ

_BENCH = Path(__file__).resolve().parent.parent / "syk-self-averaging" / \
    "fold_bench.json"


# ---------------------------------------------------------------------------
# closed-form identities
# ---------------------------------------------------------------------------
def test_prefactor_small_v():
    """C(v) ~ (pi v)^2 / 2 as v -> 0."""
    for v in (1e-3, 1e-4):
        assert prefactor_C(v) == pytest.approx((np.pi * v) ** 2 / 2, rel=1e-4)


def test_squeeze_domain():
    """theta in (0, 2pi) maps to phi in pi(1-v, 1+v)."""
    v = 0.37
    assert squeeze(0.0, v) == pytest.approx(np.pi * (1 - v))
    assert squeeze(2 * np.pi, v) == pytest.approx(np.pi * (1 + v))
    assert squeeze(np.pi, v) == pytest.approx(np.pi)


def test_f_and_g_v_to_zero_limits():
    """v -> 0: f(theta) -> theta/pi (the cot(pi v/2) ~ 2/(pi v) pole eats
    the tan ~ v zero), and g stays finite; all connected structure is then
    carried by the C(v) ~ v^2 prefactor."""
    for th in (0.4, 1.3, 2.9):
        assert complex(f_toc(th, 1e-8)).real == pytest.approx(th / np.pi,
                                                              rel=1e-5)
    assert abs(complex(g_otoc(squeeze(1.3, 1e-8)))) < 10.0


def test_G2_thermal_matches_largeq_anchor_on_real_domain():
    """The complex-safe G2 must equal g_thermal_euclidean^{1/q} at real
    theta in (0, 2pi) -- this keeps largeq_anchor's e^{g} exercised."""
    q, betaJ = 6.0, 2.0
    beta_eff = 2 * np.pi
    for th in (0.3, 1.0, np.pi, 5.0):
        ref = g_thermal_euclidean(th, beta_eff, betaJ / beta_eff) ** (1 / q)
        assert complex(G2_thermal(th, betaJ, q)).real == \
            pytest.approx(ref, rel=1e-12)
        assert abs(complex(G2_thermal(th, betaJ, q)).imag) < 1e-12


def test_lyapunov_rate_limits():
    """lambda_L/scriptJ -> 2 at T_B = infinity (betaJ -> 0): the growth
    survives at infinite temperature."""
    assert lyapunov_rate(0.0) == 2.0
    assert lyapunov_rate(1e-3) == pytest.approx(2.0, rel=1e-4)
    # conformal limit: lambda_L -> 2 pi / beta (maximal chaos)
    assert lyapunov_rate(1e4) * 1e4 / (2 * np.pi) == pytest.approx(1.0,
                                                                   rel=1e-3)


# ---------------------------------------------------------------------------
# the sign verdict (this task, closed-form side)
# ---------------------------------------------------------------------------
def test_growth_coefficient_is_real_and_negative():
    """a(v) = -1/cos(pi v/2): REAL and negative at every temperature -- the
    scrambling sign (CK's c > 0), unsuppressed at v -> 0 (T_B = infinity),
    where it tends to -1.  |a| grows with beta."""
    prev = 0.0
    for v in (0.05, 0.3, 0.6, 0.9):
        a = growth_coefficient(v)
        assert a.imag == pytest.approx(0.0, abs=1e-15)
        assert a.real == pytest.approx(-1.0 / np.cos(np.pi * v / 2.0),
                                       rel=1e-12)
        assert a.real < 0.0                      # the scrambling verdict
        assert abs(a.real) > prev                # |a| grows with v (beta)
        prev = abs(a.real)
    # the T_B = infinity hologram point: a -> -1, growth NOT suppressed
    assert growth_coefficient(1e-6).real == pytest.approx(-1.0, abs=1e-9)


@pytest.mark.parametrize("betaJ", [2.0, 0.5, 8.0])
def test_growth_coefficient_numeric_agrees(betaJ):
    """The late-time fit of otoc_ratio recovers a(v) = -1/cos(pi v/2)
    (validates the asymptotics of the transcription, not just the
    formula)."""
    v = v_of_betaJ(betaJ)
    assert growth_coefficient_numeric(betaJ) == pytest.approx(
        growth_coefficient(v).real, rel=5e-3)


# ---------------------------------------------------------------------------
# the gauge that makes the transcription correct (2026-08-24 review)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("v", [0.15, 0.4705, 0.75])
def test_otoc_is_time_translation_invariant(v):
    """script-F/F_d is a ratio of thermal correlators, so it depends only on
    time DIFFERENCES -- as the chord engine (chord_fold._segments) and ED
    (fold_bench) both manifestly do.  This pins Streicher's phi_12^+ to the
    gauge-invariant squeeze(th12+ - th34+): the gauge-fixed display
    squeeze(th12+), valid only at theta_34^+ = 0, drifts by a factor ~2
    over a 0.5 shift."""
    cfg = (1.5 * np.pi, 0.5 * np.pi, np.pi, 0.0)
    ref = complex(otoc_ratio(*cfg, v))
    for c in (-0.2, 0.1, 0.3, 0.5):
        shifted = complex(otoc_ratio(*(t + c for t in cfg), v))
        assert shifted == pytest.approx(ref, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("v", [0.2, 0.4705, 0.8])
def test_otoc_minus_toc_is_minus_two_at_coincidence(v):
    """Structural identity: at theta_2 = theta_3 the crossed and uncrossed
    branches differ by exactly -2 (the exchange term at coincidence).  The
    gauge-fixed form gives arbitrary values here; this is a second,
    configuration-by-configuration check of the corrected phi_12^+."""
    rng = np.random.default_rng(11)
    for _ in range(6):
        a, b, d = np.sort(rng.uniform(0.15, 2 * np.pi - 0.15, 3))[::-1]
        cfg = (a, b, b, d)
        diff = complex(otoc_ratio(*cfg, v)) - complex(toc_ratio(*cfg, v))
        assert diff.real == pytest.approx(-2.0, abs=1e-9)


# ---------------------------------------------------------------------------
# the HZ fold map
# ---------------------------------------------------------------------------
def test_fold_map_implements_HZ_610():
    """(6.10): tau31 and tau42 map to their KMS reflections, and
    tau43 + tau21 maps to 4 pi n - (tau43 + tau21) -- the w -> -w flip."""
    taus = (0.3, 1.1, 2.0, 2.9)
    for n in (1, 2):
        q = fold_map(taus, n)
        assert q[2] - q[0] == pytest.approx(2 * np.pi - (taus[2] - taus[0]))
        assert q[3] - q[1] == pytest.approx(2 * np.pi - (taus[3] - taus[1]))
        s = (taus[3] - taus[2]) + (taus[1] - taus[0])
        assert (q[3] - q[2]) + (q[1] - q[0]) == \
            pytest.approx(4 * np.pi * n - s)


def test_fold_negT_is_reflection_of_fold():
    """(6.12) = -(6.9): the negative-temperature variant is the H -> -H
    reflection of the fold."""
    taus = (0.3, 1.1, 2.0, 2.9)
    a = fold_map(taus, 1)
    b = fold_map_negT(taus, 1)
    for x, y in zip(a, b):
        assert x == pytest.approx(-y)


# ---------------------------------------------------------------------------
# ED ground truth (fold_bench.json)
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def bench():
    # fold_bench.json is committed: missing = broken checkout, fail not skip
    return json.loads(_BENCH.read_text())["points"]


def test_bench_fold_thetas_match_fold_map(bench):
    """The ED campaign's fold configuration is exactly fold_map of its
    crossed configuration (ties the two modules together)."""
    for pt in bench:
        expect = fold_map(pt["config_crossed"], pt["fold_n"])
        assert np.allclose(pt["fold_thetas"], expect)


def test_bench_g2_calibration_floor(bench):
    """Standard 2-pt: ED deviation-from-1 has the same SIGN as the closed
    form and sits at the O(lambda, 1/p) floor (ratio in (0.2, 0.95));
    at fixed (p, betaJ) the ratio INCREASES with N (toward the closed
    form)."""
    ratios_by_N = {}
    for pt in bench:
        for th, g in zip(pt["th_g2"], pt["g2_mean"]):
            cf = complex(G2_thermal(th, pt["betaJ"], pt["p"])).real
            r = (g - 1) / (cf - 1)
            assert 0.2 < r < 0.95, (pt["N"], pt["betaJ"], th, r)
        if pt["betaJ"] == 2.0:
            th, g = pt["th_g2"][2], pt["g2_mean"][2]
            cf = complex(G2_thermal(th, pt["betaJ"], pt["p"])).real
            ratios_by_N[pt["N"]] = (g - 1) / (cf - 1)
    Ns = sorted(ratios_by_N)
    assert len(Ns) >= 3
    assert all(ratios_by_N[a] < ratios_by_N[b]
               for a, b in zip(Ns, Ns[1:]))


def test_bench_ed_kms_reflection(bench):
    """G(2 pi - theta) = G(theta) exactly in ED (the th_g2 grid contains the
    pi/2, 3pi/2 KMS pair)."""
    for pt in bench:
        g = dict(zip(pt["th_g2"], pt["g2_mean"]))
        assert g[0.5 * np.pi] == pytest.approx(g[1.5 * np.pi], rel=1e-10)


def test_bench_2pt_fold_grows_above_one(bench):
    """The folded 2-pt (e^{+tau H} evolution, HZ fig 7c) exceeds 1 -- the
    qualitative signature of the continued closed form -- and is finite
    (boundedness of the SYK spectrum, HZ's requirement)."""
    for pt in bench:
        for g in pt["g2_fold_mean"]:
            assert 1.0 < g < 3.0


def test_bench_4pt_signs_match_closed_forms(bench):
    """TOC positive, crossed OTOC negative -- ED agrees with the closed-form
    signs at every UNFOLDED point, and now quantitatively (see
    test_bench_unfolded_is_quantitative).  The FOLDED closed form is
    deliberately not sign-asserted here: the folded angles lie far outside
    Streicher's principal window, so its branch is a choice, not a
    prediction (test_fold_angles_leave_principal_window)."""
    for pt in bench:
        v = v_of_betaJ(pt["betaJ"])
        assert pt["toc"]["re"] > 0
        assert complex(toc_ratio(*pt["config_toc"], v)).real > 0
        assert pt["otoc"]["re"] < 0
        assert complex(otoc_ratio(*pt["config_crossed"], v)).real < 0
        assert pt["fold"]["re"] < 0                       # ED, unambiguous


def test_bench_unfolded_is_quantitative(bench):
    """With the gauge-correct phi_12^+ the crossed 4-pt closed form is not
    just a sign/trend guide: ED/closed is 1.02-1.09 at every bench point,
    and at fixed betaJ = 2 it converges TOWARD 1 as N grows (1.084, 1.064,
    1.048 at N = 12, 14, 16).  Before the 2026-08-24 gauge fix the same
    ratios read 1.80-1.84 and were attributed to an O(1/p) floor."""
    for pt in bench:
        v = v_of_betaJ(pt["betaJ"])
        closed = complex(otoc_ratio(*pt["config_crossed"], v)).real
        assert 1.0 < pt["otoc"]["re"] / closed < 1.10
    b2 = sorted((pt for pt in bench if pt["betaJ"] == 2.0),
                key=lambda q: q["N"])
    v = v_of_betaJ(2.0)
    rs = [pt["otoc"]["re"] / complex(otoc_ratio(*pt["config_crossed"], v)).real
          for pt in b2]
    assert all(b < a for a, b in zip(rs, rs[1:]))     # converging toward 1
    assert rs[-1] < 1.05


def test_fold_angles_leave_principal_window():
    """The HZ fold sends the configuration far outside the window where
    Streicher's expression is the principal branch: theta_12^- and
    theta_34^- go to -3 pi and (theta_12^+ - theta_34^+) to -5 pi, against a
    principal window of (0, 2 pi).  So 'the continued semiclassical answer'
    at the fold is a BRANCH CHOICE.  This test pins that fact so the fold
    verdicts are read against the closed-form-FREE connected fractions
    (tests below), not against a branch."""
    cfg = (1.5 * np.pi, 0.5 * np.pi, np.pi, 0.0)
    for angles, inside in ((cfg, True), (fold_map(cfg, 1), False)):
        t12m = angles[0] - angles[1]
        t34m = angles[2] - angles[3]
        tp = angles[0] + angles[1] - angles[2] - angles[3]
        ok = all(0.0 < x < 2 * np.pi for x in (t12m, t34m, tp))
        assert ok is inside
    fb = fold_map(cfg, 1)
    assert fb[0] - fb[1] == pytest.approx(-3 * np.pi)
    assert fb[0] + fb[1] - fb[2] - fb[3] == pytest.approx(-5 * np.pi)


def test_bench_scrambling_sign_verdict(bench):
    """Model level: Re script-F/F_d decreases monotonically with
    Lorentzian time at EVERY (N, betaJ) point -- finite-N SYK scrambles
    (CK's c > 0, the AdS/time-delay sign), at betaJ = 2 and at the
    near-T_B-infinity point betaJ = 0.5.  A dS observer requires the
    opposite sign."""
    for pt in bench:
        seq = [pt["otoc"]["re"]] + [pt[f"growth_thL={t}"]["re"]
                                    for t in pt["thL_grid"][1:]]
        assert all(b < a for a, b in zip(seq, seq[1:])), (pt["N"],
                                                          pt["betaJ"], seq)


def test_bench_fold_lorentzian_response_tracks_continuation(bench):
    """the folded correlator's Lorentzian response
    Im script-F/F_d is nonzero with the SAME SIGN and monotone trend as the
    analytically continued closed form (both negative, growing in
    magnitude) at every point.  The unfolded ED correlator has Im = 0 by
    the symmetry of the configuration -- the fold breaks it exactly the
    way the continuation predicts."""
    for pt in bench:
        v = v_of_betaJ(pt["betaJ"])
        ims_ed, ims_cf = [], []
        for t in pt["thL_grid"][1:]:
            ims_ed.append(pt[f"fold_thL={t}"]["im"])
            th = pt["fold_thetas"]
            ims_cf.append(complex(otoc_ratio(th[0] + 1j * t, th[1] + 1j * t,
                                             th[2], th[3], v)).imag)
            assert abs(pt[f"growth_thL={t}"]["im"]) < 5e-3
        # ED side: unambiguous, and the whole robust content of this test.
        assert all(x < 0 for x in ims_ed)
        assert all(b < a for a, b in zip(ims_ed, ims_ed[1:]))
        # Continued closed form: at the folded angles its BRANCH (hence its
        # overall sign) is a choice -- see
        # test_fold_angles_leave_principal_window -- so only the MAGNITUDE
        # is a prediction.  It tracks ED well and better as N grows: at
        # betaJ = 2 the ratio |ED/closed| runs 0.69-0.73 (N=12), 0.86-0.91
        # (N=14), 1.02-1.07 (N=16).
        assert all(x != 0 for x in ims_cf)
        assert all(abs(b) > abs(a) for a, b in zip(ims_cf, ims_cf[1:]))
        if pt["betaJ"] == 2.0 and pt["N"] == 16:
            for e, c in zip(ims_ed, ims_cf):
                assert abs(e / c) == pytest.approx(1.0, abs=0.10)


# ---------------------------------------------------------------------------
# the fixed-lambda scan (fold_scan.json)
# ---------------------------------------------------------------------------
_SCAN = Path(__file__).resolve().parent.parent / "syk-self-averaging" / \
    "fold_scan.json"


@pytest.fixture(scope="module")
def scan():
    # fold_scan.json is committed: missing = broken checkout, fail not skip
    return json.loads(_SCAN.read_text())["points"]


def _fold_ratio(pt):
    v = v_of_betaJ(pt["betaJ"])
    cf = complex(otoc_ratio(*pt["fold_thetas"], v)).real
    return pt["fold"]["re"] / cf, abs(pt["fold"]["err_re"] / cf)


def test_scan_verdict_N_obstruction_not_lambda(scan):
    """the folded Euclidean drift is a genuine
    N-obstruction, not a lambda-correction.

    (a) PRIMARY, closed-form-free: the folded connected fraction
        |(F-F_d)/F_d| = |ratio|/N grows monotonically and reaches O(1) by
        N = 20 (0.365 -> 0.786), while lambda = 16/N SHRINKS: the 1/N
        hierarchy collapses under the fold, opposite to a correction that
        vanishes toward the semiclassical regime.  This needs no closed
        form and is the verdict;
    (b) CONTROL: the unfolded crossed OTOC ratio is N-stable across the
        whole p = 4 series (spread < 0.03) -- the fold is the culprit, not
        the closed form;
    (c) SECONDARY: the |folded ratio| against the continued form grows
        strictly with N by more than 5x.  Magnitudes only: the folded
        branch is a choice (test_fold_angles_leave_principal_window);
    (d) matched lambda = 2.0: the larger-N point (p=6, N=18) drifts
        further than (p=4, N=8) by many sigma, despite smaller 1/p error.
    """
    p4 = sorted((pt for pt in scan if pt["p"] == 4), key=lambda q: q["N"])
    conn = [abs(pt["fold"]["re"]) / pt["N"] for pt in p4]
    assert all(b > a for a, b in zip(conn, conn[1:]))                # (a)
    assert conn[0] < 0.4 and conn[-1] > 0.7
    v = v_of_betaJ(2.0)
    cf_otoc = {pt["N"]: complex(otoc_ratio(*pt["config_crossed"], v)).real
               for pt in p4}
    unfolded = [pt["otoc"]["re"] / cf_otoc[pt["N"]] for pt in p4]
    assert max(unfolded) - min(unfolded) < 0.03                      # (b)
    folded = [abs(_fold_ratio(pt)[0]) for pt in p4]
    assert all(b > a for a, b in zip(folded, folded[1:]))            # (c)
    assert folded[-1] > 5 * folded[0]
    anchor6 = next(pt for pt in scan if pt["p"] == 6 and pt["N"] == 18)
    anchor4 = next(pt for pt in scan if pt["p"] == 4 and pt["N"] == 8)
    r6, e6 = _fold_ratio(anchor6)
    r4, e4 = _fold_ratio(anchor4)
    assert abs(r6) - abs(r4) > 4 * np.hypot(e6, e4)                  # (d)


def test_chord_engine_converges_onto_invariant_form():
    """ANTISCRAMBLING.md's third piece of evidence for the invariant form,
    recomputed from the committed chord_fold.json rather than quoted: the
    exact N = infinity chord engine's crossed 4-pt over the closed form
    runs 0.964 -> 0.979 over lambda = 2 -> 0.5 (0.973 at lambda = 4), and
    the time-ordered control 0.683 -> 0.945 over the same leg -- both
    heading to 1 as lambda -> 0, where the double-scaled theory is the
    large-q limit."""
    with open(Path(__file__).resolve().parent / "chord_fold.json") as fh:
        d = json.load(fh)
    cc, ct = d["config_crossed"], d["config_toc"]
    ratios = {}
    for pt in d["points"]:
        v = pt["v"]
        ot = pt["otoc"]["re"] if isinstance(pt["otoc"], dict) else pt["otoc"]
        tc = pt["toc"]["re"] if isinstance(pt["toc"], dict) else pt["toc"]
        ratios[round(pt["lam"], 6)] = (
            ot / complex(otoc_ratio(*cc, v)).real,
            tc / complex(toc_ratio(*ct, v)).real)
    assert ratios[4.0][0] == pytest.approx(0.973, abs=1e-3)
    assert ratios[2.0][0] == pytest.approx(0.964, abs=1e-3)
    assert ratios[0.5][0] == pytest.approx(0.979, abs=1e-3)
    assert ratios[2.0][1] == pytest.approx(0.683, abs=1e-3)
    assert ratios[0.5][1] == pytest.approx(0.945, abs=1e-3)
    # every point within 4% of 1 for the crossed ratio; the TOC control
    # rises monotonically along lambda = 2 -> 0.5
    assert all(0.96 < r[0] < 1.0 for r in ratios.values())
    leg = [ratios[l][1] for l in sorted(ratios, reverse=True) if l <= 2.0]
    assert all(b > a for a, b in zip(leg, leg[1:]))
