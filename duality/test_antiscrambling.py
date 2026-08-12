"""Tests for antiscrambling.py (rungs 25-26): closed-form identities, the
scrambling-sign verdict, the HZ fold map, and the ED ground truth in
syk-self-averaging/fold_bench.json."""

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
# the sign verdict (rung 25, closed-form side)
# ---------------------------------------------------------------------------
def test_growth_coefficient_real_part_is_minus_one():
    """Re a(v) = -1 EXACTLY at every temperature: the scrambling sign (CK's
    c > 0), unsuppressed at v -> 0 (T_B = infinity)."""
    for v in (0.05, 0.3, 0.6, 0.9):
        assert growth_coefficient(v).real == pytest.approx(-1.0, abs=1e-12)


@pytest.mark.parametrize("betaJ", [2.0, 0.5, 8.0])
def test_growth_coefficient_numeric_agrees(betaJ):
    """The late-time fit of otoc_ratio recovers Re a = -1 (validates the
    asymptotics of the transcription, not just the formula)."""
    assert growth_coefficient_numeric(betaJ) == pytest.approx(-1.0, rel=5e-3)


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
    if not _BENCH.exists():
        pytest.skip("run syk-self-averaging/fold_bench.py first")
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
        v = v_of_betaJ(pt["betaJ"])
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
    """TOC positive, crossed OTOC negative, folded crossed negative -- ED
    agrees with the closed-form signs at every point."""
    for pt in bench:
        v = v_of_betaJ(pt["betaJ"])
        assert pt["toc"]["re"] > 0
        assert complex(toc_ratio(*pt["config_toc"], v)).real > 0
        assert pt["otoc"]["re"] < 0
        assert complex(otoc_ratio(*pt["config_crossed"], v)).real < 0
        assert pt["fold"]["re"] < 0
        assert complex(otoc_ratio(*pt["fold_thetas"], v)).real < 0


def test_bench_scrambling_sign_verdict(bench):
    """RUNG 25, model level: Re script-F/F_d decreases monotonically with
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
    """RUNG 26 (first data): the folded correlator's Lorentzian response
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
        assert all(x < 0 for x in ims_ed)
        assert all(x < 0 for x in ims_cf)
        assert all(b < a for a, b in zip(ims_ed, ims_ed[1:]))
        # quantitative tracking at the smallest lambda-correction points:
        # betaJ = 2, N = 12 tracks within 10%; document-level, not asserted
        # tighter because the Euclidean N-trend at fixed p is a genuine
        # obstruction (decided by the fixed-lambda scan below; see
        # ANTISCRAMBLING.md).
        if pt["betaJ"] == 2.0 and pt["N"] == 12:
            for e, c in zip(ims_ed, ims_cf):
                assert e / c == pytest.approx(1.0, abs=0.15)


# ---------------------------------------------------------------------------
# the fixed-lambda scan (fold_scan.json): rung-26 decider #1
# ---------------------------------------------------------------------------
_SCAN = Path(__file__).resolve().parent.parent / "syk-self-averaging" / \
    "fold_scan.json"


@pytest.fixture(scope="module")
def scan():
    if not _SCAN.exists():
        pytest.skip("run syk-self-averaging/fold_bench.py --scan first")
    return json.loads(_SCAN.read_text())["points"]


def _fold_ratio(pt):
    v = v_of_betaJ(pt["betaJ"])
    cf = complex(otoc_ratio(*pt["fold_thetas"], v)).real
    return pt["fold"]["re"] / cf, abs(pt["fold"]["err_re"] / cf)


def test_scan_verdict_N_obstruction_not_lambda(scan):
    """RUNG 26, decider #1: the folded Euclidean drift is a genuine
    N-obstruction, not a lambda-correction.

    (a) CONTROL: the unfolded crossed OTOC ratio is N-stable across the
        whole p = 4 series (spread < 0.06 around its O(1/p) offset) --
        the fold is the culprit, not the closed form;
    (b) the folded ratio grows STRICTLY with N at fixed p = 4 while
        lambda = 16/N shrinks -- opposite to a correction that vanishes
        toward the semiclassical regime;
    (c) matched lambda = 2.0: the larger-N point (p=6, N=18) drifts
        further than (p=4, N=8) by many sigma, despite smaller 1/p error;
    (d) the folded connected fraction |(F-F_d)/F_d| = |ratio|/N grows
        monotonically and reaches O(1) by N = 20: the 1/N hierarchy
        collapses under the fold (HZ's 'continuation commutes' fails at
        ED sizes)."""
    p4 = sorted((pt for pt in scan if pt["p"] == 4), key=lambda q: q["N"])
    v = v_of_betaJ(2.0)
    cf_otoc = {pt["N"]: complex(otoc_ratio(*pt["config_crossed"], v)).real
               for pt in p4}
    unfolded = [pt["otoc"]["re"] / cf_otoc[pt["N"]] for pt in p4]
    assert max(unfolded) - min(unfolded) < 0.06                      # (a)
    folded = [_fold_ratio(pt)[0] for pt in p4]
    assert all(b > a for a, b in zip(folded, folded[1:]))            # (b)
    assert folded[-1] > 5 * folded[0]
    anchor6 = next(pt for pt in scan if pt["p"] == 6 and pt["N"] == 18)
    anchor4 = next(pt for pt in scan if pt["p"] == 4 and pt["N"] == 8)
    r6, e6 = _fold_ratio(anchor6)
    r4, e4 = _fold_ratio(anchor4)
    assert r6 - r4 > 4 * np.hypot(e6, e4)                            # (c)
    conn = [abs(pt["fold"]["re"]) / pt["N"] for pt in p4]
    assert all(b > a for a, b in zip(conn, conn[1:]))                # (d)
    assert conn[-1] > 0.7
