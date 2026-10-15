"""Tests for the antipodal sign-flip probe (signflip_probe.json).

Three things are pinned.

THE CLOSED-FORM FLIP IS EXACT.  Re a goes from a(0) to -a(0) at
Delta = pi/v, at every temperature, with a(0) = -1/cos(pi v/2), which
tends to -1 at the T_B = infinity hologram point.  The rotation law holds
to 1e-6.

THE SHIFTED 2-PT FUNCTION stays finite and positive, up to a MEASURED
propagator pole.

IN THE EXACT THEORY THE FLIP IS SEMICLASSICAL.  It is absent for
lambda >= 1.8 and present for lambda <= 1.7, while the unshifted control
scrambles at every lambda.

The probe tests a necessary condition of the candidate map.  It is not a
dictionary.
"""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from antiscrambling import G2_thermal

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def data():
    with open(HERE / "signflip_probe.json") as fh:
        return json.load(fh)


def test_closed_form_flip_exact(data):
    """The antipodal half-fake-period shift flips the OTOC sign EXACTLY, at
    every temperature.  With the gauge-correct phi_12^+ (2026-08-24) the
    unshifted coefficient is real, a(0) = -1/cos(pi v/2), so the rotation
    law a(Delta) = a(0) e^{-i v Delta} makes the flip a clean negation:
    Re a(pi/v) = -a(0) > 0, same magnitude."""
    for cf in data["closed_form"]:
        a0 = cf["re_a_unshifted"]
        assert a0 == pytest.approx(-1.0 / math.cos(math.pi * cf["v"] / 2.0),
                                   rel=1e-9)
        assert a0 < 0                                   # scrambling
        assert cf["re_a_antipodal"] > 0                 # anti-scrambling
        assert cf["re_a_antipodal"] == pytest.approx(-a0, rel=1e-9)
        for row in cf["rows"]:
            assert row["rotation_law_error"] < 1e-6
    # the hologram point (betaJ -> 0, v -> 0) is the -1 -> +1 limit
    coldest = min(data["closed_form"], key=lambda c: c["v"])
    assert coldest["re_a_unshifted"] == pytest.approx(-1.0, abs=0.05)


def test_g2_stays_physical(data):
    """The shifted 2-pt function stays finite and positive all the way to the
    antipodal point, with the propagator pole a margin pi/2 beyond the
    endpoint.

    The pole location is MEASURED from G2_thermal itself (bisection on
    1/|G2|), not read back from the producer's stored constant -- the
    former assertions compared `endpoint_margin` and `pole_at` with the
    hardcoded expressions the producer wrote, so they could not fail."""
    g2 = data["g2"]
    v, betaJ = g2["v"], g2["betaJ"]
    assert g2["all_finite_positive"]
    # recompute the scan independently of the stored values
    vals = [float(np.real(G2_thermal(0.5 * math.pi + d, betaJ, 6)))
            for d in g2["deltas"]]
    assert all(np.isfinite(vals)) and all(x > 0 for x in vals)
    assert vals == pytest.approx(g2["g2"], rel=1e-12)
    # MEASURE the pole.  G2 is exactly real below it and goes complex above
    # (the branch point of base^{2/q}), so bisect on Im G2 == 0.
    endpoint = 0.5 * math.pi + g2["deltas"][-1]
    lo, hi = endpoint, math.pi + math.pi / v + 1.0
    assert complex(G2_thermal(lo, betaJ, 6)).imag == 0.0
    assert complex(G2_thermal(hi, betaJ, 6)).imag != 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if complex(G2_thermal(mid, betaJ, 6)).imag == 0.0:
            lo = mid
        else:
            hi = mid
    assert lo == pytest.approx(math.pi + math.pi / v, rel=1e-12)
    # it really is a divergence, not just a branch point: G2 ~ (theta_pole -
    # theta)^{-2/q}, so the measured log-log slope is 2/q = 1/3 at q = 6
    d1, d2 = 1e-3, 1e-9
    slope = math.log(abs(G2_thermal(lo - d2, betaJ, 6))
                     / abs(G2_thermal(lo - d1, betaJ, 6))) / math.log(d1 / d2)
    assert slope == pytest.approx(2.0 / 6.0, rel=1e-6)
    # and the endpoint margin is that MEASURED pole minus the last point
    assert lo - endpoint == pytest.approx(math.pi / 2.0, rel=1e-8)


def test_chord_flip_is_semiclassical(data):
    """Exact in lambda: unshifted always scrambles; the antipodal shift
    flips the trend for lambda <= 1.7 and not for lambda >= 1.8 (this
    configuration and window)."""
    for ch in data["chord"]:
        assert ch["truncation_drift"] < 1e-9
        rows = {round(r["delta"], 3): r for r in ch["rows"]}
        unshifted = rows[0.0]
        anti = max(ch["rows"], key=lambda r: r["delta"])
        assert unshifted["direction"] == "down"       # always scrambling
        if ch["lam"] <= 1.7:
            assert anti["direction"] == "up", ch["lam"]
        if ch["lam"] >= 1.8:
            assert anti["direction"] == "down", ch["lam"]


def test_antipodal_frame_thermal_at_tomperature(data):
    """Per-line detailed balance of the shifted Wightman weights at
    beta_eff = beta + 2 dtau; with dtau = beta_fake/2 this is
    beta_eff = beta + beta_fake -> the tomperature at the T = infinity
    hologram point.  The line law is exact by construction
    (max_line_deviation is a float-level consistency check); the
    verification with teeth is engine_check_dev: the shifted correlator
    from the contour engine (no Lehmann weights in that code path) must
    match the per-line-shifted sum."""
    for r in data["detailed_balance"]:
        assert r["max_line_deviation"] < 1e-9
        assert r["engine_check_dev"] < 1e-8
    # The hologram-point statement.  beta_eff = (1 + v) beta_fake is
    # BOOKKEEPING (true for every row by construction), so assert it as
    # such -- and get the physics from v -> 0 at the hologram point, which
    # is what makes beta_eff collapse onto beta_fake there.
    for r in data["detailed_balance"]:
        assert r["beta_eff"] == pytest.approx(
            (1.0 + r["v"]) * r["beta_fake_c"], rel=1e-9)
    small = [r for r in data["detailed_balance"] if r["betaJ"] < 0.1]
    assert small
    for r in small:
        assert r["v"] < 0.01                       # the actual content
        assert r["beta_eff"] == pytest.approx(r["beta_fake_c"], rel=0.01)


def test_phase_period_set_by_v_not_kappa(data):
    """The flip works with Delta = pi/v and fails with pi/kappa_eff where
    it works at all -- and both fail at strong coupling.  kappa_eff is a
    window-dependent effective exponent: the same estimator on the large-q
    closed form gives 0.94 > v on the same window, and the exact-chord
    value sits below that (ratio 0.85-0.95), so kappa_eff > v is NOT a
    rate decoupling.  The calibration is recomputed live for one row."""
    import math
    from antiscrambling import otoc_ratio
    from signflip_probe import CONFIG_CROSSED
    for r in data["calibrated"]:
        assert r["kappa"] > r["v_largeq"]
        assert r["kappa_closed_same_window"] > r["kappa"] > r["v_largeq"]
        assert 0.8 < r["kappa_over_closed"] < 1.0
    r = data["calibrated"][0]
    t1, t2, t3, t4 = CONFIG_CROSSED
    rs = [otoc_ratio(t1 + 1j * t, t2 + 1j * t, t3, t4, r["v_largeq"]).real
          for t in (0.5, 1.0, 1.5, 2.0)]
    d1, d2, d3 = rs[1] - rs[0], rs[2] - rs[1], rs[3] - rs[2]
    kc = 0.5 * (math.log(abs(d2 / d1)) + math.log(abs(d3 / d2))) / 0.5
    assert kc == pytest.approx(r["kappa_closed_same_window"], rel=1e-9)
    for r in data["calibrated"]:
        if r["lam"] >= 2.0:
            assert r["largeq"]["direction"] == "down"
            assert r["calibrated"]["direction"] == "down"
        if r["lam"] <= 1.0:
            assert r["largeq"]["direction"] == "up"
            assert r["calibrated"]["direction"] == "down"


def test_toc_bounded_under_shift(data):
    """The time-ordered 4-pt stays bounded under the antipodal shift
    (magnitudes essentially unchanged) while the crossed one flips."""
    for r in data["toc"]:
        assert r["bounded"]
        assert r["max_shifted"] < 2.0 * max(r["max_unshifted"], 0.1)


def test_second_configuration(data):
    """Config robustness at a GENUINELY inequivalent crossed configuration.

    The former CONFIG_CROSSED_B was the primary configuration rigidly
    translated by 0.1 pi; since the corrected otoc_ratio is translation
    invariant, that check carried provably zero information.  The current
    one differs in all three translation invariants (th12-, th34-,
    th12+ - th34+), so its a(0) is genuinely COMPLEX -- which the
    symmetric primary configuration cannot exhibit.

    Two precisions, deliberately separated: the SIGN FLIP (the physics)
    is machine-exact, because the antipodal shift negates a(0) whatever
    its phase; the full complex rotation law is verified to ~4e-5, limited
    by the finite late-time fit window at an asymmetric configuration
    (there is no closed form to compare against here)."""
    cb = data["config_b"]
    assert cb["re_a0"] < 0 < cb["re_a_antipodal"]                 # flips
    assert cb["re_a_antipodal"] == pytest.approx(-cb["re_a0"], rel=1e-12)
    assert cb["rotation_error"] < 1e-3
    assert cb["chord_lam1_unshifted"] == "down"
    assert cb["chord_lam1_antipodal"] == "up"


def test_chord_scan_reproduces_its_committed_row(data):
    """Execute the producer, do not just read its output.

    Every other test here reads signflip_probe.json, which leaves
    signflip_probe.py itself unexercised -- a corruption in the crossed
    contour assembly (the shift construction, the flavor weight, the
    direction logic) would keep the suite green until someone regenerated
    the data.  The engine is deterministic, so the cheapest committed
    chord row reproduces exactly (~0.9 s).  The producer's actual
    arguments are reconstructed from the row itself (main() overrides the
    chord_scan defaults for nmax and thL_grid).

    TOLERANCE.  Cross-platform the last two digits of this chain move, so
    the value comparisons run at rel=1e-9 -- far above the measured
    platform spread, far below any real regression (a 1e-8 perturbation
    of the producer is still caught).  The first curve is a Euclidean
    contour, where the chord correlator is exactly real: its committed
    "im" is float noise about zero (|im| <= 2.3e-15) and reproduces no
    digits at all across platforms.  That one is asserted as a BOUND."""
    import signflip_probe

    row = max(data["chord"], key=lambda r: r["lam"])
    grid = tuple(c["thL"] for c in row["rows"][0]["curve"])
    out = signflip_probe.chord_scan(row["lam"], row["betaJ"], p=row["p"],
                                    nmax=row["nmax"], thL_grid=grid)
    assert out["truncation_drift"] == pytest.approx(
        row["truncation_drift"], rel=1e-6, abs=1e-13)
    for got, want in zip(out["rows"], row["rows"]):
        assert got["direction"] == want["direction"]
        assert got["re_move"] == pytest.approx(want["re_move"],
                                               rel=1e-9, abs=1e-12)
        for gc, wc in zip(got["curve"], want["curve"]):
            assert gc["re"] == pytest.approx(wc["re"], rel=1e-9, abs=1e-12)
            if abs(wc["im"]) < 1e-12:        # analytic zero: bound it
                assert abs(gc["im"]) < 1e-12
            else:
                assert gc["im"] == pytest.approx(wc["im"],
                                                 rel=1e-9, abs=1e-12)


def test_detailed_balance_reproduces_its_committed_row(data):
    """The tomperature-thermality verdict, recomputed from scratch (~0.02 s)
    rather than read back: beta_eff = beta + beta_fake detailed balance and
    the independent-engine cross-check."""
    import signflip_probe

    want = next(r for r in data["detailed_balance"]
                if r["lam"] == 2.0 and r["betaJ"] == 2.0)
    got = signflip_probe.detailed_balance(2.0, 2.0)
    for key in ("beta_c", "dtau_c", "beta_eff", "beta_fake_c"):
        assert got[key] == pytest.approx(want[key], rel=1e-12)
    # the verdict structure itself, not just the stored residuals
    assert got["beta_eff"] == pytest.approx(got["beta_c"]
                                            + got["beta_fake_c"], rel=1e-12)
    assert got["max_line_deviation"] < 1e-12
    assert got["engine_check_dev"] < 1e-12
