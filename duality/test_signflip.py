"""Tests for the antipodal sign-flip probe (signflip_probe.json).

What is pinned: the closed-form flip is exact (Re a: -1 -> +1 at
Delta = pi/v, every temperature, rotation law to 1e-6); the shifted 2-pt
stays finite and positive; and in the exact theory the flip is
SEMICLASSICAL -- absent at lambda = 2.0, present for lambda <= 1.6, while
the unshifted control scrambles at every lambda.  The probe is a
necessary-condition test of the candidate map, not a dictionary.
"""

import json
import math
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def data():
    with open(HERE / "signflip_probe.json") as fh:
        return json.load(fh)


def test_closed_form_flip_exact(data):
    for cf in data["closed_form"]:
        assert cf["re_a_unshifted"] == pytest.approx(-1.0, abs=1e-6)
        assert cf["re_a_antipodal"] == pytest.approx(+1.0, abs=1e-6)
        for row in cf["rows"]:
            assert row["rotation_law_error"] < 1e-6


def test_g2_stays_physical(data):
    g2 = data["g2"]
    assert g2["all_finite_positive"]
    assert g2["endpoint_margin"] == pytest.approx(math.pi / 2.0)
    # pole location: pi + pi/v
    assert g2["pole_at"] == pytest.approx(math.pi + math.pi / g2["v"])


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
    # the hologram-point statement: at betaJ -> 0, beta_eff ~ beta_fake
    small = [r for r in data["detailed_balance"] if r["betaJ"] < 0.1]
    assert small
    for r in small:
        assert r["beta_eff"] == pytest.approx(r["beta_fake_c"], rel=0.01)


def test_phase_period_set_by_v_not_kappa(data):
    """The measured magnitude growth kappa exceeds the large-q v, yet the
    flip works with Delta = pi/v and fails with pi/kappa where it works at
    all -- and both fail at strong coupling."""
    for r in data["calibrated"]:
        assert r["kappa"] > r["v_largeq"]
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
    """Config robustness: the rotation law holds at a second crossed
    configuration (Re a -> -Re a to machine precision) and the chord flip
    reproduces at lambda = 1."""
    cb = data["config_b"]
    assert cb["rotation_error"] < 1e-9
    assert cb["re_a_antipodal"] == pytest.approx(-cb["re_a0"], rel=1e-9)
    assert cb["chord_lam1_unshifted"] == "down"
    assert cb["chord_lam1_antipodal"] == "up"
