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
    flips the trend for lambda <= 1.6 and not at lambda = 2.0."""
    for ch in data["chord"]:
        assert ch["truncation_drift"] < 1e-9
        rows = {round(r["delta"], 3): r for r in ch["rows"]}
        unshifted = rows[0.0]
        anti = max(ch["rows"], key=lambda r: r["delta"])
        assert unshifted["direction"] == "down"       # always scrambling
        if ch["lam"] <= 1.6:
            assert anti["direction"] == "up", ch["lam"]
            assert anti["re_move"] > 0.5
        if ch["lam"] >= 2.0:
            assert anti["direction"] == "down", ch["lam"]
