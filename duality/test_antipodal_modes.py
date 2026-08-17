"""Tests for the antipodal mode separation (antipodal_modes.json).

What is pinned: the sign structure (a flipped positive component plus a
smaller unflipped negative one wherever the flip operates), the
REFUTATION of the harmonic-tower mechanism (kappa_2/kappa_1 far from 2,
comparable harmonic weights), and the recorded fit degeneracies.
"""

import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def data():
    with open(HERE / "antipodal_modes.json") as fh:
        return json.load(fh)


def test_sign_structure_where_flip_operates(data):
    """For lambda in [0.8, 1.5] the constrained two-mode fit finds a
    positive flipped amplitude and a smaller negative unflipped one in
    the antipodal frame."""
    for pt in data["points"]:
        if 0.8 <= pt["lam"] <= 1.5:
            con = pt["antipodal"]["constrained"]
            assert con["A1"] > 0 > con["A2"]
            assert abs(con["A2"]) < abs(con["A1"])
            # and the fit actually describes the curve
            assert pt["antipodal"]["resid_over_range_con"] < 0.05


def test_lambda06_fit_degeneracy_recorded(data):
    """The lambda = 0.6 two-mode fit collapses (amplitudes blow up with
    swapped signs) -- kept as the recorded instrument failure, not
    silently dropped."""
    pt = next(p for p in data["points"] if p["lam"] == 0.6)
    con = pt["antipodal"]["constrained"]
    assert not (con["A1"] > 0 > con["A2"]) or abs(con["A1"]) > 10


def test_harmonic_tower_refuted(data):
    """The sharp signature kappa_2 = 2 kappa_1 fails: measured ratios sit
    near 1, and the harmonic weights are comparable -- the shift-response
    is not a v-harmonic tower at finite lambda."""
    checked = 0
    for h in data["harmonics"]:
        if h["kappa_1"] > 0.05:          # rate measurable on this window
            ratio = h["kappa_2"] / h["kappa_1"]
            assert ratio < 1.6, (h["lam"], ratio)
            checked += 1
        assert 0.3 < h["c2_over_c1_last"] < 3.0
    assert checked >= 1
