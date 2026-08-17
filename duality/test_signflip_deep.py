"""Tests for the deep antipodal-map campaign (signflip_deep.json).

What is pinned: the transient anti-scrambling window at the hologram
point grows as lambda falls (the semiclassical dS epoch); the MSS-evasion
mechanism (phase winding at constant magnitude); the honest negatives
(no Gibbs shape, no continued-sech shape, no fake-period recurrence --
the shifted contour is fold-adjacent and edge-dominated); and the rate
ordering that flags mode contamination.
"""

import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def data():
    with open(HERE / "signflip_deep.json") as fh:
        return json.load(fh)


def test_hologram_window_grows_semiclassically(data):
    """At beta = 0 the antipodal frame anti-scrambles in a transient
    window; the window (in Jc units) grows monotonically as lambda
    falls.  Truncation-clean."""
    rows = sorted(data["hologram"], key=lambda r: -r["lam"])
    wins = [r["window_in_Jc_units"] for r in rows]
    assert all(w > 0 for w in wins)
    assert all(wins[i + 1] > wins[i] for i in range(len(wins) - 1))
    for r in rows:
        assert r["truncation_drift"] < 1e-4
        # the unshifted control scrambles from the start
        un = r["unshifted"]
        assert un[1] < un[0]


def test_shape_negatives_edge_dominated(data):
    """Both semiclassical shape targets fail, and worse as lambda falls
    (the fold-adjacency of the shifted contour); the fake-period
    recurrence is absent.  These are the recorded costs of the map."""
    rows = sorted(data["shape"], key=lambda r: -r["lam"])
    sech = [r["max_rel_dev_sech"] for r in rows]
    gibbs = [r["max_rel_dev_gibbs"] for r in rows]
    assert all(sech[i + 1] > sech[i] for i in range(len(sech) - 1))
    assert all(gibbs[i + 1] > gibbs[i] for i in range(len(gibbs) - 1))
    for r in rows:
        assert r["fake_periodicity_dev"] > 0.5
        # the growing-along-the-shift signature: |G_obs(0)| > 1
        assert r["g_obs_abs"][0] > 1.0


def test_mss_evasion_located(data):
    """The growth coefficient's phase winds around the fake circle at
    constant magnitude: no Euclidean placement gives the definite-sign
    strip structure the chaos-bound proof requires."""
    assert data["mss"]["oscillates"]
    assert data["mss"]["abs_a_constant"]


def test_rate_ordering_flags_contamination(data):
    """The fitted growth of the shifted curve exceeds 2 pi/beta_fake >
    2 pi/beta_eff at every lambda: the flipped mode's rate cannot be
    cleanly extracted (a faster scrambling-reasserting structure
    contaminates the fit) -- recorded, not hidden."""
    for r in data["rate"]:
        assert r["kappa_obs_time"] > r["rate_fake"] > r["rate_eff"]
