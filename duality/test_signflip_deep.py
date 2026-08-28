"""Tests for the deep antipodal-map campaign (signflip_deep.json).

Pinned here: the transient anti-scrambling window at the hologram point
grows as lambda falls, which is the semiclassical dS epoch; the
MSS-evasion mechanism, which is phase winding at constant magnitude; the
honest negatives, which are no Gibbs shape, no continued-sech shape and
no fake-period recurrence, because the shifted contour is fold-adjacent
and edge-dominated; and the rate ordering that flags mode contamination.
"""

import json
import math
from pathlib import Path

import pytest

from largeq_anchor import v_of_betaJ

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
    # the quantitative form road_map quotes: the window grows like
    # ~(1/2) ln(1/lambda) (committed pairwise slopes 0.53, 0.40, 0.49) --
    # a power-law growth in 1/lambda would leave this band
    slopes = [(w2 - w1) / math.log(r1["lam"] / r2["lam"])
              for (r1, w1), (r2, w2) in zip(zip(rows, wins),
                                            zip(rows[1:], wins[1:]))]
    assert all(0.35 < s < 0.60 for s in slopes), slopes
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
    """The fitted growth of the shifted curve exceeds BOTH reference rates
    at every lambda: the flipped mode's rate cannot be cleanly extracted (a
    faster scrambling-reasserting structure contaminates the fit) --
    recorded, not hidden.

    Only kappa_obs_time > rate_fake has measured content.  The ordering
    rate_fake > rate_eff is definitional (beta_eff = (1+v) beta_fake, so
    the ratio is exactly 1 + v for any input the producer can generate);
    it is asserted below as the bookkeeping identity it is, not as a
    measurement."""
    for r in data["rate"]:
        assert r["kappa_obs_time"] > r["rate_fake"]          # measured
        assert r["rate_fake"] / r["rate_eff"] == pytest.approx(
            1.0 + v_of_betaJ(r["betaJ"]), rel=1e-9)          # bookkeeping
