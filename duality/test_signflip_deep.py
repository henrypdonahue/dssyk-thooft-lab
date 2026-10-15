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
    # no window may be censored by the grid, and the interpolated maximum
    # sits within one grid step of the grid maximum
    for r in rows:
        assert r["censored"] is False
        assert abs(r["t_window"] - r["t_window_grid"]) <= 0.25
        assert r["ts"][-1] >= 3.5
    # the growth is monotone in ln(1/lambda) but NOT a fixed (1/2)ln(1/lambda):
    # interpolated pairwise slopes 0.66, 0.38, 0.34 (the earlier 0.53/0.40/
    # 0.49 were grid-quantized, one point censored); a power law in
    # 1/lambda would leave this band
    slopes = [(w2 - w1) / math.log(r1["lam"] / r2["lam"])
              for (r1, w1), (r2, w2) in zip(zip(rows, wins),
                                            zip(rows[1:], wins[1:]))]
    assert all(0.2 < s < 0.8 for s in slopes), slopes
    for r in rows:
        assert r["truncation_drift"] < 1e-4
        # the unshifted control scrambles from the start
        un = r["unshifted"]
        assert un[1] < un[0]


def test_shape_is_the_continued_largeq_form(data):
    """For t > 0 the antipodal frame's 2-pt is the continued large-q form
    |1/cos(Jc(tau_s + i t))|^(2 Delta), with a deviation that SHRINKS toward
    semiclassics; the Euclidean 2-pt tracks the fake-period growth toward
    its singularity the same way; the naive Gibbs guess fails and fails
    worse toward semiclassics.  The t = 0 point is the sec(pi/2)
    singularity, cut off by the bounded chord spectrum (|G_obs(0)| > 1).
    The instrument is re-executed at lambda = 1 (1 s)."""
    import math
    import numpy as np
    from signflip_probe import thermal_shape
    rows = sorted(data["shape"], key=lambda r: -r["lam"])
    sech = [r["max_rel_dev_sech"] for r in rows]
    eu = [r["max_rel_dev_euclid"] for r in rows]
    gibbs = [r["max_rel_dev_gibbs"] for r in rows]
    assert all(sech[i + 1] < sech[i] for i in range(len(sech) - 1)), sech
    assert all(eu[i + 1] < eu[i] for i in range(len(eu) - 1)), eu
    assert all(gibbs[i + 1] > gibbs[i] for i in range(len(gibbs) - 1)), gibbs
    assert sech[-1] < 0.25 and eu[-1] < 0.3 and gibbs[-1] > 1.0
    for r in rows:
        assert r["g_obs_abs"][0] > 1.0
        # the stored target is the continued form, recomputed here
        Jc = math.sqrt(r["lam"]) * math.exp(r["lam"] / 8.0)
        z = Jc * (0.5 * r["beta_f"] + 1j * np.array(r["ts"][1:]))
        assert np.allclose(r["g_sech_pred"][1:],
                           np.abs(1.0 / np.cos(z)) ** 0.5, rtol=1e-12)
    live = thermal_shape(1.0)
    row = next(r for r in rows if abs(r["lam"] - 1.0) < 1e-12)
    for key in ("max_rel_dev_sech", "max_rel_dev_gibbs", "max_rel_dev_euclid"):
        assert live[key] == pytest.approx(row[key], rel=1e-9)


def test_mss_evasion_located(data):
    """The growth coefficient's phase winds around the fake circle at
    constant magnitude: no Euclidean placement gives the definite-sign
    strip structure the chaos-bound proof requires.  The verdict booleans
    are recomputed from the stored Re a(y) array and the producer is
    re-executed (0.01 s), not trusted as written."""
    import numpy as np
    from signflip_probe import mss_strip
    m = data["mss"]
    re_a = np.array(m["re_a"])
    assert bool(re_a.max() > 0.5 and re_a.min() < -0.5) == m["oscillates"] is True
    assert m["abs_a_constant"] is True
    live = mss_strip(betaJ=m["betaJ"], n_y=len(m["ys"]))
    assert np.allclose(live["re_a"], m["re_a"], rtol=1e-9, atol=1e-12)
    assert live["v"] == pytest.approx(m["v"], rel=1e-12)
    assert live["oscillates"] and live["abs_a_constant"]


def test_hologram_point_reproduces(data):
    """The lambda = 1 hologram-point row of the deep campaign is re-executed
    and its window and antipodal curve come back to 1e-9."""
    import numpy as np
    from signflip_probe import hologram_point
    row = next(r for r in data["hologram"] if abs(r["lam"] - 1.0) < 1e-12)
    live = hologram_point(1.0)
    assert live["window_in_Jc_units"] == pytest.approx(
        row["window_in_Jc_units"], rel=1e-9)
    assert live["t_window"] == pytest.approx(row["t_window"], rel=1e-9)
    assert np.allclose(live["antipodal"], row["antipodal"], rtol=1e-9, atol=1e-12)
    assert np.allclose(live["unshifted"], row["unshifted"], rtol=1e-9, atol=1e-12)


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
