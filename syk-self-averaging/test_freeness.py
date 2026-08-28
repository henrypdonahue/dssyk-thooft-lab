#!/usr/bin/env python3
"""Tests for the freeness-onset measurement (freeness.py).

Fast fresh-computation tests pin the exact t = 0 identities, which are
the classical commuting ones, and the closed-form W1(binomial, arcsine)
constant.

A quick ensemble at N = 12 asserts the measured classical -> free
crossover itself.

JSON-verdict tests assert the production N-trend.  They skip if
freeness.json has not been generated.
"""
import json
from pathlib import Path

import numpy as np
import pytest

import freeness
from freeness import (W1_BINOMIAL_VS_ARCSINE, arcsine_cdf, binomial_cdf,
                      heisenberg, instance_operators, run_point)


def test_w1_constant_closed_form():
    """W1(binomial, arcsine) = 1 - 4(sqrt(2)-1)/pi, vs fine quadrature."""
    grid = np.linspace(-2.001, 2.001, 400001)
    quad = np.trapezoid(np.abs(binomial_cdf(grid) - arcsine_cdf(grid)), grid)
    assert abs(quad - W1_BINOMIAL_VS_ARCSINE) < 1e-5


def test_t0_identities_and_invariants():
    """At t = 0 the pair (E, F) is exactly classically independent:
    m2 = tau[(EF)^2] = 1, w4 = tau[abab] = -1, and odd words vanish.
    spec(E+F) is exactly the binomial law (W1_bin = grid floor, W1_arc =
    the closed-form constant).  Heisenberg evolution preserves tau(F) = 0
    and F(t)^2 = 1 (machine precision)."""
    rng = np.random.default_rng(7)
    lam, ops = instance_operators(12, 4, rng)
    r = freeness.measure_instance(lam, ops, np.array([0.0]))
    assert abs(r["m2"][0] - 1.0) < 1e-12
    assert abs(r["w4"][0] + 1.0) < 1e-12
    assert abs(r["m1"][0]) < 1e-12
    assert abs(r["m3"][0]) < 1e-12
    assert r["W1_bin"][0] < 1e-3
    assert abs(r["W1_arc"][0] - W1_BINOMIAL_VS_ARCSINE) < 1e-3
    # invariants at a generic nonzero time
    Ft = heisenberg(ops[1], lam, 0.37)
    dim = len(lam)
    assert abs(np.trace(Ft)) / dim < 1e-12
    assert np.max(np.abs(Ft @ Ft - np.eye(dim))) < 1e-9


def test_classical_to_free_crossover_quick():
    """This test runs the measurement itself, at quick scale
    (N = 12, 6 instances).  W1_arc decays from the binomial value toward a
    small floor, while W1_bin rises toward W1(binomial, arcsine).  The
    crossover exists at t ~ 1/J.  Every alternating word (free prediction
    0) lands near 0 at late time."""
    r = run_point(N=12, p=4, n_inst=6, seed=901,
                  tJ_grid=(0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0))
    assert r["t_x_J"] is not None and 0.5 < r["t_x_J"] < 3.0
    assert r["W1_arc"][-1] < 0.1
    assert r["W1_bin"][-1] > 0.4
    assert abs(r["m2"][-1]) < 0.1
    assert abs(r["w4"][-1]) < 0.15
    assert max(abs(np.array(r["m3"]))) < 0.15


def _load_json():
    # freeness.json is committed: a missing file is a broken checkout and
    # must FAIL the verdict tests, not skip them
    path = Path(__file__).resolve().parent / "freeness.json"
    return json.loads(path.read_text())


def test_json_t0_identities():
    """Every production point starts exactly classical."""
    data = _load_json()
    for pt in data["points"]:
        assert abs(pt["m2"][0] - 1.0) < 1e-9
        assert abs(pt["w4"][0] + 1.0) < 1e-9
        assert pt["W1_bin"][0] < 1e-3


def test_json_freeness_verdict():
    """The production verdict (N = 12..24, dim 64..4096):

    (a) DEPTH: the late-time W1-to-arcsine plateau falls STRICTLY with
        every step in N, with power ~dim^-0.9.  That power is consistent
        with the ~1/dim floor of a RIGID spectrum converging to the
        arcsine law, and far below the dim^-1/2 of unstructured sampling.
        The free product is exact in the large-N limit by this trend.
    (b) The spectrum actually lands on the arcsine law: W1_bin plateaus at
        W1(binomial, arcsine) itself.  That is, the distance to the
        classical prediction saturates because spec(E+F) IS arcsine.
    (c) ONSET: the crossover sits at the DISSIPATION time t ~ 1.2/J over
        the whole range N = 12..24.  It is NOT N-independent: t_x_J runs
        1.191, 1.172, 1.233, 1.273, and the N = 16 -> 24 rise of +0.101 is
        about 4 sigma given the measured resolution (sigma(t_x_J) ~ 0.011
        at n_inst = 24, ~0.019 at n_inst = 12, from reruns at N = 12).  So
        a drift IS resolved -- but it is non-monotone (it dips at N = 16),
        the lever arm is a factor 2 in N, and ln(24/16) = 0.41, so these
        four points cannot separate a ln N scrambling time from a
        finite-size effect.  Either way the crossover sits at the
        dissipation scale, two decades of N short of where CK's statement
        lives (see README caveats).
    (d) Every alternating word lands on its free-prediction 0 within the
        (statistics-limited) floor."""
    data = _load_json()
    pts = data["points"]
    dims = np.array([pt["dim"] for pt in pts], float)
    plat = np.array([pt["plateau_W1_arc"] for pt in pts])
    assert all(b < a for a, b in zip(plat, plat[1:]))
    slope = np.polyfit(np.log(dims), np.log(plat), 1)[0]
    # two-sided: the README quotes ~dim^-0.9 as consistent with the ~1/dim
    # floor (committed slope -0.93).  A slope steeper than any 1/dim
    # mechanism (say -2) would signal a broken plateau extraction, not a
    # sharper confirmation, so it must fail too.
    assert -1.15 < slope < -0.8, slope
    for pt in pts:
        assert abs(pt["plateau_W1_bin"]
                   - data["w1_binomial_vs_arcsine"]) < 0.01
        assert 1.0 < pt["t_x_J"] < 1.5
        assert pt["plateau_m2"] < 0.005
        assert pt["plateau_w4"] < 0.005
    # the drift is real and non-monotone: pin both, so "N-independent"
    # cannot creep back into the narrative
    tx = [pt["t_x_J"] for pt in pts]
    assert tx[-1] - tx[1] > 0.05                 # resolved rise, N = 16 -> 24
    assert tx[1] < tx[0]                         # and non-monotone at N = 16
    assert max(tx) - min(tx) < 0.15              # still all at t*J ~ 1.2


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
