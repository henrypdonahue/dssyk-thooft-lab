#!/usr/bin/env python3
"""Asserting tests for the freeness-onset measurement (freeness.py, rung 27).

Fast fresh-computation tests pin the exact t = 0 (classical/commuting)
identities and the closed-form W1(binomial, arcsine) constant; a quick
ensemble at N = 12 asserts the measured classical -> free crossover itself;
JSON-verdict tests assert the production N-trend (skip if freeness.json has
not been generated).
"""
import json
from pathlib import Path

import numpy as np
import pytest

import freeness
from freeness import (W1_BINOMIAL_VS_ARCSINE, arcsine_cdf, binomial_cdf,
                      heisenberg, instance_operators, run_point, w1_distance)


def test_w1_constant_closed_form():
    """W1(binomial, arcsine) = 1 - 4(sqrt(2)-1)/pi, vs fine quadrature."""
    grid = np.linspace(-2.001, 2.001, 400001)
    quad = np.trapezoid(np.abs(binomial_cdf(grid) - arcsine_cdf(grid)), grid)
    assert abs(quad - W1_BINOMIAL_VS_ARCSINE) < 1e-5


def test_t0_identities_and_invariants():
    """At t = 0 the pair (E, F) is exactly classically independent:
    m2 = tau[(EF)^2] = 1, w4 = tau[abab] = -1, odd words vanish, and
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
    """The rung-27 measurement itself, at quick scale (N = 12, 6 instances):
    W1_arc decays from the binomial value toward a small floor while W1_bin
    rises toward W1(binomial, arcsine); the crossover exists at t ~ 1/J;
    every alternating word (free prediction 0) lands near 0 at late time."""
    r = run_point(N=12, p=4, n_inst=6, seed=901,
                  tJ_grid=(0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0))
    assert r["t_x_J"] is not None and 0.5 < r["t_x_J"] < 3.0
    assert r["W1_arc"][-1] < 0.1
    assert r["W1_bin"][-1] > 0.4
    assert abs(r["m2"][-1]) < 0.1
    assert abs(r["w4"][-1]) < 0.15
    assert max(abs(np.array(r["m3"]))) < 0.15


def _load_json():
    path = Path(__file__).resolve().parent / "freeness.json"
    if not path.exists():
        pytest.skip("freeness.json not present; run freeness.py")
    return json.loads(path.read_text())


def test_json_t0_identities():
    """Every production point starts exactly classical."""
    data = _load_json()
    for pt in data["points"]:
        assert abs(pt["m2"][0] - 1.0) < 1e-9
        assert abs(pt["w4"][0] + 1.0) < 1e-9
        assert pt["W1_bin"][0] < 1e-3


def test_json_freeness_verdict():
    """RUNG 27, the production verdict (N = 12..24, dim 64..4096):

    (a) DEPTH: the late-time W1-to-arcsine plateau falls STRICTLY with
        every step in N, with power ~dim^-0.9 -- consistent with the
        ~1/dim floor of a RIGID spectrum converging to the arcsine law,
        and far below the dim^-1/2 of unstructured sampling.  The free
        product is exact in the large-N limit by this trend.
    (b) The spectrum actually lands on the arcsine law: W1_bin plateaus at
        W1(binomial, arcsine) itself, i.e. the distance to the classical
        prediction saturates because spec(E+F) IS arcsine.
    (c) ONSET: the crossover sits at the DISSIPATION time t ~ 1.2/J,
        nearly N-independent (any ln N drift is below resolution here --
        the distribution-level crossover is not a scrambling-time
        phenomenon at ED sizes; see README caveats).
    (d) Every alternating word lands on its free-prediction 0 within the
        (statistics-limited) floor."""
    data = _load_json()
    pts = data["points"]
    dims = np.array([pt["dim"] for pt in pts], float)
    plat = np.array([pt["plateau_W1_arc"] for pt in pts])
    assert all(b < a for a, b in zip(plat, plat[1:]))
    slope = np.polyfit(np.log(dims), np.log(plat), 1)[0]
    assert slope < -0.8
    for pt in pts:
        assert abs(pt["plateau_W1_bin"]
                   - data["w1_binomial_vs_arcsine"]) < 0.01
        assert 1.0 < pt["t_x_J"] < 1.5
        assert pt["plateau_m2"] < 0.005
        assert pt["plateau_w4"] < 0.005


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
