"""Tests for the condensate / GMOR-corner / vacuum-energy module (rung 21).

Anchors: LM arXiv:2409.11324 (6.1)-(6.3) [transcribed], the Zhitnitsky
chiral value, and this repo's own Jacobi solver in the chiral corner.
"""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from condensate import (CHIRAL_VALUE, I_alpha, condensate_dimensionless,
                        duality_point_value, vacuum_energy_curve)

HERE = Path(__file__).resolve().parent


def test_chiral_value_closed_form():
    """-1/sqrt(12 pi), the Zhitnitsky condensate [LM (6.2)]."""
    assert CHIRAL_VALUE == pytest.approx(-0.16286750396764, rel=1e-10)


def test_I_alpha_independent_route():
    """I(0) via an independent midpoint/tanh-substitution quadrature
    (shares nothing with the scipy.quad splitting)."""
    # alpha = 0: I(0) = (1/2) Int_0^inf (sinh2t - 2t)/(t^3 cosh^2 t) dt
    t1 = np.linspace(1e-8, 10.0, 400001)
    f1 = (np.sinh(2 * t1) - 2 * t1) / (t1 ** 3 * np.cosh(t1) ** 2)
    t2 = np.linspace(10.0, 80.0, 200001)
    f2 = (np.sinh(2 * t2) - 2 * t2) / (t2 ** 3 * np.cosh(t2) ** 2)
    tail = 1.0 / 80.0 ** 2            # Int_80^inf 2/t^3 dt (tanh t = 1)
    val = 0.5 * (np.trapezoid(f1, t1) + np.trapezoid(f2, t2) + tail)
    # trapezoid-limited cross-check of the transcription (the sharper
    # anchors are the chiral-limit and duality-point closed forms)
    assert I_alpha(0.0) == pytest.approx(float(val), rel=2e-5)


def test_condensate_chiral_limit_with_log_rate():
    """f(-1+a) -> CHIRAL_VALUE, with the LM (6.2) log correction:
    1 - f/f_chiral ~= sqrt(3a) log(pi/a)/pi."""
    for a in (1e-3, 1e-4, 1e-5):
        f = condensate_dimensionless(-1.0 + a)
        deficit = 1.0 - f / CHIRAL_VALUE
        predicted = math.sqrt(3.0 * a) * math.log(math.pi / a) / math.pi
        assert deficit == pytest.approx(predicted, rel=0.15)


def test_duality_point_value():
    """f(0) = (log pi - 1 - gamma_E)/(2 pi sqrt(pi)); the alpha-dependent
    terms of (6.1) vanish at alpha = 0."""
    assert duality_point_value() == pytest.approx(-0.03883444, rel=1e-6)
    assert condensate_dimensionless(0.0) == pytest.approx(
        duality_point_value(), rel=1e-10)


def test_vacuum_energy_feynman_hellmann_consistency():
    """d eps/d alpha == f(alpha)/(2 sqrt(pi(1+alpha))) (the curve really is
    the Feynman-Hellmann integral of the condensate)."""
    alphas = np.array([-0.6, -0.599, -0.2, -0.199])
    eps = vacuum_energy_curve(alphas)
    for i in (0, 2):
        d_num = (eps[i + 1] - eps[i]) / (alphas[i + 1] - alphas[i])
        amid = 0.5 * (alphas[i] + alphas[i + 1])
        d_exact = condensate_dimensionless(amid) / (
            2.0 * math.sqrt(math.pi * (1.0 + amid)))
        assert d_num == pytest.approx(d_exact, rel=1e-3)


def test_committed_json_and_gmor_corner():
    """Committed campaign numbers: the vacuum-energy shift at the duality
    point is FINITE and small (the entire 2d counterpart of Sec. 5's
    'corrections'), and the solver's chiral-corner ground state approaches
    the GMOR leading order with the measured sqrt(a) NLO."""
    with open(HERE / "condensate.json") as fh:
        d = json.load(fh)
    assert d["eps_duality_point"] == pytest.approx(-0.0407, abs=0.002)
    assert d["duality_point_f"] == pytest.approx(-0.038834, abs=1e-4)
    corner = d["gmor_corner"]
    ratios = [r["ratio"] for r in sorted(corner, key=lambda r: r["a"])]
    # ratio -> 1 monotonically as a falls
    assert all(ratios[i] < ratios[i + 1] for i in range(len(ratios) - 1))
    for r in corner:
        assert 0.97 < r["nlo_coeff"] < 1.01
    # eps curve is monotone decreasing (condensate is negative throughout)
    eps = d["eps"]
    assert all(eps[i + 1] <= eps[i] + 1e-12 for i in range(len(eps) - 1))
