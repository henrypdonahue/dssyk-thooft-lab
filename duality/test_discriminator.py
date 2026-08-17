"""Tests for the KMS-axis discriminator note (rung 24).

The note's repo-side claims, pinned: the fake-temperature identity
beta_fake * lambda_L = 2 pi across the repo's v machinery (tying the
sine-dilaton (1.7) to rung 25), and the fold verdicts it quotes from the
committed JSONs.
"""

import json
import math
from pathlib import Path

import pytest

from antiscrambling import growth_coefficient, lyapunov_rate
from largeq_anchor import v_of_betaJ

HERE = Path(__file__).resolve().parent


@pytest.mark.parametrize("betaJ", [0.5, 2.0, 10.0])
def test_fake_temperature_identity(betaJ):
    """Internal consistency of the repo's own v machinery: with beta_fake
    DEFINED as beta/v, lambda_L = 2 pi v/beta makes beta_fake*lambda_L =
    2 pi an algebraic identity of the large-q parametrization -- which is
    what this test pins (post-review honesty: it is NOT an independent
    test of sine-dilaton's beta_BH = 2 pi/sin(theta); that contact is a
    units-level identification, and a naive chord-saddle check
    sin(theta*) =? v does NOT close under the beta_chord bridge --
    an open item, stated in DISCRIMINATOR.md)."""
    v = v_of_betaJ(betaJ)
    beta = betaJ                       # in units scriptJ = 1
    beta_fake = beta / v
    assert lyapunov_rate(betaJ) * beta == pytest.approx(2 * math.pi * v,
                                                        rel=1e-10)
    assert beta_fake * lyapunov_rate(betaJ) == pytest.approx(2 * math.pi,
                                                             rel=1e-10)


def test_scrambling_sign_quoted():
    """Re a(v) = -1 exactly (rung 25), the note's scrambling-side fact."""
    for v in (0.1, 0.5, 0.9):
        assert growth_coefficient(v).real == pytest.approx(-1.0, abs=1e-12)


def test_fold_verdicts_quoted():
    """The note quotes decider #2: folded connected fraction O(1) and
    growing toward small lambda; N-scaled folded ratio far from the
    continuation at lambda = 0.5."""
    with open(HERE / "chord_fold.json") as fh:
        d = json.load(fh)
    pts = {round(p["lam"], 3): p for p in d["points"]}
    assert abs(pts[0.5]["fold_connected_fraction"]["re"]) > 0.3
    assert abs(pts[0.5]["fold"]["re"] / pts[0.5]["fold"]["closed_re"]) > 10
