"""Tests for the NV-dressing identification note (NV-DRESSING.md).

These pin the two CHECKABLE halves of the identification between
Narovlansky-Verlinde's +-i beta_dS/4 dressing (2310.16994 Eqs. 33-34 and
61-62) and this repo's fake-period contour shifts, with beta_dS
identified as the fake temperature beta_fake = beta/v.

(a) THE 2-PT SECTOR IS THERMAL AT NV'S TEMPERATURE, exactly in lambda.
    beta_fake_c = beta_c/v and beta_eff = beta_c (1 + 1/v) are re-derived
    from the bridge.  Per-line detailed balance, which is NV Eq. 69 with
    beta_dS -> beta_fake, holds at float precision.  The independent
    contour-engine check holds at about 1e-14.  The hologram-point gap
    beta_eff/beta_fake - 1 is exactly v.

(b) THE DRESSING COMPOSITION ALGEBRA on the closed-form growth
    coefficient.  A single quarter-fake-period dressing is NV's "90
    degree rotation along the thermal circle", a -> -i a.  The opposed
    +- quarter dressings of the two probe pairs compose to the
    half-period antipodal flip, a -> -a.  That is identically equal to
    the one-sided Delta = pi/v shift, by translation invariance.

WHAT IS NOT TESTED HERE, by design: the constrained-sector step, which is
NV's Eq. 32 dressing plus the BRST projection.  NV-DRESSING.md states
that scope boundary.  These tests pin the contour realization only.
"""

import math

import numpy as np
import pytest

from antiscrambling import otoc_ratio
from chord import beta_chord_from_betaJ
from chord_fold import CONFIG_CROSSED
from largeq_anchor import v_of_betaJ
from signflip_probe import detailed_balance, growth_coefficient_shifted


def _coef_two_sided(v: float, d12: float, d34: float,
                    config=CONFIG_CROSSED) -> complex:
    """Growth coefficient with the (theta1, theta2) pair shifted by d12
    and the (theta3, theta4) pair by d34 -- the two-sided generalization
    of signflip_probe.growth_coefficient_shifted, needed because NV's
    O^+/O^- dressing displaces the two operator pairs in OPPOSITE
    directions (their Eq. 33)."""
    t1, t2, t3, t4 = config
    thLs = np.array([4.0, 5.0, 6.0, 7.0]) / v
    vals = np.array([complex(otoc_ratio(t1 + d12 + 1j * t, t2 + d12 + 1j * t,
                                        t3 + d34, t4 + d34, v))
                     for t in thLs])
    A = np.vstack([np.ones_like(thLs), np.exp(v * thLs),
                   np.exp(-v * thLs)]).T
    coef, *_ = np.linalg.lstsq(A.astype(complex), vals, rcond=None)
    return complex(coef[1])


@pytest.mark.parametrize("betaJ", [0.5, 2.0])
def test_quarter_dressings_compose_to_the_antipodal_flip(betaJ):
    """NV Eq. 34: O^-(tau) = O^+(tau + i beta_dS/2) -- two quarter-turns
    are one half-turn.  On the closed-form growth coefficient, with
    beta_dS/4 <-> pi/(2v):

      one-sided quarter shift:   a -> -i a   (the '90 degree rotation')
      opposed +- quarter shifts: a -> -a     (the antipodal flip)

    and the opposed-pair dressing equals the one-sided half shift
    IDENTICALLY (translation invariance of the corrected otoc_ratio),
    which is what lets signflip_probe's one-sided Delta = pi/v probe
    stand in for NV's two-sided dressing."""
    v = v_of_betaJ(betaJ)
    quarter = math.pi / (2.0 * v)
    a0 = _coef_two_sided(v, 0.0, 0.0)
    # baseline agrees with the producer's one-sided fit and closed form
    assert a0 == pytest.approx(growth_coefficient_shifted(v, 0.0), rel=1e-12)
    assert a0.real == pytest.approx(-1.0 / math.cos(math.pi * v / 2.0),
                                    rel=1e-12)
    a_quarter = _coef_two_sided(v, quarter, 0.0)
    a_nv = _coef_two_sided(v, +quarter, -quarter)     # O^+ pair vs O^- pair
    a_half = _coef_two_sided(v, 2.0 * quarter, 0.0)
    assert a_quarter == pytest.approx(-1j * a0, rel=1e-12)
    assert a_nv == pytest.approx(a_half, rel=1e-12)
    assert a_nv == pytest.approx(-a0, rel=1e-12)
    assert a_nv.real > 0 > a0.real                    # the sign flip itself


def test_dressed_frame_thermal_at_nv_temperature():
    """NV Eq. 69 detailed balance at beta_dS, with beta_dS <-> beta_fake:
    re-derive the scales from the bridge (not the stored JSON) and check
    the per-line law plus the independent contour-engine evaluation."""
    lam, betaJ = 2.0, 2.0
    v = v_of_betaJ(betaJ)
    # the defining relation behind v, asserted so the scale chain is closed
    assert math.pi * v / math.cos(math.pi * v / 2.0) == pytest.approx(
        betaJ, rel=1e-10)
    beta_c = beta_chord_from_betaJ(betaJ, lam)
    db = detailed_balance(lam, betaJ)
    assert db["dtau_c"] == pytest.approx(0.5 * beta_c / v, rel=1e-12)
    assert db["beta_fake_c"] == pytest.approx(beta_c / v, rel=1e-12)
    assert db["beta_eff"] == pytest.approx(beta_c * (1.0 + 1.0 / v),
                                           rel=1e-12)
    assert db["max_line_deviation"] < 1e-12
    assert db["engine_check_dev"] < 1e-12


def test_hologram_point_thermal_at_the_tomperature():
    """At the hologram point (beta -> 0) the dressed frame's detailed-
    balance temperature IS the tomperature: beta_eff/beta_fake - 1 = v
    exactly, and v -> 0 there.  This is the T = infinity limit of NV's
    KMS matching (their Sec. 2.3 / 3.3 choice beta_dS = 2 pi)."""
    db = detailed_balance(1.0, 0.02)
    v = v_of_betaJ(0.02)
    assert v < 0.01
    assert db["beta_eff"] / db["beta_fake_c"] - 1.0 == pytest.approx(
        v, rel=1e-9)
    assert db["max_line_deviation"] < 1e-10
    assert db["engine_check_dev"] < 1e-12
