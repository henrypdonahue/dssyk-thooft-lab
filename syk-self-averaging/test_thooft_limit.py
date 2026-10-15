"""Tests for thooft_limit.py: the fixed-p, N -> infinity shape of the
singlet channels against the lambda -> 0 ladders."""
import json
from fractions import Fraction
from pathlib import Path

import pytest

import thooft_limit as t

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def art():
    with open(HERE / "thooft_limit.json") as fh:
        return json.load(fh)


def test_ladder_moments_are_the_sech_cumulants():
    """kurtosis 3 + 1/Delta and mu6/mu2^3 = 15 + 15/Delta + 4/Delta^2,
    independent of the scale; <w^4>/<w^2> = 4 p^2 (1 + 3 Delta)."""
    for d in (2, 3, 5):
        m = t.ladder_moments(d, 7)
        assert m[4] / m[2] ** 2 == 3 + Fraction(1, d)
        assert m[6] / m[2] ** 3 == 15 + Fraction(15, d) + Fraction(4, d * d)
        assert m[4] / m[2] == 4 * 49 * (1 + 3 * d)
    assert t.ladder_moments(2, 1)[2] == 8          # 2 Delta scriptJ^2, scriptJ^2 = 2 p^2


def test_closed_forms_tend_to_the_ladders():
    big = 10 ** 6
    assert float(t.kurtosis_closed(big)) == pytest.approx(3.5, abs=1e-6)
    assert float(t.hyper_closed(big)) == pytest.approx(23.5, abs=1e-5)
    assert t.a3_ratio_closed(big) / (40 * big * big) == pytest.approx(1.0, abs=1e-5)
    assert t.omega4_closed(big) / (224 * big ** 4) == pytest.approx(1.0, abs=1e-5)


def test_tables_reproduce_from_the_stored_coefficients(art):
    """The tables are a deterministic function of the stored exact
    coefficients."""
    again = t.tables(art["coefficients"])
    assert again["a2"] == art["a2"]
    assert again["a3"] == art["a3"]


def test_a2_laws_hold_at_every_p(art):
    rows = art["a2"]
    assert [r["p"] for r in rows] == [4, 6, 8, 10, 12]
    for r in rows:
        assert r["omega2_closed_ok"] and r["omega4_closed_ok"] and r["kurtosis_closed_ok"]
        if "hyper" in r:
            assert r["hyper_closed_ok"], r["p"]
    assert sum("hyper" in r for r in rows) >= 3


def test_a2_shape_converges_to_the_delta2_ladder(art):
    rows = art["a2"]
    for r in rows:
        assert r["omega2_over_ladder"] == pytest.approx(1 - 2 / r["p"], rel=1e-14)
    k = [r["kurtosis_float"] for r in rows]
    assert k == sorted(k) and all(x < 3.5 for x in k)
    h = [r["hyper_float"] for r in rows if "hyper" in r]
    assert h == sorted(h) and all(x < 23.5 for x in h)


def test_a3_selects_the_delta3_ladder(art):
    """<w^4>/<w^2> = 4p(10p - 21): its p^2 coefficient is the Delta = 3
    ladder's 40, not the Delta = 2 ladder's 28."""
    rows = art["a3"]
    assert len(rows) >= 3
    for r in rows:
        assert r["ratio_closed_ok"], r["p"]
    fr = [r["ratio_over_ladder3"] for r in rows]
    assert fr == sorted(fr) and fr[-1] > 0.7


def test_a2_moments_recomputed_at_p4():
    r = t.a2_shape(4)
    assert r["omega2"] == t.omega2_closed(4) == 64
    assert r["omega4"] == t.omega4_closed(4) == 13312
    assert r["kurtosis"] == Fraction(13, 4)


@pytest.mark.slow
def test_a2_sixth_moment_recomputed_at_p4():
    r = t.a2_shape(4, with_mu6=True)
    assert r["hyper"] == t.hyper_closed(4) == Fraction(307, 16)


@pytest.mark.slow
def test_a3_ratio_recomputed_at_p4():
    assert t.a3_ratio(4)["ratio"] == t.a3_ratio_closed(4) == 304
