"""Tests for the exact 1/N laws (pair_moment_laws.py).

Pinned here: the degree deficits, which are the exact powers of the
leading 1/N laws; the exact leading coefficients, including the closed
forms in p certified out of sample at p = 14; the
mu_4(n=1) against mu_2(n=2) leading-coefficient coincidence; and
cross-artifact consistency with the committed float campaign in
exact_pair_moments.json.

The extraction machinery itself is re-executed live on the cheap laws.
That is producer re-execution, not just JSON reading.
"""

import json
from fractions import Fraction
from pathlib import Path

import pytest

from pair_moment_laws import _lagrange_in_p, _poly_eval, law

HERE = Path(__file__).resolve().parent

# pair_moment_laws.json is committed: missing = broken checkout, fail
with open(HERE / "pair_moment_laws.json") as fh:
    DATA = json.load(fh)
LAWS = {(r["n"], r["p"], r["k"]): r for r in DATA["laws"]}


def test_producer_reexecution_cheap_laws():
    """Execute the extractor, do not just read its output: the three
    sub-second laws, recomputed end to end (sampling, interpolation,
    series, out-of-sample certification -- law() raises internally if
    any of those fail)."""
    r = law(2, 1, 2)
    assert (r["deficit"], r["A"]) == (3, "-64")
    r = law(2, 2, 2)
    assert (r["deficit"], r["A"], r["B"]) == (4, "-1536", "-512")
    r = law(4, 2, 2)
    assert (r["deficit"], r["A"]) == (4, "-135168")
    # and each matches its committed row exactly
    assert LAWS[(2, 4, 2)]["A"] == "-135168"
    assert LAWS[(2, 4, 2)]["P_coeffs"] == r["P_coeffs"]


def test_k2_deficits_are_n_plus_2():
    """The fixed-p law mu_2 ~ N^-(n+2), exact at every computed (n, p):
    the fitted tail slopes 2.96 / 3.90 / 4.82 of the float campaign are
    the shadows of deficits 3 / 4 / 5, and the pattern holds through
    n = 4."""
    ks2 = [(n, p) for (n, p, k) in LAWS if k == 2]
    assert len(ks2) >= 17          # the campaign grid
    for (n, p) in ks2:
        assert LAWS[(n, p, 2)]["deficit"] == n + 2, (n, p)
    assert {n for (n, _) in ks2} == {1, 2, 3, 4}


def test_k4_deficits_are_n_plus_3():
    """mu_4 dies one power faster than mu_2 in the pair channel."""
    ks4 = [(n, p) for (n, p, k) in LAWS if k == 4]
    assert len(ks4) >= 5
    for (n, p) in ks4:
        assert LAWS[(n, p, 4)]["deficit"] == n + 3, (n, p)


def test_A_closed_form_n1():
    """A(n=1, p) = -16 p^2 (p-1), certified out-of-sample at p = 14 by
    the campaign; re-asserted here against every committed sample."""
    cf = next(c for c in DATA["A_of_p"] if c["n"] == 1)
    assert cf["identified"] and cf["p_check"] == 14
    coeffs = [Fraction(c) for c in cf["A_poly_coeffs"]]
    for (n, p, k), r in LAWS.items():
        if n == 1 and k == 2:
            assert Fraction(r["A"]) == -16 * p * p * (p - 1), p
            assert _poly_eval(coeffs, Fraction(p)) == Fraction(r["A"])


def test_A_closed_form_n2():
    """A(n=2, p) = -64 p^3 (p-1)(4p-5), certified out-of-sample at
    p = 14."""
    cf = next(c for c in DATA["A_of_p"] if c["n"] == 2)
    assert cf["identified"] and cf["p_check"] == 14
    for (n, p, k), r in LAWS.items():
        if n == 2 and k == 2:
            assert Fraction(r["A"]) == -64 * p ** 3 * (p - 1) * (4 * p - 5), p
    # the polynomial the campaign certified is that closed form
    coeffs = [Fraction(c) for c in cf["A_poly_coeffs"]]
    assert coeffs == [Fraction(x) for x in (0, 0, 0, -320, 576, -256)]


def test_mu4_n1_equals_mu2_n2_at_leading_order():
    """The measured coincidence: the exact leading coefficient of
    mu_4(n=1) equals that of mu_2(n=2) at every shared p -- consistent
    with ad_H-shifting between channels (ad(X_1) = X_2 + cross term,
    with the cross term entering beyond leading order).  Stated as a
    measured identity of the exact coefficients, not a derived one."""
    ps = sorted({p for (n, p, k) in LAWS if n == 1 and k == 4})
    assert ps                       # the campaign computed some
    for p in ps:
        assert LAWS[(1, p, 4)]["A"] == LAWS[(2, p, 2)]["A"], p
        assert LAWS[(1, p, 4)]["deficit"] == LAWS[(2, p, 2)]["deficit"]


def test_out_of_sample_certification_recorded():
    """Every committed law was certified against a direct exact
    evaluation at an N outside its interpolation grid (law() raises on
    mismatch, so presence in the artifact is the certificate; the
    recorded N pins which point)."""
    for key, r in LAWS.items():
        assert r["out_of_sample_N"] > 0, key
        assert len(r["P_coeffs"]) - 1 == r["deg_P"], key


def test_double_scaling_consistent_with_float_campaign():
    """Cross-artifact gate: the exact c(p) = |mu_2| p N on
    lambda_chord = 2 must reproduce the committed float campaign's
    double-scaling rows (exact_pair_moments.json) at the overlapping
    (p, N) points, to float-assembly precision."""
    with open(HERE / "exact_pair_moments.json") as fh:
        float_ds = json.load(fh)["double_scaling"]
    float_mu2 = {(r["p"], r["N"]): r["n2"]["mu2"] for r in float_ds}
    rows = {(r["p"], r["N"]): r for r in DATA["double_scaling"]["n2"]}
    overlap = set(float_mu2) & set(rows)
    assert len(overlap) == 4        # (4,16) (6,36) (8,64) (10,100)
    for key in overlap:
        assert rows[key]["mu2"] == pytest.approx(float_mu2[key], rel=1e-9)
    # and the exact scan extends beyond the float one
    assert max(p for (p, _) in rows) >= 14


def test_lagrange_helper():
    pts = [(Fraction(x), Fraction(x) ** 3 - 2) for x in (1, 2, 3, 5, 7)]
    poly = _lagrange_in_p([(int(x), y) for x, y in pts])
    assert poly == [Fraction(-2), Fraction(0), Fraction(0), Fraction(1)]
