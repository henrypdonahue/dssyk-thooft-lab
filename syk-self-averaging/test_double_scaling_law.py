"""Tests for double_scaling_law.py: the closed forms of the pair channel's
double-scaling constants, the exact sequences on three lines, the
diagonal-pair O(Delta) cancellation, and the B(2,2; p) closed form."""

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import pytest

from double_scaling_law import (B22_closed_form, P_K, _R, c_closed_form,
                                cancellation_table, delta_slope_symbolic,
                                lambda_expansion, limit_estimate)
from pair_moment_laws import mu_exact

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "duality"))


@pytest.fixture(scope="module")
def art():
    with open(HERE / "double_scaling_law.json") as fh:
        return json.load(fh)


@pytest.mark.parametrize("n,k", [(1, 0), (1, 2), (2, 2), (3, 2), (1, 6)])
def test_delta_slopes_symbolic(n, k):
    """The symbolic walker reproduces the stored P_K for K <= 8, and the
    slope depends on K = k + 2n only (the exact ad-shift)."""
    import sympy as sp
    q = sp.Symbol("q", positive=True)
    P = delta_slope_symbolic(n, k)
    K = k + 2 * n
    for qv in (sp.Rational(1, 3), sp.Rational(7, 10)):
        assert abs(float(P.subs(q, qv)) - P_K(K, float(qv))) < 1e-12


def test_P10_numerically():
    """K = 10 (too slow symbolically for the fast suite): the walker's
    finite-difference Delta-slope matches the stored polynomial."""
    import diagonal_chord as dc
    for lam in (1.0, 2.0):
        q = math.exp(-lam)
        d0 = 1.0 / 2048
        r1 = dc.channel_moments_diag(q, 2 * d0, 2, ks=(0, 2, 6))
        r2 = dc.channel_moments_diag(q, d0, 2, ks=(0, 2, 6))
        slope = 2 * (r2["mu6"] / d0) - r1["mu6"] / (2 * d0)
        assert abs(slope - lam * P_K(10, q)) < 2e-4 * lam * P_K(10, q)


def test_closed_forms_vanish_at_q_1():
    """P_K ~ (1-q)^(K/2-1): every slope dies as lambda -> 0 like lambda^(K/2)."""
    for K in (4, 6, 8, 10):
        assert P_K(K, 1.0) == 0
        assert P_K(K, 0.0) == 4 * _R[K][0]


def test_sequences_recompute(art):
    """A cheap row of every line is recomputed exactly from mu_exact."""
    for key, rows in art["sequences"].items():
        r = rows[0]
        if r.get("source"):
            continue                       # imported ladder rows
        assert str(mu_exact(r["N"], r["p"], r["n"], 2)) == r["mu2_exact"]


def test_limits_match_closed_forms(art):
    """On every line the 1/p-extrapolated constant agrees with
    lambda P_K(e^-lambda): to 1e-5 relative at K = 6 (lambda = 1, 2),
    1e-4 at K = 8 (lambda = 1, 2) and 1e-2 on the short lambda = 4 lines."""
    tol = {"lam1_n2": 2e-5, "lam2_n2": 2e-5, "lam4_n2": 5e-4,
           "lam1_n3": 2e-4, "lam2_n3": 5e-4, "lam4_n3": 1e-2}
    for key, est in art["limits"].items():
        assert est["estimated"], key
        assert est["rel_dev"] < tol[key], (key, est["c_inf_central"],
                                          est["closed_form"], est["rel_dev"])
        assert abs(est["closed_form"] - c_closed_form(est["lam"], est["K"])) < 1e-12


def test_ladder_c2_limit_is_the_closed_form(art):
    """The committed c2(p) (p <= 48) extrapolates to 4 lambda (1-q)^2 (5+3q)
    at lambda = 2 = 32.334189, and the windowed 32.334(1) contains it."""
    rows = art["sequences"]["lam2_n2"]
    assert rows[-1]["p"] >= 48
    est = limit_estimate(rows)
    assert abs(est["c_inf_central"] - 32.334189) < 3e-5
    with open(HERE / "pair_moment_ladder.json") as fh:
        lad = json.load(fh)["double_scaling"]["c2_limit"]
    assert abs(lad["c_inf_central"] - c_closed_form(2.0, 6)) < lad["spread"]


def test_cancellation_at_finite_N(art):
    """mu_2^diag p -> P_6 from above, (N-1) mu_2^pair p -> -P_6 from below,
    and the sum x p^2 stays bounded (O(Delta^2)); same for mu_0 with P_4."""
    rows = art["cancellation_lambda2"]
    P6 = c_closed_form(2.0, 6)
    P4 = c_closed_form(2.0, 4)
    d2 = [r["diag_mu2_times_p"] for r in rows]
    p2 = [r["pair_mu2_Nm1_times_p"] for r in rows]
    assert all(d2[i + 1] < d2[i] for i in range(len(d2) - 1)) and d2[-1] > P6
    assert all(p2[i + 1] < p2[i] for i in range(len(p2) - 1)) and p2[-1] > -P6
    assert abs(d2[-1] - P6) < 0.02 * P6 and abs(p2[-1] + P6) < 0.05 * P6
    s2 = [r["sum_mu2_times_p2"] for r in rows]
    assert max(s2) < 60 and min(s2) > 20
    d0 = [r["diag_mu0_times_p"] for r in rows]
    assert d0[-1] > P4 and abs(d0[-1] - P4) < 0.06 * P4
    # live recomputation of one row
    live = cancellation_table(ps=(6,))[0]
    for key in ("diag_mu2_times_p", "sum_mu2_times_p2"):
        assert abs(live[key] - rows[1][key]) < 1e-12 * abs(rows[1][key])


def test_B22_closed_form(art):
    """The N^-5 coefficient of mu_2^pair(n=2): the committed exact B at
    p = 4..12 fits 64 p^3 (p-1)(11p^3 - 36p^2 + 29p - 3), which reproduces
    the committed p = 2 value out of sample, and whose leading 704 p^7 is
    the lambda^5 coefficient of the closed form."""
    assert art["B22_check"]["all_match"]
    assert B22_closed_form(2) == -512
    exp = {e["m"]: e for e in lambda_expansion()}
    assert exp[0]["C_m_over_p_power"] == "-256"
    assert exp[1]["C_m_over_p_power"] == "704"


def test_A22_large_p_matches_lambda4_coefficient():
    """A(2,2; p) = -64 p^3 (p-1)(4p-5) -> -256 p^5 = -gamma_4 2^4 p^5."""
    p = 10 ** 6
    A = -64 * p ** 3 * (p - 1) * (4 * p - 5)
    # Keep the known finite-p correction: at p=10^6 it is 5.76e-4,
    # larger than the previous 1e-4 threshold despite the correct limit.
    correction = Fraction(A, p ** 5) + 256
    assert correction == Fraction(576, p) - Fraction(320, p ** 2)
    assert abs(correction) < Fraction(576, p)
