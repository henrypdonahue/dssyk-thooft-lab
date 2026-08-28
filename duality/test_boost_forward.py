"""Tests for the forward tower check (boost_forward.json).

Pinned here: the exact derivative moment route against brute-force
quadrature of the closed-form kernel, at machine precision; the Taylor
composition against mpmath's generic differentiator; the committed grid
cells, by re-executing the producer; and the structural facts, which are
the log-divergent mu_0 as eps -> 0, the CP = -1 coupling projection, and
controlled truncation.

And the verdict.  The tower's kurtosis band covers the SYK targets, but
it is regulator-dominated.  The difference between the full tower and the
ground state alone is at most about 5%.
"""

import json
import math
from pathlib import Path

import numpy as np
import pytest

import boost_forward as bf
from boost_mass import spectral_closed_form
from dictionary import load_reference

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def data():
    with open(HERE / "boost_forward.json") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def tower12():
    return bf.tower(12)


def test_derivative_moments_match_kernel_quadrature():
    """The K_0 Taylor-composition moments equal brute-force quadrature of
    boost_mass.spectral_closed_form with the same e^{-eps omega} regulator.
    The integrand is analytic and decaying, so the trapezoid rule is
    spectrally accurate; the two routes share nothing."""
    m, eps = 1.0, 0.5
    mus = bf.single_mass_moments(m, eps)
    ws = np.arange(-6.0, 80.0, 0.1)
    Sv = spectral_closed_form(ws, m, 1.0)
    for k in (0, 2, 4):
        quad = float(np.trapezoid(Sv * np.exp(-eps * ws) * ws ** k, ws))
        assert mus[k] == pytest.approx(quad, rel=1e-8)
    quad6 = float(np.trapezoid(Sv * np.exp(-eps * ws) * ws ** 6, ws))
    assert mus[6] == pytest.approx(quad6, rel=1e-6)


def test_series_route_matches_generic_differentiation():
    """The ODE-recurrence + composition Taylor coefficients agree with
    mpmath's generic mp.taylor of the same Wightman function."""
    import mpmath as mp
    m, eps = 0.8585, 0.4

    def W(s):
        return mp.besselk(0, m * 2.0 * mp.sin(s / 2.0)) / (2 * mp.pi)

    mine = bf.w_series(eps, m)
    with mp.workdps(40):
        ref = mp.taylor(W, mp.mpf(eps), 6)
        for a, b in zip(mine, ref):
            assert abs(float((a - b) / b)) < 1e-25


def test_committed_grid_reproduces(data, tower12):
    """Execute the producer path on two committed grid cells (one per
    coupling ansatz) and on the eps* bisection input; the engine is
    deterministic."""
    masses, w, ev, captured = tower12
    assert captured == pytest.approx(data["tower"]["parseval_captured"],
                                     rel=1e-12)
    for cell in ((1.0, 0.4, 0), (0.35, 0.05, 1)):
        ms, eps, s_pow = cell
        want = next(g for g in data["grid"]
                    if g["mscale"] == ms and g["eps"] == eps
                    and g["s_pow"] == s_pow)
        got = bf.tower_invariants(masses, w, ev, ms, eps, s_pow)
        for key in ("mu0", "mu2", "mu4", "mu6", "kurt", "r6"):
            assert got[key] == pytest.approx(want[key], rel=1e-10)
    # verdict range really is the min/max of the stored grid
    kurts = [g["kurt"] for g in data["grid"]]
    assert data["verdict"]["tower_kurt_range"] == pytest.approx(
        [min(kurts), max(kurts)], rel=1e-12)


def test_couplings_project_onto_cp_minus_tower(data, tower12):
    """The decay-constant weights select exactly the CP = -1 (symmetric-FLZ)
    levels of the certified reference spectrum, ground-state dominated."""
    masses, w, ev, captured = tower12
    sym = [r for r in load_reference() if r["CP"] == -1][:12]
    assert len(sym) == 12
    for m_i, ev_i, r in zip(masses, ev, sym):
        assert ev_i == pytest.approx(r["two_lambda"], rel=1e-9)
        assert m_i == pytest.approx(math.sqrt(r["two_lambda"]), rel=1e-9)
    assert w[0] == pytest.approx(0.888888, abs=2e-4)   # f_0 = 0.9428^2
    assert all(b < a for a, b in zip(w, w[1:]))        # falling weights
    assert captured > 0.99


def test_regulator_divergence_is_logarithmic(data):
    """mu_0(eps) grows by ~ln 2 per halving of eps (S -> 1/omega tail):
    the structural fact that makes every raw invariant regulator-dependent."""
    rows = data["divergence_eps_to_0"]
    assert [r["eps"] for r in rows] == [0.4, 0.2, 0.1, 0.05, 0.025]
    mu0 = [r["mu0"] for r in rows]
    incs = [b - a for a, b in zip(mu0, mu0[1:])]
    assert all(0.55 < d < 0.72 for d in incs), incs        # ~ln 2 = 0.693
    assert all(b > a for a, b in zip(incs, incs[1:]))      # -> ln 2
    # and kurtosis rises without bound along the same sequence
    kurts = [r["kurt"] for r in rows]
    assert all(b > a for a, b in zip(kurts, kurts[1:]))


def test_verdict_band_and_single_mass_dominance(data, tower12):
    """The two-sided verdict: the tower band COVERS the SYK targets (no
    structural exclusion), while the tower is indistinguishable from its
    ground state at the few-percent level (moments see one mass)."""
    lo, hi = data["syk_targets"]["band"]
    assert lo == pytest.approx(2.527, abs=0.01)
    assert hi == pytest.approx(4.808, abs=0.01)
    klo, khi = data["verdict"]["tower_kurt_range"]
    assert klo < lo and khi > hi
    assert data["verdict"]["syk_band_covered"]
    rel = data["verdict"]["max_tower_vs_single_mass_rel"]
    assert 0.005 < rel < 0.10, rel
    # every finite eps* is a genuine root of the reported target
    masses, w, ev, _ = tower12
    for row in data["eps_star"]:
        e = row["eps_star_mscale1"]
        if not math.isnan(e):
            k = bf.tower_invariants(masses, w, ev, 1.0, e)["kurt"]
            assert k == pytest.approx(row["target_kurt"], rel=1e-6)
    # truncation of the 12-level tower is controlled
    assert data["truncation"]["rel"] < 0.01
    assert data["truncation"]["residual_weight"] < 0.02


def test_a2_full_channel_target_is_live(data):
    """The frozen A_2 full-channel kurtosis targets trace to the exact
    finite-N machinery: recompute the cheapest point live."""
    import sys
    sys.path.insert(0, str(HERE.parent / "syk-self-averaging"))
    import exact_pair_moments as epm

    p, N = 4, 16
    pair = epm.channel_pair_moments(N, p, 2, ks=(0, 2, 4))
    diag = epm.channel_pair_moments(N, p, 2, marks=('i', 'i', 'i', 'i'),
                                    ks=(0, 2, 4))
    mu0 = N * diag["mu0"] + N * (N - 1) * pair["mu0"]
    mu2 = N * diag["mu2"] + N * (N - 1) * pair["mu2"]
    mu4 = N * diag["mu4"] + N * (N - 1) * pair["mu4"]
    kurt = mu4 * mu0 / mu2 ** 2
    frozen = next(r for r in data["syk_targets"]["a2_full"]
                  if r["p"] == p and r["N"] == N)
    assert kurt == pytest.approx(frozen["kurt"], abs=2e-4)
    # and the O_Delta targets trace to the committed chord JSON
    with open(HERE / "chord_moments.json") as fh:
        pts = json.load(fh)["points"]
    for row in data["syk_targets"]["odelta"][:4]:
        src = next(r for r in pts if r["lam"] == row["lam"])
        assert row["kurt"] == pytest.approx(
            src["channels"][row["channel"]]["kurt"], rel=1e-12)
