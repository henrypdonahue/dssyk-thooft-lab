"""Tests for the inverse conditioning study (boost_inverse.json).

Pinned here: the arbitrary-order engine against brute-force quadrature of
the closed-form kernel with the same e^{-eps omega} regulator, and against
boost_forward's fixed-order route; the dps-doubling gate's certificate;
the conditioning cells, by re-executing the matrix + SVD path live; the
separability scaling, by re-executing one fit; the node observables, by
re-finding the zeros; and the verdict numbers, recomputed from the stored
arrays.

The headline is negative and quantitative: super-geometric singular-value
decay with measured, eps-dependent exponents, saturation at K = L
moments, and a ~1 regulator leverage that node positions are exactly free
of.  Everything is conditional on the e^{-eps omega} regulator
convention, and the artifact records the eps-dependence."""

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

import boost_forward as bf
import boost_inverse as bi
from boost_mass import spectral_closed_form

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def data():
    with open(HERE / "boost_inverse.json") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# study 1: the engine
# ---------------------------------------------------------------------------
def test_engine_matches_kernel_quadrature_to_high_order():
    """The certified Taylor-composition moments equal brute-force trapezoid
    quadrature of boost_mass.spectral_closed_form with the same regulator,
    through k = 20 (the quadrature is spectrally accurate; the two routes
    share nothing).  Odd orders included: the conditioning map uses them."""
    m, eps = 1.0, 1.0
    mus, err, _ = bi.regulated_moments(m, eps, kmax=20, rel_tol=1e-40)
    assert err < 1e-40
    ws = np.arange(-8.0, 110.0, 0.1)
    Sv = spectral_closed_form(ws, m, 1.0)
    reg = Sv * np.exp(-eps * ws)
    for k in (0, 1, 2, 3, 5, 8, 12, 16, 20):
        quad = float(np.trapezoid(reg * ws ** k, ws))
        assert float(mus[k]) == pytest.approx(quad, rel=1e-10), k


def test_engine_matches_boost_forward_fixed_order():
    """At kmax <= 6 the generalized engine reproduces boost_forward's
    single_mass_moments exactly (same mathematics, independent code path
    for the composition length)."""
    for m, eps in ((0.8585230027743485, 0.4), (1.0, 1.6), (2.5, 0.05)):
        ref = bf.single_mass_moments(m, eps)
        mus, _, _ = bi.regulated_moments(m, eps, kmax=6)
        for k in (0, 2, 4, 6):
            assert float(mus[k]) == pytest.approx(ref[k], rel=1e-12)


def test_engine_gate_is_certified_and_deterministic():
    """The dps-doubling gate: the certificate beats rel_tol, a repeat call
    is bit-identical, and a run forced to higher starting precision agrees
    to well below the certificate."""
    m, eps = 3.4278, 0.1
    ws1, err1, dps1 = bi.w_series_k(eps, m, kmax=23, rel_tol=1e-60)
    ws2, err2, dps2 = bi.w_series_k(eps, m, kmax=23, rel_tol=1e-60)
    assert err1 < 1e-60 and (err1, dps1) == (err2, dps2)
    assert all(a == b for a, b in zip(ws1, ws2))
    # a run gated at tighter tolerance (hence higher precision) agrees
    ws3, _, dps3 = bi.w_series_k(eps, m, kmax=23, rel_tol=1e-80)
    assert dps3 > dps1
    with mp.workdps(90):
        for a, b in zip(ws1, ws3):
            assert abs((a - b) / b) < mp.mpf("1e-58")


def test_engine_rejects_out_of_range_regulator():
    with pytest.raises(ValueError):
        bi.w_series_k(0.0, 1.0)
    with pytest.raises(ValueError):
        bi.w_series_k(7.0, 1.0)


# ---------------------------------------------------------------------------
# study 2: conditioning
# ---------------------------------------------------------------------------
def test_conditioning_cell_reexecutes(data):
    """Rebuild the weights->moments map at the reference regulator and
    re-run the SVD path live; the stored spectra reproduce."""
    eps = data["conditioning"]["fine_scan_eps_ref"]["eps"]
    masses = bi.tower_masses(12)
    assert [float(m) for m in masses] == pytest.approx(
        data["tower"]["masses"], rel=1e-12)
    rows, err = bi.moment_matrix(masses, eps)
    assert err < data["conditioning"]["entry_rel_tol"]
    cell = next(c for c in data["conditioning"]["cells"]
                if c["eps"] == eps and c["K"] == 24)
    lg_raw = bi.singular_spectrum(rows, 24, column_normalize=False)
    lg_shape = bi.singular_spectrum(rows, 24, column_normalize=True)
    assert lg_raw == pytest.approx(cell["lg_sv_raw"], abs=1e-8)
    assert lg_shape == pytest.approx(cell["lg_sv_shape"], abs=1e-8)
    # law fit re-derives from the stored spectrum
    fit = next(r for r in data["conditioning"]["law_fit"]
               if r["eps"] == eps)
    a, b, rms = bi.decay_fit(cell["lg_sv_raw"])
    assert (a, b, rms) == pytest.approx((fit["a"], fit["b"], fit["rms"]),
                                        abs=1e-10)


def test_conditioning_headline_numbers(data):
    """The measured conditioning law: super-geometric decay (b ~ 0.05),
    worse at smaller eps (a decreasing in eps), double-precision level
    counts, and the marginal last level at the reference regulator."""
    fits = data["conditioning"]["law_fit"]
    assert [r["eps"] for r in fits] == list(bi.EPS_GRID)
    a_vals = [r["a"] for r in fits]
    assert all(x > y for x, y in zip(a_vals, a_vals[1:]))  # small eps worse
    for r in fits:
        if r["eps"] >= 0.2:
            assert 0.04 < r["b"] < 0.07, r  # accelerating decay, measured
            assert r["rms"] < 0.3
    # at a*eps = 0.1 the last singular values floor near 1e-17, the
    # quadratic law fails (recorded, not hidden)
    fit01 = next(r for r in fits if r["eps"] == 0.1)
    assert fit01["rms"] > 0.5
    c01 = next(c for c in data["conditioning"]["cells"]
               if c["eps"] == 0.1 and c["K"] == 24)
    assert c01["lg_sv_raw"][-3] - c01["lg_sv_raw"][-1] < 1.0  # the floor
    counts = {c["eps"]: c["n_above"]["1e-16"]
              for c in data["conditioning"]["cells"] if c["K"] == 24}
    assert counts == {0.1: 9, 0.2: 10, 0.4: 12, 0.8: 12, 1.6: 12}
    ref = next(c for c in data["conditioning"]["cells"]
               if c["eps"] == bi.EPS_REF and c["K"] == 24)
    # last level at the reference regulator: present but marginal
    assert 15.0 < ref["lg_sv_raw"][0] - ref["lg_sv_raw"][-1] < 16.0
    # raw and shape spectra agree closely: Boltzmann column suppression
    # is negligible against shape degeneracy (< 0.25 decades at K = 24,
    # < 0.5 in the row-starved small-K cells)
    for c in data["conditioning"]["cells"]:
        bound = 0.25 if c["K"] == 24 else 0.5
        for x, y in zip(c["lg_sv_raw"], c["lg_sv_shape"]):
            assert abs(x - y) < bound, (c["eps"], c["K"])


def test_moments_beyond_L_saturate(data):
    """The answer to 'how many moments': the double-precision level count
    saturates at K = L = 12 at every eps (one-count flicker at 0.2), the
    fine K-scan at the reference eps is monotone to saturation, and the
    residual singular-value motion from K = 12 to 24 is eps-graded --
    < 0.5 decades at a*eps <= 0.4, ~4 at a*eps = 1.6, where high moments
    still probe the mass window."""
    fine = data["conditioning"]["fine_scan_eps_ref"]
    counts = [r["n_double"] for r in fine["rows"]]
    assert all(b >= a for a, b in zip(counts, counts[1:]))
    assert fine["K_saturation_double"] == 12
    assert counts[-1] == fine["n_double_final"] == 12
    moves = {}
    for eps in bi.EPS_GRID:
        c12 = next(c for c in data["conditioning"]["cells"]
                   if c["eps"] == eps and c["K"] == 12)
        c24 = next(c for c in data["conditioning"]["cells"]
                   if c["eps"] == eps and c["K"] == 24)
        assert abs(c12["n_above"]["1e-16"] - c24["n_above"]["1e-16"]) <= 1
        moves[f"{eps}"] = max(
            abs((x - c12["lg_sv_raw"][0]) - (y - c24["lg_sv_raw"][0]))
            for x, y in zip(c12["lg_sv_raw"], c24["lg_sv_raw"]))
    assert moves == pytest.approx(
        data["conditioning"]["sv_move_K12_to_K24"], abs=1e-9)
    assert moves["0.1"] < 0.5 and moves["0.4"] < 0.5
    assert moves["1.6"] > 2.0


def test_cp_subtower_and_even_rows(data):
    """The CP = -1 6-level sub-tower is the easy case (all levels above
    double precision), and restricting to even moment orders (the
    symmetrized correlator's observables) does not rescue conditioning."""
    cp = data["conditioning"]["cp_minus_subtower_eps_ref"]
    assert len(cp["lg_sv_raw"]) == 6
    assert cp["n_above"]["1e-16"] == 6
    even = data["conditioning"]["even_only_eps_ref"]["lg_sv_raw"]
    full = next(c for c in data["conditioning"]["cells"]
                if c["eps"] == bi.EPS_REF and c["K"] == 24)["lg_sv_raw"]
    assert even[-1] - even[0] <= (full[-1] - full[0]) + 0.6


# ---------------------------------------------------------------------------
# study 3: separability
# ---------------------------------------------------------------------------
def test_separability_reexecutes_and_scales(data):
    """One live fit reproduces the artifact; the stored scan obeys the
    delta^2 law; two moments fit exactly (so 3 is the minimum); and the
    threshold delta_min = sqrt(eta/C)."""
    sep = data["separability"]
    row = next(r for r in sep["scan"] if r["delta"] == 0.1)
    live = bi.two_mass_fit(1.0, 0.1, bi.EPS_REF)
    assert live["mstar"] == pytest.approx(row["mstar"], rel=1e-10)
    assert live["res"][2:] == pytest.approx(row["res"][2:], rel=1e-6)
    # matched moments are exact: residuals 0 and 1 at solver precision
    # (the mass iterate stays an mpf; findroot tol 1e-28 on the ratio)
    assert live["res"][0] < 1e-20 and live["res"][1] < 1e-20
    assert sep["min_moments"] == 3
    assert sep["fit_exponent"] == pytest.approx(2.0, abs=0.01)
    # quarter the splitting -> 1/16 the signal, across the whole scan
    r2 = [r["res"][2] for r in sep["scan"]]
    for a, b in zip(r2, r2[1:]):
        assert a / b == pytest.approx(4.0, rel=0.03)
    assert sep["delta_min"]["1e-16"] == pytest.approx(
        math.sqrt(1e-16 / sep["C_eps_ref"]), rel=1e-9)
    # eps-dependence recorded, not hidden: C spans an order of magnitude
    Cs = {v["eps"]: v["C"] for v in sep["eps_variants"]}
    assert Cs[0.1] < sep["C_eps_ref"] < Cs[1.6]
    assert Cs[1.6] / Cs[0.1] > 5
    # the unresolved-pair mass bias is second order with an O(0.2) coeff
    assert 0.1 < sep["mstar_bias_coeff"] < 0.4


def test_higher_moments_add_little_separation(data):
    """The residual profile over k = 2..7 is flat within a factor ~1.5:
    more moments do not substitute for precision (the conditioning
    saturation, seen from the nonlinear side)."""
    for r in data["separability"]["scan"]:
        tail = r["res"][2:]
        assert max(tail) / min(tail) < 1.6, r["delta"]


# ---------------------------------------------------------------------------
# study 4: nodes against moments
# ---------------------------------------------------------------------------
def test_nodes_refind_and_match_boost_md(data):
    """Live sign-change scan + bracketed refinement reproduces the stored
    nodes, and the first nodes match BOOST.md's calibration (omega ~ 3.0
    at m = 1, ~ 4.4 at m = 2).  The regulator moves no node: e^{-eps w}
    is strictly positive, so S and S e^{-eps w} share every zero."""
    n1 = bi.scan_nodes(1.0, w_max=10.0)
    assert n1 == pytest.approx(data["nodes"]["nodes_m1"], abs=1e-8)
    assert n1[0] == pytest.approx(2.9625, abs=1e-3)
    assert data["nodes"]["nodes_m2"][0] == pytest.approx(4.4255, abs=1e-3)
    assert data["nodes"]["eps_leverage_of_nodes"] == 0.0
    # node mass-sensitivity: order one, falling with node index
    lam = data["nodes"]["lambda_node"]
    assert lam[0] == pytest.approx(bi.node_mass_sensitivity(1.0, n1[0]),
                                   abs=1e-4)
    assert all(a > b for a, b in zip(lam, lam[1:]))
    assert 0.3 < lam[2] < lam[0] < 0.7


def test_regulator_leverage_of_moments(data):
    """The decisive number: |d ln kurt/d ln eps| / |d ln kurt/d ln m| ~ 1
    across the eps grid -- one digit of regulator uncertainty costs one
    digit of mass.  One row re-executed live."""
    rows = data["nodes"]["kurt_sens"]
    for r in rows:
        assert 0.85 < r["leverage"] < 1.1, r
    ref = next(r for r in rows if r["eps"] == bi.EPS_REF)
    live = bi.kurt_sensitivities(1.0, bi.EPS_REF)
    for key in ("kurt", "dlnK_dlnm", "dlnK_dlneps", "leverage"):
        assert live[key] == pytest.approx(ref[key], rel=1e-4)
    # and kurt itself agrees with the boost_forward moment route
    mus = bf.single_mass_moments(1.0, bi.EPS_REF)
    assert ref["kurt"] == pytest.approx(
        mus[4] * mus[0] / mus[2] ** 2, rel=1e-9)


def test_digits_table_scaling_laws(data):
    """Node information transfers as sqrt(eta) (double zero of S = K^2),
    moment information linearly in eta; the 1%-eps floor is
    eta-independent.  The crossover justifies the roadmap verdict."""
    tab = data["nodes"]["digits_table"]
    assert [r["eta"] for r in tab] == [1e-4, 1e-8, 1e-12]
    for a, b in zip(tab, tab[1:]):
        assert a["dm_node"] / b["dm_node"] == pytest.approx(100.0, rel=1e-6)
        assert a["dm_kurt_eps_exact"] / b["dm_kurt_eps_exact"] == \
            pytest.approx(1e4, rel=1e-6)
        assert a["dm_floor_eps_1pct"] == b["dm_floor_eps_1pct"]
    # with eps known only to 1%, nodes beat moments at every eta shown
    for r in tab:
        assert r["dm_node"] < r["dm_floor_eps_1pct"]
    pref = data["nodes"]["node1_curvature"]["transfer_prefactor"]
    assert tab[0]["dm_node"] == pytest.approx(pref * 1e-2, rel=1e-9)
    assert 0.2 < pref < 0.8


def test_tower_node_filling(data):
    """The coupling-weighted tower fills the ground node: S_tower at the
    ground kernel's first zero is a pure excited-level signal.  Re-execute
    the cheap evaluation live and check the background-free contrast."""
    tw = data["nodes"]["tower"]
    masses, wgt, _, _ = bf.tower(12)
    w1 = tw["w1_ground"]
    with mp.workdps(30):
        g0 = bi.scaled_amplitude(w1, float(masses[0]))
        s_tow = sum(float(wgt[i])
                    * bi.scaled_amplitude(w1, float(masses[i])) ** 2
                    for i in range(12)) / math.pi
    assert abs(g0) < 1e-8                       # really the ground node
    assert s_tow == pytest.approx(tw["S_tower_at_ground_node"], rel=1e-8)
    assert tw["S_tower_min"] <= tw["S_tower_at_ground_node"]
    assert tw["min_depth"] <= tw["fill_at_ground_node"] < 0.2
    assert tw["min_depth"] > 0.01               # filled, not a zero
    assert abs(tw["node_shift_rel"]) < 0.01


# ---------------------------------------------------------------------------
# the verdict
# ---------------------------------------------------------------------------
def test_verdict_consistency(data):
    """Verdict numbers re-derive from the stored sections; the headline is
    the negative result, stated conditionally on the regulator."""
    v = data["verdict"]
    assert v["headline_negative"] is True
    ref = next(c for c in data["conditioning"]["cells"]
               if c["eps"] == bi.EPS_REF and c["K"] == 24)
    assert v["amplification_last_level_ref"] == pytest.approx(
        ref["lg_sv_raw"][0] - ref["lg_sv_raw"][-1], abs=1e-9)
    counts = {f"{c['eps']}": c["n_above"]["1e-16"]
              for c in data["conditioning"]["cells"] if c["K"] == 24}
    assert v["n_levels_double_precision"] == counts
    fit = next(r for r in data["conditioning"]["law_fit"]
               if r["eps"] == bi.EPS_REF)
    assert v["digits_per_level_ref"]["first"] == pytest.approx(
        fit["a"] + fit["b"], rel=1e-9)
    assert v["digits_per_level_ref"]["last"] > \
        2 * v["digits_per_level_ref"]["first"]
    assert "e^{-eps omega}" in v["statement"]
    assert "node-level" in v["statement"]
