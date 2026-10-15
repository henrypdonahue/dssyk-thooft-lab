"""Tests for horizon_tower.py: the orthogonalized tower, the quasinormal
ladder of each chord, and the Rindler calibration."""
import json
import math
from pathlib import Path

import numpy as np
import pytest

import horizon_tower as h
import singlet_tower as st

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def art():
    with open(HERE / "horizon_tower.json") as fh:
        return json.load(fh)


# -- the orthogonalized tower ------------------------------------------------
def test_new_chord_enters_with_coefficient_four():
    """O~_n = 4 O^_n: the top chord of A_n has constant dressing 4, and its
    weight 16 norm_n is the sector moment; for n <= 4 it is the closed form."""
    lam = 2.0
    q = math.exp(-lam)
    closed = {2: st.norm2(lam, q), 3: st.norm3(lam, q), 4: st.norm4(lam, q)}
    for r in h.orthogonal_tower(lam, nmax=4):
        assert r["top_kernel"] == {"0,0": 4.0}
        assert r["weight_sector"] == pytest.approx(r["weight_16norm"], rel=1e-12)
        assert r["norm"] == pytest.approx(closed[r["n"]], rel=1e-10)


def test_new_chord_shares_match_singlet_tower(art):
    """The weight of O~_n in A_n at lambda = 2 is the faint new-chord share
    quoted in singlet_tower.py: 14.9%, 1.6%, 0.18%."""
    shares = {r["n"]: r["share_of_dynamic"] for r in art["orthogonal_tower"]["2.0"]}
    assert shares[2] == pytest.approx(1.0)
    assert shares[3] == pytest.approx(0.149, abs=5e-4)
    assert shares[4] == pytest.approx(0.0164, abs=5e-4)
    assert shares[5] == pytest.approx(0.0018, abs=1e-4)


# -- the ladder ----------------------------------------------------------------
@pytest.mark.parametrize("delta", [2, 3, 5])
def test_ladder_density_is_normalized_with_the_sech_width(delta):
    """rho_Delta integrates to 1 and has <omega^2> = 2 Delta Jc^2, the
    curvature of sech^(2 Delta)(Jc t) at t = 0."""
    Jc = 0.7
    om = np.linspace(-80 * Jc, 80 * Jc, 40001)
    rho = h.ladder_density(om, Jc, delta)
    assert np.trapezoid(rho, om) == pytest.approx(1.0, rel=1e-8)
    assert np.trapezoid(om ** 2 * rho, om) == pytest.approx(2 * delta * Jc ** 2, rel=1e-6)


def test_ladder_poles():
    assert h.quasinormal_ladder(3, 2.0, 2) == [-6j, -8j, -10j]


@pytest.mark.parametrize("delta", [2, 3])
def test_exact_chord_density_converges_at_rate_delta_over_four(delta):
    """max |S_exact - rho_Delta| / max rho_Delta = (Delta/4) lambda + O(lambda^2)."""
    c = h.convergence(delta, lams=(0.1, 0.05), n_omega=60)
    r1, r2 = c["rows"]
    assert r2["dev"] < r1["dev"]
    assert r2["dev_over_lam"] == pytest.approx(delta / 4, rel=0.03)


def test_artifact_convergence_law(art):
    for c in art["convergence"]:
        last = c["rows"][-1]
        assert last["lam"] == 0.05
        assert last["dev_over_lam"] == pytest.approx(c["delta"] / 4, rel=0.03)
        devs = [r["dev"] for r in c["rows"]]
        assert devs == sorted(devs, reverse=True)


# -- the time domain -------------------------------------------------------------
def test_chord_decay_rate_plateaus_at_delta_a():
    for delta in (2, 3):
        r = h.decay_rates(delta, 0.03, a_taus=(6, 8))["rate_over_a"]
        assert r["8"] == pytest.approx(delta, abs=0.02)
        assert r["6"] == pytest.approx(delta, abs=0.03)


def test_massive_and_thooft_rates_grow_without_bound(art):
    for r in art["massive_rates"] + art["thooft_tower_rates"]:
        v = list(r["rate_over_a"].values())
        assert all(b > a for a, b in zip(v, v[1:]))
    for r in art["massive_rates"]:
        m = r["m_over_a"]
        # rate/a -> (M/a) sinh(a tau/2) + 1/4 (the 1/4 from K_0's z^(-1/2))
        assert r["rate_over_a"]["8"] == pytest.approx(m * math.sinh(4.0) + 0.25, rel=0.01)
    # every chord rate stays at or below its plateau
    for r in art["decay_rates"]:
        assert max(r["rate_over_a"].values()) < r["delta"] + 0.05


def test_thooft_tower_is_led_by_its_lightest_meson():
    t = h.thooft_tower_rates(0.5, n_levels=8)["rate_over_a"]
    m = h.massive_rates(0.5)["rate_over_a"]
    assert t["8"] == pytest.approx(m["8"], rel=0.02)


# -- the frequency domain ----------------------------------------------------------
def test_bessel_identity():
    for row in h.bessel_identity(((1.0, 0.7), (2.5, 3.0))):
        assert row["rel"] < 1e-12


def test_massive_density_matches_bessel_closed_form():
    """(1/pi) int_0^inf cos(omega tau) K_0((2M/a) cosh(a tau/2)) dtau
    = K_(i omega/a)(M/a)^2 / (pi a), here at a = 1."""
    import mpmath as mp
    m = 0.5
    om = np.array([0.0, 1.0, 2.5])
    S = h.massive_two_sided_density(om, m, a=1.0)
    for w, s in zip(om, S):
        closed = float(mp.re(mp.besselk(1j * w, m) ** 2)) / math.pi
        assert s == pytest.approx(closed, rel=1e-6)


def test_massive_densities_have_nodes_ladders_do_not(art):
    assert all(v >= 1 for v in art["nodes"]["massive"].values())
    assert all(v == 0 for v in art["nodes"]["ladder"].values())


def test_full_channel_fit_falls_below_every_massive_floor(art):
    fits = {(r["n"], r["lam"]): r["rel_residual"] for r in art["full_channel_fit"]}
    floor = min(r["rel_residual"] for r in art["massive_fit_control"])
    for n in (3, 4):
        series = [fits[(n, lam)] for lam in (0.4, 0.2, 0.1, 0.05)]
        assert series == sorted(series, reverse=True)
        assert fits[(n, 0.05)] < floor / 20


# -- units and the de Sitter reading -----------------------------------------------
@pytest.mark.parametrize("lam_paper", [1.0, 0.25, 0.05])
def test_ladder_spacing_is_string_scale(lam_paper):
    u = h.units_mss(lam_paper)
    assert u["scriptJ_over_JMSS"] == pytest.approx(math.sqrt(2.0), rel=1e-14)
    assert u["ladder_spacing_over_JMSS"] == pytest.approx(2 * math.sqrt(2.0), rel=1e-14)


def test_real_dimensions_above_one_are_tachyonic_in_ds():
    for r in h.ds_reading((2, 3, 5)):
        assert r["dS2_M2l2"] < 0 and r["NV_dS3_m2"] < 0
    assert h.ds_reading((2,))[0]["dS2_M2l2"] == -2
