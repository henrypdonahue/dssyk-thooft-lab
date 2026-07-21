#!/usr/bin/env python3
"""Asserting tests for the baryon-sector machinery (baryons.py, rung 23).

Every closed form is checked against an independent route: sector blocks vs
the dense builder AND vs the sparse Jordan-Wigner reference; the correlator
collapse vs full ED that knows nothing about it; the sector-m2 closed form vs
a structurally different enumeration and vs ensemble ED means; the DeltaE
statistics vs the Gaussian-sum prediction.

Run with:  pytest -q -m 'not slow'   (two heavier statistical anchors are slow)
"""

from math import comb

import numpy as np
import pytest

from dirac import (apply_monomial, coupling_variance_dirac,
                   dirac_hamiltonian, dirac_hamiltonian_reference,
                   draw_couplings)
from baryons import (baryon_correlator_ed, baryon_dagger, baryon_sign,
                     charging_curve, dense_sector_blocks,
                     deltaE_from_couplings, ensemble_run, fit_mass_models,
                     instance_observables, power_exponent_fit,
                     reassemble_dense, sector_hamiltonian, sector_m2_enumerate,
                     sector_m2_closed_form, sector_states)


# --------------------------------------------------------------------------
# 1. charge-sector machinery
# --------------------------------------------------------------------------
@pytest.mark.parametrize("ensemble", ["generic", "c_symmetric"])
def test_sector_reassembly_is_dense_h(ensemble):
    """Scattering the independently assembled sector blocks back into the
    full space must reproduce dirac_hamiltonian EXACTLY (same per-element
    accumulation order), and H must vanish identically between sectors."""
    Nc, p = 6, 4
    rng = np.random.default_rng(11)
    coup = draw_couplings(Nc, p, rng, ensemble)
    H = dirac_hamiltonian(Nc, p, couplings=coup)
    blocks = [sector_hamiltonian(Nc, coup, q) for q in range(Nc + 1)]
    assert [b.shape[0] for b in blocks] == [comb(Nc, q) for q in range(Nc + 1)]
    R = reassemble_dense(Nc, blocks)
    assert np.array_equal(R, H)          # bitwise: same terms, same order
    # block-diagonality: off-sector entries of the dense H are exact zeros
    qs = np.bitwise_count(np.arange(1 << Nc, dtype=np.int64))
    off = qs[:, None] != qs[None, :]
    assert np.abs(H[off]).max() == 0.0


def test_sector_spectrum_reassembles_dense_spectrum():
    """Concatenated sector eigenvalues = dense spectrum (independent LAPACK
    paths, so machine precision rather than bitwise)."""
    Nc, p = 7, 4
    rng = np.random.default_rng(12)
    coup = draw_couplings(Nc, p, rng, "generic")
    H = dirac_hamiltonian(Nc, p, couplings=coup)
    blocks = dense_sector_blocks(H, Nc)
    ev_b = np.sort(np.concatenate([np.linalg.eigvalsh(b) for b in blocks]))
    ev_d = np.linalg.eigvalsh(H)
    scale = max(1.0, np.abs(ev_d).max())
    assert np.abs(ev_b - ev_d).max() < 1e-11 * scale


def test_sector_blocks_match_sparse_jw_reference():
    """Sector blocks vs the sliced dirac_hamiltonian_reference (sparse
    Jordan-Wigner kron route; shares no assembly code)."""
    Nc, p = 5, 4
    rng = np.random.default_rng(13)
    coup = draw_couplings(Nc, p, rng, "generic")
    Href = dirac_hamiltonian_reference(Nc, p, coup)
    for q in range(Nc + 1):
        st = sector_states(Nc, q)
        blk = sector_hamiltonian(Nc, coup, q)
        assert np.abs(blk - Href[np.ix_(st, st)]).max() < 1e-12


@pytest.mark.parametrize("ensemble", ["generic", "c_symmetric"])
def test_edge_blocks_vanish_identically(ensemble):
    """q < p/2 blocks are identically zero (every term needs p/2 occupied
    modes); for c_symmetric (disjoint A, B) the mirror blocks q > Nc - p/2
    vanish too, so E0(Nc) = E0(Nc-1) = 0 exactly."""
    Nc, p = 6, 4
    k = p // 2
    rng = np.random.default_rng(14)
    coup = draw_couplings(Nc, p, rng, ensemble)
    blocks = [sector_hamiltonian(Nc, coup, q) for q in range(Nc + 1)]
    zero_qs = set(range(k))
    if ensemble == "c_symmetric":
        zero_qs |= set(range(Nc - k + 1, Nc + 1))
    for q in zero_qs:
        assert np.abs(blocks[q]).max() == 0.0
    for q in set(range(Nc + 1)) - zero_qs:
        assert np.abs(blocks[q]).max() > 0.0


# --------------------------------------------------------------------------
# 2. the baryon operator and the exact collapse
# --------------------------------------------------------------------------
def test_baryon_operator_structure():
    """B^dag = c^dag_0 ... c^dag_{Nc-1} has exactly one nonzero entry,
    <full|B^dag|empty> = +1 (each creation parity counts only lower-index
    modes, all empty when it acts)."""
    for Nc in (4, 5, 6):
        Bd = baryon_dagger(Nc)
        nz = np.nonzero(Bd)
        assert len(nz[0]) == 1
        assert nz[0][0] == (1 << Nc) - 1 and nz[1][0] == 0
        assert Bd[(1 << Nc) - 1, 0] == 1.0


def test_full_state_number_operator_sign():
    """<s|c^dag_A c_A|s> = (-1)^{k(k-1)/2} for A subset of occ(s): the JW
    lemma behind the DeltaE closed form, for k = 2, 3, 4."""
    for k in (2, 3, 4):
        Nc = 2 * k + 1
        full = np.array([(1 << Nc) - 1], dtype=np.int64)
        for A in [tuple(range(k)), tuple(range(1, k + 1)),
                  tuple(range(Nc - k, Nc))]:
            valid, rows, signs = apply_monomial(Nc, A, A, full)
            assert valid[0] and rows[0] == full[0]
            assert signs[0] == baryon_sign(k)


@pytest.mark.parametrize("ensemble", ["generic", "c_symmetric"])
def test_baryon_correlator_collapse(ensemble):
    """Tr(B(t)B^dag)/dim = exp(-i DeltaE t)/dim with DeltaE =
    (-1)^{k(k-1)/2} sum_A J_AA -- asserted against full ED that uses nothing
    about the collapse.  |G| = 1/dim identically: no decay, and the one
    frequency is O(Nc^{-1/2}) (generic) or exactly zero (c_symmetric); the
    T = infinity baryon 2-pt function is blind to the extensive mass."""
    Nc, p = 6, 4
    dim = 1 << Nc
    rng = np.random.default_rng(15)
    coup = draw_couplings(Nc, p, rng, ensemble)
    H = dirac_hamiltonian(Nc, p, couplings=coup)
    dE = deltaE_from_couplings(coup, p // 2)
    # the two 1-dim sectors: diagonal energies straight from the dense H
    assert H[0, 0] == 0.0                       # <empty|H|empty> exact zero
    assert abs(H[-1, -1].real - dE) < 1e-12
    if ensemble == "c_symmetric":
        assert dE == 0.0                        # no diagonal couplings at all
        assert H[-1, -1] == 0.0
    ts = np.linspace(0.0, 9.0, 31)
    G = baryon_correlator_ed(H, ts)
    assert np.abs(G - np.exp(-1j * dE * ts) / dim).max() < 1e-12
    assert np.abs(np.abs(G) * dim - 1.0).max() < 1e-12


def test_deltaE_gaussian_statistics():
    """DeltaE is a sum of C(Nc,k) i.i.d. N(0, sigma^2) diagonal couplings:
    mean 0, Var = C(Nc,k) sigma^2.  Checked against ED diagonals over 300
    draws (variance-of-variance sampling error ~ sqrt(2/(n-1)))."""
    Nc, p, n = 5, 4, 300
    k = p // 2
    var_exact = comb(Nc, k) * coupling_variance_dirac(Nc, p)
    dEs = []
    for j in range(n):
        rng = np.random.default_rng(1000 + j)
        coup = draw_couplings(Nc, p, rng, "generic")
        H = dirac_hamiltonian(Nc, p, couplings=coup)
        dE_ed = float(H[-1, -1].real - H[0, 0].real)
        assert abs(dE_ed - deltaE_from_couplings(coup, k)) < 1e-12
        dEs.append(dE_ed)
    dEs = np.array(dEs)
    assert abs(dEs.mean()) < 5 * np.sqrt(var_exact / n)
    relerr = np.sqrt(2.0 / (n - 1))
    assert abs(np.var(dEs, ddof=1) / var_exact - 1.0) < 5 * relerr


# --------------------------------------------------------------------------
# 3. sector second moments: closed form, two routes + ED ensemble
# --------------------------------------------------------------------------
@pytest.mark.parametrize("Nc,p,ensemble", [
    (4, 4, "generic"), (6, 4, "generic"), (8, 4, "generic"),
    (7, 6, "generic"), (4, 4, "c_symmetric"), (6, 4, "c_symmetric"),
    (8, 4, "c_symmetric")])
def test_sector_m2_two_routes_agree(Nc, p, ensemble):
    """Key-table enumeration vs overlap-resolved combinatorial closed form:
    independent code paths, exact agreement."""
    for q in range(Nc + 1):
        a = sector_m2_enumerate(Nc, p, q, ensemble)
        b = sector_m2_closed_form(Nc, p, q, ensemble)
        assert abs(a - b) < 1e-12 * max(1.0, abs(b))


@pytest.mark.parametrize("ensemble", ["generic", "c_symmetric"])
def test_sector_m2_matches_ed_ensemble_mean(ensemble):
    """Disorder mean of the ED Tr_q H^2/dim_q vs the closed form, within
    5 SEM per sector (cross terms vanish only in the mean, hence
    statistical); identically-zero sectors must be exactly zero."""
    Nc, p, n = 6, 4, 120
    m2 = []
    for j in range(n):
        d = instance_observables(Nc, p, 2000 + j, ensemble)
        m2.append(d["m2q"])
    m2 = np.array(m2)
    mean = m2.mean(axis=0)
    sem = m2.std(axis=0, ddof=1) / np.sqrt(n)
    for q in range(Nc + 1):
        pred = sector_m2_closed_form(Nc, p, q, ensemble)
        if pred == 0.0:
            assert mean[q] == 0.0
        else:
            assert abs(mean[q] - pred) < 5 * sem[q]


# --------------------------------------------------------------------------
# 4. charging curve: symmetries and observables
# --------------------------------------------------------------------------
def test_charging_curve_particle_hole_palindrome():
    """c_symmetric instances are exactly C-invariant, and C maps sector q to
    Nc - q, so E0(q) = E0(Nc - q) per instance; also Delta2(Nc-1) = E0(2)
    exactly (edge blocks vanish)."""
    Nc, p = 6, 4
    for seed in (21, 22, 23):
        d = instance_observables(Nc, p, seed, "c_symmetric")
        E0 = np.array(d["E0"])
        assert np.abs(E0 - E0[::-1]).max() < 1e-10
        assert E0[Nc] == 0.0 and E0[Nc - 1] == 0.0
        assert abs(d["delta2"] - E0[2]) < 1e-12
        assert d["delta2"] <= 0.0        # exactly concave edge


def test_low_charge_quarks_are_free():
    """E0(q) = 0 exactly for q < p/2 in BOTH ensembles (blocks vanish)."""
    Nc, p = 6, 4
    for ens in ("generic", "c_symmetric"):
        d = instance_observables(Nc, p, 31, ens)
        for q in range(p // 2):
            assert d["E0"][q] == 0.0


def test_mb_positive_and_extensive_trend():
    """M_B > 0 by construction (E0(Nc) >= global min); its mean grows
    strongly with Nc (the extensivity the fits quantify).  Loose statistical
    check: mean M_B at Nc = 8 exceeds that at Nc = 4 by > 8 sigma."""
    r4 = ensemble_run(4, 4, 100, "generic", seed0=900)
    r8 = ensemble_run(8, 4, 100, "generic", seed0=901)
    assert r4["MB_mean"] > 0 and r8["MB_mean"] > 0
    gap_sigma = np.hypot(r4["MB_sem"], r8["MB_sem"])
    assert r8["MB_mean"] - r4["MB_mean"] > 8 * gap_sigma


# --------------------------------------------------------------------------
# 5. fit machinery bookkeeping
# --------------------------------------------------------------------------
def test_fit_models_on_synthetic_data():
    """On exact y = 3 Nc all three models fit perfectly; AIC then ranks by
    parameter count (a*Nc wins, Delta-AIC = 2 to each 2-parameter model),
    and the power-law recovers gamma = 1."""
    import warnings
    from scipy.optimize import OptimizeWarning
    ns = np.array([4, 5, 6, 7, 8, 9, 10], float)
    ys = 3.0 * ns
    errs = 0.01 * np.ones_like(ns)
    with warnings.catch_warnings():
        # exact synthetic data -> zero residual -> degenerate covariance;
        # benign here, the test only uses the best-fit values and AICs
        warnings.simplefilter("ignore", OptimizeWarning)
        f = fit_mass_models(ns, ys, errs)
    assert abs(f["a*Nc"]["params"][0] - 3.0) < 1e-8
    assert abs(f["a*Nc^gamma"]["params"][1] - 1.0) < 1e-6
    assert f["a*Nc"]["aic"] < f["a*Nc + b"]["aic"]
    assert abs(f["a*Nc + b"]["aic"] - f["a*Nc"]["aic"] - 2.0) < 1e-6
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        p = power_exponent_fit(ns, 5.0 * ns ** -1.5, 0.001 * ns ** -1.5)
    assert abs(p["e"] + 1.5) < 1e-6


# --------------------------------------------------------------------------
# slow statistical anchors at larger Nc
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_sector_m2_ed_anchor_nc8_slow():
    """Heavier anchor: Nc = 8 ensemble means of Tr_q H^2/dim_q vs the closed
    form, both ensembles, 200 draws, 5 SEM."""
    Nc, p, n = 8, 4, 200
    for ens in ("generic", "c_symmetric"):
        m2 = []
        for j in range(n):
            m2.append(instance_observables(Nc, p, 3000 + j, ens)["m2q"])
        m2 = np.array(m2)
        mean = m2.mean(axis=0)
        sem = m2.std(axis=0, ddof=1) / np.sqrt(n)
        for q in range(Nc + 1):
            pred = sector_m2_closed_form(Nc, p, q, ens)
            if pred == 0.0:
                assert mean[q] == 0.0
            else:
                assert abs(mean[q] - pred) < 5 * sem[q]


@pytest.mark.slow
def test_deltaE_variance_anchor_nc8_slow():
    """DeltaE variance vs C(Nc,2) sigma^2 at Nc = 8, 400 draws."""
    Nc, p, n = 8, 4, 400
    var_exact = comb(Nc, 2) * coupling_variance_dirac(Nc, p)
    dEs = []
    for j in range(n):
        rng = np.random.default_rng(4000 + j)
        coup = draw_couplings(Nc, p, rng, "generic")
        dEs.append(deltaE_from_couplings(coup, 2))
    ratio = np.var(np.array(dEs), ddof=1) / var_exact
    assert abs(ratio - 1.0) < 5 * np.sqrt(2.0 / (n - 1))
