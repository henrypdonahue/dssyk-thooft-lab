#!/usr/bin/env python3
"""Asserting tests for the Dirac/QED campaign (qed_campaign.py, rung 19).

Every closed form is checked against an independent enumeration route
(sparse-operator traces + literal Isserlis Wick) and against ED disorder
ensembles; the spectral pipeline is checked against eigenbasis-free
commutator identities.  Run with:  pytest -q  (the Nc = 10 ED ensemble is
marked slow).
"""

import json
from math import comb
from pathlib import Path

import numpy as np
import pytest

from dirac import (annihilators, coupling_variance_dirac, dirac_hamiltonian,
                   draw_couplings, mn_tower)
from qed_campaign import (charging_curve, enumerated_m2_stats,
                          enumerated_mu2_charged, fit_charging,
                          instance_spectral, m2_of_hamiltonian,
                          mean_m2_exact, mean_mu2_charged_exact,
                          measure_m2_ensemble, measure_spectral,
                          monomial_matrix, mu2_charged_by_commutator,
                          relvar_m2_exact_dirac, relvar_m2_leading,
                          sector_indices, sector_spectra, t_ab_exact,
                          trace_pair_matrix, stored_keys, var_m2_exact)


# --------------------------------------------------------------------------
# 1. m2 disorder statistics: occupation counting, quadratic form, Wick
# --------------------------------------------------------------------------
def test_t_ab_occupation_counting():
    """T_AB = 2^{-|A u B|} against literal traces of sparse-operator products,
    covering disjoint, overlapping, diagonal, and k = 1 keys."""
    Nc = 6
    dim = 1 << Nc
    for (A, B) in [((0, 1), (2, 3)), ((0, 1), (1, 2)), ((0, 2), (0, 2)),
                   ((1, 3), (3, 5)), ((0,), (1,)), ((2,), (2,))]:
        O = monomial_matrix(Nc, A, B)
        t = float(np.trace(O @ O.conj().T).real) / dim
        assert abs(t - t_ab_exact(A, B)) < 1e-12


def test_instance_m2_quadratic_form_generic():
    """For a generic instance, Tr(H^2)/dim equals the full quadratic form
    sum J_a J_b T_ab machine-exactly, while the naive sum |J|^2 T_AB does
    NOT (admissible cross-pairings): the structural fact of the docstring."""
    Nc, p = 5, 4
    J = draw_couplings(Nc, p, np.random.default_rng(1), "generic")
    H = dirac_hamiltonian(Nc, p, couplings=J)
    m2 = m2_of_hamiltonian(H)
    keys = list(J.keys())
    T = trace_pair_matrix(Nc, keys)
    vec = np.array([J[key] for key in keys])
    quad = float((vec @ T @ vec).real)
    naive = sum(abs(v) ** 2 * t_ab_exact(A, B) for (A, B), v in J.items())
    assert abs(m2 - quad) < 1e-10 * m2
    assert abs(m2 - naive) > 1e-3 * m2     # the cross terms are real


def test_instance_m2_sum_of_squares_c_symmetric():
    """For the disjoint/real ensemble the naive sum IS exact: T is the pure
    pairing delta, so m2 = sum |J|^2 2^{-p} instance by instance."""
    Nc, p = 6, 4
    J = draw_couplings(Nc, p, np.random.default_rng(2), "c_symmetric")
    H = dirac_hamiltonian(Nc, p, couplings=J)
    naive = sum(abs(v) ** 2 * t_ab_exact(A, B) for (A, B), v in J.items())
    assert abs(m2_of_hamiltonian(H) - naive) < 1e-12 * naive
    # and the pairing-delta structure itself
    keys = list(J.keys())
    T = trace_pair_matrix(Nc, keys)
    idx = {key: i for i, key in enumerate(keys)}
    for i, (A, B) in enumerate(keys):
        row = T[i].copy()
        assert abs(row[idx[(B, A)]] - 2.0 ** -p) < 1e-13
        row[idx[(B, A)]] = 0.0
        assert np.abs(row).max() < 1e-13


@pytest.mark.parametrize("Nc,p", [(4, 4), (5, 4), (6, 4), (4, 2), (6, 2)])
def test_m2_closed_forms_match_enumeration_generic(Nc, p):
    """E[m2] and Var(m2) closed forms vs literal enumeration (numeric traces
    + full Isserlis Wick; no shared assembly code)."""
    me, ve = enumerated_m2_stats(Nc, p, "generic")
    assert abs(me / mean_m2_exact(Nc, p, "generic") - 1.0) < 1e-12
    assert abs(ve / var_m2_exact(Nc, p, "generic") - 1.0) < 1e-12


@pytest.mark.parametrize("Nc", [4, 5, 6])
def test_m2_closed_forms_match_enumeration_c_symmetric(Nc):
    me, ve = enumerated_m2_stats(Nc, 4, "c_symmetric")
    assert abs(me / mean_m2_exact(Nc, 4, "c_symmetric") - 1.0) < 1e-12
    assert abs(ve / var_m2_exact(Nc, 4, "c_symmetric") - 1.0) < 1e-12


def test_relvar_csym_is_two_over_n_pairs():
    """The c_symmetric relVar is exactly 2/N_pairs -- the same identity as
    the Majorana 2/C(N,p), with N_pairs = C(Nc,k)C(Nc-k,k)/2."""
    for Nc in (6, 8, 10, 12):
        n_pairs = comb(Nc, 2) * comb(Nc - 2, 2) // 2
        assert abs(relvar_m2_exact_dirac(Nc, 4, "c_symmetric")
                   - 2.0 / n_pairs) < 1e-15


def test_relvar_leading_rates():
    """Exact relVar approaches the stated large-Nc rates from the correct
    side and within 25% by Nc = 40 (p = 4)."""
    for ens in ("generic", "c_symmetric"):
        r = [relvar_m2_exact_dirac(Nc, 4, ens) / relvar_m2_leading(Nc, 4, ens)
             for Nc in (20, 40, 80, 160)]
        assert abs(r[1] - 1.0) < 0.25
        assert abs(r[-1] - 1.0) < 0.06
        # monotone approach to 1 in this window
        assert all(abs(b - 1.0) < abs(a - 1.0) for a, b in zip(r, r[1:]))


def test_m2_relvar_matches_ED_Nc8():
    """ED disorder ensembles at Nc = 8 reproduce the exact mean and relVar
    within sampling error (250 draws => ~9% on relVar)."""
    for ens in ("generic", "c_symmetric"):
        r = measure_m2_ensemble(8, 4, ens, n_inst=250, seed=17)
        assert abs(r["m2_mean"] - r["mean_exact"]) < 5 * r["m2_sem"]
        assert abs(r["relvar"] - r["relvar_exact"]) < 4 * r["relvar_err"]


@pytest.mark.slow
def test_m2_relvar_matches_ED_Nc10():
    for ens in ("generic", "c_symmetric"):
        r = measure_m2_ensemble(10, 4, ens, n_inst=300, seed=19)
        assert abs(r["m2_mean"] - r["mean_exact"]) < 5 * r["m2_sem"]
        assert abs(r["relvar"] - r["relvar_exact"]) < 4 * r["relvar_err"]


# --------------------------------------------------------------------------
# 2. charge sectors and the charging curve
# --------------------------------------------------------------------------
def test_sector_reassembly_machine_exact():
    """Sector-resolved spectra reassemble the dense spectrum machine-exactly,
    and every off-sector block of H vanishes identically."""
    Nc, p = 8, 4
    H = dirac_hamiltonian(Nc, p, rng=np.random.default_rng(3))
    sectors = sector_indices(Nc)
    assert sum(len(ix) for ix in sectors) == 1 << Nc
    eigs = np.sort(np.concatenate(sector_spectra(H, Nc, sectors)))
    dense = np.linalg.eigvalsh(H)
    assert np.abs(eigs - dense).max() < 1e-10
    for qa in range(Nc + 1):
        for qb in range(qa + 1, Nc + 1):
            block = H[np.ix_(sectors[qa], sectors[qb])]
            assert np.abs(block).max() == 0.0


def test_edge_sectors_and_particle_hole():
    """Exact sector facts: H_q = 0 for q < p/2 (so E0 = 0 there); the
    c_symmetric top sector vanishes and E0(q) = E0(Nc - q) instance-exactly;
    the generic top sector is the 1-dim sum of diagonal couplings."""
    Nc, p = 6, 4
    sectors = sector_indices(Nc)
    Jc = draw_couplings(Nc, p, np.random.default_rng(4), "c_symmetric")
    Hc = dirac_hamiltonian(Nc, p, couplings=Jc)
    spec_c = sector_spectra(Hc, Nc, sectors)
    for q in range(p // 2):
        assert np.abs(Hc[np.ix_(sectors[q], sectors[q])]).max() == 0.0
    assert np.abs(spec_c[Nc]).max() < 1e-14
    for q in range(Nc + 1):
        assert np.abs(spec_c[q] - spec_c[Nc - q]).max() < 1e-10
    Jg = draw_couplings(Nc, p, np.random.default_rng(5), "generic")
    Hg = dirac_hamiltonian(Nc, p, couplings=Jg)
    top = Hg[np.ix_(sectors[Nc], sectors[Nc])]
    assert top.shape == (1, 1)
    assert abs(top[0, 0]) > 1e-8      # J_AA break particle-hole symmetry


def test_charging_curve_shapes_and_edges():
    """Measured curve carries SEMs; the q < p/2 sectors are exactly zero for
    every instance (zero mean AND zero SEM)."""
    curve = charging_curve(6, 4, "generic", n_inst=8, seed=11)
    assert len(curve["e0_mean"]) == 7 and len(curve["e0_sem"]) == 7
    for q in range(2):
        assert curve["e0_mean"][q] == 0.0 and curve["e0_sem"][q] == 0.0
    assert all(s >= 0 for s in curve["e0_sem"])


def test_charging_fit_recovers_planted_models():
    """The fit machinery recovers a planted quadratic (capacitive) curve --
    coefficient and AIC preference -- and symmetrically for a planted |x|."""
    Nc, p = 10, 4
    q = np.arange(Nc + 1, dtype=float)
    x = q - Nc / 2.0
    for planted, expect in [(0.125 * x ** 2, "capacitive_x2"),
                            (0.7 * np.abs(x), "confining_absx")]:
        curve = dict(Nc=Nc, p=p, q=list(range(Nc + 1)),
                     e0_mean=list(-3.0 + 0.02 * x + planted),
                     e0_sem=[0.01] * (Nc + 1))
        f = fit_charging(curve)
        assert f["verdict"].startswith("best: " + expect)
        if expect == "capacitive_x2":
            assert abs(f["capacitive_x2"]["coef"][2] - 0.125) < 1e-9
            assert abs(f["capacitance"] - 4.0) < 1e-6
        else:
            assert abs(f["confining_absx"]["coef"][2] - 0.7) < 1e-9


# --------------------------------------------------------------------------
# 3. charged vs singlet spectral weight
# --------------------------------------------------------------------------
def test_mu0_charged_exact_half():
    """mu0 of the mode-averaged symmetrized c_i weight is exactly 1/2."""
    for ens in ("generic", "c_symmetric"):
        r = instance_spectral(6, 4, ens, seed=7)
        assert abs(r["charged"]["mu0"] - 0.5) < 1e-10


def test_mu2_charged_commutator_identity():
    """Eigenbasis mu2 of the charged channel equals the eigenbasis-free
    identity mu2 = mean_i ||[H, c_i]||_F^2/dim per instance."""
    Nc, p, seed = 6, 4, 8
    for ens in ("generic", "c_symmetric"):
        r = instance_spectral(Nc, p, ens, seed=seed)
        rng = np.random.default_rng(seed)
        H = dirac_hamiltonian(Nc, p,
                              couplings=draw_couplings(Nc, p, rng, ens))
        mu2_ind = mu2_charged_by_commutator(H, annihilators(Nc))
        assert abs(r["charged"]["mu2"] - mu2_ind) < 1e-8 * mu2_ind


def test_mu2_singlet_commutator_identity():
    """Eigenbasis singlet mu2 equals ||[H, M2bar]||_F^2/dim with M2bar built
    independently in the position basis via dirac.mn_tower -- validating the
    whole eigenbasis M2 assembly (cdot~ = i omega o c~)."""
    Nc, p, seed = 6, 4, 9
    dim = 1 << Nc
    for ens in ("generic", "c_symmetric"):
        r = instance_spectral(Nc, p, ens, seed=seed)
        rng = np.random.default_rng(seed)
        H = dirac_hamiltonian(Nc, p,
                              couplings=draw_couplings(Nc, p, rng, ens))
        M2 = mn_tower(H, annihilators(Nc), 2)[2]
        M2bar = M2 - (np.trace(M2) / dim) * np.eye(dim)
        assert np.abs(M2bar - M2bar.conj().T).max() < 1e-10
        comm = H @ M2bar - M2bar @ H
        mu2_ind = float(np.sum(np.abs(comm) ** 2)) / dim
        assert abs(r["singlet"]["mu2"] - mu2_ind) < 1e-8 * mu2_ind
        # mu0 consistency too: ||M2bar||_F^2/dim
        mu0_ind = float(np.sum(np.abs(M2bar) ** 2)) / dim
        assert abs(r["singlet"]["mu0"] - mu0_ind) < 1e-8 * mu0_ind


def test_w2_charged_instance_exact_c_symmetric():
    """For the disjoint/real ensemble, mu2_charged = (p/Nc) m2 holds PER
    INSTANCE (so w2_charged = 2p/Nc exactly): the disjoint-key pattern
    argument kills the cross terms of the commutator quadratic form.  The
    generic ensemble must violate it (cross terms are real)."""
    Nc, p = 6, 4
    for seed in (0, 1):
        H = dirac_hamiltonian(Nc, p, couplings=draw_couplings(
            Nc, p, np.random.default_rng(seed), "c_symmetric"))
        mu2 = mu2_charged_by_commutator(H, annihilators(Nc))
        m2 = m2_of_hamiltonian(H)
        assert abs(mu2 - (p / Nc) * m2) < 1e-12 * mu2
    Hg = dirac_hamiltonian(Nc, p, couplings=draw_couplings(
        Nc, p, np.random.default_rng(0), "generic"))
    mu2g = mu2_charged_by_commutator(Hg, annihilators(Nc))
    assert abs(mu2g - (p / Nc) * m2_of_hamiltonian(Hg)) > 0.05 * mu2g


def test_mean_mu2_charged_closed_vs_enumeration():
    """The closed form for E[mu2_charged] against the literal key-by-key
    commutator sum (independent matrices route)."""
    for (Nc, p, ens) in [(5, 4, "generic"), (4, 2, "generic"),
                         (5, 4, "c_symmetric"), (6, 4, "c_symmetric")]:
        e = mean_mu2_charged_exact(Nc, p, ens)
        n = enumerated_mu2_charged(Nc, p, ens)
        assert abs(n / e - 1.0) < 1e-12


def test_spectral_ensemble_structure():
    """Ensemble spectral run: odd moments vanish (symmetrized weights),
    kurtosis obeys Cauchy-Schwarz, mu0 is exact, and the measured E[mu2]
    of the charged channel matches the closed form within error bars."""
    r = measure_spectral(6, 4, "generic", n_inst=16, seed0=100)
    assert r["odd_moment_check"] < 1e-8
    assert abs(r["mu0_charged"] - 0.5) < 1e-10
    assert r["kurt_charged"] >= 1.0 and r["kurt_singlet"] >= 1.0
    assert r["w2_charged"] > 0 and r["w2_singlet"] > 0
    pull = abs(r["mu2_charged_mean"] - r["mu2_charged_exact"])
    assert pull < 5 * r["mu2_charged_sem"]


def test_qed_campaign_json_structure():
    """The committed campaign JSON has the advertised structure (skipped if
    the campaign has not been run in this checkout)."""
    path = Path(__file__).resolve().parent / "qed_campaign.json"
    if not path.exists():
        pytest.skip("qed_campaign.json not present; run qed_campaign.py")
    data = json.loads(path.read_text())
    for key in ("citation_context", "m2_exact_table", "m2_ensembles",
                "charging", "spectral"):
        assert key in data
    for row in data["m2_ensembles"]:
        assert abs(row["relvar"] - row["relvar_exact"]) < 4 * row["relvar_err"]
    for curve in data["charging"]:
        assert "fits" in curve and "verdict" in curve["fits"]
        assert len(curve["e0_mean"]) == curve["Nc"] + 1
    for row in data["spectral"]:
        assert row["w2_charged_sem"] > 0 and row["w2_singlet_sem"] > 0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
