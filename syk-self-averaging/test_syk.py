#!/usr/bin/env python3
"""Asserting tests for the SYK harness and the self-averaging measurement.

Run with:  pytest -q     (the RMT test is marked slow; deselect with -m 'not slow')
"""

from math import comb

import numpy as np
import pytest
import scipy.sparse as sp

from syk import (majorana_operators, majorana_hamiltonian, fermion_parity,
                 charge_operator, coupling_variance)
from pauli_strings import (majorana_string, monomial_string,
                           majorana_hamiltonian_fast, add_string_dense)
from dirac import (draw_couplings, dirac_hamiltonian,
                   dirac_hamiltonian_reference, charge_matrix, annihilators,
                   mn_operator, time_derivative, particle_hole, cp_conjugate,
                   wrong_channel_weight)
from self_averaging import relvar_m2_exact, measure_moments
import validate_syk


def test_clifford_algebra():
    """{psi_a, psi_b} = 2 delta_ab I, psi_a Hermitian."""
    N = 8
    psis = majorana_operators(N)
    dim = 1 << (N // 2)
    I = sp.identity(dim, dtype=complex)
    for a in range(N):
        assert abs((psis[a] - psis[a].getH()).toarray()).max() < 1e-12
        for b in range(N):
            anti = psis[a] @ psis[b] + psis[b] @ psis[a]
            target = (2.0 if a == b else 0.0) * I
            assert abs((anti - target).toarray()).max() < 1e-12


@pytest.mark.parametrize("p", [4, 6])
def test_hamiltonian_hermitian_and_parity(p):
    N = 8
    rng = np.random.default_rng(0)
    H = majorana_hamiltonian(N, p, rng)
    assert abs((H - H.getH()).toarray()).max() < 1e-12
    P = fermion_parity(N)
    assert abs((P @ P - sp.identity(H.shape[0])).toarray()).max() < 1e-12
    assert abs((H @ P - P @ H).toarray()).max() < 1e-12


def test_m2_equals_sum_J_squared():
    """Tr(H^2)/dim must equal sum_c J_c^2 exactly (the analytic anchor)."""
    from itertools import combinations
    N, p = 10, 4
    psis = majorana_operators(N)
    dim = 1 << (N // 2)
    rng = np.random.default_rng(3)
    sigma = np.sqrt(coupling_variance(N, p))
    combos = list(combinations(range(N), p))
    J = rng.normal(0, sigma, len(combos))
    prefac = (1j) ** (p // 2)
    H = sp.csr_matrix((dim, dim), dtype=complex)
    for c, Jc in zip(combos, J):
        t = psis[c[0]]
        for a in c[1:]:
            t = t @ psis[a]
        H = H + (prefac * Jc) * t
    m2 = (H @ H).diagonal().sum().real / dim
    assert abs(m2 - np.sum(J ** 2)) < 1e-9


# --------------------------------------------------------------------------
# Fast Pauli-string assembly (pauli_strings.py)
# --------------------------------------------------------------------------
def test_majorana_strings_match_sparse():
    """Every psi_a as a (phase, x, z) string must equal the sparse JW matrix."""
    N = 8
    psis = majorana_operators(N)
    dim = 1 << (N // 2)
    states = np.arange(dim, dtype=np.int64)
    for a in range(N):
        M = np.zeros((dim, dim), dtype=complex)
        add_string_dense(M, majorana_string(N, a), 1.0, states)
        assert np.abs(M - psis[a].toarray()).max() < 1e-14


def test_monomial_string_matches_sparse_product():
    """A full p-Majorana monomial string must equal the sparse matrix product."""
    N = 8
    psis = majorana_operators(N)
    dim = 1 << (N // 2)
    states = np.arange(dim, dtype=np.int64)
    for combo in [(0, 1, 2, 3), (0, 2, 5, 7), (1, 3, 4, 6), (2, 3, 6, 7)]:
        ref = psis[combo[0]].toarray()
        for a in combo[1:]:
            ref = ref @ psis[a].toarray()
        M = np.zeros((dim, dim), dtype=complex)
        add_string_dense(M, monomial_string(N, combo), 1.0, states)
        assert np.abs(M - ref).max() < 1e-14


@pytest.mark.parametrize("N,p", [(8, 4), (10, 4), (8, 6)])
def test_fast_hamiltonian_equals_reference(N, p):
    """Identical rng state -> the fast builder reproduces syk.majorana_hamiltonian."""
    H_ref = majorana_hamiltonian(N, p, np.random.default_rng(5)).toarray()
    H_fast = majorana_hamiltonian_fast(N, p, np.random.default_rng(5))
    assert np.abs(H_ref - H_fast).max() < 1e-13


# --------------------------------------------------------------------------
# Dirac sector (dirac.py): the U(1)/"QED" half of the duality
# --------------------------------------------------------------------------
def test_dirac_builder_matches_sparse_reference():
    """The bitwise Dirac builder must equal the independent sparse construction."""
    Nc, p = 4, 4
    couplings = draw_couplings(Nc, p, np.random.default_rng(2))
    H = dirac_hamiltonian(Nc, p, couplings=couplings)
    H_ref = dirac_hamiltonian_reference(Nc, p, couplings)
    assert np.abs(H - H_ref).max() < 1e-13


def test_dirac_charge_conserved():
    """[H, Q] = 0 for the complex SYK model (the U(1)/EM sector), with H
    genuinely built and Hermitian, and Q with integer spectrum 0..Nc."""
    Nc, p = 5, 4
    H = dirac_hamiltonian(Nc, p, rng=np.random.default_rng(3))
    Q = charge_matrix(Nc)
    assert np.abs(H - H.conj().T).max() < 1e-12
    assert np.abs(H @ Q - Q @ H).max() < 1e-12
    eigs = np.linalg.eigvalsh(Q)
    assert np.allclose(np.round(eigs), eigs, atol=1e-9)
    assert eigs.min() >= -1e-9 and eigs.max() <= Nc + 1e-9
    # and the sparse-operator Q agrees with the diagonal one
    assert np.abs(charge_operator(Nc).toarray() - Q).max() < 1e-12


def test_mn_operator_identities():
    """The corrected tower identities (README correction #3 and its family):
    M_0 = Q;  M_1 = -(p/2) i H  (paper's M_1 = H fails at O(1));
    M_2 = -sum cdot^dag cdot  (paper's + sign is wrong)."""
    Nc, p = 5, 4
    H = dirac_hamiltonian(Nc, p, rng=np.random.default_rng(7))
    cs = annihilators(Nc)
    Q = charge_matrix(Nc)

    M0 = mn_operator(H, cs, 0)
    assert np.abs(M0 - Q).max() < 1e-12

    M1 = mn_operator(H, cs, 1)
    assert np.abs(M1 + (p / 2) * 1j * H).max() < 1e-12
    assert np.abs(M1 - H).max() > 1.0          # the paper's identity fails

    M2 = mn_operator(H, cs, 2)
    cdots = [time_derivative(H, C) for C in cs]
    sum_cdd = sum(Cd.conj().T @ Cd for Cd in cdots)
    assert np.abs(M2 + sum_cdd).max() < 1e-11
    assert np.abs(M2 - sum_cdd).max() > 1.0    # the paper's + sign fails


def test_mn_hermiticity_parity():
    """M_n^dag = (-1)^n M_n holds exactly for n = 0, 1, 2 -- and genuinely
    FAILS for n = 3 (the tower is not simply +-Hermitian beyond n = 2)."""
    Nc, p = 5, 4
    H = dirac_hamiltonian(Nc, p, rng=np.random.default_rng(11))
    cs = annihilators(Nc)
    for n in range(3):
        Mn = mn_operator(H, cs, n)
        scale = max(np.abs(Mn).max(), 1e-300)
        assert np.abs(Mn.conj().T - (-1) ** n * Mn).max() / scale < 1e-12
    M3 = mn_operator(H, cs, 3)
    assert np.abs(M3.conj().T - (-1) ** 3 * M3).max() / np.abs(M3).max() > 0.05


def test_cp_photon_graviton_exact():
    """On a C-invariant instance (disjoint/real ensemble), the unitary
    particle-hole map C gives exactly CP(Q - Nc/2) = -1 (photon) and
    CP(M_1) = +1 (graviton) -- the paper's n = 0, 1 assignments."""
    Nc, p = 6, 4
    rng = np.random.default_rng(4)
    Hc = dirac_hamiltonian(Nc, p,
                           couplings=draw_couplings(Nc, p, rng, "c_symmetric"))
    U = particle_hole(Nc)
    dim = 1 << Nc
    assert np.abs(U @ U.conj().T - np.eye(dim)).max() < 1e-12
    assert np.abs(cp_conjugate(U, Hc) - Hc).max() < 1e-12   # C-invariance
    Qbar = charge_matrix(Nc) - (Nc / 2) * np.eye(dim)
    assert np.abs(cp_conjugate(U, Qbar) + Qbar).max() < 1e-12     # CP = -1
    cs = annihilators(Nc)
    M1 = mn_operator(Hc, cs, 1)
    assert np.abs(cp_conjugate(U, M1) - M1).max() < 1e-11          # CP = +1


def test_cp_higher_n_is_a_mixture():
    """The measured finding this module documents: for n >= 2 the raw M_n is a
    CP-mixture even on a C-invariant instance (total-time-derivative
    contamination), so the paper's CP = (-1)^{n+1} for the dynamical tower is
    a statement about CP-resolved correlators, not about the raw operators."""
    Nc, p = 6, 4
    rng = np.random.default_rng(4)
    Hc = dirac_hamiltonian(Nc, p,
                           couplings=draw_couplings(Nc, p, rng, "c_symmetric"))
    U = particle_hole(Nc)
    cs = annihilators(Nc)
    dim = 1 << Nc
    M2 = mn_operator(Hc, cs, 2)
    M2 = M2 - np.trace(M2) / dim * np.eye(dim)
    w = wrong_channel_weight(U, M2, (-1) ** (2 + 1))
    assert 0.05 < w < 0.95, f"n=2 wrong-channel weight {w:.3f}"


# --------------------------------------------------------------------------
# Exact-Wick adjoint variance (exact_wick.py): Eq. 3.28 beyond ED
# --------------------------------------------------------------------------
@pytest.mark.parametrize("N,p", [(8, 4), (8, 6)])
def test_wick_closed_form_matches_enumeration(N, p):
    """Var(B_jk) = 4 sigma^4 C(N-2,p-1) must equal the literal sum over all
    subset pairs of Majorana-trace contractions (verifies every sign in the
    derivation, including T_cd T_dc = +1)."""
    from exact_wick import var_offdiag_exact, var_offdiag_enumerated
    ve, vn = var_offdiag_exact(N, p), var_offdiag_enumerated(N, p)
    assert abs(vn / ve - 1.0) < 1e-12


def test_wick_mean_diag_and_ratio_artifact():
    """E[B_jj] = sigma^2 C(N,p)(1 - 2p/N): matches enumeration, and vanishes
    exactly at p = N/2 -- the derived origin of the ratio artifact."""
    from exact_wick import mean_diag_exact, mean_diag_enumerated
    for (N, p) in [(8, 4), (8, 6), (10, 4)]:
        me, mn = mean_diag_exact(N, p), mean_diag_enumerated(N, p)
        assert abs(me - mn) < 1e-10 * max(1.0, abs(me))
    assert mean_diag_exact(12, 6) == 0.0
    assert mean_diag_exact(16, 8) == 0.0


def test_wick_matches_ED_ensemble():
    """The closed-form rms must reproduce a fresh ED disorder ensemble."""
    from exact_wick import rms_offdiag_exact
    from self_averaging import measure_symmetry_violation
    N, p = 10, 4
    r = measure_symmetry_violation(N, p, n_inst=100, seed=99)
    exact = rms_offdiag_exact(N, p)
    # E[RMS] <= sqrt(E[RMS^2]) = exact (Jensen), so allow a small one-sided
    # margin plus 3 standard errors
    assert r["rms_off"] < exact * 1.005 + 3 * r["rms_off_err"]
    assert r["rms_off"] > exact * 0.95 - 3 * r["rms_off_err"]


def test_wick_double_scaling_superexponential():
    """Along p = sqrt(N), ln rms / sqrt(N) must decrease without bound:
    the Eq. 3.28 form e^{-a sqrt N} is beaten for every a (checked well past
    any ED-reachable N)."""
    from exact_wick import double_scaling_scan
    rows = double_scaling_scan(1.0, [16, 64, 100, 400, 1600])
    slopes = [r["ln_rms"] / np.sqrt(r["N"]) for r in rows]
    assert all(b < a for a, b in zip(slopes, slopes[1:]))
    assert slopes[-1] < -2.0


def test_relvar_m2_matches_ED():
    """ED relative variance of m2 matches the analytic 2/C(N,p) to ~few %."""
    N, p = 12, 4
    res = measure_moments(N, p, n_inst=400, seed=7, ks=(1,))
    relvar = res[1][0]
    analytic = relvar_m2_exact(N, p)
    assert abs(relvar - analytic) / analytic < 0.15


def test_self_averaging_decreases_with_N():
    """Singlet variance must fall as N grows (the whole point of self-averaging)."""
    rv = [relvar_m2_exact(N, 4) for N in (10, 14, 18)]
    assert rv[0] > rv[1] > rv[2]


@pytest.mark.slow
def test_bott_periodicity_two_classes():
    """A cheap slice of the RMT validation: GUE (N=10) and GSE (N=12)."""
    r10 = validate_syk.r_statistic(10, n_real=40, seed=110)
    r12 = validate_syk.r_statistic(12, n_real=40, seed=112)
    assert validate_syk._nearest_class(r10) == "GUE"
    assert validate_syk._nearest_class(r12) == "GSE"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
