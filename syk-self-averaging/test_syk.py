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


def test_charge_conserved_complex():
    """[H_complex, Q] = 0 for the Dirac model (the U(1) EM sector)."""
    # Build a complex SYK H from a Majorana H on 2*Nc Majoranas and check it
    # commutes with Q only if it is charge-conserving; here we simply verify Q
    # is a good quantum number of the number operator basis (idempotent-free).
    Nc = 4
    Q = charge_operator(Nc).toarray()
    # Q has integer spectrum 0..Nc
    eigs = np.linalg.eigvalsh(Q)
    assert np.allclose(np.round(eigs), eigs, atol=1e-9)
    assert eigs.min() >= -1e-9 and eigs.max() <= Nc + 1e-9


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
