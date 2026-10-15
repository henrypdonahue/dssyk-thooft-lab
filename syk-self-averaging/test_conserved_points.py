"""Tests for conserved_points.py: the (8, 4) per-instance conservation
explained by the spin-trace map and the grade content of H^k, plus the
two exact identities A_1 = -2 sum H_i and A_2 = 4 sum H_i^2."""

import numpy as np
import pytest

from conserved_points import (a_n_tower, binomial_assembly, chirality,
                              commutant_residual, commutator_residual,
                              dense_majoranas, flavor_parts, general_identity,
                              grade_content, identity_coefficients,
                              spin_trace, spin_trace_weight)
from pauli_strings import majorana_hamiltonian_fast


def _instance(N, p, seed=0):
    H = majorana_hamiltonian_fast(N, p, np.random.default_rng(seed))
    return H, dense_majoranas(N)


@pytest.mark.parametrize("N", [6, 8, 10])
def test_spin_trace_is_diagonal_on_grades(N):
    """Eq. (1): S(psi_J) = (-1)^g (N - 2g) psi_J on every monomial grade."""
    psis = dense_majoranas(N)
    dim = psis[0].shape[0]
    for g in range(N + 1):
        M = np.eye(dim, dtype=complex)
        for j in range(g):
            M = M @ psis[j]
        assert np.allclose(spin_trace(M, psis), spin_trace_weight(N, g) * M,
                           atol=1e-12)


@pytest.mark.parametrize("N,p", [(8, 4), (10, 4), (12, 6), (8, 2)])
def test_binomial_assembly_identity(N, p):
    """Eq. (2): A_n = sum_j (-1)^j C(n,j) S(H^{n-j}) H^j, any (N, p), n <= 4."""
    H, psis = _instance(N, p)
    A = a_n_tower(H, psis, 4)
    for n in range(1, 5):
        assert (np.linalg.norm(binomial_assembly(H, psis, n) - A[n])
                < 1e-12 * np.linalg.norm(A[n]))


@pytest.mark.parametrize("N,p", [(8, 4), (10, 4), (12, 6), (12, 4), (8, 2)])
def test_flavor_part_identities(N, p):
    """Eq. (3): sum_i H_i = pH, A_1 = -2 sum_i H_i, A_2 = 4 sum_i H_i^2 at
    every (N, p); and the pattern stops at n = 3."""
    H, psis = _instance(N, p, seed=2)
    Hi = flavor_parts(H, psis)
    A = a_n_tower(H, psis, 3)
    assert np.linalg.norm(sum(Hi) - p * H) < 1e-12 * np.linalg.norm(H)
    assert np.linalg.norm(A[1] + 2 * sum(Hi)) < 1e-12 * np.linalg.norm(A[1])
    S2 = 4 * sum(h @ h for h in Hi)
    assert np.linalg.norm(A[2] - S2) < 1e-12 * np.linalg.norm(A[2])
    # n = 3: not -8 sum H_i^3 (nor any fixed multiple of it)
    S3 = sum(h @ h @ h for h in Hi)
    c = np.vdot(S3.ravel(), A[3].ravel()) / np.vdot(S3.ravel(), S3.ravel())
    assert np.linalg.norm(A[3] - c * S3) > 0.05 * np.linalg.norm(A[3])


@pytest.mark.parametrize("N,p", [(8, 4), (12, 6)])
def test_spin_trace_kills_H_at_half_filling(N, p):
    """p = N/2: S(H) = 0 per instance (grade N/2 has weight zero)."""
    H, psis = _instance(N, p)
    assert np.linalg.norm(spin_trace(H, psis)) < 1e-13 * np.linalg.norm(H)


def test_n8_powers_of_H_live_in_grades_0_4_8():
    """At N = 8 every H^k (p = 4) has no grade-2 or grade-6 content."""
    H, psis = _instance(8, 4, seed=3)
    Hk = H
    for _ in range(3):
        Hk = Hk @ H
        gc = grade_content(Hk, psis)
        assert gc[2] < 1e-12 * gc[4] and gc[6] < 1e-12 * gc[4]
        assert gc[4] > 0.1 * gc[0]


def test_n8p4_tower_is_in_the_commutant():
    """(8, 4): A_1..A_5 lie in span{H^j, Gamma H^j} and commute with H."""
    for seed in (0, 5):
        H, psis = _instance(8, 4, seed)
        G = chirality(psis)
        A = a_n_tower(H, psis, 5)
        for n in range(1, 6):
            assert commutator_residual(A[n], H) < 1e-13
            assert commutant_residual(A[n], H, G) < 1e-10


def test_n8p4_A2_closed_form():
    """Eq. (4): A_2 = 8 H^2 + 8 (H^2)_0 - 8 (H^2)_8 at (8, 4) -- the fitted
    H^2 coefficient is exactly 8 and the remainder is central."""
    H, psis = _instance(8, 4, seed=1)
    dim = H.shape[0]
    G = chirality(psis)
    A2 = a_n_tower(H, psis, 2)[2]
    B = np.stack([np.eye(dim).ravel(), (H @ H).ravel(), G.ravel()], 1)
    c = np.linalg.lstsq(B, A2.ravel(), rcond=None)[0]
    assert abs(c[1] - 8.0) < 1e-10
    assert np.linalg.norm(B @ c - A2.ravel()) < 1e-12 * np.linalg.norm(A2)
    # the central coefficients are the stated grade projections
    a = np.trace(H @ H) / dim
    b = np.trace(G @ H @ H) / dim / (np.trace(G @ G) / dim)
    assert abs(c[0] - 8 * a) < 1e-10 * abs(a)
    assert abs(c[2] + 8 * b) < 1e-10 * max(abs(b), 1e-3)


@pytest.mark.parametrize("N,p", [(10, 4), (12, 6), (12, 4)])
def test_failure_away_from_n8p4(N, p):
    """The mechanism is N = 8 specific: elsewhere [H, A_2] is O(1)."""
    H, psis = _instance(N, p)
    A = a_n_tower(H, psis, 3)
    assert commutator_residual(A[2], H) > 0.1
    assert commutator_residual(A[3], H) > 0.1


@pytest.mark.parametrize("N,p", [(6, 2), (8, 2), (10, 2), (8, 6)])
def test_free_and_dual_corners_conserved(N, p):
    """p = 2 (free) at any N, and (8, 6) = Gamma x (8, 2): all A_n conserved."""
    H, psis = _instance(N, p)
    A = a_n_tower(H, psis, 4)
    for n in range(1, 5):
        assert commutator_residual(A[n], H) < 1e-13


@pytest.mark.parametrize("N,p", [(8, 2), (8, 4), (10, 4), (10, 6)])
def test_general_identity_every_n(N, p):
    """Eq. (3'): A_n = sum_{|T|>=2} (-2)^|T| (-1)^(n-max T) C(n-1, max T-1)
    W_T per instance, n = 2..5, any (N, p) -- the mixed-word continuation
    of A_2 = 4 sum H_i^2."""
    H, psis = _instance(N, p, seed=4)
    A = a_n_tower(H, psis, 5)
    for n in range(2, 6):
        assert (np.linalg.norm(general_identity(H, psis, n) - A[n])
                < 1e-12 * np.linalg.norm(A[n])), n


def test_general_identity_coefficients():
    """The coefficient table: n = 2 is {(1,2): 4}; n = 3 is
    4 (W_13 + W_23 - 2 W_12) - 8 W_123; the pair coefficients depend on
    the later slot only and sum to zero over it for n >= 3."""
    assert identity_coefficients(2) == {(1, 2): 4}
    assert identity_coefficients(3) == {(1, 2): -8, (1, 3): 4, (2, 3): 4,
                                        (1, 2, 3): -8}
    for n in (4, 5, 6):
        c = identity_coefficients(n)
        for b in range(2, n + 1):
            assert len({c[(a, b)] for a in range(1, b)}) == 1
        assert c[(1, n)] == 4
        assert sum(c[(1, b)] for b in range(2, n + 1)) == 4 * (-1) ** n


@pytest.mark.parametrize("N,p", [(10, 4), (12, 6), (12, 4)])
def test_antihermitian_parts_are_total_derivatives(N, p):
    """Two corollaries of Eq. (3'), per instance: A_2 is Hermitian,
    A_3 - A_3^dag = 3 [H, A_2], and A_4 - A_4^dag = 2 [H, A_3 + A_3^dag].
    The first is why the O_2 kernel of A_3 at N = infinity carries the
    commutator with coefficient exactly 3/2 (duality/singlet_tower.py)."""
    H, psis = _instance(N, p, seed=6)
    A = a_n_tower(H, psis, 4)

    def comm(X, Y):
        return X @ Y - Y @ X

    dag = lambda X: X.conj().T
    assert np.linalg.norm(A[2] - dag(A[2])) < 1e-12 * np.linalg.norm(A[2])
    assert (np.linalg.norm(A[3] - dag(A[3]) - 3 * comm(H, A[2]))
            < 1e-12 * np.linalg.norm(A[3]))
    assert (np.linalg.norm(A[4] - dag(A[4]) - 2 * comm(H, A[3] + dag(A[3])))
            < 1e-12 * np.linalg.norm(A[4]))
