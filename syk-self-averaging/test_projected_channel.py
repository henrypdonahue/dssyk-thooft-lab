"""Independent checks of the finite-N projection, including degeneracies."""
import numpy as np
import pytest
from projected_channel import spectral_split


def test_degenerate_block_is_stationary_not_just_diagonal():
    H = np.diag([0., 0., 2.])
    A = np.array([[0., 3., 1.], [3., 0., 0.], [1., 0., 0.]])
    r = spectral_split(H, A)
    assert r['mu0_stationary'] == pytest.approx(6.)
    assert r['mu0_dynamic'] == pytest.approx(2/3)
    assert r['w2_dynamic'] == pytest.approx(4.)
    theta = .713
    U = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1.]])
    rotated = spectral_split(U @ H @ U.T, U @ A @ U.T)
    for k in ('mu0_stationary', 'mu0_dynamic', 'w2_dynamic'):
        assert rotated[k] == pytest.approx(r[k], abs=1e-13)


def test_polynomial_is_removed_and_identity_shift_does_not_change_result():
    H = np.diag([-2., 1., 3.])
    A = np.ones((3, 3))
    r = spectral_split(H, A)
    shifted = spectral_split(H, A + H @ H + 7*np.eye(3))
    assert shifted['mu0_dynamic'] == pytest.approx(r['mu0_dynamic'])
    assert shifted['w2_dynamic'] == pytest.approx(r['w2_dynamic'])
    conserved = spectral_split(H, H @ H)
    assert not conserved['dynamic_resolved']
    assert conserved['w2_dynamic'] is None


def test_a1_conserved_a2_and_a3_have_resolved_weight():
    from conserved_points import a_n_tower, dense_majoranas
    from pauli_strings import majorana_hamiltonian_fast
    H = majorana_hamiltonian_fast(10, 4, np.random.default_rng(17))
    tower = a_n_tower(H, dense_majoranas(10), 3)
    assert not spectral_split(H, tower[1])['dynamic_resolved']
    for n in (2, 3):
        r = spectral_split(H, tower[n])
        assert r['dynamic_resolved']
        assert 0 < r['mu0_dynamic'] < r['mu0_traceless']
        assert r['w2_dynamic'] > 0


def test_invalid_hamiltonian_rejected():
    with pytest.raises(ValueError, match='Hermitian'):
        spectral_split(np.array([[0, 1], [0, 0]]), np.eye(2))
