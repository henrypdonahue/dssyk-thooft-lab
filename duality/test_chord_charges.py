"""The conserved-channel theorem (chord_charges.py), asserted end to end.

Pinned: the three-term close_j recursion as a matrix identity; the dressed
propagator C(m) = P_m(T) (Motzkin polynomial) and [C(m), T] = 0 for
m <= 8; the C(1), C(2) and Xhat_2 closed forms; [Xhat_n, T] = 0 for
n <= 5; and, through the INDEPENDENT chord_moments word assembly, the
theorem's punchline: every |mu_k|/gross <= 1e-12 for n <= 5, k = 1..6,
both operator orderings.  All identities are exact; residuals are float
noise, measured on the truncation-free interior block relative to scale.
"""

from math import comb

import numpy as np
import pytest

from chord import qint
from chord_charges import (close_j_map, dressed_propagator, motzkin_coeffs,
                           open_map, poly_of_matrix, transfer_formal,
                           transfer_one_matter_formal, xhat_matrix)

CASES = [(0.5, 0.25), (0.85, 0.25), (0.2, 0.7), (0.5, 1.0)]
NMAX = 24


def _interior(A, keep):
    return A[:keep, :keep]


@pytest.mark.parametrize("q,dpsi", CASES)
def test_close_recursion_identity(q, dpsi):
    """close_j T_M = q^(D+j) T close_j + [j] close_{j-1}
    + (1 - q^(2D+j)) close_{j+1}, j <= 5, as matrices (interior block)."""
    T = transfer_formal(NMAX, q)
    TM = transfer_one_matter_formal(NMAX, q, dpsi)
    closes = [close_j_map(NMAX, q, dpsi, j) for j in range(7)]
    m = NMAX + 1
    keep = NMAX - 1
    cols = [n0 * m + n1 for n0 in range(m) for n1 in range(m)
            if n0 + n1 <= keep - 1]
    for j in range(6):
        lhs = closes[j] @ TM
        rhs = q ** (dpsi + j) * (T @ closes[j]) \
            + (1.0 - q ** (2 * dpsi + j)) * closes[j + 1]
        if j >= 1:
            rhs = rhs + qint(j, q) * closes[j - 1]
        scale = max(np.abs(lhs[:keep][:, cols]).max(), 1.0)
        dev = np.abs((lhs - rhs)[:keep][:, cols]).max()
        assert dev < 1e-12 * scale, (j, dev, scale)


@pytest.mark.parametrize("q,dpsi", CASES)
def test_dressed_propagator_is_polynomial_in_T(q, dpsi):
    """C(m) = P_m(T) and [C(m), T] = 0 for m <= 8 (interior, relative)."""
    T = transfer_formal(NMAX, q)
    keep = NMAX - 8
    for m in range(0, 9):
        Cm = dressed_propagator(m, NMAX, q, dpsi)
        Pm = poly_of_matrix(motzkin_coeffs(m, q, dpsi), T)
        scale = max(np.abs(_interior(Cm, keep)).max(), 1.0)
        assert np.abs(_interior(Cm - Pm, keep)).max() < 1e-11 * scale, m
        comm = Cm @ T - T @ Cm
        assert np.abs(_interior(comm, keep)).max() < 1e-11 * scale, m


@pytest.mark.parametrize("q,dpsi", CASES)
def test_closed_forms(q, dpsi):
    """C(1) = q^D T; C(2) = q^(2D) T^2 + (1-q^(2D));
    Xhat_2 = (1-q^D)^2 T^2 + (1-q^(2D))."""
    T = transfer_formal(NMAX, q)
    keep = NMAX - 4
    one = np.eye(NMAX + 1)
    C1 = dressed_propagator(1, NMAX, q, dpsi)
    assert np.abs(_interior(C1 - q ** dpsi * T, keep)).max() < 1e-12
    C2 = dressed_propagator(2, NMAX, q, dpsi)
    ref2 = q ** (2 * dpsi) * (T @ T) + (1 - q ** (2 * dpsi)) * one
    assert np.abs(_interior(C2 - ref2, keep)).max() < 1e-12
    X2 = xhat_matrix(2, NMAX, q, dpsi)
    refx = (1 - q ** dpsi) ** 2 * (T @ T) + (1 - q ** (2 * dpsi)) * one
    assert np.abs(_interior(X2 - refx, keep)).max() < 1e-12
    # X_1 is proportional to T: the N = infinity image of A_1 = -2pH
    X1 = xhat_matrix(1, NMAX, q, dpsi)
    assert np.abs(_interior(X1 + (1 - q ** dpsi) * T, keep)).max() < 1e-12


@pytest.mark.parametrize("q,dpsi", [(0.5, 0.25), (0.2, 0.7)])
def test_xhat_commutes_with_T(q, dpsi):
    """[Xhat_n, T] = 0 for n <= 5: the conserved charges exist at every n."""
    T = transfer_formal(NMAX, q)
    keep = NMAX - 8
    for n in range(1, 6):
        X = xhat_matrix(n, NMAX, q, dpsi)
        comm = X @ T - T @ X
        scale = max(np.abs(_interior(X, keep)).max(), 1.0)
        assert np.abs(_interior(comm, keep)).max() < 1e-11 * scale, n


@pytest.mark.parametrize("q,dpsi", [(0.5, 0.25), (0.85, 0.25), (0.2, 0.7)])
def test_all_moments_vanish_word_level(q, dpsi):
    """The theorem through the INDEPENDENT engine route: the chord_moments
    word assembly gives |mu_k|/gross <= 1e-12 for n <= 5, k = 1..6, each
    operator ordering separately (odd k included -- stronger than the
    original mu_2/mu_4 discovery)."""
    from chord_moments import MomentTable, _mu_ordered, _mu_ordered_swapped
    tab = MomentTable(q, dpsi)
    for n in range(1, 6):
        for k in range(1, 7):
            for fn in (_mu_ordered, _mu_ordered_swapped):
                v, g = fn(tab, n, k)
                assert abs(v) <= 1e-12 * max(g, 1e-300), (n, k, fn.__name__)


def test_motzkin_low_orders_by_hand():
    """P_m against hand-enumerated Motzkin paths for m <= 3 at symbolic-level
    spot values (path weights: flat q^(D+j) T, down [j], up 1-q^(2D+j))."""
    q, d = 0.37, 0.6
    # m = 3: flat^3 + (up down flat@0) + (up flat@1 down) + (flat@0 up down)
    c3 = motzkin_coeffs(3, q, d)
    flat0, flat1 = q ** d, q ** (d + 1)
    up0, down1 = 1 - q ** (2 * d), qint(1, q)
    expect_T1 = up0 * down1 * flat0 * 2 + up0 * flat1 * down1
    assert c3[3] == pytest.approx(flat0 ** 3, rel=1e-14)
    assert c3[1] == pytest.approx(expect_T1, rel=1e-14)
    assert c3[0] == 0.0 and c3[2] == 0.0
