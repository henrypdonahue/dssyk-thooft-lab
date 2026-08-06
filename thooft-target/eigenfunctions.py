#!/usr/bin/env python3
r"""
Eigenfunction-level validation of the 't Hooft solver.
======================================================

The eigenVALUES of thooft_spectrum.py are certified against FLZ; until now the
eigenVECTORS were computed but never checked.  The eventual DSSYK comparison
needs matrix elements (decay constants, form factors), so the wavefunctions
must be validated too.  Three independent checks, all encoded in
test_thooft.py:

1.  Parseval / completeness of the "decay constants".  With phi_n unit-
    normalized in L^2[0,1], completeness of the eigenbasis gives exact sum
    rules for the overlaps with 1 and with x:

        sum_n ( Int_0^1 phi_n dx )^2   = ||1||^2 = 1        (even sector only)
        sum_n ( Int_0^1 x phi_n dx )^2 = ||x||^2 = 1/3      (both sectors)

    In the Chebyshev basis the overlaps are closed-form:
    Int_0^1 sqrt(1-xi^2) U_m dx = (pi/4) delta_{m0} and
    Int_0^1 xi sqrt(1-xi^2) U_m dx = (pi/8) delta_{m1}, so the first overlap
    is (pi/4) a_0^(n) -- the 't Hooft-model decay-constant column of the
    eigenvector matrix.  These are genuine tests of the eigenvectors as a
    complete orthonormal set.

2.  Independent Rayleigh-quotient recomputation.  Each reconstructed phi_n(x)
    is pushed through the *weak-form* (Gagliardo) quadratic form of the
    finite-part operator on Gauss-Legendre nodes -- the sine_solver's
    discretization-free representation, sharing no basis or closed-form matrix
    elements with the Chebyshev solver:

        2 lambda[phi] = [ (alpha+1) Int V phi^2 + Gag[phi, phi] ] / (pi^2 Int phi^2)

    (V = 1/x + 1/(1-x); the alpha -> alpha+1 shift is the endpoint self-energy
    bundled into the Gagliardo form, exactly as in sine_solver.solve).
    Recovering 2 lambda_n from the eigenvector to quadrature accuracy
    validates the wavefunction against the operator itself.

3.  Cross-solver wavefunction overlap: |<phi_n^cheb, phi_n^sine>| -> 1.  The
    sine-basis eigenvectors come from a different basis and different
    singularity handling; overlap deviating from 1 only at the sine solver's
    documented algebraic accuracy validates the shape, not just the value.
"""

from __future__ import annotations

import numpy as np

from thooft_spectrum import solve_sector

PI = np.pi


# ---------------------------------------------------------------------------
# reconstruction
# ---------------------------------------------------------------------------
def _U_all(max_degree: int, xi: np.ndarray) -> np.ndarray:
    """U_m(xi) for m = 0..max_degree, rows m, via the stable recurrence."""
    U = np.zeros((max_degree + 1, len(xi)))
    U[0] = 1.0
    if max_degree >= 1:
        U[1] = 2.0 * xi
    for m in range(2, max_degree + 1):
        U[m] = 2.0 * xi * U[m - 1] - U[m - 2]
    return U


def reconstruct(vec: np.ndarray, degrees: np.ndarray, x: np.ndarray) -> np.ndarray:
    """phi_n(x) on x in (0,1), UNIT-normalized in L^2[0,1].

    eigh(D, S) returns S-orthonormal coefficient vectors, i.e.
    Int_{-1}^{1} phi^2 dxi = 1, and with x = (1+xi)/2 that means
    Int_0^1 phi^2 dx = 1/2 -- hence the sqrt(2)."""
    xi = 2.0 * np.asarray(x) - 1.0
    U = _U_all(int(np.max(degrees)), xi)
    phi = np.sqrt(1.0 - xi ** 2) * (vec @ U[np.asarray(degrees)])
    return np.sqrt(2.0) * phi


def reconstruct_derivative(vec: np.ndarray, degrees: np.ndarray,
                           x: np.ndarray) -> np.ndarray:
    """d phi_n / dx at interior points (chain rule d xi/dx = 2), same
    normalization as reconstruct.

    Precondition: x must be interior (Gauss-Legendre nodes) -- this divides by
    sqrt(1-xi^2) and by (1-xi^2), so it blows up at the endpoints x=0,1."""
    xi = 2.0 * np.asarray(x) - 1.0
    dmax = int(np.max(degrees))
    U = _U_all(dmax + 1, xi)
    s = np.sqrt(1.0 - xi ** 2)
    total = np.zeros_like(xi)
    for a, m in zip(vec, np.asarray(degrees)):
        # d/dxi [ sqrt(1-xi^2) U_m ] = ( -xi U_m + (1-xi^2) U_m' ) / sqrt(1-xi^2)
        # with U_m'(xi) = ( (m+2) U_{m-1} - m U_{m+1} ) / ( 2 (1-xi^2) )  for m>=1
        if m == 0:
            dU = np.zeros_like(xi)
        else:
            dU = ((m + 2) * U[m - 1] - m * U[m + 1]) / (2.0 * (1.0 - xi ** 2))
        total += a * (-xi * U[m] + (1.0 - xi ** 2) * dU) / s
    return np.sqrt(2.0) * 2.0 * total


def canonicalize_signs(vecs: np.ndarray, degrees) -> np.ndarray:
    """Fix the physically-meaningless per-eigenvector sign left arbitrary by
    LAPACK: flip each column so that phi_n(x) > 0 as x -> 0+.

    Near x = 0 (xi = -1), phi ~ sqrt(1-xi^2) * sum_m a_m U_m(-1) with
    U_m(-1) = (-1)^m (m+1), so the sign at the endpoint is the sign of
    sum_m a_m (-1)^m (m+1).  Without this, quantities linear in the
    wavefunction (decay constants, form factors) carry an arbitrary
    build-dependent sign."""
    d = np.asarray(degrees)
    endpoint = ((-1.0) ** d * (d + 1.0)) @ vecs
    flip = np.where(endpoint < 0, -1.0, 1.0)
    return vecs * flip[None, :]


# ---------------------------------------------------------------------------
# 1. decay constants and Parseval sums
# ---------------------------------------------------------------------------
def decay_constants(num_basis: int = 200, num_levels: int = 60):
    """f_n = Int_0^1 phi_n dx for the symmetric sector (antisymmetric levels
    have f = 0 by parity), unit-normalized phi, in the phi_n(0+) > 0 sign
    convention (see canonicalize_signs -- without a stated convention the
    signs would be arbitrary LAPACK output).  Closed form:
    f = sqrt(2) (pi/4) a_0^(n)."""
    ev, vecs, degrees = solve_sector(num_basis, parity=0)
    vecs = canonicalize_signs(vecs, degrees)
    return np.sqrt(2.0) * (PI / 4.0) * vecs[0, :num_levels], ev[:num_levels]


def parseval_unit(num_basis: int = 200) -> float:
    """sum_n f_n^2 over the whole symmetric tower; must -> 1 (completeness)."""
    f, _ = decay_constants(num_basis, num_levels=num_basis)
    return float(np.sum(f ** 2))


def parseval_x(num_basis: int = 200) -> float:
    """sum_n ( Int x phi_n )^2 over both sectors; must -> 1/3.

    Int_0^1 x phi dx = (1/2)Int (1+xi)/2 sqrt(1-xi^2) sum a_m U_m dxi
                     = sqrt(2) [ (pi/8) a_0 + (pi/16) a_1 ]  in raw coefficients
    (a_0 from the even sector, a_1 from the odd sector)."""
    total = 0.0
    for parity in (0, 1):
        ev, vecs, degrees = solve_sector(num_basis, parity=parity)
        if parity == 0:
            overlaps = np.sqrt(2.0) * (PI / 8.0) * vecs[0, :]
        else:
            overlaps = np.sqrt(2.0) * (PI / 16.0) * vecs[0, :]
        total += float(np.sum(overlaps ** 2))
    return total


# ---------------------------------------------------------------------------
# 2. independent weak-form Rayleigh quotient
# ---------------------------------------------------------------------------
def rayleigh_weakform(vec, degrees, alpha: float = 0.0,
                      n_quad: int = 600) -> float:
    """Recompute 2*lambda from a Chebyshev eigenvector through the Gagliardo
    weak form on Gauss-Legendre nodes -- no Chebyshev closed forms involved."""
    xq, wq = np.polynomial.legendre.leggauss(n_quad)
    xq = 0.5 * (xq + 1.0)
    wq = 0.5 * wq
    phi = reconstruct(vec, degrees, xq)
    dphi = reconstruct_derivative(vec, degrees, xq)

    norm = np.sum(wq * phi ** 2)
    V = 1.0 / xq + 1.0 / (1.0 - xq)
    pot = np.sum(wq * V * phi ** 2)

    dx = xq[:, None] - xq[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        Q = np.outer(wq, wq) / (dx * dx)
    np.fill_diagonal(Q, 0.0)
    r = Q.sum(axis=1)
    off = np.sum(r * phi ** 2) - phi @ Q @ phi
    diag = np.sum((wq * dphi) ** 2)
    gag = off + 0.5 * diag

    mu = (alpha + 1.0) * pot + gag
    return float(mu / (PI ** 2 * norm))


if __name__ == "__main__":
    print("Eigenfunction-level validation (alpha = 0)\n")

    f, ev = decay_constants(200, num_levels=8)
    print("Decay constants f_n = Int phi_n dx (symmetric sector):")
    for n in range(8):
        print(f"  n_FLZ={2*n:>2}  2lam={ev[n]:>12.8f}   f = {f[n]:+.8f}")

    p1 = parseval_unit(200)
    px = parseval_x(200)
    print(f"\nParseval: sum f_n^2       = {p1:.12f}   (exact 1)")
    print(f"          sum (x-overlap)^2 = {px:.12f}   (exact 1/3 = {1/3:.12f})")

    print("\nWeak-form Rayleigh quotient vs eigenvalue (independent quadrature):")
    ev0, vecs0, deg0 = solve_sector(200, parity=0)
    ev1, vecs1, deg1 = solve_sector(200, parity=1)
    for n in range(4):
        rq = rayleigh_weakform(vecs0[:, n], deg0)
        print(f"  sym  n={n}: 2lam = {ev0[n]:.10f}   R[phi] = {rq:.10f}   "
              f"diff = {rq-ev0[n]:+.2e}")
    for n in range(2):
        rq = rayleigh_weakform(vecs1[:, n], deg1)
        print(f"  anti n={n}: 2lam = {ev1[n]:.10f}   R[phi] = {rq:.10f}   "
              f"diff = {rq-ev1[n]:+.2e}")
