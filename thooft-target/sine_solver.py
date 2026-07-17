#!/usr/bin/env python3
r"""
Independent cross-check solver for 't Hooft's equation (sine basis)
===================================================================

This is a *deliberately independent* second solver, used only to validate the
main Chebyshev-U spectral solver (thooft_spectrum.py) -- in particular the
quark-mass (potential) machinery, which is invisible at the duality point
alpha = 0.  It shares no code, no basis, and no closed-form matrix elements with
the main solver.

Method
------
Expand the wavefunction in the sine basis

    phi(x) = sum_k a_k sin(k pi x),    k = 1..N,

which satisfies the boundary conditions phi(0) = phi(1) = 0 by construction.
Unlike the Chebyshev-U basis it does NOT build in the sqrt(x(1-x)) endpoint
behaviour, so it converges only algebraically -- fine for a 3-4 digit
cross-check, useless for high precision.  That mismatch is the point: agreement
between the two methods cannot be an artefact of a shared discretization.

The finite-part (hypersingular) operator is handled through its exact weak form.
For any u, v vanishing at the endpoints,

    Int_0^1 u(x) [ FP-Int_0^1 v(y)/(y-x)^2 dy ] dx
        = - (1/2) Int_0^1 Int_0^1 [u(x)-u(y)][v(x)-v(y)] / (x-y)^2  dx dy .

The right-hand integrand is *regular* on the diagonal (it tends to u'(x)v'(x) as
y -> x), so it is evaluated by an ordinary tensor Gauss-Legendre rule.  't
Hooft's equation in FLZ normalization,

    2 pi^2 lambda phi = alpha (1/x + 1/(1-x)) phi - FP-Int phi(y)/(y-x)^2 dy,

then becomes the generalized eigenproblem

    (alpha * Pi - W) a = pi^2 lambda' G a,     2 lambda = mu / pi^2,

with G_{nm} = <sin_n, sin_m>, Pi_{nm} = <sin_n, (1/x+1/(1-x)) sin_m>, and
W_{nm} the weak-form finite-part matrix above.  (mu = 2 pi^2 lambda.)
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

PI = np.pi


def _gauss_legendre_01(n_quad: int):
    """Gauss-Legendre nodes/weights mapped from [-1,1] to [0,1]."""
    x, w = np.polynomial.legendre.leggauss(n_quad)
    return 0.5 * (x + 1.0), 0.5 * w


def build_matrices(num_basis: int = 60, n_quad: int = 400, alpha: float = 0.0):
    """Assemble (G, Pi, W) in the sine basis by Gauss-Legendre quadrature."""
    ks = np.arange(1, num_basis + 1)
    xq, wq = _gauss_legendre_01(n_quad)

    # basis and derivative values at the quadrature nodes:  Psi[k, i] = sin(k pi x_i)
    Psi = np.sin(np.outer(ks, PI * xq))          # (K, Nq)
    dPsi = (ks[:, None] * PI) * np.cos(np.outer(ks, PI * xq))  # d/dx sin(k pi x)

    # Gram matrix  G_{nm} = Int sin_n sin_m  (= 1/2 delta, but computed by quad)
    G = (Psi * wq) @ Psi.T

    # potential matrix  Pi_{nm} = Int (1/x + 1/(1-x)) sin_n sin_m
    v = 1.0 / xq + 1.0 / (1.0 - xq)
    Pi = (Psi * (wq * v)) @ Psi.T

    # weak-form finite-part matrix
    #   W = -1/2 * Int Int [sin_n(x)-sin_n(y)][sin_m(x)-sin_m(y)]/(x-y)^2
    # Off-diagonal (x_i != x_j) collapses to  Psi (diag(r) - Q) Psi^T  with
    #   Q_{ij} = w_i w_j / (x_i-x_j)^2,   r_i = sum_j Q_{ij},
    # and the removable diagonal (limit u'v') contributes  Psi' diag(w^2) Psi'^T.
    dx = xq[:, None] - xq[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        Q = np.outer(wq, wq) / (dx * dx)
    np.fill_diagonal(Q, 0.0)
    r = Q.sum(axis=1)
    L = np.diag(r) - Q                              # graph-Laplacian-like
    off = Psi @ L @ Psi.T
    diag = (dPsi * (wq ** 2)) @ dPsi.T
    # W_gag_{nm} = -1/2 Int Int [Δsin_n][Δsin_m]/(x-y)^2  (the Gagliardo form).
    # It equals <sin_n, FP sin_m> + <sin_n, (1/x+1/(1-x)) sin_m>, i.e. it already
    # bundles in the endpoint self-energy potential -- see solve().
    W = -(off + 0.5 * diag)
    return G, Pi, W


def solve(num_basis: int = 60, n_quad: int = 400, alpha: float = 0.0,
          num_levels: int = 8):
    """Lowest `num_levels` values of 2 lambda_n = M_n^2/(pi g^2), ascending.

    This mixes both CP sectors (no parity splitting); return the full ordered
    list and let the caller compare against the interleaved main spectrum.
    """
    G, Pi, W = build_matrices(num_basis, n_quad, alpha)
    # <u,FP v> = W_gag - Pi, so  alpha*Pi - <.,FP.> = (alpha+1)*Pi - W_gag.
    # The coefficient (alpha+1) = pi m^2/g^2 is the *bare* mass: the endpoint
    # self-energy from the finite part shifts alpha -> alpha+1, as it must.
    A = (alpha + 1.0) * Pi - W
    mu = eigh(A, G, eigvals_only=True)             # mu = 2 pi^2 lambda
    two_lambda = np.sort(mu) / PI ** 2
    return two_lambda[:num_levels]


if __name__ == "__main__":
    from thooft_spectrum import spectrum

    for alpha in (0.0, 0.5, 2.0):
        cheb = [L["two_lambda"] for L in spectrum(num_basis=300, num_levels=8,
                                                  alpha=alpha)]
        sine = solve(num_basis=80, n_quad=600, alpha=alpha, num_levels=8)
        print(f"\nalpha = {alpha}")
        print(f"{'n':>3} {'sine (indep.)':>16} {'chebyshev':>16} {'|diff|':>10}")
        print("-" * 50)
        for n in range(8):
            d = abs(sine[n] - cheb[n])
            print(f"{n:>3} {sine[n]:>16.9f} {cheb[n]:>16.9f} {d:>10.1e}")
