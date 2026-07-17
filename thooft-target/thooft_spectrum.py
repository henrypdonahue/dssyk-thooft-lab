#!/usr/bin/env python3
r"""
't Hooft meson spectrum at the renormalized-massless point
===========================================================

Solves 't Hooft's (1974) integral equation for the meson mass spectrum of
large-N QCD in 1+1 dimensions, at the special quark-mass point that the
DSSYK_infty <-> 't Hooft duality needs (Miyashita-Sekino-Susskind, 2607.05678).

The equation (Fateev-Lukyanov-Zamolodchikov, arXiv:0905.2280, Eq. 1.1) is

    2 pi^2 lambda  phi(x)
        = ( alpha1/x + alpha2/(1-x) ) phi(x)
          - FP-Int_0^1  phi(y) / (y - x)^2  dy ,        x in [0,1]

with the Hadamard finite-part (double-pole) integral, boundary conditions
phi(0) = phi(1) = 0, quark-mass parameters

    alpha_i = pi m_i^2 / g^2 - 1 ,

and physical meson masses  M_n^2 = 2 pi g^2 lambda_n.

The point we want:  m_q^2 (renormalized) = 0.
    In Susskind et al. the renormalized mass obeys  m_q^2 = m_bare^2 - g^2/pi,
    so alpha = pi m_q^2 / g^2 = 0.  This is FLZ's exactly-studied case
    alpha1 = alpha2 = 0  (bare mass m = g/sqrt(pi)).  It is NOT the chiral
    limit (that is alpha = -1); it is a regular theory with a mass gap.

Units.  We report the dimensionless eigenvalue
    2 lambda_n = M_n^2 / (pi g^2),
i.e. M_n^2 in units of pi g^2.  (Equivalently M_n^2 = 2 pi g^2 lambda_n.)

------------------------------------------------------------------------------
Method (the standard Chebyshev / Multhopp-type spectral basis)
------------------------------------------------------------------------------
At alpha = 0 the solution behaves as phi(x) ~ sqrt(x(1-x)) at both ends
(endpoint exponent s = 1/2, from  pi s cot(pi s) = -alpha).  Map x = (1+xi)/2,
xi in [-1,1], and expand

    phi(x) = sqrt(1 - xi^2) * sum_m a_m U_m(xi),

where U_m is the Chebyshev polynomial of the 2nd kind.  These basis functions
are exact eigenfunctions of the hypersingular (finite-part) operator:

    FP-Int_{-1}^{1}  sqrt(1-eta^2) U_m(eta) / (xi - eta)^2  d eta = -pi (m+1) U_m(xi).

Rayleigh-Ritz / Galerkin projection onto the same basis turns the problem into
a generalized symmetric-definite eigenproblem

        D a = Lambda S a ,   Lambda = 2 pi^2 lambda,

    D_{nm} = pi^2 (m+1) delta_{nm},                              (diagonal)
    S_{nm} = Int_{-1}^{1} (1-xi^2) U_n(xi) U_m(xi) d xi
           = (1/2) [ f(n-m) - f(n+m+2) ],
           f(k) = 2/(1-k^2) if k even, else 0.

S is the (symmetric, positive-definite) Gram matrix of the basis functions in
L^2[0,1]; D is diagonal.  Both are block-diagonal under xi -> -xi, i.e. under
x -> 1-x.  Even-degree m  -> phi symmetric;  odd-degree m -> phi antisymmetric.
Those two blocks are the two CP sectors.  This is a variational method, so the
eigenvalues converge monotonically from above as the basis size grows.

Then  2 lambda_n = Lambda_n / pi^2.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh


PI = np.pi


def _f(k: np.ndarray | int):
    """f(k) = Int_0^pi sin(theta) cos(k theta) d theta.

    = 2/(1-k^2) for even integer k, 0 for odd integer k.
    """
    k = np.asarray(k)
    out = np.where(k % 2 == 0, 2.0 / (1.0 - k.astype(float) ** 2), 0.0)
    return out


def _odd_harmonic(pmax: int) -> np.ndarray:
    """Cumulative odd-harmonic numbers  H[p] = sum_{k=1}^{p} 1/(2k-1),  H[0]=0."""
    H = np.zeros(pmax + 1)
    for p in range(1, pmax + 1):
        H[p] = H[p - 1] + 1.0 / (2 * p - 1)
    return H


def _I_matrix(degrees: np.ndarray) -> np.ndarray:
    """Gram-type matrix of the *unweighted* Chebyshev-U's,

        I_{nm} = Int_{-1}^{1} U_n(xi) U_m(xi) d xi
               = 2 [ H(Q) - H(R) ],   Q = (n+m)/2 + 1,  R = |n-m|/2,   (n+m even)
               = 0                                                     (n+m odd)

    with H the cumulative odd-harmonic numbers.  This is the matrix that carries
    the quark-mass (potential) term  alpha (1/x + 1/(1-x)) into the basis:
    projecting that potential onto b_n = sqrt(1-xi^2) U_n gives exactly 4*I_{nm}
    (because (1-xi^2)(1/x + 1/(1-x)) = 4).  Derivation in the README.
    """
    d = np.asarray(degrees)
    di = d[:, None]
    dj = d[None, :]
    Q = (di + dj) // 2 + 1
    R = np.abs(di - dj) // 2
    H = _odd_harmonic(int(Q.max()))
    I = 2.0 * (H[Q] - H[R])
    return np.where((di + dj) % 2 == 0, I, 0.0)


def build_sector(num_basis: int, parity: int, alpha: float = 0.0):
    """Build (D, S) for one CP sector, where D already includes the quark-mass
    (potential) term.

    parity = 0 -> even-degree basis  (phi symmetric under x -> 1-x)
    parity = 1 -> odd-degree basis   (phi antisymmetric under x -> 1-x)

    alpha = pi m_q^2 / g^2  is the (equal-mass) FLZ quark-mass parameter; the
    DSSYK duality point is alpha = 0.  The potential contributes P = 4 alpha I,
    so the RHS operator becomes  D_hyp + P  with  D_hyp = pi^2 (m+1).

    Returns degrees, D (K x K, diagonal + potential), S (K x K).
    """
    degrees = parity + 2 * np.arange(num_basis)  # e.g. even: 0,2,4,...
    D = np.diag(PI ** 2 * (degrees + 1.0))
    if alpha != 0.0:
        D = D + 4.0 * alpha * _I_matrix(degrees)

    di = degrees[:, None]
    dj = degrees[None, :]
    S = 0.5 * (_f(di - dj) - _f(di + dj + 2))
    return degrees, D, S


def solve_sector(num_basis: int, parity: int, alpha: float = 0.0):
    """Return (two_lambda, eigenvectors, degrees) for one CP sector, ascending."""
    degrees, D, S = build_sector(num_basis, parity, alpha)
    # Generalized symmetric-definite eigenproblem  D a = Lambda S a.
    Lam, vecs = eigh(D, S)
    two_lambda = Lam / PI ** 2  # 2 lambda_n = M^2 / (pi g^2)
    order = np.argsort(two_lambda)
    return two_lambda[order], vecs[:, order], degrees


def spectrum(num_basis: int = 200, num_levels: int = 40, alpha: float = 0.0):
    """Full interleaved spectrum.

    alpha = pi m_q^2 / g^2 (equal-mass); alpha = 0 is the DSSYK duality point.

    Returns a list of dicts with keys:
      n         : level index (0,1,2,... in FLZ ordering)
      two_lambda: 2 lambda_n = M_n^2 / (pi g^2)
      parity    : 'even' (symmetric) or 'odd' (antisymmetric) under x->1-x
      CP        : (-1)^(n+1)   [Susskind et al. convention]
    """
    ev_sym, _, _ = solve_sector(num_basis, parity=0, alpha=alpha)   # symmetric -> FLZ even n
    ev_anti, _, _ = solve_sector(num_basis, parity=1, alpha=alpha)  # antisym   -> FLZ odd  n

    levels = []
    for n in range(num_levels):
        if n % 2 == 0:
            val = ev_sym[n // 2]
            par = "even"
        else:
            val = ev_anti[n // 2]
            par = "odd"
        levels.append(
            dict(n=n, two_lambda=val, parity=par, CP=(-1) ** (n + 1))
        )
    return levels


if __name__ == "__main__":
    levels = spectrum(num_basis=200, num_levels=30)
    print(f"{'n':>3} {'parity':>10} {'CP':>3}   2*lambda_n = M^2/(pi g^2)")
    print("-" * 55)
    for L in levels:
        print(f"{L['n']:>3} {L['parity']:>10} {L['CP']:>+3d}   "
              f"{L['two_lambda']:.14f}")
