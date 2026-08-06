#!/usr/bin/env python3
r"""
Matched-exponent Jacobi solver for the 't Hooft equation at alpha != 0
======================================================================

Closes the last documented limitation of this module (README "Precision
notes"): the main Chebyshev-U solver (thooft_spectrum.py) hardwires the
alpha = 0 endpoint exponent 1/2 into its basis sqrt(x(1-x)) U_m, so away from
the duality point its convergence is only algebraic (measured 1e-5..2e-4
residuals in the Litvinov-Meshcheriakov sum-rule tests).  This solver builds
the TRUE endpoint exponent into the basis and restores fast convergence at
alpha != 0.

The equation (FLZ arXiv:0905.2280 Eq. 1.1, equal masses alpha1 = alpha2 =
alpha = pi m_q^2 / g^2):

    2 pi^2 lambda phi(x) = alpha (1/x + 1/(1-x)) phi(x)
                           - FP-Int_0^1 phi(y)/(y-x)^2 dy ,      x in [0,1],

phi(0) = phi(1) = 0, eigenvalue reported as 2 lambda_n = M_n^2/(pi g^2),
exactly as in thooft_spectrum.py.

Endpoint exponent ('t Hooft 1974; FLZ Sec. 1; LM arXiv:2409.11324)
------------------------------------------------------------------
phi(x) ~ x^beta as x -> 0 (and (1-x)^beta at x -> 1) with

    pi beta cot(pi beta) = -alpha ,        beta in (0, 1).

Derivation (verified numerically in test_jacobi.py): substituting
phi = x^beta into the equation, the finite-part integral produces the local
singular term x^(beta-1) * FP-Int_0^inf t^beta/(t-1)^2 dt.  With
Int_0^inf t^beta (t+z)^-2 dt = z^(beta-1) * pi beta / sin(pi beta) and the
Hadamard finite part being the mean of the z -> e^{+-i pi} boundary values,

    FP-Int_0^inf t^beta/(t-1)^2 dt = -pi beta cot(pi beta),

so the x^(beta-1) terms cancel iff alpha + pi beta cot(pi beta) = 0.  At
alpha = 0 this gives beta = 1/2 exactly (the Chebyshev case); the root is
unique in (0,1) for alpha > -1 because pi b cot(pi b) decreases
monotonically from 1 to -inf there.  It is found with scipy brentq at
machine precision; the tests assert the transcendental residual at the root
and re-derive the half-line finite-part identity itself through an
independent subtracted scipy.quad regularization.

Basis and Gram ("mass") matrix
------------------------------
    b_k(x) = C_k [x(1-x)]^beta P_k^(2beta,2beta)(2x - 1) ,   k = 0, 1, ...

with P_k^(a,a) the (ultraspherical) Jacobi polynomials.  Weight choice: the
L^2[0,1] Gram matrix of functions [x(1-x)]^beta p_k(x) is
Int_0^1 [x(1-x)]^(2beta) p_n p_m dx, i.e. the polynomials are integrated
against the weight [x(1-x)]^(2beta) -- so the Jacobi family orthogonal
w.r.t. THAT weight (parameters (2beta, 2beta)) makes the Gram matrix exactly
diagonal, and with the normalization

    C_k^-2 = ||[x(1-x)]^beta P_k||^2
           = Gamma(k+2beta+1)^2 / [ (2k+4beta+1) Gamma(k+1) Gamma(k+4beta+1) ]

(closed form from the standard Jacobi norm h_k^(a,b), re-verified against
Gauss-Jacobi quadrature in the tests) it is exactly the identity: the
generalized eigenproblem collapses to a standard symmetric one, with perfect
conditioning.  Any other classical family (e.g. P^(beta,beta), which pairs
with known Cauchy-transform closed forms) leaves a full, increasingly
ill-conditioned Gram matrix; for 2beta not an integer no polynomial family
diagonalizes the finite-part operator itself (the Tricomi/Sohngen
construction needs the weight-exponent sum at the two endpoints to be an
integer, which is what makes alpha = 0, i.e. 2beta = 1, special).

Under x -> 1-x the basis splits into the two CP sectors by the parity of k,
exactly as in the main solver (P_k^(a,a)(-xi) = (-1)^k P_k^(a,a)(xi)).

Discretization of the finite-part operator (weak form)
------------------------------------------------------
The same exact Gagliardo identity used by sine_solver.py: for u, v vanishing
at the endpoints (any Holder exponent beta > 0),

    <u, FP v> = -B[u, v] - <u, V v>,      V(x) = 1/x + 1/(1-x) = 1/(x(1-x)),
    B[u, v]   = (1/2) Int Int [u(x)-u(y)][v(x)-v(y)] / (x-y)^2 dx dy ,

so the Galerkin matrix of the full operator alpha*V - FP is

    A_nm = (alpha + 1) <b_n, V b_m> + B[b_n, b_m] ,

manifestly SYMMETRIC (asserted at assembly), with the alpha -> alpha+1 bare-
mass shift appearing exactly as in sine_solver/eigenfunctions.  Eigenvalues:
A a = mu a with mu = 2 pi^2 lambda (Gram = identity).

* Potential term: <b_n, V b_m> = C_n C_m Int [x(1-x)]^(2beta-1) P_n P_m dx
  is a polynomial against the Jacobi weight (2beta-1, 2beta-1): evaluated by
  Gauss-Jacobi quadrature with enough nodes to be EXACT (up to rounding).

* Gagliardo double integral: the integrand g(x,y) = [Du][Dv]/(x-y)^2 is
  analytic across the diagonal (g -> u'v' as y -> x) but has algebraic
  branch points at the endpoints, so tensor Gauss-Legendre stalls.  We use a
  tensor TANH-SINH (double-exponential) rule, x = (1 + tanh((pi/2) sinh u))/2,
  |u| <= 4, which handles arbitrary endpoint exponents at spectral accuracy;
  the off-diagonal part is assembled through the graph-Laplacian form
  (Q_ij = w_i w_j/(x_i-x_j)^2, exactly symmetric), and the removable diagonal
  contributes (1/2) sum_i w_i^2 b_n'(x_i) b_m'(x_i).  1-x is stored as its
  own node array (the grid is symmetric) so no precision is lost at x -> 1.

  The node count must RESOLVE the basis: the degree-2K top polynomial
  oscillates with wavelength ~1/(4K) mid-interval, where the tanh-sinh
  spacing is (pi/4)h; undersampling aliases the quadrature and produces
  spurious eigenvalues (observed, violently, at K = 220 with 801 nodes).
  Default n_ts = max(501, 12K+1), measured safe for the lowest K/2 levels.

  Measured accuracy of this scheme (see jacobi_convergence.json and the
  anchor test): at beta = 1/2 the assembled B matrix can be compared against
  the EXACT closed form B[B_n, B_m] = (pi^2/2)(m+1) delta_nm - 2 I_nm that
  follows from the Chebyshev eigenrelation of the finite-part operator plus
  the odd-harmonic closed form I_nm of thooft_spectrum._I_matrix; at
  n_ts = 501 the worst entry error over a 21-function block is ~5e-13, and
  the alpha = 0 spectrum reproduces the K-equivalent Chebyshev solver to
  ~1e-12 (same trial space, so any difference IS pure quadrature error).

What is exact, measured, conjectured
------------------------------------
Exact: the endpoint-exponent condition; the Gagliardo identity; the Gram and
potential matrix elements (Gauss-Jacobi of polynomial integrands); the
beta = 1/2 anchor closed form.  Measured (jacobi_convergence.json): the
quadrature floor (~1e-12 on low eigenvalues at default settings); at
alpha = +/-0.5 the ground state converges below 1e-11 by K = 32 while the old
solver follows the exponent-mismatch algebraic law ~K^-(4 beta) (fitted rates
match 4 beta to 2 digits: 2.34, 1.48) and still sits at 1e-7..1e-5 at K = 128;
the reference's K = 200 vs 160 self-consistency is ~4e-12, externally certified
by the LM tables (their 6-digit floor) and by the 20-digit LM sum rules
(s = 2: <= 3e-8, s >= 4: ~1e-12 relative).  Caveat, stated honestly:
the true wavefunction also carries subleading endpoint behavior NOT in this
basis (secondary exponent beta2 in (1,2) solving the same transcendental
condition, plus integer-power/log terms sourced by the regular part of the
kernel), so alpha != 0 convergence, while fast, is not guaranteed
unbounded-spectral; at the K used here that structure is below the measured
1e-11..1e-12 floor.  Conjectured: nothing load-bearing.

Conventions match thooft_spectrum.py throughout: solve_sector_jacobi returns
(two_lambda ascending, eigenvectors, degrees); spectrum_jacobi interleaves
the two CP sectors with the same explicit interleaving assertion.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.special import gammaln, roots_jacobi

from thooft_spectrum import interleave_levels

PI = np.pi


# ---------------------------------------------------------------------------
# endpoint exponent
# ---------------------------------------------------------------------------
def endpoint_exponent(alpha: float) -> float:
    """beta in (0,1) solving  pi beta cot(pi beta) = -alpha  ('t Hooft 1974).

    alpha = pi m_q^2/g^2 (FLZ convention, this repo's alpha); requires
    alpha > -1 (above the chiral point).  alpha = 0 returns exactly 1/2.
    """
    if alpha == 0.0:
        return 0.5
    if alpha <= -1.0:
        raise ValueError("endpoint exponent needs alpha > -1 (chiral limit "
                         f"is alpha = -1); got alpha = {alpha}")

    def f(b):
        return PI * b * np.cos(PI * b) / np.sin(PI * b) + alpha

    # f(0+) = 1 + alpha > 0, f(1-) -> -inf, and f is strictly decreasing.
    return brentq(f, 1e-12, 1.0 - 1e-12, xtol=1e-16, rtol=8.9e-16)


# ---------------------------------------------------------------------------
# Jacobi polynomials P_k^(a,a), values and derivatives, by stable recurrence
# ---------------------------------------------------------------------------
def _jacobi_all(max_degree: int, a: float, xi: np.ndarray) -> np.ndarray:
    """P_k^(a,a)(xi) for k = 0..max_degree (rows k), three-term recurrence.

    2k(k+2a)(2k+2a-2) P_k = (2k+2a-1)(2k+2a)(2k+2a-2) xi P_{k-1}
                            - 2(k+a-1)^2 (2k+2a) P_{k-2}
    (general Jacobi recurrence at alpha = beta = a; cross-checked against
    scipy.special.eval_jacobi in the tests)."""
    xi = np.asarray(xi, dtype=float)
    P = np.zeros((max_degree + 1, xi.size))
    P[0] = 1.0
    if max_degree >= 1:
        P[1] = (a + 1.0) * xi
    for k in range(2, max_degree + 1):
        c0 = 2.0 * k * (k + 2 * a) * (2 * k + 2 * a - 2)
        c1 = (2 * k + 2 * a - 1) * (2 * k + 2 * a) * (2 * k + 2 * a - 2)
        c2 = 2.0 * (k + a - 1) ** 2 * (2 * k + 2 * a)
        P[k] = (c1 * xi * P[k - 1] - c2 * P[k - 2]) / c0
    return P


def _jacobi_all_deriv(max_degree: int, a: float, xi: np.ndarray) -> np.ndarray:
    """d/dxi P_k^(a,a)(xi) via  P_k' = ((k+2a+1)/2) P_{k-1}^(a+1,a+1)."""
    xi = np.asarray(xi, dtype=float)
    dP = np.zeros((max_degree + 1, xi.size))
    if max_degree >= 1:
        Pup = _jacobi_all(max_degree - 1, a + 1.0, xi)
        ks = np.arange(1, max_degree + 1)
        dP[1:] = 0.5 * (ks + 2 * a + 1)[:, None] * Pup
    return dP


def _log_norm_sq(degrees: np.ndarray, beta: float) -> np.ndarray:
    """log ||b_k||^2 for b_k = [x(1-x)]^beta P_k^(2beta,2beta)(2x-1) on [0,1].

    ||b_k||^2 = Gamma(k+2b+1)^2 / [(2k+4b+1) Gamma(k+1) Gamma(k+4b+1)]
    (from the standard Jacobi norm; asserted against Gauss-Jacobi quadrature
    in test_jacobi.py)."""
    k = np.asarray(degrees, dtype=float)
    return (2.0 * gammaln(k + 2 * beta + 1) - np.log(2 * k + 4 * beta + 1)
            - gammaln(k + 1.0) - gammaln(k + 4 * beta + 1))


# ---------------------------------------------------------------------------
# tanh-sinh nodes on (0,1)
# ---------------------------------------------------------------------------
def _tanh_sinh_01(n_ts: int, u_max: float = 4.0):
    """Symmetric tanh-sinh rule: x_j, (1-x)_j (exact by symmetry), w_j.

    x(u) = (1 + tanh((pi/2) sinh u))/2 on u in [-u_max, u_max], n_ts points.
    1-x is returned as its own array (= x reversed) so that quantities near
    x = 1 keep full relative precision.  u_max = 4 puts the extreme nodes at
    x ~ 5e-38 -- far below any endpoint scale that matters here, while all
    intermediates (w^(beta-1) etc.) stay inside double range for beta > 0.05.
    """
    if n_ts % 2 == 0:
        n_ts += 1  # keep the grid symmetric about x = 1/2
    u = np.linspace(-u_max, u_max, n_ts)
    h = u[1] - u[0]
    s = np.sinh(u)
    # x = 1/(1 + exp(-pi*s)) computed stably from the small side
    x = np.empty_like(u)
    neg = s < 0
    x[neg] = np.exp(PI * s[neg]) / (1.0 + np.exp(PI * s[neg]))
    x[~neg] = 1.0 / (1.0 + np.exp(-PI * s[~neg]))
    omx = x[::-1].copy()                       # 1 - x, exact by symmetry
    # dx/du = (pi/2) cosh(u) sech^2((pi/2) sinh u); sech^2(z) = 4 e^{-2|z|}/(1+e^{-2|z|})^2
    z = (PI / 2.0) * s
    e = np.exp(-2.0 * np.abs(z))
    sech2 = 4.0 * e / (1.0 + e) ** 2
    w = h * (PI / 4.0) * np.cosh(u) * sech2
    return x, omx, w


# ---------------------------------------------------------------------------
# matrix assembly
# ---------------------------------------------------------------------------
def gagliardo_from_values(Psi: np.ndarray, dPsi: np.ndarray,
                          x: np.ndarray, omx: np.ndarray,
                          w: np.ndarray) -> np.ndarray:
    """B_nm = (1/2) IntInt [Db_n][Db_m]/(x-y)^2  from basis values Psi[n,i]
    and derivatives dPsi[n,i] on nodes (x, 1-x, w).

    Off-diagonal pairs via the graph-Laplacian form (exactly symmetric);
    the removable diagonal (limit u'v') contributes (1/2) w_i^2 u'_i v'_i.
    This is precisely the tensor quadrature of the smooth integrand g(x,y),
    so its accuracy is that of the underlying 1-d rule.

    The pairwise differences x_i - x_j are formed from whichever endpoint the
    pair is closer to (x_i - x_j vs (1-x_j) - (1-x_i)): deep tanh-sinh nodes
    near x = 1 round to exactly 1.0 in double, so the naive difference there
    is 0 and would produce inf/nan; the (1-x)-side difference keeps full
    relative precision."""
    near1 = (x[:, None] + x[None, :]) > 1.0
    dx = np.where(near1, omx[None, :] - omx[:, None], x[:, None] - x[None, :])
    np.fill_diagonal(dx, 1.0)                    # placeholder, zeroed below
    Q = np.outer(w, w) / (dx * dx)
    np.fill_diagonal(Q, 0.0)
    r = Q.sum(axis=1)
    off = (Psi * r) @ Psi.T - Psi @ Q @ Psi.T
    diag = (dPsi * (w ** 2)) @ dPsi.T
    return off + 0.5 * diag


def _auto_n_ts(num_basis: int) -> int:
    """Tanh-sinh node count that resolves the basis oscillation.

    The highest basis polynomial has degree ~2K, oscillating with wavelength
    ~1/(4K) in mid-interval where the tanh-sinh spacing is (pi/4)h; sampling
    below that aliases the Gagliardo quadrature and produces spurious
    eigenvalues.  Measured at alpha = 0 against the exact Chebyshev solver:
    n_ts = 10K is safe for the lowest K/2 levels (~1e-11); 12K adds margin."""
    return max(501, 12 * num_basis + 1)


def build_sector_jacobi(num_basis: int, parity: int, alpha: float = 0.0,
                        n_ts: int | None = None, beta: float | None = None):
    """Assemble the (symmetric) operator matrix A for one CP sector.

    A_nm = (alpha+1) <b_n, V b_m> + B[b_n, b_m],  Gram = identity, so the
    eigenproblem is standard:  A a = mu a,  mu = 2 pi^2 lambda.

    parity = 0 -> even-degree polynomials (phi symmetric under x -> 1-x),
    parity = 1 -> odd-degree (antisymmetric); identical labeling to
    thooft_spectrum.build_sector.  beta defaults to endpoint_exponent(alpha)
    (overridable for validation experiments only).  n_ts = None picks the
    resolution-matched _auto_n_ts(num_basis).

    Returns degrees, A.
    """
    if beta is None:
        beta = endpoint_exponent(alpha)
    if n_ts is None:
        n_ts = _auto_n_ts(num_basis)
    degrees = parity + 2 * np.arange(num_basis)
    dmax = int(degrees[-1])
    Cln = -0.5 * _log_norm_sq(degrees, beta)     # log normalization constants

    # --- potential term, exact by Gauss-Jacobi (2beta-1, 2beta-1) ----------
    nq = 2 * dmax + 8
    xq, wq = roots_jacobi(nq, 2 * beta - 1.0, 2 * beta - 1.0)
    Pq = _jacobi_all(dmax, 2.0 * beta, xq)[degrees]
    Pq *= np.exp(Cln)[:, None]
    pot = 2.0 ** (1.0 - 4.0 * beta) * ((Pq * wq) @ Pq.T)

    # --- Gagliardo term on tanh-sinh nodes ---------------------------------
    x, omx, w = _tanh_sinh_01(n_ts)
    xi = x - omx                                  # 2x-1, exact near both ends
    lw = np.log(x) + np.log(omx)                  # log of x(1-x)
    wbeta = np.exp(beta * lw)
    wbm1 = np.exp((beta - 1.0) * lw)
    P = _jacobi_all(dmax, 2.0 * beta, xi)[degrees]
    dP = _jacobi_all_deriv(dmax, 2.0 * beta, xi)[degrees]
    C = np.exp(Cln)[:, None]
    Psi = C * (wbeta * P)
    dPsi = C * (beta * wbm1 * (omx - x) * P + 2.0 * wbeta * dP)
    B = gagliardo_from_values(Psi, dPsi, x, omx, w)

    A = (alpha + 1.0) * pot + B
    # the weak form is symmetric; assert the discrete matrix is too, then
    # symmetrize away the last-bit rounding
    scale = np.abs(A).max()
    asym = np.abs(A - A.T).max()
    if not asym <= 1e-12 * scale:
        raise RuntimeError(f"assembled matrix not symmetric: {asym:.2e} "
                           f"vs scale {scale:.2e}")
    A = 0.5 * (A + A.T)
    return degrees, A


def solve_sector_jacobi(num_basis: int, parity: int, alpha: float = 0.0,
                        n_ts: int | None = None, beta: float | None = None,
                        eigvals_only: bool = False):
    """(two_lambda ascending, eigenvectors, degrees) for one CP sector.

    Same return convention as thooft_spectrum.solve_sector; eigenvectors are
    orthonormal coefficient columns in the unit-normalized b_k basis (the
    Gram matrix is the identity).  eigvals_only=True skips the eigenvectors
    (None in that slot, ~2x faster)."""
    degrees, A = build_sector_jacobi(num_basis, parity, alpha, n_ts, beta)
    if eigvals_only:
        mu = eigh(A, eigvals_only=True)
        return np.sort(mu / PI ** 2), None, degrees
    mu, vecs = eigh(A)
    two_lambda = mu / PI ** 2
    order = np.argsort(two_lambda)
    return two_lambda[order], vecs[:, order], degrees


def spectrum_jacobi(num_basis: int = 80, num_levels: int = 30,
                    alpha: float = 0.0, n_ts: int | None = None):
    """Full interleaved spectrum, same output schema as
    thooft_spectrum.spectrum (n, two_lambda, parity, CP), with the same
    explicit interleaving assertion instead of silent mislabeling (shared
    helper thooft_spectrum.interleave_levels)."""
    ev_sym, _, _ = solve_sector_jacobi(num_basis, 0, alpha, n_ts,
                                       eigvals_only=True)
    ev_anti, _, _ = solve_sector_jacobi(num_basis, 1, alpha, n_ts,
                                        eigvals_only=True)
    return interleave_levels(ev_sym, ev_anti, num_levels, alpha)


# ---------------------------------------------------------------------------
# spectral sum rules with a drift-corrected analytic tail
# ---------------------------------------------------------------------------
def sum_rule_G(two_lambda: np.ndarray, parity_offset: int, s: int,
               n_cut: int = 120, n_fit: int = 16,
               n_big: int = 40001) -> float:
    """G^(s) = sum_n 1/lambda_n^s over one parity sector (FLZ/LM notation;
    lambda_n = two_lambda_n / 2, so each term is 2^s / two_lambda_n^s).

    Levels with FLZ index n <= n_cut are summed from the computed spectrum;
    the tail uses the DRIFT-CORRECTED asymptotic

        2 lambda_n ~ n + 3/4 + delta(n),   delta(n) = A + B log(n + 3/4),

    with (A, B) least-squares-fitted to the next n_fit computed levels above
    the cutoff (at alpha != 0 the spectrum drifts logarithmically off the
    n + 3/4 line -- the 't Hooft/WKB log correction, clearly visible in the
    LM tables -- so the plain polygamma tail of the alpha = 0 tests is wrong
    at ~1e-4 for s = 2).  The tail is summed numerically to n_big with the
    fitted drift, then closed with a Hurwitz-zeta remainder at frozen drift.

    Measured accuracy against the 20-digit LM closed forms (K = 150 per
    sector, all five LM alpha values; table in jacobi_convergence.json):
    s = 2: <= 3e-8, s = 3: <= 2e-10, s = 4..7: <= 2e-12 relative.  The same
    machinery at alpha = 0 reproduces FLZ 7 zeta(3) and 2 to ~6e-11.
    """
    from scipy.special import zeta as _hurwitz

    evs = np.asarray(two_lambda, dtype=float)
    ns = parity_offset + 2 * np.arange(len(evs))
    head = ns <= n_cut
    total = float(np.sum(2.0 ** s / evs[head] ** s))

    sel = np.where(~head)[0][:n_fit]
    if len(sel) < 4:
        raise ValueError("not enough levels above n_cut to fit the drift; "
                         "increase num_basis or lower n_cut")
    resid = evs[sel] - (ns[sel] + 0.75)
    X = np.vstack([np.ones(len(sel)), np.log(ns[sel] + 0.75)]).T
    A, B = np.linalg.lstsq(X, resid, rcond=None)[0]

    n_start = ns[head][-1] + 2
    ntail = np.arange(n_start, n_big, 2, dtype=float)
    drift = A + B * np.log(ntail + 0.75)
    total += float(np.sum(2.0 ** s / (ntail + 0.75 + drift) ** s))
    a_frozen = (ntail[-1] + 2 + 0.75 + drift[-1]) / 2.0
    total += float(_hurwitz(s, a_frozen))
    return total


# ---------------------------------------------------------------------------
# demonstration / convergence artifact
# ---------------------------------------------------------------------------
def _convergence_scan(alphas=(0.5, -0.5), Ks=(8, 16, 32, 64, 128),
                      K_ref=200):
    """Error vs basis size K, old (Chebyshev) solver vs this solver, for
    FLZ levels n = 0 and n = 4 (symmetric-sector entries 0 and 2), against a
    high-K Jacobi reference.  Also fits the old solver's algebraic rate
    err ~ K^-p over the scan.  Returns a JSON-ready dict."""
    from thooft_spectrum import solve_sector as solve_sector_cheb

    out = {}
    for alpha in alphas:
        beta = endpoint_exponent(alpha)
        ref, _, _ = solve_sector_jacobi(K_ref, 0, alpha)
        ref2, _, _ = solve_sector_jacobi(K_ref - 40, 0, alpha)
        ref_shift = float(np.max(np.abs(ref[:6] - ref2[:6])))
        rows = []
        for K in Ks:
            old, _, _ = solve_sector_cheb(K, parity=0, alpha=alpha)
            jac, _, _ = solve_sector_jacobi(K, 0, alpha)
            rows.append(dict(
                K=K,
                old_err_n0=float(abs(old[0] - ref[0])),
                jac_err_n0=float(abs(jac[0] - ref[0])),
                old_err_n4=float(abs(old[2] - ref[2])),
                jac_err_n4=float(abs(jac[2] - ref[2])),
            ))
        lk = np.log([r["K"] for r in rows])
        le = np.log([r["old_err_n0"] for r in rows])
        p = -np.polyfit(lk, le, 1)[0]
        out[str(alpha)] = dict(
            beta=beta,
            reference_K=K_ref,
            reference_self_consistency=ref_shift,
            reference_two_lambda_n0=float(ref[0]),
            reference_two_lambda_n4=float(ref[2]),
            old_solver_algebraic_rate_n0=float(p),
            four_beta=float(4 * beta),
            rows=rows,
        )
    return out


def _sum_rule_scan(K=150):
    """Residuals of sum_rule_G against the 20-digit LM closed forms for
    s = 2..7 at all five LM alpha values -- the measurement that sets the
    tolerances asserted in test_jacobi.py.  Returns a JSON-ready dict."""
    import json
    import os

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "external_anchors.json")) as fh:
        anchors = json.load(fh)

    out = {}
    for a_str, exact in anchors["lm_sum_rules_exact"].items():
        if a_str == "comment":
            continue
        alpha = float(a_str)
        ev_s, _, _ = solve_sector_jacobi(K, 0, alpha)
        ev_a, _, _ = solve_sector_jacobi(K, 1, alpha)
        row = {}
        for s in range(2, 8):
            gp = sum_rule_G(ev_s, 0, s)
            gm = sum_rule_G(ev_a, 1, s)
            tp, tm = float(exact[f"Gp{s}"]), float(exact[f"Gm{s}"])
            row[f"s{s}_abs"] = max(abs(gp - tp), abs(gm - tm))
            row[f"s{s}_rel"] = max(abs(gp - tp) / abs(tp),
                                   abs(gm - tm) / abs(tm))
        out[a_str] = row
    return out


if __name__ == "__main__":
    import json

    print("Matched-exponent Jacobi solver for the 't Hooft equation\n")

    print("Endpoint exponents beta(alpha)  [pi b cot(pi b) = -alpha]:")
    for a in (-0.5, -0.1, 0.0, 0.1, 0.5, 1.0):
        print(f"  alpha = {a:>4}:  beta = {endpoint_exponent(a):.15f}")

    print("\nalpha = 0 sanity (same trial space as the Chebyshev solver;")
    print("difference = pure quadrature error of the Gagliardo assembly):")
    from thooft_spectrum import solve_sector as solve_sector_cheb
    for K in (20, 40, 60):
        evj, _, _ = solve_sector_jacobi(K, 0, 0.0)
        evc, _, _ = solve_sector_cheb(K, parity=0)
        d = np.max(np.abs(evj[:K // 2] - evc[:K // 2]))
        print(f"  K = {K:>3}: max |jacobi - chebyshev| over lowest K/2 "
              f"levels = {d:.2e}")

    print("\nSpectrum at alpha = 0.5 (2 lambda_n = M^2/(pi g^2)):")
    for L in spectrum_jacobi(num_basis=80, num_levels=8, alpha=0.5):
        print(f"  n = {L['n']}  {L['parity']:>5}  CP {L['CP']:+d}   "
              f"{L['two_lambda']:.12f}")

    print("\nConvergence scan (old solver vs Jacobi, alpha = +/-0.5)...")
    scan = _convergence_scan()
    print("Sum-rule residual scan (s = 2..7, five LM alphas)...")
    sums = _sum_rule_scan()
    payload = {
        "_provenance": {
            "generated_by": "jacobi_solver.py __main__",
            "reference": "high-K Jacobi solve (K=200, auto n_ts) per alpha; "
                         "reference_self_consistency = max |K=200 - K=160| "
                         "change over the lowest 6 levels of the symmetric "
                         "sector (self-convergence estimate, mildly "
                         "optimistic; externally certified by the LM tables "
                         "and 20-digit LM sum rules in test_jacobi.py)",
            "errors": "absolute errors on 2*lambda = M^2/(pi g^2); n0/n4 = "
                      "FLZ levels 0 and 4 (symmetric sector entries 0, 2)",
            "old_solver_rate": "log-log fit of old_err_n0 over the K scan; "
                               "compare four_beta: the measured algebraic "
                               "rate is ~K^-(4 beta), the exponent-mismatch "
                               "prediction",
            "sum_rules": "residuals of sum_rule_G (K=150/sector, "
                         "drift-corrected tail) vs the 20-digit LM closed "
                         "forms in external_anchors.json; these MEASURED "
                         "values set the tolerances asserted in "
                         "test_jacobi.py",
        },
        "old_vs_jacobi": scan,
        "lm_sum_rule_residuals": sums,
    }
    from pathlib import Path
    with open(Path(__file__).resolve().parent / "jacobi_convergence.json",
              "w") as fh:
        json.dump(payload, fh, indent=1)
    print("wrote jacobi_convergence.json")
    for a_str, blk in scan.items():
        print(f"\nalpha = {a_str} (beta = {blk['beta']:.6f}, reference "
              f"self-consistency {blk['reference_self_consistency']:.1e}, "
              f"old rate K^-{blk['old_solver_algebraic_rate_n0']:.2f} vs "
              f"4 beta = {blk['four_beta']:.2f}):")
        print(f"{'K':>5} {'old n0':>10} {'jacobi n0':>10} "
              f"{'old n4':>10} {'jacobi n4':>10}")
        for r in blk["rows"]:
            print(f"{r['K']:>5} {r['old_err_n0']:>10.2e} "
                  f"{r['jac_err_n0']:>10.2e} {r['old_err_n4']:>10.2e} "
                  f"{r['jac_err_n4']:>10.2e}")
    print("\nSum-rule residuals (worst of G+/G-, absolute):")
    hdr = "  ".join(f"s={s}" for s in range(2, 8))
    print(f"{'alpha':>8}   " + hdr)
    for a_str, row in sums.items():
        vals = "  ".join(f"{row[f's{s}_abs']:.0e}" for s in range(2, 8))
        print(f"{a_str:>8}   " + vals)

