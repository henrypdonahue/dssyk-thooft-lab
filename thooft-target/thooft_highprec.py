#!/usr/bin/env python3
r"""
Arbitrary-precision massless 't Hooft spectrum, in mpmath.
==========================================================

The same Rayleigh-Ritz / Chebyshev-U method as thooft_spectrum.py, but
carried out in arbitrary-precision arithmetic.  The working precision is
then not the floor; the basis truncation is.  Measured at alpha = 0 the
truncation error falls like K^-6 (2 lambda_0 moves by 5.9e-11, 5.2e-12,
4.1e-13 over K = 20 -> 30 -> 45 -> 60), so at the K = 90 the reference
builder runs the eigenvalues are good to ~3e-15 -- one digit past double
-- and 1e-18 would need K ~ 400 (tens of minutes per sector in mpmath).

There is a cheaper way to the same place, and past it.  rayleigh_polish
takes a DOUBLE-precision Ritz vector and evaluates its Rayleigh quotient
with the exact rational matrix elements in mpmath: second order in the
vector's error, O(K^2) per level.  polished_levels and extrapolate_ladder
turn that into the 40-level, 17-21 digit table
reference_spectrum_polished.json (generate_reference.py --polish).  The
full eigensolve below stays as the independent cross-check: it shares no
eigenvector with the polished route.

The matrix elements D and S are RATIONAL, times pi^2 in D.  So the only
precision limits are the working precision and the basis-truncation
error.  Both are controllable.

The generalized eigenproblem is D a = Lambda S a, with S symmetric
positive definite.

At alpha = 0, the duality point, D is diagonal and positive.  That reduces
the problem by a single O(K^2) scaling to the standard symmetric

    A = D^{-1/2} S D^{-1/2},   A = A^T,   eig(A) = {1/Lambda_n},

with no Cholesky and no matrix inversion.  Then
2 lambda_n = Lambda_n / pi^2.

For alpha != 0, D is not diagonal, and a Cholesky reduction is used
instead.
"""

from mpmath import mp, mpf, matrix, cholesky, eigsy, pi, nstr


def _f(k):
    """f(k) = 2/(1-k^2) if k even else 0  (exact rational in mpf)."""
    if k % 2 != 0:
        return mpf(0)
    return mpf(2) / (mpf(1) - mpf(k) ** 2)


def _odd_harmonic_mp(pmax):
    """H[p] = sum_{k=1}^{p} 1/(2k-1) in mpf; H[0]=0."""
    H = [mpf(0)] * (pmax + 1)
    for p in range(1, pmax + 1):
        H[p] = H[p - 1] + mpf(1) / (2 * p - 1)
    return H


def _I_mp(di, dj, H):
    """I_{nm} = Int_{-1}^1 U_n U_m dxi = 2[H(Q)-H(R)] for n+m even, else 0."""
    if (di + dj) % 2 != 0:
        return mpf(0)
    Q = (di + dj) // 2 + 1
    R = abs(di - dj) // 2
    return 2 * (H[Q] - H[R])


def build_sector_mp(num_basis, parity, alpha=0):
    """Build (degrees, D, S) for one CP sector.

    alpha = pi m_q^2 / g^2 (equal-mass FLZ parameter).  alpha = 0 is the
    duality point.  The quark-mass potential enters D as  P = 4 alpha I.
    """
    degrees = [parity + 2 * i for i in range(num_basis)]
    K = num_basis
    alpha = mpf(alpha)
    D = matrix(K, K)
    S = matrix(K, K)
    H = _odd_harmonic_mp(degrees[-1] + 1) if alpha != 0 else None
    for i in range(K):
        D[i, i] = pi ** 2 * (degrees[i] + 1)
    if alpha != 0:
        for i in range(K):
            for j in range(K):
                D[i, j] += 4 * alpha * _I_mp(degrees[i], degrees[j], H)
    for i in range(K):
        for j in range(K):
            di, dj = degrees[i], degrees[j]
            S[i, j] = (mpf(1) / 2) * (_f(di - dj) - _f(di + dj + 2))
    return degrees, D, S


def solve_sector_mp(num_basis, parity, alpha=0):
    """Return the sorted 2*lambda_n for one CP sector at current mp.dps."""
    _, D, S = build_sector_mp(num_basis, parity, alpha)
    K = num_basis
    if alpha == 0:
        # D diagonal & positive: reduce  D a = Lambda S a  by the O(K^2) scaling
        # A = D^{-1/2} S D^{-1/2}, whose eigenvalues are {1/Lambda_n}.
        d = [mp.sqrt(D[i, i]) for i in range(K)]
        A = matrix(K, K)
        for i in range(K):
            for j in range(K):
                A[i, j] = S[i, j] / (d[i] * d[j])
        E = eigsy(A, eigvals_only=True)  # skip the unused eigenvectors (~2x)
        return sorted(1 / (E[i] * pi ** 2) for i in range(K))
    # alpha != 0: D non-diagonal; reduce via Cholesky S = L L^T,
    # A = L^{-1} D L^{-T} (symmetric), eig(A) = {Lambda_n}.
    Linv = _invert_lower(cholesky(S))
    E = eigsy(Linv * D * Linv.T, eigvals_only=True)
    return sorted(E[i] / pi ** 2 for i in range(K))


def _invert_lower(L):
    """Invert a lower-triangular mpmath matrix by forward substitution."""
    n = L.rows
    Inv = matrix(n, n)
    for col in range(n):
        Inv[col, col] = mpf(1) / L[col, col]
        for i in range(col + 1, n):
            s = mpf(0)
            for k in range(col, i):
                s += L[i, k] * Inv[k, col]
            Inv[i, col] = -s / L[i, i]
    return Inv


# ---------------------------------------------------------------------------
# Rayleigh-quotient polish of double-precision Ritz vectors (alpha = 0)
# ---------------------------------------------------------------------------
def rayleigh_polish(vec, parity, dps=40):
    """2*lambda from ONE double-precision Ritz vector of the alpha = 0
    problem: the Rayleigh quotient

        2 lambda = a^T D a / (pi^2 a^T S a)

    with the exact rational matrix elements summed in mpmath.  The quotient
    is stationary at an eigenvector, so its error is SECOND order in the
    vector's: a LAPACK vector good to 1e-13 gives the eigenvalue of the
    K-truncated problem to ~1e-26, at O(K^2) cost per level instead of the
    O(K^3) of a full mpmath eigensolve.  What remains is basis truncation,
    which polished_levels removes by a K-ladder.

    vec: coefficient column from thooft_spectrum.solve_sector (any
    normalization); parity: 0 or 1.  Returns an mpf at the caller's
    precision `dps`."""
    saved = mp.dps
    mp.dps = dps
    try:
        K = len(vec)
        a = [mpf(float(v)) for v in vec]
        deg = [parity + 2 * i for i in range(K)]
        # S_ij = (1/2)[f(d_i - d_j) - f(d_i + d_j + 2)], f(k) = 2/(1 - k^2):
        # a Toeplitz piece in i - j and a Hankel piece in i + j
        toep = [mpf(1) / (1 - mpf(2 * k) ** 2) for k in range(K)]
        hank = [mpf(1) / (1 - mpf(2 * (k + 1 + parity)) ** 2)
                for k in range(2 * K - 1)]
        aSa = mpf(0)
        for i in range(K):
            row = [toep[abs(i - j)] - hank[i + j] for j in range(K)]
            aSa += a[i] * mp.fdot(row, a)
        aDa = mp.fsum((deg[i] + 1) * a[i] ** 2 for i in range(K))
        return +(aDa / aSa)            # the pi^2 of D cancels 1/pi^2
    finally:
        mp.dps = saved


def polished_levels(num_basis, num_levels, dps=40):
    """The first num_levels values of 2*lambda_n (FLZ order, sectors
    interleaved) of the K = num_basis truncation, each polished from the
    double-precision eigenvector."""
    from thooft_spectrum import solve_sector
    vecs = {par: solve_sector(num_basis, parity=par)[1] for par in (0, 1)}
    return [rayleigh_polish(vecs[n % 2][:, n // 2], n % 2, dps)
            for n in range(num_levels)]


def extrapolate_ladder(v1, v2, v3, K1, K2, K3):
    """Truncation extrapolation of one level from three basis sizes: fit
    v(K) = v_inf + c K^-r, return (v_inf, r, correction) with correction
    = v_inf - v3.  The variational values decrease with K, so the
    correction is negative; its magnitude is the honest error bar."""
    saved = mp.dps
    mp.dps = max(mp.dps, 40)
    try:
        ratio = (v1 - v2) / (v2 - v3)
        x1, x2, x3 = mp.log(K1), mp.log(K2), mp.log(K3)

        def g(r):
            return ((mp.e ** (-r * x1) - mp.e ** (-r * x2))
                    / (mp.e ** (-r * x2) - mp.e ** (-r * x3)) - ratio)
        r = mp.findroot(g, 6)
        corr = -(v2 - v3) / ((mpf(K3) / K2) ** r - 1)
        return v3 + corr, r, corr
    finally:
        mp.dps = saved


if __name__ == "__main__":
    mp.dps = 35
    K = 90
    print(f"Arbitrary precision (mp.dps={mp.dps}, basis K={K} per sector)\n")
    print("NOTE: the accuracy floor is the basis size K, not mp.dps (matrix")
    print("elements are exact rationals).  Certified table: reference_spectrum.json.\n")

    sym = solve_sector_mp(K, parity=0)
    anti = solve_sector_mp(K, parity=1)

    print(f"{'n':>3} {'parity':>6} {'CP':>3}   2*lambda_n = M^2/(pi g^2)")
    print("-" * 62)
    for n in range(20):
        if n % 2 == 0:
            v, par = sym[n // 2], "sym"
        else:
            v, par = anti[n // 2], "anti"
        cp = (-1) ** (n + 1)
        print(f"{n:>3} {par:>6} {cp:>+3d}   {nstr(v, 22)}")

    print(f"\nGround state (K={K}):")
    print(f"  2*lambda_0 = {nstr(sym[0], 20)}")
    from flz_reference import FLZ_PADE_LAMBDA0
    print(f"  FLZ Pade   = {FLZ_PADE_LAMBDA0:.15f}  (FLZ Eq. 4.37)")
    print("  reference  = reference_spectrum.json, ~15 digits (see validate.py)")
