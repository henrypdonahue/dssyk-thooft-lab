#!/usr/bin/env python3
r"""
Arbitrary-precision massless 't Hooft spectrum (mpmath).

Same Rayleigh-Ritz / Chebyshev-U method as thooft_spectrum.py, but carried out
in arbitrary-precision arithmetic so the eigenvalues can be pushed well past the
double-precision floor (~1e-13).  The matrix elements D and S are *rational*
(times pi^2 in D), so the only precision limit is the working precision and the
basis-truncation error, both controllable.

Generalized eigenproblem  D a = Lambda S a  (S symmetric positive definite).
At alpha = 0 (the duality point) D is diagonal & positive, so this reduces by a
single O(K^2) scaling to the standard symmetric
    A = D^{-1/2} S D^{-1/2},   A = A^T,   eig(A) = {1/Lambda_n},
with no Cholesky and no matrix inversion.  Then  2 lambda_n = Lambda_n / pi^2.
(For alpha != 0, D is non-diagonal and a Cholesky reduction is used instead.)
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

    alpha = pi m_q^2 / g^2 (equal-mass FLZ parameter); alpha = 0 is the duality
    point.  The quark-mass potential enters D as  P = 4 alpha I.
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


if __name__ == "__main__":
    mp.dps = 35
    K = 90
    print(f"Arbitrary precision (mp.dps={mp.dps}, basis K={K} per sector)\n")
    print("NOTE: with a truncated basis the accuracy floor is set by the basis")
    print("size K, NOT by mp.dps -- the matrix elements are exact rationals, so")
    print("arithmetic is not the bottleneck.  At this modest K the last couple of")
    print("printed digits are basis-truncation noise.  The certified table is")
    print("reference_spectrum.json (built at K=160 vs K=220 in generate_reference.py).\n")

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
    print("  FLZ Pade   = 0.737061746292690  (16 digits published)")
    print("  reference  = 0.737061746292690  (reference_spectrum.json, ~15 digits)")
