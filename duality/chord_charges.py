#!/usr/bin/env python3
r"""
The conserved-channel theorem: why the pair channel is frozen at N = inf.
=========================================================================

`chord_moments.py` found this numerically, at machine precision.  The
singlet bilinear channels A_n = sum_i psi_i (ad_H)^n psi_i carry no
propagating weight at N = infinity, in their i != j pair amplitude.  This
module proves it, and supplies the apparatus the tests assert.

THE PICTURE.  At N = infinity the bulk is a length.  |n> means the bulk
is n chords long, and T lengthens or shortens it.  A particle is a chord
stretched across the bulk; it gets dragged as the bulk moves, and that is
what puts weight at nonzero frequency.  You build A_n out of the
Hamiltonian.  The result below is that such a probe is not a particle.
It is the bulk.

SETUP (formal, unscaled chord basis; T|n> = |n+1> + [n]_q |n-1>).  In a
sequential pair word  H^a psi_i H^b psi_i H^c psi_j H^d psi_j  the i-pair
factors through the zero-chord sector as the DRESSED PROPAGATOR

    C(m) = close o T_M^m o open ,

where open |n> = |n, 0>, T_M is the one-matter transfer matrix, and
close |n0, n1> = q^(D n1) |n0 + n1>  (D = the matter weight Delta_psi).
C(m) is the whole block psi H^m psi, seen from the chord line.

LEMMA (the three-term recursion).  Define the j-fold closes

    close_j |n0, n1> = [n1][n1-1]...[n1-j+1] q^(D (n1-j)) |n0+n1-j> .

The chord rules give, with SCALAR coefficients,

    close_j o T_M = q^(D+j) T o close_j + [j]_q close_{j-1}
                    + (1 - q^(2D+j)) close_{j+1} .

The coefficients being scalars is the whole result.  The two matching
conditions collapse to [n1+1] - q^j [n1-j+1] = [j] and
q^(-D) = q^(D+j) + d_j q^(-D).  Neither depends on the chord occupation.
Had they, the family would not close, and the channel would propagate
like anything else.

Iterate from close_0.  Use close_{j>0} o open = 0 and close_0 o open = 1:

    C(m) = P_m(T) ,

a POLYNOMIAL in T: the Motzkin-path polynomial of length m.  At height j
it carries flat weight q^(D+j) T, down weight [j]_q, and up weight
1 - q^(2D+j).  First cases, verified exactly:

    C(1) = q^D T ,          C(2) = q^(2D) T^2 + (1 - q^(2D)) .

THEOREM (conserved channels).  The chord representation of the A_n
insertion is

    Xhat_n = sum_r (-1)^r C(n, r) C(n-r) T^r      [= d^n/dt^n C(t) e^{-tT}|_0]

-- a polynomial in T by the lemma.  Hence ad_T^k(Xhat_n) = 0.  EVERY
spectral moment mu_k (k >= 1, even and odd, both operator orderings) of
every A_n PAIR channel therefore vanishes at N = infinity, at every
lambda and every matter weight.  Explicitly:

    Xhat_1 = -(1 - q^D) T ,
    Xhat_2 = (1 - q^D)^2 T^2 + (1 - q^(2D)) .

Xhat_1 is proportional to T.  That is the N = infinity image of the exact
finite-N law A_1 = -2pH, so the n = 1 case lands on something known
independently.  Xhat_2 is the closed form of the A_2 chord image; note it
is conserved for the trivial reason that any polynomial in H is, so it is
not a new symmetry.  The content is that the channel reduces to such a
polynomial.

WHAT DOES PROPAGATE.  A random operator O_Delta inserts q^(Delta n-hat),
which does NOT commute with T because T changes n.  So random matter
propagates, with mu_2 = 2(1 - q^Delta).  The bulk holds particles
perfectly well.  The discriminator is not the size of the operator.  It
is whether its couplings are correlated with H's.

SCOPE, AND IT IS LOAD-BEARING.  This is the i != j pair amplitude at
strict N = infinity -- exactly the object chord_moments.py computes.  The
full channel also has a flavor-diagonal part of the SAME order, and that
part does not vanish for n >= 2.  So the tower does propagate at
N = infinity, with all of its weight diagonal.  See CHORD.md, and
diagonal_chord.py for the diagonal law.  Finite-N corrections are in
../syk-self-averaging/exact_pair_moments.py.

VERIFICATION (test_chord_charges.py): the recursion as a matrix identity;
C(m) = P_m(T) and [C(m), T] = 0 for m <= 8; the closed forms;
[Xhat_n, T] = 0 for n <= 5; and the word-level zeros |mu_k|/gross <= 1e-12
for n <= 5, even k <= 6, through the independent chord_moments assembly.

TRUNCATION.  On a finite chord lattice the identities hold exactly on the
interior block; truncation intrudes at most m rows from the top.
Residuals are measured there, relative to the block scale.
"""

from __future__ import annotations

import numpy as np

from chord import qint


# ---------------------------------------------------------------------------
# formal-basis builders (unscaled: T|n> = |n+1> + [n]|n-1>)
# ---------------------------------------------------------------------------
def transfer_formal(nmax: int, q: float) -> np.ndarray:
    T = np.zeros((nmax + 1, nmax + 1))
    for n in range(nmax + 1):
        if n + 1 <= nmax:
            T[n + 1, n] = 1.0
        if n >= 1:
            T[n - 1, n] = qint(n, q)
    return T


def transfer_one_matter_formal(nmax: int, q: float, dpsi: float) -> np.ndarray:
    """T_M on |n0, n1> (n0 below / n1 above the matter chord): open top 1,
    close from above [n1], close from below q^(D+n1)[n0]."""
    m = nmax + 1
    TM = np.zeros((m * m, m * m))
    for n0 in range(m):
        for n1 in range(m):
            i = n0 * m + n1
            if n1 + 1 <= nmax:
                TM[n0 * m + n1 + 1, i] += 1.0
            if n1 >= 1:
                TM[n0 * m + n1 - 1, i] += qint(n1, q)
            if n0 >= 1:
                TM[(n0 - 1) * m + n1, i] += q ** (dpsi + n1) * qint(n0, q)
    return TM


def open_map(nmax: int) -> np.ndarray:
    m = nmax + 1
    O = np.zeros((m * m, m))
    for n in range(m):
        O[n * m, n] = 1.0
    return O


def close_j_map(nmax: int, q: float, dpsi: float, j: int) -> np.ndarray:
    """close_j |n0,n1> = [n1]^(j)_falling q^(D(n1-j)) |n0+n1-j>.  j = 0 is
    the ordinary close."""
    m = nmax + 1
    C = np.zeros((m, m * m))
    for n0 in range(m):
        for n1 in range(m):
            tot = n0 + n1 - j
            if n1 < j or not 0 <= tot <= nmax:
                continue
            fall = 1.0
            for a in range(j):
                fall *= qint(n1 - a, q)
            C[tot, n0 * m + n1] = fall * q ** (dpsi * (n1 - j))
    return C


def dressed_propagator(m_steps: int, nmax: int, q: float,
                       dpsi: float) -> np.ndarray:
    """C(m) = close o T_M^m o open, built from the chord rules directly."""
    TM = transfer_one_matter_formal(nmax, q, dpsi)
    v = open_map(nmax)
    out = v.copy()
    for _ in range(m_steps):
        out = TM @ out
    return close_j_map(nmax, q, dpsi, 0) @ out


# ---------------------------------------------------------------------------
# the Motzkin polynomial of the lemma
# ---------------------------------------------------------------------------
def _padd(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if len(a) < len(b):
        a, b = b, a
    out = a.copy()
    out[:len(b)] += b
    return out


def motzkin_coeffs(m_steps: int, q: float, dpsi: float) -> np.ndarray:
    """Coefficients of P_m(T) (index = power of T) from the recursion:
    at height j -- flat: q^(D+j) T, down: [j]_q, up: 1 - q^(2D+j)."""
    polys = {0: np.array([1.0])}
    for _ in range(m_steps):
        new: dict[int, np.ndarray] = {}
        for j, pc in polys.items():
            flat = np.concatenate([[0.0], pc]) * q ** (dpsi + j)
            new[j] = _padd(new.get(j, np.zeros(1)), flat)
            if j >= 1:
                new[j - 1] = _padd(new.get(j - 1, np.zeros(1)),
                                   pc * qint(j, q))
            new[j + 1] = _padd(new.get(j + 1, np.zeros(1)),
                               pc * (1.0 - q ** (2 * dpsi + j)))
        polys = new
    return polys.get(0, np.array([0.0]))


def poly_of_matrix(coeffs: np.ndarray, T: np.ndarray) -> np.ndarray:
    out = np.zeros_like(T)
    P = np.eye(T.shape[0])
    for c in coeffs:
        out += c * P
        P = P @ T
    return out


def xhat_matrix(n: int, nmax: int, q: float, dpsi: float) -> np.ndarray:
    """The chord representation of the A_n insertion:
    Xhat_n = sum_r (-1)^r C(n,r) C(n-r) T^r."""
    from math import comb
    T = transfer_formal(nmax, q)
    out = np.zeros_like(T)
    Tp = np.eye(nmax + 1)
    for r in range(n + 1):
        out += (-1) ** r * comb(n, r) * (
            dressed_propagator(n - r, nmax, q, dpsi) @ Tp)
        Tp = Tp @ T
    return out


if __name__ == "__main__":
    q, dpsi, nmax = 0.5, 0.25, 20
    T = transfer_formal(nmax, q)
    keep = nmax - 8
    for m in (1, 2, 3, 4):
        Cm = dressed_propagator(m, nmax, q, dpsi)
        Pm = poly_of_matrix(motzkin_coeffs(m, q, dpsi), T)
        comm = Cm @ T - T @ Cm
        print(f"m={m}: ||C - P_m(T)|| = "
              f"{np.abs((Cm - Pm)[:keep, :keep]).max():.2e}   "
              f"||[C, T]|| = {np.abs(comm[:keep, :keep]).max():.2e}   "
              f"P_m coeffs = {np.round(motzkin_coeffs(m, q, dpsi), 6)}")
    X2 = xhat_matrix(2, nmax, q, dpsi)
    ref = (1 - q ** dpsi) ** 2 * (T @ T) + (1 - q ** (2 * dpsi)) * np.eye(nmax + 1)
    print(f"Xhat_2 closed form: dev = "
          f"{np.abs((X2 - ref)[:keep, :keep]).max():.2e}")
