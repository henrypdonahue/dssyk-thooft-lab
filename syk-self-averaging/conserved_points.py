#!/usr/bin/env python3
r"""
The (N, p) = (8, 4) per-instance conservation, explained.
=========================================================

mn_majorana_bench.py records that at (N, p) = (8, 4) every singlet
channel A_n = sum_i psi_i (ad_H)^n psi_i commutes with H per instance,
that the same is false at (12, 6) and (16, 8), and that the mechanism
was unknown.  It is a short Clifford-algebra fact.  This module states
it, and asserts every step at machine precision (test_conserved_points).

TWO EXACT IDENTITIES, ANY (N, p).  Write S(X) = sum_i psi_i X psi_i, the
spin-trace map.  On a Clifford monomial psi_J of grade g = |J| each
psi_i either lies in J (psi_i psi_J psi_i = (-1)^{g-1} psi_J) or not
((-1)^g psi_J), so

    S(psi_J) = (-1)^g (N - 2g) psi_J .                                (1)

S is diagonal on grades and annihilates grade N/2.  The nested
commutator expands binomially into

    A_n = sum_{j=0}^{n} (-1)^j C(n, j) S(H^{n-j}) H^j ,  S(H^0) = N .  (2)

Two faces of (2) are worth having by name.  With H_i = (H - psi_i H
psi_i)/2 the part of H that contains flavor i (sum_i H_i = pH),

    A_1 = -2 sum_i H_i = -2pH,        A_2 = 4 sum_i H_i^2 .          (3)

The first is the bench's exact anchor.  The second is new: the A_2
channel is, per instance, the sum of the squared single-flavor energies,
so its width is the width of sum_i H_i^2 -- the finite-N face of the
diagonal chord law, whose H_i-insertions are exactly what carries the
(1 - 2p/N) -> q^Delta weight in duality/diagonal_chord.py.

THE GENERAL IDENTITY.  A_3 is not a multiple of sum_i H_i^3 (asserted),
but the pattern does not stop: it continues with MIXED words in H and
the H_i.  For a set T of slots let W_T = sum_i prod_{t=1..n} (H_i if t
in T else H).  Then, for every n >= 2, every N and every even p,

    A_n = sum_{|T| >= 2} (-2)^|T| (-1)^(n - max T) C(n-1, max T - 1) W_T. (3')

n = 2 is (3).  n = 3 is A_3 = 4 (W_13 + W_23 - 2 W_12) - 8 W_123, i.e.
4 sum_i (H_i H H_i + H H_i^2 - 2 H_i^2 H) - 8 sum_i H_i^3.  It follows
from (1) and (2) by inclusion-exclusion: a product of m Hamiltonian
monomials has grade |c_1 xor ... xor c_m| = sum_T (-2)^(|T|-1)
|cap_T c_t|, the |T| = 1 terms cancel in the binomial sum for n >= 2,
and |cap_T c_t| = sum_i prod_{t in T} [i in c_t] is what turns H into
H_i on the slots of T.  This is the operator behind the N = infinity
tower law of duality/singlet_tower.py, where only the pairs |T| = 2
survive the double-scaled limit.

WHY (8, 4) IS SPECIAL.  For even p the Hermitian H is symmetric on each
chirality block of Cl_N^even, and so is every H^k.  At N = 8 the even
algebra is M_8(R) + M_8(R); grade-2 and grade-6 monomials act as
antisymmetric 8 x 8 matrices on each block (28 + 28 = 2 x dim so(8)),
while grades 0, 4, 8 act symmetrically (1 + 70 + 1 = 2 x 36).  Hence
every H^k has Clifford content in grades {0, 4, 8} ONLY.  By (1), S
kills grade 4 and sends grades 0 and 8 to multiples of 1 and
Gamma = psi_1 ... psi_8, both central in the even algebra.  So every
S(H^k) is central, every A_n lies in span{H^j, Gamma H^j} -- the
commutant of H -- and [H, A_n] = 0 for every n and every instance.
Explicitly, with (H^2)_g the grade-g part,

    A_2 = 8 H^2 + 8 (H^2)_0 - 8 (H^2)_8 .                            (4)

At (12, 6) the non-central grades of H^k are 4 and 8, with S-weights +4
and -4, so S(H^2) = 24 (H^2)_0 + 4 (H^2)_4 - 4 (H^2)_8 is not a
polynomial in H and the identity fails -- measured at 0.3 in the bench.
For p = N/2 the non-central content of H^k sits in the single grade N/2
only when N = 8, so (8, 4) is the unique such point, as the review found
empirically.

TWO CORNERS THE SAME ALGEBRA COVERS.  p = 2 at any N is free: A_n is the
bilinear with single-particle matrix K^n, so [H, A_n] = 0 for all n.  And
(8, 6) is (8, 2) in disguise, psi_I = +-Gamma psi_{I^c} for |I| = 6.  Both
are asserted.  Neither is in the bench grid.

Nothing here touches the N = infinity theorem (duality/chord_charges.py)
or the ladder laws; it closes the one "mechanism unknown" line in the
bench and removes the (8, 4) row from the class of unexplained data.
"""

from __future__ import annotations

import itertools
from math import comb

import numpy as np

from pauli_strings import majorana_hamiltonian_fast
from syk import majorana_operators


def spin_trace_weight(N: int, g: int) -> int:
    """Eq. (1): S(psi_J) = (-1)^g (N - 2g) psi_J for |J| = g."""
    return (-1) ** g * (N - 2 * g)


def spin_trace(X: np.ndarray, psis: list) -> np.ndarray:
    """S(X) = sum_i psi_i X psi_i."""
    return sum(P @ X @ P for P in psis)


def dense_majoranas(N: int) -> list:
    return [op.toarray().astype(complex) for op in majorana_operators(N)]


def chirality(psis: list) -> np.ndarray:
    """Gamma = psi_1 ... psi_N (phase-free product; commutes with the even
    algebra, squares to +-1)."""
    G = np.eye(psis[0].shape[0], dtype=complex)
    for P in psis:
        G = G @ P
    return G


def flavor_parts(H: np.ndarray, psis: list) -> list:
    """H_i = (H - psi_i H psi_i)/2: the terms of H that contain flavor i."""
    return [(H - P @ H @ P) / 2 for P in psis]


def grade_content(X: np.ndarray, psis: list) -> dict:
    """Frobenius norm of the grade-g Clifford component of X, per grade
    (monomials are orthonormal under Tr(M^dag M')/dim)."""
    N = len(psis)
    dim = X.shape[0]
    out = {}
    for g in range(N + 1):
        tot = 0.0
        for J in itertools.combinations(range(N), g):
            M = np.eye(dim, dtype=complex)
            for j in J:
                M = M @ psis[j]
            tot += abs(np.trace(M.conj().T @ X) / dim) ** 2
        out[g] = float(np.sqrt(tot))
    return out


def identity_coefficients(n: int) -> dict:
    """The coefficients of Eq. (3'), keyed by the slot set T (a sorted
    tuple of 1-based slots, |T| >= 2)."""
    return {T: (-2) ** len(T) * (-1) ** (n - max(T)) * comb(n - 1, max(T) - 1)
            for r in range(2, n + 1)
            for T in itertools.combinations(range(1, n + 1), r)}


def mixed_word(H: np.ndarray, Hi: list, n: int, T: tuple) -> np.ndarray:
    """W_T = sum_i prod_{t=1..n} (H_i if t in T else H), dense."""
    out = np.zeros_like(H)
    for h in Hi:
        M = np.eye(H.shape[0], dtype=complex)
        for t in range(1, n + 1):
            M = M @ (h if t in T else H)
        out += M
    return out


def general_identity(H: np.ndarray, psis: list, n: int) -> np.ndarray:
    """The right-hand side of Eq. (3'), evaluated literally."""
    Hi = flavor_parts(H, psis)
    out = np.zeros_like(H)
    for T, coef in identity_coefficients(n).items():
        out += coef * mixed_word(H, Hi, n, T)
    return out


def a_n_tower(H: np.ndarray, psis: list, n_max: int) -> dict:
    """A_n = sum_i psi_i (ad_H)^n psi_i, n = 1..n_max (dense reference)."""
    A = {n: np.zeros_like(H) for n in range(1, n_max + 1)}
    for P in psis:
        D = P
        for n in range(1, n_max + 1):
            D = H @ D - D @ H
            A[n] += P @ D
    return A


def binomial_assembly(H: np.ndarray, psis: list, n: int) -> np.ndarray:
    """Eq. (2), evaluated literally."""
    dim = H.shape[0]
    N = len(psis)
    powers = [np.eye(dim, dtype=complex)]
    for _ in range(n):
        powers.append(powers[-1] @ H)
    out = np.zeros_like(H)
    for j in range(n + 1):
        S = N * powers[0] if n - j == 0 else spin_trace(powers[n - j], psis)
        out += (-1) ** j * comb(n, j) * (S @ powers[j])
    return out


def commutant_residual(A: np.ndarray, H: np.ndarray, G: np.ndarray) -> float:
    """Relative residual of A against span{H^j, Gamma H^j}, j < dim/2 (the
    commutant of an H with simple spectrum per chirality block).  Meant
    for dim <= 16; the power basis is too ill-conditioned beyond that."""
    dim = H.shape[0]
    cols = []
    Hk = np.eye(dim, dtype=complex)
    for _ in range(dim // 2):
        cols.append(Hk.ravel())
        cols.append((G @ Hk).ravel())
        Hk = Hk @ H
    B = np.stack(cols, 1)
    c = np.linalg.lstsq(B, A.ravel(), rcond=None)[0]
    return float(np.linalg.norm(B @ c - A.ravel()) / np.linalg.norm(A))


def commutator_residual(A: np.ndarray, H: np.ndarray) -> float:
    return float(np.linalg.norm(H @ A - A @ H) / np.linalg.norm(A))


def report(N: int, p: int, seed: int = 0, n_max: int = 5) -> dict:
    rng = np.random.default_rng(seed)
    H = majorana_hamiltonian_fast(N, p, rng)
    psis = dense_majoranas(N)
    G = chirality(psis)
    A = a_n_tower(H, psis, n_max)
    Hi = flavor_parts(H, psis)
    rows = {}
    for n in range(1, n_max + 1):
        rows[n] = dict(commutator=commutator_residual(A[n], H),
                       binomial=float(np.linalg.norm(
                           binomial_assembly(H, psis, n) - A[n])
                           / np.linalg.norm(A[n])))
        if N <= 8:
            rows[n]["commutant"] = commutant_residual(A[n], H, G)
    out = dict(N=N, p=p,
               S_of_H=float(np.linalg.norm(spin_trace(H, psis))
                            / np.linalg.norm(H)),
               A2_vs_4sumHi2=float(np.linalg.norm(
                   A[2] - 4 * sum(h @ h for h in Hi)) / np.linalg.norm(A[2])),
               general_identity={n: float(np.linalg.norm(
                   general_identity(H, psis, n) - A[n]) / np.linalg.norm(A[n]))
                   for n in range(2, n_max + 1)},
               tower=rows)
    if N <= 8:
        out["grades_H2"] = grade_content(H @ H, psis)
    return out


if __name__ == "__main__":
    for (N, p) in [(8, 4), (8, 2), (8, 6), (10, 4), (12, 6)]:
        r = report(N, p)
        print(f"(N, p) = ({N}, {p}):  |S(H)|/|H| = {r['S_of_H']:.1e}   "
              f"|A_2 - 4 sum H_i^2|/|A_2| = {r['A2_vs_4sumHi2']:.1e}   "
              f"general identity n = 2..5: "
              f"{max(r['general_identity'].values()):.1e}")
        for n, row in r["tower"].items():
            extra = (f"  in commutant: {row['commutant']:.1e}"
                     if "commutant" in row else "")
            print(f"   n={n}: |[H,A_n]|/|A_n| = {row['commutator']:.1e}"
                  f"  binomial identity: {row['binomial']:.1e}{extra}")
        if "grades_H2" in r:
            live = {g: round(v, 3) for g, v in r["grades_H2"].items()
                    if v > 1e-9}
            print(f"   grade content of H^2: {live}")
