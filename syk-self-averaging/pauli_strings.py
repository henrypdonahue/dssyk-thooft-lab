#!/usr/bin/env python3
r"""
Fast Majorana-SYK Hamiltonian assembly, via symplectic Pauli strings.
=====================================================================

Every product of Jordan-Wigner Majoranas is a SINGLE Pauli string.  So
the Hamiltonian

    H = i^{p/2} sum_{i_1<...<i_p} J_{i_1...i_p} psi_{i_1} ... psi_{i_p}

is a sum of C(N,p) strings.  Each one acts on a computational-basis state
by one bit-flip pattern and one sign pattern.

Assembling H that way replaces the per-term sparse matmuls of
syk.majorana_hamiltonian, the reference implementation, which is kept for
cross-checking.  It is one to two orders of magnitude faster.  That is
what makes n_inst ~ 10^3 disorder ensembles, and N = 20-22 ED,
affordable.

REPRESENTATION.  A string is (phase, x, z), with phase in {0,1,2,3} and
x, z integer bitmasks over the nq = N/2 qubits.  It encodes

    P = i^{phase} * X^x Z^z,
    (X^x Z^z)|s> = (-1)^{popcount(z & s)} |s XOR x>

Z acts first, then X.  Multiplication is exact integer arithmetic:

    (i^{p1} X^{x1} Z^{z1}) (i^{p2} X^{x2} Z^{z2})
        = i^{p1 + p2 + 2*popcount(z1 & x2)}  X^{x1 XOR x2}  Z^{z1 XOR z2}.

BIT CONVENTION.  Qubit q of syk.majorana_operators, the leftmost kron
factor, corresponds to bit (nq-1-q) of the basis index.  That matches
scipy's kron ordering.

With Y = i X Z, the Jordan-Wigner Majoranas of syk.py are

    psi_{2q}   : phase 0, x = bit(q), z = bits(0..q-1)
    psi_{2q+1} : phase 1, x = bit(q), z = bits(0..q)

test_syk.py asserts both facts against the sparse reference.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import scipy.sparse as sp


def _qubit_bit(nq: int, q: int) -> int:
    """Bitmask of qubit q in the basis-index convention of syk.py's kron order."""
    return 1 << (nq - 1 - q)


def majorana_string(N: int, a: int):
    """(phase, x, z) for the single Majorana psi_a on nq = N/2 qubits."""
    nq = N // 2
    q = a // 2
    x = _qubit_bit(nq, q)
    z = 0
    for j in range(q):
        z |= _qubit_bit(nq, j)
    if a % 2 == 1:  # Y = i X Z carries phase i and a Z on its own qubit
        z |= _qubit_bit(nq, q)
        return 1, x, z
    return 0, x, z


def string_product(s1, s2):
    """Product of two (phase, x, z) strings.  Exact integer arithmetic."""
    p1, x1, z1 = s1
    p2, x2, z2 = s2
    phase = (p1 + p2 + 2 * bin(z1 & x2).count("1")) % 4  # py3.9: no int.bit_count
    return phase, x1 ^ x2, z1 ^ z2


def monomial_string(N: int, combo) -> tuple:
    """(phase, x, z) for the ordered product psi_{c1} psi_{c2} ... psi_{cp}."""
    s = (0, 0, 0)
    for a in combo:
        s = string_product(s, majorana_string(N, a))
    return s


_I_POW = np.array([1, 1j, -1, -1j])


def add_string_dense(H: np.ndarray, s, coef: complex, states: np.ndarray) -> None:
    """Accumulate  coef * i^{phase} X^x Z^z  into the dense matrix H in place."""
    phase, x, z = s
    signs = 1.0 - 2.0 * (np.bitwise_count(states & z) & 1)
    H[states ^ x, states] += (coef * _I_POW[phase]) * signs


# Nothing in the string table depends on the disorder draw: the (x, z, phase)
# of each monomial is fixed by (N, p), and the scatter row vector states^x and
# the Z-parity sign vector depend only on x and z.  A disorder ensemble
# therefore builds all of it once.  The (x, z, amp) list is always cached (it
# is a few thousand small ints).  The dim-length row/sign vectors are cached
# only below a dimension cap, and deduped on x and z -- there are at most `dim`
# distinct values of each, so the tables stay O(dim^2).
_STRING_TABLE: dict = {}
_ROW_CACHE: dict = {}
_SIGN_CACHE: dict = {}
_VECTOR_CACHE_MAX_DIM = 1 << 10


def _string_table(N: int, p: int, prefac_phase: int):
    """[(x, z, i^{phase + prefac_phase})] in itertools.combinations order."""
    key = (N, p, prefac_phase)
    hit = _STRING_TABLE.get(key)
    if hit is None:
        hit = []
        for combo in combinations(range(N), p):
            phase, x, z = monomial_string(N, combo)
            hit.append((x, z, _I_POW[(phase + prefac_phase) % 4]))
        _STRING_TABLE[key] = hit
    return hit


def _rows_signs(dim: int, x: int, z: int, states: np.ndarray):
    if dim > _VECTOR_CACHE_MAX_DIM:
        return states ^ x, 1.0 - 2.0 * (np.bitwise_count(states & z) & 1)
    rows = _ROW_CACHE.get((dim, x))
    if rows is None:
        rows = states ^ x
        _ROW_CACHE[(dim, x)] = rows
    signs = _SIGN_CACHE.get((dim, z))
    if signs is None:
        signs = 1.0 - 2.0 * (np.bitwise_count(states & z) & 1)
        _SIGN_CACHE[(dim, z)] = signs
    return rows, signs


def majorana_hamiltonian_fast(N: int, p: int, rng, dense: bool = True):
    """One disorder realization of the Majorana SYK Hamiltonian, assembled by
    Pauli-string bitwise arithmetic.

    Draws couplings in the same itertools.combinations order and from the same
    Gaussian as syk.majorana_hamiltonian, so identical (N, p, rng-state) give
    the identical operator (asserted in test_syk.py).  Returns a dense complex
    ndarray by default (what eigvalsh wants).  With dense=False it returns CSR.
    """
    from syk import coupling_variance
    if p % 2 != 0:
        raise ValueError("p must be even for a Hermitian Majorana Hamiltonian")
    if N % 2 != 0:
        raise ValueError("N must be even")
    dim = 1 << (N // 2)
    sigma = np.sqrt(coupling_variance(N, p))
    prefac_phase = (p // 2) % 4  # i^{p/2} folded into the string phase
    states = np.arange(dim, dtype=np.int64)
    H = np.zeros((dim, dim), dtype=complex)
    # one += per monomial, in combinations order: identical accumulation order
    # (and identical rng consumption) to the uncached route, so H is unchanged
    # to the last bit
    for (x, z, amp) in _string_table(N, p, prefac_phase):
        J = rng.normal(0.0, sigma)
        rows, signs = _rows_signs(dim, x, z, states)
        H[rows, states] += (J * amp) * signs
    if dense:
        return H
    return sp.csr_matrix(H)


if __name__ == "__main__":
    import time
    from syk import majorana_hamiltonian

    for (N, p) in [(10, 4), (12, 4), (12, 6)]:
        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(42)
        t0 = time.time()
        H_ref = majorana_hamiltonian(N, p, rng1).toarray()
        t_ref = time.time() - t0
        t0 = time.time()
        H_fast = majorana_hamiltonian_fast(N, p, rng2)
        t_fast = time.time() - t0
        err = np.abs(H_ref - H_fast).max()
        print(f"N={N:2d} p={p}:  max|H_ref - H_fast| = {err:.2e}   "
              f"reference {t_ref*1e3:7.1f} ms   fast {t_fast*1e3:6.1f} ms   "
              f"speedup x{t_ref/t_fast:5.1f}")
