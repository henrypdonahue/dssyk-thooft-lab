#!/usr/bin/env python3
r"""
Sparse SYK Hamiltonians (Majorana and complex), built by Jordan-Wigner.
======================================================================

This is the reusable core for the self-averaging tests in this module, following
the conventions of Miyashita-Sekino-Susskind (arXiv:2607.05678), Section 3.

Majorana model (their Eq. 3.2):

    H = i^{p/2}  sum_{i_1<...<i_p}  J_{i_1...i_p}  psi_{i_1} ... psi_{i_p}

with p even, N Majorana fermions psi_a (a = 0..N-1), and Gaussian couplings of
variance (their Eq. 3.8, with the overall scale J = 1)

    Var(J) = p! / N^{p-1}.

Conventions
-----------
We use Hermitian Majoranas normalized to  psi_a^2 = I,  {psi_a, psi_b} = 2 delta_ab,
so the normalized (infinite-temperature) trace of the identity is 1 and
Tr(psi_a^2)/dim = 1.  N must be even; the Hilbert space has dimension 2^{N/2},
built on N/2 qubits by the standard Jordan-Wigner strings

    psi_{2q}   = Z ... Z (q of them)  X  I ... I
    psi_{2q+1} = Z ... Z (q of them)  Y  I ... I .

The factor i^{p/2} makes H Hermitian for even p (checked in the tests).

Complex/Dirac model (their Section 3, "Dirac case"): fermions c_i = (psi_{2i} +
i psi_{2i+1})/2 carry a conserved U(1) charge Q = sum_i c_i^dagger c_i, the seed
of the paper's electromagnetic sector.  Provided here mainly so the singlet
charge operator M0 = Q (their Eq. 4.24) is available.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

# Single-qubit Pauli matrices as sparse CSR.
_I2 = sp.identity(2, format="csr", dtype=complex)
_X = sp.csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
_Y = sp.csr_matrix(np.array([[0, -1j], [1j, 0]], dtype=complex))
_Z = sp.csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))


def _kron_list(ops):
    """Sparse Kronecker product of a list of 2x2 operators (left = qubit 0)."""
    out = ops[0]
    for op in ops[1:]:
        out = sp.kron(out, op, format="csr")
    return out


def majorana_operators(N: int):
    """Return the list [psi_0, ..., psi_{N-1}] as sparse CSR matrices.

    Hermitian, psi_a^2 = I, mutually anticommuting.  N must be even.
    """
    if N % 2 != 0:
        raise ValueError("N (number of Majoranas) must be even")
    nq = N // 2
    psis = []
    for a in range(N):
        q = a // 2
        pauli = _X if (a % 2 == 0) else _Y
        ops = [_Z] * q + [pauli] + [_I2] * (nq - q - 1)
        psis.append(_kron_list(ops))
    return psis


def coupling_variance(N: int, p: int) -> float:
    """Var(J) = p! / N^{p-1}   (Eq. 3.8, overall scale set to 1)."""
    from math import factorial
    return factorial(p) / N ** (p - 1)


def majorana_hamiltonian(N: int, p: int, rng, psis=None):
    """Build one disorder realization of the Majorana SYK Hamiltonian.

    Parameters
    ----------
    N, p : ints (p even, p <= N)
    rng  : numpy Generator (the disorder draw)
    psis : optional precomputed majorana_operators(N) to avoid rebuilding.

    Returns the sparse Hermitian H (dim 2^{N/2}).
    """
    from itertools import combinations
    if p % 2 != 0:
        raise ValueError("p must be even for a Hermitian Majorana Hamiltonian")
    if psis is None:
        psis = majorana_operators(N)
    dim = 1 << (N // 2)
    sigma = np.sqrt(coupling_variance(N, p))
    prefac = (1j) ** (p // 2)

    H = sp.csr_matrix((dim, dim), dtype=complex)
    for combo in combinations(range(N), p):
        J = rng.normal(0.0, sigma)
        term = psis[combo[0]]
        for a in combo[1:]:
            term = term @ psis[a]
        H = H + (prefac * J) * term
    return H.tocsr()


def fermion_parity(N: int):
    """Total fermion parity  P = i^{N/2} psi_0 psi_1 ... psi_{N-1}  (Hermitian,
    P^2 = I).  H commutes with P, so the spectrum splits into two sectors;
    resolving them is essential for the level-statistics test (see validate_syk)."""
    psis = majorana_operators(N)
    P = psis[0]
    for a in range(1, N):
        P = P @ psis[a]
    P = ((1j) ** (N // 2)) * P
    return P.tocsr()


# --------------------------------------------------------------------------
# Complex / Dirac sector (optional; for the charged U(1) operators)
# --------------------------------------------------------------------------
def dirac_operators(N_complex: int):
    """Return (c, c_dag) lists for N_complex complex fermions, i.e. 2*N_complex
    Majoranas.  c_i = (psi_{2i} + i psi_{2i+1})/2 obeys {c_i, c_j^dag} = delta_ij."""
    psis = majorana_operators(2 * N_complex)
    c, cdag = [], []
    for i in range(N_complex):
        ci = 0.5 * (psis[2 * i] + 1j * psis[2 * i + 1])
        c.append(ci.tocsr())
        cdag.append(ci.getH().tocsr())
    return c, cdag


def charge_operator(N_complex: int):
    """Singlet charge  Q = sum_i c_i^dag c_i  (the paper's M0, Eq. 4.24)."""
    c, cdag = dirac_operators(N_complex)
    dim = 1 << N_complex
    Q = sp.csr_matrix((dim, dim), dtype=complex)
    for i in range(N_complex):
        Q = Q + cdag[i] @ c[i]
    return Q.tocsr()
