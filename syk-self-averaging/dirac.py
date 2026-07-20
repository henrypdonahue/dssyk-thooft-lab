#!/usr/bin/env python3
r"""
Charge-conserving complex (Dirac) SYK: the U(1)/"QED" sector of the duality.
============================================================================

The U(N) side of Miyashita-Sekino-Susskind (arXiv:2607.05678) is the Dirac
SYK model -- complex fermions c_i, i = 0..Nc-1, with a charge-conserving
Hamiltonian containing p/2 creation and p/2 annihilation operators per term:

    H = sum_{|A|=|B|=p/2}  J_{AB}  c^dag_{a1} ... c^dag_{ak} c_{b1} ... c_{bk},

    (A, B ascending index sets, k = p/2;  J_{BA} = conj(J_{AB}) so H = H^dag;
     diagonal couplings J_{AA} real.)

This is the sector where the paper's electromagnetic story lives: M_0 = Q is
the conserved U(1) charge ("photon"), M_1 the conserved "graviton", and the
tower  M_n = sum_i c^dag_i (d^n/dt^n) c_i  the meson candidates.  The module
provides the Hamiltonian (fast bitwise builder + independent sparse reference),
the charge, the M_n tower, and the (unitary) particle-hole/CP machinery --
everything the operator-identity and CP tests in test_syk.py exercise.

Operator identities encoded here (and asserted in the tests):

  *  [H, Q] = 0                       (charge conservation; the U(1) exists)
  *  M_0 = Q                          (paper Sec. 4.3, exact)
  *  M_1 = -(p/2) i H                 (NOT the paper's M_1 = H; forced by
                                       symmetry: Qdot = 0 makes M_1
                                       anti-Hermitian, so it cannot equal a
                                       nonzero Hermitian H.  README correction
                                       #3.)
  *  M_2 = - sum_i cdot^dag_i cdot_i  (because M_1 is conserved,
                                       d/dt M_1 = sum ci^dag cddot_i
                                       + sum cdot^dag cdot_i = 0; the paper's
                                       psi^dag psiddot = +psidot^dag psidot
                                       has the sign wrong.)
  *  M_n^dag = (-1)^n M_n  holds for n = 0, 1, 2 and FAILS for n >= 3
                                      (measured; the "n-th derivative of one
                                       Hermitian bilinear" picture survives
                                       whole only through n = 2).

CP structure (the paper's CP = (-1)^{n+1}, Sec. 4.3), as actually measured:

  *  The relevant operation is the UNITARY particle-hole map C -- see
     particle_hole() for why it cannot be antiunitary (M_1 ~ iH).
  *  On a C-invariant instance (draw_couplings ensemble="c_symmetric"),
     CP(M_0) = -1 and CP(M_1) = +1 hold EXACTLY.
  *  For n >= 2 the raw M_n operators are CP-MIXTURES (wrong-channel
     Frobenius weight 15-60%; measure with wrong_channel_weight).  This is
     structural, not numerical: M_n = d^n/dt2^n O(t1,t2)|_{t1=t2}
     differentiates only ONE time argument, so total-time-derivative
     (commutator) pieces carrying the opposite CP contaminate the operator.
     The sharp version of the paper's assignment must therefore be stated at
     the level of CP-RESOLVED spectral functions (split with cp_projections,
     then take two-point functions) -- mn_spectroscopy.py measures exactly
     that.

Conventions: occupation-number basis on Nc qubits, mode i <-> bit (Nc-1-i) of
the basis index (matching the kron order of syk.dirac_operators, against which
the fast builder is tested).  d/dt X = i[H, X].
"""

from __future__ import annotations

from itertools import combinations
from math import factorial

import numpy as np


# --------------------------------------------------------------------------
# couplings
# --------------------------------------------------------------------------
def coupling_variance_dirac(Nc: int, p: int) -> float:
    """Var(|J_AB|) chosen so <Tr H^2>/dim is extensive (~Nc), the analog of the
    Majorana convention Var = p!/N^{p-1}.  Each pair (A,B) contributes
    Var * 2^{-|A u B|} to <m2> and there are ~C(Nc,p/2)^2 ~ (Nc^{p/2}/(p/2)!)^2
    pairs, so Var = ((p/2)!)^2 2^p / Nc^{p-1} gives <m2> ~ Nc up to 1/Nc
    corrections.  The operator-identity tests are normalization-independent."""
    k = p // 2
    return factorial(k) ** 2 * 2 ** p / Nc ** (p - 1)


def draw_couplings(Nc: int, p: int, rng, ensemble: str = "generic") -> dict:
    """Hermitian Gaussian coupling table {(A, B): J_AB} over ascending p/2-index
    sets, with J_BA = conj(J_AB) and real diagonal.

    ensemble = "generic":
        all (A, B) pairs; complex off-diagonal entries with E|J|^2 = Var
        (split evenly re/im).  Generic instances break charge conjugation.

    ensemble = "c_symmetric":
        DISJOINT pairs (A ^ B = empty) with REAL symmetric couplings.  This is
        the maximal sub-ensemble that is exactly invariant under the unitary
        particle-hole map C (see particle_hole): C(c^dag_A c_B) = c^dag_B c_A
        with no leftover contraction terms only when A and B share no mode, and
        C-invariance + Hermiticity then force J real symmetric (for p/2 even;
        p/2 odd would need imaginary antisymmetric J, not implemented).
        Restricting to disjoint pairs discards an O(p^2/Nc) fraction of
        couplings -- immaterial at large Nc, and the price of an instance on
        which CP quantum numbers are sharp."""
    if p % 2 != 0:
        raise ValueError("p must be even (p/2 creators + p/2 annihilators)")
    k = p // 2
    if k > Nc:
        raise ValueError("p/2 exceeds the number of modes")
    sigma = np.sqrt(coupling_variance_dirac(Nc, p))
    subsets = list(combinations(range(Nc), k))
    J = {}
    if ensemble == "generic":
        for ia, A in enumerate(subsets):
            for B in subsets[ia:]:
                if A == B:
                    J[(A, B)] = rng.normal(0.0, sigma)
                else:
                    z = (rng.normal(0.0, sigma / np.sqrt(2))
                         + 1j * rng.normal(0.0, sigma / np.sqrt(2)))
                    J[(A, B)] = z
                    J[(B, A)] = np.conj(z)
    elif ensemble == "c_symmetric":
        if k % 2 != 0:
            raise ValueError("c_symmetric ensemble needs p/2 even (real J); "
                             "p/2 odd would require imaginary antisymmetric J")
        for ia, A in enumerate(subsets):
            for B in subsets[ia:]:
                if set(A) & set(B):
                    continue
                x = rng.normal(0.0, sigma)
                J[(A, B)] = x
                if A != B:
                    J[(B, A)] = x
    else:
        raise ValueError(f"unknown ensemble {ensemble!r}")
    return J


# --------------------------------------------------------------------------
# fast bitwise builder
# --------------------------------------------------------------------------
def _mode_bit(Nc: int, i: int) -> int:
    """Bitmask of mode i (kron ordering: mode 0 = most significant bit)."""
    return 1 << (Nc - 1 - i)


def _below_mask(Nc: int, i: int) -> int:
    """Bitmask of all modes j < i (the Jordan-Wigner string of mode i)."""
    m = 0
    for j in range(i):
        m |= _mode_bit(Nc, j)
    return m


def apply_monomial(Nc: int, A, B, states: np.ndarray):
    """Action of  c^dag_{A} c_{B}  (both ascending) on every basis state.

    Returns (valid, rows, signs): boolean mask of states not annihilated, the
    image state indices, and the fermionic signs.  Ascending monomials are
    applied right-to-left (largest annihilator first, largest creator first);
    with that order every Jordan-Wigner parity count can use the *pre-image*
    occupation for c_B and the post-annihilation occupation for c^dag_A, since
    already-cleared/set bits lie above the mode being counted."""
    Bmask = 0
    for b in B:
        Bmask |= _mode_bit(Nc, b)
    Amask = 0
    for a in A:
        Amask |= _mode_bit(Nc, a)

    occupied_B = (states & Bmask) == Bmask
    mid = states & ~Bmask
    free_A = (mid & Amask) == 0
    valid = occupied_B & free_A
    rows = (mid | Amask)

    par = np.zeros_like(states)
    for b in B:  # counts on the original occupation
        par += np.bitwise_count(states & _below_mask(Nc, b))
    for a in A:  # counts on the post-annihilation occupation
        par += np.bitwise_count(mid & _below_mask(Nc, a))
    signs = 1.0 - 2.0 * (par & 1)
    return valid, rows, signs


def dirac_hamiltonian(Nc: int, p: int, rng=None, couplings: dict = None):
    """One disorder realization of the charge-conserving complex SYK H as a
    dense (2^Nc x 2^Nc) Hermitian matrix.  Pass `couplings` (from
    draw_couplings) to reuse a table; otherwise `rng` draws a fresh one."""
    if couplings is None:
        if rng is None:
            raise ValueError("need rng or couplings")
        couplings = draw_couplings(Nc, p, rng)
    dim = 1 << Nc
    states = np.arange(dim, dtype=np.int64)
    H = np.zeros((dim, dim), dtype=complex)
    for (A, B), J in couplings.items():
        valid, rows, signs = apply_monomial(Nc, A, B, states)
        H[rows[valid], states[valid]] += J * signs[valid]
    return H


# --------------------------------------------------------------------------
# independent (slow) reference builder, for the equivalence test
# --------------------------------------------------------------------------
def dirac_hamiltonian_reference(Nc: int, p: int, couplings: dict):
    """Same H assembled from the sparse Jordan-Wigner matrices of
    syk.dirac_operators -- shares no assembly code with the fast builder."""
    from syk import dirac_operators
    c, cdag = dirac_operators(Nc)
    dim = 1 << Nc
    H = np.zeros((dim, dim), dtype=complex)
    for (A, B), J in couplings.items():
        term = np.eye(dim, dtype=complex)
        for a in A:
            term = term @ cdag[a].toarray()
        for b in B:
            term = term @ c[b].toarray()
        H += J * term
    return H


# --------------------------------------------------------------------------
# charge, tower, particle-hole / CP
# --------------------------------------------------------------------------
def charge_diagonal(Nc: int) -> np.ndarray:
    """Diagonal of Q = sum_i c^dag_i c_i in the occupation basis (= popcount)."""
    states = np.arange(1 << Nc, dtype=np.int64)
    return np.bitwise_count(states).astype(float)


def charge_matrix(Nc: int) -> np.ndarray:
    return np.diag(charge_diagonal(Nc)).astype(complex)


def annihilators(Nc: int) -> list:
    """Dense c_i matrices (from the bitwise action; O(dim) each)."""
    dim = 1 << Nc
    states = np.arange(dim, dtype=np.int64)
    out = []
    for i in range(Nc):
        valid, rows, signs = apply_monomial(Nc, (), (i,), states)
        C = np.zeros((dim, dim), dtype=complex)
        C[rows[valid], states[valid]] = signs[valid]
        out.append(C)
    return out


def time_derivative(H: np.ndarray, X: np.ndarray) -> np.ndarray:
    """d/dt X = i [H, X]."""
    return 1j * (H @ X - X @ H)


def mn_operator(H: np.ndarray, cs: list, n: int) -> np.ndarray:
    """M_n = sum_i c^dag_i (d^n/dt^n) c_i."""
    M = np.zeros_like(H)
    for C in cs:
        Cn = C
        for _ in range(n):
            Cn = time_derivative(H, Cn)
        M += C.conj().T @ Cn
    return M


def particle_hole(Nc: int) -> np.ndarray:
    """The UNITARY charge-conjugation (particle-hole) operator

        C = (c_0 + c^dag_0)(c_1 + c^dag_1) ... (c_{Nc-1} + c^dag_{Nc-1}),

    with C c_i C^dag = c^dag_i (all signs +1 in this Jordan-Wigner convention;
    asserted in the tests) and C (Q - Nc/2) C^dag = -(Q - Nc/2).

    Why the paper's "CP" must be this UNITARY map and not an antiunitary one:
    M_1 = -(p/2) i H, and any antiunitary operation that fixes H flips the
    explicit i, giving CP(M_1) = -1 -- but the paper assigns the n=1 "graviton"
    CP = +1.  Unitary C fixes i, so C M_1 C^dag = +M_1 exactly on a
    C-invariant instance.  (In the 1+1d bulk this is charge conjugation
    composed with spatial parity -- a unitary symmetry, as CP always is; the
    hologram is 0+1d, so its image here is just C.)"""
    cs = annihilators(Nc)
    U = np.eye(1 << Nc, dtype=complex)
    for C in cs:
        U = U @ (C + C.conj().T)
    return U


def cp_conjugate(U: np.ndarray, X: np.ndarray) -> np.ndarray:
    """C X C^dag (unitary conjugation)."""
    return U @ X @ U.conj().T


def cp_projections(U: np.ndarray, X: np.ndarray):
    """Split X into its CP-even and CP-odd parts, X = X_+ + X_-,
    with C X_pm C^dag = pm X_pm."""
    XC = cp_conjugate(U, X)
    return 0.5 * (X + XC), 0.5 * (X - XC)


def wrong_channel_weight(U: np.ndarray, X: np.ndarray, cp: int) -> float:
    """Fraction of the Frobenius weight of X in the OPPOSITE CP channel to
    `cp` (+1 or -1).  0 means X has the pure quantum number cp; 0.5 means an
    even mixture.  The traceless part of X should be passed in."""
    Xp, Xm = cp_projections(U, X)
    right = Xp if cp == +1 else Xm
    wrong = Xm if cp == +1 else Xp
    nr, nw = np.linalg.norm(right) ** 2, np.linalg.norm(wrong) ** 2
    return nw / (nr + nw) if nr + nw > 0 else 0.0


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    Nc, p = 5, 4
    couplings = draw_couplings(Nc, p, rng)
    H = dirac_hamiltonian(Nc, p, couplings=couplings)
    H_ref = dirac_hamiltonian_reference(Nc, p, couplings)
    Q = charge_matrix(Nc)
    cs = annihilators(Nc)

    print(f"Dirac SYK, Nc={Nc} (dim {1<<Nc}), p={p}")
    print(f"  builder vs sparse reference : {np.abs(H - H_ref).max():.2e}")
    print(f"  Hermiticity |H - H^dag|     : {np.abs(H - H.conj().T).max():.2e}")
    print(f"  charge conservation |[H,Q]| : {np.abs(H @ Q - Q @ H).max():.2e}")

    M0 = mn_operator(H, cs, 0)
    M1 = mn_operator(H, cs, 1)
    M2 = mn_operator(H, cs, 2)
    print(f"  M_0 - Q                     : {np.abs(M0 - Q).max():.2e}")
    print(f"  M_1 + (p/2) i H             : {np.abs(M1 + (p/2)*1j*H).max():.2e}")
    print(f"  M_1 - H  (paper's claim)    : {np.abs(M1 - H).max():.2e}   <- fails at O(1)")
    cdots = [time_derivative(H, C) for C in cs]
    sum_cdd = sum(Cd.conj().T @ Cd for Cd in cdots)
    print(f"  M_2 + sum cdot^dag cdot     : {np.abs(M2 + sum_cdd).max():.2e}")

    print("\nHermiticity parity M_n^dag = (-1)^n M_n (exact for n <= 2 only):")
    for n in range(5):
        Mn = mn_operator(H, cs, n)
        dev = np.abs(Mn.conj().T - (-1) ** n * Mn).max() / max(np.abs(Mn).max(), 1e-300)
        print(f"  n={n}: relative deviation = {dev:.2e}")

    print("\nCP (UNITARY particle-hole C; C-invariant disjoint/real ensemble):")
    U = particle_hole(Nc)
    dim = 1 << Nc
    print(f"  unitarity |U U^dag - 1|     : {np.abs(U @ U.conj().T - np.eye(dim)).max():.2e}")
    Qbar = Q - (Nc / 2) * np.eye(dim)
    print(f"  C Qbar C^dag + Qbar         : {np.abs(cp_conjugate(U, Qbar) + Qbar).max():.2e}"
          f"   (charge is CP-odd, exact)")
    Hc = dirac_hamiltonian(Nc, p,
                           couplings=draw_couplings(Nc, p, rng, "c_symmetric"))
    print(f"  |C Hc C^dag - Hc|           : {np.abs(cp_conjugate(U, Hc) - Hc).max():.2e}"
          f"   (C-invariant instance)")
    print(f"  |M1[Hc] + (p/2) i Hc|       : "
          f"{np.abs(mn_operator(Hc, cs, 1) + (p/2)*1j*Hc).max():.2e}")
    print("\n  wrong-CP-channel Frobenius weight of traceless M_n "
          "(0 = pure CP=(-1)^(n+1)):")
    for n in range(6):
        Mn = mn_operator(Hc, cs, n)
        Mb = Mn - np.trace(Mn) / dim * np.eye(dim)
        w = wrong_channel_weight(U, Mb, (-1) ** (n + 1))
        tag = "exact" if w < 1e-12 else "CP-mixture -> use CP-resolved correlators"
        print(f"    n={n}: {w:.4f}   ({tag})")
