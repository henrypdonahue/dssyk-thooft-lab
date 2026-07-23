#!/usr/bin/env python3
r"""
Validate the SYK construction against the textbook signature: random-matrix
level statistics that cycle through GOE / GUE / GSE with N mod 8 (the Bott
periodicity of the Majorana SYK ensemble).

Diagnostic: the mean adjacent-gap ratio  <r>,  r_n = min(s_n,s_{n-1})/max(...),
s_n = E_{n+1}-E_n.  It needs no spectral unfolding and has known universal
values:

    Poisson  0.3863     GOE  0.5307     GUE  0.5996     GSE  0.6744

Crucial subtlety (this bit trips people up): H commutes with fermion parity, so
its spectrum is a superposition of two independent blocks.  Mixed, the level
repulsion is washed out and <r> fakes the Poisson value ~0.39.  You must resolve
parity first.  In our Jordan-Wigner convention parity is diagonal --
P = (-1)^{N/2} Z^{⊗N/2} -- so a sector is just the computational-basis states of
fixed qubit-parity (even/odd popcount).

Predicted class by N mod 8 (within a parity sector):

    N mod 8 = 0 -> GOE      = 2 -> GUE      = 4 -> GSE      = 6 -> GUE
"""

from __future__ import annotations

import numpy as np

from pauli_strings import majorana_hamiltonian_fast

RMT = {"Poisson": 0.3863, "GOE": 0.5307, "GUE": 0.5996, "GSE": 0.6744}
CLASS_BY_NMOD8 = {0: "GOE", 2: "GUE", 4: "GSE", 6: "GUE"}


def _even_parity_mask(N: int) -> np.ndarray:
    """Boolean mask over the 2^{N/2} computational basis states selecting the
    even-qubit-parity (P-eigenvalue) sector."""
    nq = N // 2
    idx = np.arange(1 << nq)
    return (np.bitwise_count(idx) % 2) == 0


def mean_r(eigs: np.ndarray, trim: float = 0.15) -> float:
    """Mean adjacent-gap ratio over the central (1-2*trim) of a sorted spectrum.

    The s > 1e-12 filter below is load-bearing: in the GSE cases (N mod 8 = 4)
    every level is Kramers-doubled within a parity sector, and without dropping
    the zero gaps the r-statistic would be garbage.  It also makes the GOE/GUE
    cases robust to accidental near-degeneracies."""
    e = np.sort(eigs)
    lo, hi = int(len(e) * trim), int(len(e) * (1 - trim))
    e = e[lo:hi]
    s = np.diff(e)
    s = s[s > 1e-12]
    r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return float(np.mean(r))


def r_statistic(N: int, p: int = 4, n_real: int = 8, seed: int = 0) -> float:
    """<r> in the even-parity sector, averaged over n_real disorder realizations.

    Only valid against CLASS_BY_NMOD8 for p = 0 mod 4: the antiunitary
    classification of Majorana SYK depends on p mod 4 (p = 2 mod 4 realizes a
    different N mod 8 table), so guard rather than silently validate against
    the wrong ensemble."""
    if p % 4 != 0:
        raise ValueError(
            f"p={p}: CLASS_BY_NMOD8 is the p = 0 (mod 4) classification; "
            "p = 2 (mod 4) has a different N mod 8 table")
    mask = _even_parity_mask(N)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_real):
        H = majorana_hamiltonian_fast(N, p, rng)
        Hs = H[np.ix_(mask, mask)]
        eigs = np.linalg.eigvalsh(Hs)
        vals.append(mean_r(eigs))
    return float(np.mean(vals))


def _nearest_class(r: float) -> str:
    return min(RMT, key=lambda k: abs(RMT[k] - r))


def run(ns=(10, 12, 14, 16, 18, 20), p: int = 4):
    print("Level-statistics validation (even-parity sector, p = %d)\n" % p)
    print(f"{'N':>3} {'N%8':>4} {'<r>':>8} {'predicted':>10} {'value':>7} {'match':>6}")
    print("-" * 46)
    all_ok = True
    for N in ns:
        # more realizations for the smaller (noisier) sectors
        n_real = 40 if N <= 14 else (16 if N <= 18 else 8)
        r = r_statistic(N, p=p, n_real=n_real, seed=100 + N)
        pred = CLASS_BY_NMOD8[N % 8]
        got = _nearest_class(r)
        ok = (got == pred)
        all_ok = all_ok and ok
        print(f"{N:>3} {N % 8:>4} {r:>8.4f} {pred:>10} {RMT[pred]:>7.4f} "
              f"{'OK' if ok else 'X':>6}")
    print("-" * 46)
    print("Bott periodicity " + ("reproduced." if all_ok else "NOT fully reproduced."))
    return all_ok


if __name__ == "__main__":
    run()
