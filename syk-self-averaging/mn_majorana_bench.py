#!/usr/bin/env python3
r"""
Majorana singlet-channel spectral moments: the ED side of rung 15 proper.
=========================================================================

duality/chord_moments.py computes the EXACT N = infinity spectral moments of
the Majorana singlet bilinear channels

    A_n = sum_i psi_i (ad_H)^n psi_i ,      n = 1, 2, 3

from the chord expansion (two matter chords of dimension Delta_psi = 1/p).
This module measures the SAME channels by exact diagonalization at finite N
along fixed p = 4 (lambda_chord = 2p^2/N falling), disorder-averaged with
errors, so the chord prediction can be checked as a genuine N -> infinity
extrapolation with no Dirac-vs-Majorana channel caveat (moments_pipeline.py
measures the Dirac CP-resolved tower instead; both JSONs are kept).

Conventions (single point of truth = dirac.symmetrized_weight, mirrored
here): moments of the SYMMETRIZED T = infinity spectral function,

    mu_k = [Tr(ad_H^k(A_h) A_h) + Tr(ad_H^k(A_a) A_a^dag)]/dim
         = (1/2)[Tr(ad^k(A) A^dag) + Tr(ad^k(A^dag) A)]/dim ,

with A traceless (mu_0 subtracts the identity component; mu_{k>=1} is
insensitive to it because the trace of a commutator vanishes).  Shape
invariants reported: w2 = mu_2/(mu_0 <H^2>/dim) and kurt = mu_4 mu_0/mu_2^2.

Exact anchor asserted in the tests: A_1 = -2p H identically (psi gamma psi
summed over flavors = (N-2p) gamma for a p-Majorana string), so the n = 1
channel is conserved: mu_2 = mu_4 = 0 at EVERY N.  The chord side reproduces
this identity only up to its finite-p matter-weight error -- the measured
size of that residual is part of the rung-15 honesty budget.

Output: mn_majorana_bench.json.  Runtime ~2-4 min.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from pauli_strings import majorana_hamiltonian_fast
from syk import majorana_operators


def a_n_tower(H: np.ndarray, psis: list, n_max: int) -> dict:
    """A_n = sum_i psi_i (ad_H)^n psi_i, built incrementally: D^0 = psi_i,
    D^{k+1} = [H, D^k]."""
    dim = H.shape[0]
    A = {n: np.zeros((dim, dim), dtype=complex) for n in range(1, n_max + 1)}
    for P in psis:
        D = P.astype(complex)
        for n in range(1, n_max + 1):
            D = H @ D - D @ H
            A[n] += P @ D
    return A


def instance_moments(N: int, p: int, seed: int, n_list=(1, 2, 3)) -> dict:
    rng = np.random.default_rng(seed)
    H = majorana_hamiltonian_fast(N, p, rng)
    dim = H.shape[0]
    E, V = np.linalg.eigh(H)
    omega = E[:, None] - E[None, :]
    escale2 = float(np.mean(E ** 2))
    psis = [op.toarray() for op in majorana_operators(N)]
    A = a_n_tower(H, psis, max(n_list))
    out = {}
    for n in n_list:
        M = V.conj().T @ A[n] @ V
        M -= np.trace(M) / dim * np.eye(dim)
        Mh = 0.5 * (M + M.conj().T)
        Ma = 0.5 * (M - M.conj().T)
        W = (np.abs(Mh) ** 2 + np.abs(Ma) ** 2) / dim
        mu = {k: float(np.sum(W * omega ** k)) for k in (0, 1, 2, 4)}
        out[n] = dict(mu0=mu[0], mu1=mu[1], mu2=mu[2], mu4=mu[4],
                      escale2=escale2)
    return out


def ensemble(N: int, p: int, n_inst: int, seed0: int = 7000,
             n_list=(1, 2, 3)) -> list:
    rows = {n: [] for n in n_list}
    for j in range(n_inst):
        inst = instance_moments(N, p, seed0 + j, n_list)
        for n, d in inst.items():
            if n == 1 or abs(d["mu2"]) < 1e-9 * abs(d["mu0"]) * d["escale2"]:
                # conserved channel (n = 1 always; n >= 2 exactly at the
                # p = N/2 artifact point): kurt is 0/0, record zeros
                rows[n].append((d["mu0"], 0.0, 0.0,
                                d["mu2"], d["mu4"], d["escale2"]))
            else:
                w2 = d["mu2"] / (d["mu0"] * d["escale2"])
                kurt = d["mu4"] * d["mu0"] / d["mu2"] ** 2
                rows[n].append((d["mu0"], w2, kurt,
                                d["mu2"], d["mu4"], d["escale2"]))
    out = []
    for n, vals in rows.items():
        arr = np.array(vals)
        mean = arr.mean(axis=0)
        sem = (arr.std(axis=0, ddof=1) / np.sqrt(len(arr))
               if len(arr) > 1 else 0 * mean)
        out.append(dict(
            N=N, p=p, n=n, n_inst=n_inst,
            lambda_chord=2.0 * p * p / N,
            mu0=mean[0], mu0_err=float(sem[0]),
            w2=mean[1], w2_err=float(sem[1]),
            kurt=mean[2], kurt_err=float(sem[2]),
            mu2=mean[3], mu2_err=float(sem[3]),
            mu4=mean[4], mu4_err=float(sem[4]),
            escale2=mean[5]))
    return out


def main():
    t0 = time.time()
    rows = []
    # fixed p = 4: the 't Hooft-corner trajectory (lambda = 32/N falls)
    # plus the small-N anchors, and a fixed p = 6 series -- the pairs
    # (p=4, N=8) <-> (p=6, N=18) sit at MATCHED lambda_chord = 4.0, the
    # clean test of the chord fixed-lambda N -> infinity prediction
    # NOTE the (p=4, N=8) row sits exactly ON the p = N/2 artifact
    # (E[B_jj] ~ (1-2p/N) = 0, exact_wick.py): its n >= 2 moments vanish
    # identically -- kept as the exact exhibition of that artifact, and
    # excluded from trend fits.  (6, 22) provides the near-matched-lambda
    # anchor to (4, 10): lambda = 3.27 vs 3.20.
    for p, N, n_inst in [(4, 8, 256), (4, 10, 160), (4, 12, 96),
                         (4, 14, 64), (4, 16, 48), (4, 18, 24), (4, 20, 12),
                         (6, 14, 64), (6, 16, 48), (6, 18, 24), (6, 20, 12),
                         (6, 22, 6)]:
        print(f"p={p} N={N} ({n_inst} instances) ...", flush=True)
        rows.extend(ensemble(N, p, n_inst))
    payload = dict(
        provenance="mn_majorana_bench.py -- ED side of rung 15 proper: "
                   "Majorana singlet channels A_n = sum_i psi_i ad_H^n "
                   "psi_i, symmetrized traceless moments, fixed p = 4, "
                   "lambda_chord = 32/N; chord comparison in "
                   "duality/chord_moments.py",
        runtime_s=round(time.time() - t0, 1), points=rows)
    out = Path(__file__).resolve().parent / "mn_majorana_bench.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"wrote {out} ({payload['runtime_s']} s)")


if __name__ == "__main__":
    main()
