#!/usr/bin/env python3
r"""
ED laboratory for the anti-scrambling sign test and the Harlow-Zhao fold.
=========================================================================

Ground truth for duality/antiscrambling.py (road_map rungs 25-26).  Finite-N
Majorana SYK has a spectrum bounded above AND below, so every Euclidean-FOLDED
correlator (Harlow-Zhao 2607.14215 Sec. 6) is finite and exactly computable --
their proposed mechanism for de Sitter anti-scrambling, posed by them as an
open test "in a concrete model such as the SYK model".  This module runs that
test at ED scale, plus the scrambling-sign measurement of the standard
(unfolded) OTOC that Cui-Kolchmeyer 2607.13665 make the dS discriminator.

Correlator conventions (locked to Streicher 1911.10171 via
duality/antiscrambling.py; asserted in duality/test_antiscrambling.py):

  * theta in (0, 2pi) is Euclidean angle; complex slot times
    z = -i beta theta/(2pi)  (Im theta = Lorentzian theta-time, growth rate
    v per unit theta-time); psi(z) = e^{izH} psi e^{-izH}.
  * UNCROSSED (TOC) configurations are path-ordered as written:
        F = (1/N^2) sum_ij Tr[rho psi_i(z1) psi_i(z2) psi_j(z3) psi_j(z4)].
  * The CROSSED configuration theta1 > theta3 > theta2 > theta4 is the
    path-ordered correlator: operator sequence by descending
    standard-configuration theta, WITH the fermionic transposition sign:
        F = -(1/N^2) sum_ij Tr[rho psi_i(z1) psi_j(z3) psi_i(z2) psi_j(z4)].
    (Convention fixed empirically: this choice tracks Streicher (3.10); the
    opposite sign is off by x25.)
  * At finite dimension the ordered-trace expression is an ENTIRE function
    of the four thetas, so the FOLD -- Harlow-Zhao's continuation (6.9),
    theta^QM = (-t1, 2pi n - t2, 2pi - t3, 2pi(n+1) - t4) -- is evaluated by
    freezing the operator sequence at its standard-configuration order and
    moving the arguments.  e^{+tau H} segments appear; finite precisely
    because the spectrum is bounded (HZ's requirement).

Measurements (disorder-averaged, bootstrap errors over instances):

  g2        -- thermal 2-pt vs the large-q closed form: conventions lock +
               the O(1/p, lambda) error floor that calibrates everything
               else (the floor is LARGE at ED sizes -- deviation ratios
               ~0.5 -- so 4-pt comparisons are sign/trend-level, and every
               fold verdict is read against the same-size unfolded floor).
  g2_fold   -- the 2-pt fold (HZ fig 7c): ED literal continuation to
               theta < 0 vs the continued closed form, which has a POLE at
               theta = -pi(1/v - 1): whether the finite-N object follows
               the continued semiclassical form toward its pole is the
               2-pt version of the "continuation commutes" question.
  toc/otoc  -- Euclidean 4-pt at uncrossed/crossed configurations vs
               Streicher (3.3)/(3.10).
  growth_*  -- Lorentzian OTOC at the symmetric crossed configuration:
               script-F/F_d vs Lorentzian theta-time thL; closed form has
               Re a = -1 growth (scrambling, CK's c > 0) at EVERY
               temperature.  The ED trend is the model-level verdict.
  fold_*    -- the SAME configuration through the HZ fold map (6.9, n=1),
               Euclidean and Lorentzian: ED vs continued closed form.
               Tracking at the calibrated floor supports HZ's
               continuation-commutes assumption in SYK; refusal to follow
               the continued form (e.g. toward its poles) is a Stokes-type
               obstruction.  Either way, first data on their question.

Output: fold_bench.json.  Runtime ~15-20 min default (the crossed evaluator
is O(N^2) batched dense matmuls); `--quick` runs a reduced campaign (~2 min).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from pauli_strings import majorana_hamiltonian_fast
from syk import majorana_operators


# ---------------------------------------------------------------------------
# complex-time correlator engine (energy eigenbasis)
# ---------------------------------------------------------------------------
_PSI_CACHE = {}


def _psis_dense(N: int):
    if N not in _PSI_CACHE:
        _PSI_CACHE[N] = [op.toarray() for op in majorana_operators(N)]
    return _PSI_CACHE[N]


def eig_instance(N: int, p: int, rng):
    """One disorder instance: (E, psis) with the Majoranas rotated to the H
    eigenbasis, stacked as an (N, dim, dim) array."""
    H = majorana_hamiltonian_fast(N, p, rng)
    E, V = np.linalg.eigh(H)
    psis = np.stack([V.conj().T @ P @ V for P in _psis_dense(N)])
    return E, psis


def thetas_to_z(thetas, beta: float):
    """theta (possibly complex; Im theta = Lorentzian theta-time) -> slot
    time z = -i beta theta/(2 pi): Euclidean tau = beta Re(theta)/(2 pi),
    real time t = beta Im(theta)/(2 pi)."""
    return [-1j * beta * complex(th) / (2.0 * np.pi) for th in thetas]


def bilocal2(E, psis, beta: float, z1: complex, z2: complex) -> complex:
    """G(z1, z2) = (1/N) sum_i Tr[rho psi_i(z1) psi_i(z2)], literal order:
    sum_ab rho_a e^{i(z1-z2)(Ea-Eb)} psi_ab psi_ba."""
    N = len(psis)
    rho = np.exp(-beta * E)
    phase = np.exp(1j * (z1 - z2) * np.subtract.outer(E, E))
    val = np.einsum('a,ab,iab,iba->', rho, phase, psis, psis)
    return complex(val / (N * rho.sum()))


def bilocal4_uncrossed(E, psis, beta: float, zs) -> complex:
    """F = (1/N^2) sum_ij Tr[rho psi_i(z1) psi_i(z2) psi_j(z3) psi_j(z4)],
    literal (path) order -- factorizes into two flavor-summed pair kernels."""
    z1, z2, z3, z4 = zs
    N = len(psis)
    rho = np.exp(-beta * E)
    dE = np.subtract.outer(E, E)

    def pair_sum(w):
        d = np.exp(1j * w * E)
        return np.einsum('iab,b,ibc->ac', psis, d, psis)

    A = pair_sum(z2 - z1)
    B = pair_sum(z4 - z3)
    front = rho * np.exp(1j * (z1 - z4) * E)
    mid = np.exp(1j * (z3 - z2) * E)
    val = np.einsum('a,ab,b,ba->', front, A, mid, B)
    return complex(val / (N * N * rho.sum()))


def bilocal4_crossed(E, psis, beta: float, zs) -> complex:
    """The path-ordered CROSSED correlator (sequence frozen at the standard
    configuration's descending-theta order, fermionic sign included):

        F = -(1/N^2) sum_ij Tr[rho psi_i(z1) psi_j(z3) psi_i(z2) psi_j(z4)].

    Evaluated at arbitrary complex z's (Lorentzian continuation and HZ
    folds).  O(N^2) batched matmuls."""
    z1, z2, z3, z4 = zs
    N, dim = psis.shape[0], psis.shape[1]
    rho = np.exp(-beta * E)
    dE = np.subtract.outer(E, E)
    P1 = psis * np.exp(1j * z1 * dE)          # psi_i(z1), batched
    P2 = psis * np.exp(1j * z2 * dE)
    P3 = psis * np.exp(1j * z3 * dE)
    P4 = psis * np.exp(1j * z4 * dE)
    tot = 0.0 + 0.0j
    for i in range(N):
        #  sum_j Tr[rho P1_i P3_j P2_i P4_j]
        M = P3 @ P2[i] @ P4                    # (N, dim, dim) batched
        tot += np.einsum('a,ab,jba->', rho, P1[i], M)
    return complex(-tot / (N * N * rho.sum()))


# ---------------------------------------------------------------------------
# the measurement campaign
# ---------------------------------------------------------------------------
CONFIG_CROSSED = (1.5 * np.pi, 0.5 * np.pi, np.pi, 0.0)  # theta1 > theta3 > theta2 > theta4
CONFIG_TOC = (1.9 * np.pi, 1.2 * np.pi, 0.9 * np.pi, 0.2 * np.pi)


def hz_fold(thetas, n: int = 1):
    """HZ (6.9) in theta units."""
    t1, t2, t3, t4 = thetas
    return (-t1, 2 * np.pi * n - t2, 2 * np.pi - t3, 2 * np.pi * (n + 1) - t4)


def connected_ratio(F, G12, G34, Nflav: int, n_boot: int = 200, seed: int = 0):
    """script-F/F_d = N (<F> - <G12><G34>) / (<G12><G34>), bootstrap SEM."""
    F, G12, G34 = (np.asarray(x) for x in (F, G12, G34))

    def est(idx):
        fd = G12[idx].mean() * G34[idx].mean()
        return Nflav * (F[idx].mean() - fd) / fd

    n = len(F)
    full = est(np.arange(n))
    rng = np.random.default_rng(seed)
    boots = np.array([est(rng.integers(0, n, n)) for _ in range(n_boot)])
    return (complex(full), float(np.std(boots.real, ddof=1)),
            float(np.std(boots.imag, ddof=1)))


def run_point(N: int, p: int, betaJ: float, n_inst: int, seed: int,
              thL_grid=(0.0, 0.4, 0.8, 1.2), fold_n: int = 1) -> dict:
    scriptJ = np.sqrt(2.0) * p
    beta = betaJ / scriptJ
    rng = np.random.default_rng(seed)

    th_g2 = [0.25 * np.pi, 0.5 * np.pi, np.pi, 1.5 * np.pi]
    th_g2_fold = [-0.3 * np.pi, -0.6 * np.pi, -np.pi]
    fold_base = hz_fold(CONFIG_CROSSED, fold_n)

    configs = {"toc": ("u", CONFIG_TOC),
               "otoc": ("c", CONFIG_CROSSED),
               "fold": ("c", fold_base)}
    for thL in thL_grid[1:]:
        configs[f"growth_thL={thL}"] = ("c", (
            CONFIG_CROSSED[0] + 1j * thL, CONFIG_CROSSED[1] + 1j * thL,
            CONFIG_CROSSED[2], CONFIG_CROSSED[3]))
        configs[f"fold_thL={thL}"] = ("c", (
            fold_base[0] + 1j * thL, fold_base[1] + 1j * thL,
            fold_base[2], fold_base[3]))

    g2_rows, g2f_rows = [], []
    samples = {k: ([], [], []) for k in configs}
    t0 = time.time()
    for _ in range(n_inst):
        E, psis = eig_instance(N, p, rng)
        g2_rows.append([bilocal2(E, psis, beta,
                                 *thetas_to_z((th, 0.0), beta)).real
                        for th in th_g2])
        g2f_rows.append([bilocal2(E, psis, beta,
                                  *thetas_to_z((th, 0.0), beta)).real
                         for th in th_g2_fold])
        for key, (kind, th) in configs.items():
            zs = thetas_to_z(th, beta)
            F = (bilocal4_uncrossed if kind == "u"
                 else bilocal4_crossed)(E, psis, beta, zs)
            samples[key][0].append(F)
            samples[key][1].append(bilocal2(E, psis, beta, zs[0], zs[1]))
            samples[key][2].append(bilocal2(E, psis, beta, zs[2], zs[3]))

    g2, g2f = np.array(g2_rows), np.array(g2f_rows)
    sem = 1.0 / np.sqrt(max(n_inst - 1, 1))
    out = dict(N=N, p=p, betaJ=betaJ, beta=beta, scriptJ=scriptJ,
               n_inst=n_inst, seed=seed, thL_grid=list(thL_grid),
               fold_n=fold_n, config_crossed=list(CONFIG_CROSSED),
               config_toc=list(CONFIG_TOC),
               fold_thetas=[float(x) for x in fold_base],
               th_g2=th_g2, th_g2_fold=th_g2_fold,
               g2_mean=g2.mean(0).tolist(),
               g2_err=(g2.std(0, ddof=1) * sem).tolist(),
               g2_fold_mean=g2f.mean(0).tolist(),
               g2_fold_err=(g2f.std(0, ddof=1) * sem).tolist(),
               runtime_s=round(time.time() - t0, 1))
    for key, (F, G12, G34) in samples.items():
        val, err_re, err_im = connected_ratio(F, G12, G34, N)
        out[key] = dict(re=val.real, im=val.imag,
                        err_re=err_re, err_im=err_im)
    return out


def main(quick: bool = False):
    if quick:
        points = [dict(N=12, p=6, betaJ=2.0, n_inst=24, seed=402)]
    else:
        points = [
            dict(N=12, p=6, betaJ=2.0, n_inst=96, seed=402),
            dict(N=14, p=6, betaJ=2.0, n_inst=48, seed=400),
            dict(N=16, p=6, betaJ=2.0, n_inst=12, seed=403),
            dict(N=14, p=6, betaJ=0.5, n_inst=48, seed=401),
        ]
    results = []
    for kw in points:
        print(f"running {kw} ...", flush=True)
        results.append(run_point(**kw))
        print(f"  done in {results[-1]['runtime_s']} s", flush=True)
    payload = dict(
        provenance="fold_bench.py -- ED ground truth for "
                   "duality/antiscrambling.py (rungs 25-26); conventions in "
                   "the module docstring; N-sequence at betaJ=2 probes the "
                   "continuation-commutes question vs system size",
        points=results)
    with open(Path(__file__).resolve().parent / "fold_bench.json", "w") as fh:
        json.dump(payload, fh, indent=1)
    print("wrote fold_bench.json")


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
