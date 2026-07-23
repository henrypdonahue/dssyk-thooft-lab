#!/usr/bin/env python3
r"""
Singlet-moment pipeline: rung 15 groundwork.
============================================

The eventual N = infinity comparison ("moments vs 't Hooft sum rules",
road_map rung 15) needs the DSSYK side's spectral moments

    mu_k(n, CP) = int domega omega^k S_n^{CP}(omega)
                = Tr( ad_H^k(Mbar_n^{CP}) (Mbar_n^{CP})^dag ) / dim

of the CP-resolved, traceless tower operators.  This module measures them by
ED with disorder error bars along the FIXED-p trajectory (p = 4, Nc growing:
lambda = p^2/Nc falls toward the 't Hooft-matching corner), and reports the
dimensionless shape invariants that survive normalization ambiguities:

    w2   = mu_2 / (mu_0 * <H^2>/dim-per-mode-ish)   -- mean-square frequency
           in units of the instance's own energy scale sqrt(Tr H^2/dim),
    kurt = mu_4 mu_0 / mu_2^2                       -- shape (kurtosis-like).

Why moments and shape invariants: BOOST.md shows boost-frequency data cannot
be compared to bulk masses peak-by-peak (single masses are continua; the
hologram/bulk frequency map -- the tomperature identification -- is not yet
fixed).  Moment RATIOS of dimensionless shape parameters are the objects most
robust to that unresolved map, so this pipeline publishes exactly those, with
statistics, ready for the exact chord/Wick N = infinity computation (rung 15
proper) to be checked against.

Hermiticity subtlety (caught by this module's own test): for n >= 3 the raw
M_n is NOT +-Hermitian (dirac.py), so the ordered correlator Tr(M_n(t)M_n^dag)
has genuinely nonvanishing odd-in-omega parts -- at T_B = infinity the exact
relation is S_{AA^dag}(omega) = S_{A^dag A}(-omega), symmetric only for the
combined object.  We therefore use the SYMMETRIZED spectral function, i.e.
split each CP-projected operator into Hermitian + anti-Hermitian parts and add
their (individually omega-symmetric) weights; odd moments then vanish
identically (asserted in the tests).
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np

from dirac import (draw_couplings, dirac_hamiltonian, annihilators,
                   mn_tower, particle_hole, cp_projections,
                   symmetrized_weight)


@lru_cache(maxsize=8)
def _fixed_ops(Nc: int):
    """Disorder-independent operators, built once per Nc (not per instance)."""
    return particle_hole(Nc), annihilators(Nc)


def instance_moments(Nc: int, p: int, seed: int, n_list=(2, 3), n_max=None):
    """Moments mu_0, mu_1, mu_2, mu_4 of the CP-resolved traceless M_n spectral
    functions for ONE C-invariant disorder instance, plus the instance energy
    scale sqrt(Tr H^2/dim).  Returns {(n, 'right'/'wrong'): dict}."""
    if n_max is None:
        n_max = max(n_list)
    rng = np.random.default_rng(seed)
    H = dirac_hamiltonian(Nc, p,
                          couplings=draw_couplings(Nc, p, rng, "c_symmetric"))
    dim = 1 << Nc
    E, V = np.linalg.eigh(H)
    omega = E[:, None] - E[None, :]
    escale = float(np.sqrt(np.mean(E ** 2)))          # sqrt(Tr H^2/dim)
    U, cs = _fixed_ops(Nc)
    towers = mn_tower(H, cs, n_max)

    out = {}
    for n in n_list:
        Mb = towers[n] - np.trace(towers[n]) / dim * np.eye(dim)
        Mp, Mm = cp_projections(U, Mb)
        cp_target = (-1) ** (n + 1)
        chans = {"right": Mp if cp_target == +1 else Mm,
                 "wrong": Mm if cp_target == +1 else Mp}
        for tag, op in chans.items():
            W = symmetrized_weight(V.conj().T @ op @ V)
            mu = {k: float(np.sum(W * omega ** k)) for k in (0, 1, 2, 4)}
            out[(n, tag)] = dict(mu0=mu[0], mu1=mu[1], mu2=mu[2], mu4=mu[4],
                                 escale=escale)
    return out


def ensemble_moments(Nc: int, p: int, n_inst: int, seed0: int = 100,
                     n_list=(2, 3)):
    """Disorder mean +- standard error of the dimensionless shape invariants."""
    rows = {key: [] for n in n_list for key in [(n, "right"), (n, "wrong")]}
    for j in range(n_inst):
        inst = instance_moments(Nc, p, seed0 + j, n_list)
        for key, d in inst.items():
            w2 = d["mu2"] / (d["mu0"] * d["escale"] ** 2)
            kurt = d["mu4"] * d["mu0"] / d["mu2"] ** 2
            odd = abs(d["mu1"]) / np.sqrt(d["mu0"] * d["mu2"])
            rows[key].append((d["mu0"], w2, kurt, odd))
    result = []
    for (n, tag), vals in rows.items():
        arr = np.array(vals)
        mean = arr.mean(axis=0)
        sem = arr.std(axis=0, ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0 * mean
        result.append(dict(
            Nc=Nc, p=p, n=n, channel=tag, n_inst=n_inst,
            mu0=mean[0], mu0_err=float(sem[0]),
            w2=mean[1], w2_err=float(sem[1]),
            kurt=mean[2], kurt_err=float(sem[2]),
            odd_moment_check=float(mean[3])))
    return result


if __name__ == "__main__":
    p = 4
    all_rows = []
    print(f"CP-resolved M_n spectral moments, fixed p = {p} "
          f"(lambda = p^2/Nc falls toward the match corner)\n")
    print(f"{'Nc':>3} {'lam':>5} {'n':>2} {'chan':>6} "
          f"{'mu0':>10} {'w2 = mu2/(mu0<H^2>)':>20} {'kurt':>16}")
    for Nc, n_inst in [(6, 40), (8, 40), (10, 12)]:
        rows = ensemble_moments(Nc, p, n_inst)
        all_rows.extend(rows)
        for r in rows:
            print(f"{r['Nc']:>3} {p*p/r['Nc']:>5.2f} {r['n']:>2} "
                  f"{r['channel']:>6} {r['mu0']:>10.3e} "
                  f"{r['w2']:>12.4f} ± {r['w2_err']:<6.4f} "
                  f"{r['kurt']:>9.3f} ± {r['kurt_err']:<6.3f}")
    out = Path(__file__).resolve().parent / "moments.json"
    with open(out, "w") as fh:
        json.dump(all_rows, fh, indent=1)
    print(f"\nwrote {out}")
    print("\nThese dimensionless shape invariants (w2, kurt) per CP channel are")
    print("the map-robust targets for the exact N = infinity chord/Wick")
    print("computation of rung 15; the 't Hooft-side comparison additionally")
    print("needs the boost kernel + tomperature map (duality/BOOST.md).")
