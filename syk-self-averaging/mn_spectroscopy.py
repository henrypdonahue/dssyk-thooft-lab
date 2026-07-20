#!/usr/bin/env python3
r"""
First single-realization spectroscopy of the M_n tower (Dirac SYK, T = inf).
============================================================================

The paper's "best test of the duality" needs the spectrum of the U(N)-singlet
tower M_n = sum_i c^dag_i (d^n/dt^n) c_i.  Nobody -- the authors included --
has published even an exploratory ED look at these correlators.  This script
is that first look: honest about what it can and cannot see.

What it computes, for ONE disorder realization (legitimate because the meson
spectral function is a singlet and the measured singlet self-averaging rate is
1/C(N,p) -- see self_averaging.py; at Nc = 10, p = 4 that is 2/C(10,4) ~ 1%):

    S_n(omega) = (1/dim) sum_{a,b} |<a| Mbar_n |b>|^2
                                     delta(omega - (E_a - E_b)),

the infinite-temperature spectral weight of the traceless tower operator
Mbar_n, plus its CP-RESOLVED versions S_n^+-(omega): on a C-invariant
instance (dirac.draw_couplings ensemble="c_symmetric") the operator is split
with cp_projections into definite-CP parts BEFORE taking matrix elements.
The raw M_n for n >= 2 are CP-mixtures (see dirac.py) -- the CP-resolved
spectral functions are the objects the paper's CP = (-1)^{n+1} assignment can
be tested on.

Exact structural checks built in:
  * n = 0, 1: all spectral weight sits at omega = 0 (conserved charge and
    Hamiltonian -- the non-propagating photon/graviton).  Machine-exact.
  * n >= 2: weight moves to omega != 0 (the dynamical candidates).

What this is NOT: the duality test.  At Nc = 10-12 and p = 4 the coupling is
lambda = p^2/Nc ~ 1.3-1.6, far from the lambda -> 0 corner where the 't Hooft
match is expected (correction #2), and Nc is far from large-N.  Expect broad
spectral distributions, not sharp meson peaks; the deliverables are the
pipeline, the conserved-mode checks, the CP-channel resolution, and a first
data point for how the weight organizes -- inputs for the moments-vs-sum-rules
program (road_map.md rung 15), not a spectrum match.

Output: mn_spectra.json (histogrammed S_n^+-, summary numbers) and
mn_spectra.png.  Runtime ~seconds at Nc = 10 (default); pass --big for
Nc = 12 (~minutes, dim 4096).
"""

from __future__ import annotations

import json
import sys
from math import comb

import numpy as np

from dirac import (draw_couplings, dirac_hamiltonian, annihilators,
                   charge_matrix, mn_operator, particle_hole, cp_projections)


def spectral_weights(M: np.ndarray, V: np.ndarray, E: np.ndarray):
    """All (omega_ab, |<a|M|b>|^2/dim) pairs for one operator."""
    Mt = V.conj().T @ M @ V
    W = np.abs(Mt) ** 2 / M.shape[0]
    omega = E[:, None] - E[None, :]
    return omega.ravel(), W.ravel()


def run(Nc: int = 10, p: int = 4, n_max: int = 6, seed: int = 21,
        n_bins: int = 121, outfile: str = None):
    rng = np.random.default_rng(seed)
    couplings = draw_couplings(Nc, p, rng, ensemble="c_symmetric")
    H = dirac_hamiltonian(Nc, p, couplings=couplings)
    dim = 1 << Nc
    U = particle_hole(Nc)
    cs = annihilators(Nc)
    lam = p * p / Nc

    print(f"Dirac SYK, Nc={Nc} (dim {dim}), p={p}, C-invariant instance, "
          f"lambda = p^2/Nc = {lam:.2f}")
    print("(single realization; singlet self-averaging rate 2/C(Nc,p) = "
          f"{2/comb(Nc, p):.1e})\n")

    E, V = np.linalg.eigh(H)

    results = dict(Nc=Nc, p=p, seed=seed, lam=lam, levels=n_max + 1, spectra=[])
    print(f"{'n':>3} {'CP tgt':>7} {'w(omega=0)':>11} {'right-CP frac':>14} "
          f"{'sqrt<w^2>':>10}")
    for n in range(n_max + 1):
        Mn = mn_operator(H, cs, n)
        Mb = Mn - np.trace(Mn) / dim * np.eye(dim)
        Mp, Mm = cp_projections(U, Mb)
        cp_target = (-1) ** (n + 1)
        right, wrong = (Mp, Mm) if cp_target == +1 else (Mm, Mp)

        row = dict(n=n, cp_target=cp_target)
        total_all, static_all = 0.0, 0.0
        for tag, op in [("right", right), ("wrong", wrong)]:
            omega, w = spectral_weights(op, V, E)
            total = w.sum()
            static_mask = np.abs(omega) < 1e-8
            static = w[static_mask].sum()
            # histogram the DYNAMICAL weight only; the omega = 0 (diagonal /
            # conserved) weight is reported separately as static_weight
            hist, edges = np.histogram(
                omega[~static_mask], bins=n_bins,
                weights=w[~static_mask],
                range=(-1.05 * np.max(np.abs(E)) * 2,
                       1.05 * np.max(np.abs(E)) * 2))
            row[tag] = dict(total_weight=float(total),
                            static_weight=float(static),
                            omega_second_moment=float((w * omega ** 2).sum()
                                                      / max(total, 1e-300)),
                            hist=hist.tolist(),
                            edges=edges.tolist())
            total_all += total
            static_all += static
        frac_right = row["right"]["total_weight"] / max(total_all, 1e-300)
        rms_omega = np.sqrt(
            (row["right"]["omega_second_moment"] * row["right"]["total_weight"]
             + row["wrong"]["omega_second_moment"] * row["wrong"]["total_weight"])
            / max(total_all, 1e-300))
        row["frac_right_cp"] = float(frac_right)
        row["static_fraction"] = float(static_all / max(total_all, 1e-300))
        results["spectra"].append(row)
        print(f"{n:>3} {cp_target:>+7d} {row['static_fraction']:>11.2e} "
              f"{frac_right:>14.3f} {rms_omega:>10.3f}")

    print("\nn = 0, 1: static_fraction = 1 (conserved photon/graviton) is the")
    print("built-in exactness check; n >= 2 weight is genuinely dynamical.")
    if outfile:
        with open(outfile, "w") as fh:
            json.dump(results, fh)
        print(f"wrote {outfile}")
    return results


def plot(results, fname="mn_spectra.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [r for r in results["spectra"] if r["n"] >= 2]
    fig, axes = plt.subplots(1, len(rows), figsize=(3.4 * len(rows), 3.4))
    if len(rows) == 1:
        axes = [axes]
    for ax, r in zip(axes, rows):
        # normalize each panel by its total dynamical weight: the shape is
        # the physics; absolute weights grow ~ omega^{2n} across panels
        dyn_total = sum(np.sum(r[tag]["hist"]) for tag in ("right", "wrong"))
        for tag, color, style in [("right", "#2166ac", "-"),
                                  ("wrong", "#b2182b", "--")]:
            edges = np.array(r[tag]["edges"])
            centers = 0.5 * (edges[1:] + edges[:-1])
            width = edges[1] - edges[0]
            ax.plot(centers,
                    np.array(r[tag]["hist"]) / (width * max(dyn_total, 1e-300)),
                    style, color=color, lw=1.2,
                    label=f"CP {'=' if tag == 'right' else '!='} "
                          f"{r['cp_target']:+d}")
        ax.set_title(f"$n={r['n']}$  (target CP = {r['cp_target']:+d})",
                     fontsize=10)
        ax.text(0.03, 0.95, f"static: {r['static_fraction']:.0%}",
                transform=ax.transAxes, fontsize=8, va="top", color="0.4")
        ax.set_xlabel(r"$\omega$")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel(r"$S_n(\omega)/\int S_n$  (dynamical part)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(
        f"CP-resolved $M_n$ spectral functions, one realization "
        f"(Nc={results['Nc']}, p={results['p']}, "
        f"$\\lambda$={results['lam']:.2f} -- NOT the $\\lambda\\to0$ match regime)",
        fontsize=10)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print(f"wrote {fname}")


if __name__ == "__main__":
    big = "--big" in sys.argv
    res = run(Nc=12 if big else 10, outfile="mn_spectra.json")
    plot(res)
