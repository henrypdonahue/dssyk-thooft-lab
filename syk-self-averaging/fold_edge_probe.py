#!/usr/bin/env python3
r"""
The fold's soft-edge mechanism, measured (rung-26 epilogue).
============================================================

Decider #2 (duality/chord_fold.py) found that at matched lambda the ED
folded connected fraction sits ~1.5-1.6x above the exact N = infinity
value while unfolded observables agree to ~10% -- the fold amplifies
finite-N corrections.  Working hypothesis stated there: the folded
e^{+tau H} weights are dominated by the SPECTRAL MAXIMUM, whose finite-N
fluctuations the sharp-edged chord theory lacks.

This module measures that mechanism directly: per-instance folded and
unfolded connected fractions at the HZ (6.9) n = 1 configuration vs the
instance's spectral extremes.  Findings (N = 12, p = 4, betaJ = 2,
48 instances; fold_edge_probe.json, asserted in
duality/test_chord_fold.py::test_soft_edge_mechanism):

  * corr(E_max, folded fraction) ~ -0.86 -- the folded observable is a
    near-deterministic function of the instance's UPPER spectral edge;
  * the unfolded control shows no such dependence (corr ~ +0.15), and its
    instance-to-instance scatter is ~35x smaller;
  * the correlation with the LOWER edge is weak (~ -0.23): the (6.9) fold
    weights the upper edge specifically, as e^{+tau H} predicts.

Consequence, stated honestly: folded correlators at ED sizes measure
edge fluctuations (a Tracy-Widom-scale finite-N effect), not the bulk
double-scaled dynamics -- which is WHY their 1/N convergence to the chord
answer is slow and non-uniform, completing the decider-#2 story.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from fold_bench import (CONFIG_CROSSED, bilocal2, bilocal4_crossed,
                        eig_instance, hz_fold, thetas_to_z)


def run(N: int = 12, p: int = 4, betaJ: float = 2.0, n_inst: int = 48,
        seed: int = 777) -> dict:
    scriptJ = np.sqrt(2.0) * p
    beta = betaJ / scriptJ
    rng = np.random.default_rng(seed)
    fold_th = hz_fold(CONFIG_CROSSED, 1)
    rows = []
    for _ in range(n_inst):
        E, psis = eig_instance(N, p, rng)
        zf = thetas_to_z(fold_th, beta)
        zu = thetas_to_z(CONFIG_CROSSED, beta)
        Ff = bilocal4_crossed(E, psis, beta, zf).real
        Fu = bilocal4_crossed(E, psis, beta, zu).real
        Gf = (bilocal2(E, psis, beta, zf[0], zf[1])
              * bilocal2(E, psis, beta, zf[2], zf[3])).real
        Gu = (bilocal2(E, psis, beta, zu[0], zu[1])
              * bilocal2(E, psis, beta, zu[2], zu[3])).real
        rows.append((float(E.max()), float(E.min()),
                     float((Ff - Gf) / Gf), float((Fu - Gu) / Gu)))
    arr = np.array(rows)
    Emax, Emin, ff, fu = arr.T

    def corr(a, b):
        return float(np.corrcoef(a, b)[0, 1])

    return dict(N=N, p=p, betaJ=betaJ, n_inst=n_inst, seed=seed,
                corr_Emax_folded=corr(Emax, ff),
                corr_Emax_unfolded=corr(Emax, fu),
                corr_negEmin_folded=corr(-Emin, ff),
                folded_frac_mean=float(ff.mean()),
                folded_frac_std=float(ff.std(ddof=1)),
                unfolded_frac_mean=float(fu.mean()),
                unfolded_frac_std=float(fu.std(ddof=1)))


def main():
    t0 = time.time()
    out = run()
    out["provenance"] = ("fold_edge_probe.py -- the soft-edge mechanism "
                         "behind the folded 4-pt's finite-N behavior "
                         "(rung-26 epilogue; duality/ANTISCRAMBLING.md)")
    out["runtime_s"] = round(time.time() - t0, 1)
    path = Path(__file__).resolve().parent / "fold_edge_probe.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"corr(E_max, folded) = {out['corr_Emax_folded']:+.3f}  "
          f"(unfolded control {out['corr_Emax_unfolded']:+.3f}); "
          f"scatter ratio = "
          f"{out['folded_frac_std'] / out['unfolded_frac_std']:.0f}x")
    print(f"wrote {path} ({out['runtime_s']} s)")


if __name__ == "__main__":
    main()
