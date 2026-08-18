#!/usr/bin/env python3
r"""
The chord-side Euclidean fold at N = infinity (rung-26 decider #2).
==================================================================

Decider #1 (`syk-self-averaging/fold_bench.py --scan`, fold_scan.json)
established at ED sizes that the Harlow-Zhao-folded crossed 4-pt drifts away
from the continued large-q closed form with growing N at fixed lambda, and
that the folded connected fraction reaches O(1).  The final word was
explicitly deferred to the exact double-scaled theory: the fold is a contour
statement, chord correlators are analytic in the chord times, and the chord
engine (chord.py, validated) evaluates ANY complex contour.  This module
runs that computation -- the last open fold question.

Everything mirrors fold_bench.py exactly (same theta configurations, same
HZ (6.9) fold map n = 1, same frozen path-order convention, same connected
ratio script-F/F_d = N (F - F_d)/F_d with N = 2p^2/lambda_chord), with the
chord replacements:

  * Tr[e^{-beta H} psi(z_1) psi(z_3) psi(z_2) psi(z_4)] -> the exact
    two-matter-chord contour amplitude (Delta_psi = 1/p, fermionic
    matter-matter crossing sign = fold_bench's derived transposition minus);
  * the units bridge: beta_chord = beta_phys sigma_H.  Primary rows use the
    N = infinity form betaJ e^{-lambda/8}/sqrt(lambda); ED-comparison
    columns use the exact finite-N sigma_H(N, p) of the matched scan point.

Objects per lambda: toc / otoc / fold (Euclidean), growth_thL / fold_thL
(Lorentzian theta-time on the (theta_1, theta_2) pair), the folded and
unfolded connected fractions, the continued Streicher closed forms
(duality/antiscrambling.py, rung-25 verified transcription), and the ED
columns from fold_scan.json at matched lambda.  Truncation is doubled and
the drift stored.

Output: chord_fold.json; verdicts asserted in test_chord_fold.py and
narrated in ANTISCRAMBLING.md.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

from antiscrambling import fold_map, otoc_ratio, toc_ratio
from chord import Contour, Sector0, beta_chord_from_betaJ, sigma_H
from largeq_anchor import v_of_betaJ

HERE = Path(__file__).resolve().parent

CONFIG_CROSSED = (1.5 * np.pi, 0.5 * np.pi, np.pi, 0.0)
CONFIG_TOC = (1.9 * np.pi, 1.2 * np.pi, 0.9 * np.pi, 0.2 * np.pi)


def thetas_to_z(thetas, beta: float):
    """theta -> slot time z = -i beta theta/(2 pi)  (fold_bench convention:
    Im theta = Lorentzian theta-time)."""
    return [-1j * beta * complex(th) / (2.0 * np.pi) for th in thetas]


def _segments(zs_path, beta: float):
    """u_after values (trace order) for operators at path-ordered complex
    slot times zs_path: u_1 = beta + i(z_last - z_1) wraps; u_j =
    i(z_{j-1} - z_j)."""
    k = len(zs_path)
    us = [beta + 1j * (zs_path[-1] - zs_path[0])]
    for j in range(1, k):
        us.append(1j * (zs_path[j - 1] - zs_path[j]))
    return us


class ChordCorrelators:
    """2- and 4-pt thermal correlators at fixed (q, Delta_psi, beta_chord)."""

    def __init__(self, lam: float, p: int, beta_c: float, nmax: int):
        self.lam, self.p, self.beta_c, self.nmax = lam, p, beta_c, nmax
        self.q = math.exp(-lam)
        self.dpsi = 1.0 / p
        s = Sector0(nmax, self.q)
        self.Z = float(s.V[0] @ (np.exp(-beta_c * s.E) * s.V[0]))
        # one Contour per instance: its Sector0 eigendecomposition and
        # cached SectorM matrices are reused across every g2/f4 call
        self._c = Contour(self.q, self.nmax,
                          {'A': self.dpsi, 'B': self.dpsi},
                          fermionic={'A', 'B'})

    def g2(self, z1: complex, z2: complex) -> complex:
        us = _segments([z1, z2], self.beta_c)
        return self._c.correlator([('A', us[0]), ('A', us[1])]) / self.Z

    def f4(self, labels, zs_path) -> complex:
        """Path-ordered 4-pt: labels in path order (e.g. A B A B crossed),
        zs_path the matching slot times."""
        us = _segments(list(zs_path), self.beta_c)
        ops = [(lab, u) for lab, u in zip(labels, us)]
        return self._c.correlator(ops) / self.Z

    def f4_physical(self, thetas, kind: str) -> complex:
        """fold_bench's F at a theta configuration.  kind 'c' freezes the
        crossed path order (z1, z3, z2, z4) and carries fold_bench's
        derived transposition minus: the chord amplitude is the plain
        ordered trace <Tr[psi psi psi psi]> (whose matter-matter crossing
        sign makes it ~ -F_d), while fold_bench defines
        F = -(1/N^2) sum Tr[...] ~ +F_d -- so F = -amplitude.  kind 'u' is
        the uncrossed order (z1, z2, z3, z4), no extra sign."""
        zs = thetas_to_z(thetas, self.beta_c)
        if kind == 'c':
            return -self.f4(['A', 'B', 'A', 'B'],
                            [zs[0], zs[2], zs[1], zs[3]])
        return self.f4(['A', 'A', 'B', 'B'], zs)

    def connected_ratio(self, thetas, kind: str, Nflav: float) -> complex:
        """script-F/F_d = Nflav (F - G12 G34)/(G12 G34), fold_bench
        conventions."""
        zs = thetas_to_z(thetas, self.beta_c)
        F = self.f4_physical(thetas, kind)
        G12 = self.g2(zs[0], zs[1])
        G34 = self.g2(zs[2], zs[3])
        Fd = G12 * G34
        return Nflav * (F - Fd) / Fd


def run_point(lam: float, p: int, betaJ: float, nmax: int,
              thL_grid=(0.0, 0.4, 0.8, 1.2), fold_n: int = 1,
              beta_c: float | None = None) -> dict:
    """All fold_bench objects at one lambda, exact N = infinity."""
    if lam < 0.5:
        import warnings
        warnings.warn(
            f"lam = {lam} < 0.5: folded contours exceed float headroom "
            "here (the lam = 0.3 point fails truncation-doubling); folded "
            "outputs are OUTSIDE the validated window", stacklevel=2)
    if beta_c is None:
        beta_c = beta_chord_from_betaJ(betaJ, lam)
    Nflav = 2.0 * p * p / lam
    v = v_of_betaJ(betaJ)
    cc = ChordCorrelators(lam, p, beta_c, nmax)
    fold_base = fold_map(CONFIG_CROSSED, fold_n)

    out = dict(lam=lam, p=p, betaJ=betaJ, beta_chord=beta_c, nmax=nmax,
               Nflav=Nflav, v=v, fold_n=fold_n)

    def put(key, thetas, kind, closed):
        val = cc.connected_ratio(thetas, kind, Nflav)
        out[key] = dict(re=val.real, im=val.imag,
                        closed_re=float(np.real(closed)),
                        closed_im=float(np.imag(closed)))

    put("toc", CONFIG_TOC, 'u', toc_ratio(*CONFIG_TOC, v))
    put("otoc", CONFIG_CROSSED, 'c', otoc_ratio(*CONFIG_CROSSED, v))
    put("fold", fold_base, 'c', otoc_ratio(*fold_base, v))
    for thL in thL_grid[1:]:
        gcfg = (CONFIG_CROSSED[0] + 1j * thL, CONFIG_CROSSED[1] + 1j * thL,
                CONFIG_CROSSED[2], CONFIG_CROSSED[3])
        fcfg = (fold_base[0] + 1j * thL, fold_base[1] + 1j * thL,
                fold_base[2], fold_base[3])
        put(f"growth_thL={thL}", gcfg, 'c', otoc_ratio(*gcfg, v))
        put(f"fold_thL={thL}", fcfg, 'c', otoc_ratio(*fcfg, v))
    # raw connected fractions (no Nflav): the 1/N-hierarchy question
    for key, cfg in [("fold_connected_fraction", fold_base),
                     ("otoc_connected_fraction", CONFIG_CROSSED)]:
        zs = thetas_to_z(cfg, beta_c)
        F = cc.f4_physical(cfg, 'c')
        Fd = cc.g2(zs[0], zs[1]) * cc.g2(zs[2], zs[3])
        frac = (F - Fd) / Fd
        out[key] = dict(re=float(frac.real), im=float(frac.imag))
    return out


def main(quick: bool = False):
    t0 = time.time()
    betaJ = 2.0
    # lambda_chord = 32/N for the fold_scan p = 4 series, extended toward
    # the semiclassical corner
    # Honest cap, stated not silent: below lambda ~ 0.5 the folded contour's
    # e^{+tau E} dynamic range (e^{beta_c * band * excursion} ~ e^{20+})
    # exceeds float cancellation headroom -- the truncated evaluation stops
    # converging in nmax (measured: lam = 0.3 gives 6.6e3 -> 3.8e4 between
    # nmax 220 and 260).  The lambda-trend over [0.5, 4] is unambiguous.
    scan_lams = [4.0, 32.0 / 10, 32.0 / 12, 32.0 / 14, 2.0, 32.0 / 18,
                 1.6]
    extra_lams = [1.2, 0.8, 0.5]
    lams = scan_lams + ([] if quick else extra_lams)
    points = []
    for lam in lams:
        nmax = int(max(24, min(140, 26.0 / lam + 18)))
        t1 = time.time()
        row = run_point(lam, 4, betaJ, nmax)
        row2 = run_point(lam, 4, betaJ, 2 * nmax)   # genuine doubling
        row["truncation_drift_fold"] = abs(
            row["fold"]["re"] - row2["fold"]["re"]) / max(
            abs(row2["fold"]["re"]), 1e-12)
        row["runtime_s"] = round(time.time() - t1, 1)
        points.append(row)
        print(f"lam={lam:.3f}: otoc={row['otoc']['re']:+.3f} "
              f"(closed {row['otoc']['closed_re']:+.3f})  "
              f"fold={row['fold']['re']:+.3f} "
              f"(closed {row['fold']['closed_re']:+.3f})  "
              f"fold_frac={row['fold_connected_fraction']['re']:+.3f}  "
              f"drift={row['truncation_drift_fold']:.1e}", flush=True)
    # matched-lambda ED-bridge row: the (p=6, N=18) <-> (p=4, N=8) anchor at
    # lambda = 4 with the EXACT finite-N sigma_H bridges, for the direct
    # fold_scan comparison columns
    ed_rows = []
    for (p, N) in [(4, 8), (6, 18)]:
        lam = 2.0 * p * p / N
        beta_c = betaJ / (math.sqrt(2.0) * p) * sigma_H(N, p)
        row = run_point(lam, p, betaJ, 60, beta_c=beta_c)
        row["bridge"] = f"exact sigma_H({N},{p})"
        row["N_ed"] = N
        ed_rows.append(row)
        print(f"ED-bridge p={p} N={N}: otoc={row['otoc']['re']:+.3f} "
              f"fold={row['fold']['re']:+.3f}", flush=True)
    payload = dict(
        provenance="chord_fold.py -- rung-26 decider #2: the HZ fold "
                   "evaluated exactly at N = infinity in the chord theory "
                   "(fold_bench.py conventions mirrored; Streicher closed "
                   "forms from antiscrambling.py; ED columns from "
                   "fold_scan.json at matched lambda)",
        betaJ=betaJ, config_crossed=list(CONFIG_CROSSED),
        config_toc=list(CONFIG_TOC),
        fold_thetas=[float(x) for x in fold_map(CONFIG_CROSSED, 1)],
        runtime_s=round(time.time() - t0, 1),
        points=points, ed_bridge_points=ed_rows)
    out = HERE / "chord_fold.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"wrote {out} ({payload['runtime_s']} s)")


if __name__ == "__main__":
    main(quick="--quick" in __import__('sys').argv)
