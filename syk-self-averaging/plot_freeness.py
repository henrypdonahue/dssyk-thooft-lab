#!/usr/bin/env python3
"""Figure for the freeness-onset measurement (freeness.py, rung 27).

Reads freeness.json only -- nothing recomputed, nothing hardcoded.

(a) spec(E + F(t)) at the largest N: the exact binomial atoms at t = 0
    melt into the ARCSINE law (the free convolution) by a few 1/J;
(b) the order parameter: W1 distances of the spectrum to the arcsine
    (falling, solid) and binomial (rising, dashed) laws vs tJ, per N;
(c) the verdict: late-time plateaus vs Hilbert dimension -- a falling
    trend supports the free product in the large-N limit.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = Path(__file__).resolve().parent

try:
    data = json.loads((_HERE / "freeness.json").read_text())
except FileNotFoundError:
    sys.exit("freeness.json not found -- run  python3 freeness.py  first")

points = data["points"]
W1_REF = data["w1_binomial_vs_arcsine"]
BLUES = ["#9ecae1", "#6baed6", "#3182bd", "#08519c"]     # light -> dark with N

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15, 4.4))

# ---- (a) the spectrum becomes the arcsine law (largest N) ----------------
big = points[-1]
edges = np.array(big["hist_edges"])
mids = 0.5 * (edges[1:] + edges[:-1])
width = edges[1] - edges[0]
n_eigs = big["n_inst"] * big["dim"]
snap_cols = {0: "#999999", 1: BLUES[1], 2: BLUES[3]}
for k, (tJ, counts) in enumerate(sorted(big["spec_hists"].items(),
                                        key=lambda kv: float(kv[0]))):
    dens = np.array(counts) / (n_eigs * width)
    axA.stairs(dens, edges, color=snap_cols[k], lw=1.8,
               label=rf"$tJ = {float(tJ):g}$")
x = np.linspace(-1.999, 1.999, 400)
axA.plot(x, 1.0 / (np.pi * np.sqrt(4.0 - x ** 2)), "--", color="#b2182b",
         lw=2.0, label="arcsine (free)")
axA.set_ylim(0, 1.0)
axA.set_xlabel(r"$\mathrm{spec}(E + F(t))$")
axA.set_ylabel("spectral density")
axA.set_title(rf"(a) classical atoms $\to$ free arcsine law "
              rf"($N={big['N']}$, dim {big['dim']})")
axA.legend(frameon=False, fontsize=9)

# ---- (b) the order parameter vs time -------------------------------------
for c, pt in zip(BLUES, points):
    tJ = np.array(pt["tJ"])
    axB.errorbar(tJ, pt["W1_arc"], yerr=pt["W1_arc_sem"], fmt="-o", ms=3,
                 lw=1.8, color=c, label=rf"$N={pt['N']}$")
    axB.plot(tJ, pt["W1_bin"], "--", lw=1.2, color=c)
axB.axhline(W1_REF, color="#999999", ls=":", lw=1.2)
axB.text(0.4, W1_REF + 0.008, r"$W_1$(binomial, arcsine)", fontsize=8,
         color="#666666")
axB.set_xlabel(r"$t\,\mathcal{J}$")
axB.set_ylabel(r"$W_1$ of $\mathrm{spec}(E+F(t))$")
axB.set_title(r"(b) $\to$ arcsine (solid), $\to$ binomial (dashed)")
axB.legend(frameon=False, fontsize=9)

# ---- (c) the N-trend of the late-time plateaus ---------------------------
dims = np.array([pt["dim"] for pt in points], float)
for key, marker, label in (("plateau_W1_arc", "o", r"$W_1$ to arcsine"),
                           ("plateau_m2", "s", r"$|\tau[(EF)^2]|$"),
                           ("plateau_w4", "^", r"$|\tau[abab]|$")):
    axC.loglog(dims, [pt[key] for pt in points], marker + "-", lw=1.8,
               color={"plateau_W1_arc": BLUES[3], "plateau_m2": "#b2182b",
                      "plateau_w4": "#666666"}[key], label=label)
guide = points[0]["plateau_W1_arc"] * dims[0] / dims
axC.loglog(dims, guide, ":", color="#999999", lw=1.2,
           label=r"$\propto \mathrm{dim}^{-1}$ (rigid spectrum)")
axC.set_xlabel("Hilbert dimension")
axC.set_ylabel("late-time plateau")
axC.set_title("(c) freeness deepens with system size")
axC.legend(frameon=False, fontsize=9)

fig.tight_layout()
out = _HERE / "freeness.png"
fig.savefig(out, dpi=160)
print(f"wrote {out}")
