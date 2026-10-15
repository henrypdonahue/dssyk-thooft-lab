#!/usr/bin/env python3
"""Figures for the chord-engine campaigns (CHORD.md).

Reads chord_spectra.json, chord_fold.json and
../syk-self-averaging/mn_majorana_bench.json.  Nothing is recomputed and
nothing is hardcoded; the one curve not read from JSON is the closed form
of full_channel_chord.channel in (d).

(a) The verdict.  Exact N = infinity spectral functions S_Delta, with a
    fine window at the semiclassical scale, against the sech^{2 Delta}
    transform, as lambda falls.  Structureless continua, converging
    linearly in lambda.
(b) The length channel, and the exact-in-lambda de Sitter length growth
    <n(t)> against -(2/lambda) log sech(Jc t).
(c) The decider.  Folded against unfolded connected fractions, against
    lambda.  The fold breaks the 1/N (~lambda) hierarchy at N = infinity.
(d) The bench.  Raw ED w2 of the full A_2 channel along fixed p, plus
    the matched-lambda pair, against the exact N = infinity full channel
    (full_channel_chord.py) and its pair part, which is 0.  The rise along
    fixed p is conserved H^2 weight leaving mu_0, not a channel turning
    on (ERRATA E4; with it removed the width falls monotonically,
    syk-self-averaging/full_channel_laws.py).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from full_channel_chord import channel

_HERE = Path(__file__).resolve().parent

try:
    spectra = json.loads((_HERE / "chord_spectra.json").read_text())
    fold = json.loads((_HERE / "chord_fold.json").read_text())
    bench = json.loads((_HERE.parent / "syk-self-averaging"
                        / "mn_majorana_bench.json").read_text())
except FileNotFoundError as e:
    sys.exit(f"{e.filename} not found -- run the campaigns first")

fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))

# (a) spectral functions, Delta = 0.5 channel
ax = axes[0, 0]
cmap = plt.cm.viridis
pts = sorted(spectra["points"], key=lambda p: -p["lam"])
for i, pt in enumerate(pts):
    ch = pt["channels"]["delta=0.5"]
    om = np.array(ch["om_fine"]) / pt["Jc_bridge"]
    col = cmap(i / max(len(pts) - 1, 1))
    ax.plot(om, np.array(ch["S_fine"]) * pt["Jc_bridge"], color=col,
            label=f"$\\lambda={pt['lam']}$")
    if i == len(pts) - 1:
        ax.plot(om, np.array(ch["S_semiclassical_fine"]) * pt["Jc_bridge"],
                "k--", lw=1, label="sech$^{2\\Delta}$ transform")
ax.set_xlabel(r"$\omega/J_c$")
ax.set_ylabel(r"$J_c\,S_{\Delta=1/2}(\omega)$")
ax.set_title("(a) exact $N=\\infty$ spectra: structureless, "
             "$\\to$ semiclassical $\\propto\\lambda$")
ax.legend(fontsize=7)

# (b) length growth
ax = axes[0, 1]
for i, pt in enumerate(pts):
    g = pt["length_growth"]
    ts = np.array(g["ts"]) * pt["Jc_bridge"]
    col = cmap(i / max(len(pts) - 1, 1))
    lam = pt["lam"]
    ax.plot(ts, lam / 2 * np.array(g["n_of_t"]), color=col,
            label=f"$\\lambda={lam}$")
    if i == len(pts) - 1:
        ax.plot(ts, lam / 2 * np.array(g["n_of_t_semiclassical"]), "k--",
                lw=1, label=r"$\log\cosh(J_c t)$")
ax.set_xlabel(r"$J_c\,t$")
ax.set_ylabel(r"$(\lambda/2)\,\langle n(t)\rangle$")
ax.set_title("(b) de Sitter length growth, exact in $\\lambda$")
ax.legend(fontsize=7)

# (c) fold verdict
ax = axes[1, 0]
fpts = sorted(fold["points"], key=lambda p: p["lam"])
lams = [p["lam"] for p in fpts]
ax.plot(lams, [abs(p["fold_connected_fraction"]["re"]) for p in fpts],
        "o-", color="crimson", label="folded $|(F-F_d)/F_d|$")
ax.plot(lams, [abs(p["otoc_connected_fraction"]["re"]) for p in fpts],
        "s-", color="steelblue", label="unfolded $|(F-F_d)/F_d|$")
ax.plot(lams, [l / 10 for l in lams], "k:", lw=1,
        label=r"$\propto\lambda$ (the $1/N$ hierarchy)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$\lambda$")
ax.set_title("(c) the fold breaks the $1/N$ hierarchy at $N=\\infty$")
ax.legend(fontsize=8)

# (d) bench: raw widths only (no artifact-corrected curves;
# see mn_majorana_bench.py for why)
ax = axes[1, 1]
for p, marker in [(4, "o"), (6, "s")]:
    rows = sorted((r for r in bench["points"]
                   if r["n"] == 2 and r["p"] == p and r["N"] > 2 * p),
                  key=lambda r: r["N"])
    lam = [r["lambda_chord"] for r in rows]
    ax.errorbar(lam, [r["w2"] for r in rows],
                yerr=[r["w2_err"] for r in rows], marker=marker, ls="-",
                label=f"raw $w_2$, $p={p}$ (ED)")
# the near-matched-lambda pair: the honest negative (rises with N)
pair = [next(r for r in bench["points"]
             if r["p"] == 4 and r["N"] == 10 and r["n"] == 2),
        next(r for r in bench["points"]
             if r["p"] == 6 and r["N"] == 22 and r["n"] == 2)]
ax.plot([r["lambda_chord"] for r in pair], [r["w2"] for r in pair],
        "k:", marker="*", ms=11, lw=1,
        label="matched-$\\lambda$ pair ($N{=}10\\to22$)")
lam_grid = np.linspace(1.0, max(r["lambda_chord"] for r in bench["points"]
                               if r["n"] == 2 and r["p"] in (4, 6)), 200)
ax.plot(lam_grid, [channel(x)["w2"] for x in lam_grid], "C3-", lw=1.5,
        label="$N=\\infty$ full channel (closed form)")
ax.axhline(0.0, color="k", lw=1,
           label="$N=\\infty$ pair part (theorem): 0")
ax.set_xlabel(r"$\lambda_{\rm chord}$")
ax.set_ylabel(r"raw $w_2$ ($n=2$ channel)")
ax.set_title("(d) finite-$N$ ED against the exact $N=\\infty$ channel")
ax.legend(fontsize=7)

fig.tight_layout()
out = _HERE / "chord_summary.png"
fig.savefig(out, dpi=160)
print(f"wrote {out}")
