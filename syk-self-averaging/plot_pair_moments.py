#!/usr/bin/env python3
"""Figure for the exact finite-N pair moments (exact_pair_moments.json).
Reads the JSON only; nothing recomputed.

(a) fixed p = 4: |mu_2| and |mu_4| of the n = 1, 2, 3 pair channels vs N
    (log-log), falling to the conserved-channel theorem's N = infinity
    zero; N^-4 slope guide.
(b) the double-scaling line lambda_chord = 2: |mu_2(n=2)| * pN vs p --
    the observed c/(pN) regularity (c ~ 29-30).
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
    data = json.loads((_HERE / "exact_pair_moments.json").read_text())
except FileNotFoundError as e:
    sys.exit(f"{e.filename} not found -- run exact_pair_moments.py")

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

ax = axes[0]
rows = data["fixed_p"]
ns = np.array([r["N"] for r in rows])
for n, color in [(1, "steelblue"), (2, "crimson"), (3, "darkgreen")]:
    for k, marker, ls in [(2, "o", "-"), (4, "s", "--")]:
        vals = np.array([abs(r["channels"][f"n={n}"][f"mu{k}"])
                         for r in rows])
        ax.loglog(ns, vals, marker=marker, ls=ls, ms=4, color=color,
                  label=f"$|\\mu_{k}|$, $n={n}$")
guide = ns.astype(float) ** -4
guide *= abs(rows[-1]["channels"]["n=2"]["mu2"]) / guide[-1]
ax.loglog(ns, guide, "k:", lw=1, label=r"$N^{-4}$ guide")
ax.set_xlabel(r"$N$   (fixed $p = 4$; $\lambda_{\rm chord} = 32/N$)")
ax.set_ylabel("pair-channel moment (chord units)")
ax.set_title("(a) the measured approach to the proven $N=\\infty$ zero")
ax.legend(fontsize=7, ncol=2)

ax = axes[1]
ds = data["double_scaling"]
ps = np.array([r["p"] for r in ds])
c2 = np.array([abs(r["n2"]["mu2"]) * r["p"] * r["N"] for r in ds])
c3 = np.array([abs(r["n3"]["mu2"]) * r["p"] * r["N"] for r in ds])
ax.plot(ps, c2, "o-", color="crimson", label=r"$|\mu_2(n{=}2)|\,pN$")
ax.plot(ps, c3, "s--", color="darkgreen", label=r"$|\mu_2(n{=}3)|\,pN$")
ax.set_xlabel(r"$p$   ($\lambda_{\rm chord} = 2$, i.e. $N = p^2$)")
ax.set_ylabel(r"$|\mu_2|\, p N$")
ax.set_title(r"(b) $\lambda_{\rm chord} = 2$: the $c/(pN)$ regularity")
ax.legend(fontsize=8)

fig.tight_layout()
out = _HERE / "pair_moments.png"
fig.savefig(out, dpi=160)
print(f"wrote {out}")
