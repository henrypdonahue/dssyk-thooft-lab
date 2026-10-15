#!/usr/bin/env python3
"""Figure for the exact finite-N pair moments (exact_pair_moments.json).
Reads the JSON only.  Nothing is recomputed.

(a) Fixed p = 4.  |mu_2| and |mu_4| of the n = 1, 2, 3 pair channels
    against N, log-log, falling to the conserved-channel theorem's
    N = infinity zero.  With an N^-4 slope guide.
(b) The double-scaling line lambda_chord = 2.  |mu_2(n=2)| * pN against
    p, the c/(pN) regularity: the float campaign's p <= 10 points, and
    the exact points to p = 48 from pair_moment_laws.json /
    pair_moment_ladder.json with the closed-form limit
    4 lambda (1-q)^2 (5+3q) = 32.334 (double_scaling_law.py) as a guide.
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
# the exact points beyond the float campaign, and the closed-form limits
exact2, exact3 = {}, {}
for name in ("pair_moment_laws.json", "pair_moment_ladder.json"):
    path = _HERE / name
    if path.exists():
        with open(path) as fh:
            dd = json.load(fh).get("double_scaling", {})
        exact2.update({r["p"]: r["c"] for r in dd.get("n2", [])})
        exact3.update({r["p"]: r["c"] for r in dd.get("n3", [])})
path = _HERE / "double_scaling_law.json"
if path.exists():
    with open(path) as fh:
        dd = json.load(fh).get("sequences", {})
    exact3.update({r["p"]: r["c"] for r in dd.get("lam2_n3", [])})
if exact2:
    ax.plot(sorted(exact2), [exact2[p] for p in sorted(exact2)], ".",
            color="crimson", ms=4, label="exact, $n=2$ (laws + ladder)")
if exact3:
    ax.plot(sorted(exact3), [exact3[p] for p in sorted(exact3)], ".",
            color="darkgreen", ms=4, label="exact, $n=3$")
q = np.exp(-2.0)
ax.axhline(8 * (1 - q) ** 2 * (5 + 3 * q), color="crimson", lw=0.8, ls=":",
           label=r"$4\lambda(1-q)^2(5+3q) = 32.334$")
ax.axhline(8 * (1 - q) ** 3 * (35 + 56 * q + 35 * q * q + 10 * q ** 3),
           color="darkgreen", lw=0.8, ls=":",
           label=r"$4\lambda(1-q)^3(35+56q+35q^2+10q^3) = 223.65$")
ax.set_yscale("log")
ax.set_xlabel(r"$p$   ($\lambda_{\rm chord} = 2$, i.e. $N = p^2$)")
ax.set_ylabel(r"$|\mu_2|\, p N$")
ax.set_title(r"(b) $\lambda_{\rm chord} = 2$: $c/(pN)$ and its closed-form limits")
ax.legend(fontsize=6.5)

fig.tight_layout()
out = _HERE / "pair_moments.png"
fig.savefig(out, dpi=160)
print(f"wrote {out}")
