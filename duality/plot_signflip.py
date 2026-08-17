#!/usr/bin/env python3
"""Figure for the antipodal-map probes (signflip_probe.json,
signflip_deep.json).  Reads the JSONs only; nothing recomputed.

(a) the rotation law: Re a(Delta) around the half fake circle, closed
    form -- the exact flip at the antipode;
(b) the exact-theory onset: re_move of the antipodal trend vs lambda
    (up = anti-scrambling), with the unshifted control;
(c) the hologram point: the transient anti-scrambling window vs lambda
    (in Jc units), with its ~ half ln(1/lambda) growth;
(d) the frame's thermality: per-line detailed-balance beta_eff vs
    beta + beta_fake (exact identity), all lambda.
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
    probe = json.loads((_HERE / "signflip_probe.json").read_text())
    deep = json.loads((_HERE / "signflip_deep.json").read_text())
except FileNotFoundError as e:
    sys.exit(f"{e.filename} not found -- run signflip_probe.py [--deep]")

fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))

# (a) rotation law
ax = axes[0, 0]
for cf in probe["closed_form"]:
    ds = [r["delta"] * cf["v"] / np.pi for r in cf["rows"]]
    ax.plot(ds, [r["re_a"] for r in cf["rows"]], "o-",
            label=f"$\\beta J = {cf['betaJ']}$")
ax.axhline(-1, color="k", lw=0.8, ls=":")
ax.axhline(+1, color="k", lw=0.8, ls=":")
ax.set_xlabel(r"shift $\Delta$ in units of the half fake period $\pi/v$")
ax.set_ylabel(r"Re $a(\Delta)$")
ax.set_title("(a) closed form: the exact flip at the antipode")
ax.legend(fontsize=8)

# (b) exact-theory onset
ax = axes[0, 1]
lams = [c["lam"] for c in probe["chord"]]
anti = [max(c["rows"], key=lambda r: r["delta"])["re_move"]
        for c in probe["chord"]]
uns = [c["rows"][0]["re_move"] for c in probe["chord"]]
ax.plot(lams, anti, "o-", color="crimson", label="antipodal frame")
ax.plot(lams, uns, "s-", color="steelblue", label="unshifted control")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"late-time move of Re $\mathcal{F}/F_d$")
ax.set_title("(b) exact theory: flip on for $\\lambda \\lesssim 1.7$")
ax.legend(fontsize=8)

# (c) hologram-point window
ax = axes[1, 0]
hl = sorted(deep["hologram"], key=lambda r: r["lam"])
lams_h = np.array([r["lam"] for r in hl])
wins = np.array([r["window_in_Jc_units"] for r in hl])
ax.plot(lams_h, wins, "o-", color="darkgreen", label="measured window")
lam_fit = np.linspace(lams_h.min(), lams_h.max(), 50)
ax.plot(lam_fit, 0.5 * np.log(1.0 / lam_fit) + wins[-1]
        - 0.5 * np.log(1.0 / lams_h[-1]), "k--", lw=1,
        label=r"$\frac{1}{2}\ln(1/\lambda)$ + const")
ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"anti-scrambling window  $t^*\,J_c$")
ax.set_title("(c) hologram point: a semiclassical dS epoch")
ax.legend(fontsize=8)

# (d) detailed balance
ax = axes[1, 1]
rows = deep.get("detailed_balance", probe.get("detailed_balance", []))
if not rows:
    rows = probe["detailed_balance"]
x = np.arange(len(rows))
ax.bar(x - 0.18, [r["beta_eff"] for r in rows], 0.36,
       label=r"measured $\beta_{\rm eff}$")
ax.bar(x + 0.18, [r["beta_c"] + r["beta_fake_c"] for r in rows], 0.36,
       label=r"$\beta + \beta_{\rm fake}$")
ax.set_xticks(x)
ax.set_xticklabels([f"$\\lambda={r['lam']}$\n$\\beta J={r['betaJ']}$"
                    for r in rows], fontsize=7)
ax.set_title("(d) the frame is thermal at the tomperature (exact)")
ax.legend(fontsize=8)

fig.tight_layout()
out = _HERE / "signflip_summary.png"
fig.savefig(out, dpi=160)
print(f"wrote {out}")
