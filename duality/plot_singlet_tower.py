#!/usr/bin/env python3
"""Figure for the singlet-tower result.  Reads singlet_tower.json only;
nothing is recomputed except the two closed-form curves of panel (a).

(a) The A_2 width against lambda_chord: the identity-subtracted ratio
    mu_2/mu_0 and the dynamical width 2(1 - q^2), with the exact finite-N
    1/p limits and the dense-ED commutant-projected points.
(b) The exact dynamical spectral densities of A_2, A_3, A_4 at
    lambda_chord = 1, each normalized to unit area.
(c) Which matter chord carries the dynamical weight of A_2..A_6 at
    lambda_chord = 2.
(d) The exact finite-N cross moment <ad^2(A_4) A_2^dag> on
    lambda_chord = 1 against 1/p, with the N = infinity value with and
    without the remainder term.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
with open(HERE / "singlet_tower.json") as fh:
    D = json.load(fh)

# validated categorical order (fixed, never cycled) and neutral ink
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
SURFACE, INK, INK2, MUTED, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#c3c2b7"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "text.color": INK, "axes.titlesize": 10.5, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.labelsize": 9.5, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#e9e8e3", "grid.linewidth": 0.8, "lines.linewidth": 2.0,
    "legend.frameon": False, "legend.fontsize": 8.5,
})


def ring(ax, x, y, color, marker="o", size=7, **kw):
    """A marker with a surface-colored ring so it stays legible on a line."""
    ax.plot(x, y, marker, ms=size, mfc=color, mec=SURFACE, mew=1.6,
            ls="none", zorder=4, **kw)


fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.8))

# ---- (a) the width, raw and resolved --------------------------------------
ax = axes[0, 0]
lam = np.linspace(0.02, 5.0, 400)
q = np.exp(-lam)
raw = 16 * lam * (1 - q) * (1 - q * q + 2 * lam * q) / (
    4 * lam * (lam * (1 + q) + 2 * (1 - q)))
dyn = 2 * (1 - q * q)
ax.plot(lam, dyn, color=SERIES[0], label=r"dynamical: $2(1-q^2)$")
ax.plot(lam, raw, color=SERIES[1], label=r"identity-subtracted $\mu_2/\mu_0$")
gx = [g["lam"] for g in D["a2_gate"].values()]
gy = [g["limits"]["w2_dynamic"]["central"] for g in D["a2_gate"].values()]
ring(ax, gx, gy, SERIES[0], marker="s", label="exact finite $N$, $1/p$ limit")
ed = [r for r in D["ed_commutant"] if r["p"] == 6]
ax.errorbar([r["lam_chord"] for r in ed], [r["w2_dynamic"] for r in ed],
            yerr=[r["w2_dynamic_err"] for r in ed], fmt="none", ecolor=INK2,
            elinewidth=1, zorder=3)
ring(ax, [r["lam_chord"] for r in ed], [r["w2_dynamic"] for r in ed],
     SERIES[0], marker="o", label="dense ED, $p=6$: commutant removed")
ring(ax, [r["lam_chord"] for r in ed], [r["w2_raw"] for r in ed],
     SERIES[1], marker="o", label="dense ED, $p=6$: raw")
ax.annotate("conserved $H^2$ weight\ndilutes the raw ratio", xy=(2.6, 1.12),
            xytext=(2.75, 1.45), color=INK2, fontsize=8.5,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.set_xlabel(r"$\lambda_{\rm chord}=2p^2/N$")
ax.set_ylabel(r"width of the full $A_2$ channel (chord units)")
ax.set_title(r"(a) $A_2$ at $N=\infty$: a conserved piece plus a $\Delta=2$ chord")
ax.set_xlim(0, 5.0)
ax.set_ylim(0, 2.25)
ax.grid(True, axis="y")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2,
          columnspacing=1.6, handlelength=1.6)

# ---- (b) spectral densities ------------------------------------------------
ax = axes[0, 1]
sp = D["spectra"]
om = np.array(sp["omega"])
for i, n in enumerate(("2", "3", "4")):
    lev = sp["levels"][n]
    S = np.array(lev["total"]) / lev["mu0_exact"]
    ax.plot(om, S, color=SERIES[i], label=f"$A_{n}$")
    j = int(np.argmax(S)) if n != "2" else 12
    ax.annotate(f"$A_{n}$", xy=(om[j], S[j]), xytext=(4, 6),
                textcoords="offset points", color=INK, fontsize=9)
ax.axvline(sp["band"], color=AXIS, lw=1.0)
ax.annotate("half-width of\nthe $H$ band", xy=(sp["band"], 0.0),
            xytext=(sp["band"] + 0.08, 0.262), color=INK2, fontsize=8)
ax.set_xlabel(r"boost frequency $\omega$ (chord units)")
ax.set_ylabel(r"$S(\omega)\,/\,\mu_0^{\rm dyn}$")
ax.set_title(rf"(b) exact dynamical spectral densities, $\lambda_{{\rm chord}}={sp['lam']:g}$")
ax.set_xlim(0, om[-1])
ax.set_ylim(0, None)
ax.grid(True, axis="y")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3,
          columnspacing=1.6, handlelength=1.6)

# ---- (c) which chord carries the weight -----------------------------------
ax = axes[1, 0]
row = next(r for r in D["scan"] if r["lam"] == 2.0)
levels = sorted(row["channels"], key=int)
ypos = np.arange(len(levels))[::-1]
for y, n in zip(ypos, levels):
    left = 0.0
    shares = row["channels"][n]["sector_share"]
    for m in sorted(shares, key=int):
        w = 100 * shares[m]
        ax.barh(y, w, left=left, height=0.46, color=SERIES[int(m) - 2],
                edgecolor=SURFACE, linewidth=2.0,
                label=f"$O_{m}$" if n == levels[-1] else None)
        if w > 9:
            ax.text(left + w / 2, y, f"{w:.0f}%", ha="center", va="center",
                    color="white" if int(m) in (2, 3) else INK, fontsize=8.5)
        left += w
    new = 100 * shares[n]
    if n != "2":
        ax.text(101.5, y, f"new chord $O_{n}$: {new:.2g}%", va="center",
                color=INK2, fontsize=8.5)
ax.set_yticks(ypos)
ax.set_yticklabels([f"$A_{n}$" for n in levels])
ax.set_xlim(0, 138)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("share of the dynamical weight (%)")
ax.set_title(r"(c) what each level is made of, $\lambda_{\rm chord}=2$")
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.legend(loc="upper center", bbox_to_anchor=(0.36, -0.17), ncol=5,
          handlelength=1.2, columnspacing=1.4)

# ---- (d) the remainder gate ------------------------------------------------
ax = axes[1, 1]
g = D["gates"]["cross42_k2_lam1"]
ctl = D["remainder_control"]["lam1"]
x = [1.0 / r["p"] for r in g["rows"]]
y = [r["mu"] for r in g["rows"]]
ax.axhline(ctl["cross42_with"], color=SERIES[0], lw=2.0,
           label=f"$N=\\infty$, remainder term kept: {ctl['cross42_with']:.2f}")
ax.axhline(ctl["cross42_without"], color=SERIES[1], lw=2.0,
           label=f"$N=\\infty$, remainder dropped: {ctl['cross42_without']:.2f}")
xs = np.linspace(0.0, max(x), 200)
order = g["limit"]["fits"][-1]["order"]
pts = [(1.0 / r["p"], r["mu"]) for r in g["rows"] if r["p"] >= 6]
coef = np.polyfit([a for a, _ in pts], [b for _, b in pts], order)
ax.plot(xs, np.polyval(coef, xs), color=MUTED, lw=1.0, zorder=2,
        label=f"order-{order} fit in $1/p$: {g['limit']['central']:.2f}")
ring(ax, x, y, INK2, label="exact finite $N$ (rational arithmetic)")
ax.set_xlim(0, 0.27)
ax.set_xlabel(r"$1/p$   ($\lambda_{\rm chord}=1$, $N=2p^2$)")
ax.set_ylabel(r"$\langle{\rm ad}_H^2(A_4)\,A_2^\dagger\rangle$")
ax.set_title("(d) exact finite $N$ selects the remainder rule")
ax.grid(True, axis="y")
ax.legend(loc="lower left")

fig.tight_layout(h_pad=1.6, w_pad=2.0)
out = HERE / "singlet_tower.png"
fig.savefig(out, dpi=150)
print(f"wrote {out.name}")
