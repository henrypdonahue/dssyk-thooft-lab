#!/usr/bin/env python3
"""Plot the massless 't Hooft meson Regge trajectory, colored by CP."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from thooft_spectrum import spectrum

levels = spectrum(num_basis=400, num_levels=30)
n = [L["n"] for L in levels]
y = [L["two_lambda"] for L in levels]
cp = [L["CP"] for L in levels]

CP_MINUS = "#2166ac"   # blue
CP_PLUS = "#b2182b"    # red

fig, ax = plt.subplots(figsize=(7.5, 5.0))

# leading 't Hooft/WKB asymptote  2 lambda_n -> n + 3/4
ax.plot(n, [k + 0.75 for k in n], "-", color="0.7", lw=1.2, zorder=1,
        label=r"asymptote $n+\frac{3}{4}$")

for parity_cp, color, lab in [(-1, CP_MINUS, r"CP $=-1$ (symmetric)"),
                              (+1, CP_PLUS, r"CP $=+1$ (antisymmetric)")]:
    xs = [ni for ni, c in zip(n, cp) if c == parity_cp]
    ys = [yi for yi, c in zip(y, cp) if c == parity_cp]
    ax.plot(xs, ys, "o", color=color, ms=7, label=lab, zorder=3)

ax.set_xlabel("level $n$")
ax.set_ylabel(r"$M_n^2/(\pi g^2)\;=\;2\lambda_n$")
ax.set_title("Massless 't Hooft meson spectrum "
             r"($m_q^2=0$, i.e. $\alpha=0$)")
ax.legend(frameon=False, loc="upper left")
ax.grid(True, alpha=0.25)
fig.tight_layout()
fig.savefig("regge_trajectory.png", dpi=150)
print("wrote regge_trajectory.png")
