#!/usr/bin/env python3
"""Figures for the self-averaging test: (left) singlet relative variance vs N at
fixed p (power law) and along the double-scaling line (super-exponential);
(right) symmetry-violation ratio vs N for p = 4, 6."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from self_averaging import relvar_m2_exact, double_scaling_line

# Symmetry-violation data measured by self_averaging.py (120 realizations each,
# saved in results.txt). Read from there rather than recomputed so the figure is
# instant and uses the production-quality statistics. N=12,p=6 is a noisy
# small-N outlier (std ~ 0.21) that breaks monotonicity; shown as an open marker.
SV = {
    4: ([10, 12, 14, 16, 18],
        [3.225e-01, 1.313e-01, 6.834e-02, 4.162e-02, 2.769e-02]),
    6: ([10, 12, 14, 16, 18],
        [3.271e-01, 7.700e-01, 1.280e-01, 4.415e-02, 2.114e-02]),
}
SV_OUTLIERS = {6: {12}}

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.8))

# ---- left: singlet self-averaging ---------------------------------------
ns_fixed = np.array([10, 12, 14, 16, 18, 20, 24, 30, 40, 60, 100])
rv_fixed = [relvar_m2_exact(N, 4) for N in ns_fixed]
axL.loglog(ns_fixed, rv_fixed, "o-", color="#2166ac",
           label=r"fixed $p=4$: relVar$(m_2)=2/\binom{N}{4}\sim N^{-4}$")

ds = double_scaling_line(1.0, [10, 20, 40, 80, 160, 320, 640])
ns_ds = [r[0] for r in ds]
rv_ds = [r[2] for r in ds]
axL.loglog(ns_ds, rv_ds, "s-", color="#b2182b",
           label=r"double-scaling $p=\sqrt{N}$: super-exponential")
axL.set_xlabel("$N$")
axL.set_ylabel(r"relative variance of singlet $m_2$")
axL.set_title("Singlet self-averaging (Eq. 3.27)")
axL.legend(frameon=False, fontsize=9)
axL.grid(True, which="both", alpha=0.2)

# ---- right: symmetry violation ------------------------------------------
axR.set_title("Symmetry violation (Eq. 3.28)")
for p, color in [(4, "#2166ac"), (6, "#b2182b")]:
    ns_sv, ys = SV[p]
    ns_sv, ys = np.array(ns_sv), np.array(ys)
    good = np.array([n not in SV_OUTLIERS.get(p, set()) for n in ns_sv])
    axR.semilogy(np.sqrt(ns_sv[good]), ys[good], "o-", color=color, label=f"$p={p}$")
    if (~good).any():
        axR.semilogy(np.sqrt(ns_sv[~good]), ys[~good], "o", mfc="none",
                     color=color, label=f"$p={p}$ (small-$N$ outlier)")
axR.set_xlabel(r"$\sqrt{N}$")
axR.set_ylabel(r"RMS$_{\rm off}$ / RMS$_{\rm diag}$ of $B_{jk}$")
axR.legend(frameon=False, fontsize=9)
axR.grid(True, which="both", alpha=0.2)
axR.text(0.03, 0.03, r"straight line here $\Rightarrow e^{-a\sqrt{N}}$",
         transform=axR.transAxes, fontsize=8, color="0.4")

fig.tight_layout()
fig.savefig("self_averaging.png", dpi=150)
print("wrote self_averaging.png")
