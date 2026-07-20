#!/usr/bin/env python3
"""Figures for the self-averaging test: (left) singlet relative variance vs N at
fixed p (power law) and along the double-scaling line (super-exponential);
(right) the symmetry-violation signal rms_off vs sqrt(N) for p = 4, 6, with the
EXACT closed-form curve of exact_wick.py drawn through the ED points (no fit).

Reads results.json (written by self_averaging.py) from the module directory,
so the figure always shows the production numbers with their error bars --
nothing is hardcoded and the scripts work from any cwd.  The right panel plots
rms_off, the signal self_averaging.py designates as trustworthy; the ratio
rms_off/rms_diag (whose p ~ N/2 denominator collapse exact_wick.py now
derives) is deliberately NOT plotted.
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from self_averaging import relvar_m2_exact, double_scaling_line
from exact_wick import rms_offdiag_exact

_HERE = Path(__file__).resolve().parent

try:
    with open(_HERE / "results.json") as fh:
        results = json.load(fh)
except FileNotFoundError:
    sys.exit("results.json not found -- run  python3 self_averaging.py  first")

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.8))

# ---- left: singlet self-averaging (exact combinatorics) ------------------
ns_fixed = np.array([10, 12, 14, 16, 18, 20, 24, 30, 40, 60, 100])
rv_fixed = [relvar_m2_exact(N, 4) for N in ns_fixed]
axL.loglog(ns_fixed, rv_fixed, "o-", color="#2166ac",
           label=r"fixed $p=4$: relVar$(m_2)=2/\binom{N}{4}\sim N^{-4}$")

ds = double_scaling_line(1.0, [10, 20, 40, 80, 160, 320, 640])
axL.loglog([r[0] for r in ds], [r[2] for r in ds], "s-", color="#b2182b",
           label=r"double-scaling $p=\sqrt{N}$: super-exponential")

# ED-measured relVar(m2) with error bars, on top of the exact fixed-p curve
bm = results["B_moments"]
axL.errorbar([r["N"] for r in bm], [r["relvar_m2"] for r in bm],
             yerr=[r["relvar_m2_err"] for r in bm], fmt="x", color="#333333",
             capsize=3, label=r"ED, $200$ draws (fixed $p=4$)")

axL.set_xlabel("$N$")
axL.set_ylabel(r"relative variance of singlet $m_2$")
axL.set_title("Singlet self-averaging (Eq. 3.27)")
axL.legend(frameon=False, fontsize=9)
axL.grid(True, which="both", alpha=0.2)

# ---- right: symmetry violation (ED points + exact closed form) -----------
axR.set_title("Symmetry violation (Eq. 3.28)")
for p, color in [(4, "#2166ac"), (6, "#b2182b")]:
    rows = results["C_symmetry_violation"][f"p{p}"]["rows"]
    ns_sv = np.array([r["N"] for r in rows])
    ys = np.array([r["rms_off"] for r in rows])
    errs = np.array([r["rms_off_err"] for r in rows])
    axR.errorbar(np.sqrt(ns_sv), ys, yerr=errs, fmt="o", color=color,
                 capsize=3, label=f"$p={p}$ (ED)")
    n_dense = np.arange(ns_sv[0], ns_sv[-1] + 1, 2)
    axR.plot(np.sqrt(n_dense), [rms_offdiag_exact(int(N), p) for N in n_dense],
             "-", color=color, lw=1.0, alpha=0.7,
             label=f"$p={p}$ exact: $2\\sigma^2\\sqrt{{\\binom{{N-2}}{{{p-1}}}}}$")
axR.set_yscale("log")
axR.set_xlabel(r"$\sqrt{N}$")
axR.set_ylabel(r"RMS$_{j\neq k}\,B_{jk}$   (adjoint bilinear)")
axR.legend(frameon=False, fontsize=8)
axR.grid(True, which="both", alpha=0.2)
axR.text(0.03, 0.03,
         "curves: exact closed form (exact_wick.py), no fit.\n"
         "Fixed-$p$ asymptote is $N^{-(p-1)/2}$; along $p=\\sqrt{N}$ the exact\n"
         "form beats $e^{-a\\sqrt{N}}$ for every $a$ (Eq. 3.28, settled).",
         transform=axR.transAxes, fontsize=8, color="0.4")

fig.tight_layout()
out = _HERE / "self_averaging.png"
fig.savefig(out, dpi=150)
print(f"wrote {out}")
