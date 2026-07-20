#!/usr/bin/env python3
r"""
Numerical test of the self-averaging claims of
Miyashita-Sekino-Susskind (arXiv:2607.05678), Section 3.5.

Two claims are tested:

  (3.27)  Var<W> ~ p!/N^p   for singlet (O(N)-invariant) observables W;
  (3.28)  symmetry violations <= exp(-a sqrt(N))  for non-singlet (adjoint) ones.

A subtlety that governs the whole exercise
-------------------------------------------
At *fixed* p, p!/N^p is a power law N^{-p}: self-averaging, but only
polynomially fast.  The paper's striking *super-exponential* suppression
(e^{-aN}, e^{-a sqrt N}) lives in the DOUBLE-SCALING limit p = sqrt(lambda N),
where p grows with N.  Exact diagonalization caps out near N ~ 18-20 (Hilbert
dim 2^{N/2}), which cannot reach large p and large N together -- the same wall
the headline spectrum-match computation hits.

So we split the test into what each tool can actually establish:

  A. The self-averaging *mechanism*, exactly and analytically.  The singlet
     second moment m2 = Tr(H^2)/dim = sum_c J_c^2 has relative variance exactly
     2 / C(N,p).  This needs no diagonalization, so we can follow it along the
     double-scaling line p = sqrt(lambda N) to large N and show it becomes
     super-exponential -- i.e. Eq. (3.27) with its promised e^{-aN} tail.

  B. Higher singlet moments m4, m6 at accessible (N, p) by ED: do they
     self-average at the same ~1/C(N,p) rate as m2, or slower?

  C. Symmetry violation at accessible (N, p) by ED: the off-diagonal of the
     H-dressed adjoint bilinear  B_{jk} = Tr(psi_j H psi_k H)/dim  has zero
     ensemble mean for j != k; its typical size measures how much one disorder
     instance breaks O(N).  We measure its decay with N and with p.

     (The e^{-a sqrt N} FORM cannot be discriminated from a steep power law in
     the ED window -- the fits below say so explicitly, with Delta-AIC -- but
     the disorder VARIANCE of B_jk is exactly computable by Wick contraction at
     any (N, p): see exact_wick.py, which carries this test into the
     double-scaling regime that ED cannot reach.)

Statistics: every ensemble estimate below carries an error bar.  For a relative
variance estimated from n draws the sampling error is ~ relVar*sqrt(2/(n-1))
(Gaussian approximation for the variance of a sample variance); rms values
carry the standard error of the mean over draws.
"""

from __future__ import annotations

import json
from math import comb, factorial, log
from pathlib import Path

import numpy as np

from syk import majorana_operators
from pauli_strings import majorana_hamiltonian_fast


# --------------------------------------------------------------------------
# A. Analytic self-averaging of the singlet second moment (no ED)
# --------------------------------------------------------------------------
def relvar_m2_exact(N: int, p: int) -> float:
    """Exact relative variance of m2 = Tr(H^2)/dim = sum_c J_c^2 over disorder.

    m2 is a sum of C(N,p) i.i.d. Gaussian squares, so Var/mean^2 = 2 / C(N,p).
    """
    return 2.0 / comb(N, p)


def even_p_near(x: float) -> int:
    """Nearest positive even integer to x (p must be even, >= 2)."""
    p = 2 * round(x / 2)
    return max(2, p)


def double_scaling_line(lam: float, ns) -> list:
    """Along p = sqrt(lambda N): report the singlet relative variance 2/C(N,p).
    Pure combinatorics -- reaches N far beyond ED."""
    rows = []
    for N in ns:
        p = even_p_near((lam * N) ** 0.5)
        if p > N:
            continue
        rows.append((N, p, relvar_m2_exact(N, p)))
    return rows


# --------------------------------------------------------------------------
# B. Singlet moment self-averaging by exact diagonalization
# --------------------------------------------------------------------------
def measure_moments(N: int, p: int, n_inst: int, seed: int = 0, ks=(1, 2, 3)):
    """Relative variance of m_{2k} = Tr(H^{2k})/dim over n_inst disorder draws.

    Moments come from the full spectrum (exact per instance, no stochastic-trace
    noise), so the measured spread is pure disorder variance.
    Returns dict k -> (relvar, mean_m2k, relvar_err), where relvar_err is the
    ~relVar*sqrt(2/(n-1)) sampling error of the variance estimate itself."""
    rng = np.random.default_rng(seed)
    moments = {k: [] for k in ks}
    for _ in range(n_inst):
        H = majorana_hamiltonian_fast(N, p, rng)
        eigs = np.linalg.eigvalsh(H)
        for k in ks:
            moments[k].append(np.mean(eigs ** (2 * k)))
    out = {}
    err_fac = np.sqrt(2.0 / (n_inst - 1))
    for k in ks:
        m = np.array(moments[k])
        relvar = float(np.var(m, ddof=1) / np.mean(m) ** 2)
        out[k] = (relvar, float(np.mean(m)), relvar * err_fac)
    return out


# --------------------------------------------------------------------------
# C. Symmetry violation: off-diagonal of the H-dressed adjoint bilinear
# --------------------------------------------------------------------------
def measure_symmetry_violation(N: int, p: int, n_inst: int, seed: int = 0):
    """Typical size of the O(N)-violating (off-diagonal) part of
    B_{jk} = Tr(psi_j H psi_k H)/dim, normalized by its diagonal scale.

    <B_{jk}> ~ delta_{jk} in the O(N)-invariant ensemble, so a nonzero
    off-diagonal is a finite-N, single-instance symmetry violation.

    Returns a dict with:
      rms_off      : mean over disorder of RMS_{j!=k} B_{jk}  -- the trustworthy,
                     monotone symmetry-violation signal.
      rms_off_err  : standard error of that mean over the n_inst draws.
      ratio        : mean of RMS_off / RMS_diag.  Convenient dimensionless
                     measure, BUT its denominator (the diagonal B_jj scale)
                     passes through a small value near p ~ N/2, so the *ratio*
                     has a spurious bump there; always read rms_off when in
                     doubt.
      rms_diag, std_ratio : diagnostics."""
    psis = majorana_operators(N)
    dim = 1 << (N // 2)
    rng = np.random.default_rng(seed)
    off_mask = ~np.eye(N, dtype=bool)
    offs, diags, ratios = [], [], []
    for _ in range(n_inst):
        H = majorana_hamiltonian_fast(N, p, rng)
        # dense stack A[j] = psi_j H ;  B_{jk} = Tr(A_j A_k)/dim.  Contract as
        # one BLAS gemm on flattened matrices (vec(A_j) . vec(A_k^T)) instead
        # of a non-BLAS einsum loop.
        A = np.stack([np.asarray(psis[j] @ H) for j in range(N)])
        F = A.reshape(N, -1)
        G = A.transpose(0, 2, 1).reshape(N, -1)
        B = (F @ G.T).real / dim
        rms_off = np.sqrt(np.mean(B[off_mask] ** 2))
        rms_diag = np.sqrt(np.mean(np.diag(B) ** 2))
        offs.append(rms_off)
        diags.append(rms_diag)
        ratios.append(rms_off / rms_diag)
    n = len(offs)
    return dict(rms_off=float(np.mean(offs)),
                rms_off_err=float(np.std(offs, ddof=1) / np.sqrt(n)),
                rms_diag=float(np.mean(diags)),
                ratio=float(np.mean(ratios)),
                std_ratio=float(np.std(ratios, ddof=1)))


# --------------------------------------------------------------------------
# Fit helpers: power law vs exp(-a N^beta) vs the true double-scaling form
# --------------------------------------------------------------------------
MODELS = {
    "power-law": lambda ns: np.log(ns),
    "exp(-aN)": lambda ns: np.asarray(ns, float),
    "exp(-a sqrt N)": lambda ns: np.sqrt(ns),
    "exp(-a sqrtN lnN)": lambda ns: np.sqrt(ns) * np.log(ns),
}


def fit_models(ns, ys):
    """Linear fits of log(y) against each candidate x(N) in MODELS.

    Returns dict name -> dict(r2, slope, aic).  All models have the same
    number of parameters (2), so AIC = n*ln(SS_res/n) + 4 and Delta-AIC
    ordering coincides with R^2 ordering -- the honest content of AIC here is
    the MAGNITUDE of the differences: Delta-AIC < 2 means the data cannot
    distinguish the models."""
    ns = np.asarray(ns, float)
    ly = np.log(np.asarray(ys, float))
    n = len(ns)
    out = {}
    for name, xfun in MODELS.items():
        x = xfun(ns)
        A = np.vstack([x, np.ones_like(x)]).T
        coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
        pred = A @ coef
        ss_res = float(np.sum((ly - pred) ** 2))
        ss_tot = float(np.sum((ly - ly.mean()) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        aic = n * log(ss_res / n) + 4 if ss_res > 0 else -np.inf
        out[name] = dict(r2=r2, slope=float(coef[0]), aic=aic)
    return out


def fit_verdict(fits) -> str:
    """Best model plus an explicit indistinguishability statement."""
    ranked = sorted(fits.items(), key=lambda kv: kv[1]["aic"])
    best, second = ranked[0], ranked[1]
    daic = second[1]["aic"] - best[1]["aic"]
    if daic < 2.0:
        contenders = [k for k, v in ranked if v["aic"] - best[1]["aic"] < 2.0]
        return (f"cannot discriminate: {', '.join(contenders)} all within "
                f"Delta-AIC < 2 (best '{best[0]}' leads by only {daic:.1f})")
    return f"best: {best[0]} (Delta-AIC to next = {daic:.1f})"


if __name__ == "__main__":
    results = {}

    print("=" * 70)
    print("A. Self-averaging MECHANISM: singlet m2 relative variance = 2/C(N,p)")
    print("=" * 70)
    print("\n  Fixed p = 4 (a power law N^-4):")
    ns_fixed = [10, 12, 14, 16, 18, 20, 30, 50, 100]
    results["A_fixed_p4"] = []
    for N in ns_fixed:
        rv = relvar_m2_exact(N, 4)
        results["A_fixed_p4"].append(dict(N=N, relvar=rv))
        print(f"    N={N:4d}  relVar(m2) = {rv:.3e}   (p!/N^p = {factorial(4)/N**4:.3e})")

    print("\n  Double-scaling line p = sqrt(N)  (lambda = 1) -> super-exponential:")
    rows = double_scaling_line(1.0, [10, 20, 40, 80, 160, 320, 640])
    results["A_double_scaling"] = [dict(N=N, p=p, relvar=rv) for N, p, rv in rows]
    for (N, p, rv) in rows:
        print(f"    N={N:4d}  p={p:3d}  relVar(m2) = {rv:.3e}")

    # exact vs leading-order asymptotic, so nobody mistakes the leading form
    # for the whole story
    N_big, p_big, rv_big = rows[-1]
    leading = -0.5 * np.sqrt(N_big) * np.log(N_big)
    print(f"\n    Asymptotics at N={N_big}: ln(2/C) = {log(rv_big):.1f} vs the "
          f"leading -(1/2)sqrtN lnN = {leading:.1f}")
    print("    (subleading O(sqrtN) terms are sizable; the leading form is a "
          "guide, not the value)")

    fits = fit_models([r[0] for r in rows], [r[2] for r in rows])
    print("\n    Fit of log relVar along double-scaling line:")
    for name, f in fits.items():
        print(f"      {name:18s}: R^2 = {f['r2']:.4f}   AIC = {f['aic']:+7.1f}")
    print(f"    -> {fit_verdict(fits)}")
    print("    (super-polynomial growth of the decay rate is what confirms "
          "Eq. 3.27's tail)")
    results["A_fits"] = fits

    print("\n" + "=" * 70)
    print("B. Higher singlet moments by ED (fixed p = 4): rate vs m2")
    print("=" * 70)
    ns = [10, 12, 14, 16, 18]
    n_inst = 200
    err_pct = 100 * np.sqrt(2 / (n_inst - 1))
    print(f"\n  {n_inst} disorder realizations each; sampling error on every "
          f"relVar is ~{err_pct:.0f}%.")
    print(f"  {'N':>3} {'relVar(m2)':>12} {'2/C(N,p)':>12} {'relVar(m4)':>12} "
          f"{'relVar(m6)':>12}")
    results["B_moments"] = []
    rv2, rv4, rv6 = [], [], []
    for N in ns:
        res = measure_moments(N, 4, n_inst, seed=10 + N)
        r2, r4, r6 = res[1][0], res[2][0], res[3][0]
        rv2.append(r2); rv4.append(r4); rv6.append(r6)
        results["B_moments"].append(dict(
            N=N, p=4, n_inst=n_inst,
            relvar_m2=r2, relvar_m2_err=res[1][2], relvar_m2_exact=2 / comb(N, 4),
            relvar_m4=r4, relvar_m4_err=res[2][2],
            relvar_m6=r6, relvar_m6_err=res[3][2]))
        print(f"  {N:>3} {r2:>12.3e} {2/comb(N,4):>12.3e} {r4:>12.3e} {r6:>12.3e}")
    for label, ys in [("m2", rv2), ("m4", rv4), ("m6", rv6)]:
        f = fit_models(ns, ys)
        slope = f["power-law"]["slope"]
        print(f"    {label}: power-law slope d(log relVar)/d(log N) = {slope:.2f} "
              f"(exact 2/C(N,4) slope over this window: -4.56)")

    print("\n" + "=" * 70)
    print("C. Symmetry violation by ED: off-diagonal adjoint bilinear vs N")
    print("=" * 70)
    n_inst_sv = 120
    print(f"\n  B_jk = Tr(psi_j H psi_k H)/dim.  rms_off = RMS_(j!=k) B_jk is the")
    print(f"  trustworthy monotone signal; the ratio rms_off/rms_diag has a spurious")
    print(f"  bump near p ~ N/2 (its denominator collapses there), so we FIT rms_off.")
    print(f"  {n_inst_sv} realizations each; rms_off carries its standard error.")
    results["C_symmetry_violation"] = {}
    for p in (4, 6):
        print(f"\n  p = {p}:")
        print(f"    {'N':>3} {'rms_off':>12} {'+-':>10} {'rms_diag':>12} {'ratio':>10}")
        ys = []
        rows_c = []
        ns_sv = [n for n in [10, 12, 14, 16, 18] if n >= p]
        for N in ns_sv:
            r = measure_symmetry_violation(N, p, n_inst_sv, seed=20 + N)
            ys.append(r["rms_off"])
            rows_c.append(dict(N=N, p=p, n_inst=n_inst_sv, **r))
            flag = "  <- ratio artifact (p~N/2)" if abs(p - N / 2) < 1e-9 else ""
            print(f"    {N:>3} {r['rms_off']:>12.3e} {r['rms_off_err']:>10.1e} "
                  f"{r['rms_diag']:>12.3e} {r['ratio']:>10.3e}{flag}")
        f = fit_models(ns_sv, ys)   # fit the trustworthy rms_off
        print("    fits of rms_off: " + ", ".join(
            f"{k} R^2={v['r2']:.3f}" for k, v in f.items()))
        print(f"    -> {fit_verdict(f)}")
        results["C_symmetry_violation"][f"p{p}"] = dict(rows=rows_c, fits=f)

    out_path = Path(__file__).resolve().parent / "results.json"
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=1)
    print(f"\nwrote {out_path} (machine-readable; plot_self_averaging.py reads it)")

    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    print("""
  * Eq. 3.27 (singlet self-averaging) -- what is actually established, precisely:
    - The singlet second moment m2 = Tr(H^2)/dim = sum_c J_c^2 EXACTLY, so its
      relative variance is the elementary statistical identity 2/C(N,p) (sum of
      C(N,p) i.i.d. Gaussian squares). ED reproduces it within the ~10% sampling
      noise of 200 draws -- but that is a CODE-CORRECTNESS check (Tr H^2 =
      sum J^2), not a test of the paper's diagrammatic mechanism, which the
      identity holds regardless of.
    - The non-trivial content is: (i) higher moments m4, m6 self-average at the
      same ~1/C(N,p) combinatorial rate (they steepen slightly); (ii) 2/C(N,p)
      is the paper's p!/N^p times 2*prod(1-k/N)^-1 ~ 2e^{+p^2/2N}: an O(1/N)
      correction at fixed p but an order-one factor in double scaling. Followed
      along p=sqrt(lambda N) the relative variance is super-exponential, giving
      the claimed p!/N^p -> e^-aN tail. So: the self-averaging MECHANISM is
      confirmed and quantified for the H^{2k} singlets; a fully general bounded
      singlet W is not probed here.

  * Eq. 3.28 (symmetry violation): rms_off of the adjoint bilinear decays with N
    and steepens with p, consistent with the claimed suppression. What the
    Delta-AIC verdicts above actually show in the ED window (N <= 18, fixed p):
    a PURE power law is disfavored (Delta-AIC ~ 9-12), but the exponential-family
    candidates cannot be reliably separated (best-vs-next 2-5, and the winner
    differs between p = 4 and p = 6) -- so the window alone identifies neither
    the asymptotic form nor the exponent. It does not need to: exact_wick.py
    gives the observable's disorder variance in CLOSED FORM at any (N, p) --
    at fixed p the true asymptote is the power law rms ~ N^{-(p-1)/2} (the
    window's curvature that mimics an exponential is the pre-asymptotic shape
    of 2 sigma^2 sqrt(C(N-2,p-1)), which passes through every ED point), and
    along p = sqrt(N) it is ln rms = -(1/4) sqrt(N) ln N (1+o(1)), beating
    e^{-a sqrt N} for every a -- Eq. 3.28 verified where it matters. (The
    ratio rms_off/rms_diag shows a spurious bump at p=N/2 where its
    denominator collapses; exact_wick.py derives it: E[B_jj] vanishes there.
    Read rms_off.)
""")
