#!/usr/bin/env python3
r"""
Anti-scrambling and DSSYK_infinity: the OTOC sign test, and the fold.
=====================================================================

Two papers from July 2026 sharpen what a de Sitter static-patch hologram
must do at OTOC level.

CUI-KOLCHMEYER 2607.13665, "A de Sitter Anti-Scrambling Algebra".  The dS
observer's shockwave S-matrix has a NEGATIVE eikonal phase, a time
ADVANCE.  Their Eq. (2.74) writes the semiclassical OTOC as

    OTOC/TOC = 1 - i x h_A h_B e^{(t2+t4-t1-t3)/2}
                    / (sinh(t13/2) sinh(t24/2)),   x = c e^T,

with c > 0 for AdS (scrambling, time delay) and c < 0 for dS
(anti-scrambling).

The Maldacena-Shenker-Stanford chaos bound forces c > 0 for any unitary
quantum system in a KMS state at the dS temperature.  dS correlators
violate the KMS assumption, and CK conclude that such holograms are ruled
out.

(2607.05678, this repo's target paper, is cited there as [12].  Its
T_B = infinity hologram is NOT a KMS state at the dS temperature, so it
dodges the letter of the no-go.  Whether it can produce c < 0 through a
nontrivial dictionary is exactly the open question this module
quantifies.)

HARLOW-ZHAO 2607.14215, "Anti-scrambling and euclidean folds".  The dS
observer 2-pt function is thermal at beta_dS, and the 4-pt function
anti-scrambles, their F(w) -> F(-w).

They propose that the quantum-mechanical origin is a system with a
spectrum BOUNDED ABOVE AND BELOW, evaluated on a Euclidean FOLD.  With
beta_dS = 2pi, continue a standard thermal 4-pt correlator to

    (tau1, tau2, tau3, tau4)^QM
        = (-tau1, 2 pi n - tau2, 2 pi - tau3, 2 pi (n+1) - tau4)   (6.9)

with n an integer, or use the negative-temperature variant (6.12).  That
flips the sign of their w-variable, and so turns scrambling into
anti-scrambling.

They pose the test directly: "It would be interesting to study the
validity of this assumption in a concrete model such as the SYK model."
Finite-N SYK has a bounded spectrum, so every folded correlator is finite
and exactly computable by ED.

WHAT THIS MODULE HOLDS
----------------------
The closed forms needed to run both tests against the repo's ED
laboratory (syk-self-averaging/fold_bench.py).

* Streicher 1911.10171, large-q Majorana-SYK 4-pt functions, uniformly
  valid in beta*scriptJ.  That covers all temperatures, including the
  T_B = infinity hologram point.  The TOC is his (3.3) and the crossed
  OTOC his (3.10), with f, g, the squeezed coordinates, and the
  prefactor C(v).  Transcribed from the arXiv LaTeX source at equation
  level; provenance is below.
* The T = infinity limits.  v -> beta*scriptJ/pi, and the Lyapunov rate
  lambda_L = 2 pi v / beta -> 2 scriptJ = 2 sqrt(2) p in repo units.  The
  growth survives at infinite temperature, and its sign is the SCRAMBLING
  sign, c > 0, for every v in (0,1).  Asserted in the tests.
* The Harlow-Zhao fold map (6.9)/(6.12) in theta units, beta = 2 pi.

CONVENTIONS
-----------
theta = 2 pi tau / beta in (0, 2 pi).  Flavor-averaged bilocals are
G_hat(th_i, th_j) = (1/N) sum_k psi_k(th_i) psi_k(th_j).
F(th1..th4) = <G_hat(th1,th2) G_hat(th3,th4)> with LITERAL operator
order.  F = F_d + script-F/N, with F_d = G(th12m) G(th34m).

v solves beta scriptJ = pi v / cos(pi v / 2), in largeq_anchor.v_of_betaJ.
The repo bridge is scriptJ = sqrt(2) p, locked against ED in
test_largeq.py.
"""

from __future__ import annotations

import numpy as np

from largeq_anchor import v_of_betaJ


# ---------------------------------------------------------------------------
# Streicher closed forms (LaTeX-source transcription; complex-safe)
# ---------------------------------------------------------------------------
def prefactor_C(v: float) -> float:
    """C(v) = 2 tan(pi v/2) / (pi v/2 + cot(pi v/2))  [prefactor of (3.3),
    (3.10)].  Small v: C ~ (pi v)^2 / 2."""
    return 2.0 * np.tan(np.pi * v / 2.0) / (
        np.pi * v / 2.0 + 1.0 / np.tan(np.pi * v / 2.0))


def f_toc(theta, v: float):
    """f(theta) = 1 - (v theta/2 + cot(pi v/2)) tan(v (pi - theta)/2)
    [Streicher (3.2)]; complex-safe."""
    theta = np.asarray(theta, dtype=complex)
    return 1.0 - (v * theta / 2.0 + 1.0 / np.tan(np.pi * v / 2.0)) \
        * np.tan(v * (np.pi - theta) / 2.0)


def squeeze(theta, v: float):
    """phi with pi - phi = v (pi - theta)  [squeezed coordinates]."""
    theta = np.asarray(theta, dtype=complex)
    return np.pi - v * (np.pi - theta)


def g_otoc(phi):
    """g(phi) = 1 + (pi - phi) / (2 tan(phi/2))  [Streicher (3.9)]."""
    phi = np.asarray(phi, dtype=complex)
    return 1.0 + (np.pi - phi) / (2.0 * np.tan(phi / 2.0))


def toc_ratio(th1, th2, th3, th4, v: float):
    """script-F/F_d for the UNCROSSED (time-ordered) configuration
    theta1 > theta2 > theta3 (> theta4)  [Streicher (3.3), first branch].
    Depends only on th12m = th1-th2 and th34m = th3-th4."""
    return prefactor_C(v) * f_toc(th1 - th2, v) * f_toc(th3 - th4, v)


def otoc_ratio(th1, th2, th3, th4, v: float):
    """script-F/F_d for the CROSSED configuration theta1 > theta3 > theta2
    [Streicher (3.5)], as an analytic function of the four thetas --
    complex arguments implement Lorentzian continuation
    (tau_j -> tau_j + i t_j  <=>  theta_j -> theta_j + 2 pi i t_j / beta):

        C(v) g(ph12m) g(ph34m)
        - 2 sin(ph12p/2) / (cos(pi v/2) sin(ph12m/2) sin(ph34m/2))
        - tan(pi v/2) (pi - ph12p) / (tan(ph12m/2) tan(ph34m/2)).

    GAUGE (the 2026-08-24 review's high-severity correction).  Streicher
    writes that display after fixing the time-translation freedom:
    Appendix A, "we use time-translation invariance to set theta_34^+ = 0"
    (with theta_3 = -theta_4 = theta_34^-/2).  So the phi_12^+ appearing in
    (3.5) is the squeeze of theta_12^+ - theta_34^+, not of theta_12^+.
    The two agree in Streicher's gauge and nowhere else.  Using the
    gauge-fixed display at general angles breaks time-translation
    invariance, which script-F/F_d must have -- the chord engine's
    correlator depends only on time DIFFERENCES (chord_fold._segments), as
    does ED (fold_bench).  Evidence for this form, all in
    test_antiscrambling.py: it is translation invariant to machine
    precision (the old form drifts by a factor ~2 over a 0.5 shift); it
    satisfies OTOC - TOC = -2 exactly at the coincidence theta_2 = theta_3
    (the old form gives arbitrary values); and the exact N = infinity chord
    engine converges onto it as lambda -> 0, where the double-scaled theory
    IS the large-q limit (chord/closed 0.964 -> 0.979 over lambda = 2 ->
    0.5, the committed chord_fold.json range, tracking the TOC control;
    against the old form the same ratios run to 1.68, away from 1).
    """
    ph12m = squeeze(np.asarray(th1) - th2, v)
    ph34m = squeeze(np.asarray(th3) - th4, v)
    ph12p = squeeze(np.asarray(th1) + th2 - np.asarray(th3) - th4, v)
    return (prefactor_C(v) * g_otoc(ph12m) * g_otoc(ph34m)
            - 2.0 * np.sin(ph12p / 2.0)
            / (np.cos(np.pi * v / 2.0) * np.sin(ph12m / 2.0)
               * np.sin(ph34m / 2.0))
            - np.tan(np.pi * v / 2.0) * (np.pi - ph12p)
            / (np.tan(ph12m / 2.0) * np.tan(ph34m / 2.0)))


def G2_thermal(theta, betaJ: float, q: float):
    """Large-q Euclidean 2-pt function on the thermal circle (beta = 2 pi in
    theta units), single-valued on the principal branch of the 2/q power --
    with a branch-point divergence ~ (theta_pole - theta)^(-2/q) at
    theta = pi + pi/v (the location test_g2_stays_physical measures by
    bisection) -- valid at the complex/folded arguments the HZ map
    produces:

        G(theta) = [ cos(pi v/2) / cos(pi v (1/2 - theta/(2 pi))) ]^{2/q}.

    At real theta in (0, 2 pi) this equals g_thermal_euclidean(...)^{1/q}
    (asserted in tests, which keeps largeq_anchor's e^{g} exercised)."""
    v = v_of_betaJ(betaJ)
    theta = np.asarray(theta, dtype=complex)
    base = np.cos(np.pi * v / 2.0) / np.cos(np.pi * v * (0.5 - theta / (2.0 * np.pi)))
    return base ** (2.0 / q)


# ---------------------------------------------------------------------------
# the sign test (CK Eq. (2.74)): scrambling or anti-scrambling?
# ---------------------------------------------------------------------------
def lyapunov_rate(betaJ: float) -> float:
    """lambda_L / scriptJ = 2 pi v / (beta scriptJ).  As betaJ -> 0
    (T_B = infinity) this -> 2: the growth persists at infinite temperature
    with rate 2 scriptJ (= 2 sqrt(2) p in repo units)."""
    if betaJ == 0.0:
        return 2.0
    return 2.0 * np.pi * v_of_betaJ(betaJ) / betaJ


def growth_coefficient(v: float) -> complex:
    """Closed-form coefficient of the exponentially growing OTOC term at the
    symmetric configuration theta = (3pi/2, pi/2, pi, 0), Lorentzian
    theta-time thL added to the (theta1, theta2) flavor pair:

        script-F/F_d  ~  a(v) e^{v thL} + (decaying),
        a(v) = - 1 / cos(pi v / 2)          (REAL and negative)

    Derivation: at this configuration phi12m = phi34m = pi, so the third
    term of otoc_ratio vanishes identically and the first is the constant
    C(v).  The gauge-invariant phi12p = squeeze(th12+ - th34+) equals
    pi + 2 i v thL here, so sin(phi12p/2) = cosh(v thL) and the second term
    is exactly -2 cosh(v thL)/cos(pi v/2): a growing/decaying pair with a
    REAL coefficient, and no linear-in-thL term.

    a(v) < 0 at EVERY temperature: the SCRAMBLING sign (CK's c > 0, the
    AdS/time-delay sign).  |a| = 1/cos(pi v/2) GROWS with beta (1.030 at
    betaJ = 0.5, 1.353 at betaJ = 2, 3.194 at betaJ = 8) and tends to
    a -> -1 at the T_B = infinity hologram point (v -> 0), where the growth
    is not suppressed at all.  A dS observer needs the opposite sign (CK's
    c < 0, time-advance): under the flat (boundary time = observer time)
    dictionary DSSYK_infinity cannot anti-scramble.

    NOTE (2026-08-24): the former value a = -e^{-i pi v/2}/cos(pi v/2),
    quoted as "Re a = -1 exactly at every temperature", was an artifact of
    the phi12p gauge error documented in otoc_ratio.  The sign verdict --
    the load-bearing one -- survives and strengthens; the "exactly -1"
    number holds only in the v -> 0 limit.  The antipodal rotation law
    a(Delta) = a(0) e^{-i v Delta} (signflip_probe.py) is UNCHANGED, since
    Delta enters phi12p identically in both forms; with a(0) now real it
    reads Re a(Delta) = a(0) cos(v Delta), which flips sign at
    Delta = pi/v exactly.  Cross-checked numerically from otoc_ratio and
    against ED in test_antiscrambling.py."""
    return complex(-1.0 / np.cos(np.pi * v / 2.0))


def growth_coefficient_numeric(betaJ: float) -> float:
    """a(v) extracted numerically from otoc_ratio.  At the symmetric config
    (structure derived in growth_coefficient above)
    dev.real = C(v) - (e^{v thL} + e^{-v thL})/cos(pi v/2): fit the basis
    [1, e^{v thL}, e^{-v thL}] over a late window and return the
    e^{v thL} coefficient (= a(v), machine precision); validates
    growth_coefficient independently of its closed form."""
    v = v_of_betaJ(betaJ) if betaJ > 0 else 1e-4
    thLs = np.array([4.0, 5.0, 6.0, 7.0]) / v
    devs = np.array([complex(otoc_ratio(1.5 * np.pi + 1j * t,
                                        0.5 * np.pi + 1j * t,
                                        np.pi, 0.0, v)) for t in thLs])
    A = np.vstack([np.ones_like(thLs), np.exp(v * thLs), np.exp(-v * thLs)]).T
    coef, *_ = np.linalg.lstsq(A, devs.real, rcond=None)
    return float(coef[1])


# ---------------------------------------------------------------------------
# the Harlow-Zhao Euclidean fold (2607.14215 Eqs. (6.9), (6.12))
# ---------------------------------------------------------------------------
def fold_map(taus, n: int = 1):
    """HZ (6.9): (tau1..tau4)^QM = (-t1, 2 pi n - t2, 2 pi - t3,
    2 pi (n+1) - t4), beta = 2 pi units.  Flips the sign of the
    (tau43 + tau21) combination (their w -> -w) while preserving tau31,
    tau42 modulo the KMS reflection -- scrambling -> anti-scrambling."""
    t1, t2, t3, t4 = taus
    return (-t1,
            2.0 * np.pi * n - t2,
            2.0 * np.pi - t3,
            2.0 * np.pi * (n + 1) - t4)


def fold_map_negT(taus, n: int = 1):
    """HZ (6.12), the negative-temperature variant:
    (tau1, tau2 - 2 pi n, tau3 - 2 pi, tau4 - 2 pi (n+1)), to be used with
    beta -> -beta, H -> -H."""
    t1, t2, t3, t4 = taus
    return (t1,
            t2 - 2.0 * np.pi * n,
            t3 - 2.0 * np.pi,
            t4 - 2.0 * np.pi * (n + 1))


PROVENANCE = """
Equation-verified sources (arXiv, fetched 2026-07-23):
  1911.10171 (Streicher) LaTeX source: f (Sec. 3.1 display), (3.3) TOC
    branches, squeezed coordinates + g (Sec. 3.2 displays), (3.10) crossed
    OTOC, F = F_d + script-F/N (Sec. "Four-Point Correlators..."), (3.11)
    growth term, (3.14) lambda_L = 2 pi v / beta.
  2607.13665 (Cui-Kolchmeyer): (2.74) eikonal OTOC with x = c e^T, (2.77)
    "leading correction can only decrease the magnitude", Sec. 2.2.3
    scrambling-vs-anti-scrambling sign dictionary; Sec. 5 KMS/trace
    conclusions.
  2607.14215 (Harlow-Zhao): (6.4) large-w expansion sign, (6.8) 2-pt KMS
    invariance under the fold, (6.9)/(6.12) fold prescriptions, (6.10) the
    sign flip mechanism.
ED ground truth: syk-self-averaging/fold_bench.py (fold_bench.json),
asserted in test_antiscrambling.py.
"""

if __name__ == "__main__":
    print("Large-q SYK OTOC sign test (CK 2607.13665 Eq. 2.74 language):\n")
    print(f"{'betaJ':>7} {'v':>7} {'lambda_L/J':>11} {'Re a exact':>11} "
          f"{'Re a fit':>11}")
    for bJ in (8.0, 4.0, 2.0, 1.0, 0.5, 0.1):
        v = v_of_betaJ(bJ)
        print(f"{bJ:>7.2f} {v:>7.4f} {lyapunov_rate(bJ):>11.4f} "
              f"{growth_coefficient(v).real:>11.6f} "
              f"{growth_coefficient_numeric(bJ):>11.6f}")
    print("\na < 0 at every temperature INCLUDING T_B -> infinity:")
    print("large-q SYK scrambles (CK's c > 0, the AdS/time-delay sign).")
    print("A dS observer needs c < 0.  Under the flat dictionary DSSYK_inf")
    print("cannot anti-scramble; the fold/tomperature map has to do it.")
    print("\nHZ fold map (6.9), n=1, of (3pi/2, pi/2, pi, 0):")
    print("  ->", tuple(round(float(x), 4) for x in
                        fold_map((1.5 * np.pi, 0.5 * np.pi, np.pi, 0.0))))
    print(PROVENANCE)
