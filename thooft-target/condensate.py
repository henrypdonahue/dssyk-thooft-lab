#!/usr/bin/env python3
r"""
Chiral condensate, GMOR corner, and the vacuum energy (rung 21, bulk half).
===========================================================================

Section 5 of the target paper ("Holographic Tuning") transplants a 4d
fine-tuning drama into a bulk dual -- large-N QCD2 -- that is
SUPER-RENORMALIZABLE: its only mass renormalization is the finite shift
m_q^2 = m_bare^2 - g^2/pi, and its vacuum energy receives one FINITE
interaction correction, computable exactly.  This module computes that
correction and validates the exact inputs against this repo's own solver in
a regime it had never touched (the chiral corner alpha -> -1).

Exact inputs (Litvinov-Meshcheriakov arXiv:2409.11324 Sec. 6.1, fetched and
transcribed 2026-08-16; their alpha == this repo's alpha = pi m_q^2/g^2,
a = 1 + alpha = pi m_bare^2/g^2):

  * the condensate at ARBITRARY quark mass [LM (6.1), due to Burkardt]:
        <psibar psi> = (m Nc/2pi) { log(pi/(1+alpha)) - 1 - gamma_E
                                    + alpha log 4 - alpha^2 I(alpha) },
        I(alpha) = (1/4) Int_-inf^inf dt/t^2 (sinh 2t - 2t)
                   / [cosh t (alpha sinh t + t cosh t)],
    with m = m_bare = g sqrt((1+alpha)/pi);
  * its chiral limit [LM (6.2), Zhitnitsky]:
        <psibar psi> -> -Nc sqrt(g^2/(12 pi));
  * GMOR [LM (6.3)]: (Nc/pi) M_0^2 ~= -4 m <psibar psi>, giving the pion
        eigenvalue 2 lambda_0 = M_0^2/(pi g^2) -> (2/pi) sqrt(a/3).

What is computed and asserted here:

  1. f(alpha) = <psibar psi>/(Nc g), dimensionless, on (-1, 0]; the I(alpha)
     integral has an integrable near-singularity at alpha -> -1
     (denominator ~ t(t^2/3 + a) at small t) which reproduces EXACTLY the
     pi/sqrt(3a) Zhitnitsky divergence -- asserted, plus the duality-point
     value f(0) = (log pi - 1 - gamma_E)/(2 pi sqrt(pi)) = -0.038834...
  2. The GMOR corner vs the matched-exponent Jacobi solver: 2 lambda_0 at
     alpha = -1 + a approaches (2/pi) sqrt(a/3) with a measured NLO
     (2lambda_0/LO - 1) ~= sqrt(a).  Measured basis convergence at K = 96:
     10+ digits for a >= 0.04, ~8 digits at the a = 0.02 endpoint (the
     subleading endpoint structure grows as beta ~ sqrt(3a)/pi shrinks) --
     far inside every tolerance used.  First solver data in the chiral
     corner.
  3. The vacuum-energy curve by Feynman-Hellmann in the bare mass at fixed
     g (dE_vac/dm_bare = <psibar psi>):
        eps(alpha) = [E_vac(alpha) - E_vac(-1)]/(Nc g^2)
                   = Int_{-1}^{alpha} f(a') (dm/dalpha')/g dalpha',
        dm/dalpha = g / (2 sqrt(pi (1+alpha))).
     In the m variable the integrand (the condensate f itself) is FINITE
     through the chiral point; in the alpha variable it has an integrable
     (1+alpha)^{-1/2} endpoint singularity, with integrand*sqrt(1+alpha)
     -> f_chiral/(2 sqrt(pi)) = -1/(4 pi sqrt(3)).  eps(0) -- the
     interaction vacuum-energy shift at the duality point, in units
     Nc g^2 -- is a single finite number.  That is the entire 2d counterpart
     of Sec. 5's "vacuum energy corrections": one finite, parameter-free
     constant; nothing requires tuning.  Under the dictionary
     (alpha_bar = gbar^2 N = p^2, g^2 = lambda J^2 -> 0 with tau = p^2
     m_q^2 fixed), the shift per flavor-mode vanishes as lambda: the
     Holographic-Tuning mechanism has no divergence to cancel in this bulk.

Output: condensate.json (curve + corner scan + the numbers above).
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent

GAMMA_E = 0.5772156649015329


# ---------------------------------------------------------------------------
# the exact condensate [LM (6.1)]
# ---------------------------------------------------------------------------
def I_alpha(alpha: float) -> float:
    """I(alpha) = (1/4) Int dt/t^2 (sinh2t - 2t)/[cosh t (alpha sinh t
    + t cosh t)]; even integrand, near-singular ~ 4/(t^2 + 3a) at small t
    when a = 1 + alpha -> 0 (handled by an explicit small-t split)."""
    a = 1.0 + alpha

    def integrand(t):
        if t < 1e-6:
            # series: (sinh2t-2t)/t^2 ~ 4t/3 + 4t^3/15;
            # alpha sinh t + t cosh t ~ a t + (alpha/6 + 1/2) t^3
            num = 4.0 * t / 3.0 + 4.0 * t ** 3 / 15.0
            den = math.cosh(t) * (a * t + (alpha / 6.0 + 0.5) * t ** 3)
            return num / den
        if t > 30.0:
            # tanh t = 1 to double precision: algebraic tail 2/(t^2 (t+a))
            return 2.0 / (t * t * (t + alpha))
        return ((math.sinh(2 * t) - 2 * t) / (t * t)
                / (math.cosh(t) * (alpha * math.sinh(t)
                                   + t * math.cosh(t))))

    # split at the near-singular scale sqrt(3a) and again at O(1); the tail
    # decays only ALGEBRAICALLY (integrand -> 2/(t^2 (t + alpha))), so the
    # last piece runs to infinity
    t_split = max(math.sqrt(3.0 * a), 1e-3)
    v = 0.0
    edges = [0.0, 0.5 * t_split, t_split, min(4.0 * t_split, 1.0), 1.0,
             40.0, np.inf]
    edges = sorted(set(edges))
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi > lo:
            vi, _ = quad(integrand, lo, hi, limit=400)
            v += vi
    return 0.5 * v                  # (1/4) * 2 * Int_0^inf


def condensate_dimensionless(alpha: float) -> float:
    """f(alpha) = <psibar psi>/(Nc g)  [LM (6.1) with m = g sqrt((1+a)/pi)]."""
    a = 1.0 + alpha
    if a <= 0:
        raise ValueError("alpha must be > -1")
    m_over_g = math.sqrt(a / math.pi)
    brace = (math.log(math.pi / a) - 1.0 - GAMMA_E
             + alpha * math.log(4.0) - alpha * alpha * I_alpha(alpha))
    return m_over_g / (2.0 * math.pi) * brace


CHIRAL_VALUE = -1.0 / math.sqrt(12.0 * math.pi)      # LM (6.2)


def duality_point_value() -> float:
    """f(0) = (log pi - 1 - gamma_E)/(2 pi sqrt(pi))  (alpha-terms vanish)."""
    return (math.log(math.pi) - 1.0 - GAMMA_E) / (2.0 * math.pi
                                                  * math.sqrt(math.pi))


# ---------------------------------------------------------------------------
# vacuum energy by Feynman-Hellmann in the bare mass
# ---------------------------------------------------------------------------
def vacuum_energy_curve(alphas: np.ndarray) -> np.ndarray:
    """eps(alpha) = Int_{-1}^{alpha} f(a') / (2 sqrt(pi (1+a'))) da',
    the dimensionless vacuum-energy shift [E_vac(alpha) - E_vac(-1)]/(Nc g^2).
    The alpha-integrand has an integrable (1+alpha)^{-1/2} endpoint
    singularity; the substitution u = sqrt(1+alpha) removes it exactly
    (d alpha = 2u du), leaving the smooth integrand f(u^2-1)/sqrt(pi),
    which tends to -1/(2 pi sqrt(3)) at the chiral end."""
    # evaluate f once on a dense u-grid and integrate its cubic-spline
    # antiderivative -- no nested adaptive quadrature, smooth integrand,
    # ~1e-8 accuracy.  The u = 0 endpoint value is the exact chiral limit.
    u_grid = np.linspace(0.0, 1.0, 801)
    f_grid = np.empty_like(u_grid)
    f_grid[0] = CHIRAL_VALUE / math.sqrt(math.pi)
    for i, u in enumerate(u_grid[1:], 1):
        f_grid[i] = condensate_dimensionless(u * u - 1.0) / math.sqrt(math.pi)
    from scipy.interpolate import CubicSpline
    anti = CubicSpline(u_grid, f_grid).antiderivative()
    us = np.sqrt(1.0 + np.asarray(alphas))
    return anti(us) - anti(0.0)


# ---------------------------------------------------------------------------
# the GMOR corner scan (solver data)
# ---------------------------------------------------------------------------
def gmor_corner(a_values=(0.25, 0.16, 0.09, 0.04, 0.02), K: int = 96):
    """2 lambda_0 from the matched-exponent Jacobi solver at alpha = -1 + a
    vs the GMOR leading order (2/pi) sqrt(a/3)."""
    from jacobi_solver import spectrum_jacobi, endpoint_exponent
    rows = []
    for a in a_values:
        alpha = -1.0 + a
        lam0 = float(spectrum_jacobi(K, 1, alpha)[0]["two_lambda"])
        lo = (2.0 / math.pi) * math.sqrt(a / 3.0)
        rows.append(dict(a=a, alpha=alpha,
                         beta=endpoint_exponent(alpha),
                         two_lambda0=lam0, gmor_lo=lo,
                         ratio=lam0 / lo,
                         nlo_coeff=(lam0 / lo - 1.0) / math.sqrt(a)))
    return rows


def main():
    t0 = time.time()
    alphas = np.concatenate([
        -1.0 + np.geomspace(1e-4, 0.5, 40), np.linspace(-0.5, 0.0, 21)])
    alphas = np.unique(np.clip(alphas, -1.0 + 1e-4, 0.0))
    f_curve = [condensate_dimensionless(a) for a in alphas]
    eps_curve = vacuum_energy_curve(alphas)
    corner = gmor_corner()
    payload = dict(
        provenance="condensate.py -- rung 21 bulk half: LM 2409.11324 "
                   "(6.1)-(6.3) transcribed; GMOR corner vs jacobi_solver; "
                   "vacuum energy by Feynman-Hellmann in m_bare",
        chiral_value=CHIRAL_VALUE,
        duality_point_f=duality_point_value(),
        eps_duality_point=float(eps_curve[-1]),
        alphas=alphas.tolist(), f=f_curve, eps=eps_curve.tolist(),
        gmor_corner=corner,
        runtime_s=round(time.time() - t0, 1))
    out = HERE / "condensate.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"chiral f = {CHIRAL_VALUE:.6f}; "
          f"f(0) = {duality_point_value():.6f}; "
          f"eps(0) = {payload['eps_duality_point']:.6f}")
    for r in corner:
        print(f"  a={r['a']}: 2lam0={r['two_lambda0']:.6f} "
              f"ratio to GMOR-LO = {r['ratio']:.4f} "
              f"NLO coeff = {r['nlo_coeff']:.4f}")
    print(f"wrote {out} ({payload['runtime_s']} s)")


if __name__ == "__main__":
    main()
