#!/usr/bin/env python3
r"""
The forward tower check: 't Hooft masses through the boost kernel.
==================================================================

BOOST.md names three observables that could discriminate the headline
test.  Its item (b) is "moment ratios against the kernel-transformed
't Hooft tower".  This module RUNS that forward direction.

No inversion, no new theory.  Sum the single-mass boost spectral function
of boost_mass.py over the certified tower, with decay-constant couplings.
Take moments.  Compare the scale-free invariants against the repo's exact
SYK-side values.

TOWER AND COUPLINGS
-------------------
S_tower(omega) = sum_n c_n S(omega; M_n, a),
S(omega; m, a) = (1/(pi a)) e^{pi omega/a} K_{i omega/a}(m/a)^2.

M_n = sqrt(two_lambda_n) in units sqrt(pi g^2) = 1, from thooft-target.

The coupling ansatz is c_n = f_n^2 * (M_n^2)^s.  f_n are the measured
decay constants (eigenfunctions.decay_constants), and s in {0, 1}
parametrizes the operator-normalization ambiguity, field-like against
current-like.

The f_n^2 weighting projects onto the CP = -1 (symmetric-FLZ) tower
automatically: antisymmetric levels have f = 0 by parity, which matches
the CP = (-1)^(n+1) assignment of the even-n_paper channels
(dictionary.py).  The CP = +1 channel would need the x-overlap couplings
instead.  That is a second ansatz axis, and it is not implemented.
Stated, not hidden.

THE STRUCTURAL FACT THIS ESTABLISHES FIRST
------------------------------------------
RAW BOOST MOMENTS DIVERGE.

The worldline 2-pt function of a free field behaves as
W ~ -(1/2pi) ln|tau| at short times.  So S(omega) -> 1/omega at large
omega, independent of the mass.  Measured: the mean of omega*S over
[60, 80] is 0.994 at m = a = 1.

Hence mu_0 diverges logarithmically, and every higher moment diverges as
a power.  So all moments here carry the standard Wightman regulator
e^{-eps omega}, the same eps as boost_mass.spectral_from_wightman.  Every
invariant is then a function of (eps, a) -- and with the masses fixed by
the yardstick, of the two dimensionless combinations (m/a, a*eps).

The regulator is not a numerical convenience.  It stands in for the
operator smearing that the missing tomperature/state map must supply.  So
the eps-dependence of the invariants measures how much of the "spectral
shape" is dictionary rather than physics.

EXACT MOMENT ROUTE (no quadrature)
----------------------------------
mu_k(eps) = int domega omega^k S(omega) e^{-eps omega}
          = 2 pi (-1)^k  d^k/ds^k [ (1/2pi) K_0( (2m/a) sin(a s/2) ) ]
            at s = eps                                  (0 < a eps < 2 pi)

These are derivatives of the SAME Euclidean-section Wightman function
that boost_mass.wightman_tau continues.

For even k they equal the even moments of the symmetrized spectral
function, up to the regulator asymmetry (e^{-eps omega} against a
symmetric cutoff).  That asymmetry is O(eps), and it belongs to the same
ambiguity budget.  The KMS mirror weight e^{-(2pi/a - eps) omega} on
omega < 0 is what makes the negative-frequency side converge.

The derivatives are computed by exact Taylor composition.  K_0
derivatives come from the Bessel ODE recurrence
K_0^{(j+2)} = K_0^{(j)} - d^j(K_0'/z), composed with the sine worldline,
in mpmath scalars.  That is two Bessel evaluations per (mass, eps) in
total.  It is verified against brute-force quadrature of
spectral_closed_form to machine precision (test_boost_forward.py).

THE COMPARISON, AND THE VERDICT
-------------------------------
Numbers are in boost_forward.json.

SYK targets, exact, from this repo.  The chord O_Delta channel kurtosis
kurt = mu_4 mu_0 / mu_2^2 spans 2.53 -> 4.81 over the committed lambda
grid (the Delta = 1 and 2/p closed forms, chord_moments.json).  The full
A_2 channel at finite N along lambda_chord = 2 gives 4.59 -> 4.11 over
(p, N) = (4,16)..(10,100), from the exact_pair_moments machinery.

Tower result.  kurt(mscale, eps) sweeps a range that COVERS the entire
SYK band, and it diverges logarithmically as eps -> 0, which is the
regulator-dominated regime.

So the boost kernel does not structurally exclude the match.  But the
match point carries no discriminating power by itself, because eps is a
free dial that can land kurt on any target.  The discriminating content
is sharper, and it is what the JSON records:

  1. The eps* required to hit each SYK target at fixed a.  That is the
     SPEC the tomperature/state map must reproduce.  One map has to give
     one eps*, not one per channel.
  2. The difference between the full tower and the ground state alone, at
     the invariant level, is at most about 5% across the grid.  Boltzmann
     suppression plus the f_n^2 falloff put 89% of the coupling weight on
     the ground state.  So even with the map fixed, moment invariants at
     these eps see essentially ONE mass.  Node structure, not moments, is
     where the tower lives.  That is BOOST.md fact 2, now quantified at
     the invariant level.

kurt(eps) is not monotone.  It turns back up past a*eps ~ 1.6, where the
KMS mirror weight starts to dominate.  So eps* is a bracketed root, not
an inverse function, and the grid records the full shape.

Scope: the free-particle kernel (BOOST.md's calibration scope), CP = -1
couplings, and symmetric-regulator convention differences at O(eps) left
unresolved.  This bounds what moment data can show.  It does not fix the
map.

Output: boost_forward.json.  Runtime: seconds.
"""

from __future__ import annotations

import json
import sys
import time
from math import comb, factorial
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "thooft-target"))

from eigenfunctions import decay_constants  # noqa: E402

ORDERS = (0, 2, 4, 6)
_SERIES_LEN = 7          # Taylor order 0..6


# ---------------------------------------------------------------------------
# exact Taylor derivatives of W(s) = K_0(z(s))/2pi, z = (2m/a) sin(a s/2)
# ---------------------------------------------------------------------------
def _k0_derivs(z0, n: int):
    """K_0^{(j)}(z0), j = 0..n, from the Bessel ODE recurrence
    K_0'' = K_0 - K_0'/z (differentiated with Leibniz on K_0'/z).
    mpmath scalars; exactly two Bessel evaluations."""
    c = [mp.besselk(0, z0), -mp.besselk(1, z0)]
    for j in range(n - 1):
        # c_{j+2} = c_j - sum_i C(j,i) c_{i+1} (-1)^{j-i} (j-i)! z0^{-(j-i+1)}
        acc = mp.mpf(0)
        for i in range(j + 1):
            acc += (comb(j, i) * c[i + 1] * (-1) ** (j - i)
                    * factorial(j - i) / z0 ** (j - i + 1))
        c.append(c[j] - acc)
    return c


def _series_mul(a, b):
    out = [mp.mpf(0)] * _SERIES_LEN
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if i + j < _SERIES_LEN:
                out[i + j] += ai * bj
    return out


def w_series(eps: float, m: float, a: float = 1.0):
    """Taylor coefficients w_k of W(eps + t) = K_0((2m/a) sin(a(eps+t)/2))/2pi
    to order 6: exact composition of the K_0 series (ODE recurrence) with the
    sine worldline series.  Runs at 40 working digits: the recurrence
    coefficients grow like j!/z0^j at small z0 = (2m/a) sin(a eps/2), so
    the default 15 digits would shed accuracy exactly in the small-eps
    corner the divergence scan probes."""
    with mp.workdps(40):
        return _w_series_impl(eps, m, a)


def _w_series_impl(eps: float, m: float, a: float):
    eps_, m_, a_ = mp.mpf(eps), mp.mpf(m), mp.mpf(a)
    if not 0 < a_ * eps_ < 2 * mp.pi:
        raise ValueError("need 0 < a*eps < 2*pi for convergence")
    # z(eps + t) series: (2m/a) sin(a eps/2 + a t/2)
    half = a_ * eps_ / 2
    s0, c0 = mp.sin(half), mp.cos(half)
    z = []
    for k in range(_SERIES_LEN):
        trig = (s0, c0, -s0, -c0)[k % 4]
        z.append((2 * m_ / a_) * trig * (a_ / 2) ** k / factorial(k))
    z0, dz = z[0], [mp.mpf(0)] + z[1:]
    cder = _k0_derivs(z0, _SERIES_LEN - 1)
    # Horner composition of P(w) = sum_j cder_j w^j / j! with w = dz (w_0 = 0)
    res = [mp.mpf(0)] * _SERIES_LEN
    res[0] = cder[-1] / factorial(_SERIES_LEN - 1)
    for j in range(_SERIES_LEN - 2, -1, -1):
        res = _series_mul(res, dz)
        res[0] += cder[j] / factorial(j)
    return [r / (2 * mp.pi) for r in res]


def single_mass_moments(m: float, eps: float, a: float = 1.0) -> dict:
    """mu_k(eps) = int domega omega^k S(omega; m, a) e^{-eps omega} for
    k = 0, 2, 4, 6: exact via 2 pi (-1)^k k! [t^k] W-series."""
    ws = w_series(eps, m, a)
    return {k: float(2 * mp.pi * (-1) ** k * factorial(k) * ws[k])
            for k in ORDERS}


# ---------------------------------------------------------------------------
# the tower
# ---------------------------------------------------------------------------
def tower(num_levels: int = 12, num_basis: int = 200):
    """(masses, weights_f2, two_lambda, parseval_captured): the CP = -1
    (symmetric-FLZ) levels in units sqrt(pi g^2) = 1, decay-constant
    weights f_n^2 (Parseval-normalized: the full sum is exactly 1)."""
    f, ev = decay_constants(num_basis, num_levels=num_levels)
    masses = np.sqrt(ev)
    w = f ** 2
    return masses, w, ev, float(np.sum(w))


def tower_invariants(masses, weights, two_lambda, mscale: float, eps: float,
                     s_pow: int = 0) -> dict:
    """Scale-free invariants of the coupling-weighted tower at kernel scale
    m_n/a = mscale * M_n and regulator a*eps = eps (computed in a = 1
    units; only these two combinations matter).  c_n = f_n^2 (M_n^2)^s."""
    tot = {k: 0.0 for k in ORDERS}
    for m_i, w_i, ev_i in zip(masses, weights, two_lambda):
        c = w_i * (ev_i ** s_pow)
        mu = single_mass_moments(m_i * mscale, eps)
        for k in ORDERS:
            tot[k] += c * mu[k]
    return dict(
        mu0=tot[0], mu2=tot[2], mu4=tot[4], mu6=tot[6],
        kurt=tot[4] * tot[0] / tot[2] ** 2,
        r6=tot[6] * tot[0] ** 2 / tot[2] ** 3)


def eps_star_for_kurt(masses, weights, two_lambda, mscale: float,
                      target: float, s_pow: int = 0,
                      lo: float = 0.02, hi: float = 3.0) -> float:
    """A regulator eps* at which the tower kurtosis equals `target`
    (bracketed bisection on [lo, hi]; kurt is continuous but not monotone,
    so this returns one valid root; NaN if [lo, hi] does not bracket)."""
    def kurt(e):
        return tower_invariants(masses, weights, two_lambda,
                                mscale, e, s_pow)["kurt"]
    klo, khi = kurt(lo), kurt(hi)
    if not (khi <= target <= klo):
        return float("nan")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if kurt(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# SYK-side targets
# ---------------------------------------------------------------------------
def syk_targets() -> dict:
    """The exact SYK-side kurtosis values this comparison aims at: the chord
    O_Delta closed-form channels over the committed lambda grid, and the
    finite-N full A_2 channel along lambda_chord = 2 (values recomputed and
    pinned in the tests; frozen here for the campaign record)."""
    with open(HERE / "chord_moments.json") as fh:
        pts = json.load(fh)["points"]
    odelta = [dict(lam=row["lam"], channel=name,
                   kurt=row["channels"][name]["kurt"],
                   w2=row["channels"][name]["w2"])
              for row in pts
              for name in ("O_delta=2/p", "O_delta=1")]
    # full A_2 channel, exact finite N (recomputed live in the tests from
    # exact_pair_moments.channel_pair_moments; 7 s, so frozen here)
    a2_full = [dict(p=4, N=16, kurt=4.5875), dict(p=6, N=36, kurt=4.3151),
               dict(p=8, N=64, kurt=4.1857), dict(p=10, N=100, kurt=4.1119)]
    kurts = [r["kurt"] for r in odelta] + [r["kurt"] for r in a2_full]
    return dict(odelta=odelta, a2_full=a2_full,
                band=[min(kurts), max(kurts)])


# ---------------------------------------------------------------------------
# campaign
# ---------------------------------------------------------------------------
MSCALES = (0.35, 0.7, 1.0, 1.5, 2.5)
EPSES = (0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.0)


def run(out_path=None) -> dict:
    t0 = time.time()
    masses, w, ev, captured = tower(12)
    grid = []
    for s_pow in (0, 1):
        for ms in MSCALES:
            for eps in EPSES:
                inv = tower_invariants(masses, w, ev, ms, eps, s_pow)
                grid.append(dict(s_pow=s_pow, mscale=ms, eps=eps, **inv))
    kurts = [g["kurt"] for g in grid]

    # the eps -> 0 regulator divergence (log signature: equal increments
    # per halving) and the tower-vs-single-mass discrimination
    div = [dict(eps=e, **tower_invariants(masses, w, ev, 1.0, e))
           for e in (0.4, 0.2, 0.1, 0.05, 0.025)]
    single = []
    for ms in MSCALES:
        for eps in (0.1, 0.4, 1.6):
            kt = tower_invariants(masses, w, ev, ms, eps)["kurt"]
            mu = single_mass_moments(masses[0] * ms, eps)
            k1 = mu[4] * mu[0] / mu[2] ** 2
            single.append(dict(mscale=ms, eps=eps, kurt_tower=kt,
                               kurt_ground_only=k1,
                               rel_diff=abs(kt / k1 - 1.0)))
    max_rel = max(r["rel_diff"] for r in single)

    # truncation control: 12 vs 8 levels at a central cell, plus the
    # conservative tail bound residual_weight * mu_k(m_last)
    inv12 = tower_invariants(masses, w, ev, 1.0, 0.4)
    m8, w8, ev8, _ = tower(8)
    inv8 = tower_invariants(m8, w8, ev8, 1.0, 0.4)
    trunc = dict(
        kurt_L12=inv12["kurt"], kurt_L8=inv8["kurt"],
        rel=abs(inv8["kurt"] / inv12["kurt"] - 1.0),
        residual_weight=1.0 - captured)

    targets = syk_targets()
    lo, hi = targets["band"]
    eps_star = []
    for label, target in (
            [("odelta_min", lo), ("odelta_delta1_lam0", 3.8573),
             ("a2_full_p10", 4.1119), ("odelta_max", hi)]):
        eps_star.append(dict(target_label=label, target_kurt=target,
                             eps_star_mscale1=eps_star_val(
                                 masses, w, ev, 1.0, target)))
    verdict = dict(
        structurally_excluded=False,
        tower_kurt_range=[min(kurts), max(kurts)],
        syk_band_covered=bool(min(kurts) < lo and max(kurts) > hi),
        regulator_dominated=True,
        max_tower_vs_single_mass_rel=max_rel,
        statement=(
            "The boost kernel does not exclude the match: the tower "
            "kurtosis sweeps the whole SYK band as the regulator eps "
            "varies.  But eps is exactly the freedom the missing "
            "tomperature/state map must fix (mu_0 diverges as eps -> 0), "
            "and at any fixed (eps, a) the tower's invariants differ from "
            f"the ground-state-only ones by at most {max_rel:.1%}: moment "
            "invariants see one mass.  The spec for the map is the "
            "eps_star table; the tower itself lives in the node "
            "structure, not the moments."))

    payload = dict(
        provenance="boost_forward.py -- the forward tower check of BOOST.md "
                   "item (b): 't Hooft masses + measured decay constants "
                   "through the exact boost kernel (Wightman-regulated "
                   "moments via K_0 Taylor composition, no quadrature); "
                   "SYK targets from chord_moments.json (O_Delta closed "
                   "forms) and the exact finite-N full A_2 channel",
        units="masses in sqrt(pi g^2) = 1; kernel depends on (mscale = 1/a, "
              "a*eps) only; coupling c_n = f_n^2 (M_n^2)^s_pow",
        tower=dict(two_lambda=[float(x) for x in ev],
                   masses=[float(x) for x in masses],
                   f_squared=[float(x) for x in w],
                   parseval_captured=captured),
        grid=grid, divergence_eps_to_0=div, single_mass_comparison=single,
        truncation=trunc, syk_targets=targets, eps_star=eps_star,
        verdict=verdict, runtime_s=round(time.time() - t0, 1))
    out = Path(out_path) if out_path else HERE / "boost_forward.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"tower kurt range over grid: {min(kurts):.3f} .. {max(kurts):.3f} "
          f"(SYK band {lo:.3f} .. {hi:.3f})")
    for r in eps_star:
        print(f"  eps* ({r['target_label']}, kurt = {r['target_kurt']:.3f}) "
              f"= {r['eps_star_mscale1']:.4f}")
    print(f"tower vs ground-state-only: max rel diff {max_rel:.2%}")
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


def eps_star_val(masses, w, ev, mscale, target):
    return eps_star_for_kurt(masses, w, ev, mscale, target)


if __name__ == "__main__":
    run()
