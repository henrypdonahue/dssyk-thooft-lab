#!/usr/bin/env python3
r"""
The inverse tower check: conditioning of moments -> spectrum on the boost kernel.
=================================================================================

road_map.md Part III leaves the moments -> spectrum direction open: "how
many moments a Regge trajectory needs" and "which finite-N object carries
the tower".  boost_forward.py ran the forward direction and found the
moment invariants regulator-dominated.  This module makes the inverse
direction QUANTITATIVE on the exactly solvable kernel of BOOST.md: given
regulated boost-kernel moments of a 't Hooft tower, how much of the tower
is recoverable, at what precision cost, and how does that compare with
node-level data?

Four bounded studies share one engine.  All statements are conditional on
the e^{-eps omega} regulator convention of boost_forward.py (the standard
Wightman i-eps): eps stands in for the smearing the missing
tomperature/state map must supply, so every number below records its
eps-dependence instead of hiding it.

ENGINE (study 1)
----------------
boost_forward's exact moment route, generalized to arbitrary order:

  mu_k(eps) = int domega omega^k S(omega; m, a) e^{-eps omega}
            = 2 pi (-1)^k d^k/ds^k [ K_0((2m/a) sin(a s/2)) / 2pi ]|_{s=eps}

computed to k <= 23 by the same K_0 ODE-recurrence Taylor composition
(reusing boost_forward._k0_derivs), now with an adaptive-precision
dps-doubling gate: the series is recomputed at doubled working precision
until successive runs agree coefficient-wise to rel_tol (default 1e-30;
the conditioning study runs at 1e-60).  The certified relative error is
returned with the coefficients.  Measured: at kmax = 23 the gate
converges at dps ~ 300 with certified error < 1e-140 across the whole
(m, eps) range used here, i.e. the recurrence has no catastrophic
cancellation in this range.  Validation (test_boost_inverse.py): against
trapezoid quadrature of boost_mass.spectral_closed_form with the same
regulator -- agreement at k = 20 is ~2e-15, quadrature-limited -- and
against boost_forward.single_mass_moments at k <= 6 to 1e-12.

WEIGHTS-FROM-MOMENTS CONDITIONING (study 2)
-------------------------------------------
The linear map V[k, n] = mu_k(M_n; eps) from level weights to regulated
moments, for the first L = 12 levels of the certified reference spectrum
(both CP sectors; masses in sqrt(pi g^2) = 1).  Rows are scaled to unit
2-norm, which is the "each moment measured to relative precision"
noise model; singular values are computed in mpmath at 220 digits from
entries certified to 1e-60, so every reported value is far above the
arithmetic floor.  Column-normalized ("shape") and raw-column (physical
weight) spectra agree to < 0.25 decades at K = 24 (< 0.5 in the
row-starved small-K cells): after row equilibration the residual
Boltzmann column-norm spread is < 0.05 decades, negligible against the
shape degeneracy.

Findings (the numbers are in boost_inverse.json, and the headline is
NEGATIVE in the expected direction):

* The singular values fall super-geometrically:
  log10(sigma_1/sigma_j) ~= a(eps) (j-1) + b(eps) (j-1)^2 with
  (a, b) = (1.06, 0.049) at a*eps = 0.2, (0.78, 0.057) at 0.4
  (rms 0.08), (0.44, 0.050) at 1.6.  Recovering level j costs 0.8
  digits of moment precision more than level j-1 at the start and 2.0
  by level 12 (a*eps = 0.4): the cost per level GROWS.  At a*eps = 0.1
  the decay is steeper still (a = 1.54) and the last singular values
  floor near 1e-17, so the quadratic fit degrades (rms 0.68).
* "How many moments does a 12-level trajectory need": 12.  The
  double-precision level count saturates at K = 12 at every eps on the
  grid (fine K-scan at a*eps = 0.4; a one-count threshold flicker at
  a*eps = 0.2).  Beyond K = L the singular values themselves move by
  < 0.5 decades at a*eps <= 0.4, 1.0 at 0.8, and 4.0 at 1.6: high
  moments stay informative only while k/eps is inside the mass window,
  and no number of further moments buys back the levels that precision
  has already lost.
* With 16-digit moment data: 9 recoverable weight directions at
  a*eps = 0.1, 10 at 0.2, 12 at 0.4 (the last marginal, amplification
  10^15.5), 12 at 1.6.  With 8-digit data: 6 / 6 / 7 / 10.  Smaller
  eps is WORSE -- it weights the mass-blind S -> 1/omega tail -- so
  the regulator dial the state map must fix also sets how many levels
  exist numerically.
* Restricting to even orders (the symmetrized correlator's observable
  set, 12 rows) is ~2 decades worse than all 24 orders; the CP = -1
  6-level sub-tower (the physical CP-resolved channel) is the easy
  case: all 6 levels above double precision, 5.8 decades total.

TWO-MASS SEPARABILITY (study 3)
-------------------------------
The nonlinear identifiability core: an equal-weight pair at
m0 (1 -+ delta/2) against the best single-mass fit (c, m*) matching
mu_0 and mu_1 exactly.  With 2 moments the fit is exact (2 parameters),
so THREE moments is the minimum that can distinguish at all, and the
third carries a relative residual measured to scale as

  r_2(delta) = C(eps) delta^2,  exponent 1.999 (fit over delta = 0.4
  .. 0.025), C = 0.008 (a*eps = 0.1), 0.046 (0.4), 0.090 (1.6).

So the resolvable splitting at moment precision eta is
delta_min = sqrt(eta / C): ~5e-8 at double precision, ~5e-4 at 8 digits
(a*eps = 0.4).  Higher moments do not help much: the residual profile
over k = 2..7 is flat within a factor ~1.4.  The fitted single mass is
biased low by 0.21 delta^2 (relative) -- an unresolved pair silently
drags the recovered mass below the mean.

NODE-LEVEL AGAINST MOMENT-LEVEL DATA (study 4)
----------------------------------------------
The roadmap verdict to test: "the inversion needs either node-level
data, or the state map fixed first."  Quantified here as digits of the
mass per digit of data.

* Node positions (zeros of K_{i omega}(m/a)) are EXACTLY independent of
  the e^{-eps omega} regulator: leverage 0 by construction.  Measured
  log-sensitivity d ln omega_j / d ln m = 0.54 / 0.44 / 0.40 for the
  first three nodes (omega = 2.963, 4.534, 5.880) at m/a = 1.  But
  S ~ K^2 has DOUBLE zeros, so locating a node from S-data of relative
  precision eta transfers as a square root: dm/m ~= 0.46 sqrt(eta).
* Moment invariants (kurt at fixed eps) transfer linearly,
  dm/m ~= 1.24 eta -- better than nodes for all eta < ~0.1 IF eps is
  known exactly.  But the measured regulator leverage
  |d ln kurt / d ln eps| / |d ln kurt / d ln m| = 0.93 .. 1.02 over
  a*eps = 0.2 .. 0.8: an unknown relative error u in the smearing scale
  costs the mass ~u, one for one, at ANY data precision.  So nodes win
  whenever u > ~0.46 sqrt(eta), and with the map entirely unfixed
  (u = O(1), the current situation) moments carry ZERO digits of mass
  while nodes carry half the data digits.
* For the coupling-weighted tower the ground node fills in
  (S_tower > 0 where S_ground = 0): the filled-node depth, 5.7e-2 of
  the local maximum, is a background-free, eps-independent measure of
  excited-level weight -- the ground contribution vanishes exactly
  there -- against the <= 5% shift the same tower leaves in moment
  invariants (boost_forward.json), which rides on a baseline the
  regulator sweeps by a factor ~10.

VERDICT.  The moments -> spectrum inversion on the solvable kernel is
exponentially ill-conditioned with a measured, accelerating exponent
(0.8 digits per level at the start, 2.0 by level 12, at the reference
regulator a*eps = 0.4; steeper at smaller eps), saturates at K = L
moments, and -- decisive for the roadmap question -- inherits a
one-for-one mass error from any uncertainty in the regulator scale that
the missing state map is supposed to fix.  Node-level data is free of
exactly that error.  This answers "how many moments a Regge trajectory
needs" on the solvable kernel: as many as there are levels, each ~1-2
digits more expensive than the last, and all of them hostage to the
state map; the roadmap's "the inversion needs node-level data" verdict
stands, now with numbers.

Scope: the free-particle kernel of BOOST.md (single particle, sharp
worldline).  Composite-operator smearing and interactions would dress
the kernel; the e^{-eps omega} regulator is a convention standing in
for the unfixed state map, and every eps-dependent number above is
recorded on an eps grid rather than at a preferred point.

Output: boost_inverse.json.  Runtime: ~1-3 min on 2 cores.
"""

from __future__ import annotations

import json
import sys
import time
from math import factorial
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "thooft-target"))

from boost_forward import _k0_derivs  # noqa: E402  (shared single source)
import dictionary  # noqa: E402

K_MAX = 23           # moment orders 0..K_MAX
SVD_DPS = 220        # working digits for singular values (entries at 1e-60)
COND_TOL = 1e-60     # entry certification for the conditioning study
L_LEVELS = 12


# ---------------------------------------------------------------------------
# study 1: the engine -- arbitrary-order regulated moments, certified
# ---------------------------------------------------------------------------
def _series_mul(a, b, length):
    out = [mp.mpf(0)] * length
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if i + j < length:
                out[i + j] += ai * bj
    return out


def _w_series_impl(eps, m, a, length):
    """Taylor coefficients of W(eps + t) = K_0((2m/a) sin(a(eps+t)/2))/2pi
    to order length-1, at the CURRENT working precision: boost_forward's
    exact composition (K_0 ODE recurrence o sine worldline), general order."""
    eps_, m_, a_ = mp.mpf(eps), mp.mpf(m), mp.mpf(a)
    if not 0 < a_ * eps_ < 2 * mp.pi:
        raise ValueError("need 0 < a*eps < 2*pi for convergence")
    half = a_ * eps_ / 2
    s0, c0 = mp.sin(half), mp.cos(half)
    z = []
    for k in range(length):
        trig = (s0, c0, -s0, -c0)[k % 4]
        z.append((2 * m_ / a_) * trig * (a_ / 2) ** k / factorial(k))
    z0, dz = z[0], [mp.mpf(0)] + z[1:]
    cder = _k0_derivs(z0, length - 1)
    res = [mp.mpf(0)] * length
    res[0] = cder[-1] / factorial(length - 1)
    for j in range(length - 2, -1, -1):
        res = _series_mul(res, dz, length)
        res[0] += cder[j] / factorial(j)
    return [r / (2 * mp.pi) for r in res]


def w_series_k(eps, m, a=1.0, kmax=K_MAX, rel_tol=1e-30, dps0=None,
               max_doublings=8):
    """Certified Taylor coefficients w_0..w_kmax of the regulated Wightman
    function, by a dps-doubling convergence gate: recompute at doubled
    working precision until successive runs agree coefficient-wise to
    rel_tol.  Returns (coeffs, certified_rel_err, dps_used)."""
    length = kmax + 1
    dps = dps0 or max(50, 40 + 2 * kmax + int(-mp.log10(rel_tol)))
    prev = None
    worst = mp.inf
    for _ in range(max_doublings + 1):
        with mp.workdps(dps):
            cur = _w_series_impl(eps, m, a, length)
        if prev is not None:
            worst = mp.mpf(0)
            for p_, c_ in zip(prev, cur):
                d = abs(p_) if c_ == 0 else abs((p_ - c_) / c_)
                worst = max(worst, d)
            if worst < rel_tol:
                return cur, float(worst), dps
        prev = cur
        dps *= 2
    raise RuntimeError(
        f"w_series_k gate failed at dps={dps // 2} (worst={float(worst)})")


def regulated_moments(m, eps, a=1.0, kmax=K_MAX, rel_tol=1e-30):
    """mu_k(eps) = int domega omega^k S(omega; m, a) e^{-eps omega} for
    k = 0..kmax, as mpmath scalars, with the gate certificate.
    Returns (mus, certified_rel_err, dps_used)."""
    ws, err, dps = w_series_k(eps, m, a, kmax, rel_tol)
    mus = [2 * mp.pi * (-1) ** k * factorial(k) * ws[k]
           for k in range(kmax + 1)]
    return mus, err, dps


def regulated_moments_float(m, eps, a=1.0, kmax=K_MAX):
    """Convenience float view of regulated_moments."""
    mus, _, _ = regulated_moments(m, eps, a, kmax)
    return [float(x) for x in mus]


# ---------------------------------------------------------------------------
# study 2: weights-from-moments conditioning for the certified tower
# ---------------------------------------------------------------------------
def tower_masses(num_levels=L_LEVELS):
    """The first num_levels certified 't Hooft masses, BOTH CP sectors
    (the full tower the dictionary must extract), as mpmath scalars in
    sqrt(pi g^2) = 1 units."""
    rows = dictionary.load_reference()
    return [mp.sqrt(mp.mpf(repr(r["two_lambda"])))
            for r in rows[:num_levels]]


def moment_matrix(masses, eps, kmax=K_MAX, rel_tol=COND_TOL):
    """V[k][n] = mu_k(M_n; eps): the linear map from level weights to
    regulated moments.  Entries certified to rel_tol.  Returns
    (rows, max_certified_err)."""
    cols, worst = [], 0.0
    for m_n in masses:
        mus, err, _ = regulated_moments(float(m_n), eps, kmax=kmax,
                                        rel_tol=rel_tol)
        cols.append(mus)
        worst = max(worst, err)
    rows = [[cols[n][k] for n in range(len(masses))]
            for k in range(kmax + 1)]
    return rows, worst


def singular_spectrum(rows, num_moments, column_normalize=True,
                      dps=SVD_DPS):
    """log10 singular values (descending) of the first num_moments rows,
    row-equilibrated to unit 2-norm (the 'relative precision per moment'
    noise model), optionally column-normalized ('shape' spectrum)."""
    ncols = len(rows[0])
    with mp.workdps(dps):
        A = mp.matrix(num_moments, ncols)
        for k in range(num_moments):
            rn = mp.sqrt(mp.fsum([x * x for x in rows[k]]))
            for n in range(ncols):
                A[k, n] = rows[k][n] / rn
        if column_normalize:
            for n in range(ncols):
                cn = mp.sqrt(mp.fsum([A[k, n] ** 2
                                      for k in range(num_moments)]))
                for k in range(num_moments):
                    A[k, n] /= cn
        s = sorted((x for x in mp.svd_r(A, compute_uv=False)), reverse=True)
        return [float(mp.log10(x)) for x in s]


def column_lognorms(rows, num_moments):
    """log10 column 2-norms of the row-equilibrated matrix: the residual
    Boltzmann suppression of heavy levels after the noise model."""
    ncols = len(rows[0])
    with mp.workdps(SVD_DPS):
        rn = [mp.sqrt(mp.fsum([x * x for x in rows[k]]))
              for k in range(num_moments)]
        return [float(mp.log10(mp.sqrt(mp.fsum(
            [(rows[k][n] / rn[k]) ** 2 for k in range(num_moments)]))))
            for n in range(ncols)]


def decay_fit(lg_sv):
    """Least-squares fit log10(sigma_1/sigma_j) = a (j-1) + b (j-1)^2.
    Returns (a, b, rms)."""
    d = np.array([lg_sv[0] - x for x in lg_sv])
    x = np.arange(len(lg_sv), dtype=float)
    X = np.column_stack([x, x ** 2])
    coef, *_ = np.linalg.lstsq(X, d, rcond=None)
    rms = float(np.sqrt(np.mean((X @ coef - d) ** 2)))
    return float(coef[0]), float(coef[1]), rms


def informative_counts(lg_sv, etas=(1e-8, 1e-16, 1e-34)):
    """#{j : sigma_j / sigma_1 > eta}: weight directions above the noise
    floor for data of relative precision eta."""
    return {f"{eta:.0e}": sum(1 for x in lg_sv
                              if lg_sv[0] - x < -mp.log10(eta))
            for eta in etas}


# ---------------------------------------------------------------------------
# study 3: two-mass separability (nonlinear identifiability)
# ---------------------------------------------------------------------------
def two_mass_fit(m0, delta, eps, kmax=7, rel_tol=1e-40):
    """Equal-weight pair at m0 (1 -+ delta/2) against the best single-mass
    fit: (c, m*) matches mu_0 and mu_1 exactly (2 parameters = 2 moments,
    always solvable), residuals live in k >= 2.  Returns dict with m*,
    the relative m* bias, and the relative residual profile r_k."""
    with mp.workdps(60):
        ma = mp.mpf(m0) * (1 - mp.mpf(delta) / 2)
        mb = mp.mpf(m0) * (1 + mp.mpf(delta) / 2)
        A, _, _ = regulated_moments(ma, eps, kmax=kmax, rel_tol=rel_tol)
        B, _, _ = regulated_moments(mb, eps, kmax=kmax, rel_tol=rel_tol)
        tgt = [(x + y) / 2 for x, y in zip(A, B)]
        ratio = tgt[1] / tgt[0]

        def f(mass):
            mm, _, _ = regulated_moments(mass, eps, kmax=1,
                                         rel_tol=rel_tol)
            return mm[1] / mm[0] - ratio

        # the mass iterate stays an mpf end to end, so the matched
        # moments close at solver precision, not float precision
        mstar = mp.findroot(f, mp.mpf(m0), tol=1e-28)
        fit, _, _ = regulated_moments(mstar, eps, kmax=kmax,
                                      rel_tol=rel_tol)
        c = tgt[0] / fit[0]
        res = [float(abs(c * fit[k] / tgt[k] - 1)) for k in range(kmax + 1)]
        return dict(mstar=float(mstar),
                    mstar_rel_bias=float(1 - mstar / mp.mpf(m0)),
                    res=res)


# ---------------------------------------------------------------------------
# study 4: node-level data against moment-level data
# ---------------------------------------------------------------------------
def scaled_amplitude(w, m):
    """g(omega) = e^{pi omega/2} K_{i omega}(m): real, O(1)-amplitude
    oscillation whose zeros are the nodes of S = g^2/pi (a = 1 units).
    The e^{-eps omega} regulator multiplies S and moves NO zero."""
    return float(mp.e ** (mp.pi * mp.mpf(w) / 2)
                 * mp.besselk(1j * mp.mpf(w), mp.mpf(m)).real)


def scan_nodes(m, w_max=10.0, step=0.25):
    """All nodes of S(omega; m, 1) in (0.5, w_max]: sign-change scan of the
    scaled amplitude plus bracketed (anderson) refinement.  Deterministic;
    no unbracketed root-finding (zeros crowd at large omega)."""
    with mp.workdps(30):
        grid = [0.5 + i * step for i in range(int((w_max - 0.5) / step) + 1)]
        vals = [scaled_amplitude(w, m) for w in grid]
        nodes = []
        for a_, b_, va, vb in zip(grid, grid[1:], vals, vals[1:]):
            if va * vb < 0:
                nodes.append(float(mp.findroot(
                    lambda w: scaled_amplitude(w, m), (a_, b_),
                    solver="anderson", tol=1e-24)))
        return nodes


def node_mass_sensitivity(m, node, h=1e-6):
    """Lambda = d ln omega_node / d ln m by centered difference with
    bracketed re-refinement around the given node."""
    with mp.workdps(30):
        out = []
        for mm in (m * (1 + h), m * (1 - h)):
            out.append(mp.findroot(lambda w: scaled_amplitude(w, mm),
                                   (node - 0.2, node + 0.2),
                                   solver="anderson", tol=1e-24))
        return float((out[0] - out[1]) / (2 * h * m) * m / mp.mpf(node))


def node_transfer(m, node, half_window=1.3, n_grid=53):
    """The S-domain location transfer at a (double) zero: S'' at the node
    and the local maximum of S, giving
    d omega = sqrt(2 eta S_ref / S'') for data of relative precision eta."""
    with mp.workdps(30):
        def S(w):
            g = mp.e ** (mp.pi * w / 2) * mp.besselk(1j * w, mp.mpf(m)).real
            return g * g / mp.pi
        spp = float(mp.diff(S, mp.mpf(node), 2))
        grid = np.linspace(node - half_window, node + half_window, n_grid)
        sref = max(float(S(mp.mpf(w))) for w in grid)
        return spp, sref


def kurt_sensitivities(m, eps, h=1e-6):
    """kurt = mu_4 mu_0 / mu_2^2 at fixed eps, with log-derivatives in m
    and in eps, and the regulator leverage |d ln K/d ln eps| /
    |d ln K/d ln m| -- mass digits lost per digit of eps uncertainty."""
    def kurt(mass, e):
        mm, _, _ = regulated_moments(mass, e, kmax=4, rel_tol=1e-30)
        return mm[4] * mm[0] / mm[2] ** 2
    with mp.workdps(40):
        k0 = kurt(m, eps)
        dm = (kurt(m * (1 + h), eps) - kurt(m * (1 - h), eps)) / (2 * h) / k0
        de = (kurt(m, eps * (1 + h)) - kurt(m, eps * (1 - h))) / (2 * h) / k0
        return dict(kurt=float(k0), dlnK_dlnm=float(dm),
                    dlnK_dlneps=float(de),
                    leverage=float(abs(de / dm)))


def tower_node_observables():
    """The coupling-weighted CP = -1 tower (boost_forward.tower) at node
    level: the ground kernel's first node fills in, and the filled depth /
    shifted minimum measure excited weight background-free."""
    import boost_forward as bf
    masses, wgt, _, _ = bf.tower(12)
    with mp.workdps(30):
        w1 = scan_nodes(float(masses[0]), w_max=4.0, step=0.2)[0]

        def S_tower(w):
            return sum(float(wgt[i])
                       * scaled_amplitude(w, float(masses[i])) ** 2
                       for i in range(len(masses))) / float(mp.pi)

        grid = np.linspace(w1 - 1.5, w1 + 1.5, 61)
        vals = [S_tower(float(w)) for w in grid]
        local_max = max(vals)
        i_min = min(range(len(grid)),
                    key=lambda i: vals[i] if abs(grid[i] - w1) < 0.7
                    else np.inf)
        # refine the tower minimum by golden-section on the bracket
        lo, hi = float(grid[max(i_min - 1, 0)]), float(
            grid[min(i_min + 1, len(grid) - 1)])
        phi = (np.sqrt(5.0) - 1) / 2
        a_, b_ = lo, hi
        for _ in range(40):
            c_, d_ = b_ - phi * (b_ - a_), a_ + phi * (b_ - a_)
            if S_tower(c_) < S_tower(d_):
                b_ = d_
            else:
                a_ = c_
        w_min = 0.5 * (a_ + b_)
        return dict(w1_ground=float(w1),
                    S_tower_at_ground_node=S_tower(float(w1)),
                    w_min=w_min, S_tower_min=S_tower(w_min),
                    local_max=local_max,
                    fill_at_ground_node=S_tower(float(w1)) / local_max,
                    min_depth=S_tower(w_min) / local_max,
                    node_shift_rel=(w_min - float(w1)) / float(w1))


# ---------------------------------------------------------------------------
# campaign
# ---------------------------------------------------------------------------
EPS_GRID = (0.1, 0.2, 0.4, 0.8, 1.6)
K_GRID = (4, 8, 12, 16, 20, 24)
EPS_REF = 0.4
DELTAS = (0.4, 0.2, 0.1, 0.05, 0.025)


def run(out_path=None) -> dict:
    t0 = time.time()
    masses = tower_masses(L_LEVELS)

    # ---- study 2: conditioning ------------------------------------------
    cells, law_fit, mats = [], [], {}
    worst_entry_err = 0.0
    for eps in EPS_GRID:
        rows, err = moment_matrix(masses, eps)
        mats[eps] = rows
        worst_entry_err = max(worst_entry_err, err)
        for K in K_GRID:
            lg_shape = singular_spectrum(rows, K, column_normalize=True)
            lg_raw = singular_spectrum(rows, K, column_normalize=False)
            cells.append(dict(eps=eps, K=K, lg_sv_shape=lg_shape,
                              lg_sv_raw=lg_raw,
                              n_above=informative_counts(lg_raw)))
        full = next(c for c in cells if c["eps"] == eps and c["K"] == 24)
        a_c, b_c, rms = decay_fit(full["lg_sv_raw"])
        law_fit.append(dict(eps=eps, a=a_c, b=b_c, rms=rms))
        print(f"eps={eps}: sv12/sv1 = 1e{full['lg_sv_raw'][-1]:.1f}, "
              f"law a={a_c:.2f} b={b_c:.3f} (rms {rms:.2f}), "
              f"n_double={full['n_above']['1e-16']} "
              f"[{time.time() - t0:.0f}s]")

    # fine K-scan at the reference eps: where does information saturate?
    fine = []
    for K in range(3, 25):
        lg = singular_spectrum(mats[EPS_REF], K, column_normalize=False)
        fine.append(dict(K=K, lg_sv_min=lg[-1] - lg[0],
                         n_double=informative_counts(lg)["1e-16"]))
    n_final = fine[-1]["n_double"]
    k_sat = next(r["K"] for r in fine if r["n_double"] == n_final)

    # eps-dependence of the saturation: how far the (sigma_1-relative)
    # spectrum still moves from K = 12 to K = 24 -- the high orders stay
    # mass-informative only while k/eps is inside the mass window
    sv_move = {}
    for eps in EPS_GRID:
        c12 = next(c for c in cells if c["eps"] == eps and c["K"] == 12)
        c24 = next(c for c in cells if c["eps"] == eps and c["K"] == 24)
        sv_move[f"{eps}"] = max(
            abs((x - c12["lg_sv_raw"][0]) - (y - c24["lg_sv_raw"][0]))
            for x, y in zip(c12["lg_sv_raw"], c24["lg_sv_raw"]))

    # variants at the reference eps: even-only rows (the symmetric
    # correlator's observable set) and the CP = -1 6-level sub-tower
    rows_ref = mats[EPS_REF]
    even_rows = [rows_ref[k] for k in range(0, 24, 2)]
    lg_even = singular_spectrum(even_rows, 12, column_normalize=False)
    cp_cols = [0, 2, 4, 6, 8, 10]     # FLZ-even = CP -1 levels
    cp_rows = [[rows_ref[k][n] for n in cp_cols] for k in range(24)]
    lg_cp = singular_spectrum(cp_rows, 24, column_normalize=False)

    conditioning = dict(
        noise_model="rows scaled to unit 2-norm: every moment known to "
                    "the same relative precision; 'shape' additionally "
                    "column-normalized",
        eps_grid=list(EPS_GRID), K_grid=list(K_GRID),
        entry_rel_tol=COND_TOL, svd_dps=SVD_DPS,
        max_certified_entry_err=worst_entry_err,
        trust_floor_log10=-55.0,
        cells=cells, law_fit=law_fit,
        fine_scan_eps_ref=dict(eps=EPS_REF, rows=fine,
                               K_saturation_double=k_sat,
                               n_double_final=n_final),
        sv_move_K12_to_K24=sv_move,
        even_only_eps_ref=dict(orders=list(range(0, 24, 2)),
                               lg_sv_raw=lg_even),
        cp_minus_subtower_eps_ref=dict(
            level_indices=cp_cols, lg_sv_raw=lg_cp,
            n_above=informative_counts(lg_cp)),
        column_lognorms_eps_ref=column_lognorms(rows_ref, 24),
        statement=(
            "Super-geometric ill-conditioning: log10(sv1/svj) ~= "
            "a(eps)(j-1) + b(eps)(j-1)^2 with the fitted (a, b) above "
            "(at a*eps = 0.1 the mid-range is steeper than the fit and "
            "the last singular values floor near 1e-17, hence the rms); "
            "the cost per level grows, smaller eps is worse, the level "
            "count saturates at K = L moments, and further moments help "
            "only at large eps (sv_move_K12_to_K24).  All conditional "
            "on the e^{-eps omega} regulator convention."))

    # ---- study 3: two-mass separability ---------------------------------
    scan = [dict(delta=d, eps=EPS_REF, **two_mass_fit(1.0, d, EPS_REF))
            for d in DELTAS]
    lr = np.log10([r["res"][2] for r in scan])
    ld = np.log10(DELTAS)
    slope = float(np.polyfit(ld, lr, 1)[0])
    C_ref = float(scan[-1]["res"][2] / DELTAS[-1] ** 2)
    eps_var = []
    for eps in (0.1, 1.6):
        r = two_mass_fit(1.0, 0.1, eps)
        eps_var.append(dict(eps=eps, delta=0.1, res2=r["res"][2],
                            C=r["res"][2] / 0.01))
    bias_c = float(scan[-1]["mstar_rel_bias"] / DELTAS[-1] ** 2)
    separability = dict(
        protocol="equal-weight pair at m0(1 -+ delta/2), m0 = 1; single-"
                 "mass fit (c, m*) matches mu_0, mu_1 exactly; residuals "
                 "r_k = |c mu_k(m*)/mu_k^pair - 1| for k >= 2",
        scan=scan, fit_exponent=slope, C_eps_ref=C_ref,
        eps_variants=eps_var, mstar_bias_coeff=bias_c,
        min_moments=3,
        delta_min=dict((f"{eta:.0e}", float(np.sqrt(eta / C_ref)))
                       for eta in (1e-16, 1e-8, 1e-4)),
        statement=(
            "Two moments always fit a pair exactly, so 3 is the minimum "
            "that can distinguish; the third moment's signal is "
            f"C(eps) delta^2 (exponent {slope:.3f}), C = {C_ref:.3f} at "
            f"a*eps = {EPS_REF}, so delta_min = sqrt(eta/C) and the "
            f"recovered single mass is biased low by {bias_c:.2f} "
            "delta^2."))
    print(f"separability: exponent {slope:.3f}, C({EPS_REF}) = {C_ref:.4f} "
          f"[{time.time() - t0:.0f}s]")

    # ---- study 4: nodes against moments ---------------------------------
    nodes_m1 = scan_nodes(1.0, w_max=10.0)
    nodes_m2 = scan_nodes(2.0, w_max=12.0)
    lam = [node_mass_sensitivity(1.0, w) for w in nodes_m1[:3]]
    spp, sref = node_transfer(1.0, nodes_m1[0])
    node_pref = float(np.sqrt(2 * sref / spp) / nodes_m1[0] / lam[0])
    ksens = [dict(eps=e, **kurt_sensitivities(1.0, e))
             for e in (0.2, 0.4, 0.8)]
    lam_kurt = abs(next(r for r in ksens if r["eps"] == EPS_REF)
                   ["dlnK_dlnm"])
    digits = []
    for eta in (1e-4, 1e-8, 1e-12):
        digits.append(dict(
            eta=eta,
            dm_node=node_pref * float(np.sqrt(eta)),
            dm_kurt_eps_exact=eta / lam_kurt,
            dm_floor_eps_1pct=0.01 * next(
                r for r in ksens if r["eps"] == EPS_REF)["leverage"]))
    tower_nodes = tower_node_observables()
    nodes = dict(
        nodes_m1=nodes_m1, nodes_m2=nodes_m2,
        lambda_node=lam,
        node1_curvature=dict(Spp=spp, S_ref=sref,
                             transfer_prefactor=node_pref),
        kurt_sens=ksens, digits_table=digits,
        eps_leverage_of_nodes=0.0,
        crossover=(
            "nodes carry more mass digits than kurt whenever the relative "
            f"eps uncertainty u > {node_pref:.2f} sqrt(eta); with the "
            "state map unfixed (u = O(1)) moments carry zero digits"),
        tower=tower_nodes,
        statement=(
            "Node positions are exactly regulator-independent (the "
            "e^{-eps omega} factor moves no zero) but S has double "
            f"zeros, so location transfers as {node_pref:.2f} sqrt(eta); "
            f"kurt transfers linearly (eta/{lam_kurt:.2f}) but carries "
            "regulator leverage ~1: one digit of eps buys one digit of "
            "mass, at any data precision."))
    print(f"nodes: omega_1(m=1) = {nodes_m1[0]:.4f}, Lambda = {lam[0]:.3f}, "
          f"leverage = {[round(r['leverage'], 2) for r in ksens]} "
          f"[{time.time() - t0:.0f}s]")

    # ---- verdict ---------------------------------------------------------
    ref_fit = next(r for r in law_fit if r["eps"] == EPS_REF)
    ref_cell = next(c for c in cells
                    if c["eps"] == EPS_REF and c["K"] == 24)
    verdict = dict(
        headline_negative=True,
        conditioning_law=(
            f"log10(sv1/svj) ~= a(j-1) + b(j-1)^2, a = "
            f"{law_fit[0]['a']:.2f}..{law_fit[-1]['a']:.2f} (eps = "
            f"{EPS_GRID[0]}..{EPS_GRID[-1]}), b ~= "
            f"{np.mean([r['b'] for r in law_fit]):.3f}"),
        digits_per_level_ref=dict(
            eps=EPS_REF, first=ref_fit["a"] + ref_fit["b"],
            last=ref_fit["a"] + ref_fit["b"] * (2 * (L_LEVELS - 1) - 1)),
        moments_needed=(
            f"the double-precision level count saturates at K = {k_sat} "
            f"moments for {L_LEVELS} levels at every eps on the grid; "
            f"beyond K = L the singular values themselves move by "
            f"{sv_move[f'{EPS_GRID[0]}']:.1f} decades (a*eps = "
            f"{EPS_GRID[0]}) up to {sv_move[f'{EPS_GRID[-1]}']:.1f} "
            f"(a*eps = {EPS_GRID[-1]}): high orders stay informative "
            "only while k/eps is inside the mass window"),
        n_levels_double_precision={
            f"{c['eps']}": c["n_above"]["1e-16"]
            for c in cells if c["K"] == 24},
        amplification_last_level_ref=(ref_cell["lg_sv_raw"][0]
                                      - ref_cell["lg_sv_raw"][-1]),
        separability_scaling=(
            f"r = C(eps) delta^2, C = {eps_var[0]['C']:.3f}.."
            f"{eps_var[1]['C']:.3f} over a*eps = 0.1..1.6; 3 moments "
            "minimum"),
        node_vs_moment=(
            f"nodes: eps-leverage 0, dm/m = {node_pref:.2f} sqrt(eta); "
            f"moments: dm/m = eta/{lam_kurt:.2f} conditional on exact "
            "eps, plus ~1.0 x (relative eps uncertainty) always"),
        statement=(
            "On the solvable kernel the moments -> spectrum inversion is "
            "super-geometrically ill-conditioned (measured exponents "
            "above), needs exactly as many moments as levels, and every "
            "moment-side mass digit is hostage to the regulator scale "
            "the missing tomperature/state map must fix.  Node-level "
            "data has zero regulator leverage.  The roadmap's 'the "
            "inversion needs node-level data' verdict stands, "
            "quantified.  All numbers are conditional on the "
            "e^{-eps omega} regulator convention; composite-operator "
            "smearing and interactions would dress the kernel."))

    payload = dict(
        provenance="boost_inverse.py -- conditioning of the moments -> "
                   "spectrum inverse problem on the exactly solvable "
                   "boost kernel (BOOST.md), answering road_map.md Part "
                   "III's 'how many moments' question: certified "
                   "arbitrary-order regulated moments (K_0 Taylor "
                   "composition, dps-doubling gate), mpmath SVD "
                   "conditioning of the weights->moments map for the "
                   "certified 12-level tower, two-mass identifiability, "
                   "and the node-against-moment information comparison",
        units="masses in sqrt(pi g^2) = 1, a = 1: only (m/a, a*eps) "
              "enter; regulator e^{-eps omega} as in boost_forward.py",
        engine=dict(k_max=K_MAX, gate="dps-doubling to rel_tol",
                    default_rel_tol=1e-30,
                    conditioning_rel_tol=COND_TOL),
        tower=dict(masses=[float(m) for m in masses],
                   source="thooft-target/reference_spectrum.json, first "
                          f"{L_LEVELS} levels, both CP sectors"),
        conditioning=conditioning, separability=separability,
        nodes=nodes, verdict=verdict,
        runtime_s=round(time.time() - t0, 1))
    out = Path(out_path) if out_path else HERE / "boost_inverse.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"verdict: {verdict['moments_needed']}")
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


if __name__ == "__main__":
    run()
