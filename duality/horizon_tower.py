#!/usr/bin/env python3
r"""
The singlet tower at the stretched horizon: quasinormal ladders, not mesons.
===========================================================================

singlet_tower.py writes every Majorana singlet channel
A_n = sum_i psi_i (ad_H)^n psi_i (n <= 6) at N = infinity as a conserved
polynomial in H plus matter chords of dimension 2..n.  This module asks
what that content looks like in the frame where the proposal says the
mesons live, and whether a massive meson could be hiding in it.

THE ORTHOGONALIZED TOWER
------------------------
The decomposition is triangular, and the chord that is new at level n
enters A_n with a constant dressing, coefficient exactly 4 (the a = 1,
b = n term of kappa_ab, see singlet_tower.py).  So projecting A_n off the
commutant of H and off every dressed lower chord f(H) O_m g(H), m < n,
leaves

    O~_n = 4 O^_n ,     <O~_n O~_n> = 16 norm_n ,

a pure matter chord of dimension Delta = n.  orthogonal_tower() checks
the coefficient and the weight against sector_moment.

THE SEMICLASSICAL FORM: A QUASINORMAL LADDER
--------------------------------------------
As lambda -> 0 the T = infinity two-point function of a matter chord of
dimension Delta tends to sech^(2 Delta)(Jc t) (Jc the chord-units scale
of chord.py), whose spectral density is

    rho_Delta(omega) = 2^(2 Delta - 1) |Gamma(Delta + i omega/(2 Jc))|^2
                       / (2 pi Jc Gamma(2 Delta)) .

It is meromorphic, with poles at omega = -2 i Jc (Delta + k), k = 0, 1, ...:
a quasinormal ladder with spacing a = 2 Jc and lowest decay rate
Delta * a.  Measured on the exact N = infinity density (chord.py,
s_omega_smooth) for Delta = 2..6, the maximum deviation is (Delta/4)
lambda + O(lambda^2) (convergence()).  In the time domain the effective
decay rate -d ln G/dt plateaus at Delta * a (decay_rates()).

Every dressed sector of a full channel inherits the same poles: at small
lambda the dressing contributes a polynomial in omega.  full_channel_fit()
fits the exact A_3 and A_4 densities to span{omega^(2j) rho_m}; the
residual falls steeply with lambda (about lambda^3).  That fit is a
consistency check, not a discriminator on its own: a massive-field
density also fits the same basis to a few 1e-3 on a finite window
(massive_fit_control), but at a fixed floor with nothing to send to zero.
The sharp discriminator is the time domain.

UNITS: THE LADDER IS STRING-SCALE
---------------------------------
Miyashita-Sekino-Susskind 2506.18054 (Sec. 2) measure energies in string
units, J_MSS = p J with this repo's coupling Var(J) = p!/N^(p-1) J^2,
J_MSS the inverse string length; the flat-space limit keeps energies
fixed in those units as N -> infinity at fixed lambda.  With
<H^2> = C(N,p) p!/N^(p-1) J^2 -> J_MSS^2 e^(-lambda/4) / lambda_paper,

    scriptJ = Jc sqrt(<H^2>) = sqrt(2) J_MSS ,   a = 2 scriptJ = 2 sqrt(2) J_MSS,

exactly at the semiclassical bridge (units_mss()).  So the ladder of the
singlet tower sits at the string scale: the scale of the meson masses,
and of the stretched horizon of 2506.18054 Sec. 6, which sits at a string
distance from the Rindler horizon with a string-scale local Unruh
temperature.  (The fundamental fermion, Delta = 1/p, decays at
2 scriptJ/p, the cosmic scale of Rahman-Susskind 2401.08555.)

THE CALIBRATION: WHAT A MASS LOOKS LIKE AT A RINDLER HORIZON
------------------------------------------------------------
The T = infinity correlator Tr(O(t) O)/dim is regular and real at t = 0:
it is a two-sided (thermofield-double) correlator, and so is
sech^(2 Delta).  For a probe at Rindler radius 1/a, two-sided free-field
correlators are:

  * a scale-invariant operator of dimension Delta:
        cosh^(-2 Delta)(a tau/2)        <->  |Gamma(Delta + i omega/a)|^2 ,
    the quasinormal ladder -i a (Delta + k);
  * a free field of mass M:
        K_0((2M/a) cosh(a tau/2))       <->  K_(i omega/a)(M/a)^2 ,
    by DLMF 10.32.17, K_mu K_nu = 2 int K_(mu+nu)(2z cosh t)
    cosh((mu-nu)t) dt.  That is ENTIRE in omega: no quasinormal poles.
    Its effective decay rate grows like M sinh(a tau/2), without bound,
    and the density has real zeros (nodes) at large omega.

The one-sided (same-wedge) responses differ by the KMS factor e^(pi omega/a):
e^(pi omega/a) |Gamma(Delta + i omega/a)|^2 against e^(pi omega/a)
K_(i omega/a)(M/a)^2 (duality/BOOST.md).  The factor is entire, so the
dichotomy, poles against no poles, does not depend on that choice; the
two-sided form is used because the hologram correlator is omega-symmetric.

A gapped tower of mesons coupled with positive weights, like the 't Hooft
model at m_q^2 = 0 (gapped, FLZ 0905.2280), is a positive sum of the
second kind: entire, super-exponentially decaying (thooft_tower_rates(),
with the repo's own masses and decay constants).

With a = 2 Jc the N = infinity singlet channels are, as lambda -> 0,
exactly the scale-invariant responses with Delta = 2, 3, 5.  The argument
is the exact decomposition (every A_n is a polynomial in H plus matter
chords, nothing else) plus the measured convergence of every chord's
density to |Gamma(m + i omega/a)|^2; so the channel densities are sums of
polynomial multiples of those ladders, with no entire K^2 component for
any M > 0.  The time domain shows it directly: the exact decay rates
plateau at Delta * a, while any massive field, and the 't Hooft tower,
has a rate that grows without bound.  In de Sitter
language the same statement reads: a real Delta > 1 corresponds to
M^2 l^2 = Delta (1 - Delta) < 0 in dS_2 (4 Delta(1 - Delta) in the dS_3
convention of Narovlansky-Verlinde 2310.16994, Eq. 3), while a massive
particle with M l > 1/2 has complex Delta = 1/2 + i sqrt(M^2 l^2 - 1/4)
and quasinormal modes with nonzero real part (ds_reading()).

WHAT THIS DOES AND DOES NOT SAY
-------------------------------
It says that at N = infinity and small lambda the singlet channels carry
the analytic structure of scale-invariant operators of dimension 2, 3, 5
at the string-scale ladder spacing, and no massive component at all.
That is consistent with the picture of 2506.18054 Sec. 6, in which the
stretched horizon is the deconfined, effectively conformal region and
mesons are what escapes it; it means the mesons, if the duality holds,
are not visible in these channels at N = infinity.

It does not exclude mesons at finite N or finite lambda (at finite lambda
the spectrum is bounded and every correlator eventually decays as a power
law).  The calibration is for free fields probed at a fixed distance
1/a from the horizon; a dictionary that smears operators down to the
horizon (radius << 1/M) sees a conformal window up to a tau ~
(2/a) ln(1/(M rho)).  The identification a = 2 Jc (ladder spacing =
acceleration) is the natural one but is not supplied by the proposal.
Majorana A_n only.

Output: horizon_tower.json.  `python horizon_tower.py` prints the
headline numbers and writes the artifact (about two minutes).
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import gammaln, k0e, loggamma

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# the semiclassical ladder
# ---------------------------------------------------------------------------
def jc_bridge(lam: float) -> float:
    """The semiclassical frequency scale in chord units (<H^2> = 1):
    G_Delta(t) -> sech^(2 Delta)(Jc t), Jc = sqrt(lambda) e^(lambda/8)
    (chord.py, run_spectra)."""
    return math.sqrt(lam) * math.exp(lam / 8.0)


def ladder_density(omegas, Jc: float, delta: float) -> np.ndarray:
    """Spectral density of sech^(2 Delta)(Jc t), normalized to unit
    integral over omega:  2^(2D-1) |Gamma(D + i w/2Jc)|^2 / (2 pi Jc Gamma(2D))."""
    x = np.asarray(omegas, dtype=float) / (2.0 * Jc)
    lg = (2.0 * np.real(loggamma(delta + 1j * x)) + (2 * delta - 1) * math.log(2.0)
          - gammaln(2 * delta) - math.log(Jc))
    return np.exp(lg) / (2.0 * math.pi)


def quasinormal_ladder(delta: float, a: float, kmax: int = 4) -> list:
    """The poles of |Gamma(Delta + i omega/a)|^2 in the lower half plane:
    omega_k = -i a (Delta + k)."""
    return [complex(0.0, -a * (delta + k)) for k in range(kmax + 1)]


# ---------------------------------------------------------------------------
# the orthogonalized tower
# ---------------------------------------------------------------------------
def orthogonal_tower(lam: float, nmax: int = 5) -> list:
    """For n = 2..nmax: the coefficient of the new chord O^_n in A_n (with
    its dressing), its norm, and the weight 16 norm_n of O~_n = 4 O^_n,
    checked against sector_moment(n, 0, n); with its share of A_n's
    dynamical weight."""
    from singlet_tower import TowerChord
    q = math.exp(-lam)
    T = TowerChord(lam, q)
    T.decompose(nmax - 2, check_len=0)
    rows = []
    for n in range(2, nmax + 1):
        kern = T.a_kernels(n)
        top = kern[n]
        weight = float(T.sector_moment(n, 0, n))
        total = float(sum(T.sector_moment(n, 0, m) for m in kern))
        rows.append(dict(n=n, top_kernel={f"{i},{j}": float(c) for (i, j), c in top.items()},
                         norm=float(T.norm[n]), weight_16norm=16.0 * float(T.norm[n]),
                         weight_sector=weight, share_of_dynamic=weight / total))
    return rows


# ---------------------------------------------------------------------------
# convergence of the exact chord density to the ladder
# ---------------------------------------------------------------------------
def convergence(delta: int, lams=(0.4, 0.2, 0.1, 0.05), n_omega: int = 160) -> dict:
    """Max |S_exact - rho_Delta| / max rho_Delta on [0, 6 Jc sqrt(Delta)],
    S_exact the exact N = infinity, T = infinity density of a Delta-chord
    (chord.s_omega_smooth), and the ratio dev/lambda."""
    from chord import q_of_lambda, s_omega_smooth
    rows = []
    for lam in lams:
        q = q_of_lambda(lam)
        Jc = jc_bridge(lam)
        band = 2.0 / math.sqrt(1.0 - q)
        om = np.linspace(0.0, min(6.0 * Jc * math.sqrt(delta), 1.998 * band), n_omega)
        n_theta = int(min(4000, max(800, 60.0 / lam)))
        S = s_omega_smooth(om, q, float(delta), n_theta)
        G = ladder_density(om, Jc, delta)
        dev = float(np.max(np.abs(S - G)) / G.max())
        rows.append(dict(lam=lam, dev=dev, dev_over_lam=dev / lam))
    return dict(delta=delta, rows=rows)


# ---------------------------------------------------------------------------
# time domain: the decay-rate plateau
# ---------------------------------------------------------------------------
def decay_rates(delta: int, lam: float, a_taus=(2, 4, 6, 8), n_t: int = 321) -> dict:
    """Effective decay rate -d ln G/dt of the exact chord correlator
    (chord.g2_engine), in units of a = 2 Jc, at the given a*t."""
    from chord import g2_engine, q_of_lambda
    q = q_of_lambda(lam)
    a = 2.0 * jc_bridge(lam)
    nmax = int(min(2500, 40.0 / lam + 100))
    ts = np.linspace(0.0, max(a_taus) / a, n_t)
    G = g2_engine(ts, q, float(delta), nmax).real
    rate = -np.gradient(np.log(np.abs(G)), ts) / a
    return dict(delta=delta, lam=lam, nmax=nmax,
                rate_over_a={str(x): float(np.interp(x / a, ts, rate)) for x in a_taus})


def massive_rates(m_over_a: float, a_taus=(2, 4, 6, 8)) -> dict:
    """Effective decay rate of the two-sided massive correlator
    K_0((2M/a) cosh(a tau/2)), in units of a."""
    x = np.linspace(0.0, max(a_taus), 801)
    z = 2.0 * m_over_a * np.cosh(x / 2.0)
    lg = np.log(k0e(z)) - z
    rate = -np.gradient(lg, x)
    return dict(m_over_a=m_over_a,
                rate_over_a={str(t): float(np.interp(t, x, rate)) for t in a_taus})


def thooft_tower_rates(m0_over_a: float, n_levels: int = 12, a_taus=(2, 4, 6, 8)) -> dict:
    """The same for the 't Hooft tower at the duality point: masses
    M_n = M_0 sqrt(2 lambda_n / 2 lambda_0) and weights f_n^2 from the
    symmetric-sector decay constants (thooft-target/eigenfunctions.py),
    G = sum f_n^2 K_0((2 M_n/a) cosh(a tau/2)).  Positive and gapped, so
    its rate grows without bound like the lightest meson's."""
    sys.path.insert(0, str(HERE.parent / "thooft-target"))
    from eigenfunctions import decay_constants
    f, ev = decay_constants(200, n_levels)
    masses = m0_over_a * np.sqrt(np.asarray(ev[:n_levels]) / ev[0])
    x = np.linspace(0.0, max(a_taus), 801)
    ch = np.cosh(x / 2.0)
    # sum f^2 K0(2 m cosh) in log form, relative to the lightest exponent
    z0 = 2.0 * masses[0] * ch
    terms = [f[i] ** 2 * k0e(2.0 * masses[i] * ch) * np.exp(-2.0 * (masses[i] - masses[0]) * ch)
             for i in range(n_levels)]
    lg = np.log(np.sum(terms, axis=0)) - z0
    rate = -np.gradient(lg, x)
    return dict(m0_over_a=m0_over_a, n_levels=n_levels,
                rate_over_a={str(t): float(np.interp(t, x, rate)) for t in a_taus})


# ---------------------------------------------------------------------------
# the calibration in frequency: entire against meromorphic
# ---------------------------------------------------------------------------
def bessel_identity(points=((1.0, 0.7), (2.5, 3.0), (0.3, 1.5))) -> list:
    """K_(i nu)(z)^2 = 2 int_0^inf K_0(2 z cosh t) cos(2 nu t) dt
    (DLMF 10.32.17 with mu = -nu = i nu), checked in mpmath."""
    import mpmath as mp
    out = []
    with mp.workdps(20):
        for z, nu in points:
            lhs = mp.re(mp.besselk(1j * nu, z) ** 2)
            rhs = 2 * mp.quad(lambda t: mp.besselk(0, 2 * z * mp.cosh(t))
                              * mp.cos(2 * nu * t), [0, 1, 2, 4])
            out.append(dict(z=z, nu=nu, lhs=float(lhs), rhs=float(rhs),
                            rel=float(abs(lhs - rhs) / abs(lhs))))
    return out


def massive_two_sided_density(omegas, m_over_a: float, a: float = 1.0) -> np.ndarray:
    """S(omega) = (1/pi) int_0^inf cos(omega tau) K_0((2M/a) cosh(a tau/2)) dtau
    by quadrature: an independent route to the closed form
    K_(i omega/a)(M/a)^2 / (pi a) of DLMF 10.32.17."""
    tmax = 2.0 / a * math.acosh(max(40.0 / (2.0 * m_over_a), 1.0)) + 2.0 / a
    t = np.linspace(0.0, tmax, 20001)
    z = 2.0 * m_over_a * np.cosh(a * t / 2.0)
    G = k0e(z) * np.exp(-z)
    return np.array([np.trapezoid(G * np.cos(w * t), t) / math.pi
                     for w in np.atleast_1d(omegas)])


def node_count(values) -> int:
    v = np.asarray(values)
    return int(np.sum(np.diff(np.sign(v[np.abs(v) > 0])) != 0))


# ---------------------------------------------------------------------------
# full channels: polynomial times ladder
# ---------------------------------------------------------------------------
def full_channel_fit(n: int, lam: float, n_omega: int = 120, deg: int = 2) -> dict:
    """Fit the exact dynamical density of A_n (singlet_tower.spectral_density)
    on [0, 8 Jc] to span{omega^(2j) rho_m(omega): m in its sectors,
    j <= deg}, least squares; report the relative residual."""
    from singlet_tower import TowerChord, spectral_density
    Jc = jc_bridge(lam)
    q = math.exp(-lam)
    om = np.linspace(0.0, 8.0 * Jc, n_omega)
    n_theta = int(min(3000, max(800, 60.0 / lam)))
    S = np.array(spectral_density(n, om, lam, n_theta=n_theta)["total"])
    sectors = sorted(TowerChord(lam, q).a_kernels(n))
    cols = [(om / (2 * Jc)) ** (2 * j) * ladder_density(om, Jc, m)
            for m in sectors for j in range(deg + 1)]
    A = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(A, S, rcond=None)
    resid = float(np.max(np.abs(A @ c - S)) / np.max(np.abs(S)))
    return dict(n=n, lam=lam, sectors=sectors, deg=deg, rel_residual=resid)


def massive_fit_control(m_over_a: float, sectors=(2, 3), deg: int = 2,
                        n_omega: int = 120) -> dict:
    """The same fit applied to the two-sided massive density on the same
    window (a = 2 Jc = 1): a finite-window least-squares fit is not by
    itself a test of analyticity, and this records the floor it reaches
    for a field that has no quasinormal modes."""
    Jc = 0.5
    om = np.linspace(0.0, 8.0 * Jc, n_omega)
    S = massive_two_sided_density(om, m_over_a, a=2.0 * Jc)
    cols = [(om / (2 * Jc)) ** (2 * j) * ladder_density(om, Jc, m)
            for m in sectors for j in range(deg + 1)]
    A = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(A, S, rcond=None)
    return dict(m_over_a=m_over_a, sectors=list(sectors), deg=deg,
                rel_residual=float(np.max(np.abs(A @ c - S)) / np.max(np.abs(S))))


# ---------------------------------------------------------------------------
# units and the de Sitter reading
# ---------------------------------------------------------------------------
def units_mss(lam_paper: float) -> dict:
    """Chord units against MSS string units (2506.18054 Sec. 2,
    J_MSS = p J): <H^2> -> e^(-lambda/4) J_MSS^2 / lambda_paper with
    lambda = 2 lambda_paper, so scriptJ = Jc sqrt(<H^2>) = sqrt(2) J_MSS."""
    lam = 2.0 * lam_paper
    E_c = math.exp(-lam / 8.0) / math.sqrt(lam_paper)      # in units of J_MSS
    scriptJ = jc_bridge(lam) * E_c
    return dict(lambda_paper=lam_paper, lambda_chord=lam,
                chord_unit_over_JMSS=E_c, scriptJ_over_JMSS=scriptJ,
                ladder_spacing_over_JMSS=2.0 * scriptJ)


def ds_reading(deltas=(2, 3, 5)) -> list:
    """M^2 l^2 = Delta(1 - Delta) (dS_2) and 4 Delta(1 - Delta) (the dS_3
    convention of NV 2310.16994 Eq. 3): negative for every Delta > 1."""
    return [dict(delta=d, dS2_M2l2=d * (1 - d), NV_dS3_m2=4 * d * (1 - d)) for d in deltas]


# ---------------------------------------------------------------------------
def run(out_path=None) -> dict:
    t0 = time.time()
    payload = dict(
        provenance="horizon_tower.py -- the orthogonalized singlet tower, the "
                   "lambda -> 0 quasinormal ladder of each chord, and the "
                   "Rindler calibration against massive fields and the "
                   "'t Hooft tower",
        orthogonal_tower={str(lam): orthogonal_tower(lam) for lam in (2.0, 0.5)},
        convergence=[convergence(d) for d in (2, 3, 4, 5, 6)],
        decay_rates=[decay_rates(d, lam) for d in (2, 3) for lam in (0.5, 0.1, 0.03)],
        massive_rates=[massive_rates(m) for m in (0.1, 0.5, 1.0)],
        thooft_tower_rates=[thooft_tower_rates(m) for m in (0.1, 0.5, 1.0)],
        bessel_identity=bessel_identity(),
        full_channel_fit=[full_channel_fit(n, lam) for n in (3, 4)
                          for lam in (0.4, 0.2, 0.1, 0.05)],
        massive_fit_control=[massive_fit_control(m, s) for m in (0.1, 0.5, 1.0)
                             for s in ((2, 3), (2, 3, 4))],
        units_mss=[units_mss(lp) for lp in (1.0, 0.25, 0.05)],
        ds_reading=ds_reading())
    om = np.linspace(0.0, 12.0, 241)
    payload["nodes"] = dict(
        massive={str(m): node_count(massive_two_sided_density(om, m)) for m in (0.5, 1.0)},
        ladder={str(d): node_count(ladder_density(om, 0.5, d)) for d in (2, 3, 5)},
        omega_window=[0.0, 12.0], a=1.0)
    payload["runtime_s"] = round(time.time() - t0, 1)
    out = Path(out_path) if out_path else HERE / "horizon_tower.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    return payload


if __name__ == "__main__":
    p = run()
    print("Orthogonalized tower, lambda_chord = 2:  O~_n = 4 O^_n")
    for r in p["orthogonal_tower"]["2.0"]:
        print(f"  n={r['n']}  top kernel {r['top_kernel']}  16 norm = {r['weight_16norm']:.6f}"
              f"  (sector {r['weight_sector']:.6f})  share of A_n {r['share_of_dynamic']:.4f}")
    print("\nExact chord density against the ladder |Gamma(D + i w/2Jc)|^2:")
    for c in p["convergence"]:
        print(f"  Delta={c['delta']}: dev/lambda = "
              + ", ".join(f"{r['dev_over_lam']:.3f}@{r['lam']}" for r in c["rows"])
              + f"   (Delta/4 = {c['delta'] / 4:.3f})")
    print("\nDecay rate / a at a*t = 2, 4, 6, 8:")
    for r in p["decay_rates"]:
        print(f"  chord Delta={r['delta']} lambda={r['lam']}: "
              + " ".join(f"{v:5.2f}" for v in r["rate_over_a"].values()))
    for r in p["massive_rates"]:
        print(f"  free field M/a={r['m_over_a']}:      "
              + " ".join(f"{v:5.2f}" for v in r["rate_over_a"].values()))
    for r in p["thooft_tower_rates"]:
        print(f"  't Hooft tower M0/a={r['m0_over_a']}:  "
              + " ".join(f"{v:5.2f}" for v in r["rate_over_a"].values()))
    print("\nFull channels as polynomial x ladder (relative residual):")
    for r in p["full_channel_fit"]:
        print(f"  A_{r['n']} lambda={r['lam']}: sectors {r['sectors']}  residual {r['rel_residual']:.2e}")
    for r in p["massive_fit_control"]:
        print(f"  control: free field M/a={r['m_over_a']} on sectors {r['sectors']}: "
              f"residual {r['rel_residual']:.2e} (fixed floor)")
    print("\nNodes on [0, 12a]:", p["nodes"])
    print("Units (MSS string units):", p["units_mss"][0])
    print(f"\nwrote horizon_tower.json ({p['runtime_s']} s)")
