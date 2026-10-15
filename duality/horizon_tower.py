#!/usr/bin/env python3
r"""
The singlet tower at the stretched horizon: a gapless, scale-invariant spectrum.
==============================================================================

singlet_tower.py writes every Majorana singlet channel
A_n = sum_i psi_i (ad_H)^n psi_i (n <= 6) at N = infinity as a conserved
polynomial in H plus matter chords of dimension 2..n.  This module asks
what that content looks like in the frame where the proposal says the
mesons live, and what mass spectrum it corresponds to there.

THE ORTHOGONALIZED TOWER
------------------------
The decomposition is triangular, and the chord that is new at level n
enters A_n with a constant dressing, coefficient exactly 4 (the T = {1, n}
term (-2)^2 C(n-1, n-1) of the A_n identity).  So projecting A_n off the
commutant of H and off every dressed lower chord f(H) O_m g(H), m < n,
leaves O~_n = 4 O^_n, of weight 16 norm_n: a pure matter chord of
dimension Delta = n (orthogonal_tower()).

THE SEMICLASSICAL FORM: A QUASINORMAL LADDER
--------------------------------------------
As lambda -> 0 the T = infinity two-point function of a matter chord of
dimension Delta tends to sech^(2 Delta)(Jc t), Jc = sqrt(lambda) e^(lambda/8)
in chord units, with spectral density

    rho_Delta(omega) = 2^(2 Delta - 1) |Gamma(Delta + i omega/(2 Jc))|^2
                       / (2 pi Jc Gamma(2 Delta)) ,

poles at omega = -2 i Jc (Delta + k): a quasinormal ladder with spacing
a = 2 Jc.  Against the exact N = infinity density (chord.s_omega_smooth),
Delta = 2..6, the deviation at the fixed scale Jc is (Delta/4) lambda +
O(lambda^2) (convergence()).  Almost all of it is a scale shift: at the
best-fit scale J_eff = Jc (1 - (Delta/4) lambda + ...) the shape deviates
by about 0.045 lambda for every Delta (scale_fit()).  So one common
spacing a for all channels holds at leading order only.  In time, the
exact decay rate -d ln G/dt plateaus at Delta * a (decay_rates()).

A full channel's density is a sum of polynomial multiples of its sectors'
ladders, since a dressing is a polynomial in omega at small lambda.
full_channel_fit() checks that; it is a consistency check only, because
|Gamma(m + 1 + ix)|^2 = (m^2 + x^2) |Gamma(m + ix)|^2 makes the basis
rank-deficient (an even polynomial times rho_2), and a massive-field
density also fits it on a finite window to ~5e-4 (massive_fit_control).

UNITS
-----
2506.18054 measures energies in string units, J_MSS = p J, and keeps them
fixed in the flat-space limit at fixed lambda.  In the Majorana
translation used here (Var(J) = p! J^2 / N^(p-1), lambda_chord =
2 lambda_paper), <H^2> -> J_MSS^2 e^(-lambda/4)/lambda_paper and
scriptJ = Jc sqrt(<H^2>) = sqrt(2) J_MSS (units_mss(); exact because the
e^(lambda/8) in Jc cancels the e^(-lambda/4)).  The O(1) factor is
convention-dependent: the paper's complex model (q = e^(-4 lambda) in
2006.13983's convention) changes it.  Either way the ladder sits at the
string scale, the scale of the meson masses and of the stretched horizon
of 2506.18054 Sec. 6.

THE CALIBRATION: MASS SPECTRA AT A RINDLER HORIZON
--------------------------------------------------
The T = infinity correlator Tr(O(t) O)/dim is real, even and regular at
t = 0, the structure of a two-sided correlator (sech^(2 Delta)(scriptJ t)
first becomes singular at t = i beta_R/2, beta_R = 2 pi/a).  For a probe
at Rindler radius 1/a, a free field of mass M gives

    K_0((2M/a) cosh(a tau/2))  <->  S_M(omega) = K_(i omega/a)(M/a)^2 / (pi a)

(DLMF 10.32.17; bessel_identity()), entire in omega, with nodes, and a
decay rate (M/a) sinh(a tau/2) + 1/4 tanh(a tau/2) that grows without
bound.  A field theory with Kallen-Lehmann mass density rho(M) gives
S(omega) = int rho(M) S_M(omega) dM, and the ladder is such a
superposition:

    int_0^inf z^(2 Delta - 1) K_(ix)(z)^2 dz
        = 2^(2 Delta - 3) Gamma(Delta)^2 |Gamma(Delta + ix)|^2 / Gamma(2 Delta)

(the Mellin transform of K_nu^2 at nu = ix; mass_spectrum_identity()
checks it).  So a dimension-Delta ladder
is the response of the GAPLESS, scale-invariant mass density
rho(M) ~ M^(2 Delta - 1), down to M = 0; its poles come from the M -> 0
end.  A gapped spectrum (support M >= M_0 > 0), like the discrete 't Hooft
tower at m_q^2 = 0 (FLZ), gives an entire S(omega): no poles.  The
transform is invertible (the cosine representation turns S into the
K_0-transform of rho, which is injective), so the exact channel densities,
sums of |Gamma(m + i omega/a)|^2 with polynomial dressings, fix rho(M) as
a power-law density M^(2m-1)(c_0 + c_1 M^2 + ...) with no gap and no
discrete masses.

The decay rates illustrate the same thing but do not prove it: the
plateau at Delta * a is the slowest, M -> 0 part, and a gapped admixture
would decay faster and hide underneath it (thooft_tower_rates() and
massive_rates() show the gapped rates).

WHAT THIS DOES AND DOES NOT SAY
-------------------------------
At N = infinity and small lambda, the Majorana singlet channels have a
gapless, scale-invariant Rindler mass spectrum, dimensions 2, 3, 5, with
no gap and no discrete masses: the spectrum of a conformal near-horizon
region, consistent with the deconfined stretched horizon of 2506.18054
Sec. 6, and unlike the gapped, discrete 't Hooft tower.  It does not
exclude mesons at finite N or finite lambda (at finite lambda the
spectrum is bounded and correlators end in a t^-3 power law,
late_time_tail()).  The calibration is for free fields probed at a fixed
distance 1/a from the horizon; a = 2 Jc is the natural identification,
not one the proposal supplies.  The Dirac channels are not computed.

Output: horizon_tower.json.  `python horizon_tower.py` prints the tables
and writes the artifact (a few minutes).
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


def scale_fit(delta: int, lam: float, n_omega: int = 160) -> dict:
    """The exact chord density against the ladder at the best-fit scale:
    most of the (Delta/4) lambda deviation at the fixed scale Jc is a
    shift of the effective scale, J_eff/Jc - 1 = -(Delta/4) lambda + ...;
    what is left is the shape deviation, about 0.045 lambda for every
    Delta."""
    from scipy.optimize import minimize_scalar
    from chord import q_of_lambda, s_omega_smooth
    q = q_of_lambda(lam)
    Jc = jc_bridge(lam)
    om = np.linspace(0.0, 6.0 * Jc * math.sqrt(delta), n_omega)
    S = s_omega_smooth(om, q, float(delta), int(min(4000, max(800, 60.0 / lam))))

    def dev(J):
        return float(np.max(np.abs(S - ladder_density(om, J, delta))) / S.max())
    best = minimize_scalar(dev, bounds=(0.7 * Jc, 1.1 * Jc), method="bounded")
    return dict(delta=delta, lam=lam, dev_at_Jc=dev(Jc), J_best_over_Jc=float(best.x / Jc),
                scale_shift_over_lam=float((best.x / Jc - 1.0) / lam),
                dev_at_best=float(best.fun), dev_at_best_over_lam=float(best.fun / lam))


def mass_spectrum_identity(points=((2, 0.0), (2, 0.7), (3, 2.5), (5, 1.0))) -> list:
    """The ladder is a gapless superposition of massive responses:
    int_0^inf z^(2D-1) K_(ix)(z)^2 dz = 2^(2D-3) Gamma(D)^2 |Gamma(D+ix)|^2 / Gamma(2D)
    (the Mellin transform of K_nu^2 at nu = ix), i.e. a Rindler
    Kallen-Lehmann mass density
    rho(M) ~ M^(2 Delta - 1) reaching down to M = 0."""
    import mpmath as mp
    out = []
    with mp.workdps(20):
        for d, x in points:
            lhs = mp.quad(lambda z: z ** (2 * d - 1) * mp.re(mp.besselk(1j * x, z)) ** 2,
                          [0, 1, 5, 20, mp.inf])
            rhs = (2 ** (2 * d - 3) * mp.gamma(d) ** 2 * abs(mp.gamma(d + 1j * x)) ** 2
                   / mp.gamma(2 * d))
            out.append(dict(delta=d, x=x, lhs=float(lhs), rhs=float(rhs),
                            rel=float(abs(lhs - rhs) / abs(rhs))))
    return out


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


def late_time_tail(delta: int, lam: float, nmaxes=(1500, 3000), n_t: int = 500) -> dict:
    """Where the ladder regime ends at finite lambda.  The exact correlator
    (g2_engine) on a log grid up to t = 0.4 nmax/band, for two
    truncations: the time and height of the decay-rate peak (rate/a), and
    the exponent of the late power-law tail (fit over the last decade of
    t).  The density of states has square-root band edges; pairs of
    states near the SAME edge make the spectral density non-analytic at
    omega = 0 (an omega^2 ln|omega| term), which gives a monotone t^-3
    tail.  Pairs at opposite edges add a smaller part oscillating at
    2 E_max."""
    from chord import g2_engine, q_of_lambda
    q = q_of_lambda(lam)
    a = 2.0 * jc_bridge(lam)
    band = 2.0 / math.sqrt(1.0 - q)
    rows = []
    for nmax in nmaxes:
        tmax = 0.4 * nmax / band
        ts = np.geomspace(0.05 / a, tmax, n_t)
        G = np.abs(g2_engine(ts, q, float(delta), nmax).real)
        lg = np.log(G)
        rate = -np.gradient(lg, ts) / a
        i = int(np.argmax(rate))
        sel = ts > tmax / 10.0
        slope = float(np.polyfit(np.log(ts[sel]), lg[sel], 1)[0])
        rows.append(dict(nmax=nmax, a_t_peak=float(a * ts[i]), peak_rate_over_a=float(rate[i]),
                         tail_exponent=slope, min_G=float(G.min())))
    return dict(delta=delta, lam=lam, truncations=rows,
                tail_exponent_drift=abs(rows[0]["tail_exponent"] - rows[-1]["tail_exponent"]))


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
                              * mp.cos(2 * nu * t), [0, 1, 2, 4, 6, 8])
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
        late_time_tail=[late_time_tail(d, lam) for d in (2, 3) for lam in (1.0, 0.5, 0.3)],
        massive_rates=[massive_rates(m) for m in (0.1, 0.5, 1.0)],
        thooft_tower_rates=[thooft_tower_rates(m) for m in (0.1, 0.5, 1.0)],
        bessel_identity=bessel_identity(),
        full_channel_fit=[full_channel_fit(n, lam) for n in (3, 4)
                          for lam in (0.4, 0.2, 0.1, 0.05)],
        massive_fit_control=[massive_fit_control(m, s) for m in (0.1, 0.5, 1.0, 2.0, 3.5, 6.0)
                             for s in ((2, 3), (2, 3, 4))],
        scale_fit=[scale_fit(d, lam) for d in (2, 3, 5) for lam in (0.1, 0.05)],
        mass_spectrum_identity=mass_spectrum_identity(),
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
    print("\nFinite lambda: decay-rate peak and the late power-law tail:")
    for r in p["late_time_tail"]:
        t = r["truncations"][-1]
        print(f"  Delta={r['delta']} lambda={r['lam']}: peak rate/a {t['peak_rate_over_a']:.3f} "
              f"at a*t={t['a_t_peak']:.2f}; tail exponent {t['tail_exponent']:.2f} "
              f"(truncation drift {r['tail_exponent_drift']:.2f})")
    print("\nFull channels as polynomial x ladder (relative residual):")
    for r in p["full_channel_fit"]:
        print(f"  A_{r['n']} lambda={r['lam']}: sectors {r['sectors']}  residual {r['rel_residual']:.2e}")
    for r in p["massive_fit_control"]:
        print(f"  control: free field M/a={r['m_over_a']} on sectors {r['sectors']}: "
              f"residual {r['rel_residual']:.2e} (fixed floor)")
    print("\nNodes on [0, 12a]:", p["nodes"])
    print("Units (MSS string units):", p["units_mss"][0])
    print(f"\nwrote horizon_tower.json ({p['runtime_s']} s)")
