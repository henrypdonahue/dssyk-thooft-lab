#!/usr/bin/env python3
r"""
Boost-to-mass dictionary calibration (road_map rung 13).
========================================================

The DSSYK Hamiltonian generates BOOSTS (paper Sec. 2: H = i d/dt, Rindler time
in the flat-space limit), and every spectrum-match proposal silently assumes
"peaks in M_n boost-frequency correlators = 't Hooft meson masses".  This
module computes what a mass actually looks like in a boost-frequency spectral
function, exactly, in the simplest calibration case: a free scalar particle of
mass m in 1+1d Minkowski vacuum, probed along a uniformly accelerated
worldline of proper acceleration a (boost frequency omega = frequency
conjugate to the observer's proper time; for the trajectory at Rindler radius
1/a this is boost time in units where the horizon generator is rescaled by a).

Closed form (derived by expanding the Wightman function in plane waves,
k = m sinh(theta), and using int du e^{i nu u - i z sinh u} =
2 e^{pi nu / 2} K_{i nu}(z) twice):

    S(omega) = (1/(pi a)) e^{pi omega / a} [ K_{i omega/a}(m/a) ]^2 .

K with imaginary order is real, so S > 0; detailed balance
S(omega) = e^{2 pi omega/a} S(-omega) (the Unruh temperature a/2pi) is
manifest.  The module verifies this closed form against an INDEPENDENT route:
the position-space Wightman function W(tau) = (1/2pi) K_0(m d(tau)) with the
exact worldline geodesic distance continued through the Euclidean section
(d_E = (2/a) sin(a s/2), s = eps + i tau), numerically Fourier transformed --
the epsilon regulator enters the comparison exactly as e^{-|omega| ... } no:
as the known factor e^{-omega eps}.

Calibration facts this establishes (tested in test_boost.py):

1.  A SINGLE discrete mass gives a smooth CONTINUUM in boost frequency -- no
    pole, support on all omega.  Consequence for the duality test: broad
    continua in M_n spectroscopy do NOT falsify a discrete meson tower, and
    peak-hunting at fixed setup is the wrong observable.
2.  The mass enters only through the scaling variable m/a (with one overall
    1/a): S(omega; m, a) = (1/a) s(omega/a, m/a).  Masses must be extracted
    from the SHAPE (the K_{i omega/a}(m/a) kernel -- e.g. via moments or
    kernel inversion), not positions.  A 't Hooft tower {M_n} appears as a
    SUM of such continua: S_tot = sum_n c_n S(omega; M_n, a).
3.  Large mass is Boltzmann-suppressed, S ~ e^{-2m/a} at omega = 0 (the
    K_{0}(x)^2 ~ (pi/2x) e^{-2x} tail): heavy tower members are exponentially
    faint at fixed a.
4.  DOCUMENTED TENSION: the bulk response is KMS-asymmetric at the Unruh
    temperature, but the T_B = infinity hologram correlator
    Tr(M(t)M)/dim is exactly omega-SYMMETRIC (measured in mn_spectroscopy).
    The two can only match through a nontrivial state/frequency map (the
    "tomperature" story of Rahman-Susskind 2401.08555) -- this unresolved map
    is precisely what rung 13 flags as the gating dictionary question, now
    stated sharply.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# closed form
# ---------------------------------------------------------------------------
def besselk_imag_order(nu: float, x: float) -> float:
    """K_{i nu}(x) for real nu, x > 0 -- real-valued; via mpmath."""
    import mpmath as mp
    v = mp.besselk(1j * nu, mp.mpf(x))
    return float(v.real)   # imaginary part is zero to precision


def spectral_closed_form(omega, m: float, a: float = 1.0) -> np.ndarray:
    """S(omega) = (1/(pi a)) e^{pi omega/a} K_{i omega/a}(m/a)^2, vectorized."""
    omega = np.atleast_1d(np.asarray(omega, float))
    out = np.empty_like(omega)
    for i, w in enumerate(omega):
        K = besselk_imag_order(w / a, m / a)
        out[i] = np.exp(np.pi * w / a) * K * K / (np.pi * a)
    return out


# ---------------------------------------------------------------------------
# independent verification route: numeric FT of the K0 Wightman function
# ---------------------------------------------------------------------------
def wightman_tau(tau, m: float, a: float = 1.0, eps: float = 0.05):
    """W_eps(tau) = (1/2pi) K_0( m * (2/a) sin( a (eps + i tau)/2 ) ):
    the Euclidean-section geodesic distance continued to s = eps + i tau.
    Equals int (domega/2pi) S(omega) e^{-omega eps} e^{-i omega tau}."""
    from scipy.special import kv
    s = eps + 1j * np.asarray(tau, float)
    arg = m * (2.0 / a) * np.sin(a * s / 2.0)
    return kv(0, arg) / (2.0 * np.pi)


def spectral_from_wightman(omegas, m: float, a: float = 1.0,
                           eps: float = 0.05, tau_max: float = None,
                           n_tau: int = 2 ** 15) -> np.ndarray:
    """S(omega) recovered as e^{+omega eps} * int dtau e^{+i omega tau}
    W_eps(tau) by direct quadrature -- shares nothing with the closed form."""
    if tau_max is None:
        tau_max = 14.0 / a
    tau = np.linspace(-tau_max, tau_max, n_tau)
    W = wightman_tau(tau, m, a, eps)
    dt = tau[1] - tau[0]
    out = np.empty(len(np.atleast_1d(omegas)))
    for i, w in enumerate(np.atleast_1d(omegas)):
        out[i] = np.real(np.sum(W * np.exp(1j * w * tau)) * dt) * np.exp(w * eps)
    return out


if __name__ == "__main__":
    m, a = 1.0, 1.0
    print("Boost-frequency spectral function of ONE mass m (a = 1):\n")
    print("verification: closed form vs independent K0-Wightman Fourier transform")
    ws = np.array([-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0])
    Sc = spectral_closed_form(ws, m, a)
    Sn = spectral_from_wightman(ws, m, a)
    print(f"{'omega':>7} {'closed form':>13} {'numeric FT':>13} {'rel.diff':>10}")
    for w, sc, sn in zip(ws, Sc, Sn):
        print(f"{w:>7.2f} {sc:>13.6e} {sn:>13.6e} {abs(sn/sc-1):>10.2e}")

    print("\nKMS / detailed balance  S(w) = e^{2 pi w} S(-w):")
    for w in (0.5, 1.0, 2.0):
        lhs = spectral_closed_form([w], m, a)[0]
        rhs = np.exp(2 * np.pi * w / a) * spectral_closed_form([-w], m, a)[0]
        print(f"  w={w}:  ratio = {lhs/rhs:.12f}")

    print("\nA single mass is a smooth continuum -- S at omega = 0..3 for two")
    print("masses (no pole anywhere; mass only reshapes/suppresses):")
    for mm in (1.0, 2.0):
        Ss = spectral_closed_form([0.0, 1.0, 2.0, 3.0], mm, a)
        print(f"  m={mm}: " + "  ".join(f"{s:.3e}" for s in Ss))
    print(f"\n  Boltzmann suppression at omega=0: S(m=2)/S(m=1) = "
          f"{spectral_closed_form([0.0],2,a)[0]/spectral_closed_form([0.0],1,a)[0]:.4f}"
          f"   (vs e^-2(m2-m1)/a = {np.exp(-2.0):.4f} asymptotically)")

    print("\nScaling collapse S(w; m, a) = (1/a) s(w/a, m/a):")
    for (mm, aa, w) in [(1.0, 1.0, 1.0), (2.0, 2.0, 2.0), (3.0, 3.0, 3.0)]:
        print(f"  a={aa}: a*S(a*w0) = "
              f"{aa*spectral_closed_form([w], mm, aa)[0]:.10e}")
