#!/usr/bin/env python3
"""Tests for the boost-to-mass dictionary calibration (rung 13)."""

import numpy as np
import pytest

from boost_mass import (spectral_closed_form, spectral_from_wightman,
                        besselk_imag_order)


def test_closed_form_vs_independent_ft():
    """The derived S(w) = (1/pi a) e^{pi w/a} K_{iw/a}(m/a)^2 must match the
    independent numeric Fourier transform of the K0 Wightman function."""
    ws = np.array([-1.0, 0.0, 0.7, 2.0])
    Sc = spectral_closed_form(ws, m=1.0, a=1.0)
    Sn = spectral_from_wightman(ws, m=1.0, a=1.0)
    assert np.all(np.abs(Sn / Sc - 1.0) < 1e-8)


def test_kms_detailed_balance():
    """S(w) = e^{2 pi w / a} S(-w): the Unruh temperature a/2pi, exact."""
    for w in (0.3, 1.0, 2.5):
        lhs = spectral_closed_form([w], 1.0, 1.0)[0]
        rhs = np.exp(2 * np.pi * w) * spectral_closed_form([-w], 1.0, 1.0)[0]
        assert abs(lhs / rhs - 1.0) < 1e-10


def test_scaling_collapse():
    """S(w; m, a) = (1/a) s(w/a, m/a): mass only enters via m/a."""
    ref = 1.0 * spectral_closed_form([1.0], 1.0, 1.0)[0]
    for aa in (2.0, 3.0):
        val = aa * spectral_closed_form([aa * 1.0], aa * 1.0, aa)[0]
        assert abs(val / ref - 1.0) < 1e-8


def test_single_mass_is_a_continuum_with_nodes():
    """Calibration fact 1: one discrete mass -> positive continuum, no
    isolated support (every unit interval carries weight), but with
    oscillatory near-nodes at large w (zeros of K_{iw}(m))."""
    ws = np.linspace(0.0, 8.0, 161)
    S = spectral_closed_form(ws, 1.0, 1.0)
    assert np.all(S >= 0)
    for lo in range(0, 8):
        seg = S[(ws >= lo) & (ws <= lo + 1)]
        assert seg.max() > 1e-3          # support everywhere
    assert S.min() < 1e-2 * S.max()      # ...but genuine near-nodes exist


def test_nodes_encode_the_mass():
    """Calibration fact 2: the first spectral node moves with the mass
    (zeros of K_{iw}(m/a) in w) -- the mass lives in the node/oscillation
    pattern, not in a peak position."""
    def first_node(m):
        ws = np.linspace(1.5, 6.0, 226)
        S = spectral_closed_form(ws, m, 1.0)
        for i in range(1, len(ws) - 1):
            if S[i] < S[i - 1] and S[i] < S[i + 1] and S[i] < 1e-2:
                return ws[i]
        raise AssertionError("no node found")
    n1, n2 = first_node(1.0), first_node(2.0)
    assert n2 > n1 + 0.5                  # m=1: ~3.0; m=2: ~4.4


def test_boltzmann_suppression_of_heavy_masses():
    """Calibration fact 3: heavier masses are exponentially fainter at fixed
    a (K_0(x)^2 ~ e^{-2x} tail): S(0; 2m) << S(0; m)."""
    s1 = spectral_closed_form([0.0], 1.0, 1.0)[0]
    s2 = spectral_closed_form([0.0], 2.0, 1.0)[0]
    s3 = spectral_closed_form([0.0], 3.0, 1.0)[0]
    assert s2 < 0.2 * s1
    assert s3 < 0.2 * s2


def test_besselk_imaginary_order_is_real():
    """K_{i nu}(x) is real for real nu, x > 0 (sanity of the mpmath path)."""
    import mpmath as mp
    for nu, x in [(0.5, 1.0), (2.0, 0.3), (5.0, 2.0)]:
        v = mp.besselk(1j * nu, mp.mpf(x))
        assert abs(float(v.imag)) < 1e-15 * max(1.0, abs(float(v.real)))
        assert besselk_imag_order(nu, x) == pytest.approx(float(v.real))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
