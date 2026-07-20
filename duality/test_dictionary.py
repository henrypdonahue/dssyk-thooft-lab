#!/usr/bin/env python3
"""Tests for the comparison dictionary: every convention that could silently
wreck a DSSYK <-> 't Hooft comparison is pinned here."""

import math

import numpy as np
import pytest

from dictionary import (M1_COEFFICIENT, lambda_paper, lambda_chord,
                        chord_from_paper, paper_from_chord, q_parameter,
                        alpha_bar, gbar_squared, string_tension,
                        n_paper_from_flz, n_flz_from_paper, cp_paper,
                        flz_parity, load_reference, mass_squared_ratios,
                        interleave_splittings)


def test_lambda_conventions():
    """The factor-2 landmine: lambda_chord = 2 * lambda_paper, q = e^-lambda_chord."""
    N, p = 640, 26
    assert lambda_paper(N, p) == pytest.approx(26 * 26 / 640)
    assert lambda_chord(N, p) == pytest.approx(2 * 26 * 26 / 640)
    assert chord_from_paper(lambda_paper(N, p)) == pytest.approx(lambda_chord(N, p))
    assert paper_from_chord(lambda_chord(N, p)) == pytest.approx(lambda_paper(N, p))
    assert q_parameter(lambda_chord(N, p)) == pytest.approx(
        math.exp(-2 * 26 * 26 / 640))
    # q -> 1 exactly as lambda -> 0 (the corner where the match should be exact)
    assert q_parameter(chord_from_paper(1e-9)) == pytest.approx(1.0, abs=1e-8)


def test_coupling_map():
    """alpha_bar = gbar^2 N = p^2 and gbar^2 = lambda (paper Sec. 4.1)."""
    N, p = 100, 8
    assert alpha_bar(p) == 64
    assert gbar_squared(N, p) == pytest.approx(0.64)
    assert gbar_squared(N, p) * N == pytest.approx(alpha_bar(p))
    # string tension tau = p^2 m_q^2; the duality point is m_q^2 -> 0
    assert string_tension(0.0, p) == 0.0
    assert string_tension(2.5, p) == pytest.approx(64 * 2.5)


def test_index_offset():
    """Correction #1: n_paper = n_FLZ + 2, and n_paper < 2 has no FLZ level."""
    assert n_paper_from_flz(0) == 2
    assert n_paper_from_flz(7) == 9
    assert n_flz_from_paper(2) == 0
    for bad in (0, 1):
        with pytest.raises(ValueError):
            n_flz_from_paper(bad)


def test_cp_assignments():
    """CP = (-1)^(n+1): photon -1, graviton +1, eta (n=2) -1, f1 (n=3) +1;
    and the FLZ-parity <-> CP correspondence is consistent through the offset."""
    assert cp_paper(0) == -1   # photon
    assert cp_paper(1) == +1   # graviton
    assert cp_paper(2) == -1   # eta
    assert cp_paper(3) == +1   # f1
    for n_flz in range(20):
        # symmetric FLZ levels must be CP = -1, antisymmetric CP = +1
        expected = -1 if flz_parity(n_flz) == "sym" else +1
        assert cp_paper(n_paper_from_flz(n_flz)) == expected


def test_m1_coefficient():
    """Correction #3 as a constant: M_1 = -(p/2) i H."""
    assert M1_COEFFICIENT(4) == -2.0
    assert M1_COEFFICIENT(8) == -4.0


def test_reference_loads_and_matches_flz():
    """The decorated reference table must reproduce the FLZ anchors."""
    rows = load_reference()
    assert rows[0]["two_lambda"] == pytest.approx(0.737061746292690, abs=1e-14)
    assert rows[1]["two_lambda"] == pytest.approx(1.7537313369175, abs=1e-12)
    for r in rows:
        assert r["n_paper"] == r["n_flz"] + 2
        assert r["CP"] == cp_paper(r["n_paper"])
        assert r["parity"] == flz_parity(r["n_flz"])


def test_mass_ratios_parameter_free():
    """Fingerprint #1: M_n^2/M_eta^2 ratios; first excited/ground ~ 2.3794."""
    ratios = mass_squared_ratios()
    assert ratios[0]["ratio"] == pytest.approx(1.0)
    assert ratios[1]["ratio"] == pytest.approx(2.379354709, abs=1e-8)
    # monotone increasing along the trajectory
    vals = [r["ratio"] for r in ratios]
    assert all(b > a for a, b in zip(vals, vals[1:]))


def test_interleave_fingerprint():
    """Fingerprint #2: strictly interleaved alternating-CP tower with
    near-unit gaps converging to 1 (2lam ~ n + 3/4), NOT degenerate doublets."""
    gaps = interleave_splittings()
    for g in gaps:
        assert g["CP_pair"][0] == -g["CP_pair"][1]      # CP alternates
        assert 0.9 < g["gap"] < 1.1                      # near-unit, never ~0
    assert gaps[0]["gap"] == pytest.approx(1.0166695906, abs=1e-9)
    # |gap - 1| shrinks with n (both subsequences approach unit spacing)
    devs = [abs(g["gap"] - 1.0) for g in gaps]
    assert devs[-1] < devs[0]
    assert max(devs[-4:]) < 0.001


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
