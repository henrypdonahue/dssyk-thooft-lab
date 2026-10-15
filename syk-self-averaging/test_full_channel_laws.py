"""Tests for full_channel_laws.py: the exact fixed-p full-channel scans,
the peak, and the certified law N w2 -> 8p(p-2)."""

import json
from fractions import Fraction
from pathlib import Path

import pytest

from full_channel_laws import (closed_form_A0, closed_form_A2, closed_form_L,
                               full_law, full_moments, sigma_h2, width_law)

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def art():
    with open(HERE / "full_channel_laws.json") as fh:
        return json.load(fh)


def test_scan_rows_recompute(art):
    """Cheap rows of every scan are recomputed exactly."""
    for key, s in art["scans"].items():
        p, n = s["p"], s["n"]
        for r in s["rows"][:3]:
            f = full_moments(r["N"], p, n)
            assert abs(float(f[2] / f[0]) - r["w2"]) < 1e-13 * abs(r["w2"])
            assert abs(float(f[2] / f[0] * sigma_h2(r["N"], p))
                       - r["w2_physical"]) < 1e-12 * abs(r["w2_physical"])


def test_peak_then_decay(art):
    """The full width peaks inside lambda_chord in [1.2, 1.7] and then
    falls monotonically; N w2 approaches 8p(p-2) from below."""
    for key, s in art["scans"].items():
        rows = sorted(s["rows"], key=lambda r: r["N"])
        w2 = [r["w2"] for r in rows]
        ipk = max(range(len(w2)), key=lambda i: w2[i])
        assert 1.2 <= rows[ipk]["lam_chord"] <= 1.7, (key, rows[ipk])
        assert all(w2[i + 1] < w2[i] for i in range(ipk, len(w2) - 1)), key
        if s["n"] == 2:
            L = closed_form_L(s["p"])
            Nw = [r["N_w2"] for r in rows[ipk:]]
            assert all(Nw[i + 1] > Nw[i] for i in range(len(Nw) - 1))
            assert Nw[-1] < L and Nw[-1] > 0.9 * L, (key, Nw[-1], L)


def test_bench_window_is_the_rising_side(art):
    """The ED bench's p = 4 window (N <= 20) sits entirely before the peak."""
    s = art["scans"]["p4_n2"]
    assert s["peak_N"] == 24
    rows = sorted((r for r in s["rows"] if r["N"] <= 24), key=lambda r: r["N"])
    w2 = [r["w2"] for r in rows]
    assert all(w2[i + 1] > w2[i] for i in range(len(w2) - 1))


def test_width_laws_certified(art):
    """N w2 -> 8p(p-2), mu0 -> 64p^2(p-1)^2/N^2, mu2 -> 512p^3(p-1)^2(p-2)/N^3:
    the stored laws match the closed forms at every p, the closed forms
    were guessed on p = 2, 4, 6 and checked at 8, 10, 12, and the p = 4
    law is recomputed live."""
    assert art["L_identified"] and art["coefficients_identified"]
    for w in art["width_laws"]:
        p = w["p"]
        assert Fraction(w["L"]) == closed_form_L(p)
        assert Fraction(w["mu0_law"]["A"]) == closed_form_A0(p)
        assert Fraction(w["mu2_law"]["A"]) == closed_form_A2(p)
        assert w["mu0_law"]["deficit"] == 2
        assert w["mu2_law"]["deficit"] == (3 if p > 2 else 6)
    assert art["L_checks"] == [8, 10, 12]
    live = width_law(4, 2)
    assert Fraction(live["L"]) == 64
    assert live["mu2_law"]["A"] == art["width_laws"][1]["mu2_law"]["A"]


def test_physical_width_is_N_independent(art):
    """w2 x <H^2>/dim -> 8p(p-2) in the paper's units (J = 1): the finite
    O(J) width at fixed p."""
    for key in ("p4_n2", "p6_n2"):
        s = art["scans"][key]
        p = s["p"]
        last = sorted(s["rows"], key=lambda r: r["N"])[-1]
        assert abs(last["w2_physical"] - 8 * p * (p - 2)) < 0.12 * 8 * p * (p - 2)
        assert last["w2_physical"] < 8 * p * (p - 2)


def test_p2_is_the_free_point():
    """p = 2: the full A_2 is central (conserved_points.py), so mu_2 vanishes
    at every N and the law's deficit jumps."""
    l2 = full_law(2, 2, 2)
    assert Fraction(l2["A"]) == 0
    for N in (6, 8, 10):
        assert full_moments(N, 2, 2)[2] == 0
