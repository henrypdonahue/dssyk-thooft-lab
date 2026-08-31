#!/usr/bin/env python3
"""Tests for the derived large-q ladder-kernel skeleton.

Four things: the derived-vs-transcribed gate (the point of the module),
the even/odd classification and its limits, the special poles m = 0 and
m = +-v, and the CMS closed-form cross-checks.  Every headline number in
largeq_kernel.json is pinned here.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from largeq_anchor import quantization_sets
from largeq_kernel import (M_MAX, cms_norm, cms_psi_e, cms_psi_o,
                           compare_at_v, derived_quantization_sets,
                           doubling_shift, eigenfunction_check, ode_residual,
                           scramblon_mode, zero_mode)

_JSON = Path(__file__).resolve().parent / "largeq_kernel.json"


# ---------------------------------------------------------------------------
# the gate: derived reproduces transcribed
# ---------------------------------------------------------------------------
def test_derived_matches_transcribed():
    """The derived k = 1 spectrum must reproduce the CMS (3.3)
    transcription at <= 1e-8 relative (the gate; the achieved level,
    ~1e-13 at these K, is pinned via the json test below)."""
    for v, K in ((0.05, 160), (0.35, 160), (0.65, 160), (0.95, 256)):
        row = compare_at_v(v, K=K)
        assert row["max_rel_dev_e"] < 1e-8, (v, row["max_rel_dev_e"])
        assert row["max_rel_dev_o"] < 1e-8, (v, row["max_rel_dev_o"])


def test_even_odd_classification():
    """The parity blocks land on the right transcribed sets: swapping the
    labels must FAIL by O(1) while the correct pairing agrees to 1e-8.
    (M^e and M^o interleave, so the swap error is about one unit.)"""
    v = 0.5
    d_e, d_o = derived_quantization_sets(v, K=160)
    t_e, t_o = quantization_sets(v, m_max=M_MAX)
    assert max(abs(a - b) for a, b in zip(d_e, t_e)) < 1e-8
    assert max(abs(a - b) for a, b in zip(d_o, t_o)) < 1e-8
    swapped = min(abs(a - b) for a, b in zip(d_e, t_o))
    assert swapped > 0.3          # the labels are not interchangeable


def test_infinite_T_degeneration_of_derived_sets():
    """v -> 0: the DERIVED M^e -> odd integers and M^o -> even integers
    (the transcription-side version of this is in test_largeq.py; here it
    holds for the eigensolver output with no root-finding involved)."""
    d_e, d_o = derived_quantization_sets(0.05, m_max=8.0, K=160)
    for k, m in enumerate(d_e[:4]):
        assert m == pytest.approx(2 * k + 1, abs=0.02)
    for k, m in enumerate(d_o[:3]):
        assert m == pytest.approx(2 * k + 2, abs=0.02)


def test_derived_sets_interleave():
    """The CP-alternation slot survives the derivation: e/o/e/o... at
    intermediate coupling."""
    d_e, d_o = derived_quantization_sets(0.5, K=160)
    merged = sorted([(m, "e") for m in d_e] + [(m, "o") for m in d_o])
    kinds = [k for _, k in merged]
    assert kinds[0] == "e"
    assert all(a != b for a, b in zip(kinds, kinds[1:]))


def test_grid_doubling_convergence():
    """Spectral convergence: doubling K moves no root by more than 1e-9
    at the worst scan point v = 0.95 (the K -> 2K shifts across the whole
    scan are recorded per v in the json)."""
    assert doubling_shift(0.95, K=320) < 1e-9


# ---------------------------------------------------------------------------
# the special poles m = 0 and m = +-v [CMS (3.13)]
# ---------------------------------------------------------------------------
def test_special_modes_solve_the_ode():
    """tan(z/2) and sec(z/2) are EXACT solutions at eigenvalue 0 (m = 0)
    and 1/4 (m = +-v): checked against the ODE by finite differences, not
    against their own hand-derived derivatives."""
    for v in (0.3, 0.6, 0.9):
        assert ode_residual(zero_mode, 0.0, v) < 1e-7
        assert ode_residual(scramblon_mode, 0.25, v) < 1e-7


def test_special_modes_violate_dirichlet_and_stay_off_the_sets():
    """The m = 0 and m = +-v slots are Regge poles, not tower members:
    both closed-form solutions are nonzero at the boundary, and every
    derived root sits strictly above v."""
    for v in (0.3, 0.6, 0.9):
        L = v * np.pi
        assert abs(zero_mode(-L)) > 0.4        # tan(v pi / 2) != 0
        assert abs(scramblon_mode(-L)) > 1.0   # sec(v pi / 2) > 1
        d_e, d_o = derived_quantization_sets(v, K=160)
        assert min(d_e + d_o) > v + 0.05


def test_m_equals_v_is_a_formal_root_with_trivial_eigenfunction():
    """m = v solves the M^o root condition EXACTLY (fo(v) = 0), yet the
    CMS (3.2) odd eigenfunction is identically zero there -- the mode
    carries nothing, which is why the transcription's m > v restriction
    is a derivation, not a convention.  The even branch's trivial slot is
    m = 0: psi^e_0 == 0."""
    v = 0.6
    c = v * np.tan(np.pi * v / 2.0)
    fo = lambda m: m * np.tan(np.pi * m / 2.0) - c      # the M^o condition
    assert abs(fo(v)) < 1e-14                            # formal root at m = v
    # ... yet the nearest ACTUAL root sits a finite distance above v
    t_e, t_o = quantization_sets(v, m_max=M_MAX)
    assert min(t_e + t_o) - v > 0.5
    x = np.linspace(0.3, 5.9, 41)
    assert np.max(np.abs(cms_psi_o(v, v, x))) < 1e-13
    assert np.max(np.abs(cms_psi_e(0.0, v, x))) < 1e-13
    # ... and m = 0 is NOT a root of the even condition: the m -> 0 limit
    # of m cot(m pi/2) is 2/pi, so fe(0) = 2/pi + v tan(v pi/2) > 0
    assert 2 / np.pi + v * np.tan(np.pi * v / 2) > 0.6


# ---------------------------------------------------------------------------
# CMS closed-form cross-checks
# ---------------------------------------------------------------------------
def test_eigenfunctions_match_cms_closed_forms():
    """The numerical eigenvectors reproduce the CMS (3.2) closed forms,
    and their x-space norms reproduce the CMS (3.4) closed form."""
    rows = eigenfunction_check(0.6, K=160, n_roots=2)
    assert len(rows) == 4
    for r in rows:
        assert r["max_abs_dev"] < 1e-6, r
        assert r["norm_rel_dev"] < 1e-6, r


def test_cms_printed_example_v_3_over_5():
    """The one table CMS print (below their (3.3)), at v = 3/5, from the
    DERIVED spectrum: M^e = {1.35, 3.16, 5.10, 7.07, 9.06},
    M^o = {2.23, 4.13, 6.09, 8.06, 10.05}."""
    d_e, d_o = derived_quantization_sets(0.6, m_max=10.5, K=160)
    assert [round(m, 2) for m in d_e[:5]] == [1.35, 3.16, 5.10, 7.07, 9.06]
    assert [round(m, 2) for m in d_o[:5]] == [2.23, 4.13, 6.09, 8.06, 10.05]


def test_cms_norm_closed_form_at_a_root():
    """CMS (3.4) against direct quadrature at the first even root of
    v = 0.6 (root from the transcription; the norm formula only holds ON
    the quantization set)."""
    v = 0.6
    m = quantization_sets(v, m_max=3.0)[0][0]
    x = np.linspace(0.0, 2 * np.pi, 20001)
    num = np.trapezoid(cms_psi_e(m, v, x) ** 2, x)
    assert num == pytest.approx(cms_norm(m, v, "e"), rel=1e-6)


# ---------------------------------------------------------------------------
# the campaign record
# ---------------------------------------------------------------------------
# largeq_kernel.json is committed: missing = broken checkout, fail not skip
def test_json_headline_numbers():
    """The recorded scan achieved what the module claims: every v in
    {0.05, ..., 0.95}, max relative deviation below 1e-11 (three decades
    under the 1e-8 gate), doubling shifts below 1e-10, and per-row
    consistency of the recorded maxima."""
    payload = json.loads(_JSON.read_text())
    assert payload["v_scan"] == [round(0.05 * i, 2) for i in range(1, 20)]
    ach = payload["achieved"]
    assert ach["max_rel_dev_overall"] < 1e-11
    assert ach["max_doubling_shift"] < 1e-10
    per_v = payload["per_v"]
    assert len(per_v) == 19
    assert ach["max_rel_dev_overall"] == max(
        max(r["max_rel_dev_e"], r["max_rel_dev_o"]) for r in per_v)
    for r in per_v:
        for br in ("e", "o"):
            assert len(r[f"derived_{br}"]) == len(r[f"transcribed_{br}"])
    sp = payload["special_poles"]
    assert sp["min_root_over_v"] > 1.0
    assert sp["residual_zero_mode"] < 1e-7
    assert sp["residual_scramblon"] < 1e-7
    assert sp["psi_e_trivial_at_m0"] < 1e-13
    assert sp["psi_o_trivial_at_mv"] < 1e-13
    for r in payload["eigenfunction_check"]:
        assert r["max_abs_dev"] < 1e-7      # achieved 1.35e-8 (v=0.9 top root)
        assert r["norm_rel_dev"] < 1e-11    # achieved 5.8e-13
    ex = payload["cms_printed_example"]
    assert ex["Me_2dp"] == ex["cms_printed_Me"]
    assert ex["Mo_2dp"] == ex["cms_printed_Mo"]


def test_json_reproducible_at_one_v():
    """Determinism: re-deriving v = 0.5 at the producer's K reproduces the
    recorded roots to 1e-10 (the slack is BLAS thread-count rounding, not
    physics; the discretization error itself is another decade down)."""
    payload = json.loads(_JSON.read_text())
    K = payload["method"]["K"]
    row = next(r for r in payload["per_v"] if r["v"] == 0.5)
    d_e, d_o = derived_quantization_sets(0.5, K=K)
    cut = M_MAX - 1e-6
    d_e = [m for m in d_e if m < cut]
    d_o = [m for m in d_o if m < cut]
    assert np.allclose(d_e, row["derived_e"], rtol=0, atol=1e-10)
    assert np.allclose(d_o, row["derived_o"], rtol=0, atol=1e-10)


def test_derived_matches_transcribed_at_producer_precision():
    """The producer's K = 400 gate, re-derived live at three v: the
    derived sets match the transcribed CMS roots to 1e-11 relative (the
    stored full-scan number is 3e-12), so the headline digits are not a
    stored scalar the suite merely reads back."""
    for v in (0.2, 0.5, 0.8):
        row = compare_at_v(v, K=400)
        assert row["max_rel_dev_e"] < 1e-11, (v, row["max_rel_dev_e"])
        assert row["max_rel_dev_o"] < 1e-11, (v, row["max_rel_dev_o"])


def test_v_outside_open_interval_is_rejected():
    for v in (0.0, 1.0, -0.2, 1.5):
        with pytest.raises(ValueError):
            derived_quantization_sets(v)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
