#!/usr/bin/env python3
"""Regression tests for the massless 't Hooft spectrum.

Run with: pytest -q, or python3 -m pytest -q.

Every check that validate.py PRINTS, this file ASSERTS, with an explicit
tolerance.  So the package fails loudly if a change breaks the physics.

The whole suite is double precision and runs in a few seconds.  The sine
cross-check quadrature dominates the cost.
"""

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import zeta, polygamma

from thooft_spectrum import solve_sector, spectrum, _I_matrix
import sine_solver

PI = np.pi

# FLZ Tables 1 & 2 anchors: one provenance-checked copy (flz_reference.py)
from flz_reference import FLZ, FLZ_PADE_LAMBDA0, sum_inv_sq_with_tail


@pytest.fixture(scope="module")
def sectors():
    """Both CP sectors at a fixed, generous basis size (alpha = 0)."""
    ev_sym, _, _ = solve_sector(300, parity=0)
    ev_anti, _, _ = solve_sector(300, parity=1)
    return ev_sym, ev_anti


def test_flz_tables(sectors):
    """Reproduce FLZ Tables 1 & 2 to high precision for n <= 29."""
    ev_sym, ev_anti = sectors
    worst = 0.0
    for n, ref in FLZ.items():
        mine = ev_sym[n // 2] if n % 2 == 0 else ev_anti[n // 2]
        worst = max(worst, abs(mine - ref))
    assert worst < 1e-10, f"worst FLZ deviation {worst:.2e} exceeds 1e-10"


def test_ground_state_matches_pade(sectors):
    ev_sym, _ = sectors
    assert abs(ev_sym[0] - FLZ_PADE_LAMBDA0) < 1e-12


def test_sum_rules():
    """FLZ Eq. (1.8):  sum 1/lambda_2m^2 = 7 zeta(3),  sum 1/lambda_2m+1^2 = 2."""
    ev_sym, _, _ = solve_sector(400, parity=0)
    ev_anti, _, _ = solve_sector(400, parity=1)
    n_cut = 120

    assert abs(sum_inv_sq_with_tail(ev_sym, 0, n_cut) - 7 * zeta(3)) < 1e-8
    assert abs(sum_inv_sq_with_tail(ev_anti, 1, n_cut) - 2.0) < 1e-8


def test_variational_from_above():
    """Rayleigh-Ritz: the ground state must decrease monotonically toward, and
    stay above, the true value as the basis grows."""
    prev = None
    for K in [4, 8, 16, 32, 64, 128]:
        ev, _, _ = solve_sector(K, parity=0)
        v = ev[0]
        assert v >= FLZ_PADE_LAMBDA0 - 1e-12, "eigenvalue fell below the truth"
        if prev is not None:
            assert v <= prev + 1e-14, "eigenvalue not monotonically decreasing"
        prev = v


def test_asymptotics():
    """2*lambda_n -> n + 3/4 ('t Hooft/WKB); the gap must shrink with n."""
    levels = spectrum(num_basis=300, num_levels=30)
    gaps = [abs(L["two_lambda"] - (L["n"] + 0.75)) for L in levels[8:]]
    assert gaps[-1] < 3e-5
    assert gaps[-1] < gaps[0]  # shrinking toward the asymptote


def test_cp_assignment():
    """Symmetric (even-degree) levels are FLZ even n with CP = -1, and the
    two sectors strictly interleave.  interleave_levels sets CP and parity
    from n by definition, so the content is checked against the sector
    solves themselves: level n of the merged list IS the (n//2)-th
    eigenvalue of the parity-(n%2) sector, and the sectors alternate."""
    levels = spectrum(num_basis=100, num_levels=10)
    ev_sym, _, _ = solve_sector(100, parity=0)
    ev_anti, _, _ = solve_sector(100, parity=1)
    for L in levels:
        n = L["n"]
        src = ev_sym if n % 2 == 0 else ev_anti
        assert L["two_lambda"] == pytest.approx(src[n // 2], rel=1e-14)
        assert L["CP"] == (-1) ** (n + 1)
        assert L["parity"] == ("even" if n % 2 == 0 else "odd")
    vals = [L["two_lambda"] for L in levels]
    assert all(b > a for a, b in zip(vals, vals[1:]))


def _U(n):
    """Chebyshev U_n(x) as a numpy-callable, stable at |x| -> 1."""
    def f(x):
        x = np.clip(x, -1.0, 1.0)
        th = np.arccos(x)
        s = np.sin(th)
        near = s < 1e-9
        out = np.where(near, (n + 1) * np.sign(x) ** n, 0.0)
        safe = np.where(near, 1.0, s)
        return np.where(near, out, np.sin((n + 1) * th) / safe)
    return f


@pytest.mark.parametrize("degrees", [[0, 2, 4, 6], [1, 3, 5, 7]])
def test_I_matrix_vs_quadrature(degrees):
    """Closed-form I_{nm} = Int_{-1}^1 U_n U_m dxi must match numerical quad."""
    degrees = np.array(degrees)
    I = _I_matrix(degrees)
    for i, n in enumerate(degrees):
        for j, m in enumerate(degrees):
            ref, _ = quad(lambda x: _U(n)(x) * _U(m)(x), -1, 1, limit=200)
            assert abs(I[i, j] - ref) < 1e-10


def test_potential_matrix_vs_quadrature():
    """The full potential matrix element  <b_n,(1/x+1/(1-x)) b_m>  (b = sqrt(1-xi^2)U)
    must equal 4 I_{nm}.  An independent 1/x-form quadrature checks it."""
    degrees = np.array([0, 2, 4])
    P = 4.0 * _I_matrix(degrees)  # what the solver adds (per unit alpha)
    for i, n in enumerate(degrees):
        for j, m in enumerate(degrees):
            def integrand(xi):
                x = 0.5 * (1 + xi)
                b = (1 - xi ** 2)  # sqrt(1-xi^2)^2, the two basis weights
                return b * _U(n)(xi) * _U(m)(xi) * (1.0 / x + 1.0 / (1.0 - x))
            ref, _ = quad(integrand, -1, 1, limit=400, points=[0.0])
            assert abs(P[i, j] - ref) < 1e-8


# --------------------------------------------------------------------------
# External anchors (external_anchors.json): FLZ higher sum rules and the
# Litvinov-Meshcheriakov alpha != 0 tables / exact alpha-dependent sum rules
# --------------------------------------------------------------------------
import math
import os as _os
from math import factorial as _factorial

with open(_os.path.join(_os.path.dirname(__file__), "external_anchors.json")) as _fh:
    ANCHORS = json.load(_fh)


def _G_sum(evs_two_lambda, offset, s, n_cut=120):
    """G^(s) = sum 1/lambda_n^s over one parity sector, with the analytic
    polygamma tail from the leading asymptotic 2*lambda_n -> n + 3/4."""
    total, n = 0.0, offset
    for tl in evs_two_lambda:
        if n > n_cut:
            break
        total += 2.0 ** s / tl ** s
        n += 2
    a = (n + 0.75) / 2.0
    total += (-1) ** s * polygamma(s - 1, a) / _factorial(s - 1)
    return total


def test_flz_higher_sum_rules():
    """FLZ Eq. (1.8) beyond s = 2: G3+- and G4+- against the exact zeta-value
    closed forms.  These are four more independent whole-spectrum
    certificates."""
    ev_sym, _, _ = solve_sector(400, parity=0)
    ev_anti, _, _ = solve_sector(400, parity=1)
    targets = {
        (3, 0): -4 * PI ** 2 / 3 + 28 * zeta(3),
        (3, 1): -8.0 / 3.0 + 4 * PI ** 2 / 9,
        (4, 0): (-2 * PI ** 2 + 42 * zeta(3) - (7 / 3) * PI ** 2 * zeta(3)
                 + (49 / 2) * zeta(3) ** 2 + (31 / 2) * zeta(5)),
        (4, 1): 11.0 / 3.0 - (7 / 9) * PI ** 2 + (7 / 6) * PI ** 2 * zeta(3)
                - (31 / 4) * zeta(5),
    }
    for (s, parity), tgt in targets.items():
        evs = ev_sym if parity == 0 else ev_anti
        got = _G_sum(evs, parity, s)
        assert abs(got - tgt) < 1e-10, f"G{s} parity {parity}: {got} vs {tgt}"
    # and the tabulated FLZ Table 3 values agree with those closed forms
    assert abs(targets[(3, 0)] - ANCHORS["flz_table3_alpha0"]["3"][0]) < 1e-12
    assert abs(targets[(4, 1)] - ANCHORS["flz_table3_alpha0"]["4"][1]) < 1e-12


def test_flz_table3_all_sixteen_anchors():
    """Every FLZ Table 3 value at alpha = 0 -- G^(s)+- for s = 1..8, sixteen
    numbers, ten of them previously stored but unconsumed -- is reproduced
    from the computed spectrum with the Eq. (1.6) tail, at the level of
    the table's own last printed digit (~1e-12 absolute; the s = 7, 8
    G+ entries are printed to 1e-10).  The s = 1 case is FLZ's
    regularized sum (their Eq. 4.18), whose closed forms log(8 pi) - 1 and
    log(8 pi) - 3 are checked against the table too."""
    from flz_reference import flz_sum_rule
    ev_sym, _, _ = solve_sector(400, parity=0)
    ev_anti, _, _ = solve_sector(400, parity=1)
    table = ANCHORS["flz_table3_alpha0"]
    for s in range(1, 9):
        gp, gm = table[str(s)]
        # tolerance = one unit in the table's last printed digit (14
        # significant figures; the s = 8 entry is truncated, not rounded)
        # plus the solver's ~2e-12 floor
        tol_p = 10.0 ** (math.floor(math.log10(abs(gp))) - 13) + 3e-12
        got_p = flz_sum_rule(ev_sym, 0, s, tail="eq16")
        got_m = flz_sum_rule(ev_anti, 1, s, tail="eq16")
        assert abs(got_p - gp) < tol_p, (s, "+", got_p, gp, tol_p)
        assert abs(got_m - gm) < 5e-12, (s, "-", got_m, gm)
    # the two s = 1 closed forms of Eq. 1.8
    assert abs(math.log(8 * PI) - 1 - table["1"][0]) < 1e-13
    assert abs(math.log(8 * PI) - 3 - table["1"][1]) < 1e-13
    # the Eq. (1.6) tail is what buys the last two digits at s = 2 and the
    # last three at s = 1
    assert abs(flz_sum_rule(ev_sym, 0, 2, tail="leading") - 7 * zeta(3)) < 2e-10
    assert abs(flz_sum_rule(ev_sym, 0, 2, tail="eq16") - 7 * zeta(3)) < 5e-12
    assert abs(flz_sum_rule(ev_sym, 0, 1, tail="leading") - table["1"][0]) < 1e-8


def test_lm_alpha_nonzero_eigenvalues():
    """External anchor at alpha != 0: the spectrum must match
    the published Litvinov-Meshcheriakov tables (arXiv:2409.11324).  Their 6
    significant digits in lambda give ~1e-5 on 2*lambda for the low levels."""
    for a_str, tabs in ANCHORS["lm_two_lambda_tables"].items():
        if a_str == "comment":
            continue
        alpha = float(a_str)
        levels = spectrum(num_basis=300, num_levels=10, alpha=alpha)
        for n, L in enumerate(levels):
            ref = tabs["even"][n // 2] if n % 2 == 0 else tabs["odd"][n // 2]
            assert abs(L["two_lambda"] - ref) < 3e-5, (
                f"alpha={alpha} n={n}: {L['two_lambda']:.6f} vs LM {ref}")


def test_lm_exact_alpha_sum_rules():
    """Sharper alpha != 0 anchor: the DIFFERENCE G2(alpha) - G2(0) converges
    absolutely (no alpha-dependent tail needed).  The exact 20-digit
    Litvinov-Meshcheriakov sum-rule values therefore test the whole alpha != 0
    spectrum against exact analytic results.  The basis' hardwired alpha = 0
    endpoint exponent limits the accuracy here (algebraic convergence off the
    duality point).  The tolerance below IS that measured limitation."""
    ev0_s, _, _ = solve_sector(300, parity=0)
    ev0_a, _, _ = solve_sector(300, parity=1)
    G2p0, G2m0 = 7 * zeta(3), 2.0
    for a_str in ("-0.5", "-0.1", "0.1", "0.5"):
        alpha = float(a_str)
        eva_s, _, _ = solve_sector(300, parity=0, alpha=alpha)
        eva_a, _, _ = solve_sector(300, parity=1, alpha=alpha)
        dGp = sum(4 / eva_s[m] ** 2 - 4 / ev0_s[m] ** 2 for m in range(100))
        dGm = sum(4 / eva_a[m] ** 2 - 4 / ev0_a[m] ** 2 for m in range(100))
        tp = float(ANCHORS["lm_sum_rules_exact"][a_str]["Gp2"])
        tm = float(ANCHORS["lm_sum_rules_exact"][a_str]["Gm2"])
        assert abs(G2p0 + dGp - tp) < 5e-4
        assert abs(G2m0 + dGm - tm) < 5e-4


# --------------------------------------------------------------------------
# Eigenfunction-level validation (eigenfunctions.py)
# --------------------------------------------------------------------------
def test_parseval_completeness():
    """Completeness of the eigenvectors: sum_n (Int phi_n)^2 = 1 (symmetric
    tower) and sum_n (Int x phi_n)^2 = 1/3 (both sectors), up to the
    basis-truncation tail.  This checks the S-orthonormal basis is complete
    (the sums equal (pi^2/8)(S^-1)_00 and its x-analog for ANY complete
    S-orthonormal set), so it is a normalization check, not a check of the
    individual wavefunctions -- the Rayleigh test below is that."""
    from eigenfunctions import parseval_unit, parseval_x
    assert abs(parseval_unit(200) - 1.0) < 5e-5
    assert abs(parseval_x(200) - 1.0 / 3.0) < 5e-6


def test_decay_constant_asymptotics():
    """f_n (n + 3/4) -> 2/pi along the symmetric tower, monotonically from
    above, with the fitted constant within 1e-5 of 2/pi and a positive
    1/nu correction; f_0 to 13 digits."""
    from eigenfunctions import decay_constant_asymptotics
    d = decay_constant_asymptotics()
    assert abs(d["c0"] - 2.0 / PI) < 1e-5
    assert 0.05 < d["a"] < 0.08
    vals = [v for _, v in d["f_nu"]]
    assert all(vals[i + 1] < vals[i] for i in range(len(vals) - 1))
    assert vals[-1] > 2.0 / PI
    assert abs(d["f0"] - 0.9428101504585) < 2e-13


def test_rayleigh_weakform_recovers_eigenvalues():
    """Each eigenvector, pushed through the INDEPENDENT Gagliardo weak form,
    must reproduce its eigenvalue.  That validates the wavefunctions
    against the operator, not just against FLZ numbers.  On Gauss-Legendre
    nodes the endpoint sqrt(x) behavior caps this at ~1e-6; on the
    tanh-sinh rule the residual is at the 1e-12 level, i.e. the
    wavefunctions are right to the precision of the eigenvalues."""
    from eigenfunctions import rayleigh_weakform
    ev0, vecs0, deg0 = solve_sector(200, parity=0)
    ev1, vecs1, deg1 = solve_sector(200, parity=1)
    for n in range(3):
        assert abs(rayleigh_weakform(vecs0[:, n], deg0) - ev0[n]) < 1e-4
    assert abs(rayleigh_weakform(vecs1[:, 0], deg1) - ev1[0]) < 1e-4
    for n in range(6):
        assert abs(rayleigh_weakform(vecs0[:, n], deg0, n_quad=1201,
                                     rule="tanh-sinh") - ev0[n]) < 5e-12
    for n in range(3):
        assert abs(rayleigh_weakform(vecs1[:, n], deg1, n_quad=1201,
                                     rule="tanh-sinh") - ev1[n]) < 5e-12


def test_wavefunction_crosscheck_sine():
    """|<phi_0^cheb, phi_0^sine>| -> 1: the two solvers agree on the SHAPE of
    the ground state, at the sine solver's documented algebraic accuracy."""
    from eigenfunctions import reconstruct
    ev0, vecs0, deg0 = solve_sector(200, parity=0)
    _, svecs = sine_solver.solve(num_basis=80, n_quad=600, num_levels=1,
                                 return_vectors=True)
    xq, wq = np.polynomial.legendre.leggauss(400)
    xq = 0.5 * (xq + 1.0)
    wq = 0.5 * wq
    phi_cheb = reconstruct(vecs0[:, 0], deg0, xq)
    ks = np.arange(1, 81)
    phi_sine = np.sin(np.outer(xq, np.pi * ks)) @ svecs[:, 0]
    overlap = abs(np.sum(wq * phi_cheb * phi_sine))
    norms = np.sqrt(np.sum(wq * phi_cheb ** 2) * np.sum(wq * phi_sine ** 2))
    assert overlap / norms > 0.9999


def test_decay_constants_structure():
    """In the physical phi_n(0+) > 0 convention (canonicalize_signs), the
    decay constants are ALL POSITIVE, monotonically decreasing, with f_0
    dominant.  (Without a stated convention the signs are arbitrary
    LAPACK output.)"""
    from eigenfunctions import decay_constants
    f, _ = decay_constants(200, num_levels=8)
    assert f[0] > 0.9
    assert all(fi > 0 for fi in f)
    assert all(f[i + 1] < f[i] for i in range(7))


@pytest.mark.parametrize("alpha,tol", [(0.0, 1.5e-2), (0.5, 5e-3), (2.0, 5e-4)])
def test_independent_sine_crosscheck(alpha, tol):
    """The Chebyshev spectrum must agree with the *independent* sine-basis
    solver (different basis, different discretization, different singularity
    handling).  This is the real check on the quark-mass machinery, which is
    invisible at alpha = 0.  Tolerances reflect the sine basis' slow (algebraic)
    convergence, not the Chebyshev precision."""
    cheb = [L["two_lambda"]
            for L in spectrum(num_basis=300, num_levels=6, alpha=alpha)]
    sine = sine_solver.solve(num_basis=90, n_quad=700, alpha=alpha, num_levels=6)
    for n in range(6):
        assert abs(cheb[n] - sine[n]) < tol, (
            f"alpha={alpha} n={n}: cheb={cheb[n]:.6f} sine={sine[n]:.6f}")


def test_sine_converges_to_chebyshev():
    """Sanity: the sine-basis error must shrink as its basis grows (proving the
    small alpha=0 discrepancy is resolution, not a systematic disagreement)."""
    cheb0 = spectrum(300, 1, 0.0)[0]["two_lambda"]
    errs = [abs(sine_solver.solve(N, 4 * N + 200, 0.0, 1)[0] - cheb0)
            for N in (40, 80, 160)]
    assert errs[0] > errs[1] > errs[2]
    assert errs[-1] < 1e-3


def test_highprec_alpha_nonzero_path():
    """The mpmath Cholesky path (alpha != 0) agrees with the double-precision
    Chebyshev solver.  This is the only coverage of that branch (it has no
    other caller)."""
    from mpmath import mp
    from thooft_highprec import solve_sector_mp
    mp.dps = 30
    K, alpha = 16, 0.5
    mp_ev = [float(v) for v in solve_sector_mp(K, parity=0, alpha=alpha)]
    ev, _, _ = solve_sector(K, parity=0, alpha=alpha)
    assert max(abs(mp_ev[n] - ev[n]) for n in range(5)) < 1e-12


# --------------------------------------------------------------------------
# the committed artifacts, against the live solver
# --------------------------------------------------------------------------
def _reference_rows():
    path = Path(__file__).resolve().parent / "reference_spectrum.json"
    with open(path) as fh:
        return json.load(fh)


def _polished_rows():
    path = Path(__file__).resolve().parent / "reference_spectrum_polished.json"
    with open(path) as fh:
        return json.load(fh)["rows"]


def test_reference_spectrum_matches_live_solver():
    """reference_spectrum.json is the module's flagship output and is consumed
    cross-module by duality/dictionary.py, but nothing asserted its values
    against a live solve: its producer (generate_reference.py, mpmath) is a
    reproduce-slow target CI never runs, and the only other check covered
    row 0 and row 1 only.  Every level is now pinned."""
    rows = _reference_rows()
    assert len(rows) >= 24
    live = {d["n"]: d for d in spectrum(num_basis=300, num_levels=len(rows))}
    for r in rows:
        d = live[r["n"]]
        assert d["two_lambda"] == pytest.approx(float(r["value"]), rel=1e-10), \
            r["n"]
        # bookkeeping the dictionary layer relies on
        assert r["CP"] == (-1) ** (r["n"] + 1) == d["CP"]
        assert r["parity"] == ("sym" if r["n"] % 2 == 0 else "anti")
        assert d["parity"] == ("even" if r["n"] % 2 == 0 else "odd")
    vals = [float(r["value"]) for r in rows]
    assert all(b > a for a, b in zip(vals, vals[1:]))     # strictly ordered


def test_double_precision_csv_matches_reference():
    """spectrum_double_precision.csv is the 40-row published table; it must
    agree with the certified reference over their common range at the ~1e-11
    the README advertises, and stay interleaved beyond it."""
    import csv
    path = Path(__file__).resolve().parent / "spectrum_double_precision.csv"
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) >= 40
    ref = {r["n"]: float(r["value"]) for r in _reference_rows()}
    key = "two_lambda_n = M^2/(pi g^2)"
    for row in rows:
        n = int(row["n"])
        assert row["CP"] == ("-1" if n % 2 == 0 else "+1")
        assert row["parity_under_x->1-x"] == ("even" if n % 2 == 0 else "odd")
        if n in ref:
            assert float(row[key]) == pytest.approx(ref[n], abs=5e-12), n
    # against the 17-digit polished table, all 40 rows: the CSV is the
    # K = 800 reciprocal solve, good to 13 digits (measured 9e-14)
    pol = {r["n"]: float(r["value"]) for r in _polished_rows()}
    for row in rows:
        assert float(row[key]) == pytest.approx(pol[int(row["n"])], abs=4e-13)
    vals = [float(r[key]) for r in rows]
    assert all(b > a for a, b in zip(vals, vals[1:]))
    # the asymptote 2 lambda_n -> n + 3/4
    assert vals[-1] - (len(vals) - 1) == pytest.approx(0.75, abs=0.01)


# --------------------------------------------------------------------------
# FLZ spectral-determinant constants d+- (external_anchors.json, FLZ Eq. 5.1)
# --------------------------------------------------------------------------
def test_flz_d_constants_anchor(sectors):
    """CONSUMES the flz_d_constants_alpha0 anchor.  d+- are the constants in
    the large-(-lambda) asymptotics of FLZ's spectral determinants (their
    Eq. 4.16), computable from the eigenvalues alone through the convergent
    products of their Eq. (4.21) -- see flz_dconstants.py.  FLZ have no
    analytic expression for them individually; their Eq. (5.1) numbers are
    12-digit numerics, so the anchor's own rounding floor is +-5e-13.
    Measured here (K = 300, n_cut = 80, Eq. 1.6 tail + residual fit):
    d+ dev -1.5e-12, d- dev -7.6e-13 -- basis truncation of the levels
    near n = 90 at this K, not a floor: the ladder converges, to
    +6.7e-13 / +2.0e-13 by K = 600 (test_flz_d_constants_converged).  The
    exact-ratio check d-/d+ = sqrt(2)/pi (their Eq. 4.20) is
    convention-free end to end; measured 1.9e-13 with this 80-level head.
    Tolerances carry ~7x margin; a 1-2e-11 relative error in one low head
    level trips the d+ check, d- is gated at ~2.3e-11 relative (1e-11
    absolute on 0.43), and inside the 30..90 fit window the residual fit
    absorbs up to about half of a perturbation."""
    from flz_dconstants import d_constant_sector
    ev_sym, ev_anti = sectors
    dp = d_constant_sector(ev_sym, 0)
    dm = d_constant_sector(ev_anti, 1)
    ref = ANCHORS["flz_d_constants_alpha0"]
    dp_flz, dm_flz = float(ref["d+"]), float(ref["d-"])
    # the transcribed pair itself must satisfy the exact Eq. (4.20) ratio
    assert ref["ratio_exact"] == "sqrt(2)/pi"
    assert abs(dm_flz / dp_flz - np.sqrt(2.0) / PI) < 2e-13
    # normalization-free: our ratio against the exact sqrt(2)/pi
    assert abs(dm / dp - np.sqrt(2.0) / PI) < 5e-12
    # absolute: both constants at the anchor's own published precision
    assert abs(dp - dp_flz) < 1e-11, f"d+ dev {dp - dp_flz:+.2e}"
    assert abs(dm - dm_flz) < 1e-11, f"d- dev {dm - dm_flz:+.2e}"


def test_flz_d_constants_tail_ladder(sectors):
    """The Eq. (4.21) product truncated at n_cut = 80 is ~7e-7 off (that IS
    the size of the discarded tail, ~1/(pi^4 nu_cut^2)); appending the
    transcribed Eq. (1.6) tail recovers ~1e-10 (its own O(log^2 n/n^4)
    error); the measured-residual fit recovers the last digits.  Each rung
    must beat the previous by well over an order of magnitude."""
    from flz_dconstants import d_constant_sector
    ev_sym, _ = sectors
    dp_flz = float(ANCHORS["flz_d_constants_alpha0"]["d+"])
    bare = abs(d_constant_sector(ev_sym, 0, tail=False, resid_fit=False)
               - dp_flz)
    tail = abs(d_constant_sector(ev_sym, 0, resid_fit=False) - dp_flz)
    fit = abs(d_constant_sector(ev_sym, 0) - dp_flz)
    assert 1e-7 < bare < 1e-5          # measured 7.3e-7
    assert tail < 1e-9                 # measured 1.5e-10
    assert fit < 1e-11                 # measured 1.5e-12 at K = 300
    assert bare > 30 * tail > 30 * fit


def test_flz_d_constants_cross_solver_jacobi():
    """The same products from the Jacobi solver at alpha = 0: identical
    trial space, independent Gagliardo-quadrature assembly, so agreement is
    limited only by its documented quadrature floor.  Measured at K = 60:
    d+ dev -1.5e-8, d- dev -6.6e-9, ratio dev 2.1e-10."""
    from flz_dconstants import d_constants
    dp, dm = d_constants(num_basis=60, n_cut=50, fit_window=(30, 60),
                         solver="jacobi")
    ref = ANCHORS["flz_d_constants_alpha0"]
    assert abs(dp - float(ref["d+"])) < 2e-7
    assert abs(dm - float(ref["d-"])) < 2e-7
    assert abs(dm / dp - np.sqrt(2.0) / PI) < 5e-9


@pytest.mark.parametrize("parity", [0, 1])
def test_reciprocal_solver_preserves_eigenvector_normalization(parity):
    """The stable alpha=0 solve still returns S-orthonormal Ritz vectors."""
    from thooft_spectrum import build_sector
    _, D, S = build_sector(160, parity)
    values, vectors, _ = solve_sector(160, parity)
    V = vectors[:, :20]
    assert np.max(np.abs(V.T @ S @ V - np.eye(20))) < 2e-12
    residual = D @ V - (S @ V) * (PI ** 2 * values[:20])[None, :]
    assert np.linalg.norm(residual) / np.linalg.norm(D @ V) < 2e-12
    only, _, _ = solve_sector(160, parity, eigvals_only=True)
    assert np.allclose(only[:20], values[:20], rtol=2e-13, atol=1e-13)


# --------------------------------------------------------------------------
# the polished reference table (generate_reference.py --polish)
# --------------------------------------------------------------------------
def test_polished_table_certifies_the_older_tables():
    """reference_spectrum_polished.json: 40 levels, 17-21 digits, from
    exact-rational Rayleigh quotients of double-precision Ritz vectors on
    a K-ladder with the measured K^-6 truncation law extrapolated away.

    It certifies the two older anchors rather than replacing them: every
    row of the 24-level mpmath-eigensolver table agrees with it to at
    least that table's own claimed stable digits, and FLZ's printed
    Tables 1-2 agree to within one unit of their last printed digit
    (1.0e-12 worst, at n = 19; FLZ print 14 significant figures)."""
    from mpmath import mp, mpf
    pol = _polished_rows()
    assert [r["n"] for r in pol] == list(range(40))
    saved = mp.dps
    mp.dps = 30
    try:
        P = {r["n"]: mpf(r["value"]) for r in pol}
        for r in pol:
            assert r["stable_digits"] >= 17
            assert 5.99 < r["truncation_rate"] < 6.01          # K^-6
            assert r["richardson_correction"] < 0                # from above
            assert r["CP"] == (-1) ** (r["n"] + 1)
            assert r["parity"] == ("sym" if r["n"] % 2 == 0 else "anti")
        for r in _reference_rows():
            rel = abs(mpf(r["value"]) - P[r["n"]]) / P[r["n"]]
            assert rel < mpf(10) ** (-r["stable_digits"]), r["n"]
            assert mpf(r["value"]) > P[r["n"]]                   # variational
        worst = max(abs(mpf(repr(v)) - P[n]) for n, v in FLZ.items())
        assert worst < mpf("1.1e-12")
        assert abs(mpf(repr(FLZ_PADE_LAMBDA0)) - P[0]) < mpf("5e-16")
    finally:
        mp.dps = saved


def test_live_solver_against_polished_table():
    """The double-precision reciprocal solve at K = 800 reproduces all 40
    polished levels to 2.5e-13 (measured 9e-14): truncation is below
    rounding there, and there is no roundoff growth with K."""
    pol = {r["n"]: float(r["value"]) for r in _polished_rows()}
    for K, tol in ((300, 1.5e-12), (800, 2.5e-13)):
        live = spectrum(num_basis=K, num_levels=40)
        assert max(abs(d["two_lambda"] - pol[d["n"]]) for d in live) < tol, K


def test_rayleigh_polish_and_ladder_extrapolation():
    """The producer's two steps, re-executed.  The exact-rational Rayleigh
    quotient of the K = 300 double eigenvector gives the K = 300 Ritz value
    far past double precision: it sits ABOVE the stored infinite-basis
    value (variational) by the K^-6 truncation, 5.7e-18 at n = 0.  And the
    three-rung extrapolation recovers a planted limit and rate."""
    from mpmath import mp, mpf
    from thooft_highprec import extrapolate_ladder, rayleigh_polish
    pol = _polished_rows()
    saved = mp.dps
    mp.dps = 30
    try:
        for parity in (0, 1):
            _, vecs, _ = solve_sector(300, parity=parity)
            v = rayleigh_polish(vecs[:, 0], parity)
            gap = v - mpf(pol[parity]["value"])
            assert 0 < gap < mpf("1e-16"), (parity, gap)
        _, vecs, _ = solve_sector(300, parity=0)
        gap0 = rayleigh_polish(vecs[:, 0], 0) - mpf(pol[0]["value"])
        assert mpf("4e-18") < gap0 < mpf("8e-18")
        f = lambda K: mpf(2) + mpf(3) / mpf(K) ** 6
        v, r, corr = extrapolate_ladder(f(600), f(900), f(1350), 600, 900, 1350)
        assert abs(v - 2) < mpf("1e-25") and abs(r - 6) < mpf("1e-10")
        assert corr < 0
    finally:
        mp.dps = saved


def test_flz_d_constants_converged():
    """With the reciprocal solve the d+- ladder converges in K.  At
    K = 800 with a 160-level head the exact ratio sqrt(2)/pi is met to
    ~2e-15, and both constants sit inside the +-5e-13 rounding of FLZ's
    12 printed digits: d+ = 0.9631784563984(1), d- = 0.43358263983325(1)."""
    from flz_dconstants import d_constant_sector
    ev = {par: solve_sector(800, parity=par, eigvals_only=True)[0]
          for par in (0, 1)}
    ref = ANCHORS["flz_d_constants_alpha0"]
    vals = []
    for n_cut, win in ((160, (60, 200)), (200, (80, 260))):
        dp = d_constant_sector(ev[0], 0, n_cut=n_cut, fit_window=win)
        dm = d_constant_sector(ev[1], 1, n_cut=n_cut, fit_window=win)
        assert abs(dm / dp - np.sqrt(2.0) / PI) < 3e-14
        assert 0 < dp - float(ref["d+"]) < 5e-13
        assert 0 < dm - float(ref["d-"]) < 5e-13
        vals.append((dp, dm))
    assert abs(vals[0][0] - vals[1][0]) < 2e-13
    assert abs(vals[0][0] - 0.9631784563984) < 1.5e-13
    assert abs(vals[0][1] - 0.43358263983325) < 5e-14
    # the default 80-level head: converged in K, limited by its tail model
    dp300 = d_constant_sector(solve_sector(300, parity=0, eigvals_only=True)[0], 0)
    dp600 = d_constant_sector(solve_sector(600, parity=0, eigvals_only=True)[0], 0)
    dp800 = d_constant_sector(ev[0], 0)
    assert abs(dp600 - dp800) < 1e-13 < abs(dp300 - dp800)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
