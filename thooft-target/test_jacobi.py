#!/usr/bin/env python3
"""Asserting tests for the matched-exponent Jacobi solver (jacobi_solver.py).

Every tolerance below is tied to a measurement recorded in
jacobi_convergence.json (regenerate with `python3 jacobi_solver.py`), with a
safety margin of roughly 10-100x over the measured residual -- not to hope.
The three correctness gates of the workstream:

  (a) alpha = 0: FLZ Tables 1-2 to ~1e-8 at moderate K   -> measured 5.3e-9
      at K = 60 (test_alpha0_flz_tables);
  (b) alpha != 0: LM anchors at least as well as the old solver, plus a
      demonstrated convergence-rate improvement                -> measured
      factor ~1.5e5 (alpha=0.5) and ~1.8e7 (alpha=-0.5) at K = 32
      (test_convergence_beats_old_solver);
  (c) higher-s LM sum rules s = 2..7 at all five LM alpha values including
      alpha = 1 (no published eigenvalue table there -- a pure prediction
      against their exact closed forms)   -> measured s2 <= 3e-8,
      s3 <= 2e-10, s4..7 <= ~2e-12 relative (test_lm_sum_rules_s2_to_s7).
"""

import json
import os

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_jacobi, roots_jacobi, zeta

from jacobi_solver import (
    endpoint_exponent, _jacobi_all, _jacobi_all_deriv, _log_norm_sq,
    _tanh_sinh_01, gagliardo_from_values, build_sector_jacobi,
    solve_sector_jacobi, spectrum_jacobi, sum_rule_G,
)

PI = np.pi
HERE = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(HERE, "external_anchors.json")) as _fh:
    ANCHORS = json.load(_fh)

# FLZ Tables 1 & 2 anchors: one provenance-checked copy (flz_reference.py)
from flz_reference import FLZ

LM_ALPHAS = ("-0.5", "-0.1", "0.1", "0.5")          # eigenvalue tables
SUM_RULE_ALPHAS = ("-0.5", "-0.1", "0.1", "0.5", "1")  # exact sum rules


@pytest.fixture(scope="module")
def lm_sectors():
    """Both CP sectors at K = 150 for every alpha with LM anchors."""
    out = {}
    for a_str in SUM_RULE_ALPHAS:
        alpha = float(a_str)
        ev_s, _, _ = solve_sector_jacobi(150, 0, alpha)
        ev_a, _, _ = solve_sector_jacobi(150, 1, alpha)
        out[a_str] = (ev_s, ev_a)
    return out


# ---------------------------------------------------------------------------
# endpoint exponent
# ---------------------------------------------------------------------------
def test_endpoint_exponent_alpha0_exact():
    assert endpoint_exponent(0.0) == 0.5


def test_endpoint_exponent_solves_transcendental():
    """beta solves pi b cot(pi b) = -alpha at machine precision, is monotone
    in alpha, and refuses alpha <= -1 (at/below the chiral point)."""
    prev = 0.0
    for alpha in (-0.9, -0.5, -0.1, 0.1, 0.5, 1.0, 5.0):
        b = endpoint_exponent(alpha)
        assert 0.0 < b < 1.0
        assert abs(PI * b / np.tan(PI * b) + alpha) < 5e-14 * max(1, abs(alpha))
        assert b > prev
        prev = b
    with pytest.raises(ValueError):
        endpoint_exponent(-1.0)


def test_endpoint_condition_source_identity():
    """The transcendental condition comes from the half-line finite part
    FP-Int_0^inf t^b/(t-1)^2 dt = -pi b cot(pi b)  ('t Hooft 1974; FLZ
    Sec. 1).  Verify that identity by an INDEPENDENT regularization:
    subtract 1 + b(t-1) on [0,2] (making the integrand regular at t = 1;
    the subtracted pieces have exact finite parts -2 and PV 0) and integrate
    the remainder with scipy.quad."""
    for beta in (0.371009648203552, 0.5, 0.645773676543405):
        i1, _ = quad(lambda t: (t ** beta - 1 - beta * (t - 1)) / (t - 1) ** 2,
                     0.0, 2.0, points=[1.0], limit=200)
        i2, _ = quad(lambda t: t ** beta / (t - 1) ** 2, 2.0, np.inf,
                     limit=200)
        fp = i1 + i2 - 2.0
        target = (0.0 if beta == 0.5
                  else -PI * beta * np.cos(PI * beta) / np.sin(PI * beta))
        assert abs(fp - target) < 1e-8


# ---------------------------------------------------------------------------
# basis machinery
# ---------------------------------------------------------------------------
def test_jacobi_recurrence_vs_scipy():
    """The in-house ultraspherical recurrence (values and derivatives) must
    match scipy.special.eval_jacobi -- independent implementation."""
    rng = np.random.default_rng(7)
    xi = rng.uniform(-1, 1, 9)
    for a in (0.742019296407, 1.0, 1.169214093401):   # 2*beta at LM alphas
        P = _jacobi_all(15, a, xi)
        dP = _jacobi_all_deriv(15, a, xi)
        for k in range(16):
            assert np.max(np.abs(P[k] - eval_jacobi(k, a, a, xi))) < 1e-11
            if k >= 1:
                ref = 0.5 * (k + 2 * a + 1) * eval_jacobi(k - 1, a + 1,
                                                          a + 1, xi)
                assert np.max(np.abs(dP[k] - ref)) < 1e-10


def test_norm_closed_form_gram_is_identity():
    """The Gamma-function normalization must make the Gauss-Jacobi Gram
    matrix exactly the identity (this is the 'mass-matrix conditioning'
    claim: weight [x(1-x)]^(2beta) pairs with the P^(2beta,2beta) family)."""
    beta = endpoint_exponent(0.5)
    for parity in (0, 1):
        deg = parity + 2 * np.arange(20)
        xq, wq = roots_jacobi(2 * int(deg[-1]) + 8, 2 * beta, 2 * beta)
        P = _jacobi_all(int(deg[-1]), 2 * beta, xq)[deg]
        P = P * np.exp(-0.5 * _log_norm_sq(deg, beta))[:, None]
        G = 2.0 ** (-4 * beta - 1.0) * ((P * wq) @ P.T)
        assert np.max(np.abs(G - np.eye(len(deg)))) < 1e-11


def test_potential_matrix_vs_independent_quadrature():
    """The solver's Gauss-Jacobi potential block <b_n, V b_m> (V = 1/x +
    1/(1-x) = 1/(x(1-x)), exact for polynomial integrands) against direct
    adaptive scipy.quad of [x(1-x)]^(2beta-1) p_n p_m -- an independent
    route through the endpoint-singular integrand."""
    alpha = 0.5
    beta = endpoint_exponent(alpha)
    degrees = 2 * np.arange(4)                    # even sector, K = 4
    Cln = -0.5 * _log_norm_sq(degrees, beta)
    # the solver's rule, replicated at the same size
    xq, wq = roots_jacobi(2 * int(degrees[-1]) + 8, 2 * beta - 1,
                          2 * beta - 1)
    Pq = _jacobi_all(int(degrees[-1]), 2 * beta, xq)[degrees]
    Pq = Pq * np.exp(Cln)[:, None]
    pot = 2.0 ** (1 - 4 * beta) * ((Pq * wq) @ Pq.T)
    for i in range(3):
        for j in range(3):
            di, dj = int(degrees[i]), int(degrees[j])
            ci, cj = np.exp(Cln[i]), np.exp(Cln[j])

            def integrand(x):
                xi = 2 * x - 1
                pi_ = eval_jacobi(di, 2 * beta, 2 * beta, xi)
                pj_ = eval_jacobi(dj, 2 * beta, 2 * beta, xi)
                return (x * (1 - x)) ** (2 * beta - 1) * ci * cj * pi_ * pj_

            ref, _ = quad(integrand, 0, 1, limit=400)
            assert abs(pot[i, j] - ref) < 1e-8


def test_gagliardo_exact_anchor_beta_half():
    """The quadrature engine against an EXACT closed form: at beta = 1/2,
    B[B_n, B_m] = (pi^2/2)(m+1) delta_nm - 2 I_nm for the Chebyshev-weighted
    functions B_n = sqrt(1-xi^2) U_n (from the exact eigenrelation of the
    finite-part operator + the Gagliardo identity + the odd-harmonic closed
    form _I_matrix of thooft_spectrum -- none of which share code with the
    tanh-sinh assembly).  Measured worst entry error 4.7e-13 at n_ts=501."""
    from thooft_spectrum import _I_matrix

    x, omx, w = _tanh_sinh_01(501)
    xi = x - omx
    s = np.sqrt(4 * x * omx)                    # sqrt(1-xi^2)
    md = 20
    U = np.zeros((md + 2, x.size))
    U[0] = 1.0
    U[1] = 2 * xi
    for m in range(2, md + 2):
        U[m] = 2 * xi * U[m - 1] - U[m - 2]
    Psi = s * U[:md + 1]
    om2 = 4 * x * omx                           # 1 - xi^2
    dPsi = np.zeros_like(Psi)
    for m in range(md + 1):
        dU = (np.zeros_like(xi) if m == 0
              else ((m + 2) * U[m - 1] - m * U[m + 1]) / (2 * om2))
        dPsi[m] = 2 * (-xi * U[m] + om2 * dU) / s
    B = gagliardo_from_values(Psi, dPsi, x, omx, w)
    deg = np.arange(md + 1)
    exact = (PI ** 2 / 2) * np.diag(deg + 1.0) - 2 * _I_matrix(deg)
    assert np.max(np.abs(B - exact)) < 1e-11


def test_assembled_matrix_is_hermitian():
    """The weak form guarantees a symmetric matrix; the builder asserts the
    discrete asymmetry is at rounding level and then symmetrizes exactly."""
    for alpha in (0.0, 0.5, -0.5):
        _, A = build_sector_jacobi(30, 0, alpha)
        assert np.array_equal(A, A.T)


# ---------------------------------------------------------------------------
# gate (a): alpha = 0
# ---------------------------------------------------------------------------
def test_alpha0_flz_tables():
    """Reproduce the FLZ 14-digit tables to < 1e-8 at moderate K = 60
    (measured worst deviation 5.3e-9; the old solver needs its exact
    closed-form matrix elements plus K = 300 for its 1e-10)."""
    ev_s, _, _ = solve_sector_jacobi(60, 0, 0.0)
    ev_a, _, _ = solve_sector_jacobi(60, 1, 0.0)
    worst = 0.0
    for n, ref in FLZ.items():
        mine = ev_s[n // 2] if n % 2 == 0 else ev_a[n // 2]
        worst = max(worst, abs(mine - ref))
    assert worst < 1e-8, f"worst FLZ deviation {worst:.2e}"


def test_alpha0_equals_chebyshev_trial_space():
    """At alpha = 0 (beta = 1/2 exactly) the Jacobi trial space IS the
    Chebyshev solver's trial space, so the two spectra must agree to the
    quadrature floor (measured 6.5e-13 at K = 40 over the lowest 20
    levels) -- any difference is pure Gagliardo-quadrature error."""
    from thooft_spectrum import solve_sector
    for parity in (0, 1):
        evj, _, _ = solve_sector_jacobi(40, parity, 0.0)
        evc, _, _ = solve_sector(40, parity=parity)
        assert np.max(np.abs(evj[:20] - evc[:20])) < 1e-11


# ---------------------------------------------------------------------------
# gate (b): alpha != 0 anchors and convergence-rate improvement
# ---------------------------------------------------------------------------
def test_lm_tables_alpha_nonzero(lm_sectors):
    """Litvinov-Meshcheriakov eigenvalue tables: over the lowest 10 levels
    the deviation must be < 3e-5 (the tolerance the existing solver's test
    passes at, i.e. 'at least as well'), and over ALL 30 tabulated levels
    < 2e-4 (the LM tables' own 6-significant-digit rounding floor -- the
    solver is far more accurate than the anchor here)."""
    for a_str in LM_ALPHAS:
        ev_s, ev_a = lm_sectors[a_str]
        tab = ANCHORS["lm_two_lambda_tables"][a_str]
        worst_low, worst_all = 0.0, 0.0
        for m in range(15):
            de = abs(ev_s[m] - tab["even"][m])
            do = abs(ev_a[m] - tab["odd"][m])
            worst_all = max(worst_all, de, do)
            if 2 * m < 10:
                worst_low = max(worst_low, de)
            if 2 * m + 1 < 10:
                worst_low = max(worst_low, do)
        assert worst_low < 3e-5, f"alpha={a_str}: low levels {worst_low:.2e}"
        assert worst_all < 2e-4, f"alpha={a_str}: all levels {worst_all:.2e}"


def test_convergence_beats_old_solver():
    """THE rate gate: at K = 32 the Jacobi ground-state error must beat the
    old solver's by > 1e3 (measured: 1.5e5 at alpha = 0.5, 1.8e7 at
    alpha = -0.5) and sit below 1e-10 absolute (measured ~4e-12 with
    reference self-consistency ~4e-12, itself certified externally by the
    LM tables and 20-digit sum rules).  Also: the Jacobi error must fall
    monotonically over K = 8 -> 16 -> 32 (spectral-like), while the old
    solver's decreases only algebraically (~K^(-4 beta), measured rates
    2.34 and 1.48 vs 4 beta = 2.34, 1.48 -- see jacobi_convergence.json)."""
    from thooft_spectrum import solve_sector
    for alpha in (0.5, -0.5):
        ref, _, _ = solve_sector_jacobi(200, 0, alpha)
        jac_errs = []
        for K in (8, 16, 32):
            jac, _, _ = solve_sector_jacobi(K, 0, alpha)
            jac_errs.append(abs(jac[0] - ref[0]))
        old, _, _ = solve_sector(32, parity=0, alpha=alpha)
        old_err = abs(old[0] - ref[0])
        assert jac_errs[0] > jac_errs[1] > jac_errs[2]
        assert jac_errs[2] < 1e-10
        assert old_err > 1e3 * jac_errs[2], (
            f"alpha={alpha}: old {old_err:.2e} vs jacobi {jac_errs[2]:.2e}")


def test_convergence_json_matches_live_solver():
    """The committed jacobi_convergence.json must describe THIS code: replay
    the K = 16 row at alpha = 0.5 against the stored reference values."""
    from thooft_spectrum import solve_sector
    with open(os.path.join(HERE, "jacobi_convergence.json")) as fh:
        blk = json.load(fh)["old_vs_jacobi"]["0.5"]
    ref_n0 = blk["reference_two_lambda_n0"]
    row = next(r for r in blk["rows"] if r["K"] == 16)
    jac, _, _ = solve_sector_jacobi(16, 0, 0.5)
    old, _, _ = solve_sector(16, parity=0, alpha=0.5)
    assert abs(abs(jac[0] - ref_n0) - row["jac_err_n0"]) < 1e-10
    assert abs(abs(old[0] - ref_n0) - row["old_err_n0"]) < 1e-10
    assert abs(blk["beta"] - endpoint_exponent(0.5)) < 1e-14


def test_spectrum_interleaving_and_cp():
    """spectrum_jacobi keeps the main solver's labeling contract at
    alpha != 0: strict CP interleaving, CP = (-1)^(n+1)."""
    levels = spectrum_jacobi(num_basis=40, num_levels=12, alpha=0.5)
    for L in levels:
        assert L["CP"] == (-1) ** (L["n"] + 1)
        assert L["parity"] == ("even" if L["n"] % 2 == 0 else "odd")
    vals = [L["two_lambda"] for L in levels]
    assert all(v1 < v2 for v1, v2 in zip(vals, vals[1:]))


# ---------------------------------------------------------------------------
# gate (c): the higher-s LM sum rules
# ---------------------------------------------------------------------------
def test_sum_rule_machinery_alpha0():
    """The drift-corrected tail machinery, validated where the answer is
    exactly known (FLZ Eq. 1.8 and the s = 3 closed forms): measured
    residuals ~6e-11 at K = 150, n_cut = 120."""
    ev_s, _, _ = solve_sector_jacobi(150, 0, 0.0)
    ev_a, _, _ = solve_sector_jacobi(150, 1, 0.0)
    targets = {
        (2, 0): 7 * zeta(3.0), (2, 1): 2.0,
        (3, 0): -4 * PI ** 2 / 3 + 28 * zeta(3.0),
        (3, 1): -8.0 / 3.0 + 4 * PI ** 2 / 9,
    }
    for (s, parity), tgt in targets.items():
        got = sum_rule_G(ev_s if parity == 0 else ev_a, parity, s)
        assert abs(got - tgt) < 1e-9, f"s={s} parity={parity}"


def test_lm_sum_rules_s2_to_s7(lm_sectors):
    """All LM exact alpha-dependent sum rules G^(s)+- for s = 2..7 at all
    five alpha values -- including alpha = 1, where LM publish no eigenvalue
    table, so this is a pure prediction against their closed forms.
    Tolerances are the MEASURED residuals (jacobi_convergence.json,
    lm_sum_rule_residuals: s2 <= 3e-8, s3 <= 2e-10, s4..7 <= ~2e-12
    relative) with ~10-50x margin.  The old solver's documented limitation
    was 1e-5..2e-4 at s = 2 only; this beats it by > 3 orders of magnitude
    and extends to s = 7."""
    tol_abs = {2: 2e-7, 3: 2e-9}
    tol_rel = 1e-10  # s = 4..7 (large G^(s)+ values: relative is the honest metric)
    for a_str in SUM_RULE_ALPHAS:
        ev_s, ev_a = lm_sectors[a_str]
        exact = ANCHORS["lm_sum_rules_exact"][a_str]
        for s in range(2, 8):
            gp = sum_rule_G(ev_s, 0, s)
            gm = sum_rule_G(ev_a, 1, s)
            tp, tm = float(exact[f"Gp{s}"]), float(exact[f"Gm{s}"])
            if s in tol_abs:
                assert abs(gp - tp) < tol_abs[s], f"alpha={a_str} Gp{s}"
                assert abs(gm - tm) < tol_abs[s], f"alpha={a_str} Gm{s}"
            else:
                assert abs(gp - tp) < tol_rel * abs(tp) + 1e-11, \
                    f"alpha={a_str} Gp{s}"
                assert abs(gm - tm) < tol_rel * abs(tm) + 1e-11, \
                    f"alpha={a_str} Gm{s}"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
