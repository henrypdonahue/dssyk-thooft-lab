#!/usr/bin/env python3
"""Asserting regression tests for the massless 't Hooft spectrum.

Run with:  pytest -q        (or  python3 -m pytest -q)

Every check that validate.py *prints*, this file *asserts* with an explicit
tolerance, so the package fails loudly if a change breaks the physics.  The whole
suite is double precision and runs in a few seconds (the sine cross-check
quadrature dominates the cost).
"""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import zeta, polygamma

from thooft_spectrum import solve_sector, spectrum, _I_matrix
import sine_solver

PI = np.pi

# FLZ Tables 1 & 2 (arXiv:0905.2280), 2*lambda_n, 14 significant digits.
FLZ = {
    0: 0.73706174629269, 1: 1.7537313369175, 2: 2.7481609123706,
    3: 3.7510575817054, 4: 4.7492953810375, 5: 5.7504926236487,
    6: 6.7496294196488, 7: 7.7502843971925, 8: 8.7497715807892,
    9: 9.7501851352539, 10: 10.749845089160, 11: 11.750130142515,
    12: 12.749888008416, 13: 13.750096503972, 14: 14.749915244446,
    15: 15.750074428438, 16: 16.749933611057, 17: 17.750059159035,
    18: 18.749946584034, 19: 19.750048157169, 20: 20.749956088173,
    21: 21.750039967130, 22: 22.749963259761, 23: 23.750033705317,
    24: 24.749968804883, 25: 25.750028810060, 26: 26.749973181145,
    27: 27.750024910394, 28: 28.749976695731, 29: 29.750021753287,
}
FLZ_PADE_LAMBDA0 = 0.737061746292690  # FLZ Eq. (4.37)


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

    def sum_inv_sq(evs, offset):
        total, n = 0.0, offset
        for two_lam in evs:
            if n > n_cut:
                break
            total += 4.0 / two_lam ** 2
            n += 2
        total += polygamma(1, (n + 0.75) / 2.0)  # analytic tail, 2lam ~ n+3/4
        return total

    assert abs(sum_inv_sq(ev_sym, 0) - 7 * zeta(3)) < 1e-8
    assert abs(sum_inv_sq(ev_anti, 1) - 2.0) < 1e-8


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
    """Symmetric (even-degree) levels are FLZ even n with CP = -1; the
    interleaving alternates."""
    levels = spectrum(num_basis=100, num_levels=10)
    for L in levels:
        assert L["CP"] == (-1) ** (L["n"] + 1)
        assert L["parity"] == ("even" if L["n"] % 2 == 0 else "odd")


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
    must equal 4 I_{nm}, checked by an independent 1/x-form quadrature."""
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
import json as _json
import os as _os
from math import factorial as _factorial

with open(_os.path.join(_os.path.dirname(__file__), "external_anchors.json")) as _fh:
    ANCHORS = _json.load(_fh)


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
    closed forms -- four more independent whole-spectrum certificates."""
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


def test_lm_alpha_nonzero_eigenvalues():
    """The previously-open external gap: at alpha != 0 the spectrum must match
    the published Litvinov-Meshcheriakov tables (arXiv:2409.11324, 6
    significant digits in lambda -> ~1e-5 on 2*lambda for the low levels)."""
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
    absolutely (no alpha-dependent tail needed), so the exact 20-digit
    Litvinov-Meshcheriakov sum-rule values test the whole alpha != 0 spectrum
    against exact analytic results.  Accuracy here is limited by the basis'
    hardwired alpha = 0 endpoint exponent (algebraic convergence off the
    duality point) -- the tolerance below IS that measured limitation."""
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
    basis-truncation tail."""
    from eigenfunctions import parseval_unit, parseval_x
    assert abs(parseval_unit(200) - 1.0) < 5e-5
    assert abs(parseval_x(200) - 1.0 / 3.0) < 5e-6


def test_rayleigh_weakform_recovers_eigenvalues():
    """Each eigenvector, pushed through the INDEPENDENT Gagliardo weak form on
    Gauss-Legendre nodes, must reproduce its eigenvalue (validates the
    wavefunctions against the operator, not just against FLZ numbers)."""
    from eigenfunctions import rayleigh_weakform
    ev0, vecs0, deg0 = solve_sector(200, parity=0)
    ev1, vecs1, deg1 = solve_sector(200, parity=1)
    for n in range(3):
        assert abs(rayleigh_weakform(vecs0[:, n], deg0) - ev0[n]) < 1e-4
    assert abs(rayleigh_weakform(vecs1[:, 0], deg1) - ev1[0]) < 1e-4


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
    dominant.  (Without a stated convention the signs are arbitrary LAPACK
    output -- an earlier draft asserted an 'alternation' that was exactly
    that artifact.)"""
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


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
