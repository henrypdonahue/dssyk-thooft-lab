#!/usr/bin/env python3
"""Asserting tests for exact_wick_dirac.py.

Every closed form is checked against a literal enumeration route that shares
no assembly code with it (T-matrix traces via dirac.apply_monomial chains vs
pure binomial formulas; Isserlis pair sums vs the moments_pipeline ED means),
plus instance-level identities against the dense ED-side algebra.

Run with:  pytest -q -m 'not slow'   (ED ensembles / Nc=8 / p=8 are slow)
"""

import json
from math import comb, factorial, log
from pathlib import Path

import numpy as np
import pytest

from dirac import (annihilators, coupling_variance_dirac, dirac_hamiltonian,
                   draw_couplings, mn_tower, particle_hole)
from exact_wick_dirac import (chain_trace, instance_b, ln_rms_offdiag_csym,
                              m2_from_pairs, m2_mean_exact, mean_diag_csym,
                              mu0_exact_csym, mu2_exact_csym,
                              offdiag_census, ordered_keys,
                              rms_offdiag_csym, rms_offdiag_large_nc,
                              stats_enumerated, t_matrix, unordered_pairs,
                              var_offdiag_csym)

HERE = Path(__file__).resolve().parent


# --------------------------------------------------------------------------
# Part 1: Var(B_ij) and E[B_ii] closed forms vs literal enumeration
# --------------------------------------------------------------------------
@pytest.mark.parametrize("Nc,p", [(5, 4), (6, 4)])
def test_offdiag_closed_form_vs_enumeration(Nc, p):
    """E[B_ij] = 0 and Var(B_ij) = 4 sigma^4 C(Nc-2,k-1) C(Nc-k-1,k)/4^(p+1),
    against the full quadratic-form enumeration over ordered coupling keys."""
    mean, var = stats_enumerated(Nc, p, 0, 1)
    assert mean == 0.0
    assert var == pytest.approx(var_offdiag_csym(Nc, p), rel=1e-12)


def test_offdiag_census_p4():
    """The derived admissible-pair census: 2 L nonzero T entries (L Type I +
    L Type II), every value -2^-(p+1), and T symmetric (the measured
    T_{alpha beta} = T_{beta alpha} fact the Wick pairing sum relies on)."""
    Nc, p = 6, 4
    k = p // 2
    L = comb(Nc - 2, k - 1) * comb(Nc - k - 1, k)
    n_nonzero, values, symmetric = offdiag_census(Nc, p, 0, 1)
    assert n_nonzero == 2 * L
    assert values == [pytest.approx(-2.0 ** -(p + 1))]
    assert symmetric


@pytest.mark.parametrize("Nc,p,i", [(5, 4, 0), (6, 4, 0), (6, 4, 2)])
def test_mean_diag_closed_form_vs_enumeration(Nc, p, i):
    """E[B_ii] = sigma^2 C(Nc-1,k) C(Nc-k-1,k) 2^-(p+1), any i."""
    mean, _ = stats_enumerated(Nc, p, i, i)
    assert mean == pytest.approx(mean_diag_csym(Nc, p), rel=1e-12)


def test_mean_diag_vs_m2_relation():
    """E[B_ii] = E[m2] (Nc - p)/(2 Nc), with E[m2] recomputed numerically
    from the pair operators (independently reproducing the qed_campaign.py
    c_symmetric closed form)."""
    for (Nc, p) in [(5, 4), (6, 4), (8, 4)]:
        m2 = m2_mean_exact(Nc, p)
        k = p // 2
        m2_closed = (coupling_variance_dirac(Nc, p) * comb(Nc, k)
                     * comb(Nc - k, k) * 4.0 ** -k)
        assert m2 == pytest.approx(m2_closed, rel=1e-12)
        assert mean_diag_csym(Nc, p) == pytest.approx(
            m2 * (Nc - p) / (2.0 * Nc), rel=1e-12)


def test_mean_vanishing_artifact_at_nc_equals_p():
    """At Nc = p every coupling touches every mode: E[B_ii] = 0 exactly AND
    B_ij vanishes identically (Var = 0) -- the U(N) analog of the Majorana
    p = N/2 ratio artifact, derived rather than flagged."""
    Nc, p = 4, 4
    assert mean_diag_csym(Nc, p) == 0.0
    assert var_offdiag_csym(Nc, p) == 0.0
    mean_d, _ = stats_enumerated(Nc, p, 0, 0)
    mean_o, var_o = stats_enumerated(Nc, p, 0, 1)
    assert mean_d == 0.0 and mean_o == 0.0 and var_o == 0.0


def test_instance_quadratic_form_identity():
    """For one c_symmetric instance, B_ij from the enumerated T-matrix
    quadratic form must equal the dense-algebra Tr(c_i H c^dag_j H)/dim."""
    Nc, p = 6, 4
    k = p // 2
    rng = np.random.default_rng(21)
    coup = draw_couplings(Nc, p, rng, "c_symmetric")
    keys = ordered_keys(Nc, k)
    J = np.array([coup[key] for key in keys])
    for (i, j) in [(0, 1), (2, 4), (0, 0)]:
        T = t_matrix(Nc, k, i, j, keys)
        b_quad = float(J @ T @ J)
        assert b_quad == pytest.approx(instance_b(Nc, p, coup, i, j),
                                       rel=1e-10, abs=1e-12)


def test_large_nc_rate_and_stable_log():
    """Fixed-p law: Var -> k (k!)^2 / Nc^(p-1); and the lgamma form of
    ln rms agrees with the direct closed form."""
    p, k = 4, 2
    for Nc in [10, 20, 40]:
        assert ln_rms_offdiag_csym(Nc, p) == pytest.approx(
            log(rms_offdiag_csym(Nc, p)), rel=1e-12)
    ratio = (var_offdiag_csym(10 ** 4, p)
             / (k * factorial(k) ** 2 / (10 ** 4) ** (p - 1)))
    assert ratio == pytest.approx(1.0, rel=5e-3)     # 1 + O(p^2/Nc)
    assert rms_offdiag_large_nc(10 ** 4, p) == pytest.approx(
        rms_offdiag_csym(10 ** 4, p), rel=3e-3)


@pytest.mark.slow
def test_offdiag_census_p8():
    """p = 8 (k = 4, the next c_symmetric-admissible p) at the minimal
    Nc = 9: census count 2 L, uniform value -2^-9, T symmetric, and the
    closed form against the full enumeration."""
    Nc, p = 9, 8
    k = p // 2
    L = comb(Nc - 2, k - 1) * comb(Nc - k - 1, k)
    n_nonzero, values, symmetric = offdiag_census(Nc, p, 0, 1)
    assert n_nonzero == 2 * L
    assert values == [pytest.approx(-2.0 ** -(p + 1))]
    assert symmetric
    mean, var = stats_enumerated(Nc, p, 0, 1)
    assert mean == 0.0
    assert var == pytest.approx(var_offdiag_csym(Nc, p), rel=1e-12)


@pytest.mark.slow
def test_var_offdiag_ed_ensemble():
    """Sample Var(B_01) and mean B_00 over an ED disorder ensemble at
    Nc = 6 against the closed forms (moment-based error bars, ddof=1)."""
    Nc, p, n_inst = 6, 4, 400
    rng = np.random.default_rng(2026)
    b01, b00 = [], []
    ann = annihilators(Nc)
    dim = 1 << Nc
    c0, c1 = ann[0], ann[1]
    for _ in range(n_inst):
        coup = draw_couplings(Nc, p, rng, "c_symmetric")
        H = dirac_hamiltonian(Nc, p, couplings=coup)
        b01.append(float(np.real(np.trace(c0 @ H @ c1.conj().T @ H))) / dim)
        b00.append(float(np.real(np.trace(c0 @ H @ c0.conj().T @ H))) / dim)
    b01 = np.array(b01)
    b00 = np.array(b00)
    # mean diagonal: straight 3-sigma SEM check
    sem = b00.std(ddof=1) / np.sqrt(n_inst)
    assert abs(b00.mean() - mean_diag_csym(Nc, p)) < 3 * sem
    # off-diagonal variance: delta-method error on the sample variance
    v = b01.var(ddof=1)
    m4 = np.mean((b01 - b01.mean()) ** 4)
    sem_v = np.sqrt(max(m4 - v ** 2, 0.0) / n_inst)
    assert abs(v - var_offdiag_csym(Nc, p)) < 4 * sem_v
    # off-diagonal mean is exactly zero in expectation
    assert abs(b01.mean()) < 4 * b01.std(ddof=1) / np.sqrt(n_inst)


# --------------------------------------------------------------------------
# Part 2: M_2 pair decomposition and exact mu0 / mu2
# --------------------------------------------------------------------------
def test_m2_pair_decomposition_identity():
    """M_2 = sum_uv x_u x_v G_uv must reproduce dirac.mn_tower's M_2 for a
    c_symmetric instance (independent assembly: per-pair sparse commutators
    vs derivative chains on the dense H)."""
    Nc, p = 6, 4
    rng = np.random.default_rng(3)
    coup = draw_couplings(Nc, p, rng, "c_symmetric")
    H = dirac_hamiltonian(Nc, p, couplings=coup)
    M2_ref = mn_tower(H, annihilators(Nc), 2)[2]
    M2_pairs = m2_from_pairs(Nc, p, coup)
    assert np.abs(M2_ref - M2_pairs).max() < 1e-12
    assert np.abs(M2_ref.imag).max() < 1e-12       # M_2 is real symmetric


def test_mu0_channel_conventions_vs_pipeline_instance():
    """Per-instance mu0 of the CP-resolved traceless M_2, assembled from the
    pair decomposition + particle_hole projection, must match
    moments_pipeline.instance_moments (n = 2: right = CP-odd, wrong =
    CP-even) for the same seed."""
    from exact_wick_dirac import instance_mu_ed
    Nc, p, seed = 6, 4, 100
    dim = 1 << Nc
    rng = np.random.default_rng(seed)          # instance_moments' convention
    coup = draw_couplings(Nc, p, rng, "c_symmetric")
    M2 = m2_from_pairs(Nc, p, coup)
    M2b = M2 - np.trace(M2) / dim * np.eye(dim)
    U = particle_hole(Nc).real
    M2c = U @ M2b @ U.T
    mu0_odd = float(np.sum((0.5 * (M2b - M2c)) ** 2)) / dim
    mu0_even = float(np.sum((0.5 * (M2b + M2c)) ** 2)) / dim
    ed = instance_mu_ed(Nc, p, seed)
    assert mu0_odd == pytest.approx(ed["odd"][0], rel=1e-10)
    assert mu0_even == pytest.approx(ed["even"][0], rel=1e-10)


def _moments_rows(Nc):
    rows = json.loads((HERE / "moments.json").read_text())
    return {r["channel"]: r for r in rows if r["Nc"] == Nc and r["n"] == 2}


def test_mu0_exact_vs_moments_json_nc6():
    """Exact Isserlis E[mu0] vs the stored ED disorder means (40 instances):
    the ED-free cross-check of the moments pipeline at Nc = 6."""
    exact = mu0_exact_csym(6, 4)
    ed = _moments_rows(6)
    for ch, tag in (("odd", "right"), ("even", "wrong")):
        z = (exact[ch] - ed[tag]["mu0"]) / ed[tag]["mu0_err"]
        assert abs(z) < 3.0, f"channel {ch}: z = {z:.2f}"


@pytest.mark.slow
def test_mu0_exact_vs_moments_json_nc8():
    """Same at Nc = 8 (210 pairs, ~30 s)."""
    exact = mu0_exact_csym(8, 4)
    ed = _moments_rows(8)
    for ch, tag in (("odd", "right"), ("even", "wrong")):
        z = (exact[ch] - ed[tag]["mu0"]) / ed[tag]["mu0_err"]
        assert abs(z) < 3.0, f"channel {ch}: z = {z:.2f}"


@pytest.mark.slow
def test_mu0_mu2_exact_vs_fresh_ed_ensemble_nc6():
    """Exact E[mu0] AND E[mu2] (sextic Isserlis) vs a fresh 300-instance ED
    ensemble at Nc = 6 -- the direct calibrated check of the mu2 pairing
    algebra (moments.json stores no plain mu2 mean)."""
    from exact_wick_dirac import instance_mu_ed
    Nc, p, n_inst = 6, 4, 300
    vals = {("odd", 0): [], ("odd", 2): [], ("even", 0): [], ("even", 2): []}
    for j in range(n_inst):
        inst = instance_mu_ed(Nc, p, 7000 + j)
        for ch in ("odd", "even"):
            vals[(ch, 0)].append(inst[ch][0])
            vals[(ch, 2)].append(inst[ch][1])
    mu0_ex = mu0_exact_csym(Nc, p)
    mu2_ex = mu2_exact_csym(Nc, p, progress=False)
    for ch in ("odd", "even"):
        for k_mom, ex in ((0, mu0_ex[ch]), (2, mu2_ex[ch])):
            arr = np.array(vals[(ch, k_mom)])
            sem = arr.std(ddof=1) / np.sqrt(n_inst)
            z = (arr.mean() - ex) / sem
            assert abs(z) < 3.5, f"mu{k_mom}[{ch}]: z = {z:.2f}"
