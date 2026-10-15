"""Tests for full_channel_chord.py: the exact N = infinity law of the FULL
A_2 singlet channel, its closed forms, its gates, and its limits."""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from full_channel_chord import (channel, finite_n_full, inverse_p_limit,
                                mu0_closed, mu2_closed, mu4_closed, mu_full,
                                symbolic)

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def art():
    with open(HERE / "full_channel_chord.json") as fh:
        return json.load(fh)


@pytest.mark.parametrize("lam", [0.25, 1.0, 2.0, 3.0])
def test_matching_sum_equals_closed_forms(lam):
    """The enumerated matching sum reproduces the closed forms numerically."""
    q = math.exp(-lam)
    for k, f in ((0, mu0_closed), (2, mu2_closed), (4, mu4_closed)):
        assert abs(mu_full(k, lam, q) - f(lam, q)) < 1e-10 * abs(f(lam, q))


def test_symbolic_closed_forms_and_semiclassics():
    """sympy: the matching sum IS the closed form, and at small lambda
    w2 = 4 lam - 4 lam^2 + ..., kurt -> 7/2 (a Delta = 2 matter chord)."""
    s = symbolic()
    assert all(s["closed_form_checks"].values())
    assert s["w2_small_lambda"].startswith("4*lam - 4*lam**2")
    assert s["kurt_small_lambda"].startswith("7/2")


def test_lambda2_headline_numbers():
    """lambda_chord = 2: mu0 = 32 exactly, w2 = 1.316906, kurt = 3.852."""
    c = channel(2.0)
    assert abs(c["mu0"] - 32.0) < 1e-12
    assert abs(c["mu2"] - 42.141005) < 1e-6
    assert abs(c["w2"] - 1.316906) < 1e-6
    assert abs(c["kurt"] - 3.85217) < 1e-5


def test_finite_n_gate_recomputed_and_extrapolates(art):
    """A slice of the lambda_chord = 2 gate series is recomputed from the
    exact engine, and the stored series extrapolates in 1/p to the closed
    forms to better than 1e-4 relative (mu0, mu2) and 2e-5 (w2)."""
    rows = art["gates"]["lam2"]["rows"]
    for r in rows:
        if r["p"] in (4, 8):
            f = finite_n_full(r["N"], r["p"])
            assert str(f[0]) == r["mu0_exact"] and str(f[2]) == r["mu2_exact"]
    cf = art["gates"]["lam2"]["closed_form"]
    for key, tol in (("mu0", 1e-4), ("mu2", 1e-4), ("w2", 2e-5)):
        lim = inverse_p_limit(rows, key)
        assert abs(lim["central"] - cf[key]) < tol * abs(cf[key]), (key, lim, cf)
    assert abs(cf["w2"] - channel(2.0)["w2"]) < 1e-12


def test_second_line_gate(art):
    """lambda_chord = 1 (p <= 18): the same closed forms, to 1e-3."""
    rows = art["gates"]["lam1"]["rows"]
    cf = art["gates"]["lam1"]["closed_form"]
    for key in ("mu0", "mu2", "w2"):
        lim = inverse_p_limit(rows, key, p_min=6)
        assert abs(lim["central"] - cf[key]) < 1e-3 * abs(cf[key]), (key, lim, cf)


def test_w2_is_nonmonotone_in_lambda(art):
    """The N = infinity full width peaks near lambda_chord ~ 1 and falls on
    both sides -- the chord image of the fixed-p peak at lambda ~ 4/3."""
    scan = sorted(art["scan"], key=lambda r: r["lam"])
    w = [r["w2"] for r in scan]
    peak = max(range(len(w)), key=lambda i: w[i])
    assert 0.5 <= scan[peak]["lam"] <= 1.6
    assert w[0] < w[peak] and w[-1] < w[peak]


def test_operator_identity_dense():
    """A_2 = 4 sum_{c,c'} |c cap c'| J_c J_c' psi_c psi_c' per instance."""
    import sys
    from itertools import combinations
    sys.path.insert(0, str(HERE.parent / "syk-self-averaging"))
    from pauli_strings import majorana_hamiltonian_fast
    from syk import coupling_variance, majorana_operators
    for (N, p) in ((8, 4), (10, 4)):
        H = majorana_hamiltonian_fast(N, p, np.random.default_rng(3))
        rng = np.random.default_rng(3)
        sigma = np.sqrt(coupling_variance(N, p))
        combos = list(combinations(range(N), p))
        J = {c: rng.normal(0.0, sigma) for c in combos}
        psis = [op.toarray() for op in majorana_operators(N)]
        pref = (1j) ** (p // 2)
        G = {}
        for c in combos:
            M = np.eye(H.shape[0], dtype=complex)
            for a in c:
                M = M @ psis[a]
            G[c] = pref * M
        assert np.allclose(sum(J[c] * G[c] for c in combos), H)
        A2 = np.zeros_like(H)
        for P in psis:
            D = P.astype(complex)
            D = H @ D - D @ H
            D = H @ D - D @ H
            A2 += P @ D
        A2id = sum(4 * len(set(c) & set(d)) * J[c] * J[d] * (G[c] @ G[d])
                   for c in combos for d in combos)
        assert np.linalg.norm(A2id - A2) < 1e-12 * np.linalg.norm(A2)
