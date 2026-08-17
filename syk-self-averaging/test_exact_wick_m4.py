#!/usr/bin/env python3
"""Asserting tests for the m4 extension of exact_wick.py: exact disorder
statistics of the singlet fourth moment m4 = Tr(H^4)/dim.

Verification chain (each closed form against an independent route sharing
no assembly code):
  * exchange sign and support facts: literal Pauli-string traces;
  * Z, W combinatorics: literal set enumeration (frozensets, no binomials);
  * E[m4], Var(m4): the collected coupling polynomial + Gaussian moment
    formula (no Wick pairing taxonomy) at (6,2), (8,2), (6,4) fast and
    (8,4) slow-marked;
  * relVar(m4): fresh ED disorder ensembles at N = 12, 14;
  * micro-anchor (N=2, p=2): hand values E[m4] = 3 sigma^4, Var = 96 sigma^8.

Run with:  pytest -q -m 'not slow'   (the (8,4) enumeration is ~25 s)
"""

from itertools import combinations, product
from math import comb, exp, log, sqrt

import numpy as np
import pytest

from syk import coupling_variance
from exact_wick import (crossing_sum, crossing_q, mean_m4_exact,
                        quad_null_count, quad_null_alternating,
                        var_m4_exact, relvar_m4_exact, ln_relvar_m4,
                        relvar_m4_double_scaling_scan,
                        _mean_var_m4_enumerated, _trace_of_product)


# --------------------------------------------------------------------------
# the two structural facts of the derivation, witnessed by string traces
# --------------------------------------------------------------------------
@pytest.mark.parametrize("N,p", [(6, 2), (8, 4)])
def test_exchange_sign_witness(N, p):
    """Tr(Gamma_c Gamma_d Gamma_c Gamma_d)/dim = (-1)^{|c n d|} for EVERY pair
    of p-subsets (the four i^{p/2} prefactors multiply to +1 for even p)."""
    for c in combinations(range(N), p):
        for d in combinations(range(N), p):
            t = _trace_of_product(N, [c, d, c, d])
            assert abs(t - (-1) ** len(set(c) & set(d))) < 1e-12


def test_support_fact_witness():
    """A Gamma product has nonzero trace iff the symmetric difference of the
    subsets is empty -- and then the trace is exactly +-1."""
    N, p = 6, 4
    rng = np.random.default_rng(0)
    subsets = list(combinations(range(N), p))
    for _ in range(400):
        tup = [subsets[i] for i in rng.integers(0, len(subsets), 4)]
        sym = frozenset()
        for s in tup:
            sym = sym.symmetric_difference(s)
        t = _trace_of_product(N, tup)
        if sym:
            assert t == 0.0
        else:
            assert t == 1 or t == -1     # exactly +-1, never +-i


# --------------------------------------------------------------------------
# Z and W combinatorics vs literal set enumeration
# --------------------------------------------------------------------------
@pytest.mark.parametrize("N,p", [(6, 2), (6, 4), (8, 2)])
def test_quad_null_count_and_alternating(N, p):
    """Z = #4-tuples with empty symdiff and W = sum (-1)^{|x1 n x2|} over them,
    counted literally with frozensets (no binomials)."""
    subsets = [frozenset(c) for c in combinations(range(N), p)]
    z = w = 0
    for x1, x2 in product(subsets, repeat=2):
        s12 = x1.symmetric_difference(x2)
        inner = 0
        for x3, x4 in product(subsets, repeat=2):
            if x3.symmetric_difference(x4) == s12:
                inner += 1
        z += inner
        w += inner * (-1) ** len(x1 & x2)
    assert quad_null_count(N, p) == z
    assert quad_null_alternating(N, p) == w


# --------------------------------------------------------------------------
# closed forms vs the enumerated coupling polynomial
# --------------------------------------------------------------------------
def test_m4_micro_anchor_n2():
    """N=2, p=2: H = J Gamma with Gamma^2 = 1, so m4 = J^4 exactly:
    E[m4] = 3 sigma^4 (the degenerate c=d=e=f diagonal), Var = 96 sigma^8
    (E[J^8] = 105 sigma^8).  Hand values, no code shared with anything."""
    s2 = coupling_variance(2, 2)
    assert abs(mean_m4_exact(2, 2) - 3 * s2 ** 2) < 1e-12
    assert abs(var_m4_exact(2, 2) - 96 * s2 ** 4) < 1e-12
    mb, vb = _mean_var_m4_enumerated(2, 2)
    assert abs(mb - 3 * s2 ** 2) < 1e-12
    assert abs(vb - 96 * s2 ** 4) < 1e-12


@pytest.mark.parametrize("N,p", [(6, 2), (8, 2), (6, 4)])
def test_m4_closed_forms_match_enumeration(N, p):
    """E[m4] and Var(m4) closed forms vs the literal polynomial route."""
    mb, vb = _mean_var_m4_enumerated(N, p)
    assert abs(mean_m4_exact(N, p) / mb - 1.0) < 1e-9
    assert abs(var_m4_exact(N, p) / vb - 1.0) < 1e-9


def test_xor_filter_is_faithful():
    """The (8,4) speedup skips tuples with nonzero total support; those traces
    vanish identically (sampled witness) and the filtered enumeration agrees
    with the unfiltered one exactly."""
    N, p = 6, 4
    subsets = list(combinations(range(N), p))
    masks = [sum(1 << a for a in c) for c in subsets]
    rng = np.random.default_rng(1)
    seen_nonzero_xor = 0
    for _ in range(300):
        idx = rng.integers(0, len(subsets), 4)
        if masks[idx[0]] ^ masks[idx[1]] ^ masks[idx[2]] ^ masks[idx[3]]:
            seen_nonzero_xor += 1
            assert _trace_of_product(N, [subsets[i] for i in idx]) == 0.0
    assert seen_nonzero_xor > 100
    m1, v1 = _mean_var_m4_enumerated(N, p, filter_xor=False)
    m2, v2 = _mean_var_m4_enumerated(N, p, filter_xor=True)
    assert m1 == m2 and v1 == v2


@pytest.mark.slow
def test_m4_closed_forms_match_enumeration_8_4():
    """(8,4): the largest affordable enumeration (~25 s with the XOR filter)."""
    mb, vb = _mean_var_m4_enumerated(8, 4, filter_xor=True)
    assert abs(mean_m4_exact(8, 4) / mb - 1.0) < 1e-9
    assert abs(var_m4_exact(8, 4) / vb - 1.0) < 1e-9


# --------------------------------------------------------------------------
# ED disorder-ensemble cross-check (fresh seeds, no tuning)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("N,n_inst", [(12, 250), (14, 200)])
def test_relvar_m4_matches_ED_ensemble(N, n_inst):
    """Measured relVar(m4) over a fresh ED ensemble must agree with the exact
    closed form within 3 sampling standard errors."""
    from self_averaging import measure_moments
    res = measure_moments(N, 4, n_inst, seed=40 + N, ks=(2,))
    relvar, _, err = res[2]
    exact = relvar_m4_exact(N, 4)
    assert abs(relvar - exact) < 3 * err
    # and the ensemble-mean m4 agrees with E[m4] (SEM ~ sqrt(relvar/n))
    mean_m4 = res[2][1]
    sem = mean_m4 * sqrt(relvar / n_inst)
    assert abs(mean_m4 - mean_m4_exact(N, 4)) < 4 * sem


# --------------------------------------------------------------------------
# the crossing sum is the chord q-parameter (factor-2 landmine check)
# --------------------------------------------------------------------------
def test_crossing_q_chord_limit_fixed_p():
    """At fixed p, K/C -> exp(-2 p^2/N) with O(N^-2) error -- and is far from
    exp(-p^2/N), nailing lambda_chord = 2p^2/N (not p^2/N)."""
    p = 4
    errs = []
    for N in (100, 1000, 10000):
        s = crossing_q(N, p)
        errs.append(abs(s - exp(-2.0 * p * p / N)))
    assert errs[0] > errs[1] > errs[2]
    assert errs[1] < 3e-4 and errs[2] < 3e-6
    # scaling ~ N^-2: two decades in N shrink the error by ~1e4
    assert errs[0] / errs[2] > 1e3
    # the factor-2 landmine: exp(-p^2/N) is orders of magnitude farther away
    s = crossing_q(10000, p)
    assert abs(s - exp(-p * p / 10000)) > 100 * abs(s - exp(-2 * p * p / 10000))


def test_crossing_q_chord_limit_double_scaling():
    """Along p = sqrt(N) (lambda = 1): K/C -> exp(-2 lambda) = e^-2, with the
    O(1/p) error decreasing monotonically."""
    errs = []
    for (N, p) in [(100, 10), (400, 20), (1600, 40), (6400, 80)]:
        errs.append(abs(crossing_q(N, p) - exp(-2.0)))
    assert errs[0] > errs[1] > errs[2] > errs[3]
    assert errs[-1] / exp(-2.0) < 0.06


# --------------------------------------------------------------------------
# asymptotics: fixed-p rate and the double-scaling Eq. (3.27) law for m4
# --------------------------------------------------------------------------
def test_relvar_m4_leading_term_fixed_p():
    """relVar(m4) = 8/C(N,p) (1 + o(1)) at fixed p -- 4x relVar(m2), a pure
    N^-p power law at large N; the correction dies like ~1/C."""
    from self_averaging import relvar_m2_exact
    ratios = [relvar_m4_exact(N, 4) * comb(N, 4) / 8.0 for N in (12, 20, 40, 60)]
    assert all(r > 1.0 for r in ratios)          # Z + 2W correction is positive
    assert ratios[0] > ratios[1] > ratios[2] > ratios[3]
    assert abs(ratios[-1] - 1.0) < 2e-4          # measured 1.00013 at (60,4)
    assert abs(relvar_m4_exact(60, 4) / (4 * relvar_m2_exact(60, 4)) - 1) < 2e-4


def test_relvar_m4_double_scaling_superexponential():
    """Along p = sqrt(N): ln relVar(m4)/sqrt(N) decreases without bound (the
    same -(1/2) sqrt N ln N law as m2), and relVar(m4)/relVar(m2) -> 4."""
    rows = relvar_m4_double_scaling_scan(1.0, [16, 64, 100, 400, 1600])
    slopes = [r["ln_relvar_m4"] / sqrt(r["N"]) for r in rows]
    assert all(b < a for a, b in zip(slopes, slopes[1:]))
    assert slopes[-1] < -4.0
    for r in rows[1:]:
        ratio = exp(r["ln_relvar_m4"] - r["ln_relvar_m2"])
        assert abs(ratio - 4.0) < 1e-3
    # leading-coefficient check: ln relVar(m4) / (-(1/2) sqrt N ln N) -> 1
    r = rows[-1]
    assert abs(r["ln_relvar_m4"] / (-0.5 * sqrt(r["N"]) * log(r["N"])) - 1) < 0.25


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
