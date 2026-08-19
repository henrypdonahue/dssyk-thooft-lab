"""Exact finite-N pair moments (exact_pair_moments.py), validated four ways.

Pinned: the profile-sum evaluator equals the literal Isserlis enumeration
(all subset tuples, Pauli-string traces -- zero shared code) exactly; the
A_1 = -2pH flavor-sum identity ((N-1) mu2_pair + mu2_diag = 0) holds at
machine precision, gating the whole pipeline end to end; disorder-averaged
dense ED agrees within errors; t4bar converges to the chord t4 at matched
q = e^(-2p^2/N); and the committed campaign shows mu_2 falling to the
conserved-channel theorem's N = infinity zero.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from exact_pair_moments import (channel_pair_moments, literal_t4bar, t2bar,
                                t4bar)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "duality"))


@pytest.mark.parametrize("word", [(1, 1, 1, 1), (2, 0, 1, 1), (0, 2, 2, 0),
                                  (2, 2, 1, 1), (3, 1, 1, 1)])
def test_profile_sum_equals_literal_enumeration(word):
    """(N,p) = (6,2): the fast evaluator vs brute force over all subset
    tuples with Pauli-string traces (K = 4 and K = 6 words)."""
    assert t4bar(6, 2, *word) == pytest.approx(literal_t4bar(6, 2, *word),
                                               abs=1e-13)


def test_profile_sum_literal_diagonal_and_p4():
    """Diagonal marks at (6,2); pair marks at (8,4) (K = 4)."""
    for word in [(1, 1, 1, 1), (2, 1, 2, 1)]:
        lit = literal_t4bar(6, 2, *word, marks=('i', 'i', 'i', 'i'))
        assert t4bar(6, 2, *word, marks=('i', 'i', 'i', 'i')) == \
            pytest.approx(lit, abs=1e-13)
    for word in [(1, 1, 1, 1), (2, 0, 1, 1)]:
        assert t4bar(8, 4, *word) == pytest.approx(
            literal_t4bar(8, 4, *word), abs=1e-13)


def test_interleaved_and_invalid_marks():
    """Referee finding, fixed: interleaved distinct labels carry the
    psi-psi crossing sign (now included -- matches the literal trace);
    odd-multiplicity mark patterns have no sign law and raise."""
    for word in [(1, 1, 1, 1), (2, 0, 2, 0)]:
        lit = literal_t4bar(6, 2, *word, marks=('i', 'j', 'i', 'j'))
        assert t4bar(6, 2, *word, marks=('i', 'j', 'i', 'j')) == \
            pytest.approx(lit, abs=1e-13)
    with pytest.raises(ValueError):
        t4bar(6, 2, 1, 1, 1, 1, marks=('i', 'j', 'k', 'k'))


def test_t2bar_against_crossing_q():
    """t2bar(1,1) = E[(-1)^(|c&d|)]-free single-chord word: for the word
    H psi H psi the one H-chord crosses the psi-chord once, so
    t2bar(1,1) = E[(-1)^([i in c])] = 1 - 2p/N exactly."""
    for (N, p) in [(8, 4), (12, 4), (20, 6), (100, 4)]:
        assert t2bar(N, p, 1, 1) == pytest.approx(1.0 - 2.0 * p / N,
                                                  rel=1e-12)
        assert t2bar(N, p, 2, 0) == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize("N,p", [(8, 4), (10, 4), (12, 6), (30, 4), (64, 4)])
def test_a1_flavor_sum_identity(N, p):
    """A_1 = -2pH per instance forces (N-1) mu2_pair(n=1) + mu2_diag(n=1)
    = 0 exactly at every finite N -- an end-to-end gate spanning both
    matter structures and the imported assembly."""
    pair = channel_pair_moments(N, p, 1)
    diag = channel_pair_moments(N, p, 1, marks=('i', 'i', 'i', 'i'))
    scale = max(abs(diag["mu2"]), 1e-300)
    assert abs((N - 1) * pair["mu2"] + diag["mu2"]) < 1e-11 * scale


@pytest.mark.slow
def test_ed_cross_check():
    """Disorder-averaged dense ED at (10,4): every word within 4 sigma
    (independent route incl. the sigma_H normalization)."""
    from math import comb
    from syk import coupling_variance, majorana_hamiltonian, majorana_operators
    N, p, M = 10, 4, 200
    rng = np.random.default_rng(4242)
    psis = [m.toarray() for m in majorana_operators(N)]
    dim = psis[0].shape[0]
    sigH = math.sqrt(comb(N, p) * coupling_variance(N, p))
    words = [(1, 1, 1, 1), (2, 0, 1, 1), (2, 2, 1, 1)]
    acc = {w: [] for w in words}
    for _ in range(M):
        H = majorana_hamiltonian(N, p, rng).toarray()
        pw = {0: np.eye(dim), 1: H, 2: H @ H}
        i_, j_ = psis[0], psis[1]
        for (a, b, c, d) in words:
            val = np.trace(pw[a] @ i_ @ pw[b] @ i_ @ pw[c] @ j_
                           @ pw[d] @ j_).real / dim
            acc[(a, b, c, d)].append(val / sigH ** (a + b + c + d))
    for w in words:
        arr = np.array(acc[w])
        sem = arr.std(ddof=1) / np.sqrt(M)
        assert abs(arr.mean() - t4bar(N, p, *w)) < 4.0 * sem, w


def test_chord_limit():
    """t4bar -> chord t4 at matched q = e^(-2p^2/N), with the deviation
    falling ~ N^-2 (the Krawtchouk-microscopy rate)."""
    from chord_moments import MomentTable
    word = (1, 1, 1, 1)
    devs = []
    for N in (64, 256, 1024):
        q = math.exp(-32.0 / N)
        devs.append(abs(t4bar(N, 4, *word) - MomentTable(q, 0.25).t4(*word)))
    assert devs[0] > 8.0 * devs[1] > 8.0 * 8.0 * devs[2]


def test_committed_campaign():
    """The committed JSON: the A_1 gate on every fixed-p row; mu_2(n = 2)
    falls monotonically in magnitude toward the N = infinity zero along
    fixed p = 4 AND along the double-scaling line; the fixed-p tail is
    consistent with a power law (local slope steeper than 1/N^3)."""
    with open(HERE / "exact_pair_moments.json") as fh:
        data = json.load(fh)
    rows = data["fixed_p"]
    assert len(rows) >= 10
    for r in rows:
        assert r["a1_identity_residual"] < 1e-10
    mu2 = [abs(r["channels"]["n=2"]["mu2"]) for r in rows]
    ns = [r["N"] for r in rows]
    assert all(a > b for a, b in zip(mu2, mu2[1:]))
    # local log-slope over the tail: |mu2| ~ N^-s with s > 3
    s = (math.log(mu2[-2] / mu2[-1])
         / math.log(ns[-1] / ns[-2]))
    assert s > 3.0
    ds = [abs(r["n2"]["mu2"]) for r in data["double_scaling"]]
    assert all(a > b for a, b in zip(ds, ds[1:]))
