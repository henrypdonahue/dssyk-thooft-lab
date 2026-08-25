"""Exact finite-N pair moments (exact_pair_moments.py), validated five ways.

Pinned: the profile-sum evaluator equals the literal Isserlis enumeration
exactly (all subset tuples, Pauli-string traces -- zero shared code).  The
A_1 = -2pH flavor-sum identity ((N-1) mu2_pair + mu2_diag = 0) holds at
machine precision, and gates the whole pipeline end to end.
Disorder-averaged dense ED agrees within errors.  t4bar converges to the
chord t4 at matched q = e^(-2p^2/N).  The committed campaign shows mu_2
falling to the conserved-channel theorem's N = infinity zero.
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
    psi-psi crossing sign (now included -- matches the literal trace).
    Odd-multiplicity mark patterns have no sign law and raise."""
    for word in [(1, 1, 1, 1), (2, 0, 2, 0)]:
        lit = literal_t4bar(6, 2, *word, marks=('i', 'j', 'i', 'j'))
        assert t4bar(6, 2, *word, marks=('i', 'j', 'i', 'j')) == \
            pytest.approx(lit, abs=1e-13)
    with pytest.raises(ValueError):
        t4bar(6, 2, 1, 1, 1, 1, marks=('i', 'j', 'k', 'k'))


def test_t2bar_against_crossing_q():
    """t2bar(1,1) = E[(-1)^(|c&d|)]-free single-chord word.  For the word
    H psi H psi the one H-chord crosses the psi-chord once, so
    t2bar(1,1) = E[(-1)^([i in c])] = 1 - 2p/N exactly."""
    for (N, p) in [(8, 4), (12, 4), (20, 6), (100, 4)]:
        assert t2bar(N, p, 1, 1) == pytest.approx(1.0 - 2.0 * p / N,
                                                  rel=1e-12)
        assert t2bar(N, p, 2, 0) == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize("N,p", [(8, 4), (10, 4), (12, 6), (30, 4), (64, 4)])
def test_a1_flavor_sum_identity(N, p):
    """A_1 = -2pH per instance forces (N-1) mu2_pair(n=1) + mu2_diag(n=1)
    = 0 exactly at every finite N.  This is an end-to-end gate spanning
    both matter structures and the imported assembly."""
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


def test_dp_tables_equal_profile_sum():
    """The CRT residue-table DP equals the closed Venn-profile sum on
    every m <= 3 signature (all S, up to two mark vectors), with exact
    Fractions.  The DP is the only path for m >= 4, so this equality on
    the shared domain is its correctness anchor."""
    from itertools import combinations, product
    from exact_pair_moments import (_expectation_core,
                                    _expectation_core_profile)
    cases = 0
    for m in (1, 2, 3):
        pairs = list(combinations(range(m), 2))
        Ss = [frozenset(c) for r in range(len(pairs) + 1)
              for c in combinations(pairs, r)]
        for S in Ss:
            for nbits in range(3):
                for bits in product(*[list(product((0, 1), repeat=m))]
                                    * nbits):
                    bits = tuple(b for b in bits if any(b))
                    assert _expectation_core(10, 4, m, S, bits) == \
                        _expectation_core_profile(10, 4, m, S, bits)
                    cases += 1
    assert cases > 500


def test_k8_word_literal():
    """A K = 8 word (four subsets -- the DP-only regime) against the
    literal Isserlis enumeration at (4,2)."""
    for word in [(2, 2, 2, 2), (3, 1, 2, 2)]:
        assert t4bar(4, 2, *word) == pytest.approx(
            literal_t4bar(4, 2, *word), abs=1e-12)


@pytest.mark.parametrize("N,p", [(8, 4), (16, 4)])
def test_a1_flavor_sum_identity_mu4(N, p):
    """The mu_4 version of the A_1 gate: (N-1) mu4_pair + mu4_diag = 0
    exactly (K = 6 words, so this spans the m <= 3 DP signatures).
    test_k8_word_literal covers m = 4."""
    pair = channel_pair_moments(N, p, 1, ks=(0, 2, 4))
    diag = channel_pair_moments(N, p, 1, marks=('i', 'i', 'i', 'i'),
                                ks=(0, 2, 4))
    scale = max(abs(diag["mu4"]), 1e-300)
    assert abs((N - 1) * pair["mu4"] + diag["mu4"]) < 1e-10 * scale


def test_committed_campaign():
    """The committed JSON: the A_1 gate holds on every fixed-p row.
    mu_2(n = 2) falls monotonically in magnitude toward the N = infinity
    zero along fixed p = 4 AND along the double-scaling line.  The fixed-p
    tail is consistent with a power law (local slope steeper than
    1/N^3)."""
    with open(HERE / "exact_pair_moments.json") as fh:
        data = json.load(fh)
    rows = data["fixed_p"]
    assert len(rows) >= 10
    for r in rows:
        assert r["a1_identity_residual"] < 1e-10
    mu2 = [abs(r["channels"]["n=2"]["mu2"]) for r in rows]
    ns = [r["N"] for r in rows]
    assert all(a > b for a, b in zip(mu2, mu2[1:]))
    for key in ("n=2", "n=3"):
        m4 = [abs(r["channels"][key]["mu4"]) for r in rows[4:]]
        assert all(a > b for a, b in zip(m4, m4[1:])), key
        assert all(r["a1_identity_residual_mu4"] < 1e-9 for r in rows)
    # local log-slope over the tail: |mu2| ~ N^-s with s > 3
    s = (math.log(mu2[-2] / mu2[-1])
         / math.log(ns[-1] / ns[-2]))
    assert s > 3.0
    ds = [abs(r["n2"]["mu2"]) for r in data["double_scaling"]]
    assert all(a > b for a, b in zip(ds, ds[1:]))


def test_committed_fingerprints_are_the_ones_the_docs_quote():
    """Pin the actual numbers README.md, duality/CHORD.md and road_map.md
    quote as the target laws for the analytic 1/N chord correction.

    test_committed_campaign only asks for monotone decrease and a tail
    slope steeper than 3, which a regression from N^-4 to N^-3.2, or from
    c/(pN) to c/(pN)^1.5, would survive.  These are the fingerprints
    themselves:

      * fixed p = 4 tail slope of |mu_2| is n + 2 for n = 1, 2, 3
        (measured 2.96, 3.90, 4.82) -- so the n = 2 channel is the N^-4
        one the top-level README quotes;
      * on the double-scaling line lambda_chord = 2, |mu_2| * p * N is a
        constant: c ~ 29 for n = 2 (drifting slowly upward) and c ~ 200-205
        for n = 3 (spread 2.1% end to end, and non-monotone).
    """
    with open(HERE / "exact_pair_moments.json") as fh:
        data = json.load(fh)
    rows = data["fixed_p"]
    n1, n2 = rows[-2], rows[-1]
    for n, expected in ((1, 3.0), (2, 4.0), (3, 5.0)):
        a = abs(n1["channels"][f"n={n}"]["mu2"])
        b = abs(n2["channels"][f"n={n}"]["mu2"])
        slope = math.log(a / b) / math.log(n2["N"] / n1["N"])
        assert slope == pytest.approx(expected, abs=0.25), (n, slope)

    cs2 = [abs(r["n2"]["mu2"]) * r["p"] * r["N"]
           for r in data["double_scaling"]]
    cs3 = [abs(r["n3"]["mu2"]) * r["p"] * r["N"]
           for r in data["double_scaling"]]
    assert all(28.0 < c < 31.0 for c in cs2), cs2
    assert all(195.0 < c < 210.0 for c in cs3), cs3
    # n = 2 drifts slowly upward; n = 3 is flat and NOT monotone
    assert all(b > a for a, b in zip(cs2, cs2[1:])), cs2
    assert max(cs3) / min(cs3) < 1.03, cs3
    assert not all(b > a for a, b in zip(cs3, cs3[1:])), cs3


def test_mu0_of_the_n2_channel_changes_sign():
    """The one non-monotone feature, stated in the docs: the n = 2 pair
    mu_0 crosses zero near N = 12 at p = 4, so |mu_0| is NOT monotone
    there.  This is the exact version of an effect first seen in ED."""
    with open(HERE / "exact_pair_moments.json") as fh:
        rows = json.load(fh)["fixed_p"]
    mu0 = [(r["N"], r["channels"]["n=2"]["mu0"]) for r in rows]
    signs = [s for _, s in mu0]
    assert signs[0] > 0 and signs[-1] < 0
    crossings = [n for (n, a), (_, b) in zip(mu0, mu0[1:]) if a * b < 0]
    assert crossings == [10], crossings          # crosses between 10 and 12


def test_full_channel_is_not_conserved_only_the_pair_part_is():
    """The scope of the chord theorem, pinned.

    duality/chord_charges.py proves the i != j PAIR amplitude of every A_n
    is conserved at N = infinity.  The physical channel is

        <A_n A_n> = N <X_i X_i>  +  N(N-1) <X_i X_j>,   X_i = psi_i ad_H^n psi_i

    and the two terms are the SAME order: the per-pair amplitude is smaller
    than the diagonal one by exactly the 1/N the multiplicity restores.  So
    conserving the pair part does not conserve the channel.

    Control: at n = 1 the two parts must cancel identically, because
    A_1 = -2pH is exactly conserved.  At n = 2 they must not.  Without this
    test the docs could drift back to claiming the tower is invisible at
    N = infinity, which it is not."""
    for (p, N) in [(4, 16), (6, 36)]:
        def full(n):
            pair = channel_pair_moments(N, p, n, ks=(0, 2))
            diag = channel_pair_moments(N, p, n, marks=('i', 'i', 'i', 'i'),
                                        ks=(0, 2))
            mu0 = N * diag["mu0"] + N * (N - 1) * pair["mu0"]
            mu2 = N * diag["mu2"] + N * (N - 1) * pair["mu2"]
            return mu0, mu2, pair

        mu0_1, mu2_1, pair1 = full(1)
        assert abs(mu2_1 / mu0_1) < 1e-12, (p, N, mu2_1 / mu0_1)   # A_1 = -2pH
        # the pair part alone is conserved at N = infinity but NOT at finite N,
        # so at these sizes it is nonzero -- the cancellation above is genuine
        assert abs(pair1["mu2"]) > 1e-6

        mu0_2, mu2_2, _ = full(2)
        w2 = mu2_2 / mu0_2
        assert w2 > 0.9, (p, N, w2)      # the full A_2 channel propagates
    # and the width RISES along the double-scaling line rather than falling
    w2s = []
    for (p, N) in [(4, 16), (6, 36)]:
        pair = channel_pair_moments(N, p, 2, ks=(0, 2))
        diag = channel_pair_moments(N, p, 2, marks=('i', 'i', 'i', 'i'),
                                    ks=(0, 2))
        w2s.append((N * diag["mu2"] + N * (N - 1) * pair["mu2"])
                   / (N * diag["mu0"] + N * (N - 1) * pair["mu0"]))
    assert w2s[1] > w2s[0], w2s


def test_expectation_core_m5_vs_literal_enumeration():
    """The m = 5 signature path, against a literal enumeration.

    test_dp_tables_equal_profile_sum anchors the DP only where the closed
    Venn-profile reference exists (m <= 3), and test_k8_word_literal
    reaches m = 4.  But the committed n = 3 / mu_4 column is produced
    entirely by K = 10 words, i.e. m = 5, so a regression in the m >= 4
    factorization or canonicalization would leave every other test green
    while silently moving numbers the docs quote.

    Here the expectation is recomputed by brute force over all
    C(N,p)^5 subset tuples at (N, p) = (5, 2) -- 10^5 tuples per signature.
    """
    from itertools import combinations, product
    from fractions import Fraction
    from exact_pair_moments import _expectation_core

    N, p, m = 5, 2, 5
    subsets = [frozenset(c) for c in combinations(range(N), p)]

    def literal(S, marks_bits):
        n_dist = len(marks_bits)
        dist = list(range(n_dist))          # the distinguished elements
        tot = 0
        for combo in product(subsets, repeat=m):
            par = 0
            for (a, b) in S:
                par ^= len(combo[a] & combo[b]) & 1
            for e, bits in enumerate(marks_bits):
                for al in range(m):
                    if bits[al] and dist[e] in combo[al]:
                        par ^= 1
            tot += -1 if par else 1
        return Fraction(tot, len(subsets) ** m)

    rng = np.random.default_rng(4)
    pairs = list(combinations(range(m), 2))
    checked = 0
    for _ in range(4):
        k = int(rng.integers(1, 4))
        S = frozenset(pairs[i] for i in
                      rng.choice(len(pairs), size=k, replace=False))
        nbits = int(rng.integers(0, 3))
        bits = tuple(tuple(int(b) for b in rng.integers(0, 2, size=m))
                     for _ in range(nbits))
        bits = tuple(b for b in bits if any(b))
        assert _expectation_core(N, p, m, S, bits) == literal(S, bits), \
            (sorted(S), bits)
        checked += 1
    assert checked == 4
