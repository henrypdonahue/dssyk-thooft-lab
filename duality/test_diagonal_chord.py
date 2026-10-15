"""Tests for the N = infinity flavor-DIAGONAL channel (diagonal_chord.py).

The parity law is derived from the finite-N sign law, and gated four
independent ways.

1. INTERNAL AGREEMENT, at machine precision.  The flux walker against the
   literal matching enumeration, against the validated transfer-matrix
   engine on the 2-mark restriction.  Plus the exact psi^2 = 1 collapses,
   which refute the pairing-weight alternative.

2. THE FINITE-N MARK FACTOR is EXACTLY (1 - 2p/N).  That is the
   microscopic origin of q^{Delta}.  The crossed-signature factorization
   converges at the O(1/p) line rate.

3. CONVERGENCE.  The A_1 = -2pH identity's chord shadow, and the exact
   finite-N diagonal moments, both converge onto the chord values along
   lambda_chord = 2 at the O(1/p) rate.  At fixed p, at the word level,
   that sharpens to about 1/N^2.

4. THE G3 NUMBERS.  The diagonal width w2(Delta -> 0, lambda = 2) is the
   closed form (1-q)(5+3q) = 4.674383, which the Richardson of the walker
   values (4.67 +- 0.09) brackets.  The FULL-channel width is a large
   diagonal-against-pair cancellation: its exact N = infinity value is
   full_channel_chord.py's closed form, 1.316906, which the short
   finite-N Richardson estimate kept here (~1.31) brackets.

5. THE ARTIFACT.  diagonal_chord.json is consumed: slices of every block
   are recomputed and compared.
"""

import math
import sys
from fractions import Fraction
from pathlib import Path

import pytest

from chord_moments import MomentTable
from diagonal_chord import (DiagonalMomentTable, bruteforce_parity_word,
                            channel_moments_diag, delta_slope_closed,
                            diag_word_amplitude, finite_n_diag, finite_n_pair,
                            w2_delta0, w2_delta0_closed, word_diag)

HERE = Path(__file__).resolve().parent
Q, DP = math.exp(-1.0), 0.3


def test_walker_matches_bruteforce_enumeration():
    """The flux walker and the literal (K-1)!! matching enumeration share
    only the microscopic parity law; they must agree to machine precision
    on every word shape, including odd-H words (identically zero)."""
    for word in [word_diag(1, 1, 1, 1), word_diag(2, 0, 1, 3),
                 word_diag(0, 2, 2, 2), word_diag(3, 1, 0, 2),
                 word_diag(2, 2, 2, 2),
                 ['H'] * 2 + ['P', 'P'] + ['H'] * 4,
                 ['P'] + ['H'] * 3 + ['P'] + ['H'] * 3]:
        w = diag_word_amplitude(Q, DP, word)
        b = bruteforce_parity_word(Q, DP, word)
        assert w == pytest.approx(b, rel=1e-12, abs=1e-14), word
    assert diag_word_amplitude(Q, DP, word_diag(1, 1, 1, 0)) == 0.0


def test_two_mark_restriction_equals_the_engine():
    """Two marks of one Majorana are the ordinary single-matter-chord
    object.  The walker's t2 must reproduce the validated transfer-matrix
    engine (chord_moments.MomentTable.t2) to machine precision -- this
    anchors the walker's crossing bookkeeping to the engine the module
    chain trusts."""
    eng = MomentTable(Q, DP)
    diag = DiagonalMomentTable(Q, DP)
    for a in range(0, 7):
        for b in range(0, 7 - a):
            if (a + b) % 2:
                continue
            assert diag.t2(a, b) == pytest.approx(eng.t2(a, b), rel=1e-12,
                                                  abs=1e-14), (a, b)


def test_psi_squared_collapses_exactly():
    """psi^2 = 1 is automatic in the parity law: adjacent mark pairs are
    the identity.  These exact collapses are also what refute the
    'q^{2 Delta} for double-crossing chords' pairing-weight alternative,
    which would give q^{2 Delta} != 1 for t4(1,0,1,0)."""
    tab = DiagonalMomentTable(Q, DP)
    # Tr[H psi psi H psi psi] = Tr[H^2] -> 1 in chord units
    assert abs(tab.t4(1, 0, 1, 0) - 1.0) < 1e-14
    # Tr[H^2 psi psi H^2 psi psi] = Tr[H^4] -> 2 + q (Touchard-Riordan)
    assert tab.t4(2, 0, 2, 0) == pytest.approx(2.0 + Q, rel=1e-13)
    assert tab.t4(1, 0, 3, 0) == pytest.approx(2.0 + Q, rel=1e-13)
    # general collapse: marks drop out entirely
    assert tab.t4(3, 0, 3, 0) == pytest.approx(
        diag_word_amplitude(Q, DP, ['H'] * 6), rel=1e-13)


def test_cyclic_invariance():
    """The amplitude is a trace: rotating the word around the circle must
    not change it (the walker evaluates a cut line, so this is a real
    check of the cut-independence argument)."""
    tab = DiagonalMomentTable(Q, DP)
    ref = tab.t4(1, 2, 0, 3)
    assert tab.t4(2, 0, 3, 1) == pytest.approx(ref, rel=1e-12)
    assert tab.t4(0, 3, 1, 2) == pytest.approx(ref, rel=1e-12)
    assert tab.t4(3, 1, 2, 0) == pytest.approx(ref, rel=1e-12)


def test_delta_zero_channel_vanishes():
    """At Delta_psi = 0 the marks decouple, every word in the ad-expansion
    of a channel is the same pure-H moment, and the alternating binomial
    assembly cancels identically -- the diagonal analog of the
    odelta mu_k(qt = 1) = 0 pin."""
    for n in (1, 2):
        ch = channel_moments_diag(Q, 0.0, n, ks=(0, 2, 4))
        assert abs(ch["mu0"]) < 1e-12 * max(ch["gross2"], 1.0)
        assert abs(ch["mu2"]) < 1e-12 * max(ch["gross2"], 1.0)
        assert abs(ch["mu4"]) < 1e-12 * max(ch["gross4"], 1.0)


def test_mark_factor_is_exact_and_factorization_converges():
    """The microscopic origin of the q^{Delta} mark weight, at the
    signature level of the finite-N machinery: (1) the lone-chord mark
    factor is EXACTLY (1 - 2p/N) at every (N, p) -- the parity law's
    weight is its N -> infinity image with O(1/p^2) deviation on the
    lambda_chord = 2 line; (2) for a crossed pair with one marked chord
    the full expectation converges onto the factorized q * q^{Delta} at
    the O(1/p) line rate (the crossing factor dominates the correction)."""
    sys.path.insert(0, str(HERE.parent / "syk-self-averaging"))
    import exact_pair_moments as epm

    rels = []
    for pp in (4, 8, 16):
        N = pp * pp
        e1 = epm.signature_expectation(N, pp, 1, frozenset(), ((1,),))
        assert e1 == Fraction(N - 2 * pp, N)          # exact, every (N, p)
        qd = math.exp(-2.0 / pp)
        assert abs(float(e1) - qd) < 2.5 / pp ** 2
        e2 = float(epm.signature_expectation(
            N, pp, 2, frozenset({(0, 1)}), ((1, 0),)))
        tgt = math.exp(-2.0) * qd
        rels.append(abs(e2 - tgt) / tgt)
        epm._EXP_CACHE.clear()
        epm._DP_CACHE.clear()
    assert rels[0] > rels[1] > rels[2]
    assert 3.0 < 16 * rels[2] < 5.0                   # ~4.2/p, measured


def test_g1_a1_identity_chord_shadow():
    """A_1 = -2pH forces (N-1) mu_k^pair + mu_k^diag = 0 exactly at every
    finite N (machine-gated in test_exact_pair_moments.py).  Its strict
    N = infinity shadow: the chord diagonal must equal the limit of
    -(N-1) mu_k^pair along the line, with the residual falling at the
    O(1/p) line rate.  Measured residuals (mu2): 0.35, 0.23, 0.17."""
    resid = {2: [], 4: []}
    for (pp, N) in [(4, 16), (6, 36), (8, 64)]:
        ch = channel_moments_diag(math.exp(-2.0), 1.0 / pp, 1, ks=(0, 2, 4))
        pr = finite_n_pair(N, pp, 1, ks=(0, 2, 4))
        for k in (2, 4):
            resid[k].append(abs(ch[f"mu{k}"] + (N - 1) * pr[f"mu{k}"])
                            / abs(ch[f"mu{k}"]))
    for k in (2, 4):
        assert resid[k][0] > resid[k][1] > resid[k][2], resid[k]
    assert resid[2][-1] < 0.20
    assert resid[4][-1] < 0.22


def test_g2_finite_n_converges_onto_the_chord_diagonal():
    """The exact finite-N diagonal moments approach the chord values along
    lambda_chord = 2 at the O(1/p) line rate (p * rel_dev ~ 1.3-1.5,
    falling), and at fixed p with the lattice lambda matched exactly the
    WORD-level deviation falls like ~1/N^2 -- the assembled mu_2 at fixed
    p is itself O(lambda^2) ~ 1/N^2, so the sharp gate lives on the O(1)
    words."""
    rels = []
    for (pp, N) in [(4, 16), (6, 36), (8, 64), (10, 100)]:
        ch = channel_moments_diag(math.exp(-2.0), 1.0 / pp, 2, ks=(0, 2))
        fN = finite_n_diag(N, pp, 2, ks=(0, 2))
        rels.append(abs(ch["mu2"] - fN["mu2"]) / abs(ch["mu2"]))
    assert all(b < a for a, b in zip(rels, rels[1:])), rels
    assert rels[-1] < 0.15
    assert all(1.0 < pp * r < 1.7
               for pp, r in zip((4, 6, 8, 10), rels)), rels

    sys.path.insert(0, str(HERE.parent / "syk-self-averaging"))
    import exact_pair_moments as epm
    devs = []
    for N in (64, 256, 1024):
        lamN = 32.0 / N
        ch = diag_word_amplitude(math.exp(-lamN), 0.25, word_diag(1, 1, 1, 1))
        fN = epm.t4bar(N, 4, 1, 1, 1, 1, marks=('i', 'i', 'i', 'i'))
        devs.append(abs(ch - fN))
        epm._DP_CACHE.clear()
        epm._EXP_CACHE.clear()
    assert devs[1] < devs[0] / 8.0 and devs[2] < devs[1] / 8.0, devs
    assert devs[-1] < 5e-4


def test_g3_diagonal_width_and_full_channel_context():
    """The payoff numbers, kept honestly apart.  The DIAGONAL width at
    lambda_chord = 2 rises 3.46 -> 4.08 along Delta = 1/p and
    extrapolates to w2 = 4.67 +- 0.09 at Delta -> 0 -- the exact
    N = infinity chord value this module supplies.  The FULL-channel
    width is a large diagonal-vs-pair cancellation; the exact finite-N
    series 1.044 -> 1.194 Richardson-extrapolates to ~1.31 in 1/p (data,
    not a chord computation -- the missing chord piece is the pair
    channel's 1/N coefficient)."""
    w2s = [channel_moments_diag(math.exp(-2.0), 1.0 / pp, 2,
                                ks=(0, 2))["w2"]
           for pp in (4, 6, 8, 10)]
    assert all(b > a for a, b in zip(w2s, w2s[1:])), w2s
    assert w2s[0] == pytest.approx(3.462, abs=0.01)
    assert w2s[-1] == pytest.approx(4.075, abs=0.01)

    d0 = w2_delta0(2.0)
    assert d0["w2_delta0"] == pytest.approx(4.67, abs=0.15)
    assert d0["extrap_spread"] < 0.15
    # the exact limit, and the Richardson brackets it within its spread
    q = math.exp(-2.0)
    assert d0["w2_delta0_closed"] == pytest.approx((1 - q) * (5 + 3 * q), rel=1e-14)
    assert abs(d0["w2_delta0"] - d0["w2_delta0_closed"]) < d0["extrap_spread"]

    full = []
    for (pp, N) in [(4, 16), (6, 36), (8, 64), (10, 100)]:
        pr = finite_n_pair(N, pp, 2, ks=(0, 2))
        dg = finite_n_diag(N, pp, 2, ks=(0, 2))
        full.append([pp, (N * dg["mu2"] + N * (N - 1) * pr["mu2"])
                     / (N * dg["mu0"] + N * (N - 1) * pr["mu0"])])
    # the pinned series (also guarded in test_exact_pair_moments.py)
    for (pp, w2), want in zip(full, (1.044, 1.118, 1.164, 1.194)):
        assert w2 == pytest.approx(want, abs=0.01), (pp, w2)
    rich = [w2 + (w2 - w1) * (1.0 / pb) / (1.0 / pa - 1.0 / pb)
            for (pa, w1), (pb, w2) in zip(full, full[1:])]
    assert 1.25 < rich[-1] < 1.40, rich
    # and the exact N = infinity value it was estimating (full_channel_chord)
    from full_channel_chord import channel
    assert abs(rich[-1] - channel(2.0)["w2"]) < max(rich) - min(rich) + 0.01


def test_delta_slopes_are_closed_forms():
    """mu_k^diag = lambda P_K(q) Delta + O(Delta^2) with P_K depending on
    K = k + 2n only: the walker's finite-difference slopes at two lambdas
    match delta_slope_closed for K = 4, 6, 8 across three (n, k) each."""
    for lam in (1.0, 2.0):
        q = math.exp(-lam)
        d0 = 1.0 / 2048
        for (n, k) in [(1, 2), (2, 0), (2, 2), (1, 4), (3, 0), (3, 2), (2, 4),
                       (1, 6)]:
            K = k + 2 * n
            ks = tuple(sorted({0, 2, k}))
            r1 = channel_moments_diag(q, 2 * d0, n, ks=ks)
            r2 = channel_moments_diag(q, d0, n, ks=ks)
            slope = 2 * (r2[f"mu{k}"] / d0) - r1[f"mu{k}"] / (2 * d0)
            want = lam * delta_slope_closed(K, q)
            assert abs(slope - want) < 5e-5 * abs(want), (n, k, slope, want)
    assert w2_delta0_closed(2.0) == pytest.approx(
        delta_slope_closed(6, math.exp(-2.0)) / delta_slope_closed(4, math.exp(-2.0)))


def test_channel_moments_diag_accepts_any_ks():
    """Regression: asking for a single k used to raise KeyError."""
    r = channel_moments_diag(math.exp(-2.0), 0.1, 2, ks=(4,))
    assert set(r) == {"mu4", "gross4"}
    r = channel_moments_diag(math.exp(-2.0), 0.1, 2, ks=(0, 2))
    assert "w2" in r and "kurt" not in r


@pytest.fixture(scope="module")
def art():
    import json
    with open(HERE / "diagonal_chord.json") as fh:
        return json.load(fh)


def test_artifact_points_recompute(art):
    """diagonal_chord.json: the campaign points at two lambdas, the line
    rows, the G1 residuals and the G3 numbers are recomputed exactly."""
    for pt in art["points"]:
        if pt["lam"] not in (2.0, 0.5):
            continue
        for n in (1, 2, 3):
            live = channel_moments_diag(pt["q"], pt["delta_psi"], n,
                                        ks=(0, 2, 4, 6))
            for key, val in pt["channels"][f"n={n}"].items():
                if isinstance(val, float) and math.isfinite(val):
                    assert live[key] == pytest.approx(val, rel=1e-12, abs=1e-14), \
                        (pt["lam"], n, key)
    row = art["line_lambda2"][0]
    assert row["p"] == 4 and row["N"] == 16
    dg = finite_n_diag(16, 4, 2, ks=(0, 2))
    assert dg["mu2"] == pytest.approx(row["finiteN_mu2"], rel=1e-12)
    ch = channel_moments_diag(math.exp(-2.0), 0.25, 2, ks=(0, 2))
    assert ch["mu2"] == pytest.approx(row["chord_mu2"], rel=1e-12)
    g3 = art["g3"]
    assert g3["delta0"]["w2_delta0_closed"] == pytest.approx(
        w2_delta0_closed(2.0), rel=1e-14)
    assert abs(g3["delta0"]["w2_delta0"] - g3["delta0"]["w2_delta0_closed"]) \
        < g3["delta0"]["extrap_spread"]
    from full_channel_chord import channel
    assert g3["full_channel_w2_inf_closed_form"] == pytest.approx(
        channel(2.0)["w2"], rel=1e-14)
    assert abs(g3["full_channel_w2_inf_estimate"]
               - g3["full_channel_w2_inf_closed_form"]) < 0.01
