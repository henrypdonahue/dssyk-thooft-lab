"""Tests for singlet_tower.py: the full A_n channels at N = infinity as
matter chords of dimension 2..n.

The n = 2 statement is a theorem and is tested as one: symbolically, by
its proof's two matching classes, against the published O_Delta closed
forms, against the exact finite-N full channel and against dense ED.
The k = 1, 2, 3 decompositions are exact identities checked in computer
algebra."""

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

import singlet_tower as st
from chord import bruteforce_word
from chord_moments import odelta_moments
from full_channel_chord import (dynamic_channel, mu0_closed, mu2_closed,
                                mu4_closed, word_amplitude)

HERE = Path(__file__).resolve().parent
POINTS = [(Fraction(2), Fraction(1, 7)), (Fraction(1, 3), Fraction(3, 5)),
          (Fraction(5, 2), Fraction(2, 9))]


@pytest.fixture(scope="module")
def art():
    with open(HERE / "singlet_tower.json") as fh:
        return json.load(fh)


# --------------------------------------------------------------------------
# the engine against independent code
# --------------------------------------------------------------------------
def test_omega0_is_the_full_channel_word():
    """16 x <H^A Omega_0 H^B Omega_0> is full_channel_chord.word_amplitude,
    an independent enumeration of the same chord rule."""
    for lam, q in POINTS[:2]:
        T = st.TowerChord(lam, q)
        for (A, B) in [(0, 0), (1, 1), (2, 0), (2, 2), (3, 1), (4, 0), (5, 1)]:
            assert 16 * T.omega_word(0, 0, A, B) == word_amplitude(A, B, lam, q)


def test_matter_and_h_words_match_the_chord_engine():
    """G_m(A, B) and <H^k> against chord.bruteforce_word (which shares
    nothing with the count tables here)."""
    q = 0.37
    T = st.TowerChord(1.0, q)
    for m in (2, 3, 4):
        for (A, B) in [(0, 0), (1, 1), (2, 0), (3, 1), (2, 2), (4, 2)]:
            word = ['H'] * A + ['M'] + ['H'] * B + ['M']
            ref = bruteforce_word(q, word, {'M': float(m)})
            assert T.matter_word(A, B, m) == pytest.approx(ref, rel=1e-13)
    for k in (2, 4, 6, 8):
        assert T.h_moment(k) == pytest.approx(
            bruteforce_word(q, ['H'] * k, {}), rel=1e-13)
    # q-Hermite orthogonality for the chord measure: <h_i h_j> = [i]_q! d_ij
    for i in range(5):
        for j in range(i % 2, 5, 2):
            hi, hj = T.hermite(i), T.hermite(j)
            ip = sum(a * b * T.h_moment(x + y) for x, a in hi.items()
                     for y, b in hj.items())
            want = T.hermite_norm(i) if i == j else 0.0
            assert ip == pytest.approx(want, abs=1e-12)


def test_identity_coefficients_agree_across_modules():
    """The exact operator identity behind the chord rule is asserted
    against dense matrices in syk-self-averaging/test_conserved_points.py;
    the coefficient table used here is the same one."""
    sys.path.insert(0, str(HERE.parent / "syk-self-averaging"))
    import conserved_points as cp
    for n in (2, 3, 4, 5):
        assert st.identity_coefficients(n) == cp.identity_coefficients(n)
        assert st.kappa(n) == {T: c for T, c in
                               st.identity_coefficients(n).items()
                               if len(T) == 2}
    assert st.a_terms(3) == [(-8, 0, 0, 1), (4, 0, 1, 0), (4, 1, 0, 0)]


def test_enclosed_remainder_is_the_krawtchouk_expansion():
    """The O(1) remainder rule.  A chord of p indices that encloses a
    product of grade s = 2p - 2x (two enclosed chords sharing x indices)
    averages to f(x) = E[(-1)^|c cap S|], an exact Krawtchouk ratio, and
    p [f(x)/f(0) - 1] -> 2 lambda x on a double-scaling line.  That limit
    is the 2 lambda X-weight of _omega_terms."""
    from math import comb

    def f(N, p, s):
        return Fraction(sum((-1) ** j * comb(s, j) * comb(N - s, p - j)
                            for j in range(p + 1)), comb(N, p))

    for lam in (1.0, 2.0):
        for x in (1, 2):
            devs = []
            for p in (20, 40, 80, 160):
                N = int(round(2 * p * p / lam))
                r = float(p * (f(N, p, 2 * p - 2 * x) / f(N, p, 2 * p) - 1))
                devs.append(abs(r - 2 * lam * x))
            # an O(1/p) approach: the deviation halves when p doubles
            assert all(b < 0.6 * a for a, b in zip(devs, devs[1:])), devs
            assert devs[-1] < 0.1 * 2 * lam * x, devs
    # and the leading term itself: f(0) -> q^2
    assert float(f(2 * 160 * 160, 160, 320)) == pytest.approx(
        math.exp(-2.0), rel=0.02)


# --------------------------------------------------------------------------
# n = 2: the theorem
# --------------------------------------------------------------------------
def test_a2_theorem_symbolic():
    """<H^A Omega_0 H^B Omega_0> = [mu(1-q)/(1+q)]^2 <H^(A+B)(H^2-1)^2>
    + norm_2 G_2(A, B), identically in (lambda, q), every word through
    six H's."""
    import sympy as sp
    lam, q = sp.symbols("lam q", positive=True)
    T = st.TowerChord(lam, q)
    mu = lam / 2
    kap = mu * (1 - q) / (1 + q)
    for L in (0, 2, 4, 6):
        for A in range(L + 1):
            B = L - A
            hh = T.h_moment(L + 4) - 2 * T.h_moment(L + 2) + T.h_moment(L)
            rhs = kap ** 2 * hh + st.norm2(lam, q) * T.matter_word(A, B, 2)
            assert sp.cancel(T.omega_word(0, 0, A, B) - rhs) == 0, (A, B)
    # the projections: Omega_0 overlaps h_2 only
    assert sp.cancel(T.projection(0, 2) - mu * (1 - q)) == 0
    assert sp.cancel(T.projection(0, 0)) == 0
    assert sp.cancel(st.norm2(lam, q) - (mu * (1 - q) + mu ** 2 * (1 + q)
                                         - mu ** 2 * (1 - q) ** 2 / (1 + q))) == 0


def test_a2_proof_classes():
    """The two classes of the proof, enumerated separately.  Class I (both
    chords of one insertion end on the other) is [mu(1-q) + mu^2(1+q)]
    G_2; class II is [mu(1-q)/(1+q)]^2 times class II of the plain word
    <H^A (H^2-1) H^B (H^2-1)>, by the slot-swap involution."""
    lam, q = Fraction(3, 2), Fraction(2, 7)
    mu = lam / 2
    T = st.TowerChord(lam, q)
    for (A, B) in [(0, 0), (1, 1), (2, 0), (2, 2), (3, 1), (4, 2)]:
        K = A + B + 4
        s1, s2 = (A, A + 1), (A + 2 + B, A + 3 + B)
        cls1 = cls2 = plain2 = 0
        for of, cross in st._diagrams(K):
            c1 = [of[s] for s in s1]
            c2 = [of[s] for s in s2]
            if c1[0] == c1[1] or c2[0] == c2[1]:
                continue                      # self-contraction: the -1 of H^2-1
            e, f = st._pair(*c1), st._pair(*c2)
            nx = len(cross)
            if e == f:
                cls1 += q ** nx * ((mu * mu - mu) if e in cross else (mu + mu * mu))
            else:
                sign = (-1) ** ((e in cross) + (f in cross))
                cls2 += sign * mu * mu * q ** nx
                plain2 += q ** nx
        assert cls1 == (mu * (1 - q) + mu * mu * (1 + q)) * T.matter_word(A, B, 2)
        assert cls2 == (mu * (1 - q) / (1 + q)) ** 2 * plain2
        assert cls1 + cls2 == T.omega_word(0, 0, A, B)


@pytest.mark.parametrize("lam", [0.1, 0.5, 1.0, 2.0, 4.0])
def test_a2_dynamical_channel_is_a_delta2_matter_chord(lam):
    """mu_0 splits as stated, w2_dyn = 2(1 - q^2), and the kurtosis and
    the sixth moment are the published O_Delta closed forms at Delta = 2
    (chord_moments.odelta_moments)."""
    q = math.exp(-lam)
    s = st.a2_split(lam)
    assert s["mu0"] == pytest.approx(mu0_closed(lam, q), rel=1e-13)
    assert s["mu2"] == pytest.approx(mu2_closed(lam, q), rel=1e-13)
    assert s["w2_dynamic"] == pytest.approx(mu2_closed(lam, q) / s["mu0_dynamic"],
                                            rel=1e-13)
    od = odelta_moments(q, 2.0)
    assert s["w2_dynamic"] == pytest.approx(od["mu2"], rel=1e-13)
    assert mu4_closed(lam, q) / s["mu0_dynamic"] == pytest.approx(od["mu4"],
                                                                 rel=1e-12)
    T = st.TowerChord(lam, q)
    assert T.a_moment(2, 6) / s["mu0_dynamic"] == pytest.approx(od["mu6"],
                                                                rel=1e-10)
    d = dynamic_channel(lam)
    assert d["w2_dynamic"] == pytest.approx(2 * (1 - q * q), rel=1e-13)
    assert d["kurt_dynamic"] == pytest.approx(od["kurt"], rel=1e-12)
    assert d["stationary_fraction"] == pytest.approx(s["stationary_fraction"],
                                                     rel=1e-13)


def test_a2_headline_numbers():
    """lambda_chord = 2: 32.9 % of mu_0 is conserved H^2 weight; the
    dynamical width is 1.963369 where the raw one is 1.316906.  The raw
    width's maximum near lambda ~ 1 is the dilution: the dynamical width
    is monotone and saturates at 2."""
    s = st.a2_split(2.0)
    assert s["stationary_fraction"] == pytest.approx(0.329262, abs=1e-6)
    assert s["w2_dynamic"] == pytest.approx(1.963369, abs=1e-6)
    assert s["w2_raw"] == pytest.approx(1.316906, abs=1e-6)
    lams = [0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
    dyn = [st.a2_split(x)["w2_dynamic"] for x in lams]
    raw = [st.a2_split(x)["w2_raw"] for x in lams]
    assert all(b > a for a, b in zip(dyn, dyn[1:])) and dyn[-1] < 2.0
    assert max(raw) == raw[4] and raw[-1] < raw[4] > raw[0]
    # semiclassics: 2(1 - q^2) = 4 lam - 4 lam^2 + ..., the Delta = 2 slope
    assert st.a2_split(1e-4)["w2_dynamic"] == pytest.approx(4e-4, rel=1e-3)


# --------------------------------------------------------------------------
# the decomposition of the building blocks
# --------------------------------------------------------------------------
def _symbolic_decomposition(kmax, check_len):
    import sympy as sp
    lam, q = sp.symbols("lam q", positive=True)
    T = st.TowerChord(lam, q)
    d = T.decompose(kmax, check_len=check_len)
    assert d["max_residual"] == 0
    closed = st.omega_kernels_closed(lam, q)
    for k in range(kmax + 1):
        assert set(d["kernels"][k]) == set(closed[k])
        for m in closed[k]:
            assert set(d["kernels"][k][m]) == set(closed[k][m])
            for key, v in closed[k][m].items():
                assert sp.cancel(d["kernels"][k][m][key] - v) == 0, (k, m, key)
    for m, f in ((2, st.norm2), (3, st.norm3), (4, st.norm4)):
        if m <= kmax + 2:
            assert sp.cancel(d["norm"][m] - f(lam, q)) == 0, m
    return d


def test_decomposition_symbolic_k1():
    """Omega_0 and Omega_1, identically in (lambda, q): Omega_1 =
    [q/(1+q^2)] {H, O^_2} + O^_3 modulo polynomials in H, with zero
    residual on every two-point word through four H's."""
    d = _symbolic_decomposition(1, 4)
    assert d["n_checked"] >= 20


@pytest.mark.slow
def test_decomposition_symbolic_k2():
    """Omega_2 as well, identically in (lambda, q): the kernels are solved
    linearly and every two-point word through four H's then vanishes;
    kernels and norms equal the closed forms.  (About a minute of sympy;
    the fast suite checks the same statement in exact rational arithmetic
    at three generic points, test_decomposition_exact_k3.)"""
    d = _symbolic_decomposition(2, 4)
    assert d["n_checked"] >= 60


@pytest.mark.parametrize("lam,q", POINTS)
def test_decomposition_exact_k3(lam, q):
    """Through Omega_3 in exact rational arithmetic at generic points
    (lambda and q independent): zero residual on every checked word, a
    nonzero Delta = 5 chord, and the one-strand-crossed coefficient of
    the subleading sector, q [k]_q / (1 + q^(k+1)), at k = 1, 2, 3."""
    T = st.TowerChord(lam, q)
    d = T.decompose(3, check_len=3)
    assert d["max_residual"] == 0 and d["n_checked"] >= 80
    assert all(d["norm"][m] > 0 for m in (2, 3, 4, 5))
    for k in (1, 2, 3):
        qk = sum(q ** j for j in range(1, k + 1))
        f = d["kernels"][k][k + 1]
        assert f == {(1, 0): qk / (1 + q ** (k + 1)),
                     (0, 1): qk / (1 + q ** (k + 1))}
    # every kernel is symmetric (the blocks are Hermitian) and has the
    # block's H-parity
    for k in range(4):
        for m, f in d["kernels"][k].items():
            for (i, j), c in f.items():
                assert f[(j, i)] == c
                assert (i + j - (k + 2 - m)) % 2 == 0 and i + j <= k + 2 - m


@pytest.mark.slow
def test_decomposition_exact_k4():
    """Omega_4 and the Delta = 6 chord, same statement one level up
    (about 40 s: the K = 14 word tables)."""
    lam, q = POINTS[0]
    T = st.TowerChord(lam, q)
    d = T.decompose(4, check_len=2)
    assert d["max_residual"] == 0 and d["n_checked"] >= 110
    assert d["norm"][6] > 0
    q4 = sum(q ** j for j in range(1, 5)) / (1 + q ** 5)
    assert d["kernels"][4][5] == {(1, 0): q4, (0, 1): q4}


def test_block_projections_stop_at_k_plus_2():
    """Omega_k has k + 2 chord ends, so it overlaps h_j only for
    j <= k + 2 with j = k mod 2: the stationary part is a polynomial of
    degree k + 2."""
    lam, q = POINTS[0]
    T = st.TowerChord(lam, q)
    for k in range(3):
        for j in range(0, k + 7):
            c = T.projection(k, j)
            if j > k + 2 or (j - k) % 2:
                assert c == 0, (k, j)
    assert T.projection(0, 2) == lam * (1 - q) / 2
    assert T.projection(1, 1) == lam


@pytest.mark.parametrize("lam,q", POINTS[:2])
def test_channel_moments_from_blocks_and_from_sectors(lam, q):
    """mu_k of the full A_n channels three ways: assembled from the Omega
    words, from the matter-chord sectors, and the closed forms."""
    T = st.TowerChord(lam, q)
    assert T.a_moment(2, 2) == mu2_closed(lam, q)
    assert T.a_moment(2, 4) == mu4_closed(lam, q)
    for n, k, f in ((3, 2, st.a3_mu2), (3, 4, st.a3_mu4), (4, 2, st.a4_mu2)):
        assert T.a_moment(n, k) == f(lam, q) == T.sector_moment(n, k), (n, k)
    assert T.a_moment(5, 2) == T.sector_moment(5, 2)
    assert T.sector_moment(2, 0) == 16 * st.norm2(lam, q)
    with pytest.raises(ValueError):
        T.a_moment(3, 0)


def test_a3_kernel_closed_form():
    """A_3 = (poly in H) + nu0 [6 [H, O2] - 2 (1-q)^2/(1+q^2) {H, O2}]
    + 4 O^_3 in units of the un-normalized chords: the commutator
    coefficient is the exact identity A_3 - A_3^dag = 3 [H, A_2], the
    anticommutator admixture vanishes like lambda^2."""
    for lam, q in POINTS:
        T = st.TowerChord(lam, q)
        k = T.a_kernels(3)
        alpha = (k[2][(1, 0)] - k[2][(0, 1)]) / 2
        beta = (k[2][(1, 0)] + k[2][(0, 1)]) / 2
        assert alpha == 6
        assert beta == -2 * (1 - q) ** 2 / (1 + q * q)
        assert k[3] == {(0, 0): 4}
    assert st.TowerChord(Fraction(2), Fraction(1, 7)).a_kernels(2) == \
        {2: {(0, 0): 4}}


def test_norm_scalings_at_small_lambda():
    """norm_2 ~ lam^2, norm_3 ~ 2 lam^2, and the higher chords are
    suppressed: norm_4 ~ (52/5) lam^4 and norm_5 ~ lam^4 too."""
    lam = 1e-3
    q = math.exp(-lam)
    assert st.norm2(lam, q) / lam ** 2 == pytest.approx(1.0, abs=1e-3)
    assert st.norm3(lam, q) / lam ** 2 == pytest.approx(2.0, abs=5e-3)
    lam = 1e-2
    q = math.exp(-lam)
    assert st.norm4(lam, q) / lam ** 4 == pytest.approx(52 / 5, rel=0.05)
    n5 = []
    for lam in (Fraction(1, 20), Fraction(1, 40)):
        q = Fraction(math.exp(-float(lam))).limit_denominator(10 ** 9)
        T = st.TowerChord(lam, q)
        T.decompose(3, check_len=1)
        n5.append(float(T.norm[5]) / float(lam) ** 4)
    assert 0.5 < n5[1] / n5[0] < 1.6 and n5[0] > 0      # lam^4, not lam^2


def test_spectral_density_reproduces_the_moments():
    """The sector decomposition pushed through the closed-form vertex of
    review (2.30): int S = dynamical mu_0 and int omega^2 S = mu_2, for
    A_2 and A_3.  A_2's density is monotone (no interior maximum); A_3's
    has the single broad maximum of a derivative kernel, and no nodes."""
    lam = 1.0
    q = math.exp(-lam)
    band = 2.0 / math.sqrt(1.0 - q)
    om = np.linspace(0.0, 2.0 * band, 241)
    T = st.TowerChord(lam, q)
    for n in (2, 3):
        S = np.array(st.spectral_density(n, om, lam, n_theta=400)["total"])
        m0 = 2.0 * np.trapezoid(S, om)
        m2 = 2.0 * np.trapezoid(S * om ** 2, om)
        assert m0 == pytest.approx(float(T.sector_moment(n, 0)), rel=2e-4)
        assert m2 == pytest.approx(float(T.a_moment(n, 2)), rel=5e-4)
        assert np.all(S >= -1e-12)
        peak = int(np.argmax(S))
        if n == 2:
            assert peak == 0 and np.all(np.diff(S) <= 1e-9)
        else:
            assert 0.5 < om[peak] < band


# --------------------------------------------------------------------------
# the committed artifact
# --------------------------------------------------------------------------
def test_artifact_exact_checks(art):
    """Three generic rational points through Omega_4 (levels n <= 6), and
    one pushed to words of six H's for the blocks with closed forms."""
    for row in art["exact_checks"]:
        assert row["kmax"] == 4
        assert row["max_residual"] == "0" and row["n_checked"] >= 110
        assert row["closed_kernels_match"] and row["closed_norms_match"]
        assert row["moments_match"]
    assert len(art["exact_checks"]) == 3
    deep = art["exact_check_deep"]
    assert deep["max_residual"] == "0" and deep["check_len"] == 6
    assert deep["n_checked"] >= 120


def test_artifact_norm_scaling(art):
    """norm_m / lambda^(2 floor(m/2)) settles toward a constant as
    lambda -> 0: 1, 2 and 52/5 for m = 2, 3, 4 (closed forms), and finite
    limits for m = 5, 6 (successive changes shrink).  One row is
    recomputed through Omega_3."""
    rows = sorted(art["norm_scaling"], key=lambda r: -r["lam"])
    last = rows[-1]["norm_over_power"]
    assert last["2"] == pytest.approx(1.0, abs=0.01)
    assert last["3"] == pytest.approx(2.0, abs=0.05)
    assert last["4"] == pytest.approx(52 / 5, rel=0.05)
    for m in ("5", "6"):
        seq = [r["norm_over_power"][m] for r in rows]
        steps = [b - a for a, b in zip(seq, seq[1:])]
        assert all(s > 0 for s in steps)
        assert all(b < a for a, b in zip(steps, steps[1:])), (m, steps)
    live = st.norm_scaling(lams=(Fraction(1, 10),), kmax=3)[0]
    stored = next(r for r in rows if r["lam"] == 0.1)
    for m in ("2", "3", "4", "5"):
        assert live["norm_over_power"][m] == pytest.approx(
            stored["norm_over_power"][m], rel=1e-12)


def test_artifact_spectra(art):
    """The stored spectral densities close on the exact moments; A_2's is
    monotone and the higher levels have one interior maximum."""
    sp = art["spectra"]
    for n, lev in sp["levels"].items():
        assert lev["mu0_integral"] == pytest.approx(lev["mu0_exact"], rel=1e-4)
        assert lev["mu2_integral"] == pytest.approx(lev["mu2_exact"], rel=1e-4)
        assert lev["min_value"] >= 0.0
        tot = np.array(lev["total"])
        rising = np.diff(tot) > 0
        flips = int(np.sum(rising[1:] != rising[:-1]))
        if n == "2":
            assert lev["omega_at_max"] == 0.0 and not rising.any()
        else:
            assert 0.5 < lev["omega_at_max"] < 1.2 * sp["band"]
            assert flips == 1                    # up, then down: one maximum


def test_artifact_scan_recomputes(art):
    """One scan row is recomputed live; the sector shares say how little
    of A_4 and A_5 is new."""
    row = next(r for r in art["scan"] if r["lam"] == 2.0)
    live = st.scan_row(2.0)
    for n in ("2", "3", "4", "5"):
        for key in ("mu0_dynamic", "mu2", "w2_dynamic", "kurt_dynamic"):
            assert live["channels"][n][key] == pytest.approx(
                row["channels"][n][key], rel=1e-9), (n, key)
        assert sum(row["channels"][n]["sector_share"].values()) == \
            pytest.approx(1.0, abs=1e-12)
    assert row["channels"]["2"]["sector_share"] == {"2": 1.0}
    assert row["channels"]["3"]["sector_share"]["3"] == pytest.approx(0.1489, abs=1e-3)
    assert row["channels"]["4"]["sector_share"]["4"] == pytest.approx(0.0164, abs=1e-3)
    assert row["channels"]["5"]["sector_share"]["5"] == pytest.approx(0.0018, abs=2e-4)
    assert row["channels"]["6"]["sector_share"]["6"] == pytest.approx(0.0002, abs=5e-5)
    # outside O_2 and O_3: between 0.9 % and 5.4 % of A_4..A_6 at every
    # lambda scanned (0.05 to 4)
    for r in art["scan"]:
        for n in ("4", "5", "6"):
            sh = r["channels"][n]["sector_share"]
            assert 0.008 < 1 - sh["2"] - sh["3"] < 0.055, (r["lam"], n, sh)


def test_artifact_share_scaling(art):
    """Toward small lambda the even chords fade out of A_3..A_6 linearly
    in lambda and the odd ones stay: A_3, A_4 -> O_3; A_5 -> 99.5 % O_3 +
    0.53 % O_5; A_6 -> 99.2 % O_3 + 0.80 % O_5.  One row is recomputed
    through Omega_3."""
    rows = sorted(art["share_scaling"], key=lambda r: -r["lam"])
    a, b = rows[-3], rows[-2]            # lambda = 1/160 and 1/320
    assert a["lam"] == pytest.approx(2 * b["lam"])
    for n, m in (("3", "2"), ("4", "2"), ("4", "4"), ("5", "2"), ("5", "4"),
                 ("6", "2"), ("6", "4")):
        ratio = b["shares"][n][m] / a["shares"][n][m]
        assert ratio == pytest.approx(0.5, abs=0.03), (n, m, ratio)   # ~ lambda
    last = rows[-1]["shares"]
    assert last["3"]["3"] > 0.995 and last["4"]["3"] > 0.995
    assert last["5"]["5"] == pytest.approx(0.0053, abs=2e-4)
    assert last["6"]["5"] == pytest.approx(0.0080, abs=3e-4)
    assert last["6"]["6"] < 1e-4
    # the odd chord's share RISES toward its limit as lambda falls
    o5 = [r["shares"]["5"]["5"] for r in rows]
    assert all(y > x for x, y in zip(o5, o5[1:]))
    live = st.share_scaling(lams=(Fraction(1, 10),), kmax=3)[0]
    stored = next(r for r in rows if r["lam"] == 0.1)
    for n in ("3", "4", "5"):
        for m, v in live["shares"][n].items():
            assert v == pytest.approx(stored["shares"][n][m], rel=1e-10)


def test_a2_gate_against_exact_finite_n(art):
    """The exact finite-N full channel with its H^2 component removed
    (exact rationals, rows of full_channel_chord.json) extrapolates in 1/p
    to the theorem on both double-scaling lines; one row is recomputed."""
    for key, tol in (("lam2", 5e-5), ("lam1", 5e-4)):
        g = art["a2_gate"][key]
        for name in ("mu0_stationary", "mu0_dynamic", "w2_dynamic"):
            lim, cf = g["limits"][name]["central"], g["closed_form"][name]
            assert abs(lim - cf) < tol * abs(cf), (key, name, lim, cf)
    with open(HERE / "full_channel_chord.json") as fh:
        rows = json.load(fh)["gates"]["lam2"]["rows"]
    r = rows[2]
    live = st.finite_n_a2_split(r["N"], r["p"], Fraction(r["mu0_exact"]),
                                Fraction(r["mu2_exact"]))
    stored = art["a2_gate"]["lam2"]["rows"][2]
    assert float(live["w2_dynamic"]) == pytest.approx(stored["w2_dynamic"],
                                                      rel=1e-13)
    # the removal matters: raw and dynamical widths differ by ~50 % there
    assert stored["w2_dynamic"] > 1.4 * stored["w2_raw"]


def test_higher_gates_against_exact_finite_n(art):
    """mu_4(A_2), mu_2 and mu_4 of A_3, mu_2(A_4): the exact finite-N
    series approach the N = infinity closed forms, and their 1/p limits
    land within the stated tolerance (the tolerance tracks how deep each
    series could be pushed; the spreads are sensitivity estimates)."""
    tol = {"n2_k4_lam2": 3e-4, "n3_k2_lam2": 3e-4, "n3_k4_lam2": 1e-2,
           "n4_k2_lam2": 1e-2, "cross42_k2_lam2": 1e-4,
           "cross42_k2_lam1": 1e-4}
    assert set(art["gates"]) == set(tol)
    for key, g in art["gates"].items():
        vals = [r["mu"] for r in g["rows"]]
        assert all(b > a for a, b in zip(vals, vals[1:])), key   # rising
        assert vals[-1] < g["closed_form"], key                   # from below
        assert g["rel_dev"] < tol[key], (key, g["rel_dev"])
        assert g["closed_form"] == pytest.approx(
            st.closed_form_target(g["n"], g["k"], g["lam"]), rel=1e-12)
    # three cheap points recomputed exactly
    for key, p in (("n2_k4_lam2", 4), ("n3_k2_lam2", 4)):
        g = art["gates"][key]
        r = next(x for x in g["rows"] if x["p"] == p)
        assert str(st.finite_n_full(r["N"], p, g["n"], g["k"])) == r["mu_exact"]
    g = art["gates"]["cross42_k2_lam1"]
    r = next(x for x in g["rows"] if x["p"] == 4)
    assert str(st.finite_n_cross(r["N"], 4, 4, 2, 2)) == r["mu_exact"]
    # the cross moment with equal levels is the ordered half of mu_k
    assert st.finite_n_cross(16, 4, 2, 2, 2) == st.finite_n_full(16, 4, 2, 2)
    T = st.TowerChord(Fraction(2), Fraction(1, 7))
    assert T.a_cross_moment(3, 3, 2) == T.a_moment(3, 2)


def test_enclosed_remainder_is_needed_at_n4(art):
    """The O(1) remainder of the self-contracted Omega_2 term is not
    optional.  Switched off, <ad^2(A_4) A_2^dag> would be 150.647 instead
    of 155.029 at lambda_chord = 1 (2.8 % apart) and mu_2(A_4) 6531.9
    instead of 6593.9 at lambda_chord = 2.  The exact finite-N cross
    moment on lambda_chord = 1 (K = 8 words, p <= 20) lands on the value
    WITH the remainder and excludes the one without by many spreads.  The
    control values are recomputed live."""
    ctl = art["remainder_control"]
    live = st.remainder_control()
    for key in ctl:
        for name, v in ctl[key].items():
            assert live[key][name] == pytest.approx(v, rel=1e-12), (key, name)
    assert ctl["lam1"]["cross42_with"] == pytest.approx(155.0285, abs=1e-3)
    assert ctl["lam1"]["cross42_without"] == pytest.approx(150.6473, abs=1e-3)
    assert ctl["lam2"]["mu2_a4_without"] == pytest.approx(6531.893, abs=1e-2)
    for key in ("lam1", "lam2"):
        g = art["gates"][f"cross42_k2_{key}"]
        lim, spread = g["limit"]["central"], g["limit"]["spread"]
        with_, without = ctl[key]["cross42_with"], ctl[key]["cross42_without"]
        assert g["closed_form"] == pytest.approx(with_, rel=1e-12)
        assert abs(lim - with_) < 0.15 * abs(with_ - without), (key, lim)
        assert abs(lim - without) > 0.8 * abs(with_ - without), (key, lim)
        assert spread < 0.5 * abs(with_ - without), (key, spread)


def test_ed_commutant_projection(art):
    """Dense ED, per instance, removing the whole commutant of H: the
    stationary weight is the H^2 projection (exact hypergeometric
    formula), and at p = 6 the dynamical width sits within 2 % of the
    N = infinity 2(1 - q^2) while the raw width is 2-7 times smaller.  The
    smallest row is re-executed."""
    rows = art["ed_commutant"]
    for r in rows:
        assert abs(r["mu0_stationary"] - r["h2_projection_exact"]) < \
            max(4 * r["mu0_stationary_err"], 0.05 * r["h2_projection_exact"]), r
    p6 = [r for r in rows if r["p"] == 6]
    assert len(p6) >= 3
    for r in p6:
        assert abs(r["w2_dynamic"] / r["w2_dynamic_ninf"] - 1) < 0.02, r
        assert r["w2_raw"] < 0.5 * r["w2_dynamic"]
    r0 = rows[0]
    live = st.ed_commutant_split(r0["N"], r0["p"], r0["n_inst"], r0["seed0"])
    for key in ("mu0", "mu0_stationary", "mu0_dynamic", "w2_dynamic"):
        assert live[key] == pytest.approx(r0[key], rel=1e-9), key


def test_fixed_p_dynamical_width(art):
    """Along fixed p the identity-subtracted width has a maximum; with the
    H^2 component removed the width falls monotonically from N = 12 on,
    and in units of J the mean-square frequency rises monotonically to
    8p(p-2).  Two rows are recomputed exactly."""
    for p_str, rows in art["fixed_p_dynamic"].items():
        p = int(p_str)
        rows = sorted(rows, key=lambda r: r["N"])
        raw = [r["w2_raw"] for r in rows]
        dyn = [r["w2_dynamic"] for r in rows]
        phys = [r["omega2_dynamic_physical"] for r in rows]
        ipk = raw.index(max(raw))
        assert 0 < ipk < len(raw) - 1                         # raw: a maximum
        tail = [d for r, d in zip(rows, dyn) if r["N"] >= 16]
        assert all(b < a for a, b in zip(tail, tail[1:]))     # dynamical: none
        assert all(b > a for a, b in zip(phys, phys[1:]))
        assert 0.9 * 8 * p * (p - 2) < phys[-1] < 8 * p * (p - 2)
        # the small-N end sits on the N = infinity law
        assert abs(dyn[2] / rows[2]["w2_dynamic_ninf"] - 1) < 0.06
    live = st.fixed_p_dynamic(4, (16, 32))
    stored = {r["N"]: r for r in art["fixed_p_dynamic"]["4"]}
    for r in live:
        assert r["w2_dynamic"] == pytest.approx(stored[r["N"]]["w2_dynamic"],
                                                rel=1e-13)
