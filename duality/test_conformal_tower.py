"""Tests for conformal_tower.py: the lambda -> 0 dressings of the singlet
tower against the tanh-frame SL(2) OPE and one Liouville rung."""
import json
from ast import literal_eval
from fractions import Fraction
from pathlib import Path

import pytest

import conformal_tower as ct

HERE = Path(__file__).resolve().parent
F = Fraction


@pytest.fixture(scope="module")
def limits():
    """n <= 5 (kmax = 3), both chord rules: a few seconds."""
    return dict(on=ct.chord_limits(3, remainder=True),
                off=ct.chord_limits(3, remainder=False))


@pytest.fixture(scope="module")
def art():
    with open(HERE / "conformal_tower.json") as fh:
        return json.load(fh)


# -- the OPE side, by itself -------------------------------------------------
def test_frame_top_coefficient_is_block_coefficient():
    """The w^(n-m) coefficient of the frame prediction is the plain 1F1
    block coefficient 4 C(n,m) (m)_(n-m)/(2m)_(n-m), for every n <= 8."""
    frame = ct.frame_dressings(8)
    for (n, m), P in frame.items():
        assert P[n - m] == ct.block_coefficient(n, m)
        assert all((n - m - e) % 2 == 0 for e in P)


def test_frame_known_values():
    """Hand-checkable values: A_3 > 6 w O_2 (the exact identity
    A_3 - A_3^dag = 3[H, A_2] at leading order) and the frame corrections
    of A_4 and A_6 in the O_2 sector."""
    frame = ct.frame_dressings(6)
    assert frame[(3, 2)] == {1: 6}
    assert frame[(4, 2)] == {2: F(36, 5), 0: F(16, 5)}
    assert frame[(6, 2)] == {4: F(60, 7), 2: F(96, 7), 0: F(64, 7)}
    assert frame[(5, 3)] == {2: F(80, 7), 0: F(80, 7)}


def test_ladder_moments():
    """<w^k> of sech^(2m): 2m, 12m^2 + 4m."""
    for m in range(1, 7):
        assert ct.ladder_moment(m, 2) == 2 * m
        assert ct.ladder_moment(m, 4) == 12 * m * m + 4 * m
        assert ct.ladder_moment(m, 3) == 0


def test_gr_ratio_matches_largeq_anchor():
    """The exact GR couplings equal the float tower of largeq_anchor (its
    Gamma-function transcription of GR (2.15)) at large q, in absolute
    normalization (q^2 c_n^2, which fixes the h = 3 check) and in ratio."""
    import largeq_anchor as la
    q = 1e8
    rows = la.gr_tower(q, n_max=3)
    for n in (1, 2, 3):
        assert float(ct.gr_c2_q2(n)) == pytest.approx(q * q * rows[n - 1]["c2"], rel=1e-6)
    assert ct.gr_c2_q2(1) == 1
    for n in (2, 3):
        assert float(ct.gr_ratio(n)) == pytest.approx(rows[n - 1]["c2"] / rows[0]["c2"],
                                                       rel=1e-6)
    assert ct.gr_ratio(2) == F(44, 5103)
    assert ct.gr_ratio(3) == F(1, 4000)


# -- the chord engine's limit -----------------------------------------------
def test_no_remainder_rule_is_consistent():
    """The chord rule with the remainder off still decomposes exactly into
    dressed matter chords (zero residual on every checked word, exact
    rationals at two generic points), so its dressings are well defined."""
    import singlet_tower as st
    for lam, q in st.EXACT_POINTS[:2]:
        d = ct.NoRemainder(lam, q).decompose(3, check_len=3)
        assert d["max_residual"] == 0 and d["n_checked"] == 84


def test_limits_are_clean(limits):
    """Rational reconstruction is exact to roundoff of the extrapolation,
    no limit moves when an end point is dropped, and every sigma term and
    wrong-parity w term vanishes in the limit (the unsquared scaled size)."""
    for d in limits.values():
        assert d["max_residual"] < 1e-24
        assert d["max_drift"] < 1e-16
        assert d["sub_max"] < 1e-26


def test_limit_guard_catches_divergence():
    """A stray c/lambda extrapolates to a clean-looking finite number
    (c * sum 1/lambda_i, 390144 c on this grid); the end-point drift
    exposes it."""
    xs = list(ct.LAMS)
    finite = ct._limit(xs, [F(2) + F(3, 7) * x for x in xs])
    assert finite == (F(2), 0.0)
    L, drift = ct._limit(xs, [F(2) + F(3, 7) / x for x in xs])
    assert sum(1 / x for x in xs) == 390144
    assert L == 2 + F(3, 7) * 390144
    assert drift > 0.1 > ct.DRIFT_TOL


def test_norms(limits):
    """nu_m = norm_m / lambda^(2 floor(m/2)) at lambda -> 0: 1, 2, 52/5
    (the closed forms) and 88/7; identical with and without the remainder."""
    assert limits["on"]["norms"] == {2: 1, 3: 2, 4: F(52, 5), 5: F(88, 7)}
    assert limits["off"]["norms"] == limits["on"]["norms"]


def test_no_remainder_is_the_ope(limits):
    """Result 1: with the grade remainder off, every dressing through
    n = 5 is exactly the tanh-frame OPE (both ways: no predicted sector is
    missing and none is extra)."""
    frame = ct.frame_dressings(5)
    assert set(limits["off"]["dressings"]) == set(frame)
    for key, P in frame.items():
        assert limits["off"]["dressings"][key] == P, key


def test_remainder_is_one_rung(limits):
    """Result 2: with it on, the dressings through n = 5 are the OPE plus
    one Liouville rung, -2 a^2 sech^2 (W = 1 suffices at this depth, and
    iterating the rung changes nothing before n - m = 4)."""
    frame = ct.frame_dressings(5)
    rung = ct.ladder_dressings(frame, 5)
    assert set(limits["on"]["dressings"]) == set(rung)
    for key, P in rung.items():
        assert limits["on"]["dressings"][key] == P, key
    assert rung == ct.ladder_dressings(frame, 5, rung=(1,)) \
        == ct.ladder_dressings(frame, 5, iterations=5)
    assert limits["on"]["dressings"][(4, 2)] == {2: F(36, 5), 0: F(56, 5)}


def test_rung_series():
    """sech^2 tau = 1 - tau^2 + 2/3 tau^4 - 17/45 tau^6 + ..., and the
    rung is iterated only on request: the first difference is at
    n - m = 4, in the three coefficients that see tau^2."""
    assert ct.sech2_series(6) == (1, 0, -1, 0, F(2, 3), 0, F(-17, 45))
    frame = ct.frame_dressings(7)
    one, it = ct.ladder_dressings(frame, 7), ct.ladder_dressings(frame, 7, iterations=7)
    assert [k for k in sorted(one) if one[k] != it[k]] == [(6, 2), (7, 2), (7, 3)]


def test_shares(limits):
    """Result 3: exact lambda -> 0 shares, and the engine's own sector
    weights agree with the ones rebuilt from the limits."""
    sh = ct.limiting_shares(limits["on"]["dressings"], limits["on"]["norms"], 5)
    assert sh[3]["shares"] == {3: 1} and sh[3]["fading"] == [2]
    assert sh[4]["shares"] == {3: 1} and sh[4]["fading"] == [2, 4]
    assert sh[5]["shares"] == {3: F(2047, 2058), 5: F(11, 2058)}
    assert sh[5]["fading"] == [2, 4]
    # the engine's sector weights / lambda^(n-1), extrapolated: zero for the
    # fading sectors, the same split for the surviving ones
    mu = {m: v for (n, m), v in limits["on"]["mu0"].items() if n == 5}
    tot = sum(mu.values())
    assert float(mu[5] / tot) == pytest.approx(11 / 2058, rel=1e-10)
    assert abs(float(mu[2] / tot)) < 1e-10 and abs(float(mu[4] / tot)) < 1e-10


def test_survival_follows_the_norm_law(limits):
    """Result 3 from the data: every sector of A_n has a nonzero scaled
    dressing (F_{n,m} ~ a^(n-m)), so with the norm law the weight of O_m
    goes as lambda^(n - m + 2 floor(m/2)), and the engine's own sector
    weights / lambda^(n-1) confirm it: nonzero for odd m, zero for even m,
    at every n <= 5."""
    d = limits["on"]
    for n in range(3, 6):
        for m in range(2, n + 1):
            assert d["dressings"][(n, m)], (n, m)
            mu = float(d["mu0"][(n, m)])
            assert (mu > 1e-3) if m % 2 else (abs(mu) < 1e-9), (n, m, mu)
    sh = ct.limiting_shares(d["dressings"], d["norms"], 5)
    for n, r in sh.items():
        assert r["lead_power"] == n - 1
        assert set(r["shares"]) == {m for m in range(3, n + 1, 2)}
        assert r["fading"] == list(range(2, n + 1, 2))


def test_odd_ope_ratio(limits):
    """C_5^2/C_3^2 = 11/700 at T = infinity, against GR's 44/5103 in the
    IR: different, by (27/20)^2."""
    c2 = ct.ope_coefficients(limits["on"]["norms"])
    assert c2[3] == F(8, 9) and c2[5] == F(22, 1575)
    assert c2[5] / c2[3] == F(11, 700)
    assert c2[5] / c2[3] / ct.gr_ratio(2) == F(27, 20) ** 2


# -- the artifact ------------------------------------------------------------
def _parse(d):
    return {literal_eval(k): {int(e): F(v) for e, v in P.items()} for k, P in d.items()}


def test_artifact_matches_current_code(art):
    """The committed dressings against the OPE and ladder recomputed now
    (so a stale artifact or a changed prediction cannot hide behind the
    stored booleans), and its limit guards."""
    nmax = art["nmax"]
    frame = ct.frame_dressings(nmax)
    assert _parse(art["frame_dressings"]) == frame
    assert _parse(art["chord_off"]) == frame
    rung = ct.ladder_dressings(frame, nmax)
    on = _parse(art["chord_on"])
    assert set(on) == set(rung)
    assert [k for k in sorted(rung) if on[k] != rung[k]] == []
    # with W = 1, or with the rung iterated, the three tau^2-sensitive miss
    for other in (ct.ladder_dressings(frame, nmax, rung=(1,)),
                  ct.ladder_dressings(frame, nmax, iterations=nmax)):
        assert [k for k in sorted(other) if on[k] != other[k]] == [(6, 2), (7, 2), (7, 3)]
    norms = {int(m): F(v) for m, v in art["norms"].items()}
    shares = ct.limiting_shares(on, norms, nmax)
    for n, r in shares.items():
        assert {int(m): F(v) for m, v in art["shares"][str(n)]["shares"].items()} \
            == r["shares"]
    assert art["max_residual"] < 1e-20 and art["max_drift"] < 1e-12
    assert art["sub_max"] < 1e-20


def test_artifact_through_n6(art):
    """The n <= 6 part of the committed run: every check true, the A_6
    shares exact."""
    assert all(art["top_is_block_coefficient"].values())
    assert all(art["frame_matches_no_remainder"].values())
    assert art["norms_off_equal"]
    assert art["norms"]["6"] == "4120/7"
    assert art["shares"]["6"]["shares"] == {"3": "20443/20608", "5": "165/20608"}


def test_artifact_n7(art):
    """The committed run reaches n = 7: OPE exact there too, N_7 -> 720,
    the GR relation at h = 3, 5, 7, one sech^2 rung exact everywhere (the
    three tau^2-sensitive coefficients each need Liouville's -1), and the
    n = 8 and h = 9 predictions."""
    assert art["nmax"] == 7
    assert art["norms"]["7"] == "720"
    assert all(k in art["frame_matches_no_remainder"] for k in ("(7, 2)", "(7, 6)"))
    assert {h: r["equal"] for h, r in art["odd_ope_vs_gr"].items()} == \
        {"3": True, "5": True, "7": True}
    assert art["odd_ratio_7_3"]["chord_over_gr"] == "100/49"
    misses = {(r["n"], r["m"], r["w_power"]): r["w2"] for r in art["rung_tau2_fit"]}
    assert misses == {(6, 2, 0): "-1", (7, 2, 1): "-1", (7, 3, 0): "-1"}
    assert all(art["rung_matches_full"].values())
    assert [k for k, v in art["rung_w1_matches_full"].items() if not v] == \
        ["(6, 2)", "(7, 2)", "(7, 3)"]
    assert art["iterated_rung_misses"] == ["(6, 2)", "(7, 2)", "(7, 3)"]
    nxt = art["prediction_next_level"]
    assert nxt["n"] == 8
    assert nxt["chord_on"]["2"] == {"0": "76352/35", "2": "93488/105",
                                    "4": "9892/105", "6": "28/3"}
    assert nxt["chord_on"]["4"] == {"0": "6064/11", "2": "8744/33", "4": "980/33"}
    assert art["shares"]["7"]["shares"] == {"3": "540648/548261", "5": "7578/548261",
                                             "7": "5/78323"}
    assert F(art["prediction"]["nu"]) == F(21482496, 143)
    # the engine's own extrapolated sector weights agree with the rebuilt shares
    eng = art["shares_engine_mu0"]["7"]
    assert eng["7"] == pytest.approx(5 / 78323, rel=1e-8)


def test_gff_coefficients():
    """The generalized-free-fermion coefficients behind gff_c2: the block
    expansion of F_0 = -chi^(2D) + (chi/(1-chi))^(2D) at D = 0.1, and
    c_h^2/D^2 -> gff_c2(h) as D -> 0; and the GFF form of the relation is
    the GR form exactly."""
    import mpmath as mp
    with mp.workdps(30):
        def c2(D, n):
            h = 2 * D + 2 * n + 1
            return (2 * mp.gamma(h) ** 2 * mp.gamma(h + 2 * D - 1)
                    / (mp.gamma(2 * D) ** 2 * mp.gamma(2 * h - 1)
                       * mp.gamma(h - 2 * D + 1))), h
        D, x = mp.mpf("0.1"), mp.mpf("0.4")
        f0 = -x ** (2 * D) + (x / (1 - x)) ** (2 * D)
        blocks = sum(c * x ** h * mp.hyp2f1(h, h, 2 * h, x)
                     for c, h in (c2(D, n) for n in range(60)))
        assert abs(blocks - f0) < mp.mpf(10) ** -25
        D = mp.mpf(10) ** -12
        for n in (1, 2, 3, 4):
            c, _ = c2(D, n)
            g = ct.gff_c2(2 * n + 1)
            assert abs(c / D ** 2 / (g.numerator / mp.mpf(g.denominator)) - 1) < 1e-10
    for h in (3, 5, 7, 9, 11):
        assert ct.odd_ope_from_gff(h) == ct.odd_ope_from_gr(h)


def test_gr_relation_closed_form():
    """C_h^2 = N (1 - k_c)^2 c_h^2(GR) in closed form: 2 (s+1)/(s(s-1)) *
    4 (2n)!^2/(4n)!, s = h(h-1)/2, h = 2n + 1."""
    from math import factorial
    for n in range(1, 6):
        h, s = 2 * n + 1, n * (2 * n + 1)
        closed = 2 * F(s + 1, s * (s - 1)) * F(4 * factorial(2 * n) ** 2, factorial(4 * n))
        assert ct.odd_ope_from_gr(h) == closed


@pytest.mark.slow
def test_n6_live():
    """Recompute n <= 6 (about a minute)."""
    on = ct.chord_limits(4, remainder=True)
    off = ct.chord_limits(4, remainder=False)
    frame = ct.frame_dressings(6)
    assert off["dressings"] == frame
    assert on["norms"][6] == F(4120, 7)
    rung = ct.ladder_dressings(frame, 6)
    assert all(on["dressings"][k] == P for k, P in rung.items())
    w1 = ct.ladder_dressings(frame, 6, rung=(1,))
    assert [k for k, P in w1.items() if on["dressings"][k] != P] == [(6, 2)]
