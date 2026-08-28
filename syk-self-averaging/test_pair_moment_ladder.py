"""Tests for the mu_6 / K = 12 ladder campaign (pair_moment_ladder.py).

Pinned here: the deficit rule deficit = K/2 + 1 (K = k + 2n) at every
computed law, committed and new; the exact n = 1 -> 2 coefficient
identity A(n=1, k) = A(n=2, k-2) at every rung and every shared p, and
the MEASURED failure of every deeper step; the K = 12 frontier values;
the exact k = 10 flavor-sum gate; the closed forms (including the
A(n=3) form the committed campaign reported as unidentified); and the
double-scaling extension with its cross-artifact overlap.

The cheapest new law is re-executed live (producer re-execution, not
just JSON reading); the K = 12 gate is re-executed in the slow tier.
"""

import json
from fractions import Fraction
from pathlib import Path

import pytest

from pair_moment_laws import _poly_eval, law
from pair_moment_ladder import a1_identity_k10, rung_table

HERE = Path(__file__).resolve().parent

with open(HERE / "pair_moment_ladder.json") as fh:
    DATA = json.load(fh)
with open(HERE / "pair_moment_laws.json") as fh:
    COMMITTED = json.load(fh)
LAWS = {(r["n"], r["k"], r["p"]): r for r in DATA["laws"]}
OLD = {(r["n"], r["k"], r["p"]): r for r in COMMITTED["laws"]}
ALL = {**OLD, **LAWS}


def test_producer_reexecution_cheapest_new_law():
    """law(p=2, n=1, k=6): the first mu_6 law, recomputed end to end
    (sampling, interpolation, series, out-of-sample certification)."""
    r = law(2, 1, 6)
    assert (r["deficit"], r["A"], r["B"]) == (5, "-57344", "-28672")
    assert LAWS[(1, 6, 2)]["A"] == r["A"]
    assert LAWS[(1, 6, 2)]["P_coeffs"] == r["P_coeffs"]
    # and it coincides with the committed mu_4(n=2) law: the ladder
    assert OLD[(2, 4, 2)]["A"] == r["A"]


def test_k12_gate_recorded_exact_zero():
    """The K = 12 machinery is gated by the exact identity
    (N-1) mu_10^pair(n=1) + mu_10^diag(n=1) = 0, in Fractions."""
    g = DATA["a1_gate_k10"]
    assert g["exact_zero"] is True
    assert Fraction(g["residual"]) == 0
    assert (g["N"], g["p"], g["k"]) == (8, 2, 10)
    assert (7 * Fraction(g["pair"]) + Fraction(g["diag"])) == 0


@pytest.mark.slow
def test_k12_gate_reexecution():
    """Re-run the k = 10 flavor-sum gate (both channels, K = 12 words)."""
    g = a1_identity_k10(N=8, p=2)
    assert g["exact_zero"] is True
    assert g["pair"] == DATA["a1_gate_k10"]["pair"]


def test_deficit_rule_K_half_plus_1_everywhere():
    """deficit = K/2 + 1 with K = k + 2n, at every law in both artifacts:
    n + 2 for k = 2, n + 3 for k = 4, n + 4 for k = 6, ... through the
    K = 12 frontier where every deficit is 7."""
    assert len(ALL) >= 45
    for (n, k, p), r in ALL.items():
        K = k + 2 * n
        assert r["deficit"] == K // 2 + 1, (n, k, p)
    frontier = [(n, k) for (n, k, p) in LAWS if k + 2 * n == 12]
    assert set(frontier) >= {(1, 10), (2, 8), (3, 6), (4, 4), (5, 2)}
    assert DATA["ladder"]["deficit_rule_K_half_plus_1"] is True


def test_mu6_laws_present():
    """mu_6 (k = 6): n = 1 over p = 2..12, n = 2 at p = 2, 4, n = 3 at
    p = 2 -- the extension road_map.md Part III named as unlocked."""
    assert {p for (n, k, p) in LAWS if (n, k) == (1, 6)} >= {2, 4, 6, 8, 10, 12}
    assert {p for (n, k, p) in LAWS if (n, k) == (2, 6)} >= {2, 4}
    assert (3, 6, 2) in LAWS


def test_n1_to_n2_step_exact_at_every_rung():
    """A(n=1, k) = A(n=2, k-2) exactly, at every shared p, for
    k = 4 (committed), 6, 8, 10 -- rungs K = 6, 8, 10, 12."""
    pairs = 0
    for (n, k, p), r in ALL.items():
        if n == 1 and k >= 4 and (2, k - 2, p) in ALL:
            assert Fraction(r["A"]) == Fraction(ALL[(2, k - 2, p)]["A"]), (k, p)
            assert r["deficit"] == ALL[(2, k - 2, p)]["deficit"]
            pairs += 1
    assert pairs >= 12          # k=4: p=2,4,6; k=6: p=2..12(+14); k=8,10
    # the rung table records the same verdict
    for chk in DATA["ladder"]["adjacent_step_checks"]:
        if chk["lhs"].startswith("A(n=1,"):
            assert chk["equal"] is True, chk


def test_deeper_steps_are_broken():
    """The measured inequalities: every step with n >= 2 on the left
    differs, at every rung and every shared p.  Stated as measured
    exact facts of the coefficients."""
    broken = 0
    for chk in DATA["ladder"]["adjacent_step_checks"]:
        if not chk["lhs"].startswith("A(n=1,"):
            assert chk["equal"] is False, chk
            assert Fraction(chk["lhs_A"]) != Fraction(chk["rhs_A"])
            broken += 1
    assert broken >= 6
    # the K = 10 rung at p = 2 in one line: 21 : 21 : 14 : 9 times -2^17
    unit = -(1 << 17)
    assert [Fraction(ALL[(n, 10 - 2 * n, 2)]["A"]) / unit
            for n in (1, 2, 3, 4)] == [21, 21, 14, 9]


def test_k12_frontier_values():
    """The K = 12 frontier at p = 2, exact leading coefficients."""
    want = {(1, 10): "-155713536", (2, 8): "-155713536",
            (3, 6): "-97255424", (4, 4): "-54001664", (5, 2): "-32505856"}
    for (n, k), A in want.items():
        assert LAWS[(n, k, 2)]["A"] == A, (n, k)
        assert LAWS[(n, k, 2)]["deficit"] == 7
    # the k = 2 deficit law mu_2 ~ N^-(n+2) now certified through n = 5
    assert LAWS[(5, 2, 2)]["deficit"] == 5 + 2


def test_subleading_shared_at_p2_only():
    """Measured, not generalized: the subleading B coefficient is also
    shared across the n = 1 -> 2 step at p = 2 (K = 10 and 12), but not
    at p = 4."""
    assert LAWS[(1, 8, 2)]["B"] == LAWS[(2, 6, 2)]["B"] == "-2064384"
    assert LAWS[(1, 10, 2)]["B"] == LAWS[(2, 8, 2)]["B"] == "-173015040"
    assert LAWS[(1, 8, 4)]["B"] != LAWS[(2, 6, 4)]["B"]


def test_rung_table_helper_is_consistent():
    """rung_table over both artifacts reproduces the stored verdicts."""
    t = rung_table(COMMITTED["laws"] + DATA["laws"])
    assert t["deficit_rule_K_half_plus_1"] is True
    stored = {(c["lhs"], c["rhs"], c["p"]): c["equal"]
              for c in DATA["ladder"]["adjacent_step_checks"]}
    fresh = {(c["lhs"], c["rhs"], c["p"]): c["equal"]
             for c in t["adjacent_step_checks"]}
    assert fresh == stored


def test_p14_and_p16_n3_laws_match_committed_cross_check():
    """The committed campaign's failed degree-5 fit recorded the direct
    p = 14 value of A(n=3, k=2); the pooled path reproduces it."""
    assert LAWS[(3, 2, 14)]["A"] == "-670309412864"
    assert LAWS[(3, 2, 16)]["deficit"] == 5


def test_out_of_sample_certification_recorded():
    for key, r in LAWS.items():
        assert r["out_of_sample_N"] > 0, key
        assert len(r["P_coeffs"]) - 1 == r["deg_P"], key


def _form(n, k):
    return next(c for c in DATA["A_closed_forms"]
                if c["n"] == n and c["k"] == k)


def test_A_closed_form_n3_identified():
    """The closed form the committed campaign reported as unidentified
    (degree > 5 at p = 14): with the p grid deepened to 16 it is the
    degree-7 polynomial

        A(n=3, k=2; p) = -256 p^4 (p-1) (34 p^2 - 108 p + 91),

    certified out of sample at p = 18 and re-asserted here at every
    sampled p in both artifacts."""
    cf = _form(3, 2)
    assert cf["identified"] and cf["degree"] == 7 and cf["p_check"] == 18
    coeffs = [Fraction(c) for c in cf["A_poly_coeffs"]]
    for (n, k, p), r in ALL.items():
        if (n, k) == (3, 2):
            want = -256 * p ** 4 * (p - 1) * (34 * p * p - 108 * p + 91)
            assert Fraction(r["A"]) == want, p
            assert _poly_eval(coeffs, Fraction(p)) == want, p
    # the k = 2 family: A(n) = -2^(2n+2) p^(n+1) (p-1) Q_n(p), deg Q_n = n-1
    for (n, k, p), r in ALL.items():
        if k == 2 and n <= 3:
            q = Fraction(r["A"]) / (-(4 ** (n + 1)) * p ** (n + 1) * (p - 1))
            assert q.denominator == 1, (n, p)


def test_double_scaling_overlap_and_extension():
    """Cross-artifact gate: the exact c_2(p) = |mu_2| p N on
    lambda_chord = 2 reproduces the committed rows at p = 14, 16, then
    extends them; the sequence is monotone rising through the deepest
    point, and the n = 3 line continues past its committed p = 12."""
    rows = {r["p"]: r for r in DATA["double_scaling"]["n2"]}
    old = {r["p"]: r for r in COMMITTED["double_scaling"]["n2"]}
    for p in (14, 16):
        assert rows[p]["mu2"] == pytest.approx(old[p]["mu2"], rel=1e-12)
        assert rows[p]["N"] == p * p
    ps = sorted(rows)
    assert ps[-1] >= 28
    cs = [rows[p]["c"] for p in ps]
    assert all(b > a for a, b in zip(cs, cs[1:]))
    assert rows[28]["c"] == pytest.approx(31.218, abs=5e-4)
    n3 = {r["p"]: r["c"] for r in DATA["double_scaling"]["n3"]}
    assert n3[14] == pytest.approx(208.333, abs=5e-4)
    assert n3[16] == pytest.approx(209.759, abs=5e-4)


def test_c2_limit_is_an_estimate_with_a_spread():
    """The limit extrapolation is labeled an estimate, carries its
    window spread, and its windows agree with each other at the 1e-2
    level; the exact points stay in the artifact for anyone to refit."""
    lim = DATA["double_scaling"]["c2_limit"]
    assert lim["estimated"] is True
    assert "estimate" in lim["note"]
    vals = [f["c_inf"] for f in lim["fits"]]
    assert max(vals) - min(vals) == pytest.approx(lim["spread"])
    assert lim["spread"] < 2e-2
    assert 32.0 < lim["c_inf_central"] < 32.6
    # and it lies above every exact point, since the sequence rises
    assert lim["c_inf_central"] > max(r["c"] for r in DATA["double_scaling"]["n2"])


def test_K8_family_closed_form():
    """The other K = 8 family, A(1,6) = A(2,4): the campaign's six-sample
    fit failed at degree 5 (recorded for A(2,4), whose entry is left as
    the campaign wrote it); the --refit over eight samples is degree 7
    and certifies out of sample at p = 18:

        A(1,6; p) = A(2,4; p) = -256 p^4 (p-1) (34 p^2 - 93 p + 64),

    the same -2^8 p^4 (p-1) prefactor and the same leading 34 p^2 as
    A(3,2)'s quadratic.  Re-asserted at every sample of BOTH families."""
    cf = _form(1, 6)
    assert cf["identified"] and cf["degree"] == 7 and cf["p_check"] == 18
    coeffs = [Fraction(c) for c in cf["A_poly_coeffs"]]
    checked = 0
    for (n, k, p), r in ALL.items():
        if (n, k) in ((1, 6), (2, 4)):
            want = -256 * p ** 4 * (p - 1) * (34 * p * p - 93 * p + 64)
            assert Fraction(r["A"]) == want, (n, k, p)
            assert _poly_eval(coeffs, Fraction(p)) == want, (n, k, p)
            checked += 1
    assert checked >= 15
    old = _form(2, 4)
    assert old["p_check"] == 14 and not old["identified"] and old["degree"] == 5
    assert LAWS[(1, 6, 14)]["A"] == LAWS[(2, 4, 14)]["A"]
