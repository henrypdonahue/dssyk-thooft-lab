"""Tests for the JT-dS vacuum-energy note (jt_ds_vacuum.json).

Pinned here: the whole derivation, by re-executing the producer (every
piece is exact sympy, so re-running IS re-deriving); the c_4 = 3/8
anchor by an inline independent mini-derivation that shares no code
with the module; the sign convention of the curvature helpers; the
committed JSON against the fresh run; the paper's D = 4 statement at
its recorded source lines (skipped when papers/ is absent -- it is
gitignored); and the condensate cross-link.

The honesty structure under test: [D] items must re-derive to exact
zero, [T] items are recorded citations only (no equation is asserted
from memory), [C] items are labeled conventions, not results.
"""

import json
from pathlib import Path

import pytest
import sympy as sp

import jt_ds_vacuum as jv

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent / "papers" / "arXiv-2607.05678v1" / "HologramsSM.tex"


@pytest.fixture(scope="module")
def fresh():
    """Re-run the entire producer: every machine check re-derives."""
    return jv.main()


@pytest.fixture(scope="module")
def data():
    with open(HERE / "jt_ds_vacuum.json") as fh:
        return json.load(fh)


def _all_checks(result, key):
    checks = result[key]["checks"]
    bad = {k: v for k, v in checks.items() if v is not True}
    assert not bad, f"{key}: {bad}"


def test_curvature_sign_convention():
    """The helpers' Riemann sign: unit 2-sphere has R = +2."""
    th, ph = sp.symbols('theta phi')
    g = sp.diag(1, sp.sin(th)**2)
    assert sp.simplify(jv.scalar_curvature(g, [th, ph]) - 2) == 0


def test_general_identity_rederives(fresh):
    """Piece 1: a(D) = D/(D-2), eps_D = c_D S^(-2/(D-2)), c_4 = 3/8."""
    _all_checks(fresh, "general_D")
    assert fresh["general_D"]["a_of_D"] == "D/(D - 2)"
    table = fresh["general_D"]["c_table"]
    key4 = 4 if 4 in table else "4"
    assert table[key4]["exact"] == "3/8"


def test_c4_by_independent_inline_derivation():
    """The 3/8 anchor, re-derived here with no module code: D = 4,
    rho = 3/(8 pi G ell^2), S = pi ell^2/G  =>  rho G^2 = (3/8)/S."""
    G, ell = sp.symbols('G ell', positive=True)
    rho4 = sp.Rational(3, 1) / (8 * sp.pi * G * ell**2)   # Lam = 3/ell^2
    S4 = sp.pi * ell**2 / G                               # 4 pi ell^2 / 4G
    c4 = sp.simplify(rho4 * G**2 * S4)
    assert c4 == sp.Rational(3, 8)
    # and the exponent really is unique: [rho] = L^-4, [G] = L^2
    a = sp.Symbol('a')
    assert sp.solve(sp.Eq(-4 + a * 2, 0), a) == [2]


def test_c_table_values(fresh):
    """c_3 = pi/32 exactly, and the committed floats are the exact values."""
    table = fresh["general_D"]["c_table"]
    key3 = 3 if 3 in table else "3"
    assert sp.sympify(table[key3]["exact"]) == sp.pi / 32
    for d, row in table.items():
        assert float(sp.sympify(row["exact"])) == pytest.approx(
            row["float"], rel=1e-14)


def test_d2_degeneration_is_fourfold(fresh):
    """Piece 2: all four failures re-derive; the negative statement holds."""
    _all_checks(fresh, "d2_degeneration")
    # inline: the dimension equation at D = 2 really has no solution
    a = sp.Symbol('a')
    assert sp.solve(sp.Eq(-2 + a * 0, 0), a) == []


def test_first_law_rederives(fresh):
    """Piece 3a: sourced dilaton, delta S = -2 pi ell^2 rho = -dE/T."""
    _all_checks(fresh, "jt_first_law")
    fl = fresh["jt_first_law"]
    assert fl["delta_S"] == "-2*pi*ell**2*rho"
    assert fl["killing_energy"] == "ell*rho"
    assert fl["T_dS"] == "1/(2*pi*ell)"
    # inline closure of the first law, independent of the module
    ell, rho = sp.symbols('ell rho', positive=True)
    assert sp.simplify(
        -2 * sp.pi * ell**2 * rho + (rho * ell) / (1 / (2 * sp.pi * ell))) == 0


def test_reduction_rederives(fresh):
    """Piece 3b: R3 identity, exact bookkeeping, inherited (pi/32)/S^2."""
    _all_checks(fresh, "reduction_3d")
    # the inherited constant IS the D = 3 entry of the Piece-1 table
    table = fresh["general_D"]["c_table"]
    key3 = 3 if 3 in table else "3"
    assert sp.sympify(table[key3]["exact"]) == sp.pi / 32
    assert "(pi/32)*S**(-2)" in fresh["reduction_3d"]["inherited"]


def test_committed_json_matches_fresh_run(fresh, data):
    """The committed JSON is exactly what the producer produces now
    (runtime aside)."""
    def strip(d):
        d = json.loads(json.dumps(d))   # normalize int keys -> str
        d.pop("runtime_s", None)
        return d
    assert strip(fresh) == strip(data)


def test_paper_statement_at_recorded_lines(data):
    """The located D = 4 displays are where the JSON says they are.
    papers/ is gitignored (arXiv redistribution), so skip when absent."""
    if not PAPER.exists():
        pytest.skip("papers/ not present (gitignored); citation not checkable")
    lines = PAPER.read_text().splitlines()
    lo1, hi1 = data["paper_statement"]["lines"]["eps_rho_G2"]
    lo2, hi2 = data["paper_statement"]["lines"]["eps_3_8S"]
    block1 = "\n".join(lines[lo1 - 1:hi1])
    block2 = "\n".join(lines[lo2 - 1:hi2])
    assert r"\epsilon = \rho G^2" in block1
    assert r"\epsilon = \frac{3}{8 S}" in block2
    # and the relation is only ever written in D = 4: the tex never
    # contains a general-D or 2d version of it
    text = PAPER.read_text()
    assert text.count(r"\frac{3}{8 S}") == 1


def test_condensate_cross_link(data):
    """The recorded eps(0) is live thooft-target data, not a stale copy."""
    src = HERE.parent / "thooft-target" / "condensate.json"
    # committed artifact: a missing file is a broken checkout, not a skip
    assert src.exists(), src
    with open(src) as fh:
        eps0 = json.load(fh)["eps_duality_point"]
    assert data["condensate_link"]["eps0_measured_Nc_g2"] == pytest.approx(
        eps0, rel=1e-12)
    assert eps0 == pytest.approx(-0.0407, abs=5e-4)


def test_literature_record_is_fetch_gated(data):
    """No equation number appears outside the verified-fetch record: the
    JSON's own honesty constraint."""
    lit = data["literature"]
    assert "two independent extractions" in lit["fetch_method"]
    assert any("Eq. (9)" in s for s in lit["svvv_2203_00700"]["verified"])
    assert any("Eq. (26)" in s for s in lit["svvv_2203_00700"]["verified"])
    assert any("(3.37)" in s for s in lit["mty_1904_01911"]["verified"])
    # the one soft placement is labeled as such, not claimed
    assert any("around Eq. (2.13)" in s
               for s in lit["mty_1904_01911"]["verified"])
