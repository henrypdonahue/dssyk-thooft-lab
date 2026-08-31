#!/usr/bin/env python3
r"""
The JT-de Sitter analog of epsilon = 3/(8S): derived where derivable,
and non-unique where it is not.
=====================================================================

THE PAPER'S STATEMENT, LOCATED.  Miyashita-Sekino-Susskind
arXiv:2607.05678, section "Fine Tuning?", writes the vacuum-energy /
entropy relation only in D = 4 (HologramsSM.tex, lines 1014-1016 and
1018-1020; both are unnumbered \be...\ee displays):

    \epsilon = \rho G^2          and          \epsilon = \frac{3}{8 S} .

papers/ is gitignored (arXiv redistribution), so the source is
local-only; the test skips the line check when it is absent.  The
paper's holographic-tuning argument, though, lives in JT-de Sitter.
road_map.md Part III: "Deriving the JT-dS analog is the still-open
theory half of the vacuum-energy story."  This module closes what is
closable and states precisely what is not.  Every equation below is
re-derived by sympy when the test runs.  Provenance labels, used
throughout:

    [D]  derived here and machine-checked (sympy, exact),
    [T]  transcribed from fetched literature and cross-verified,
    [C]  conditional on a stated convention choice.

PIECE 1 -- THE GENERAL-D IDENTITY [D]
-------------------------------------
Conventions: signature (-,+,...,+), Riemann sign such that spheres and
dS have R > 0, Einstein equation G_mn + Lam g_mn = 8 pi G T_mn, vacuum
stress T_mn = -rho g_mn (so a positive rho is a positive energy
density, and moving Lam to the source side gives rho = Lam/(8 pi G)).
dS_D fixes Lam = (D-1)(D-2)/(2 ell^2) -- verified componentwise from
the explicit static-patch metric at D = 3 and D = 4, and by the
maximally-symmetric algebra for symbolic D.  Static-patch horizon
entropy S = Omega_{D-2} ell^{D-2}/(4G), the sphere r = ell.

Dimensional analysis ([rho] = L^-D, [G] = L^{D-2} from S dimensionless)
fixes the UNIQUE exponent making rho G^a dimensionless:

    a(D) = D/(D-2).

Eliminating ell between rho(ell) and S(ell):

    eps_D := rho G^{D/(D-2)} = c_D S^{-2/(D-2)},
    c_D = ((D-1)(D-2)/(16 pi)) (Omega_{D-2}/4)^{2/(D-2)}.

c_4 = 3/8 EXACTLY and a(4) = 2: the paper's pair of displays is the
D = 4 instance, reproduced with no discrepancy.  c_3 = pi/32.  Note
that the 1/S power itself is dimension-specific: the exponent is
-2/(D-2), and equals -1 only in D = 4.

PIECE 2 -- THE D -> 2 CONTINUATION DEGENERATES, PRECISELY [D]
-------------------------------------------------------------
This is a negative statement, and it is fourfold.  Each failure is a
separate machine check.

  (i)  No exponent exists.  The dimension equation -D + a(D-2) = 0 has
       no solution at D = 2 (it reads -2 = 0): continuing the D >= 3
       definition of G through S = Area/4G makes G dimensionless while
       [rho] = L^-2, so NO power rho G^a is dimensionless.  Formally
       a(D) = D/(D-2) poles at D = 2.
  (ii) The Einstein equation has no 2d ancestor.  G_mn == 0
       IDENTICALLY for an arbitrary 2d metric (checked with all three
       components free functions of both coordinates), so the relation
       rho = Lam/(8 pi G) has no Einstein-equation derivation, and the
       vacuum equation would force Lam = 0.  Relatedly the 2d EH
       action is topological: int sqrt(g) R over the round sphere is
       8 pi at every radius.
  (iii) The continued Lam dies while JT's does not.  Lam(D, ell) -> 0
       as D -> 2 at fixed ell, and c_D -> 0.  JT's 2/ell^2 is finite
       because it enters through the DILATON POTENTIAL, not through
       the continued Einstein Lam.
  (iv) The elimination step dies.  At D = 2 the "horizon sphere" is
       Omega_0 = 2 points with ell-independent measure, so
       S = Omega_0 ell^0/(4G) cannot be solved for ell, and in real 2d
       dilaton gravity S = phi_0/4G + phi_h/4G is set by the dilaton,
       not by ell.

So anything called "the JT-dS analog" must come from new input --
first-law bookkeeping or a reduction from a parent -- not from
continuation.  That is what Piece 3 supplies.

PIECE 3 -- WHAT JT-dS ACTUALLY SAYS
-----------------------------------
Stated convention [C], chosen to match Svesko-Verheijden-Verlinde-
Visser arXiv:2203.00700 ("Quasi-local energy and microcanonical
entropy in two-dimensional nearly de Sitter gravity"), their Eq. (26):

    I = (1/16 pi G2) int d^2x sqrt(-g) [ (phi0 + phi) R - 2 phi/ell^2 ]
        + boundary,

dS2 static patch ds^2 = -(1 - r^2/ell^2) dt^2 + dr^2/(1 - r^2/ell^2),
vacuum dilaton phi = phi_r r/ell (SVVV, text after their Eq. (15)),
horizon entropy S = phi0/4G2 + phi_h/4G2 (SVVV Eq. (29)), temperature
T_dS = 1/(2 pi ell) (SVVV Eq. (16); machine-checked here as
kappa/2 pi from the surface gravity).

(a) The intrinsic first law [D].  Add matter vacuum energy
T_mn = -rho g_mn.  The geometry cannot respond: the phi variation
forces R = 2/ell^2 regardless of matter, so ell is RIGID.  What
responds is the dilaton.  At fixed boundary data phi_r (the
coefficient of the unique static homogeneous mode r/ell; a constant is
NOT a homogeneous solution, which is machine-checked), the static
solution shifts by a constant:

    phi = phi_r r/ell - 8 pi G2 rho ell^2 ,

verified componentwise against the metric equation of motion
-grad_m grad_n phi + g_mn Box phi + (phi/ell^2) g_mn = 8 pi G2 T_mn.
(That EOM is the standard JT variation [T]: in 2d the phi G_mn piece
vanishes identically by (ii) above, which is why only Hessian terms
survive.)  Hence, with E = rho ell the Killing energy of the patch
(E = int_Sigma T_mn n^m xi^n dV; the redshift factors cancel exactly,
machine-checked):

    delta S = delta phi_h / 4 G2 = -2 pi ell^2 rho = -delta E / T_dS .

Sign and content agree with the Gibbons-Hawking patch first law as
stated in SVVV Eq. (9), -delta H_xi = T_GH delta S_GH [T].  The
structural contrast with D >= 3 is the actual content: there, rho and
S are LOCKED algebraically (eps_D = c_D S^{-2/(D-2)}); in JT-dS, rho
moves S only differentially, through the dilaton, at fixed ell, and S
carries free additive data (phi0, phi_r) that no vacuum-energy
relation fixes.  In the unique intrinsic dimensionless pairing
eps_JT := rho ell^2 (ell, not G: no power of G2 can absorb dimensions
in 2d),

    delta S = -2 pi eps_JT .

(b) The reduction-inherited relation [D for the D = 3 parent; T for
D = 4].  Reduce the dS3 static patch on its horizon circle,
ds^2_3 = g_ab dx^a dx^b + Phi^2 dtheta^2.  The exact identity
R_3 = R_2 - 2 Box_2 Phi / Phi is machine-checked for an arbitrary
diagonal base with full (t, r) dependence -- and since any 2d metric
is locally diagonalizable while both sides are base scalars, that
loses no generality.  Discarding the total derivative, the 3d EH
action lands exactly on the action above with

    phi = Phi (= r on the solution),  G2 = G3/(2 pi),
    ell_2 = ell_3,  phi0 = 0,          S_JT = S_3   (exactly).

The parent relation of Piece 1 then pulls back to a genuine JT-dS
statement:

    eps_3 = rho_3 G3^3 = (pi/32) S^{-2},

where rho_3 is the PARENT vacuum energy -- in JT variables, the
dilaton-potential coefficient in disguise, NOT 2d matter.  From a
D = 4 parent, the JT sector instead arises only in the near-Nariai
limit of Schwarzschild-dS4 (MTY arXiv:1904.01911 Eqs. (3.37), (3.39)
[T]) and inherits the D = 4 exponent -1 with the parent's S.

THE VERDICT, SHARPLY.  "The JT-dS analog of eps = 3/(8S)" is not one
equation.  It is (a) the exact intrinsic first law
delta S = -2 pi rho ell^2, PLUS (b) the reduction-inherited family
eps_D = c_D S^{-2/(D-2)}, whose exponent depends on the parent
dimension -- S^-2 from dS3, S^-1 from near-Nariai dS4 -- and the
intrinsic JT data (phi0, phi_r, ell) cannot choose between parents.
That choice is a dictionary input, not a gravity output.  Alternative
packagings [C] relabel the same ambiguity: SVVV take phi and G2
dimensionless while the reduction hands us both with length dimension;
a field-dependent G_eff = G2/phi_h (so S = 1/4G_eff mimics area/4G)
changes which symbol carries the dimension, never the conclusion that
"eps" in 2d needs ell.

WHAT THIS SAYS ABOUT THE MEASURED CONDENSATE
--------------------------------------------
thooft-target/condensate.py measures the interaction vacuum-energy
shift of the 2d bulk matter at the duality point:
eps(0) = -0.0407 Nc g^2, relative to the free chiral point.  If that
matter is minimally coupled to JT-dS in the convention above -- which
is the paper's own assumption 1, not a result -- the first law
converts it to a FINITE horizon-entropy shift
delta S = -2 pi ell^2 rho = +2 pi (0.0407) Nc (g ell)^2, positive
because the measured shift is negative.  The number is conditional on
(g ell)^2, a dimensionless knob the dictionary does not currently fix
(under the paper's flat-space limit g^2 = lambda J^2 -> 0 at fixed p
it vanishes per mode), and on the subtraction point E_vac(-1).  What
is NOT conditional: the shift is first-order, finite, and
divergence-free, so on the JT-dS side, exactly as on the flat-space
side, there is nothing for Holographic Tuning to cancel.  The 4d
drama has no 2d counterpart -- now on both sides of the first law.

LITERATURE RECORD (what was fetched against what is merely cited).
Both papers were fetched 2026-08-27 through ar5iv HTML with two
independent extractions each, agreeing; equation numbers below are
cited only where both extractions returned the equation verbatim.
SVVV 2203.00700: Eqs. (9), (16), (26), (29), and the statement after
Eq. (15) that Phi = Phi_r r/L solves the EOM.  MTY 1904.01911:
Eq. (2.1) (action, ell = 1 units, no 1/4G -- their phi absorbs it),
Eqs. (3.37) and (3.39) (SdS4 near-Nariai reduction), and
S_{0,dS} = 8 pi phi_0 in the discussion around their Eq. (2.13) (the
two extractions placed it at/just after (2.13); cited as "around").
Everything else here is derived, not transcribed.

Output: jt_ds_vacuum.json.  Runtime: seconds.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# curvature helpers (Riemann sign: unit 2-sphere has R = +2, checked in test)
# ---------------------------------------------------------------------------
def christoffel(g, coords):
    """Gamma^a_{bc} for metric matrix g in the given coordinates."""
    n = len(coords)
    ginv = g.inv()
    Gam = [[[sp.S(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                expr = sum(ginv[a, d] * (sp.diff(g[d, b], coords[c])
                                         + sp.diff(g[d, c], coords[b])
                                         - sp.diff(g[b, c], coords[d]))
                           for d in range(n)) / 2
                Gam[a][b][c] = sp.together(expr)
    return Gam


def ricci(g, coords):
    """R_{bc} = R^a_{b a c} = d_a Gam^a_{cb} - d_c Gam^a_{ab} + GG - GG."""
    n = len(coords)
    Gam = christoffel(g, coords)
    Ric = sp.zeros(n)
    for b in range(n):
        for c in range(n):
            expr = sp.S(0)
            for a in range(n):
                expr += (sp.diff(Gam[a][c][b], coords[a])
                         - sp.diff(Gam[a][a][b], coords[c]))
                for d in range(n):
                    expr += (Gam[a][a][d] * Gam[d][c][b]
                             - Gam[a][c][d] * Gam[d][a][b])
            Ric[b, c] = sp.cancel(sp.together(expr))
    return Ric


def scalar_curvature(g, coords):
    Ric = ricci(g, coords)
    ginv = g.inv()
    n = len(coords)
    return sp.cancel(sum(ginv[b, c] * Ric[b, c]
                         for b in range(n) for c in range(n)))


def hessian_and_box(g, coords, phi):
    """(grad_b grad_c phi, Box phi) for a scalar phi."""
    n = len(coords)
    Gam = christoffel(g, coords)
    ginv = g.inv()
    H = sp.zeros(n)
    for b in range(n):
        for c in range(n):
            H[b, c] = sp.diff(phi, coords[b], coords[c]) - sum(
                Gam[d][b][c] * sp.diff(phi, coords[d]) for d in range(n))
    box = sum(ginv[b, c] * H[b, c] for b in range(n) for c in range(n))
    return H, box


def omega_sphere(d: int):
    """Omega_{d-2}: volume of the unit (d-2)-sphere, exact."""
    return 2 * sp.pi**(sp.Rational(d - 1, 2)) / sp.gamma(sp.Rational(d - 1, 2))


# ---------------------------------------------------------------------------
# Piece 1: the general-D identity  [D]
# ---------------------------------------------------------------------------
def piece1_general_identity() -> dict:
    """eps_D := rho G^{D/(D-2)} = c_D S^{-2/(D-2)}, c_4 = 3/8 exactly."""
    D = sp.Symbol('D', positive=True)          # D > 2 throughout this block
    G, S, ell = sp.symbols('G S ell', positive=True)
    Om = sp.Symbol('Omega', positive=True)     # Omega_{D-2}, symbolic first
    a = sp.Symbol('a')

    # dimensional analysis: [rho] = L^-D, [G] = L^{D-2}; unique exponent
    a_D = sp.solve(sp.Eq(-D + a * (D - 2), 0), a)[0]
    check_aD = bool(sp.simplify(a_D - D / (D - 2)) == 0)
    check_a4 = a_D.subs(D, 4) == 2             # the paper's eps = rho G^2

    # maximally symmetric input, symbolic D: R_mn = k g_mn with
    # k = (D-1)/ell^2 (spot-verified componentwise at D = 3, 4 below);
    # Einstein G_mn + Lam g_mn = 0 then fixes Lam algebraically.
    k = (D - 1) / ell**2
    Lam_sol = sp.solve(sp.Eq(k - sp.Rational(1, 2) * D * k + sp.Symbol('L'), 0),
                       sp.Symbol('L'))[0]
    Lam = (D - 1) * (D - 2) / (2 * ell**2)
    check_lambda = bool(sp.simplify(Lam_sol - Lam) == 0)

    # explicit-metric verification of the maximally symmetric input, D = 3, 4
    t, r, th, ph = sp.symbols('t r theta phi_c')
    f = 1 - r**2 / ell**2
    g4 = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th)**2)
    Ric4 = ricci(g4, [t, r, th, ph])
    check_ric4 = sp.simplify(Ric4 - (3 / ell**2) * g4) == sp.zeros(4)
    R4 = sp.simplify(scalar_curvature(g4, [t, r, th, ph]))
    check_einstein4 = sp.simplify(
        Ric4 - sp.Rational(1, 2) * R4 * g4 + (3 / ell**2) * g4) == sp.zeros(4)
    g3 = sp.diag(-f, 1 / f, r**2)
    Ric3 = ricci(g3, [t, r, th])
    check_ric3 = sp.simplify(Ric3 - (2 / ell**2) * g3) == sp.zeros(3)
    R3 = sp.simplify(scalar_curvature(g3, [t, r, th]))
    check_einstein3 = sp.simplify(
        Ric3 - sp.Rational(1, 2) * R3 * g3 + (1 / ell**2) * g3) == sp.zeros(3)

    # rho = Lam/(8 pi G) (T_mn = -rho g_mn convention), S = Om ell^{D-2}/4G;
    # eliminate ell
    rho = Lam / (8 * sp.pi * G)
    ell_from_S = (4 * G * S / Om)**(sp.S(1) / (D - 2))
    eps = rho.subs(ell, ell_from_S) * G**a_D
    c_D = (D - 1) * (D - 2) / (16 * sp.pi) * (Om / 4)**(sp.S(2) / (D - 2))
    ratio = sp.simplify(sp.powsimp(sp.expand_power_base(
        eps / (c_D * S**(-sp.S(2) / (D - 2))), force=True), force=True))
    check_identity_symbolic = ratio == 1
    check_identity_subs = all(
        sp.simplify(eps.subs(D, d) / (c_D * S**(-sp.S(2) / (D - 2))).subs(D, d)) == 1
        for d in range(3, 9))

    # the c_D table with the real sphere volumes, exact
    table = {}
    for d in range(3, 9):
        cd = sp.simplify(c_D.subs({D: d, Om: omega_sphere(d)}))
        table[d] = {"exact": sp.sstr(cd), "float": float(cd)}
    c4 = sp.simplify(c_D.subs({D: 4, Om: omega_sphere(4)}))
    check_c4 = c4 == sp.Rational(3, 8)
    c3 = sp.simplify(c_D.subs({D: 3, Om: omega_sphere(3)}))
    check_c3 = c3 == sp.pi / 32

    # the paper's D = 4 pair reproduced verbatim: eps = rho G^2 = (3/8)/S
    rho4 = (3 / ell**2) / (8 * sp.pi * G)
    S4 = (4 * sp.pi) * ell**2 / (4 * G)
    check_paper_pair = sp.simplify(rho4 * G**2 - sp.Rational(3, 8) / S4) == 0

    checks = dict(
        exponent_is_D_over_Dm2=check_aD, exponent_at_4_is_2=bool(check_a4),
        lambda_from_einstein_symbolic_D=check_lambda,
        ricci_dS4_explicit=bool(check_ric4), einstein_dS4=bool(check_einstein4),
        ricci_dS3_explicit=bool(check_ric3), einstein_dS3=bool(check_einstein3),
        identity_symbolic_D=bool(check_identity_symbolic),
        identity_substitution_D3_to_8=bool(check_identity_subs),
        c4_equals_3_over_8=bool(check_c4), c3_equals_pi_over_32=bool(check_c3),
        paper_D4_pair_reproduced=bool(check_paper_pair))
    return {
        "a_of_D": sp.sstr(a_D),
        "c_of_D": "(D-1)*(D-2)/(16*pi) * (Omega_{D-2}/4)**(2/(D-2))",
        "identity": "rho*G**(D/(D-2)) == c_D * S**(-2/(D-2))",
        "c_table": table,
        "checks": checks,
    }


# ---------------------------------------------------------------------------
# Piece 2: the fourfold D -> 2 degeneration  [D]
# ---------------------------------------------------------------------------
def piece2_d2_degeneration() -> dict:
    """Four separate, machine-checked failures of the direct continuation."""
    D = sp.Symbol('D', positive=True)
    G, S, ell = sp.symbols('G S ell', positive=True)
    a = sp.Symbol('a')

    # (i) no exponent: -D + a(D-2) = 0 has no solution at D = 2
    check_no_exponent = sp.solve(sp.Eq(-2 + a * 0, 0), a) == []

    # (ii) G_mn == 0 identically: arbitrary 2d metric, all three
    # components free functions of BOTH coordinates
    t, r = sp.symbols('t r')
    A = sp.Function('A')(t, r)
    B = sp.Function('B')(t, r)
    C = sp.Function('C')(t, r)
    g2 = sp.Matrix([[-A, C], [C, B]])
    Ric2 = ricci(g2, [t, r])
    R2 = scalar_curvature(g2, [t, r])
    check_einstein_2d_zero = sp.simplify(
        Ric2 - sp.Rational(1, 2) * R2 * g2) == sp.zeros(2)

    # ...and the 2d EH action is topological: int sqrt(g) R over the round
    # sphere of ANY radius is 8 pi (= 4 pi chi), no ell-dependence to vary
    th, ph = sp.symbols('theta phi_c')
    gs = sp.diag(ell**2, ell**2 * sp.sin(th)**2)
    Rs = sp.simplify(scalar_curvature(gs, [th, ph]))
    integral = sp.integrate(sp.integrate(
        Rs * sp.sqrt(gs.det()), (th, 0, sp.pi)), (ph, 0, 2 * sp.pi))
    check_gauss_bonnet = sp.simplify(integral - 8 * sp.pi) == 0

    # (iii) the continued Lam and c_D both die at fixed ell
    Lam = (D - 1) * (D - 2) / (2 * ell**2)
    check_lambda_limit = sp.limit(Lam, D, 2) == 0
    c_D_real = (D - 1) * (D - 2) / (16 * sp.pi) * (sp.S(2) / 4)**(sp.S(2) / (D - 2))
    check_cD_limit = sp.limit(c_D_real, D, 2, '+') == 0   # Omega_0 = 2

    # (iv) the elimination dies: S(ell) is ell-independent at D = 2
    S_of_ell = 2 * ell**(D - 2) / (4 * G)                  # Omega_0 = 2
    check_dS_dell_zero = sp.simplify(sp.diff(S_of_ell.subs(D, 2), ell)) == 0
    check_no_ell_solution = sp.solve(
        sp.Eq(S, S_of_ell.subs(D, 2)), ell) == []

    checks = dict(
        no_dimensionless_exponent=bool(check_no_exponent),
        einstein_tensor_2d_identically_zero=bool(check_einstein_2d_zero),
        eh_action_topological_gauss_bonnet=bool(check_gauss_bonnet),
        lambda_continuation_vanishes=bool(check_lambda_limit),
        c_D_limit_vanishes=bool(check_cD_limit),
        entropy_ell_independent_at_D2=bool(check_dS_dell_zero),
        elimination_step_has_no_solution=bool(check_no_ell_solution))
    return {
        "statement": ("The D >= 3 derivation has NO direct D = 2 continuation: "
                      "(i) no dimensionless rho*G**a exists ([G]=L**0, "
                      "[rho]=L**-2); (ii) G_mn == 0 identically, so "
                      "rho = Lam/8piG has no Einstein ancestor and the EH term "
                      "is topological; (iii) Lam(D,ell) -> 0 and c_D -> 0 at "
                      "fixed ell, while JT's 2/ell**2 lives in the dilaton "
                      "potential; (iv) S = Omega_0/4G is ell-independent, so "
                      "the elimination cannot run, and the real 2d entropy "
                      "S0 + phi_h/4G is dilaton data."),
        "checks": checks,
    }


# ---------------------------------------------------------------------------
# Piece 3a: the intrinsic JT-dS first law  [D], convention [C] = SVVV Eq. (26)
# ---------------------------------------------------------------------------
def piece3_first_law() -> dict:
    """delta S = -2 pi ell^2 rho = -delta E / T_dS, all machine-checked."""
    t, r = sp.symbols('t r')
    ell = sp.Symbol('ell', positive=True)
    G2 = sp.Symbol('G_2', positive=True)
    phi_r, rho, c0 = sp.symbols('phi_r rho c_0', real=True)
    f = 1 - r**2 / ell**2
    g = sp.diag(-f, 1 / f)
    coords = [t, r]

    def eom_residual(phi, Tmn):
        # -grad grad phi + g Box phi + (phi/ell^2) g = 8 pi G2 T_mn
        H, box = hessian_and_box(g, coords, phi)
        return sp.simplify(-H + g * box + (phi / ell**2) * g
                           - 8 * sp.pi * G2 * Tmn)

    # vacuum: phi = phi_r r/ell solves; a constant does NOT (no static
    # homogeneous constant mode, so phi_r fixes the solution)
    check_vacuum = eom_residual(phi_r * r / ell, sp.zeros(2)) == sp.zeros(2)
    const_resid = eom_residual(c0, sp.zeros(2))
    check_constant_not_homogeneous = sp.simplify(
        const_resid - (c0 / ell**2) * g) == sp.zeros(2) and c0 != 0

    # sourced: T_mn = -rho g_mn shifts the dilaton by -8 pi G2 rho ell^2
    Tmn = -rho * g
    phi_sourced = phi_r * r / ell - 8 * sp.pi * G2 * rho * ell**2
    check_sourced = eom_residual(phi_sourced, Tmn) == sp.zeros(2)

    # geometry rigid: the phi-variation gives R = 2/ell^2, matter-blind
    check_R = sp.simplify(scalar_curvature(g, coords) - 2 / ell**2) == 0

    # temperature from surface gravity at the horizon r = ell
    kappa = sp.simplify(sp.Abs(sp.diff(f, r).subs(r, ell)) / 2)
    T_dS = sp.simplify(kappa / (2 * sp.pi))
    check_T = sp.simplify(T_dS - 1 / (2 * sp.pi * ell)) == 0

    # Killing energy E = int_Sigma T_mn n^m xi^n sqrt(h) dr: redshifts cancel
    n_up = sp.Matrix([1 / sp.sqrt(f), 0])       # future unit normal n^mu
    xi_up = sp.Matrix([1, 0])                   # xi = d_t
    Tnxi = sp.simplify(sum(Tmn[m, n] * n_up[m] * xi_up[n]
                           for m in range(2) for n in range(2)))
    integrand = sp.simplify(Tnxi / sp.sqrt(f))  # sqrt(h) = 1/sqrt(f)
    check_redshift_cancel = sp.simplify(integrand - rho) == 0
    E = sp.integrate(integrand, (r, 0, ell))
    check_E = sp.simplify(E - rho * ell) == 0

    # the first law: delta S = delta phi_h/4G2 = -2 pi ell^2 rho = -E/T
    delta_phi_h = sp.simplify(phi_sourced.subs(r, ell) - phi_r)
    delta_S = sp.simplify(delta_phi_h / (4 * G2))
    check_delta_S = sp.simplify(delta_S + 2 * sp.pi * ell**2 * rho) == 0
    check_first_law = sp.simplify(delta_S + E / T_dS) == 0
    eps_JT = rho * ell**2
    check_eps_JT = sp.simplify(delta_S + 2 * sp.pi * eps_JT) == 0

    checks = dict(
        vacuum_solution=bool(check_vacuum),
        constant_is_not_homogeneous=bool(check_constant_not_homogeneous),
        sourced_solution=bool(check_sourced),
        geometry_rigid_R_2_over_ell2=bool(check_R),
        temperature_1_over_2pi_ell=bool(check_T),
        redshifts_cancel=bool(check_redshift_cancel),
        killing_energy_rho_ell=bool(check_E),
        delta_S_minus_2pi_ell2_rho=bool(check_delta_S),
        first_law_delta_S_eq_minus_dE_over_T=bool(check_first_law),
        delta_S_eq_minus_2pi_epsJT=bool(check_eps_JT))
    return {
        "convention": ("SVVV 2203.00700 Eq. (26): I = (1/16 pi G2) int "
                       "sqrt(-g) [(phi0+phi) R - 2 phi/ell**2] + bdy; "
                       "S = phi0/4G2 + phi_h/4G2 (Eq. 29); T = 1/2pi ell "
                       "(Eq. 16); vacuum phi = phi_r r/ell (after Eq. 15)"),
        "sourced_dilaton": sp.sstr(phi_sourced),
        "delta_S": sp.sstr(delta_S),
        "killing_energy": sp.sstr(E),
        "T_dS": sp.sstr(T_dS),
        "first_law": "delta_S = -delta_E/T_dS  (= SVVV Eq. (9) sign)",
        "eps_JT": "rho*ell**2",
        "delta_S_in_eps_JT": "-2*pi*eps_JT",
        "checks": checks,
    }


# ---------------------------------------------------------------------------
# Piece 3b: the dS3 -> JT reduction and the inherited relation  [D]
# ---------------------------------------------------------------------------
def piece3_reduction() -> dict:
    """R3 = R2 - 2 Box Phi/Phi exactly; the JT-dS it lands on; eps_3 pullback."""
    t, r, th = sp.symbols('t r theta')
    ell = sp.Symbol('ell', positive=True)
    G3, S = sp.symbols('G_3 S', positive=True)

    # the reduction identity, arbitrary diagonal base with full (t, r)
    # dependence (any 2d metric is locally diagonalizable; both sides are
    # base scalars, so this is fully general at the level used)
    A = sp.Function('A')(t, r)
    B = sp.Function('B')(t, r)
    Phi = sp.Function('Phi')(t, r)
    base = sp.Matrix([[-A, 0], [0, B]])
    g3 = sp.diag(-A, B, Phi**2)
    R3 = scalar_curvature(g3, [t, r, th])
    R2 = scalar_curvature(base, [t, r])
    _, boxPhi = hessian_and_box(base, [t, r], Phi)
    check_identity = sp.simplify(R3 - (R2 - 2 * boxPhi / Phi)) == 0

    # on the dS3 static patch the base is the dS2 static patch and Phi = r,
    # i.e. the JT vacuum solution with phi_r = ell
    check_dilaton_is_r = sp.simplify(
        (sp.Symbol('phi_r', real=True) * r / ell).subs(
            sp.Symbol('phi_r', real=True), ell) - r) == 0

    # action bookkeeping: (1/16 pi G3) int d^3x sqrt(-g3)(R3 - 2 Lam3)
    #   = (2 pi/16 pi G3) int d^2x sqrt(-g2) Phi (R2 - 2/ell^2) + total deriv
    # => SVVV Eq. (26) with phi = Phi, G2 = G3/2pi, phi0 = 0, ell_2 = ell_3.
    # sqrt(-g3) = sqrt(-g2) * Phi (Phi > 0) and int dtheta = 2 pi:
    check_measure = sp.simplify(g3.det() - base.det() * Phi**2) == 0
    # Phi * (2 Box Phi / Phi) = 2 Box Phi is the discarded total derivative;
    # what remains is exactly Phi*(R2 - 2 Lam3): an algebraic identity
    Lam3 = 1 / ell**2
    check_action = sp.simplify(
        Phi * (R3 - 2 * Lam3) - (Phi * (R2 - 2 * Lam3) - 2 * boxPhi)) == 0

    # entropies match exactly: S3 = 2 pi ell/4G3 = phi_h/4G2 with
    # phi_h = ell, G2 = G3/(2 pi)
    G2 = G3 / (2 * sp.pi)
    S3 = 2 * sp.pi * ell / (4 * G3)
    S_JT = ell / (4 * G2)
    check_entropy = sp.simplify(S3 - S_JT) == 0

    # the inherited relation: rho3 = Lam3/(8 pi G3); eps_3 = rho3 G3^3
    # equals (pi/32) S^-2 with S = S3.  This is Piece 1 at D = 3, now a
    # genuine JT-dS statement about the dilaton-potential coefficient.
    rho3 = Lam3 / (8 * sp.pi * G3)
    check_inherited = sp.simplify(
        rho3 * G3**3 - (sp.pi / 32) / S3**2) == 0

    checks = dict(
        reduction_identity_R3=bool(check_identity),
        dilaton_is_r_on_solution=bool(check_dilaton_is_r),
        measure_factor_Phi=bool(check_measure),
        action_lands_on_svvv_26=bool(check_action),
        entropies_match_exactly=bool(check_entropy),
        inherited_eps3_pi_over_32_S2=bool(check_inherited))
    return {
        "identity": "R_3 == R_2 - 2*Box_2(Phi)/Phi   (exact, no kinetic term)",
        "mapping": {"phi": "Phi (= r on the dS3 static patch)",
                    "G_2": "G_3/(2*pi)", "ell_2": "ell_3", "phi_0": 0,
                    "S_JT": "S_3 exactly"},
        "inherited": "rho_3*G_3**3 == (pi/32)*S**(-2)",
        "rho3_meaning": ("rho_3 = Lam_3/(8*pi*G_3) is the PARENT vacuum "
                         "energy; in JT variables it is the dilaton-potential "
                         "coefficient, not 2d matter"),
        "d4_route": ("[T] near-Nariai SdS4 only: MTY 1904.01911 Eqs. (3.37) "
                     "ds^2 = rho_e^2 ds^2_dS2 + rho_e^2 (1+delta)^2 dOmega_2^2 "
                     "and (3.39) phi_0 = R_dS4^2 rho_e^2/4G_N, phi = 2 phi_0 "
                     "delta; inherits the D = 4 exponent -1 with the parent's "
                     "S; the empty-dS4 sphere reduction is NOT JT (it "
                     "generates dilaton kinetic terms)"),
        "checks": checks,
    }


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------
def condensate_link() -> dict:
    """The measured 2d vacuum-energy shift, pushed through the first law."""
    path = HERE.parent / "thooft-target" / "condensate.json"
    eps0 = None
    if path.exists():
        with open(path) as fh:
            eps0 = json.load(fh)["eps_duality_point"]
    return {
        "eps0_measured_Nc_g2": eps0,
        "prediction": ("delta_S = -2*pi*ell**2*rho = "
                       "+2*pi*|eps0|*Nc*(g*ell)**2 > 0; finite, first-order, "
                       "divergence-free"),
        "hedges": ["conditional on JT-dS being the bulk (the paper's "
                   "assumption 1, unproved)",
                   "conditional on minimal coupling to g only (matter blind "
                   "to phi)",
                   "(g*ell)**2 is a dictionary knob, currently unfixed; "
                   "under the paper's flat-space limit it vanishes per mode",
                   "the zero point E_vac(-1) is a subtraction convention"],
        "not_conditional": ("the shift is finite and first-order: nothing "
                            "for Holographic Tuning to cancel, on either "
                            "side of the first law"),
    }


def main() -> dict:
    t0 = time.time()
    p1 = piece1_general_identity()
    p2 = piece2_d2_degeneration()
    p3a = piece3_first_law()
    p3b = piece3_reduction()
    out = {
        "provenance": ("jt_ds_vacuum.py -- the JT-dS analog of eps = 3/(8S): "
                       "general-D identity derived and certified; D -> 2 "
                       "degeneration as a fourfold negative statement; JT-dS "
                       "first law + dS3 reduction machine-checked; labels "
                       "[D]/[T]/[C] per docstring"),
        "paper_statement": {
            "arxiv": "2607.05678",
            "file": "papers/arXiv-2607.05678v1/HologramsSM.tex",
            "section": "Fine Tuning?",
            "lines": {"eps_rho_G2": [1014, 1016], "eps_3_8S": [1018, 1020]},
            "displays": ["\\epsilon = \\rho G^2",
                         "\\epsilon = \\frac{3}{8 S}"],
            "note": ("both displays unnumbered \\be...\\ee; papers/ is "
                     "gitignored, source local-only, test skips if absent"),
        },
        "conventions": {
            "signature": "(-,+,...,+), spheres and dS have R > 0",
            "einstein": "G_mn + Lam*g_mn = 8*pi*G*T_mn",
            "vacuum_stress": "T_mn = -rho*g_mn  =>  rho = Lam/(8*pi*G)",
            "entropy_D": "S = Omega_{D-2}*ell**(D-2)/(4*G)",
            "jt": "SVVV 2203.00700 Eq. (26); S: Eq. (29); T: Eq. (16)",
        },
        "general_D": p1,
        "d2_degeneration": p2,
        "jt_first_law": p3a,
        "reduction_3d": p3b,
        "nonuniqueness": {
            "statement": ("'the JT-dS analog' is (a) the intrinsic first law "
                          "delta_S = -2*pi*rho*ell**2 PLUS (b) the "
                          "reduction-inherited family eps_D = c_D * "
                          "S**(-2/(D-2)) whose exponent the intrinsic JT data "
                          "(phi_0, phi_r, ell) cannot choose: S**-2 from dS3, "
                          "S**-1 from near-Nariai dS4.  The parent choice is "
                          "a dictionary input."),
            "alternatives": [
                "[C] parent D = 3: eps_3 = rho_3*G_3**3 = (pi/32)/S**2, "
                "machine-checked end to end",
                "[C] parent D = 4: exponent -1 with the parent's S, via "
                "near-Nariai only (MTY Eqs. 3.37, 3.39)",
                "[C] G_eff = G_2/phi_h so S = 1/(4*G_eff) mimics area/4G: "
                "relabels which symbol carries dimension, changes nothing",
                "[C] phi, G_2 dimensionless (SVVV) vs length-dimension "
                "(reduction): normalization only",
            ],
        },
        "literature": {
            "fetch_method": ("ar5iv HTML, 2026-08-27, two independent "
                             "extractions per paper, agreeing; equation "
                             "numbers cited only where both returned the "
                             "equation verbatim"),
            "svvv_2203_00700": {
                "title": ("Quasi-local energy and microcanonical entropy in "
                          "two-dimensional nearly de Sitter gravity"),
                "verified": ["Eq. (9): -delta H_xi = T_GH delta S_GH "
                             "(minus-sign patch first law)",
                             "Eq. (16): T_GH = 1/(2*pi*L)",
                             "Eq. (26): the JT action used here",
                             "Eq. (29): S = phi_0/4G_2 + phi_h/4G_2",
                             "after Eq. (15): Phi = Phi_r r/L solves the EOM"],
            },
            "mty_1904_01911": {
                "title": "Two dimensional Nearly de Sitter gravity",
                "verified": ["Eq. (2.1): action, ell = 1, 1/4G absorbed",
                             "Eq. (3.37): near-Nariai reduction ansatz",
                             "Eq. (3.39): phi_0 = R_dS4^2 rho_e^2/4G_N, "
                             "phi = 2 phi_0 delta",
                             "around Eq. (2.13): S_{0,dS} = 8*pi*phi_0 (both "
                             "extractions place it at/just after (2.13))"],
            },
            "background_only": ["none: every citation above was fetched"],
        },
        "condensate_link": condensate_link(),
        "runtime_s": round(time.time() - t0, 2),
    }
    return out


if __name__ == "__main__":
    result = main()
    all_checks = {}
    for key in ("general_D", "d2_degeneration", "jt_first_law", "reduction_3d"):
        all_checks.update({f"{key}.{k}": v
                           for k, v in result[key]["checks"].items()})
    assert all(all_checks.values()), {k: v for k, v in all_checks.items()
                                      if not v}
    out_path = HERE / "jt_ds_vacuum.json"
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=1)
        fh.write("\n")
    print(f"all {len(all_checks)} machine checks pass")
    print(f"a(D) = {result['general_D']['a_of_D']},  "
          f"c_4 = {result['general_D']['c_table'][4]['exact']},  "
          f"c_3 = {result['general_D']['c_table'][3]['exact']}")
    print(f"delta_S = {result['jt_first_law']['delta_S']}  "
          f"(= -delta_E/T_dS, E = {result['jt_first_law']['killing_energy']})")
    print(f"inherited: {result['reduction_3d']['inherited']}")
    print(f"wrote {out_path}  ({result['runtime_s']}s)")
