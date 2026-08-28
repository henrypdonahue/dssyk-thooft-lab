#!/usr/bin/env python3
r"""
The charged (Dirac / complex-SYK) chord layer over the Majorana engine.
=======================================================================

Double-scaled COMPLEX SYK is the U(1)-charged Dirac model of
syk-self-averaging/dirac.py.  Berkooz-Narovlansky-Raj solve it by chord
combinatorics in "Complex Sachdev-Ye-Kitaev model in the double scaling
limit", arXiv:2006.13983 (JHEP 02 (2021) 113).  Called BNR below.

Every rule this module encodes is equation-checked against that paper.

* Model [BNR (2.2), (2.4)]: H = sum J_I^{I'} psibar_{I'} psi_I, with
  |I| = |I'| = p_BNR creation and p_BNR annihilation operators per term,
  and lambda_BNR = p_BNR^2/N fixed [BNR (2.8)].  dirac.py's conventions
  (Nc modes, k = p/2 creators and k annihilators, {c, c^dag} = 1) map by
  p_BNR = p/2, N = Nc, psi = sqrt(2) c.

* mu = 0 moments [BNR (4.3)]:
      m_k = (J e^{lambda_BNR/2})^k  sum_{chord diagrams} q^{crossings},
      q = e^{-4 lambda_BNR}.
  The chord-diagram sum is IDENTICAL to the Majorana one.  Same transfer
  matrix T|n> = |n+1> + [n]_q |n-1> [review 2407.09396 (2.15)-(2.16)].
  Only the value of q changes.  So this module builds on duality/chord.py
  instead of duplicating it.

* Fixed specific charge Qs = (Q - Nc/2)/Nc [BNR (5.16)-(5.17)]:
      q(Qs) = exp(-4 lambda_BNR / (1 - 4 Qs^2)),
  with the same chord sum per sector.  Corroborated by arXiv:2512.07715,
  which writes the grand-canonical transfer matrix and shows it is block
  diagonal with the canonical sector matrices.

* Charged random probes [BNR (7.1)-(7.2)]: a creators and b annihilators,
  with couplings INDEPENDENT of H's.  An H-chord crossing the matter
  chord weighs qtilde = q(Qs)^{(a+b)/(2 p_BNR)} [BNR (7.17)].  That is
  the Majorana matter rule with

      Delta = (a + b)/p            (repo p = 2 p_BNR).

  The fixed-charge 2-pt function [BNR (7.18), (7.23)] is the Majorana
  theta-integral [review (2.30)] at (q(Qs), Delta), with per-side energy
  rescalings e^{B_1}, e^{B_2} [BNR (7.21)-(7.22)] that vanish at Qs = 0
  or a = b.

THE CONVENTIONS BRIDGE (the Dirac analog of the factor-2 landmine)
------------------------------------------------------------------
For dirac.py's ensemble, with repo p = total operators per term,

    q_Dirac = exp(-4 (p/2)^2 / Nc) = exp(-p^2/Nc) = exp(-lambda_paper).

So on the DIRAC side the paper's lambda = p^2/N (2607.05678 Sec. 4.1) IS
the chord decay constant.  There is no factor-2 clash here.

Equivalently, q_Dirac equals the Majorana-literature q at (N, p) =
(2 Nc, p): 2 Nc Majorana modes make Nc complex modes, and
lambda_chord(2 Nc, p) = 2 p^2/(2 Nc) = p^2/Nc.  Both statements are
tested against dictionary.py's converters.

VERIFICATION (test_dirac_chord.py; the measured numbers are in the
docstrings there)
-----------------------------------------------------------------------
* Microscopic leg, the Dirac Krawtchouk analog.  The exact finite-Nc
  crossed/uncrossed trace average over random monomial 4-words tends to
  e^{-p^2/Nc}.  Seeded compressed-register sampling; no chords and no
  transfer matrix are involved.
* ED leg, gate G1.  Disorder-averaged Tr H^{2k}/dim ratios of
  dirac.dirac_hamiltonian ensembles converge onto <0|T^{2k}|0> at
  q = q_Dirac(Nc) as Nc grows (c_symmetric ensemble).  The generic
  ensemble shares the limit but converges more slowly, because index
  coincidences inside a term are O(k^2/Nc) frequent.
* Conserved leg, gate G2.  M_1 = -(p/2) i H exactly (dirac.py), so its
  chord image is -(p/2) i sigma_H T, an insertion of the transfer matrix
  itself.  Its 2-pt function is time-independent, and its Lehmann weight
  away from omega = 0 vanishes identically.
* Census leg, gate G3 pair half.  The lone-crossing factor of a BARE
  c_i / c^dag_i chord against a c_symmetric H-chord is exactly 1 - p/Nc
  (derivation below), measured per sample at Nc = 16, 32, 64 with the
  compressed-register sampler, converging on q_Dirac^{1/p} at the census
  rate p^2/(2 Nc^2).
* Theorem-transfer leg, gate G3 pair half.  chord_charges' recursion,
  dressed-propagator and charge identities re-verified at
  (q_Dirac(Nc, p), Delta = 1/p), plus the word-level moment zeros through
  the independent chord_moments assembly at the Dirac q.

GATE G3, SPLIT: the H-correlated tower M_{n>=2}, pair against diagonal
----------------------------------------------------------------------
BNR's matter rules cover probes whose couplings are drawn INDEPENDENTLY
of the Hamiltonian's [BNR Sec. 7.1].  The repo's tower operators
M_n = sum_i c^dag_i (d^n/dt^n) c_i for n >= 2 are polynomials in the
SAME couplings as H (M_2 = sum x_u x_v G_uv, exact_wick_dirac.py), and
no published source provides H-correlated chord rules.  Checked: BNR
2006.13983; the review 2407.09396, whose full text mentions complex SYK
only in its reference list; the U(M)-symmetric chords of 2302.11882; and
the transfer-matrix construction of 2512.07715.  So
MN_TOWER_CHORD_RULES_PUBLISHED stays False.

THE PAIR HALF NEVER NEEDED SUCH RULES.  On the Majorana side,
chord_moments.py expands A_n = sum_i psi_i (ad_H)^n psi_i binomially
into explicit words H^a psi H^b psi ... whose matter chords are BARE
single fermions, and the conserved-channel theorem of chord_charges.py
is a statement about (q, Delta) alone.  The Dirac tower has the
identical word structure -- (ad_H)^n c = sum_r (-1)^r C(n,r) H^{n-r} c H^r,
with the i^n from d/dt = i ad_H an overall constant -- so the i != j
pair moments reduce to words Tr[H^a c^dag_i H^b c_i H^c c^dag_j H^d c_j]
with two bare, sequential, mutually non-crossing c-chords.  The one
microscopic input the transfer needs is the crossing weight of a bare
c chord against an H-chord.

THE CENSUS (derived here, in the style of exact_wick_dirac.py Part 1;
measured per sample in the tests).  The c_symmetric H-chord endpoint is
h_u = c^dag_A c_B + c^dag_B c_A = X + X^dag, A n B = empty,
|A| = |B| = k = p/2.  h_u^2 = prod_A n prod_B (1-n) +
prod_B n prod_A (1-n) is diagonal, so the uncrossed word has
tr(c_i c^dag_i h_u h_u)/dim = 2^{-2k} for EVERY u, i included or not.
The crossed word tr(c_i h_u c^dag_i h_u) splits into two cases.
i outside A u B: c^dag_i moves through the even h_u with no
Jordan-Wigner sign, so crossed = uncrossed exactly.  i inside A u B:
all four cross terms die -- (X, X) and (X^dag, X^dag) are unbalanced in
mode i's counts, and in (X, X^dag) and (X^dag, X) an occupation clash
kills the word (for i in A, the outer c^dag_A finds mode i refilled by
the inner c^dag_i, or c^dag_i finds it already filled by c^dag_A; i in B
mirrors under h.c.).  So the crossed word vanishes IDENTICALLY there,
and the ensemble-averaged lone-crossing factor is EXACTLY the avoidance
probability

    R(Nc) = 1 - p/Nc = e^{-p/Nc} + O(p^2/Nc^2)  ->  q_Dirac^{1/p}.

Consistency: the Majorana census is the sign average
E[(-1)^{|{i} n I|}] = 1 - 2p/N, and at N = 2 Nc that equals 1 - p/Nc.
The Dirac ensemble KILLS the crossed word where the Majorana one
sign-flips it, and the two censuses agree under the same N = 2 Nc bridge
as q itself.  So the Dirac single-fermion chord is BARE with
Delta = 1/p at q = q_Dirac: the same (q, Delta) evaluation map as BNR's
independent probes [BNR (7.17)], now established for the bare fermion
the tower words actually contain, where BNR's independence assumption is
unavailable.

CONSEQUENCE, the Dirac conserved pair channel.  chord_charges' theorem
applies verbatim at (q_Dirac(Nc, p), Delta = 1/p): the chord
representation of every Dirac pair-channel insertion is a polynomial in
the transfer matrix (xhat_dirac below), with

    Xhat_1 = -(1 - q^{1/p}) T ,
    Xhat_2 = (1 - q^{1/p})^2 T^2 + (1 - q^{2/p}) ,

so every spectral moment of every M_n PAIR (i != j) channel vanishes at
Nc = infinity, at every n, both orderings, odd k included.  n = 1 lands
on gate G2's exact law M_1 = -(p/2) i H.  The per-word normalization
tr(c_i c^dag_i)/dim = 1/2 per matter chord (psi = sqrt(2) c) is uniform
across a channel's words, so it cannot disturb the zeros.  The theorem
side holds at ANY q, so it covers the fixed-charge sectors q(Qs) too;
the census itself is the unprojected trace, i.e. the Qs = 0 statement.
And the census is derived in the c_symmetric ensemble -- the CP-sharp
ensemble the moments pipeline actually uses; the generic ensemble shares
the Nc = infinity limit (gate G1's pattern) but its census is not
written here.

STILL NOT established, and said plainly: the flavor-diagonal (i = j)
half of the channel -- which on the Majorana side carries ALL of the
tower's N = infinity weight (diagonal_chord.py) -- has no Dirac
derivation here, and the CP RESOLUTION of the Dirac M_n (the
cp_projections split that the CP-resolved data is actually waiting on)
is untouched.  Gate G3's comparison to the moments_pipeline data therefore
remains open; this module closes its pair half only.  The
neutral-bilinear channel below stays what it was: BNR's INDEPENDENT
random bilinear (a = b = 1, Delta = 2/p), a legitimate published object
that is not the M_n tower.
"""

from __future__ import annotations

import math

import numpy as np

from chord import (Contour, Sector0, g2_closed_form, rho_theta,
                   energy_theta, s_omega_smooth, vertex2_matrix)
from chord_charges import transfer_formal, xhat_matrix

# what gate G3 would need from the LITERATURE and does not have; tested so
# the claim is pinned
MN_TOWER_CHORD_RULES_PUBLISHED = False
# the PAIR (i != j) half of gate G3 is derived in-repo instead: the
# bare-fermion census (module docstring) + the conserved-channel theorem
# (chord_charges.py) evaluated at (q_dirac, Delta = 1/p).  The
# flavor-diagonal half and the CP resolution remain open.
MN_TOWER_PAIR_CHANNEL_DERIVED = True


# ---------------------------------------------------------------------------
# conventions bridge
# ---------------------------------------------------------------------------
def lambda_dirac(Nc: int, p: int) -> float:
    """The Dirac chord decay constant lambda with q = e^{-lambda}:
    lambda = 4 lambda_BNR = 4 (p/2)^2/Nc = p^2/Nc  [BNR (2.8), (4.3)].
    Numerically equal to dictionary.lambda_paper(Nc, p) and to
    dictionary.lambda_chord(2 Nc, p); both identities are tested."""
    return p * p / float(Nc)


def q_dirac(Nc: int, p: int) -> float:
    """q = exp(-p^2/Nc): the mu = 0 (equivalently Qs = 0 sector) H-chord
    crossing weight of the Dirac model [BNR (4.3)]."""
    return math.exp(-lambda_dirac(Nc, p))


def q_dirac_charge(Nc: int, p: int, qs: float) -> float:
    """Crossing weight in the fixed specific-charge sector
    q(Qs) = exp(-4 lambda_BNR/(1 - 4 Qs^2))  [BNR (5.17), (7.17)].
    Qs = (Q_repo - Nc/2)/Nc, Qs in (-1/2, 1/2); Qs = 0 is half filling."""
    if not -0.5 < qs < 0.5:
        raise ValueError("specific charge must satisfy |Qs| < 1/2")
    return math.exp(-lambda_dirac(Nc, p) / (1.0 - 4.0 * qs * qs))


def delta_charged(a: int, b: int, p: int) -> float:
    """Chord dimension of a random probe with a creators + b annihilators:
    Delta = (a + b)/p (repo p), from qtilde = q^{(p_M + pbar_M)/(2 p_BNR)}
    [BNR (7.17)].  Neutral bilinear (a = b = 1): Delta = 2/p."""
    if a < 0 or b < 0 or a + b == 0:
        raise ValueError("need a, b >= 0 with a + b > 0")
    return (a + b) / float(p)


def probe_is_fermionic(a: int, b: int) -> bool:
    """Matter-matter crossings carry (-1)^{(a1+b1)(a2+b2)} [BNR (3.10) with
    p_ij = 0 at large N], i.e. the transposition sign iff both probes have
    odd total fermion number.  Feed the odd ones to Contour's `fermionic`."""
    return (a + b) % 2 == 1


def charge_sector_asymmetry(Nc: int, p: int, a: int, b: int, qs: float):
    """(e^{B1}, e^{B2}): per-H-insertion energy rescalings on the two sides
    of a charged probe in the fixed-charge 2-pt function [BNR (7.21)-(7.22)],
    mapped to repo conventions (p_BNR = p/2, p_M = b, pbar_M = a -- BNR's
    psi is the annihilator -- and n = -Qs Nc):

        B1 =  (p/2) p_M/Nc * 2n/(Nc-2n) - (p/2) pbar_M/Nc * 2n/(Nc+2n)
        B2 = -(p/2) p_M/Nc * 2n/(Nc+2n) + (p/2) pbar_M/Nc * 2n/(Nc-2n)

    Both are 1 at Qs = 0, and swap under a <-> b together with n -> -n
    (charge conjugation).  These are BNR's O(N^0) saddle factors; treat
    them as the stated-order asymmetry, not an all-orders finite-Nc law."""
    n = -qs * Nc
    half_p = 0.5 * p
    pm, pbm = b, a
    b1 = (half_p * pm / Nc * (2.0 * n / (Nc - 2.0 * n))
          - half_p * pbm / Nc * (2.0 * n / (Nc + 2.0 * n)))
    b2 = (-half_p * pm / Nc * (2.0 * n / (Nc + 2.0 * n))
          + half_p * pbm / Nc * (2.0 * n / (Nc - 2.0 * n)))
    return math.exp(b1), math.exp(b2)


# ---------------------------------------------------------------------------
# units bridge: exact ensemble means of Tr H^2/dim for dirac.py's ensembles
# ---------------------------------------------------------------------------
def m2_mean_dirac(Nc: int, p: int, ensemble: str = "c_symmetric") -> float:
    """E[Tr H^2/dim] closed forms for dirac.draw_couplings ensembles, with
    sigma^2 = coupling_variance_dirac = (k!)^2 2^p / Nc^{p-1}, k = p/2.

    c_symmetric (disjoint real pairs; the moments-pipeline ensemble):
        sigma^2 C(Nc,k) C(Nc-k,k) 4^{-k}
    (the qed_campaign/exact_wick_dirac closed form, restated).

    generic (all ordered pairs (A,B), complex Hermitian):
        sigma^2 C(Nc,k) sum_{m=0}^{k} C(k,m) C(Nc-k,k-m) 2^{m-2k},
    m = |A cap B|; each ordered pair contributes sigma^2 2^{-|A u B|}
    (dirac.coupling_variance_dirac docstring), |A u B| = 2k - m.
    Both forms are verified against seeded ED means in the tests."""
    k = p // 2
    sigma2 = (math.factorial(k) ** 2) * 2 ** p / float(Nc) ** (p - 1)
    if ensemble == "c_symmetric":
        return sigma2 * math.comb(Nc, k) * math.comb(Nc - k, k) / 4.0 ** k
    if ensemble == "generic":
        tot = sum(math.comb(k, m) * math.comb(Nc - k, k - m) * 2.0 ** (m - 2 * k)
                  for m in range(k + 1))
        return sigma2 * math.comb(Nc, k) * tot
    raise ValueError(f"unknown ensemble {ensemble!r}")


def sigma_H_dirac(Nc: int, p: int, ensemble: str = "c_symmetric") -> float:
    """sqrt(E[Tr H^2/dim]): the repo-units <-> chord-units energy scale,
    t_chord = t_repo * sigma_H (the analog of chord.sigma_H)."""
    return math.sqrt(m2_mean_dirac(Nc, p, ensemble))


# ---------------------------------------------------------------------------
# the charged transfer matrix and Hamiltonian moments (gate G1 machinery)
# ---------------------------------------------------------------------------
def dirac_sector0(nmax: int, Nc: int, p: int, qs: float = 0.0) -> Sector0:
    """The Dirac chord transfer matrix in the Qs sector: chord.Sector0 at
    q = q_dirac_charge(Nc, p, Qs).  Identical structure to the Majorana
    T|n> = |n+1> + [n]_q |n-1> because the chord-diagram sum is identical
    [BNR (4.3), (5.13); blocks of the grand-canonical matrix per
    2512.07715] -- only q changes."""
    return Sector0(nmax, q_dirac_charge(Nc, p, qs))


def hamiltonian_moment_ratios(q: float, kmax: int = 3, nmax: int = 24) -> dict:
    """Normalized chord moments r_{2k} = m_{2k}/m_2^k = <0|T^{2k}|0> for
    k <= kmax, by exact integer T-powers (no exponentials, no closed
    forms).  In chord units m_2 = 1; repo units multiply m_{2k} by
    sigma_H_dirac^{2k}.  Truncation nmax >= 2*kmax is exact for these
    moments (T^k reaches chord number k at most)."""
    s = Sector0(nmax, q)
    v = s.vac()
    out = {}
    for j in range(1, 2 * kmax + 1):
        v = s.apply_T(v)
        if j % 2 == 0:
            out[j] = float(v[0].real)
    return out


def dirac_contour(nmax: int, Nc: int, p: int, probes: dict,
                  qs: float = 0.0) -> Contour:
    """chord.Contour over the Dirac q, with matter dimensions derived from
    the probes' fermion content.  probes = {label: (a, b)} with a creators
    and b annihilators per label [BNR (7.1)].  Odd-total labels get the
    fermionic crossing sign.  All Contour topologies (2-pt, 4-pt, folds,
    complex contours) are then available at the Dirac parameters."""
    deltas = {lab: delta_charged(a, b, p) for lab, (a, b) in probes.items()}
    ferm = {lab for lab, (a, b) in probes.items() if probe_is_fermionic(a, b)}
    return Contour(q_dirac_charge(Nc, p, qs), nmax, deltas, fermionic=ferm)


# ---------------------------------------------------------------------------
# gate G2: the conserved M_1 channel
# ---------------------------------------------------------------------------
def m1_channel_correlator(ts: np.ndarray, Nc: int, p: int, nmax: int = 60,
                          beta: float = 0.0, qs: float = 0.0) -> np.ndarray:
    """G_{M1}(t)/|c_1|^2 sigma_H^2 = <T(t) T(0)>_beta in chord units, where
    the chord image of M_1 = -(p/2) i H (exact operator identity,
    dirac.py) is c_1 T with c_1 = -(p/2) i sigma_H: an insertion of the
    transfer matrix itself, NOT a matter chord.  Computed by explicit
    evolution (eigen-decomposition; commutativity is never invoked), so
    time-independence is a measured outcome, not an input:

        G(t) = <0| e^{-(beta/2)T} e^{+itT} T e^{-itT} T e^{-(beta/2)T} |0>
               / <0| e^{-beta T} |0>.

    The M_1 channel is conserved: G(t) = G(0) for all t, so its spectral
    weight sits entirely at omega = 0 (gate G2)."""
    s = dirac_sector0(nmax, Nc, p, qs)
    v0 = s.vac()                               # chord basis throughout
    if beta != 0.0:
        v0 = s.evolve(v0, 0.5 * beta)
    z = complex(np.vdot(v0, v0))               # <0|e^{-beta T}|0>
    ts = np.atleast_1d(ts)
    out = np.empty(len(ts), dtype=complex)
    for i, t in enumerate(ts):
        w = s.apply_T(v0)                      # T |v0>
        w = s.evolve(w, +1j * t)               # e^{-itT} T |v0>   (u = +it)
        w = s.apply_T(w)                       # T e^{-itT} T |v0>
        w = s.evolve(w, -1j * t)               # e^{+itT} T e^{-itT} T |v0>
        out[i] = complex(np.vdot(v0, w)) / z
    return out


def m1_spectral_lines(Nc: int, p: int, nmax: int = 60, beta: float = 0.0,
                      qs: float = 0.0):
    """Lehmann lines of the M_1 channel through the SAME generic machinery
    used for matter channels (chord.spectral_lines pattern with the matter
    matrix replaced by T): (omega_jk, w_jk) with omega_jk = E_j - E_k.
    Because T is diagonal in its own eigenbasis, every off-diagonal weight
    vanishes identically -- the gate-G2 statement in spectral form."""
    s = dirac_sector0(nmax, Nc, p, qs)
    c = s.V.T @ s.vac().real
    if beta != 0.0:
        c = c * np.exp(-0.5 * beta * s.E)
    # rotate the chord-space tridiagonal T into the theta basis through the
    # same algebra a matter channel uses; <theta_j|T|theta_k> = E_j delta_jk
    # is then a RESULT of the diagonalization, not an input
    Tchord = np.zeros((nmax + 1, nmax + 1))
    Tchord[np.arange(1, nmax + 1), np.arange(nmax)] = s.off
    Tchord[np.arange(nmax), np.arange(1, nmax + 1)] = s.off
    M = s.V.T @ Tchord @ s.V
    W = np.outer(c, c) * M / float(c @ c)
    om = np.subtract.outer(s.E, s.E)
    return om.ravel(), W.ravel()


def m1_offdiagonal_weight(Nc: int, p: int, nmax: int = 60,
                          beta: float = 0.0) -> float:
    """Fraction of the M_1-channel Lehmann weight at omega != 0.  Exactly 0
    (to matrix-algebra precision): the conserved-channel gate number."""
    om, w = m1_spectral_lines(Nc, p, nmax, beta)
    off = np.abs(w[np.abs(om) > 1e-9]).sum()
    tot = np.abs(w).sum()
    return float(off / tot) if tot > 0 else 0.0


# ---------------------------------------------------------------------------
# gate G3, pair half: the Dirac conserved pair channel
# ---------------------------------------------------------------------------
def bare_fermion_crossing(Nc: int, p: int) -> float:
    """Exact finite-Nc lone-crossing factor of a BARE c_i / c^dag_i chord
    against one c_symmetric H-chord,

        E[tr(c_i h_u c^dag_i h_u)] / E[tr(c_i c^dag_i h_u h_u)] = 1 - p/Nc :

    the crossed word vanishes identically when i lies in the H-chord's
    support and equals the uncrossed word exactly when it does not
    (census, module docstring).  DERIVED; measured per sample in the
    tests.  Nc -> infinity limit e^{-p/Nc} = q_dirac^{1/p}, deviation
    p^2/(2 Nc^2) + O(p^3/Nc^3).  Equals the Majorana sign census
    1 - 2p/N at N = 2 Nc."""
    if not 0 < p <= Nc:
        raise ValueError("need 0 < p <= Nc")
    return 1.0 - p / float(Nc)


def xhat_dirac(n: int, Nc: int, p: int, nmax: int = 24) -> np.ndarray:
    """The chord representation of the Dirac pair-channel insertion of
    M_n's i != j amplitude: chord_charges.xhat_matrix at
    (q_dirac(Nc, p), Delta = 1/p), in the formal chord basis
    T|m> = |m+1> + [m]_q |m-1>.  A polynomial in T by the
    conserved-channel theorem, so every M_n pair moment vanishes at
    Nc = infinity.  Closed forms: xhat_dirac(1) = -(1 - q^{1/p}) T (the
    Nc = infinity image of the exact law M_1 = -(p/2) i H, gate G2) and
    xhat_dirac(2) = (1 - q^{1/p})^2 T^2 + (1 - q^{2/p})."""
    return xhat_matrix(n, nmax, q_dirac(Nc, p), 1.0 / p)


# ---------------------------------------------------------------------------
# charged-probe two-point functions [BNR (7.18), (7.23)]
# ---------------------------------------------------------------------------
def g2_charged(beta1: complex, beta2: complex, Nc: int, p: int, a: int,
               b: int, qs: float = 0.0, n_theta: int = 200) -> complex:
    """Unnormalized fixed-charge-sector Euclidean 2-pt amplitude of a random
    charged probe (a creators, b annihilators), in chord units:

        G = int dth1 dth2 rho rho e^{-beta1 s1 E1 - beta2 s2 E2} K(th1,th2),

    with rho and E at q = q(Qs), K the (2.30)-type vertex at
    Delta = (a+b)/p -- identical to BNR (7.23)'s
    (qtilde^2; q)_inf / (qtilde e^{i(+-th1 +- th2)}; q)_inf -- and
    (s1, s2) = charge_sector_asymmetry the per-side BNR (7.21)-(7.22)
    rescalings.  At Qs = 0 this reduces EXACTLY to chord.g2_closed_form
    at (q_dirac, Delta) (asserted in the tests).  Divide by the sector
    partition function amplitude for normalized correlators; the
    sector-uniform scale factors [BNR (7.19)-(7.20)] then drop."""
    q = q_dirac_charge(Nc, p, qs)
    delta = delta_charged(a, b, p)
    s1, s2 = charge_sector_asymmetry(Nc, p, a, b, qs)
    x, wq = np.polynomial.legendre.leggauss(n_theta)
    th = 0.5 * np.pi * (x + 1.0)
    w = 0.5 * np.pi * wq
    rho = rho_theta(th, q)
    E = energy_theta(th, q)
    V = vertex2_matrix(th, q, delta)
    f1 = w * rho * np.exp(-np.asarray(beta1) * s1 * E)
    f2 = w * rho * np.exp(-np.asarray(beta2) * s2 * E)
    return complex(f1 @ V @ f2)


def neutral_bilinear_spectrum(omegas: np.ndarray, Nc: int, p: int,
                              n_theta: int = 2000) -> np.ndarray:
    """S(omega) of BNR's INDEPENDENT random neutral bilinear probe
    (a = b = 1, Delta = 2/p) at T_B = infinity in the Qs = 0 sector: the
    exact smooth spectral function chord.s_omega_smooth at
    (q_dirac, 2/p).  Moment closures (tested): int S = 1,
    int omega^2 S = 2(1 - q^{2/p}).

    NOT the M_n tower: the tower's couplings are H-correlated and its
    chord rules are unpublished (module docstring, gate G3)."""
    return s_omega_smooth(np.atleast_1d(omegas), q_dirac(Nc, p),
                          delta_charged(1, 1, p), n_theta)


def neutral_bilinear_moments(Nc: int, p: int) -> dict:
    """Closed-form chord-unit moments of the independent neutral bilinear
    channel: mu0 = 1, mu2 = 2(1 - q^Delta), Delta = 2/p (the review's
    tested normalization, at the Dirac q)."""
    q = q_dirac(Nc, p)
    d = delta_charged(1, 1, p)
    return dict(delta=d, q=q, mu0=1.0, mu2=2.0 * (1.0 - q ** d))


if __name__ == "__main__":
    Nc, p = 12, 4
    q = q_dirac(Nc, p)
    print(f"Dirac chord layer: Nc={Nc}, p={p} "
          f"(k={p//2} creators + {p//2} annihilators per term)")
    print(f"  lambda_dirac = p^2/Nc = {lambda_dirac(Nc, p):.6f}  "
          f"q = {q:.6f}")
    r = hamiltonian_moment_ratios(q)
    print(f"  chord moment ratios: r4 = {r[4]:.6f} (= 2+q = {2+q:.6f}), "
          f"r6 = {r[6]:.6f} (= 5+6q+3q^2+q^3 = "
          f"{5+6*q+3*q*q+q**3:.6f})")
    ts = np.linspace(0.0, 4.0, 9)
    g = m1_channel_correlator(ts, Nc, p, nmax=40)
    print(f"  M_1 channel |G(t) - G(0)| max = "
          f"{np.abs(g - g[0]).max():.2e}  (conserved; gate G2)")
    print(f"  M_1 off-omega=0 Lehmann weight = "
          f"{m1_offdiagonal_weight(Nc, p, nmax=40):.2e}")
    g0 = g2_charged(0.4, 0.3, Nc, p, 1, 1)
    gm = g2_closed_form(0.4, 0.3, q, delta_charged(1, 1, p))
    print(f"  g2 charged (Qs=0) vs Majorana closed form: "
          f"{g0.real:.10f} vs {gm.real:.10f}")
    print(f"  neutral bilinear moments: {neutral_bilinear_moments(Nc, p)}")
    print(f"  bare-fermion lone crossing (census): 1 - p/Nc = "
          f"{bare_fermion_crossing(Nc, p):.6f} -> q^(1/p) = "
          f"{q ** (1.0 / p):.6f}")
    X1 = xhat_dirac(1, Nc, p, nmax=16)
    dev = np.abs((X1 + (1 - q ** (1.0 / p))
                  * transfer_formal(16, q))[:10, :10]).max()
    print(f"  Xhat_1 + (1 - q^(1/p)) T dev = {dev:.2e}  "
          f"(pair half of gate G3 closed; diagonal half open)")
    print(f"  M_n tower chord rules published: "
          f"{MN_TOWER_CHORD_RULES_PUBLISHED}; pair half derived in-repo: "
          f"{MN_TOWER_PAIR_CHANNEL_DERIVED} (see docstring)")
