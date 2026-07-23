#!/usr/bin/env python3
r"""
Baryon-sector first contact (road_map rung 23): ED data for the paper's
baryon claims in the one-flavor U(N)/Dirac-SYK dictionary.
=======================================================================

What the paper actually says about baryons -- quoted in full, because it is
short (arXiv:2607.05678, HologramsSM.tex):

  line 273:  "The SM has features closely resembling the real Standard
      Model: A U(1) x SU(N) gauge theory, aka QED x QCD; electric forces;
      quarks; mesons; baryons; Regge trajectories; confinement."

  line 871:  "In a similar way there are forces between heavy states such as
      Baryons in the 't Hooft model.  These forces may be identified as
      Gravitational.  They are too weak to affect mesons unless the number
      of mesons is of order N^{1/2}, but they affect Baryons."

That is the paper's ENTIRE baryon discussion.  Nothing is derived: the
"baryon mass ~ N" is the background large-N fact (Witten 1979; 't Hooft-model
baryons) the sentences lean on, and "gravitational-strength inter-baryon
forces" is asserted, not computed.  This module produces the first ED data on
the hologram side.

The baryon operator and an EXACT collapse (finding 1)
-----------------------------------------------------
In the one-flavor dictionary the quarks are the c^dag_i, and the baryon is
the epsilon-tensor contraction of all Nc of them:

    B^dag = c^dag_0 c^dag_1 ... c^dag_{Nc-1},

the unique maximal-charge operator (it fills every mode).  The charge sectors
Q = 0 and Q = Nc are ONE-DIMENSIONAL (|empty>, |full>), and every term of the
charge-conserving H carries p/2 >= 2 annihilators, so

    <empty|H|empty> = 0                                  (exact: c_B|empty>=0),
    <full |H|full > = (-1)^{k(k-1)/2} sum_A J_AA,   k = p/2   (exact, generic),

where the sign is a Jordan-Wigner lemma for this module's ascending-monomial
convention: <s|c^dag_A c_A|s> = (-1)^{k(k-1)/2} for A a subset of occ(s) (the
annihilation parities are counted on s, the re-creation parities on s minus A;
the difference is #{pairs in A} = k(k-1)/2).  Both statements are asserted
against dense ED in test_baryons.py.  Consequently the T = infinity baryon
two-point function collapses to a single undamped oscillation,

    G(t) = Tr(B(t) B^dag)/dim = exp(-i DeltaE t)/dim,
    DeltaE = <full|H|full> - <empty|H|empty> = (-1)^{k(k-1)/2} sum_A J_AA,

with |G(t)| = 1/dim identically.  DeltaE is a sum of C(Nc,k) i.i.d. real
Gaussians, Var(DeltaE) = C(Nc,k) sigma^2 with sigma^2 = ((k)!)^2 2^p/Nc^{p-1},
so the ONLY frequency the correlator carries is

    DeltaE_typ = sqrt(C(Nc,k)) sigma  ~  Nc^{-1/2}   (p = 4)  -->  0,

and in the c_symmetric (CP-sharp) ensemble, which has no diagonal couplings,
DeltaE = 0 EXACTLY.  So: the naive baryon two-point function at T_B = infinity
carries NO information about the extensive baryon mass -- its one frequency
VANISHES at large Nc instead of growing like Nc.  Worse, B^dag annihilates
every state of nonzero charge, so <psi|B(t)B^dag|psi> = 0 for ANY state
without an empty-sector component (including the interacting ground state,
which sits at Q ~ Nc/2): in the one-flavor model there is no baryon
spectroscopy to be had from two-point functions at all.  This exact collapse
is itself the first result of this module.

The physical mass proxy: sector energetics (finding 2, measured)
----------------------------------------------------------------
The baryon mass must therefore be read from CHARGE-SECTOR ENERGETICS.  With
E0(q) the ground-state energy of the charge-q sector (dimension C(Nc,q)),

    M_B(Nc) = E0(Q = Nc) - E0(global ground state)

is the cost of charging the equilibrium state up to the unique baryon
(= fully filled) state -- the 0+1d analog of the baryon rest mass.  Honest
framing: in this one-flavor model E0(Nc) is O(Nc^{-1/2}) (exactly 0 for
c_symmetric), so M_B is dominated by the DEPTH of the charging curve, i.e.
"baryon mass ~ N" here reduces to extensivity of the SYK ground-state energy.
The module measures M_B(Nc) with disorder statistics (mean +- SEM, ddof=1)
for Nc = 4..10 at p = 4 in both ensembles and fits a*Nc + b vs a*Nc^gamma vs
a*Nc (weighted chi^2; AIC = chi2 + 2K, same-noise comparison as in
self_averaging.py -- Delta-AIC < 2 means indistinguishable).

Baryon-baryon interaction proxy (finding 3, measured)
-----------------------------------------------------
There is no spatial separation in 0+1d, so no force in the literal sense.
The accessible analog is the CURVATURE of the charging curve at maximal
charge,

    Delta2(Nc-1) = E0(Nc) + E0(Nc-2) - 2 E0(Nc-1)
                 = [cost of the 2nd baryon-hole] - [cost of the 1st],

the interaction energy of two holes in the baryon.  The sharp testable
version of "inter-baryon forces are gravitational-strength" is an Nc-scaling
HIERARCHY: is Delta2 parametrically suppressed relative to a meson reference
scale?  The module measures Delta2(Nc-1) statistics and fits the power law
of |mean| and of the ratio |Delta2|/escale.  This is a charging-curve
statement, NOT a spatial force measurement -- stated wherever the numbers
are reported.

Denominator caveat (load-bearing for the headline exponent): the reference
scale used here, escale = sqrt(Tr H^2/dim), is the TOTAL many-body spectral
width -- it grows like sqrt(Nc) by the coupling normalization -- and NOT the
O(1) excitation energy of the M_n meson oscillations that this repo measures
in mn_spectroscopy.py.  The ratio exponent therefore depends on the
denominator choice: against the many-body width the asymptotic edge estimate
reads ratio ~ 1/Nc, while against an O(1) meson excitation energy the same
estimate reads |Delta2| ~ Nc^{-1/2}.  The QUALITATIVE conclusion -- the edge
curvature is parametrically suppressed at large Nc under either denominator
-- is robust to the choice; the specific power is not.

Exact edge structure (derived here, asserted in tests): every H term needs
B a subset of occ (k modes) and A inside the complement-plus-B, so blocks
with q < k vanish IDENTICALLY (free quarks below charge p/2); in the
c_symmetric ensemble (A, B disjoint) the mirror blocks q > Nc - k vanish too,
so E0(Nc) = E0(Nc-1) = 0 exactly and Delta2(Nc-1) = E0(Nc-2) = E0(2) < 0
(particle-hole symmetry E0(q) = E0(Nc-q) holds per instance): the c_symmetric
charging curve is exactly CONCAVE at the edge -- the second hole is cheaper
than the first (which is free).  The generic ensemble has no such identity;
its edge blocks are fed by the diagonal and near-diagonal couplings.

Sector second moment, closed form (cross-validation anchor)
-----------------------------------------------------------
For both ensembles the disorder mean of the sector-restricted second moment
is (Wick: only the pairing J_AB with J_BA = conj(J_AB) survives, and
T_AB T_BA = T_AB T_AB^dag is 0/1 on the diagonal),

    E[Tr_q H^2]/dim_q = sigma^2 sum_{(A,B) in table} C(Nc-|AuB|, q-k) / C(Nc,q),

with the sum over the SAME ordered-pair table draw_couplings uses (all pairs
incl. diagonal for "generic"; ordered disjoint pairs for "c_symmetric").
Grouping by overlap j = |A n B| gives the equivalent combinatorial form
sum_j n_pairs(j) C(Nc-2k+j, q-k) with n_pairs(j) = C(Nc,k)C(k,j)C(Nc-k,k-j)
(generic; c_symmetric keeps only j = 0).  Both routes are computed
independently and asserted equal, and the ensemble ED means are checked
against them within error bars.  Summing over q recovers the
coupling_variance_dirac normalization statement Tr H^2/dim ~ Nc.

Measured results (p = 4, Nc = 4..10, production run -> baryons.json)
--------------------------------------------------------------------
  *  M_B grows strongly and monotonically in both ensembles.  Weighted fits:
     generic:      best a*Nc^gamma with gamma = 0.795 +- 0.021 (Delta-AIC 4.0
                   over a*Nc + b with a = 0.73 +- 0.02, b = 1.19 +- 0.13);
     c_symmetric:  best a*Nc + b with a = 0.619 +- 0.004, b = -1.04 +- 0.03
                   (Delta-AIC 115 over the power law; gamma alone = 1.30
                   because of the negative offset).
     Honest verdict: the window is CONSISTENT with "baryon mass ~ N" (both
     ensembles are well described by a linear law with an O(1) offset; pure
     a*Nc is strongly rejected in-window by finite-size offsets) but Nc <= 10
     cannot pin the asymptotic exponent; the two ensembles' preferred
     effective forms differ at this size.  c_symmetric chi2/dof ~ 15 for the
     best model: with SEMs this small the fits RESOLVE finite-size curvature
     beyond any 2-parameter form.
  *  The hierarchy (effective in-window exponents): |Delta2(Nc-1)|/escale
     falls like Nc^(-1.79 +- 0.03) (generic: good fit, chi2/dof = 5.2/5;
     |Delta2| ~ Nc^(-1.31 +- 0.03), Delta2 > 0, convex edge) and like
     Nc^(-0.85 +- 0.04) (c_symmetric: the power-law form is itself REJECTED,
     chi2/dof = 194/5, error quoted after sqrt(chi2/dof) inflation;
     Delta2 < 0 exactly = E0(2), concave edge, |Delta2| saturating in-window
     with its own rejected fit Nc^(+0.34 +- 0.12), chi2/dof = 629/5).
     What these numbers DEMONSTRATE is the qualitative hierarchy only: the
     edge curvature falls relative to the many-body width, by many sigma, in
     both ensembles.  They do NOT confirm the specific asymptotic edge
     estimate |Delta2| ~ Nc^{-1/2}, escale ~ Nc^{+1/2} => ratio ~ 1/Nc:
     the generic in-window ratio exponent -1.79 +- 0.03 sits ~25 fit-sigma
     from -1, generic |Delta2| goes as Nc^{-1.31} not Nc^{-0.5}, and the
     c_symmetric exponents fail their own fits (pre-asymptotic window).
     The estimate is supported at the level of VALUES, not window exponents:
     c_symmetric Delta2 = E0(2) exactly, and the heuristic sparse-GOE edge
     -2 sqrt(C(Nc-2,2) sigma^2) (csym_edge_estimate) gives -2.68 at Nc = 10
     vs measured -2.52 +- 0.01, with the measured/estimate ratio rising
     0.87 -> 0.94 over the production window Nc = 6..10 (0.96 at Nc = 12 in
     an out-of-window check; approaching 1 from below) and the
     estimate's pre-asymptotic peak at Nc = 5 + sqrt(7) ~ 7.6 matching the
     observed in-window saturation of |Delta2|.  So ratio ~ 1/Nc is the
     plausible ASYMPTOTIC reading, consistent with but NOT established by
     this window.  In the only sense accessible in 0+1d the baryon-baryon
     term is parametrically suppressed relative to the many-body width --
     compatible with the "gravitational-strength" reading (see the
     denominator caveat above for how the power depends on the reference
     scale) -- with the sign (attractive vs repulsive analog)
     ensemble-DEPENDENT and therefore not a robust prediction of the
     construction.
  *  Scale check: generic escale ~ Nc^(0.48 +- 0.01) (the designed sqrt(Nc);
     chi2/dof = 7.6/5); c_symmetric escale shows Nc^(1.16 +- 0.08) in-window
     (power law REJECTED, chi2/dof = 589/5, error sqrt(chi2/dof)-inflated)
     because the disjointness constraint discards an O(p^2/Nc) coupling
     fraction that is order one at these sizes -- ALL c_symmetric window
     exponents are effective, not asymptotic.
  *  The collapse frequency: measured std(DeltaE) falls from 2.36 (Nc = 4)
     to 1.83 (Nc = 10), tracking sqrt(C(Nc,2)) sigma exactly, while M_B grows
     4.0 -> 8.4: the naive correlator's only frequency and the physical mass
     proxy move in OPPOSITE directions with Nc.

What is exact vs measured vs conjectured
----------------------------------------
  EXACT (asserted): sector reassembly = dense H; the correlator collapse;
    DeltaE closed form; DeltaE = 0 (c_symmetric); vanishing edge blocks;
    Delta2 = E0(2) (c_symmetric); PH palindromy; sector-m2 closed form
    (two routes).
  MEASURED (with error bars): M_B(Nc) statistics and its extensivity fits;
    Delta2 statistics and scaling; the Delta2/escale hierarchy.  Where a
    power-law fit is statistically rejected (chi2 > 5 dof) the exponent is
    reported as an effective in-window slope with sqrt(chi2/dof)-inflated
    errors, and flagged as rejected in the JSON, the prints, and the figure.
  HEURISTIC (checked numerically, not proved): the sparse-GOE semicircle
    edge estimate csym_edge_estimate for the c_symmetric E0(2); the
    ratio ~ 1/Nc asymptotic reading of the hierarchy built on it.
  CONJECTURED / not tested: any identification of Delta2 with a bulk
    gravitational force; the O(sqrt N)-meson mechanism (needs multi-meson
    states, out of ED reach here); multi-baryon physics (needs >= 2 flavors,
    out of scope).

Conventions: mode i <-> bit (Nc-1-i) (dirac.py); H from dirac_hamiltonian;
d/dt X = i[H, X]; ensembles "generic" and "c_symmetric" as in draw_couplings.
Run from inside syk-self-averaging/.  __main__ writes baryons.json and
baryons.png in ~2-4 min.
"""

from __future__ import annotations

import json
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

from dirac import (apply_monomial, coupling_variance_dirac, dirac_hamiltonian,
                   draw_couplings)
from self_averaging import fit_verdict


# --------------------------------------------------------------------------
# charge-sector machinery
# --------------------------------------------------------------------------
def sector_states(Nc: int, q: int) -> np.ndarray:
    """Basis indices of the charge-q sector (popcount = q), ascending."""
    states = np.arange(1 << Nc, dtype=np.int64)
    return states[np.bitwise_count(states) == q]


def sector_hamiltonian(Nc: int, couplings: dict, q: int) -> np.ndarray:
    """H restricted to the charge-q sector, assembled DIRECTLY on the sector
    basis (dimension C(Nc,q)) from the coupling table -- shares only the
    per-monomial physics (apply_monomial) with dirac_hamiltonian, not the
    dense full-space assembly.  Raises if any term leaves the sector, which
    doubles as a structural charge-conservation check."""
    st = sector_states(Nc, q)
    d = len(st)
    pos = np.full(1 << Nc, -1, dtype=np.int64)
    pos[st] = np.arange(d)
    Hq = np.zeros((d, d), dtype=complex)
    cols = np.arange(d)
    for (A, B), J in couplings.items():
        valid, rows, signs = apply_monomial(Nc, A, B, st)
        r = pos[rows[valid]]
        if r.size and r.min() < 0:
            raise AssertionError("term left the charge sector: H not "
                                 "charge-conserving?")
        Hq[r, cols[valid]] += J * signs[valid]
    return Hq


def reassemble_dense(Nc: int, blocks: list) -> np.ndarray:
    """Scatter the Nc+1 sector blocks back into the full 2^Nc space (for the
    machine-exactness test against dirac_hamiltonian)."""
    dim = 1 << Nc
    H = np.zeros((dim, dim), dtype=complex)
    for q, blk in enumerate(blocks):
        st = sector_states(Nc, q)
        H[np.ix_(st, st)] = blk
    return H


def dense_sector_blocks(H: np.ndarray, Nc: int) -> list:
    """Slice a dense charge-conserving H into its Nc+1 sector blocks (the
    fast production route; sector_hamiltonian is the independent one)."""
    states = np.arange(1 << Nc, dtype=np.int64)
    qs = np.bitwise_count(states)
    return [H[np.ix_(np.flatnonzero(qs == q), np.flatnonzero(qs == q))]
            for q in range(Nc + 1)]


# --------------------------------------------------------------------------
# the baryon operator and its exact two-point collapse
# --------------------------------------------------------------------------
def baryon_dagger(Nc: int) -> np.ndarray:
    """B^dag = c^dag_0 c^dag_1 ... c^dag_{Nc-1} as a dense matrix.  With the
    module's JW convention every creation parity counts only lower-index
    modes, all empty at the moment of application, so B^dag|empty> = +|full>
    with sign exactly +1 (asserted in the tests)."""
    dim = 1 << Nc
    states = np.arange(dim, dtype=np.int64)
    valid, rows, signs = apply_monomial(Nc, tuple(range(Nc)), (), states)
    B = np.zeros((dim, dim), dtype=complex)
    B[rows[valid], states[valid]] = signs[valid]
    return B


def baryon_sign(k: int) -> int:
    """(-1)^{k(k-1)/2}: the JW sign of <s|c^dag_A c_A|s> for |A| = k."""
    return -1 if (k * (k - 1) // 2) % 2 else 1


def deltaE_from_couplings(couplings: dict, k: int) -> float:
    """DeltaE = <full|H|full> - <empty|H|empty> = (-1)^{k(k-1)/2} sum_A J_AA
    (diagonal couplings are real by construction; 0 for c_symmetric, which
    has no diagonal couplings)."""
    tot = sum(J.real if isinstance(J, complex) else J
              for (A, B), J in couplings.items() if A == B)
    return baryon_sign(k) * float(tot)


def baryon_correlator_ed(H: np.ndarray, ts: np.ndarray) -> np.ndarray:
    """G(t) = Tr(B(t) B^dag)/dim by dense ED, using NOTHING about the
    collapse: eigendecompose H, rotate B and B^dag, phase-sum.  The test
    asserts this equals exp(-i DeltaE t)/dim to machine precision."""
    dim = H.shape[0]
    Nc = int(np.log2(dim))
    Bd = baryon_dagger(Nc)
    Bop = Bd.conj().T
    E, V = np.linalg.eigh(H)
    Bt = V.conj().T @ Bop @ V         # B in the eigenbasis
    Bdt = V.conj().T @ Bd @ V
    G = np.empty(len(ts), dtype=complex)
    for i, t in enumerate(ts):
        ph = np.exp(1j * E * t)
        # Tr(e^{iHt} B e^{-iHt} B^dag) = sum_{mn} ph_m Bt[m,n] conj(ph_n) Bdt[n,m]
        G[i] = np.einsum("m,mn,n,nm->", ph, Bt, ph.conj(), Bdt) / dim
    return G


# --------------------------------------------------------------------------
# sector second moment: closed form, two independent routes
# --------------------------------------------------------------------------
def sector_m2_enumerate(Nc: int, p: int, q: int, ensemble: str) -> float:
    """E[Tr_q H^2]/dim_q by direct enumeration of the SAME ordered-pair table
    draw_couplings builds (structure only; no random numbers)."""
    k = p // 2
    sigma2 = coupling_variance_dirac(Nc, p)
    subsets = list(combinations(range(Nc), k))
    tot = 0
    for A in subsets:
        for B in subsets:
            if ensemble == "c_symmetric" and (set(A) & set(B)):
                continue
            u = len(set(A) | set(B))
            tot += comb(Nc - u, q - k) if q - k >= 0 else 0
    return sigma2 * tot / comb(Nc, q)


def sector_m2_closed_form(Nc: int, p: int, q: int, ensemble: str) -> float:
    """Same quantity via the overlap-resolved combinatorial closed form
    sum_j C(Nc,k)C(k,j)C(Nc-k,k-j) C(Nc-2k+j, q-k)  (j = 0 only for
    c_symmetric) -- shares no code with sector_m2_enumerate."""
    k = p // 2
    sigma2 = coupling_variance_dirac(Nc, p)
    if q < k:
        return 0.0
    js = [0] if ensemble == "c_symmetric" else range(k + 1)
    tot = 0
    for j in js:
        npairs = comb(Nc, k) * comb(k, j) * comb(Nc - k, k - j)
        tot += npairs * comb(Nc - 2 * k + j, q - k)
    return sigma2 * tot / comb(Nc, q)


def csym_edge_estimate(Nc: int, p: int) -> float:
    """HEURISTIC sparse-GOE estimate of the c_symmetric edge curvature,
    which equals E0(k) exactly (k = p/2; for p = 4 that is Delta2(Nc-1) =
    E0(2)).  The q = k sector Hamiltonian is a C(Nc,k) x C(Nc,k) real
    symmetric matrix with zero diagonal whose row (state S, |S| = k) has one
    independent N(0, sigma^2) entry for each k-subset A disjoint from S,
    i.e. row variance R = C(Nc-k, k) sigma^2; the semicircle edge is then

        E0(k) ~ -2 sqrt(R) = -2 sqrt(C(Nc-k, k) sigma^2).

    This is an asymptotic heuristic (finite degree, sparse), NOT exact: the
    measured mean E0(2) approaches it from below, ratio 0.87 -> 0.94 over
    the production window Nc = 6..10 and 0.96 at Nc = 12
    (test_baryons.py checks the value at Nc = 8).  For p = 4 the estimate
    peaks at Nc = 5 + sqrt(7) ~ 7.65 before its asymptotic Nc^{-1/2} decay
    -- matching the observed in-window saturation of |Delta2| -- because
    sigma^2 = 4 * 2^4 / Nc^3 while C(Nc-2, 2) ~ Nc^2 / 2."""
    k = p // 2
    return -2.0 * float(np.sqrt(comb(Nc - k, k)
                                * coupling_variance_dirac(Nc, p)))


# --------------------------------------------------------------------------
# charging curve, per-instance observables, ensemble statistics
# --------------------------------------------------------------------------
def instance_charging_curve(H: np.ndarray, Nc: int):
    """E0(q) for q = 0..Nc and the per-sector second moments Tr_q H^2/dim_q,
    from the dense charge-conserving H."""
    blocks = dense_sector_blocks(H, Nc)
    E0 = np.empty(Nc + 1)
    m2q = np.empty(Nc + 1)
    for q, blk in enumerate(blocks):
        E0[q] = blk[0, 0].real if blk.shape[0] == 1 else \
            np.linalg.eigvalsh(blk)[0]
        m2q[q] = float(np.linalg.norm(blk) ** 2) / blk.shape[0]
    return E0, m2q


def instance_observables(Nc: int, p: int, seed: int, ensemble: str) -> dict:
    """One disorder instance: charging curve E0(q), baryon mass proxy
    M_B = E0(Nc) - min_q E0(q), edge curvature Delta2(Nc-1), many-body
    spectral width escale = sqrt(Tr H^2/dim) (NOT the O(1) meson excitation
    scale -- see the denominator caveat in the module docstring), and the
    collapse frequency DeltaE."""
    rng = np.random.default_rng(seed)
    couplings = draw_couplings(Nc, p, rng, ensemble)
    H = dirac_hamiltonian(Nc, p, couplings=couplings)
    E0, m2q = instance_charging_curve(H, Nc)
    qstar = int(np.argmin(E0))
    dims = np.array([comb(Nc, q) for q in range(Nc + 1)], dtype=float)
    escale = float(np.sqrt(np.sum(m2q * dims) / (1 << Nc)))
    return dict(
        E0=E0, m2q=m2q, qstar=qstar,
        MB=float(E0[Nc] - E0[qstar]),
        delta2=float(E0[Nc] + E0[Nc - 2] - 2.0 * E0[Nc - 1]),
        escale=escale,
        deltaE=float(H[-1, -1].real - H[0, 0].real),
        deltaE_closed=deltaE_from_couplings(couplings, p // 2),
    )


def _mean_sem(a: np.ndarray):
    a = np.asarray(a, dtype=float)
    m = float(a.mean())
    s = float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else 0.0
    return m, s


def ensemble_run(Nc: int, p: int, n_inst: int, ensemble: str,
                 seed0: int = None) -> dict:
    """Disorder statistics (mean, SEM with ddof=1) of the per-instance
    observables.  Seeds are seed0 + j, deterministic per (Nc, ensemble)."""
    if seed0 is None:
        seed0 = 10000 * Nc + (500000 if ensemble == "c_symmetric" else 0)
    E0s, MBs, d2s, escs, qss, dEs, ratios = [], [], [], [], [], [], []
    for j in range(n_inst):
        d = instance_observables(Nc, p, seed0 + j, ensemble)
        E0s.append(d["E0"])
        MBs.append(d["MB"])
        d2s.append(d["delta2"])
        escs.append(d["escale"])
        qss.append(d["qstar"])
        dEs.append(d["deltaE"])
        ratios.append(d["delta2"] / d["escale"])
    E0s = np.array(E0s)
    mb_m, mb_s = _mean_sem(MBs)
    d2_m, d2_s = _mean_sem(d2s)
    ra_m, ra_s = _mean_sem(ratios)
    es_m, es_s = _mean_sem(escs)
    qs_m, qs_s = _mean_sem(qss)
    return dict(
        Nc=Nc, p=p, ensemble=ensemble, n_inst=n_inst, seed0=seed0,
        E0_mean=[float(x) for x in E0s.mean(axis=0)],
        E0_sem=[float(x) for x in E0s.std(axis=0, ddof=1) / np.sqrt(n_inst)],
        MB_mean=mb_m, MB_sem=mb_s, MB_std=float(np.std(MBs, ddof=1)),
        delta2_mean=d2_m, delta2_sem=d2_s,
        ratio_mean=ra_m, ratio_sem=ra_s,
        escale_mean=es_m, escale_sem=es_s,
        qstar_mean=qs_m, qstar_sem=qs_s,
        deltaE_std=float(np.std(dEs, ddof=1)),
        deltaE_std_exact=float(
            np.sqrt(comb(Nc, p // 2) * coupling_variance_dirac(Nc, p))
            if ensemble == "generic" else 0.0),
    )


# --------------------------------------------------------------------------
# fits: extensivity of M_B and the Delta2 scaling hierarchy
# --------------------------------------------------------------------------
def fit_mass_models(ns, ys, errs) -> dict:
    """Weighted least-squares fits of the M_B(Nc) means:  a*Nc + b  vs
    a*Nc^gamma  vs  a*Nc.  With fixed per-point noise the AIC comparison is
    AIC = chi2_min + 2K (additive constants common to all models dropped) --
    self_averaging.fit_verdict then applies unchanged.  Delta-AIC < 2 means
    the window cannot distinguish the models; say so, do not oversell."""
    ns = np.asarray(ns, float)
    ys = np.asarray(ys, float)
    errs = np.asarray(errs, float)
    models = {
        "a*Nc + b": (lambda n, a, b: a * n + b, (1.0, 0.0)),
        "a*Nc^gamma": (lambda n, a, g: a * n ** g, (ys[-1] / ns[-1], 1.0)),
        "a*Nc": (lambda n, a: a * n, (ys[-1] / ns[-1],)),
    }
    out = {}
    for name, (f, p0) in models.items():
        popt, pcov = curve_fit(f, ns, ys, p0=p0, sigma=errs,
                               absolute_sigma=True, maxfev=20000)
        perr = np.sqrt(np.diag(pcov))
        chi2 = float(np.sum(((ys - f(ns, *popt)) / errs) ** 2))
        K = len(popt)
        out[name] = dict(params=[float(x) for x in popt],
                         param_errs=[float(x) for x in perr],
                         chi2=chi2, dof=len(ns) - K, K=K,
                         aic=chi2 + 2 * K)
    return out


def power_exponent_fit(ns, ys, errs) -> dict:
    """Power-law a*Nc^e fit (weighted); used for the Delta2 / escale / ratio
    scalings.  ys must be positive -- callers pass |mean| and report the sign
    separately.

    Honesty bookkeeping: the raw covariance errors (a_err, e_err) assume the
    power-law form is CORRECT.  When the fit is statistically rejected
    (chi2 >> dof) they understate the model uncertainty by ~sqrt(chi2/dof),
    so the dict also carries a_err_scaled / e_err_scaled -- the standard
    sqrt(chi2/dof) inflation, floored at 1 -- and a boolean 'rejected'
    (chi2 > 5 dof, the same criterion the figure uses).  Every consumer that
    QUOTES an exponent must quote the scaled error and flag rejection; a
    rejected fit's exponent is an effective in-window slope, not a law."""
    ns = np.asarray(ns, float)
    ys = np.asarray(ys, float)
    errs = np.asarray(errs, float)
    f = lambda n, a, e: a * n ** e
    popt, pcov = curve_fit(f, ns, ys, p0=(ys[-1], -1.0), sigma=errs,
                           absolute_sigma=True, maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    chi2 = float(np.sum(((ys - f(ns, *popt)) / errs) ** 2))
    dof = len(ns) - 2
    scale = max(1.0, float(np.sqrt(chi2 / dof))) if dof > 0 else 1.0
    return dict(a=float(popt[0]), e=float(popt[1]),
                a_err=float(perr[0]), e_err=float(perr[1]),
                a_err_scaled=float(perr[0]) * scale,
                e_err_scaled=float(perr[1]) * scale,
                chi2=chi2, dof=dof, rejected=bool(chi2 > 5.0 * dof))


# --------------------------------------------------------------------------
# production run
# --------------------------------------------------------------------------
N_INST = {4: 400, 5: 400, 6: 400, 7: 300, 8: 250, 9: 150, 10: 100}
COLORS = {"generic": "#2166ac", "c_symmetric": "#b2182b"}   # validated pair


def make_figure(rows: dict, fits: dict, scal: dict, path: Path, p: int):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15, 4.6))
    ncs = sorted({r["Nc"] for rr in rows.values() for r in rr})
    nmax = max(ncs)

    # ---- (a) mean charging curves at the largest Nc ----------------------
    for ens, rr in rows.items():
        r = [x for x in rr if x["Nc"] == nmax][0]
        qs = np.arange(nmax + 1)
        axA.errorbar(qs, r["E0_mean"], yerr=r["E0_sem"], fmt="o-", ms=4,
                     lw=1.5, capsize=2, color=COLORS[ens], label=ens)
    rg = [x for x in rows["generic"] if x["Nc"] == nmax][0]
    qstar = int(np.argmin(rg["E0_mean"]))
    axA.annotate("", xy=(nmax, rg["E0_mean"][nmax]),
                 xytext=(nmax, rg["E0_mean"][qstar]),
                 arrowprops=dict(arrowstyle="<->", color="#555555", lw=1.2))
    axA.text(nmax - 0.4, 0.5 * (rg["E0_mean"][qstar] + rg["E0_mean"][nmax]),
             r"$M_B$", ha="right", va="center", fontsize=11, color="#333333")
    axA.axhline(0.0, color="#999999", lw=0.8, alpha=0.5)
    axA.set_xlabel("charge sector $q$")
    axA.set_ylabel(r"mean $E_0(q)$")
    axA.set_title(f"Charging curve, $N_c={nmax}$, $p={p}$")
    axA.legend(frameon=False, fontsize=9)
    axA.grid(True, alpha=0.2)

    # ---- (b) M_B vs Nc with fits -----------------------------------------
    for ens, rr in rows.items():
        ns = np.array([r["Nc"] for r in rr], float)
        ys = np.array([r["MB_mean"] for r in rr])
        es = np.array([r["MB_sem"] for r in rr])
        axB.errorbar(ns, ys, yerr=es, fmt="o", ms=5, capsize=2,
                     color=COLORS[ens], label=f"{ens} (ED)")
        lin = fits[ens]["a*Nc + b"]["params"]
        pw = fits[ens]["a*Nc^gamma"]["params"]
        nn = np.linspace(ns[0], ns[-1], 100)
        axB.plot(nn, lin[0] * nn + lin[1], "-", lw=1.2, color=COLORS[ens],
                 alpha=0.7)
        axB.plot(nn, pw[0] * nn ** pw[1], "--", lw=1.2, color=COLORS[ens],
                 alpha=0.7)
        ge = fits[ens]["a*Nc^gamma"]
        # honesty: sqrt(chi2/dof)-inflate the quoted gamma error and flag
        # statistically rejected fits (same policy as power_exponent_fit)
        gsc = max(1.0, np.sqrt(ge["chi2"] / ge["dof"]))
        gbad = ge["chi2"] > 5 * ge["dof"]
        axB.text(0.03, 0.95 - 0.09 * (ens == "c_symmetric"),
                 rf"{ens}: $\gamma = {ge['params'][1]:.3f} \pm "
                 rf"{ge['param_errs'][1] * gsc:.3f}$"
                 + (rf"  (REJECTED, $\chi^2$/dof = {ge['chi2']:.0f}/"
                    rf"{ge['dof']})" if gbad else ""),
                 transform=axB.transAxes, fontsize=9, va="top",
                 color=COLORS[ens])
    axB.set_xlabel("$N_c$")
    axB.set_ylabel(r"$M_B = E_0(N_c) - E_0(\mathrm{ground})$")
    axB.set_title(r"Baryon mass proxy (solid $aN_c\!+\!b$, dashed "
                  r"$aN_c^\gamma$)")
    axB.legend(frameon=False, fontsize=9, loc="lower right")
    axB.grid(True, alpha=0.2)

    # ---- (c) the hierarchy: |Delta2| vs meson scale (log-log) ------------
    for ens, rr in rows.items():
        ns = np.array([r["Nc"] for r in rr], float)
        d2 = np.array([abs(r["delta2_mean"]) for r in rr])
        d2e = np.array([r["delta2_sem"] for r in rr])
        d2m = np.array([r["delta2_mean"] for r in rr])
        sgn = (r"$\Delta_2>0$" if (d2m > 0).all() else
               r"$\Delta_2<0$" if (d2m < 0).all() else "mixed sign")
        axC.errorbar(ns, d2, yerr=d2e, fmt="o", ms=5, capsize=2,
                     color=COLORS[ens],
                     label=rf"$|\Delta_2(N_c\!-\!1)|$ {ens} ({sgn})")
        s = scal[ens]["delta2"]
        nn = np.linspace(ns[0], ns[-1], 50)
        axC.plot(nn, s["a"] * nn ** s["e"], "-", lw=1.1, color=COLORS[ens],
                 alpha=0.6)
        # honesty: scaled errors everywhere; flag statistically rejected fits
        axC.text(0.03, 0.16 - 0.09 * (ens == "c_symmetric"),
                 rf"{ens}: $|\Delta_2| \sim N_c^{{{s['e']:.2f} \pm "
                 rf"{s['e_err_scaled']:.2f}}}$"
                 + (rf"  (REJECTED, $\chi^2$/dof = {s['chi2']:.0f}/"
                    rf"{s['dof']})" if s["rejected"] else ""),
                 transform=axC.transAxes, fontsize=9, color=COLORS[ens])
    ns = np.array([r["Nc"] for r in rows["generic"]], float)
    esc = np.array([r["escale_mean"] for r in rows["generic"]])
    esce = np.array([r["escale_sem"] for r in rows["generic"]])
    axC.errorbar(ns, esc, yerr=esce, fmt="s", ms=5, capsize=2,
                 color="#777777",
                 label=r"many-body width $\sqrt{Tr H^2/dim}$ (generic)")
    se = scal["generic"]["escale"]
    axC.text(0.03, 0.25, rf"width $\sim N_c^{{{se['e']:.2f} \pm "
             rf"{se['e_err_scaled']:.2f}}}$"
             + (rf"  (REJECTED, $\chi^2$/dof = {se['chi2']:.0f}/"
                rf"{se['dof']})" if se["rejected"] else ""),
             transform=axC.transAxes, fontsize=9, color="#777777")
    axC.set_xscale("log")
    axC.set_yscale("log")
    axC.set_xlabel("$N_c$")
    axC.set_ylabel("energy")
    axC.set_title("Edge curvature vs many-body spectral width")
    axC.legend(frameon=False, fontsize=8, loc="upper right")
    axC.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    p = 4
    ncs = sorted(N_INST)
    rows = {"generic": [], "c_symmetric": []}

    print("=" * 74)
    print("Rung 23: baryon sector, p = 4, Nc = 4..10, both ensembles")
    print("=" * 74)

    print("\n1. The exact collapse of the naive baryon two-point function")
    print("-" * 74)
    Nc_demo = 6
    rng = np.random.default_rng(2607)
    coup = draw_couplings(Nc_demo, p, rng, "generic")
    H = dirac_hamiltonian(Nc_demo, p, couplings=coup)
    ts = np.linspace(0.0, 8.0, 33)
    G = baryon_correlator_ed(H, ts)
    dE = deltaE_from_couplings(coup, p // 2)
    pred = np.exp(-1j * dE * ts) / (1 << Nc_demo)
    dev = float(np.abs(G - pred).max())
    print(f"  Nc={Nc_demo} generic: max|G(t) - exp(-i dE t)/dim| = {dev:.2e}"
          f"   (dE = {dE:+.4f})")
    print(f"  |G(t)|*dim - 1: max = {float(np.abs(np.abs(G)*(1<<Nc_demo)-1).max()):.2e}"
          "   (NO decay, a single undamped oscillation)")
    print("  -> Tr(B(t)B^dag)/dim = exp(-i DeltaE t)/dim EXACTLY;"
          " the T=inf baryon")
    print("     2-pt function carries NO information about the extensive"
          " baryon mass.")
    collapse = dict(Nc=Nc_demo, p=p, max_dev=dev, deltaE=dE)

    print("\n2. Charging curves, M_B, Delta2 (disorder mean +- SEM, ddof=1)")
    print("-" * 74)
    print(f"{'ens':>12} {'Nc':>3} {'n':>4} {'M_B':>16} {'q*':>5} "
          f"{'Delta2(Nc-1)':>18} {'escale':>14} {'dE_std':>7}")
    for ens in rows:
        for Nc in ncs:
            r = ensemble_run(Nc, p, N_INST[Nc], ens)
            rows[ens].append(r)
            print(f"{ens:>12} {Nc:>3} {r['n_inst']:>4} "
                  f"{r['MB_mean']:>9.4f} ± {r['MB_sem']:<5.4f}"
                  f"{r['qstar_mean']:>5.2f} "
                  f"{r['delta2_mean']:>10.4f} ± {r['delta2_sem']:<6.4f}"
                  f"{r['escale_mean']:>8.4f} ± {r['escale_sem']:<5.4f}"
                  f"{r['deltaE_std']:>7.3f}")
        print()

    print("3. Extensivity fits of M_B(Nc)  (the 'baryon mass ~ N' test)")
    print("-" * 74)
    fits = {}
    for ens in rows:
        ns = [r["Nc"] for r in rows[ens]]
        ys = [r["MB_mean"] for r in rows[ens]]
        es = [r["MB_sem"] for r in rows[ens]]
        f = fit_mass_models(ns, ys, es)
        fits[ens] = f
        print(f"  {ens}:")
        for name, d in f.items():
            ps = ", ".join(f"{v:.4f}±{e:.4f}"
                           for v, e in zip(d["params"], d["param_errs"]))
            print(f"    {name:12s}: [{ps}]  chi2/dof = "
                  f"{d['chi2']:.1f}/{d['dof']}  AIC = {d['aic']:7.1f}")
        print(f"    -> {fit_verdict(f)}")

    print("\n4. The hierarchy: Delta2 vs the many-body spectral width")
    print("-" * 74)
    print("  (errors on exponents are sqrt(chi2/dof)-inflated where the fit")
    print("   fails; REJECTED = chi2 > 5 dof: exponent is an effective")
    print("   in-window slope, not a law)")
    scal = {}
    for ens in rows:
        ns = [r["Nc"] for r in rows[ens]]
        d2m = np.array([r["delta2_mean"] for r in rows[ens]])
        d2e = np.array([r["delta2_sem"] for r in rows[ens]])
        ram = np.array([r["ratio_mean"] for r in rows[ens]])
        rae = np.array([r["ratio_sem"] for r in rows[ens]])
        sgn = "negative" if (d2m < 0).all() else \
              ("positive" if (d2m > 0).all() else "MIXED SIGN")
        insignif = int(np.sum(np.abs(d2m) < 2 * d2e))
        scal[ens] = dict(
            delta2=power_exponent_fit(ns, np.abs(d2m), d2e),
            ratio=power_exponent_fit(ns, np.abs(ram), rae),
            escale=power_exponent_fit(
                ns, [r["escale_mean"] for r in rows[ens]],
                [r["escale_sem"] for r in rows[ens]]),
            sign=sgn, n_insignificant=insignif)
        s = scal[ens]
        print(f"  {ens}: Delta2 mean is {sgn}"
              + (f" ({insignif} pts < 2 sigma from 0)" if insignif else ""))
        for label, key in (("|Delta2|       ", "delta2"),
                           ("escale (width) ", "escale"),
                           ("|Delta2|/escale", "ratio")):
            d = s[key]
            print(f"    {label} ~ Nc^({d['e']:+.2f} ± "
                  f"{d['e_err_scaled']:.2f})   chi2/dof = "
                  f"{d['chi2']:.1f}/{d['dof']}"
                  + ("   REJECTED (in-window slope only)"
                     if d["rejected"] else ""))
    print("\n  Reading (charging-curve proxy, NOT a spatial force): the edge")
    print("  curvature falls relative to the many-body width by many sigma in")
    print("  both ensembles -- that QUALITATIVE hierarchy is the demonstrated")
    print("  result.  The specific asymptotic estimate ratio ~ 1/Nc is NOT")
    print("  confirmed in-window: the generic ratio exponent "
          f"({scal['generic']['ratio']['e']:+.2f}) sits "
          f"{abs(scal['generic']['ratio']['e'] + 1.0) / scal['generic']['ratio']['e_err_scaled']:.0f} "
          "fit-sigma from -1,")
    print("  and the c_symmetric fits are rejected (pre-asymptotic window).")

    print("\n5. The c_symmetric edge VALUE vs the heuristic sparse-GOE estimate")
    print("-" * 74)
    print("  Delta2 = E0(2) exactly; estimate = -2 sqrt(C(Nc-2,2) sigma^2);")
    print(f"  estimate peaks at Nc = 5 + sqrt(7) = {5 + np.sqrt(7):.2f} "
          "(the in-window saturation).")
    edge_rows = []
    for r in rows["c_symmetric"]:
        est = csym_edge_estimate(r["Nc"], p)
        edge_rows.append(dict(Nc=r["Nc"], estimate=est,
                              measured=r["delta2_mean"],
                              measured_sem=r["delta2_sem"],
                              ratio=r["delta2_mean"] / est))
        print(f"    Nc={r['Nc']:>2}: estimate {est:+.4f}   measured "
              f"{r['delta2_mean']:+.4f} ± {r['delta2_sem']:.4f}   "
              f"ratio {r['delta2_mean'] / est:.3f}")
    print("  -> values consistent (ratio -> 1 from below); the window's")
    print("     effective exponents are NOT (see section 4): ratio ~ 1/Nc is")
    print("     an asymptotic reading, unestablished at Nc <= 10.")

    out = dict(meta=dict(p=p, n_inst=N_INST, ensembles=list(rows),
                         description="rung 23 baryon-sector ED, see "
                                     "baryons.py docstring"),
               collapse_demo=collapse, rows=rows, fits_MB=fits,
               scalings=scal,
               csym_edge_estimate=dict(
                   formula="-2*sqrt(C(Nc-2,2)*sigma^2), heuristic "
                           "sparse-GOE edge; peak at Nc = 5 + sqrt(7)",
                   peak_Nc=float(5 + np.sqrt(7)), rows=edge_rows))
    here = Path(__file__).resolve().parent
    with open(here / "baryons.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {here / 'baryons.json'}")
    make_figure(rows, fits, scal, here / "baryons.png", p)
    print(f"wrote {here / 'baryons.png'}")
