#!/usr/bin/env python3
r"""
The Dirac/QED campaign: U(N) self-averaging, charge-sector energetics, and
charged-vs-singlet spectral weight (road_map rung 19).
===========================================================================

Why this module exists.  Every executed self-averaging numeric in this repo is
Majorana/O(N), while the headline claims of Miyashita-Sekino-Susskind
(arXiv:2607.05678) are Dirac/U(N): "SM dual of DSSYK_infinity is 't Hooft
model, QCD(2) + QED", with the SYK fermions as electrically charged quarks.
The paper's U(1) story is anchored in one citation, at the M_0 (n = 0
"photon") discussion (HologramsSM.tex line 918):

    "M_0 is the total conserved U(1) charge.  It has the quantum numbers of
     the time component of the photon but because it is conserved all its
     matrix elements are time independent.  There is no propagating photon.
     But photon exchange gives rise to long range confining Coulomb forces.
     (See [Susskind:2020fiu].)  This parallels the 't Hooft model."

[Susskind:2020fiu] = L. Susskind, "Electric Forces in the Charged SYK Model",
JHAP 5, no.4 (2025), arXiv:2012.12326 -- cited exactly once, for exactly this
claim.  The sharp 0+1d observables that claim points at are (i) the ENERGETICS
of the U(1) charge sectors (what does it cost to charge the system?) and (ii)
the frequency scales of charged vs neutral (singlet) operators.  This module
measures both, and first nails the U(N) ensembles' self-averaging statistics
in closed form so every measured number sits on an exact anchor.

Model and ensembles (dirac.py conventions throughout):

    H = sum_{|A|=|B|=k} J_AB c^dag_A c_B,   k = p/2,   J_BA = conj(J_AB),
    sigma^2 = Var|J_AB| = coupling_variance_dirac(Nc, p) = (k!)^2 2^p/Nc^(p-1)

with the stored coupling table containing BOTH ordered keys (A,B) and (B,A):
"generic" = all pairs, complex off-diagonal (E z^2 = 0, E|z|^2 = sigma^2),
real diagonal J_AA; "c_symmetric" = disjoint pairs only, real symmetric --
the CP-sharp sub-ensemble.

--------------------------------------------------------------------------
1. Exact disorder statistics of m2 = Tr(H^2)/dim  (the U(N) analog of the
   Majorana relVar(m2) = 2/C(N,p))
--------------------------------------------------------------------------
Write m2 as the quadratic form

    m2 = sum_{alpha,beta} J_alpha J_beta T_{alpha beta},
    T_{alpha beta} = Tr(O_alpha O_beta)/dim,   O_(A,B) = c^dag_A c_B,

over stored keys.  Occupation counting gives every needed trace:

  * T_{alpha, alphabar} = Tr(O O^dag)/dim = 2^{-|A u B|}   (alphabar = (B,A)):
    O^dag|s> is nonzero iff s contains A and avoids B\A, i.e. |A u B| modes
    are forced and the rest are free -- hence 2^{Nc-|AuB|}/2^{Nc}.
  * General T_{alpha beta}, alpha = (A,B), beta = (C,D): nonzero iff for
    EVERY mode i the membership pattern (a,b,c,d) in {A,B,C,D} satisfies
    a + c = b + d (the monomial is diagonal), each admissible mode pattern
    being one of  (0000), (1100), (1001), (0110), (0011), (1111).  Every mode
    is then acted an even number of times, so the Jordan-Wigner sign is the
    same for all valid states and |T| = 2^{-u}, u = |A u B u C u D|.

STRUCTURAL FACT (verified by enumeration in the tests, and a correction to
the naive picture): the instance-level identity

    m2 = sum_stored |J_AB|^2 T_AB           (T_AB := 2^{-|AuB|})

holds EXACTLY for the c_symmetric ensemble -- disjointness (A n B = C n D =
empty) kills every admissible pattern except C = B, D = A, so T is the pure
pairing delta -- but FAILS for the generic ensemble, where overlapping keys
admit cross-pairings ((C,D) != (B,A) with T != 0).  For generic instances
only the full quadratic form reproduces Tr(H^2)/dim (machine-checked).

Disorder statistics by Wick/Isserlis on the jointly-Gaussian J (this
automatically carries the real-vs-complex fourth moments E|z|^4 = 2 sigma^4,
E x^4 = 3 sigma^4, and the (A,B)/(B,A) conjugate pairing):

    E[m2]  = sum_{alpha beta} C_{alpha beta} T_{alpha beta},
    Var m2 = 2 tr[(C T)^2],          C_{alpha beta} = E[J_alpha J_beta].

generic:  C = sigma^2 P  (P = transpose permutation (A,B)->(B,A)), hence

    E[m2]  = sigma^2 sum_m N_m 2^{-(2k-m)},
             N_m = C(Nc,m) C(Nc-m,k-m) C(Nc-k,k-m)   (m = |A n B|),
    Var m2 = 2 sigma^4 sum_{alpha beta} T_{alpha beta}^2
           = 2 sigma^4 S(Nc,k),

and the admissible-pattern generating function
f = 1 + (xy + xw + yz + zw + xyzw)/4 gives the closed form

    S(Nc,k) = sum_{j=0}^{k} sum_{m=0}^{k-j}
              Nc! / [ ((k-j-m)!)^2 (m!)^2 j! (Nc-2k+j)! ] * 4^{-(2k-j)}.

c_symmetric:  C = sigma^2 (I + P) and T = 2^{-2k} * P exactly, so m2 is a sum
of N_pairs = C(Nc,k) C(Nc-k,k)/2 iid Gaussian squares and

    E[m2]  = sigma^2 C(Nc,k) C(Nc-k,k) 4^{-k},
    Var m2 = 4 sigma^4 C(Nc,k) C(Nc-k,k) 16^{-k},
    relVar = 4 / [C(Nc,k) C(Nc-k,k)] = 2 / N_pairs        (EXACT)

-- literally the same statistical identity as the Majorana relVar = 2/C(N,p)
(sum of #couplings iid squares), with C(N,p) -> N_pairs.

Large-Nc rates (leading order; both verified against the exact forms):

    generic:      relVar(m2) -> 2 p! / Nc^p     -- the SAME leading
                  coefficient as the Majorana 2/C(N,p) ~ 2 p!/N^p;
    c_symmetric:  relVar(m2) -> 4 ((p/2)!)^2 / Nc^p
                  = [2/C(p,p/2)] * (2 p!/Nc^p).

So both U(N) ensembles self-average at the Majorana power Nc^{-p}; the
CP-sharp ensemble is asymptotically C(p,p/2)/2 times SMALLER in variance
(3x at p = 4) -- its trace is a clean sum of squares, while the generic
ensemble's admissible cross-pairings enhance Var by C(p,p/2) at leading
order.  (At small Nc the finite-size corrections invert the ranking:
exact relVar at Nc = 8, p = 4 is 7.803e-3 generic vs 9.524e-3 c_symmetric.)

Verified two independent ways: (a) closed forms vs literal enumeration --
numeric monomial traces from the sparse syk.dirac_operators route plus
full-Isserlis Wick, sharing no assembly code with the formulas; (b) vs ED
disorder ensembles at Nc = 8 and 10 (test_qed_campaign.py; Nc = 10 marked
slow).

--------------------------------------------------------------------------
2. Charge-resolved spectra and the charging curve E0(Q)
--------------------------------------------------------------------------
[H, Q] = 0, so H block-diagonalizes over popcount sectors of dim C(Nc,q).
Reassembling the sector spectra reproduces the dense spectrum machine-
exactly (tested), and off-sector blocks vanish identically.

Exact sector facts (all tested):
  * H_q = 0 identically for q < k (every term annihilates k particles), so
    E0(q) = 0 exactly there;
  * the top sector q = Nc is 1-dimensional with H = sum_A J_AA * (sign): it
    vanishes for c_symmetric (no diagonal couplings) but NOT for generic --
    the generic ensemble is NOT particle-hole symmetric even in
    distribution (the J_AA break it), while c_symmetric instances obey
    E0(q) = E0(Nc - q) exactly, sector by sector.

The charging curve is the disorder mean E0(q) +- SEM at Nc = 8, 10, p = 4,
both ensembles.  Around half filling, x = q - Nc/2, we compare on the
interior sectors k <= q <= Nc-k (where H_q is nontrivial):

    capacitive:  E0 = a + mu x + x^2/(2C)     ("Coulomb-blockade-like")
    confining :  E0 = a + mu x + s |x|        ("linear in the charge")

STATISTICS (the part that was wrong once).  The sector means E0(q) are
measured on the SAME disorder instances for every q, and each instance's
overall coupling scale moves all its sectors together, so the components
of the mean vector are strongly positively correlated: measured
interior-window correlations run from ~ +0.25 to +0.98 (generic) and up
to exactly +1 (c_symmetric mirror pairs, see below).  A chi^2 built from per-sector SEMs
as if the points were independent is therefore NOT a calibrated
likelihood statistic.  The fits here are GLS with the FULL estimated
covariance of the interior mean vector (eigenvalue-thresholded
pseudo-inverse; Hartlap factor (n - d - 2)/(n - 1) correcting the
inverse-of-estimated-covariance bias), so chi2_gls is approximately
calibrated -- up to the O(d/n_inst) noise of the estimated covariance
itself (d <= 7 fit points from n_inst >= 60 instances).  The diagonal
chi^2 is still emitted as chi2_diag, labeled a reference number only;
GLS makes the capacitive-vs-confining preference STRONGER, not weaker.

For c_symmetric the mirror E0(q) = E0(Nc - q) is instance-EXACT (tested),
so mirrored sectors carry ZERO independent information: the covariance of
the full interior window is exactly singular (mirror pairs have
correlation +1), and an independence fit would double-count every
off-center deviation.  The fit folds to q <= Nc/2 -- 3 unique interior
points at Nc = 8, not 5 -- and drops the linear term (mu = 0 exactly by
the same symmetry), leaving K = 2 parameters per model.  Generic curves
keep all interior points and K = 3 (the J_AA particle-hole breaking makes
mu physical).

Within each window both models have the same K, so AIC = chi2_gls + 2K
gives Delta-AIC = Delta chi2_gls; Delta-AIC < 2 means the window cannot
distinguish them, exactly as in self_averaging.py.  The comparison is
RELATIVE, not an absolute goodness-of-fit claim: the capacitive model
wins in every window (GLS Delta-AIC ~ 2.1-2.5e3 generic, ~ 1.5-2.4e3 in
the short folded c_symmetric windows, vs diagonal reference Delta-chi2
of only ~ 90-690 -- GLS strengthens the preference because the shared
noise mode cannot absorb a kink-vs-parabola shape difference; the OLD
independence-fit c_symmetric Delta-AIC values 530/1027 double-counted
mirrored sectors and are retired), but its own chi2_gls/dof is >> 1 for
the generic curves (the data resolve corrections beyond a pure parabola
-- quartic terms and edge effects); chi2_gls/dof is reported next to
every fit so nobody mistakes "much better than linear" for "exactly
quadratic".

INTERPRETATION (what a huge Delta-AIC buys and what it does not).
Decisively rejecting the kinked a + mu x + s|x| alternative is mostly a
test of the SMOOTHNESS/analyticity of E0(Q) at half filling: any chaotic
charge-conserving Hamiltonian whose ground-energy density is a smooth
function of the filling would reject the |x| kink just as hard, so the
model comparison by itself is WEAK evidence FOR a specific
Coulomb/capacitive bulk reading.  The positive physical content is the
measured curvature (the capacitance C and its Nc trend), not the win
over the strawman.  This is the first U(1)-sector energetics data for
the duality.

--------------------------------------------------------------------------
3. Charged vs singlet spectral weight
--------------------------------------------------------------------------
T = infinity SYMMETRIZED spectral functions (moments_pipeline conventions:
split O into Hermitian + anti-Hermitian parts, add the two omega-symmetric
weights; odd moments vanish identically):

  charged:  O = c_i, averaged over i and disorder.  Exact anchors:
      mu0 = (Tr c c^dag + Tr c^dag c)/(2 dim) = 1/2   per mode, exactly;
      mu2 = ||[H, c_i]||_F^2 / dim                    per instance, exactly
            (the Hermitian/anti-Hermitian cross term cancels; tested), and
      E[mu2] = (sigma^2/Nc) sum_m N_m 2^{-(2k-m)} (2k - m)     (generic)
      E[mu2] = (p/Nc) E[m2]                                    (c_symmetric)
    from ||[O_alpha, c_i]||_F^2/dim = 2^{-|AuB|} * (2 if i in A\B; 1 if
    i in A n B; 0 else) -- occupation counting again -- summed over keys.
  singlet:  O = M2bar = traceless M_2 = -sum_i cdot^dag_i cdot_i (the first
    dynamical singlet: M_0 = Q and M_1 = -(p/2)iH are conserved), assembled
    in the H eigenbasis as cdot~ = i omega o c~ and checked against
    dirac.mn_tower in the tests.  M_2 is Hermitian so symmetrization is
    trivial; we use the raw (not CP-resolved) M2bar so both ensembles are
    treated identically -- on c_symmetric instances its weight is the sum of
    the two CP channels of moments_pipeline.py.

Reported scale: w2 = mu2/(mu0 * Tr(H^2)/dim), the mean-square frequency in
units of the instance's own energy scale (the moments_pipeline invariant),
plus the charged/singlet ratio.  The exact anchors give the ensemble-mean
prediction w2_charged ~= 2p/Nc (1 + O(1/Nc)) for both ensembles.  For
c_symmetric the statement is far stronger and INSTANCE-EXACT:

    mu2_charged = (p/Nc) m2   per instance  =>  w2_charged = 2p/Nc exactly

-- the same disjoint-key pattern argument that makes T a pairing delta also
kills the cross terms of the commutator quadratic form
sum_i ||[H, c_i]||^2, locking the charged frequency scale rigidly to the
instance's total energy scale (verified to machine precision in the tests;
for generic instances the cross terms push E[mu2_charged] BELOW
(p/Nc) E[m2] by exactly 23.1% (Nc=6), 20.6% (Nc=7), 18.6% (Nc=8) at
p = 4 -- a deficit that shrinks with Nc, from the exact formulas above).

HONESTY BOX -- what this does and does not test.  In 0+1d there is no
spatial confinement: no Wilson loop, no linear potential in a distance, no
"single quark removed to infinity".  The 0+1d shadow of the paper's
"confinement of U(1) charge" consists of exactly two sharp statements:
(i) the sector energetics E0(Q) of deliverable 2 -- how the ground energy
rises as the hologram is charged (the paper's Coulomb-force claim via
arXiv:2012.12326), and (ii) the frequency-scale separation between charged
(quark-like c_i) and singlet (meson-like M_n) channels measured here.
Whether charged weight sits at higher frequency than singlet weight is a
measured statement about those scales, NOT a demonstration of confinement;
the bulk-side meaning of the comparison further depends on the unresolved
hologram/bulk frequency map (duality/BOOST.md).  We say so wherever the
numbers are reported.

Statistics: every ensemble number carries an SEM with ddof=1; the charging
curves additionally carry the full covariance of the mean vector (the
sectors are measured on the same instances and are strongly correlated --
see part 2), and their model fits are GLS against that covariance;
disorder instances are seeded reproducibly.  The default __main__ takes
~2-3 min and writes qed_campaign.json + qed_campaign.png; --big adds the
Nc = 10 spectral run and larger ensembles.
"""

from __future__ import annotations

import argparse
import json
from functools import lru_cache
from math import comb, factorial
from pathlib import Path

import numpy as np

from dirac import (annihilators, coupling_variance_dirac, dirac_hamiltonian,
                   draw_couplings, symmetrized_weight, time_derivative)


# --------------------------------------------------------------------------
# 1a. closed forms for the m2 disorder statistics
# --------------------------------------------------------------------------
def t_ab_exact(A, B) -> float:
    """T_AB = Tr(c^dag_A c_B c^dag_B c_A)/dim = 2^{-|A u B|} (occupation
    counting: the |A u B| modes are forced, the rest free)."""
    return 2.0 ** -len(set(A) | set(B))


def _n_pairs_overlap(Nc: int, k: int, m: int) -> int:
    """Number of ordered pairs (A, B) of k-subsets of [Nc] with |A n B| = m."""
    return comb(Nc, m) * comb(Nc - m, k - m) * comb(Nc - k, k - m)


def s_gen_exact(Nc: int, p: int) -> float:
    """S(Nc,k) = sum_{alpha beta} T_{alpha beta}^2 over ALL ordered key pairs
    of the generic ensemble: the admissible-pattern double sum (see module
    docstring).  j = #modes with pattern (1111), m = #modes with (1001)=(0110)
    partner count; u = |A u B u C u D| = 2k - j."""
    k = p // 2
    total = 0.0
    for j in range(k + 1):
        for m in range(k - j + 1):
            n1 = k - j - m
            if Nc - 2 * k + j < 0:
                continue
            count = (factorial(Nc)
                     // (factorial(n1) ** 2 * factorial(m) ** 2
                         * factorial(j) * factorial(Nc - 2 * k + j)))
            total += count * 4.0 ** -(2 * k - j)
    return total


def mean_m2_exact(Nc: int, p: int, ensemble: str = "generic") -> float:
    """E[m2] in closed form (module docstring, part 1)."""
    k = p // 2
    sigma2 = coupling_variance_dirac(Nc, p)
    if ensemble == "c_symmetric":
        return sigma2 * comb(Nc, k) * comb(Nc - k, k) * 4.0 ** -k
    if ensemble != "generic":
        raise ValueError(f"unknown ensemble {ensemble!r}")
    return sigma2 * sum(_n_pairs_overlap(Nc, k, m) * 2.0 ** -(2 * k - m)
                        for m in range(k + 1))


def var_m2_exact(Nc: int, p: int, ensemble: str = "generic") -> float:
    """Var(m2) in closed form: 2 sigma^4 S(Nc,k) (generic);
    4 sigma^4 C(Nc,k) C(Nc-k,k) 16^{-k} (c_symmetric)."""
    k = p // 2
    sigma2 = coupling_variance_dirac(Nc, p)
    if ensemble == "c_symmetric":
        return 4.0 * sigma2 ** 2 * comb(Nc, k) * comb(Nc - k, k) * 16.0 ** -k
    if ensemble != "generic":
        raise ValueError(f"unknown ensemble {ensemble!r}")
    return 2.0 * sigma2 ** 2 * s_gen_exact(Nc, p)


def relvar_m2_exact_dirac(Nc: int, p: int, ensemble: str = "generic") -> float:
    """Exact relative variance of m2 = Tr(H^2)/dim over disorder."""
    return var_m2_exact(Nc, p, ensemble) / mean_m2_exact(Nc, p, ensemble) ** 2


def relvar_m2_leading(Nc: int, p: int, ensemble: str = "generic") -> float:
    """Leading large-Nc rate: 2 p!/Nc^p (generic; the Majorana coefficient),
    4 ((p/2)!)^2/Nc^p (c_symmetric)."""
    if ensemble == "c_symmetric":
        return 4.0 * factorial(p // 2) ** 2 / Nc ** p
    return 2.0 * factorial(p) / Nc ** p


def mean_mu2_charged_exact(Nc: int, p: int, ensemble: str = "generic") -> float:
    r"""E[mu2] of the mode-averaged symmetrized c_i spectral function:
    E[ ||[H, c_i]||_F^2/dim ] = (sigma^2/Nc) sum_alpha 2^{-|AuB|} (2k - |AnB|)
    (each stored key contributes 2 per i in A\B and 1 per i in A n B)."""
    k = p // 2
    sigma2 = coupling_variance_dirac(Nc, p)
    if ensemble == "c_symmetric":
        return (2.0 * k / Nc) * mean_m2_exact(Nc, p, "c_symmetric")
    if ensemble != "generic":
        raise ValueError(f"unknown ensemble {ensemble!r}")
    return (sigma2 / Nc) * sum(
        _n_pairs_overlap(Nc, k, m) * 2.0 ** -(2 * k - m) * (2 * k - m)
        for m in range(k + 1))


# --------------------------------------------------------------------------
# 1b. literal enumeration (verification route -- shares no assembly code
#     with the closed forms: sparse-operator traces + full Isserlis Wick)
# --------------------------------------------------------------------------
def monomial_matrix(Nc: int, A, B) -> np.ndarray:
    """Dense O_(A,B) = c^dag_A c_B from the sparse Jordan-Wigner operators of
    syk.dirac_operators (independent of dirac.apply_monomial)."""
    from syk import dirac_operators
    c, cdag = dirac_operators(Nc)
    dim = 1 << Nc
    O = np.eye(dim, dtype=complex)
    for a in A:
        O = O @ cdag[a].toarray()
    for b in B:
        O = O @ c[b].toarray()
    return O


def stored_keys(Nc: int, p: int, ensemble: str) -> list:
    """The ordered key set of dirac.draw_couplings for this ensemble."""
    return list(draw_couplings(Nc, p, np.random.default_rng(0), ensemble))


def trace_pair_matrix(Nc: int, keys: list) -> np.ndarray:
    """T_{alpha beta} = Tr(O_alpha O_beta)/dim for every ordered key pair,
    via one flattened gemm over the dense monomial matrices."""
    dim = 1 << Nc
    Os = np.stack([monomial_matrix(Nc, A, B) for (A, B) in keys])
    F = Os.reshape(len(keys), -1)
    G = Os.transpose(0, 2, 1).reshape(len(keys), -1)
    T = (F @ G.T) / dim
    if np.abs(T.imag).max() > 1e-12:
        raise AssertionError("trace pair matrix must be real")
    return np.ascontiguousarray(T.real)


def coupling_covariance_matrix(keys: list, sigma2: float,
                               ensemble: str) -> np.ndarray:
    """C_{alpha beta} = E[J_alpha J_beta] over the stored keys.
    generic: sigma^2 * P (transpose permutation; E z^2 = 0 off-diagonal,
    real diagonal).  c_symmetric: sigma^2 (I + P) (J real, J_AB = J_BA)."""
    idx = {key: i for i, key in enumerate(keys)}
    n = len(keys)
    C = np.zeros((n, n))
    for i, (A, B) in enumerate(keys):
        C[i, idx[(B, A)]] += sigma2
        if ensemble == "c_symmetric":
            C[i, i] += sigma2
    return C


def enumerated_m2_stats(Nc: int, p: int, ensemble: str):
    """(E[m2], Var m2) by literal Isserlis over numeric traces:
    E = sum C o T,  E[m2^2] = E^2 + 2 tr[(C T)^2]  (the two nontrivial
    pairings are equal by symmetry of T and C).  No closed-form input."""
    keys = stored_keys(Nc, p, ensemble)
    T = trace_pair_matrix(Nc, keys)
    C = coupling_covariance_matrix(keys, coupling_variance_dirac(Nc, p),
                                   ensemble)
    mean = float(np.sum(C * T))
    CT = C @ T
    var = 2.0 * float(np.trace(CT @ CT))
    return mean, var


def enumerated_mu2_charged(Nc: int, p: int, ensemble: str) -> float:
    """E[mu2_charged] = (sigma^2/Nc) sum_i sum_alpha ||[O_alpha, c_i]||_F^2/dim
    with literal commutators of numeric matrices (verification route)."""
    from syk import dirac_operators
    c, _ = dirac_operators(Nc)
    keys = stored_keys(Nc, p, ensemble)
    dim = 1 << Nc
    sigma2 = coupling_variance_dirac(Nc, p)
    total = 0.0
    for (A, B) in keys:
        O = monomial_matrix(Nc, A, B)
        for i in range(Nc):
            X = O @ c[i].toarray() - c[i].toarray() @ O
            total += float(np.sum(np.abs(X) ** 2)) / dim
    return sigma2 * total / Nc


# --------------------------------------------------------------------------
# 1c. ED disorder ensemble for m2
# --------------------------------------------------------------------------
def m2_of_hamiltonian(H: np.ndarray) -> float:
    """Tr(H^2)/dim = ||H||_F^2/dim (H Hermitian)."""
    return float(np.sum(np.abs(H) ** 2)) / H.shape[0]


def measure_m2_ensemble(Nc: int, p: int, ensemble: str, n_inst: int,
                        seed: int = 0) -> dict:
    """Sample mean and relative variance of m2 over n_inst draws, with the
    repo's error conventions (SEM ddof=1; relVar error ~ relVar sqrt(2/(n-1)))."""
    rng = np.random.default_rng(seed)
    vals = np.empty(n_inst)
    for j in range(n_inst):
        coup = draw_couplings(Nc, p, rng, ensemble)
        vals[j] = m2_of_hamiltonian(dirac_hamiltonian(Nc, p, couplings=coup))
    mean = float(vals.mean())
    relvar = float(vals.var(ddof=1) / mean ** 2)
    return dict(Nc=Nc, p=p, ensemble=ensemble, n_inst=n_inst,
                m2_mean=mean, m2_sem=float(vals.std(ddof=1) / np.sqrt(n_inst)),
                relvar=relvar,
                relvar_err=relvar * float(np.sqrt(2.0 / (n_inst - 1))),
                relvar_exact=relvar_m2_exact_dirac(Nc, p, ensemble),
                mean_exact=mean_m2_exact(Nc, p, ensemble))


# --------------------------------------------------------------------------
# 2. charge sectors and the charging curve
# --------------------------------------------------------------------------
def sector_indices(Nc: int) -> list:
    """Basis-state indices of each charge sector q = 0..Nc (popcount)."""
    states = np.arange(1 << Nc, dtype=np.int64)
    pop = np.bitwise_count(states)
    return [np.flatnonzero(pop == q) for q in range(Nc + 1)]


def sector_spectra(H: np.ndarray, Nc: int, sectors: list = None) -> list:
    """Eigenvalues of every charge block H_q (ascending within each sector)."""
    if sectors is None:
        sectors = sector_indices(Nc)
    return [np.linalg.eigvalsh(H[np.ix_(ix, ix)]) for ix in sectors]


def charging_curve(Nc: int, p: int, ensemble: str, n_inst: int,
                   seed: int = 0) -> dict:
    """Disorder-averaged ground energy per charge sector, E0(q) +- SEM,
    PLUS the full covariance of the mean vector (Cov(E0)/n_inst, ddof=1):
    the sectors of one instance are strongly correlated across q, so the
    SEMs alone under-specify the error model of the mean curve."""
    rng = np.random.default_rng(seed)
    sectors = sector_indices(Nc)
    e0 = np.empty((n_inst, Nc + 1))
    for j in range(n_inst):
        coup = draw_couplings(Nc, p, rng, ensemble)
        H = dirac_hamiltonian(Nc, p, couplings=coup)
        for q, ix in enumerate(sectors):
            Hq = H[np.ix_(ix, ix)]
            e0[j, q] = np.linalg.eigvalsh(Hq)[0].real if len(ix) > 1 \
                else Hq[0, 0].real
    cov_mean = np.cov(e0, rowvar=False) / n_inst      # ddof=1, matches SEM^2
    return dict(Nc=Nc, p=p, ensemble=ensemble, n_inst=n_inst,
                q=list(range(Nc + 1)),
                e0_mean=[float(x) for x in e0.mean(axis=0)],
                e0_sem=[float(x) for x in e0.std(axis=0, ddof=1)
                        / np.sqrt(n_inst)],
                e0_cov_mean=[[float(v) for v in row] for row in cov_mean])


def _gls_fit(A: np.ndarray, y: np.ndarray, cov: np.ndarray,
             n_cov: int = None):
    """Generalized least squares: coef minimizing r^T Cinv r, with Cinv the
    eigenvalue-thresholded pseudo-inverse of the (estimated) covariance of
    y.  When cov was estimated from n_cov samples, the Hartlap factor
    (n_cov - rank - 2)/(n_cov - 1) corrects the mean bias of the inverse
    (it rescales chi2 only; the coefficients are scale-invariant).
    Returns (coef, chi2, rank, cond, hartlap)."""
    w, U = np.linalg.eigh(0.5 * (cov + cov.T))
    tol = w.max() * 1e-10
    good = w > tol
    rank = int(good.sum())
    winv = np.where(good, 1.0 / np.where(good, w, 1.0), 0.0)
    Cinv = (U * winv) @ U.T
    hartlap = 1.0
    if n_cov is not None and n_cov > rank + 2:
        hartlap = (n_cov - rank - 2) / (n_cov - 1.0)
    AtC = A.T @ Cinv
    coef = np.linalg.solve(AtC @ A, AtC @ y)
    r = y - A @ coef
    chi2 = hartlap * float(r @ Cinv @ r)
    cond = float(w.max() / w[good].min())
    return coef, chi2, rank, cond, hartlap


def fit_charging(curve: dict) -> dict:
    """GLS fit of the interior charging curve (k <= q <= Nc - k) against the
    capacitive and confining models (module docstring, part 2).

    The E0(q) means are measured on the same instances and are strongly
    correlated across q, so the fit uses the FULL covariance of the mean
    vector (curve["e0_cov_mean"]; falls back to diagonal SEM^2 if absent,
    e.g. for synthetic curves).  For c_symmetric the mirror
    E0(q) = E0(Nc-q) is instance-exact, making mirrored sectors literal
    duplicates (covariance exactly singular): the fit folds to q <= Nc/2
    and drops the linear term (mu = 0 exactly), K = 2; generic keeps the
    full window and K = 3.  Both models share K, so AIC = chi2_gls + 2K
    gives Delta-AIC = Delta chi2_gls; |Delta-AIC| < 2 = indistinguishable.
    chi2_diag (per-sector SEMs, independence assumed) is emitted as a
    reference number only -- it is NOT a calibrated statistic."""
    Nc, k = curve["Nc"], curve["p"] // 2
    q = np.array(curve["q"], dtype=float)
    y_all = np.array(curve["e0_mean"])
    s_all = np.array(curve["e0_sem"])
    cov_all = (np.array(curve["e0_cov_mean"]) if "e0_cov_mean" in curve
               else np.diag(s_all ** 2))
    keep = (q >= k) & (q <= Nc - k)
    folded = curve.get("ensemble") == "c_symmetric"
    mirror_max_dev = None
    if folded:
        interior = np.flatnonzero(keep)
        mirror_max_dev = float(max(abs(y_all[i] - y_all[Nc - i])
                                   for i in interior))
        keep &= (q <= Nc / 2.0)
    idx = np.flatnonzero(keep)
    x = q[idx] - Nc / 2.0
    y = y_all[idx]
    s = s_all[idx]
    C = cov_all[np.ix_(idx, idx)]
    n_pts = len(idx)
    K = 2 if folded else 3
    dof = n_pts - K
    n_cov = curve.get("n_inst")
    corr = C / np.sqrt(np.outer(np.diag(C), np.diag(C)))
    off = corr[~np.eye(n_pts, dtype=bool)]
    out = dict(folded=bool(folded), n_points=n_pts, K=K, dof=dof,
               corr_offdiag_range=([float(off.min()), float(off.max())]
                                   if off.size else None))
    if mirror_max_dev is not None:
        out["mirror_max_dev"] = mirror_max_dev
    for name, col in [("capacitive_x2", x ** 2), ("confining_absx", np.abs(x))]:
        cols = ([np.ones_like(x), col] if folded
                else [np.ones_like(x), x, col])
        A = np.vstack(cols).T
        coef, chi2, rank, cond, hartlap = _gls_fit(A, y, C, n_cov)
        coefd, *_ = np.linalg.lstsq(A / s[:, None], y / s, rcond=None)
        chi2_diag = float(np.sum(((y - A @ coefd) / s) ** 2))
        out[name] = dict(coef=[float(c) for c in coef],
                         chi2_gls=chi2, aic=chi2 + 2 * K,
                         chi2_gls_per_dof=chi2 / dof if dof > 0 else np.nan,
                         chi2_diag=chi2_diag,
                         chi2_diag_per_dof=(chi2_diag / dof if dof > 0
                                            else np.nan),
                         cov_rank=rank, cov_cond=cond, hartlap=hartlap)
    quad_c = out["capacitive_x2"]["coef"][-1]
    out["capacitance"] = float(1.0 / (2.0 * quad_c)) if quad_c != 0 else np.inf
    daic = out["confining_absx"]["aic"] - out["capacitive_x2"]["aic"]
    out["delta_aic_conf_minus_cap"] = daic
    out["delta_chi2_diag_conf_minus_cap"] = (
        out["confining_absx"]["chi2_diag"]
        - out["capacitive_x2"]["chi2_diag"])
    if abs(daic) < 2.0:
        out["verdict"] = ("indistinguishable in this window "
                          f"(|Delta-AIC| = {abs(daic):.1f} < 2)")
    else:
        best = "capacitive_x2" if daic > 0 else "confining_absx"
        out["verdict"] = f"best: {best} (Delta-AIC = {abs(daic):.1f})"
    return out


# --------------------------------------------------------------------------
# 3. charged vs singlet spectral weight
# --------------------------------------------------------------------------
def _sym_moments(W: np.ndarray, omega: np.ndarray) -> dict:
    return {f"mu{n}": float(np.sum(W * omega ** n)) for n in (0, 1, 2, 4)}


_annihilators_cached = lru_cache(maxsize=8)(annihilators)   # per-Nc, not per-instance


def instance_spectral(Nc: int, p: int, ensemble: str, seed: int) -> dict:
    """Symmetrized T = infinity spectral moments of the charged operator c_i
    (mode-averaged) and of the singlet M2bar, for one disorder instance.

    Everything is assembled in the H eigenbasis: c~ = V^dag c V, the
    commutator is the elementwise product [H, c]~ = omega o c~, and
    M2 = -sum_i cdot^dag_i cdot_i.  (The eigenbasis assembly is asserted
    against dirac.mn_tower in the tests.)"""
    rng = np.random.default_rng(seed)
    coup = draw_couplings(Nc, p, rng, ensemble)
    H = dirac_hamiltonian(Nc, p, couplings=coup)
    dim = 1 << Nc
    E, V = np.linalg.eigh(H)
    omega = E[:, None] - E[None, :]
    escale2 = float(np.mean(E ** 2))                 # Tr(H^2)/dim
    cs = _annihilators_cached(Nc)

    Wc = np.zeros((dim, dim))
    M2t = np.zeros((dim, dim), dtype=complex)
    for c in cs:
        ct = V.conj().T @ c @ V
        Wc += symmetrized_weight(ct)
        cdot = 1j * omega * ct
        M2t -= cdot.conj().T @ cdot
    Wc /= Nc                                         # mode-averaged, mu0 = 1/2
    charged = _sym_moments(Wc, omega)

    M2t -= (np.trace(M2t) / dim) * np.eye(dim)       # traceless M2bar
    herm_dev = float(np.abs(M2t - M2t.conj().T).max())
    Ws = np.abs(M2t) ** 2 / dim                      # Hermitian: symmetric
    singlet = _sym_moments(Ws, omega)

    return dict(charged=charged, singlet=singlet, escale2=escale2,
                m2_herm_dev=herm_dev,
                w2_charged=charged["mu2"] / (charged["mu0"] * escale2),
                w2_singlet=singlet["mu2"] / (singlet["mu0"] * escale2),
                kurt_charged=charged["mu4"] * charged["mu0"] / charged["mu2"] ** 2,
                kurt_singlet=singlet["mu4"] * singlet["mu0"] / singlet["mu2"] ** 2)


def mu2_charged_by_commutator(H: np.ndarray, cs: list) -> float:
    """Mode-averaged mu2 of the symmetrized c_i spectral function via the
    exact identity mu2 = ||[H, c_i]||_F^2/dim -- an eigenbasis-free route
    used to cross-validate instance_spectral."""
    dim = H.shape[0]
    total = 0.0
    for c in cs:
        X = time_derivative(H, c)
        total += float(np.sum(np.abs(X) ** 2)) / dim
    return total / len(cs)


def measure_spectral(Nc: int, p: int, ensemble: str, n_inst: int,
                     seed0: int = 500) -> dict:
    """Disorder mean +- SEM of the charged/singlet frequency-scale
    observables."""
    fields = ("w2_charged", "w2_singlet", "kurt_charged", "kurt_singlet")
    rows = {f: [] for f in fields}
    rows["ratio"] = []
    mu2c, mu0c, odd = [], [], []
    for j in range(n_inst):
        r = instance_spectral(Nc, p, ensemble, seed0 + j)
        for f in fields:
            rows[f].append(r[f])
        rows["ratio"].append(r["w2_charged"] / r["w2_singlet"])
        mu2c.append(r["charged"]["mu2"])
        mu0c.append(r["charged"]["mu0"])
        odd.append(abs(r["charged"]["mu1"]) + abs(r["singlet"]["mu1"]))
    out = dict(Nc=Nc, p=p, ensemble=ensemble, n_inst=n_inst,
               mu0_charged=float(np.mean(mu0c)),      # exact 1/2
               mu2_charged_mean=float(np.mean(mu2c)),
               mu2_charged_sem=float(np.std(mu2c, ddof=1) / np.sqrt(n_inst)),
               mu2_charged_exact=mean_mu2_charged_exact(Nc, p, ensemble),
               odd_moment_check=float(np.max(odd)))
    for f in list(fields) + ["ratio"]:
        arr = np.asarray(rows[f])
        out[f] = float(arr.mean())
        out[f + "_sem"] = float(arr.std(ddof=1) / np.sqrt(n_inst))
    return out


# --------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------
_COLORS = {"generic": "#2166ac", "c_symmetric": "#b2182b"}


def make_figure(results: dict, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axL, axM, axR) = plt.subplots(1, 3, figsize=(15, 4.8))

    # ---- left: relVar(m2) exact vs Nc, both ensembles, with ED points ----
    ns = np.arange(6, 21, 2)
    for ens, color in _COLORS.items():
        axL.loglog(ns, [relvar_m2_exact_dirac(int(N), 4, ens) for N in ns],
                   "-", color=color, label=f"{ens}: exact")
        axL.loglog(ns, [relvar_m2_leading(int(N), 4, ens) for N in ns],
                   "--", color=color, lw=1.0, alpha=0.6,
                   label=(r"$2\,p!/N_c^p$" if ens == "generic"
                          else r"$4\,((p/2)!)^2/N_c^p$"))
    for row in results["m2_ensembles"]:
        color = _COLORS[row["ensemble"]]
        axL.errorbar([row["Nc"]], [row["relvar"]], yerr=[row["relvar_err"]],
                     fmt="x", color="#333333", capsize=3, zorder=5)
    axL.errorbar([], [], fmt="x", color="#333333", capsize=3,
                 label="ED disorder ensembles")
    axL.set_xlabel("$N_c$")
    axL.set_ylabel(r"relVar$(m_2)$, $m_2 = \mathrm{Tr}(H^2)/\mathrm{dim}$")
    axL.set_title("U(N) self-averaging ($p=4$): exact vs ED")
    axL.set_xlim(5.5, 21.5)
    axL.set_xticks([6, 8, 10, 14, 20])
    axL.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axL.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axL.legend(frameon=False, fontsize=8, loc="lower left")
    axL.grid(True, which="both", alpha=0.2)

    # ---- middle: charging curves E0(q) with fits -------------------------
    for curve in results["charging"]:
        Nc, ens = curve["Nc"], curve["ensemble"]
        color = _COLORS[ens]
        x = np.array(curve["q"], dtype=float) - Nc / 2.0
        fmt = "o" if Nc == 10 else "s"
        alpha = 1.0 if Nc == 10 else 0.45
        axM.errorbar(x, curve["e0_mean"], yerr=curve["e0_sem"], fmt=fmt,
                     ms=4.5, color=color, alpha=alpha, capsize=2,
                     label=f"{ens}, $N_c={Nc}$")
        if Nc == 10:
            fit = curve["fits"]["capacitive_x2"]
            k = curve["p"] // 2
            xf = np.linspace(k - Nc / 2.0, Nc / 2.0 - k, 120)
            if len(fit["coef"]) == 2:      # folded c_symmetric: mu = 0 exact
                a, cq = fit["coef"]
                mu = 0.0
            else:
                a, mu, cq = fit["coef"]
            axM.plot(xf, a + mu * xf + cq * xf ** 2, "-", lw=1.0,
                     color=color, alpha=0.7)
    axM.set_xlabel(r"$x = Q - N_c/2$")
    axM.set_ylabel(r"$\overline{E_0(Q)}$   (ground energy per sector)")
    axM.set_title("Charging curves (curves: capacitive fit, interior $Q$)")
    axM.legend(frameon=False, fontsize=8)
    axM.grid(True, alpha=0.2)
    axM.text(0.5, 0.9,
             "generic top sector is 1-dim ($J_{AA}$ only):\n"
             "particle-hole asymmetry is real, not noise",
             transform=axM.transAxes, fontsize=7.5, color="0.4",
             ha="center", va="top")

    # ---- right: charged vs singlet frequency scales ----------------------
    xpos, labels = [], []
    for i, row in enumerate(results["spectral"]):
        color = _COLORS[row["ensemble"]]
        base = i * 1.0
        axR.errorbar([base - 0.12], [row["w2_charged"]],
                     yerr=[row["w2_charged_sem"]], fmt="o", color=color,
                     capsize=3)
        axR.errorbar([base + 0.12], [row["w2_singlet"]],
                     yerr=[row["w2_singlet_sem"]], fmt="s", color=color,
                     capsize=3, mfc="none")
        axR.annotate(f"ratio {row['ratio']:.2f}",
                     (base, max(row["w2_charged"], row["w2_singlet"]) * 1.06),
                     ha="center", fontsize=8, color="0.3")
        xpos.append(base)
        labels.append(f"{row['ensemble']}\n$N_c={row['Nc']}$")
    axR.errorbar([], [], fmt="o", color="0.3", label=r"charged $c_i$")
    axR.errorbar([], [], fmt="s", mfc="none", color="0.3",
                 label=r"singlet $\bar M_2$")
    axR.set_xticks(xpos)
    axR.set_xticklabels(labels, fontsize=8)
    axR.set_xlim(min(xpos) - 0.6, max(xpos) + 0.6)
    axR.set_ylim(0.0, 1.18)
    axR.set_ylabel(r"$w_2 = \mu_2/(\mu_0\,\mathrm{Tr}(H^2)/\mathrm{dim})$")
    axR.set_title("Charged vs singlet frequency scale")
    axR.legend(frameon=False, fontsize=8, loc="lower right")
    axR.grid(True, alpha=0.2)
    axR.text(0.03, 0.03,
             "frequency-scale separation, not spatial\n"
             "confinement (0+1d) -- see module docstring",
             transform=axR.transAxes, fontsize=7.5, color="0.4")

    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"wrote {out_png}")


# --------------------------------------------------------------------------
# campaign driver
# --------------------------------------------------------------------------
def main(big: bool = False) -> dict:
    here = Path(__file__).resolve().parent
    p = 4
    results = {"p": p, "citation_context": (
        "arXiv:2012.12326 (Susskind, 'Electric Forces in the Charged SYK "
        "Model') is cited once in arXiv:2607.05678, at the M_0 discussion: "
        "'photon exchange gives rise to long range confining Coulomb "
        "forces' -- the claim behind the charging-curve measurement.")}

    print("=" * 72)
    print("1. U(N) self-averaging of m2 = Tr(H^2)/dim: exact vs ED")
    print("=" * 72)
    print(f"\n  exact relVar(m2), p = {p} "
          "(generic leading 2p!/Nc^p, c_symmetric leading 4((p/2)!)^2/Nc^p):")
    results["m2_exact_table"] = []
    print(f"  {'Nc':>4} {'generic':>12} {'leading':>12} "
          f"{'c_symmetric':>12} {'leading':>12}")
    for Nc in (6, 8, 10, 12, 14, 16):
        row = dict(Nc=Nc,
                   generic=relvar_m2_exact_dirac(Nc, p, "generic"),
                   generic_leading=relvar_m2_leading(Nc, p, "generic"),
                   c_symmetric=relvar_m2_exact_dirac(Nc, p, "c_symmetric"),
                   c_symmetric_leading=relvar_m2_leading(Nc, p, "c_symmetric"))
        results["m2_exact_table"].append(row)
        print(f"  {Nc:>4} {row['generic']:>12.4e} {row['generic_leading']:>12.4e}"
              f" {row['c_symmetric']:>12.4e} {row['c_symmetric_leading']:>12.4e}")

    results["m2_ensembles"] = []
    plan = [(8, 300), (10, 100 if not big else 300)]
    print("\n  ED disorder ensembles (sampling error on relVar "
          "~ relVar*sqrt(2/(n-1))):")
    for Nc, n_inst in plan:
        for ens in ("generic", "c_symmetric"):
            r = measure_m2_ensemble(Nc, p, ens, n_inst, seed=40 + Nc)
            results["m2_ensembles"].append(r)
            pull = (r["relvar"] - r["relvar_exact"]) / r["relvar_err"]
            print(f"    Nc={Nc:>2} {ens:>12} n={n_inst:>3}: "
                  f"relVar = {r['relvar']:.3e} +- {r['relvar_err']:.1e} "
                  f"(exact {r['relvar_exact']:.3e}, pull {pull:+.1f}s); "
                  f"<m2> = {r['m2_mean']:.4f} +- {r['m2_sem']:.4f} "
                  f"(exact {r['mean_exact']:.4f})")

    print("\n" + "=" * 72)
    print("2. Charging curves E0(Q): the U(1)-sector energetics")
    print("=" * 72)
    results["charging"] = []
    for Nc, n_inst in [(8, 200 if not big else 400),
                       (10, 60 if not big else 150)]:
        for ens in ("generic", "c_symmetric"):
            curve = charging_curve(Nc, p, ens, n_inst, seed=60 + Nc)
            curve["fits"] = fit_charging(curve)
            results["charging"].append(curve)
            f = curve["fits"]
            print(f"\n  Nc={Nc}, {ens}, {n_inst} draws "
                  f"(interior fit window {p//2} <= Q <= {Nc - p//2}):")
            e0 = curve["e0_mean"]
            sem = curve["e0_sem"]
            print("    Q     : " + " ".join(f"{q:>8d}" for q in curve["q"]))
            print("    E0    : " + " ".join(f"{v:>8.3f}" for v in e0))
            print("    +-SEM : " + " ".join(f"{v:>8.3f}" for v in sem))
            lo, hi = f["corr_offdiag_range"]
            print(f"    [mean-vector cross-Q correlations in [{lo:+.2f}, "
                  f"{hi:+.2f}] -- fits are GLS with the full covariance]")
            if f["folded"]:
                print("    [c_symmetric mirror E0(q) = E0(Nc-q) is "
                      f"instance-exact (max dev {f['mirror_max_dev']:.1e}): "
                      f"folded to Q <= Nc/2, {f['n_points']} unique pts, "
                      "mu = 0 exact, K = 2]")
            cap = f["capacitive_x2"]
            con = f["confining_absx"]
            if f["folded"]:
                cap_str = (f"E0 = {cap['coef'][0]:.3f} "
                           f"{cap['coef'][1]:+.3f} x^2")
                con_str = (f"E0 = {con['coef'][0]:.3f} "
                           f"{con['coef'][1]:+.3f} |x|")
            else:
                cap_str = (f"E0 = {cap['coef'][0]:.3f} "
                           f"{cap['coef'][1]:+.3f} x {cap['coef'][2]:+.3f} x^2")
                con_str = (f"E0 = {con['coef'][0]:.3f} "
                           f"{con['coef'][1]:+.3f} x {con['coef'][2]:+.3f} |x|")
            print(f"    capacitive:  {cap_str}   "
                  f"chi2_gls/dof = {cap['chi2_gls_per_dof']:.1f} "
                  f"({f['n_points']} pts, dof {f['dof']}); "
                  f"C = {f['capacitance']:.3f}")
            print(f"    confining :  {con_str}   "
                  f"chi2_gls/dof = {con['chi2_gls_per_dof']:.1f}")
            print(f"    -> {f['verdict']}  [GLS Delta-chi2; diagonal-"
                  f"independence reference Delta-chi2 = "
                  f"{f['delta_chi2_diag_conf_minus_cap']:.1f} is NOT "
                  "calibrated (correlated sectors)]")
            print("       [relative comparison only: chi2_gls/dof > 1 means "
                  "resolvable non-parabolic corrections]")

    print("\n  (Rejecting the |x| kink mostly tests smoothness/analyticity "
          "of E0(Q) at")
    print("   half filling -- any chaotic charge-conserving H with a smooth "
          "ground-")
    print("   energy density does the same; weak evidence FOR a specific "
          "Coulomb/")
    print("   capacitive bulk reading.  The physics is C and its Nc trend.)")

    print("\n" + "=" * 72)
    print("3. Charged (c_i) vs singlet (M2bar) spectral weight, T = infinity")
    print("=" * 72)
    results["spectral"] = []
    spec_plan = [(8, 32 if not big else 80)]
    if big:
        spec_plan.append((10, 6))
    for Nc, n_inst in spec_plan:
        for ens in ("generic", "c_symmetric"):
            r = measure_spectral(Nc, p, ens, n_inst, seed0=500 + Nc)
            results["spectral"].append(r)
            print(f"\n  Nc={Nc}, {ens}, {n_inst} draws:")
            print(f"    mu0_charged = {r['mu0_charged']:.6f} (exact 0.5); "
                  f"mu2_charged = {r['mu2_charged_mean']:.4f} "
                  f"+- {r['mu2_charged_sem']:.4f} "
                  f"(exact E[mu2] = {r['mu2_charged_exact']:.4f})")
            if ens == "c_symmetric":
                anchor_label = "instance-EXACT identity 2p/Nc"
                anchor_val = 2 * p / Nc
            else:
                # generic: 2p/Nc is the c_symmetric value; the generic
                # ensemble-mean anchor is 2 E[mu2_charged]/E[m2] (< 2p/Nc by
                # the cross-term deficit, see docstring)
                anchor_label = "generic anchor 2*E[mu2]/E[m2]"
                anchor_val = (2 * mean_mu2_charged_exact(Nc, p, ens)
                              / mean_m2_exact(Nc, p, ens))
            print(f"    w2_charged = {r['w2_charged']:.4f} "
                  f"+- {r['w2_charged_sem']:.4f}   "
                  f"({anchor_label} = {anchor_val:.3f});  "
                  f"kurt = {r['kurt_charged']:.2f}")
            print(f"    w2_singlet = {r['w2_singlet']:.4f} "
                  f"+- {r['w2_singlet_sem']:.4f};  "
                  f"kurt = {r['kurt_singlet']:.2f}")
            print(f"    charged/singlet frequency-scale ratio = "
                  f"{r['ratio']:.3f} +- {r['ratio_sem']:.3f}")
    print("\n  (Frequency-scale separation in 0+1d, NOT spatial confinement;")
    print("   the confinement-adjacent statement is the E0(Q) curve above.)")

    out_json = here / "qed_campaign.json"
    with open(out_json, "w") as fh:
        json.dump(results, fh, indent=1)
    print(f"\nwrote {out_json}")
    make_figure(results, here / "qed_campaign.png")
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--big", action="store_true",
                    help="larger ensembles + the Nc=10 spectral run")
    main(big=ap.parse_args().big)
