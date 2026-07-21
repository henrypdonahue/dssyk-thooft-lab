#!/usr/bin/env python3
r"""
Exact Wick statistics on the U(N)/Dirac side: adjoint bilinear + singlet mu0.
=============================================================================

Two exact-disorder computations for the charge-conserving complex (Dirac) SYK
model of dirac.py, both in the C-invariant "c_symmetric" ensemble (disjoint
index pairs, real symmetric couplings -- the CP-sharp sub-ensemble):

    H = sum_u x_u h_u,   h_u = c^dag_A c_B + c^dag_B c_A   (u = {A,B},
    A n B = empty, |A| = |B| = k = p/2),   x_u iid N(0, sigma^2),
    sigma^2 = coupling_variance_dirac(Nc, p) = (k!)^2 2^p / Nc^(p-1).

Part 1: Var(B_ij) -- the U(N)-emergence analog of Majorana Var(B_jk)
--------------------------------------------------------------------
The symmetry-violation observable is the adjoint bilinear

    B_ij = Tr(c_i H c^dag_j H)/dim   (i != j),  and the diagonal B_ii.

B_ij is QUADRATIC in the Gaussian x_u, so its disorder statistics reduce by
Wick to trace combinatorics, exactly as in exact_wick.py -- but the
disjointness constraint A n B = empty changes the census.  With

    T_{alpha beta} = Tr(c_i O_alpha c^dag_j O_beta)/dim,   O_(A,B) = c^dag_A c_B

over ordered keys alpha = (A,B), the trace is nonzero only if every mode's
creation and annihilation counts balance.  Per mode the admissible patterns
are: mode i in A only or in C only; mode j in B only or in D only; any other
mode untouched, in A n D, or in B n C (beta = (C,D)).  Solving these
constraints kills the (i in A, j in B) and (i in C, j in D) cases (they force
|D| = k-1) and leaves exactly two families:

    Type I : alpha = (A,B) with i in A, j not in A u B,
             beta  = (B, (A\{i}) u {j});
    Type II: alpha = (A,B) with j in B, i not in A u B,
             beta  = ((B\{j}) u {i}, A).

Each admissible product touches 2k+1 = p+1 modes, every one with FORCED
occupation on the surviving basis states, so |T| = 2^-(p+1); the sign is
uniformly -1 and T_{alpha beta} = T_{beta alpha} (both asserted by literal
enumeration at p = 4 and p = 8 in the tests -- only T^2 and the symmetry
enter the variance).  No admissible beta equals alpha or (B,A), so

    E[B_ij] = 0   exactly (i != j).

In the unordered-pair variables x_u the quadratic form B_ij = x^T M x has M a
perfect matching between the L = C(Nc-2, k-1) C(Nc-k-1, k) pairs containing i
(and avoiding j) and the L pairs containing j (avoiding i), every nonzero
entry -2^-(p+1) in both orders.  Var(x^T M x) = 2 sigma^4 tr(M^2) gives,
EXACTLY,

    Var(B_ij) = 4 sigma^4 C(Nc-2, k-1) C(Nc-k-1, k) / 4^(p+1)
              = (k!)^4 C(Nc-2, k-1) C(Nc-k-1, k) / Nc^(2p-2).

Fixed p, large Nc:  Var -> k (k!)^2 / Nc^(p-1), i.e.

    rms(B_ij) = sqrt(k) k! Nc^-(p-1)/2  (1 + O(1/Nc)),

the SAME power law N^-(p-1)/2 as the Majorana adjoint bilinear
rms(B_jk) ~ (2 p!/sqrt((p-1)!)) N^-(p-1)/2 (exact_wick.py): U(N) symmetry
emergence is polynomially slow at fixed p on the Dirac side too, just with a
smaller coefficient (p = 4: 2 sqrt(2) ~ 2.83 vs 48/sqrt(6) ~ 19.6 -- per
complex mode vs per Majorana mode).

Diagonal: only beta = (B,A) survives, and its trace vanishes unless
i not in A u B (where it is +2^-(p+1)), so

    E[B_ii] = sigma^2 C(Nc-1, k) C(Nc-k-1, k) 2^-(p+1)
            = E[m2] (Nc - p) / (2 Nc),

with E[m2] = E[Tr H^2/dim] = sigma^2 C(Nc,k) C(Nc-k,k) 4^-k the c_symmetric
mean of qed_campaign.py (recomputed here numerically from the h_u, not
imported).  MEAN-VANISHING ARTIFACT: E[B_ii] = 0 exactly at Nc = p (every
coupling then touches every mode; indeed B_ij vanishes identically there,
Var = 0 too).  In Majorana units N_maj = 2 Nc this is the same (1 - 2p/N)
factor as the Majorana p = N/2 artifact of exact_wick.py, with an extra 1/2
from the c-vs-psi normalization (c c^dag + c^dag c = 1 vs psi^2 = 1).

Part 2: exact disorder-averaged mu0 (and mu2 at Nc = 6) of traceless M_2
------------------------------------------------------------------------
The first dynamical singlet is M_2 = -sum_i cdot^dag_i cdot_i,
cdot_i = i[H, c_i] (operator identity, dirac.py).  It is an explicit
quadratic form in the couplings:

    M_2 = sum_{u,v} x_u x_v G_uv,    G_uv = -sum_i [h_u, c_i]^T [h_v, c_i]

(real matrices; G_uv^T = G_vu).  Traceless and CP-projected per COEFFICIENT
(both operations are linear):  Gbar_uv = G_uv - (Tr G_uv/dim) 1, and with the
unitary particle-hole C of dirac.particle_hole,

    G^pm_uv = (Gbar_uv pm C Gbar_uv C^T)/2,

so  M2bar^ch = sum x_u x_v G^ch_uv with ch = CP-odd ("right" channel for
n = 2, target CP = (-1)^(n+1) = -1) or CP-even ("wrong").  Then

    mu0(ch) = Tr((M2bar^ch)^2)/dim

(M_2 is Hermitian, so the moments_pipeline symmetrized-spectral mu0 is
exactly this Frobenius weight) is QUARTIC in x, and Isserlis' three pairings
give the exact ensemble mean

    E[mu0^ch] = sigma^4/dim [ Tr(S_ch^2)
                + sum_{uv ordered} ( Tr(G^ch_uv G^ch_uv) + ||G^ch_uv||_F^2 ) ],
    S_ch = sum_u G^ch_uu,

a sum over n_pairs^2 coupling pairs with every trace evaluated numerically
(n_pairs = C(Nc,k) C(Nc-k,k)/2 = 45, 210, 630 at Nc = 6, 8, 10 for p = 4).
This is rung 15's exact-finite-N half for the n = 2 singlet: the moments
pipeline's ED means reproduced with NO diagonalization and NO sampling.
(The symbolic N = infinity chord computation stays out of scope.)

mu2(ch) = ||[H, M2bar^ch]||_F^2/dim is SEXTIC in x: with
C_{a,uv} = [h_a, Gs^ch_uv] (Gs = the u<->v symmetrized Gbar^ch), Isserlis'
15 pairings collapse by the u<->v symmetry to

    E[mu2^ch] = sigma^6/dim [ 2 S_AA + 4 S_AB + T_PP + 4 T_PQ + 4 T_QQ ],
    S_AA = sum_{a,u,v} Tr(C_{a,uv}^T C_{a,uv}),
    S_AB = sum_{a,u,v} Tr(C_{a,uv}^T C_{u,av}),
    P_a  = sum_u C_{a,uu},   Q_v = sum_a C_{a,av},
    T_PP = sum_a Tr(P_a^T P_a),  T_PQ = sum_a Tr(P_a^T Q_a),
    T_QQ = sum_v Tr(Q_v^T Q_v).

S_AB is an irreducible triple sum: ~n_pairs^3 operator-product traces.  At
Nc = 6 (n = 45, dim 64) that is ~10^5 small dense products -- minutes.  At
Nc = 8 (n = 210, dim 256) it is ~9.3e6 products of 256x256 operators --
order days on this machine -- so mu2 is delivered EXACTLY at Nc = 6 only,
with the Nc = 8 cost documented instead of half-done.

Honesty box
-----------
* The uniform sign -2^-(p+1) and the symmetry T_{alpha beta} = T_{beta alpha}
  are measured facts (enumeration), not derived here; the magnitude and the
  admissible-pair census are derived.  Var and E depend only on the measured
  facts squared/symmetrized, and every closed form is asserted against a
  literal enumeration sharing no assembly code.
* E[mu0], E[mu2] are exact ensemble means; the moments.json comparison
  numbers are ED disorder means with SEM (ddof=1), so agreement is a
  statistical statement (z-scores stored in exact_wick_dirac.json).
* moments.json stores w2 = mu2/(mu0 <H^2>/dim) as a mean of per-instance
  RATIOS; the exact ratio-of-means E[mu2]/(E[mu0] E[m2]) differs at
  O(relVar) ~ 1e-2.  The JSON stores the exact ratio-of-means labeled as
  such; the calibrated mu2 check is the direct slow-test ED ensemble mean.
"""

from __future__ import annotations

from itertools import combinations
from math import comb, factorial, log, sqrt

import numpy as np
import scipy.sparse as sp

from dirac import (annihilators, apply_monomial, coupling_variance_dirac,
                   dirac_hamiltonian, draw_couplings, mn_tower, particle_hole)


# ---------------------------------------------------------------------------
# Part 1: closed forms
# ---------------------------------------------------------------------------
def var_offdiag_csym(Nc: int, p: int) -> float:
    """Var(B_ij) over c_symmetric disorder, i != j:
    4 sigma^4 C(Nc-2, k-1) C(Nc-k-1, k) / 4^(p+1), exact."""
    k = p // 2
    sigma2 = coupling_variance_dirac(Nc, p)
    return 4.0 * sigma2 ** 2 * comb(Nc - 2, k - 1) * comb(Nc - k - 1, k) \
        / 4.0 ** (p + 1)


def rms_offdiag_csym(Nc: int, p: int) -> float:
    """Typical U(N)-symmetry violation sqrt(Var(B_ij))."""
    return sqrt(var_offdiag_csym(Nc, p))


def ln_rms_offdiag_csym(Nc: int, p: int) -> float:
    """ln rms(B_ij) evaluated stably at large (Nc, p) via lgamma."""
    from math import lgamma
    k = p // 2
    ln_sigma2 = 2.0 * lgamma(k + 1) + p * log(2.0) - (p - 1) * log(Nc)
    ln_L = (lgamma(Nc - 1) - lgamma(k) - lgamma(Nc - k)          # C(Nc-2,k-1)
            + lgamma(Nc - k) - lgamma(k + 1) - lgamma(Nc - 2 * k))
    return 0.5 * log(4.0) + ln_sigma2 + 0.5 * ln_L - 0.5 * (p + 1) * log(4.0)


def mean_diag_csym(Nc: int, p: int) -> float:
    """E[B_ii] = sigma^2 C(Nc-1,k) C(Nc-k-1,k) 2^-(p+1)
    = E[m2] (Nc-p)/(2 Nc); vanishes exactly at Nc = p (the U(N) analog of
    the Majorana p = N/2 mean-vanishing artifact)."""
    k = p // 2
    sigma2 = coupling_variance_dirac(Nc, p)
    return sigma2 * comb(Nc - 1, k) * comb(Nc - k - 1, k) / 2.0 ** (p + 1)


def rms_offdiag_large_nc(Nc: int, p: int) -> float:
    """Leading fixed-p law sqrt(k) k! Nc^-(p-1)/2 -- the Dirac-side
    statement to put next to the Majorana 2 p!/sqrt((p-1)!) N^-(p-1)/2."""
    k = p // 2
    return sqrt(k) * factorial(k) * Nc ** (-(p - 1) / 2.0)


# ---------------------------------------------------------------------------
# Part 1: literal enumeration (no assembly code shared with the closed forms)
# ---------------------------------------------------------------------------
def ordered_keys(Nc: int, k: int) -> list:
    """All ordered disjoint (A, B) with |A| = |B| = k -- the stored key set of
    the c_symmetric coupling table."""
    subsets = list(combinations(range(Nc), k))
    return [(A, B) for A in subsets for B in subsets
            if A != B and not set(A) & set(B)]


def chain_trace(Nc: int, ops) -> float:
    """Tr(ops[0] ops[1] ... ops[-1])/dim, each op a (create_set,
    annihilate_set) monomial, composed right-to-left with
    dirac.apply_monomial on all basis states at once."""
    dim = 1 << Nc
    states0 = np.arange(dim, dtype=np.int64)
    cur = states0.copy()
    sign = np.ones(dim)
    valid = np.ones(dim, dtype=bool)
    for (A, B) in reversed(ops):
        v, rows, s = apply_monomial(Nc, A, B, cur)
        valid &= v
        sign = sign * s
        cur = np.where(v, rows, cur)
    good = valid & (cur == states0)
    return float(sign[good].sum()) / dim


def t_matrix(Nc: int, k: int, i: int, j: int, keys=None) -> np.ndarray:
    """T_{alpha beta} = Tr(c_i O_alpha c^dag_j O_beta)/dim over ordered keys
    (j = i gives the diagonal observable's matrix)."""
    if keys is None:
        keys = ordered_keys(Nc, k)
    n = len(keys)
    T = np.zeros((n, n))
    for a, (A, B) in enumerate(keys):
        for b, (C, D) in enumerate(keys):
            T[a, b] = chain_trace(
                Nc, [((), (i,)), (A, B), ((j,), ()), (C, D)])
    return T


def _coupling_covariance(Nc: int, k: int, keys) -> np.ndarray:
    """Covariance E[J_alpha J_beta] = sigma^2 (delta_{alpha beta} +
    delta_{alpha, transpose(beta)}) of the stored ordered-key vector."""
    sigma2 = coupling_variance_dirac(Nc, 2 * k)
    idx = {key: a for a, key in enumerate(keys)}
    n = len(keys)
    C = sigma2 * np.eye(n)
    for a, (A, B) in enumerate(keys):
        C[a, idx[(B, A)]] += sigma2
    return C


def stats_enumerated(Nc: int, p: int, i: int, j: int):
    """(E[B_ij], Var(B_ij)) by the full quadratic form over ordered stored
    keys: B = J^T T J, E = tr(C T), Var = 2 tr((C Ts)^2), Ts = (T + T^T)/2.
    Exponentially expensive; small (Nc, p) only."""
    k = p // 2
    keys = ordered_keys(Nc, k)
    T = t_matrix(Nc, k, i, j, keys)
    C = _coupling_covariance(Nc, k, keys)
    Ts = 0.5 * (T + T.T)
    mean = float(np.trace(C @ T))
    var = 2.0 * float(np.trace(C @ Ts @ C @ Ts))
    return mean, var


def offdiag_census(Nc: int, p: int, i: int, j: int):
    """(n_nonzero, values, symmetric) of the enumerated T matrix -- checks the
    derived admissible-pair count 2 L, the uniform value -2^-(p+1), and
    T_{alpha beta} = T_{beta alpha}."""
    k = p // 2
    T = t_matrix(Nc, k, i, j)
    nz = np.abs(T) > 1e-14
    values = sorted(set(np.round(T[nz], 12)))
    return int(nz.sum()), values, bool(np.allclose(T, T.T, atol=1e-14))


def instance_b(Nc: int, p: int, couplings: dict, i: int, j: int) -> float:
    """Single-instance B_ij = Tr(c_i H c^dag_j H)/dim by dense ED-side
    algebra (independent of the T-matrix route)."""
    H = dirac_hamiltonian(Nc, p, couplings=couplings)
    cs = annihilators(Nc)
    dim = 1 << Nc
    val = np.trace(cs[i] @ H @ cs[j].conj().T @ H) / dim
    assert abs(val.imag) < 1e-12
    return float(val.real)


# ---------------------------------------------------------------------------
# Part 2: pair operators and the M_2 quadratic form
# ---------------------------------------------------------------------------
def unordered_pairs(Nc: int, k: int) -> list:
    """Unordered disjoint pairs (A, B) in the iteration order of
    dirac.draw_couplings (A before B in lex subset order), so that
    couplings[(A, B)] is exactly the basis variable x_u of pair u."""
    subsets = list(combinations(range(Nc), k))
    pairs = []
    for ia, A in enumerate(subsets):
        for B in subsets[ia:]:
            if not set(A) & set(B):
                pairs.append((A, B))
    return pairs


def monomial_sparse(Nc: int, A, B) -> sp.csr_matrix:
    """Sparse real matrix of c^dag_A c_B from dirac.apply_monomial."""
    dim = 1 << Nc
    states = np.arange(dim, dtype=np.int64)
    valid, rows, signs = apply_monomial(Nc, A, B, states)
    return sp.csr_matrix((signs[valid], (rows[valid], states[valid])),
                         shape=(dim, dim))


def pair_operator(Nc: int, A, B) -> sp.csr_matrix:
    """h_u = c^dag_A c_B + c^dag_B c_A (real symmetric)."""
    return monomial_sparse(Nc, A, B) + monomial_sparse(Nc, B, A)


def _annihilators_sparse(Nc: int) -> list:
    return [monomial_sparse(Nc, (), (i,)) for i in range(Nc)]


def _ph_perm_sign(Nc: int):
    """(perm, sign) of the unitary particle-hole C = dirac.particle_hole as a
    signed permutation: C e_s = sign[s] e_perm[s]."""
    U = particle_hole(Nc)
    Ur = U.real
    perm = np.argmax(np.abs(Ur), axis=0)
    sign = Ur[perm, np.arange(Ur.shape[0])]
    assert np.abs(np.abs(sign) - 1.0).max() < 1e-12
    return perm, sign


def _ph_conj_sparse(X: sp.spmatrix, perm, sign) -> sp.csr_matrix:
    """C X C^T for sparse real X: (CXC^T)[perm[a], perm[b]] = s_a s_b X[a,b]."""
    coo = X.tocoo()
    return sp.csr_matrix((coo.data * sign[coo.row] * sign[coo.col],
                          (perm[coo.row], perm[coo.col])), shape=X.shape)


def pair_commutators(Nc: int, k: int, pairs=None):
    """D[u][i] = [h_u, c_i] (sparse real) for every unordered pair u."""
    if pairs is None:
        pairs = unordered_pairs(Nc, k)
    cs = _annihilators_sparse(Nc)
    return [[(h @ c - c @ h).tocsr() for c in cs]
            for h in (pair_operator(Nc, A, B) for (A, B) in pairs)]


def g_matrix(D_u, D_v) -> sp.csr_matrix:
    """G_uv = -sum_i D_ui^T D_vi (real; G_uv^T = G_vu)."""
    G = None
    for Du_i, Dv_i in zip(D_u, D_v):
        term = Du_i.T @ Dv_i
        G = term if G is None else G + term
    return (-G).tocsr()


def m2_from_pairs(Nc: int, p: int, couplings: dict) -> np.ndarray:
    """Instance M_2 = sum_uv x_u x_v G_uv (dense real), for the identity test
    against dirac.mn_tower."""
    k = p // 2
    pairs = unordered_pairs(Nc, k)
    D = pair_commutators(Nc, k, pairs)
    x = np.array([couplings[pair] for pair in pairs])
    dim = 1 << Nc
    M2 = np.zeros((dim, dim))
    for u in range(len(pairs)):
        for v in range(len(pairs)):
            M2 += x[u] * x[v] * g_matrix(D[u], D[v]).toarray()
    return M2


# ---------------------------------------------------------------------------
# Part 2: exact E[mu0] of traceless CP-resolved M_2 by Isserlis
# ---------------------------------------------------------------------------
def mu0_exact_csym(Nc: int, p: int, progress: bool = False) -> dict:
    """Exact disorder mean of mu0 = Tr((M2bar^ch)^2)/dim for the CP-odd
    ("right" for n = 2) and CP-even ("wrong") channels and their sum, by the
    three Isserlis pairings over all coupling pairs.  Cost ~ n_pairs^2 sparse
    operator products."""
    k = p // 2
    dim = 1 << Nc
    pairs = unordered_pairs(Nc, k)
    n = len(pairs)
    D = pair_commutators(Nc, k, pairs)
    perm, sign = _ph_perm_sign(Nc)
    sigma2 = coupling_variance_dirac(Nc, p)
    eye = sp.identity(dim, format="csr")

    # channel accumulators: [odd, even]
    t23 = [0.0, 0.0]                       # sum of Tr(G G) + ||G||_F^2 terms
    S_ch = [np.zeros((dim, dim)), np.zeros((dim, dim))]
    for u in range(n):
        if progress and u % 20 == 0:
            print(f"    mu0 pair row {u}/{n}", flush=True)
        for v in range(u, n):
            G = g_matrix(D[u], D[v])
            tr_uv = G.diagonal().sum() / dim
            Gb = (G - tr_uv * eye).tocsr()
            Gb2 = _ph_conj_sparse(Gb, perm, sign)
            w = 1.0 if u == v else 2.0     # ordered-pair multiplicity
            for ch, Gch in enumerate((0.5 * (Gb - Gb2), 0.5 * (Gb + Gb2))):
                # Tr(Gch_uv Gch_uv) + Tr(Gch_uv Gch_vu); Gch_vu = Gch_uv^T
                tr_gg = Gch.multiply(Gch.T).sum()
                tr_ggT = Gch.multiply(Gch).sum()
                t23[ch] += w * (tr_gg + tr_ggT)
                if u == v:
                    S_ch[ch] += Gch.toarray()
    out = {}
    for ch, tag in enumerate(("odd", "even")):
        t1 = float(np.sum(S_ch[ch] * S_ch[ch].T))     # Tr(S^2)
        out[tag] = sigma2 ** 2 * (t1 + t23[ch]) / dim
    out["total"] = out["odd"] + out["even"]
    out["n_pairs"] = n
    return out


def m2_mean_exact(Nc: int, p: int) -> float:
    """E[Tr H^2/dim] = sigma^2 sum_u Tr(h_u^2)/dim, numerically from the pair
    operators (matches the qed_campaign.py closed form
    sigma^2 C(Nc,k) C(Nc-k,k) 4^-k; recomputed, not imported)."""
    k = p // 2
    sigma2 = coupling_variance_dirac(Nc, p)
    dim = 1 << Nc
    tot = 0.0
    for (A, B) in unordered_pairs(Nc, k):
        h = pair_operator(Nc, A, B)
        tot += h.multiply(h.T).sum()
    return sigma2 * tot / dim


# ---------------------------------------------------------------------------
# Part 2: exact E[mu2] at small Nc (sextic Isserlis, ~n_pairs^3 traces)
# ---------------------------------------------------------------------------
def _gs_channel_dense(Nc: int, k: int, pairs, D, perm, sign):
    """Dense symmetrized channel coefficients Gs^ch_uv = (Gbar_uv +
    Gbar_vu)/2 projected on the two CP channels; returned as two dicts
    {(u, v): array} with u <= v (odd first, even second)."""
    dim = 1 << Nc
    n = len(pairs)
    eye = sp.identity(dim, format="csr")
    gs = [{}, {}]
    for u in range(n):
        for v in range(u, n):
            G = g_matrix(D[u], D[v])
            tr_uv = G.diagonal().sum() / dim
            Gb = (G - tr_uv * eye).tocsr()
            Gs = 0.5 * (Gb + Gb.T)               # u<->v symmetrized
            Gs2 = _ph_conj_sparse(Gs, perm, sign)
            gs[0][(u, v)] = (0.5 * (Gs - Gs2)).toarray()
            gs[1][(u, v)] = (0.5 * (Gs + Gs2)).toarray()
    return gs


def mu2_exact_csym(Nc: int, p: int, progress: bool = True) -> dict:
    """Exact disorder mean of mu2(ch) = ||[H, M2bar^ch]||_F^2/dim by the 15
    Isserlis pairings (collapsed by u<->v symmetry; see module docstring).
    Feasible at Nc = 6, p = 4 (~10^5 dense 64x64 products, ~1 min); the
    Nc = 8 cost (~9.3e6 products of 256x256) is documented as out of reach,
    not attempted."""
    k = p // 2
    dim = 1 << Nc
    pairs = unordered_pairs(Nc, k)
    n = len(pairs)
    D = pair_commutators(Nc, k, pairs)
    perm, sign = _ph_perm_sign(Nc)
    sigma2 = coupling_variance_dirac(Nc, p)
    hs = [pair_operator(Nc, A, B).toarray() for (A, B) in pairs]
    gs_ch = _gs_channel_dense(Nc, k, pairs, D, perm, sign)

    def gs_get(gs, u, v):
        return gs[(u, v)] if u <= v else gs[(v, u)]

    out = {}
    for ch, tag in ((0, "odd"), (1, "even")):
        gs = gs_ch[ch]
        # cheap terms: P_a = [h_a, sum_u Gs_uu],  Q_v = sum_a [h_a, Gs_av]
        Sdiag = np.zeros((dim, dim))
        for u in range(n):
            Sdiag += gs[(u, u)]
        P = [hs[a] @ Sdiag - Sdiag @ hs[a] for a in range(n)]
        Q = []
        for v in range(n):
            acc = np.zeros((dim, dim))
            for a in range(n):
                X = gs_get(gs, a, v)
                acc += hs[a] @ X - X @ hs[a]
            Q.append(acc)
        T_PP = sum(float(np.sum(P[a] * P[a])) for a in range(n))
        T_PQ = sum(float(np.sum(P[a] * Q[a])) for a in range(n))
        T_QQ = sum(float(np.sum(Q[v] * Q[v])) for v in range(n))
        # NB Tr(X^T Y) = sum(X * Y) elementwise for real X, Y.

        # S_AA = sum_{a,u,v} ||C_{a,uv}||_F^2  (C = [h_a, Gs_uv])
        S_AA = 0.0
        for a in range(n):
            if progress and a % 10 == 0:
                print(f"    mu2[{tag}] S_AA row {a}/{n}", flush=True)
            ha = hs[a]
            for u in range(n):
                for v in range(u, n):
                    X = gs[(u, v)]
                    Cm = ha @ X - X @ ha
                    w = 1.0 if u == v else 2.0
                    S_AA += w * float(np.sum(Cm * Cm))
        # S_AB = sum_{a,u,v} Tr(C_{a,uv}^T C_{u,av})
        S_AB = 0.0
        for a in range(n):
            if progress and a % 10 == 0:
                print(f"    mu2[{tag}] S_AB row {a}/{n}", flush=True)
            ha = hs[a]
            for u in range(n):
                hu = hs[u]
                for v in range(n):
                    X = gs_get(gs, u, v)
                    Y = gs_get(gs, a, v)
                    C1 = ha @ X - X @ ha
                    C2 = hu @ Y - Y @ hu
                    S_AB += float(np.sum(C1 * C2))
        out[tag] = sigma2 ** 3 * (2 * S_AA + 4 * S_AB
                                  + T_PP + 4 * T_PQ + 4 * T_QQ) / dim
    out["n_pairs"] = n
    return out


# ---------------------------------------------------------------------------
# ED-side instance moments (independent eigenbasis route, for cross-checks)
# ---------------------------------------------------------------------------
def instance_mu_ed(Nc: int, p: int, seed: int):
    """Per-instance (mu0, mu2) of the CP-resolved traceless M_2 channels via
    the moments_pipeline eigenbasis conventions (delegated to
    moments_pipeline.instance_moments -- the established ED route this
    module's Wick sums are checked against)."""
    from moments_pipeline import instance_moments
    inst = instance_moments(Nc, p, seed, n_list=(2,))
    return {"odd": (inst[(2, "right")]["mu0"], inst[(2, "right")]["mu2"]),
            "even": (inst[(2, "wrong")]["mu0"], inst[(2, "wrong")]["mu2"]),
            "escale2": inst[(2, "right")]["escale"] ** 2}
