#!/usr/bin/env python3
r"""
The singlet tower at N = infinity: every A_n is a stack of matter chords.
=========================================================================

full_channel_chord.py gives the exact moments of the FULL A_2 channel at
N = infinity, and leaves two things open: the conserved H^2 admixture
that mean subtraction does not remove, and every n >= 3.  This module
settles the first, and the second through n = 8 (a general-n proof is
open).  The answer has one shape at every level computed (n = 2..8):

    A_n  =  (a polynomial in H)  +  sum_{m=2..n}  F_{n,m}(H_L, H_R) . O_m ,

where O_m is an ordinary matter chord of dimension Delta = m -- an
H-chord that crosses it weighs q^m -- the O_m of different m are
independent, and F_{n,m} is a polynomial of degree n - m that acts by
left and right multiplication with H.  One new matter chord enters at
each level: O_n first appears in A_n.

So the paper's singlet tower (their Sec. 4.3; here its Majorana form) is,
at N = infinity and fixed lambda, a triangular change of basis of a tower
of matter chords of INTEGER dimension 2, 3, 4, ...  Everything about
those is known in closed form: their two-point functions are the
standard DSSYK matter correlators, so the complete spectral density of
each A_n channel follows.

THE EXACT OPERATOR IDENTITY (any N, any even p, per instance)
-------------------------------------------------------------
Write H_i for the part of H that contains flavor i, and for a set T of
slots let W_T be the word H^n with H replaced by H_i on the slots in T,
summed over i:

    W_T = sum_i  prod_{t=1..n} ( H_i if t in T else H ).

Then

    A_n = sum_{T subset [n], |T| >= 2}
              (-2)^|T| (-1)^(n - max T) C(n-1, max T - 1)  W_T .        (1)

n = 2 is A_2 = 4 sum_i H_i^2 (conserved_points.py).  n = 3 is
A_3 = 4 (W_13 + W_23 - 2 W_12) - 8 W_123.  The derivation is two lines:
the spin-trace map S(X) = sum_i psi_i X psi_i multiplies a Clifford
monomial of grade g by N - 2g (p even), the grade of a product of m
Hamiltonian monomials is |c_1 xor ... xor c_m| = sum_T (-2)^(|T|-1)
|cap_T c_t| by inclusion-exclusion, and A_n = sum_j (-1)^j C(n,j)
S(H^(n-j)) H^j.  The identity is asserted against dense matrices for
n <= 5 (syk-self-averaging/test_conserved_points.py).

THE BUILDING BLOCKS
-------------------
Only the pairs survive the double-scaled limit.  A triple intersection
of three distinct random p-subsets is O(1/p).  With

    Omega_k = sum_i H_i H^k H_i ,     kappa_ab = 4 (-1)^(n-b) C(n-1, b-1),

(1) becomes

    A_n  ~  sum_{a<b} kappa_ab  H^(a-1) Omega_(b-a-1) H^(n-b)            (2)

up to exactly conserved operators and o(1).  Omega_0 = A_2 / 4.

THE CHORD RULE FOR Omega_k
--------------------------
Omega_k is the word H^(k+2) with the first and last Hamiltonian
monomials weighted by the number X of Majorana indices they share.  In
the double-scaled limit the pairwise intersections of distinct chords
are independent Poisson(lambda/2) variables, mu = lambda/2, and a
crossing of two chords carries (-1)^X.  So with q = e^(-lambda):

  * every ordinary crossing weighs q;
  * the weighted pair carries E[X (-1)^(sX)] = mu or -mu q (s = 1 iff
    its two chords cross);
  * if two insertions weight the SAME pair of chords, the pair carries
    E[X^2 (-1)^(sX)] = mu + mu^2 or q (mu^2 - mu);
  * if the first and last slots are contracted WITH EACH OTHER, X = p.
    The enclosing chord then crosses every enclosed chord, and the exact
    crossing average depends on the grade of the enclosed product:
    p (1 - 2p/N)^grade = p q^k [1 + (2 lambda / p) sum_pairs X + ...].
    The O(p) term is p q^k times a polynomial in H: conserved, dropped.
    The O(1) remainder is 2 lambda times an X-weight on each pair of
    enclosed chords.  It is not optional.  It first matters at k = 2,
    and the exact finite-N cross moment <ad^2(A_4) A_2> tells the two
    apart by a wide margin (gates below).

THE DECOMPOSITION
-----------------
Modulo polynomials in H, and writing O^_m for the un-normalized chord
(norm_m = <O^_m O^_m>),

    Omega_0 = O^_2
    Omega_1 = [q / (1+q^2)] {H, O^_2}  +  O^_3
    Omega_2 = a (H^2 O^_2 + O^_2 H^2) + b H O^_2 H + c O^_2
              + [q / (1-q+q^2)] {H, O^_3}  +  O^_4
    Omega_3 .. Omega_6: computed by decompose(); new chords O^_5 .. O^_8

with [5]_q = 1+q+q^2+q^3+q^4 and

    a = q^2 (1+q+q^2) / ((1+q^2) [5]_q),      b = q^2 (1+q)^2 / ((1+q^2) [5]_q),
    c = [1 + q - q^2 - q^3 + 2 lambda q^2 [5]_q] / [5]_q ,

    norm_2 = lambda (1 - q^2 + 2 lambda q) / (2 (1+q))
    norm_3 = lambda (1+q) [1 - q^6 + 2 lambda q^2 (2 - q + 2q^2)]
             / (2 (1+q^2) (1+q+q^2))
    norm_4 = lambda (1-q)^2 (1+q)^2 [ (1 - q^6)(1+q+2q^2+q^3+q^4)
             + 2 lambda q^3 (2 + q + 2q^2 + q^3 + 2q^4) ]
             / (2 (1+q^2) (1-q+q^2) [5]_q) .

The first sub-leading coefficient has a one-line reading, and it holds
at every level computed.  A plain chord end that sits inside a bundle of
m strands crosses only the strands on its way out.  q/(1+q^2) solves
(q, q) = x [(1, q^2) + (q^2, 1)]: one strand crossed either way, against
none or both for a chord end placed to the left or right of a pointlike
O_2.  In general the O_(k+1) sector of Omega_k is

    [ (q + q^2 + ... + q^k) / (1 + q^(k+1)) ] {H, O^_(k+1)} ,   k = 1, ..., 6.

WHAT IS PROVEN AND WHAT IS VERIFIED
-----------------------------------
k = 0 is a theorem.  In <H^A Omega_0 H^B Omega_0> sort the matchings by
whether both chords of one insertion end on the other insertion.

  (I)  They do.  The two chords form a bundle.  Every H-chord that
       separates the insertions crosses both strands and weighs q^2.
       The bundle's own weight is (mu + mu^2) + q (mu^2 - mu), nested
       plus crossed.  Total: [mu (1-q) + mu^2 (1+q)] G_2(A, B), with G_2
       the two-point word of a Delta = 2 matter chord.
  (II) They do not.  Then at each insertion, swapping which of its two
       adjacent slots carries which chord is an involution.  It flips
       the mutual crossing of the two chords and changes nothing else,
       because no node sits between adjacent slots.  It turns the weight
       mu into -mu q.  The same involution on a plain H^2 turns 1 into
       q.  So class II equals [mu (1-q)/(1+q)]^2 times the same class of
       <H^A (H^2 - 1) H^B (H^2 - 1)>.
  Class I of that plain word is (1+q) G_2.  Hence, exactly,

    <H^A Omega_0 H^B Omega_0> = [mu(1-q)/(1+q)]^2 <H^(A+B) (H^2-1)^2>
                                + norm_2 G_2(A, B),
    norm_2 = mu(1-q) + mu^2(1+q) - mu^2 (1-q)^2/(1+q) .

For the full channel A_2 = 4 Omega_0 this is

    A_2 - <A_2> = [2 lambda (1-q)/(1+q)] (H^2 - 1) + nu O_2 ,
    nu^2 = 8 lambda (1 - q^2 + 2 lambda q) / (1+q) ,

    mu_0 = 4 lambda^2 (1-q)^2/(1+q)  [stationary]  +  nu^2  [dynamical],
    w2_dyn = mu_2 / nu^2 = 2 (1 - q^2)       exactly, at every lambda.

That is the matter-chord law mu_2 = 2(1 - q^Delta) at Delta = 2.  The
raw width full_channel_chord.py quotes, 1.316906 at lambda_chord = 2, is
this 1.963369 diluted by the conserved H^2 weight (32.9 % of mu_0 there).

k = 1..6 are exact identities verified by computer algebra, not proofs.
decompose() solves each sector's kernel LINEARLY from the
cross-correlators with the lower blocks and then checks every remaining
two-point word.  The residuals are identically zero: symbolically in
(lambda, q) for k <= 2, and in exact rational arithmetic at generic
points for k = 3..6 (test_singlet_tower.py; lambda and q are
independent there, as the identities are polynomial in both; k = 5 and 6,
the levels n = 7 and 8, are singlet_tower_deep.json).  n_checked counts
each word (A, B) and its mirror (B, A), which give the same equation, and
words already used in the sector solves: the independent checks are 37
at n <= 6, 55 at n = 7 and 78 at n = 8 per point.  A zero residual tests
the STRUCTURE (a polynomial in H plus dressed matter chords), not the
values: it stays zero under several wrong chord rules (remainder off or
rescaled, other pair weights), because the strands between the two
insertions always form a matter chord.  The values (norms, dressings,
shares) rest on the derivation of the rule, the n = 2 theorem and the
finite-N gates at n = 3, 4 (ERRATA E26).  A proof for general k is
open.

WHAT IT SAYS ABOUT THE TOWER
----------------------------
  * The channel content of A_n at N = infinity is a conserved polynomial
    in H plus matter chords of dimension 2..n.  No discrete lines: each
    sector is the smooth Askey-Wilson continuum of chord.py, weighted by
    a polynomial.  A_2's density is monotone; A_3's and A_4's have the
    single broad maximum of a derivative kernel and no nodes.
  * The new chord at each level is faint.  At lambda_chord = 2 the
    dynamical weight of A_3 is 85 % O_2 and 15 % O_3; O_4 is 1.6 % of
    A_4, O_5 0.2 % of A_5, O_6 0.02 % of A_6, O_7 0.002 % of A_7, O_8
    0.0002 % of A_8.  Over lambda_chord in [0.05, 4] the weight of
    A_4..A_6 outside O_2 and O_3 stays between 0.9 % and 5.4 % (at
    lambda_chord = 2 it is 6.9 % for A_7 and 9.1 % for A_8).
  * The norms are graded: norm_m scales as lambda^(2 floor(m/2)),

        norm_2 -> lambda^2,  norm_3 -> 2 lambda^2,  norm_4 -> (52/5) lambda^4,
        norm_5 -> (88/7) lambda^4,  norm_6 -> (4120/7) lambda^6,
        norm_7 -> 720 lambda^6,     norm_8 ~ 1e5 lambda^8

    (the first three from the closed forms, m = 5..7 the exact limits of
    conformal_tower.py, m = 8 measured; m = 7, 8 were predicted by the law
    before they were computed).
  * Norms are not shares, because the dressings carry powers of the
    frequency.  Toward small lambda (share_scaling, exact arithmetic
    down to lambda = 1/1280) the EVEN chords fade out of every level
    above the second, linearly in lambda, and the odd ones stay:

        A_3 -> 99.65 % O_3;   A_4 -> 99.61 % O_3;
        A_5 -> 99.19 % O_3 + 0.53 % O_5;   A_6 -> 98.96 % O_3 + 0.80 % O_5;
        A_7 -> 98.40 % O_3 + 1.38 % O_5 + 0.006 % O_7;
        A_8 -> 97.80 % O_3 + 1.99 % O_5 + 0.014 % O_7

    (values at lambda = 1/1280, the rest in the fading even chords; A_3
    and A_4 tend to O_3 alone, and O_8 fades linearly like the other even
    chords).  A_2 is O_2 at every lambda.  So the dimensions that survive
    the semiclassical limit are 2, 3, 5, 7 -- the pattern of the large-q
    SYK bilinear spectrum h = 2, 3, 5, 7, ... (Maldacena-Stanford
    1604.07818).  O_7 surviving and O_8 fading was the pattern's
    prediction for n = 7, 8 before they were computed.  conformal_tower.py
    takes the limit through n = 7: with the chord rule's grade remainder
    off the chords are weight-m primaries of a conformal OPE of the
    bilocal (with it on, plus one Liouville rung), the survival follows from
    the norm law, and the shares above have exact rational limits
    (A_5 -> 2047/2058 O_3 + 11/2058 O_5).
  * The O_2 kernel of A_3 is (3/2)(E - E') plus an anticommutator
    admixture -(1-q)^2 / (3 (1+q^2)) of it: A_3 = (3/2) [H, A_2] + 4 O^_3
    + ..., a total derivative plus one new chord.  (A_3 - A_3^dag =
    3 [H, A_2] is an exact operator identity.)  The derivative is why
    O_2's share of A_3 dies at small lambda.

GATES (test_singlet_tower.py, singlet_tower.json)
-------------------------------------------------
  * The n = 2 theorem against the exact finite-N full channel with its
    H^2 component removed (exact rationals, the committed
    full_channel_chord.json rows): on lambda_chord = 2 the 1/p limits
    are w2_dyn = 1.96336, stationary 10.5365, dynamical 21.4638
    against 1.963369, 10.536378, 21.463622; on lambda_chord = 1,
    1.7293 against 1.729329.
  * The same split in dense ED, per instance, removing the full
    commutant of H (projected_channel.spectral_split): the stationary
    weight is the H^2 projection, and at p = 6, N = 16..20 the dynamical
    width is 2.03, 2.03, 2.01 against the N = infinity 2.00, while the
    raw width sits at 0.29-0.41.
  * Along fixed p the identity-subtracted width has a maximum (N = 24 at
    p = 4, full_channel_laws.py).  With the H^2 component removed it has
    none: that maximum is conserved weight leaving the denominator.
  * mu_4 of A_2, mu_2 and mu_4 of A_3, mu_2 of A_4, and the cross moment
    <ad^2(A_4) A_2^dag> against the exact finite-N engine on
    double-scaling lines.  The cross moment is the sharp one for the
    remainder rule: 155.03 with it and 150.65 without at
    lambda_chord = 1, and the exact series (p <= 20) extrapolates to
    155.027 (the fits of order 3 to 5 agree to 0.06).
  * The spectral density built from the decomposition and the
    closed-form vertex of review 2407.09396 (2.30) integrates back to
    the moments.

SCOPE.  Majorana A_n, strict N = infinity at fixed lambda_chord.  The
Dirac M_n tower and its CP resolution are not covered: there
c_i c_i^dag = 1 - n_i replaces psi_i^2 = 1, and the flavor sum leaves an
occupation factor between the two halves of M_2 instead of removing the
fermion marks.  A matter-chord continuum in boost frequency
is not a mass spectrum; duality/BOOST.md is the calibration for that
step.

Output: singlet_tower.json.  `python singlet_tower.py` prints the
headline numbers in seconds; `--run` regenerates the artifact (the
finite-N gates dominate: about three CPU-hours, about half of it in the
ten K = 10 points (n3_k4 and n4_k2, p <= 12); `--resume` keeps the
points already computed).
`--deep 5 6` writes singlet_tower_deep.json, the levels n = 7 and 8: the
K = 16 and 18 word tables are built in parallel (build_word_tables, about
15 minutes on 6 cores, nearly all of it K = 18), then the same exact
checks and scalings (about 10 more minutes).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# chord diagrams
# ---------------------------------------------------------------------------
def _matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for k in range(1, len(items)):
        rest = items[1:k] + items[k + 1:]
        for rem in _matchings(rest):
            yield [(a, items[k])] + rem


def _diagrams(K: int, first: int = None):
    """Every perfect matching of K nodes on a line, as (of, cross): the
    chord index of each node, and the set of crossing chord-index pairs.
    With `first`, only the matchings that pair node 0 with node `first`
    (the unit of work of build_word_tables)."""
    if first is None:
        gen = _matchings(list(range(K)))
    else:
        gen = ([(0, first)] + rem for rem in
               _matchings([x for x in range(1, K) if x != first]))
    for m in gen:
        of = [0] * K
        for ci, (x, y) in enumerate(m):
            of[x] = ci
            of[y] = ci
        cross = set()
        for i in range(len(m)):
            a, b = m[i]
            for j in range(i + 1, len(m)):
                c, d = m[j]
                if (a < c < b) != (a < d < b):
                    cross.add((i, j))
        yield of, cross


def _pair(x: int, y: int) -> tuple:
    return (x, y) if x < y else (y, x)


def _omega_terms(sc) -> tuple:
    """The X-weighted terms of one Omega_k insertion whose slots sit on
    the chords `sc` (length k + 2).  Each term is (r, pair): an X-weight
    on that pair of distinct chords, times (2 lambda)^r.

    First and last slot on distinct chords: the pair itself, r = 0.
    First and last contracted with each other: the O(1) remainder of the
    X = p term, an X-weight on every pair of enclosed slots that sit on
    two distinct chords leaving the enclosure, r = 1."""
    if sc[0] != sc[-1]:
        return ((0, _pair(sc[0], sc[-1])),)
    inner = sc[1:-1]
    seen: dict = {}
    for c in inner:
        seen[c] = seen.get(c, 0) + 1
    ext = [c for c in inner if seen[c] == 1]
    return tuple((1, _pair(ext[i], ext[j]))
                 for i in range(len(ext)) for j in range(i + 1, len(ext)))


# The tables below are pure combinatorics: they do not depend on lambda or
# q, so one enumeration serves floats, exact rationals and sympy symbols.

@lru_cache(maxsize=None)
def omega_word_table(k1: int, k2: int, A: int, B: int) -> tuple:
    """<H^A Omega_k1 H^B Omega_k2> as a count table.  Keys (kind, nx, r):
    nx crossings in the diagram, r powers of 2 lambda, and
    kind 0 = both insertions weight the same pair, not crossed
         1 = same pair, crossed
         2 = two different pairs, an even number of them crossed
         3 = two different pairs, an odd number crossed."""
    key = (k1, k2, A, B)
    if key in _PRELOADED:
        return _PRELOADED[key]
    if (A + B + k1 + k2) % 2:
        return ()
    acc: dict = {}
    spec = _word_slots(key)
    for of, cross in _diagrams(spec[0]):
        _count_word(acc, of, cross, spec)
    return tuple(sorted(acc.items()))


# Tables computed ahead of time by build_word_tables (identical to what
# omega_word_table would compute; the deep levels need them in bulk).
_PRELOADED: dict = {}


def _word_slots(key) -> tuple:
    """(K, slots of the first insertion, slots of the second)."""
    k1, k2, A, B = key
    n1, n2 = k1 + 2, k2 + 2
    K = A + B + n1 + n2
    return K, range(A, A + n1), range(A + n1 + B, K)


def _count_word(acc: dict, of, cross, spec) -> None:
    """Add one diagram's terms to the count table of one word."""
    _, s1, s2 = spec
    t1 = _omega_terms([of[s] for s in s1])
    if not t1:
        return
    t2 = _omega_terms([of[s] for s in s2])
    if not t2:
        return
    nx = len(cross)
    for r1, e in t1:
        ce = e in cross
        for r2, f in t2:
            if e == f:
                key = (1 if ce else 0, nx, r1 + r2)
            else:
                key = (3 if ce != (f in cross) else 2, nx, r1 + r2)
            acc[key] = acc.get(key, 0) + 1


def _branch_tables(task) -> list:
    """Worker: the partial tables of several words of one size K from
    the diagrams that pair node 0 with node `first`."""
    K, first, keys = task
    specs = [_word_slots(k) for k in keys]
    accs = [dict() for _ in keys]
    for of, cross in _diagrams(K, first):
        for acc, spec in zip(accs, specs):
            _count_word(acc, of, cross, spec)
    return accs


def needed_word_keys(kmax: int, check_len: int) -> dict:
    """Every omega_word_table key that decompose(kmax, check_len) asks
    for (its loops, mirrored), grouped by diagram size K."""
    words = TowerChord._words
    keys = set()
    for k in range(kmax + 1):
        for m in range(2, k + 2):
            k1 = m - 2
            for (A, B) in words((k1 + k) % 2, max(check_len, k + 4 - m)):
                keys.add((k1, k, A, B))
        for (A, B) in words(0, check_len):
            keys.add((k, k, A, B))
        for k1 in range(k):
            for (A, B) in words((k1 + k) % 2, check_len):
                keys.add((k1, k, A, B))
    by_K: dict = {}
    for key in keys:
        by_K.setdefault(_word_slots(key)[0], []).append(key)
    return {K: sorted(v) for K, v in sorted(by_K.items())}


def build_word_tables(keys_by_K: dict, workers: int = None,
                      progress: bool = True) -> dict:
    """Compute omega_word_table for many keys at once: one enumeration of
    the K-node diagrams per size, split over node 0's partner and run in
    parallel.  The result is installed in _PRELOADED and returned; it is
    identical to the one-at-a-time tables (test_singlet_tower.py)."""
    import multiprocessing as mp
    workers = workers or max(1, mp.cpu_count() - 1)
    out = {}
    ctx = mp.get_context("spawn")
    for K, keys in sorted(keys_by_K.items()):
        keys = [k for k in keys if k not in _PRELOADED]
        if not keys:
            continue
        t0 = time.time()
        tasks = [(K, first, keys) for first in range(1, K)]
        with ctx.Pool(processes=min(workers, len(tasks))) as pool:
            parts = pool.map(_branch_tables, tasks)
        for i, key in enumerate(keys):
            acc: dict = {}
            for part in parts:
                for kk, v in part[i].items():
                    acc[kk] = acc.get(kk, 0) + v
            out[key] = tuple(sorted(acc.items()))
        _PRELOADED.update(out)
        if progress:
            print(f"  K = {K}: {len(keys)} word tables in "
                  f"{time.time() - t0:.1f} s", flush=True)
    return out


@lru_cache(maxsize=None)
def omega_one_table(k: int, A: int) -> tuple:
    """<H^A Omega_k>.  Keys (crossed, nx, r)."""
    n1 = k + 2
    K = A + n1
    if K % 2:
        return ()
    acc: dict = {}
    for of, cross in _diagrams(K):
        t1 = _omega_terms([of[s] for s in range(A, K)])
        nx = len(cross)
        for r1, e in t1:
            key = (1 if e in cross else 0, nx, r1)
            acc[key] = acc.get(key, 0) + 1
    return tuple(sorted(acc.items()))


@lru_cache(maxsize=None)
def h_table(k: int) -> tuple:
    """<H^k>: counts by crossing number (Touchard-Riordan)."""
    if k % 2:
        return ()
    acc: dict = {}
    for _, cross in _diagrams(k):
        acc[len(cross)] = acc.get(len(cross), 0) + 1
    return tuple(sorted(acc.items()))


@lru_cache(maxsize=None)
def matter_table(A: int, B: int) -> tuple:
    """<H^A O H^B O> for one matter chord: counts by (H-H crossings,
    H-chords that cross the matter chord)."""
    if (A + B) % 2:
        return ()
    acc: dict = {}
    for m in _matchings(list(range(A + B))):
        nx = 0
        nm = 0
        for i in range(len(m)):
            a, b = m[i]
            if (a < A) != (b < A):
                nm += 1
            for j in range(i + 1, len(m)):
                c, d = m[j]
                if (a < c < b) != (a < d < b):
                    nx += 1
        acc[(nx, nm)] = acc.get((nx, nm), 0) + 1
    return tuple(sorted(acc.items()))


# ---------------------------------------------------------------------------
# small exact linear algebra (floats, Fractions or sympy)
# ---------------------------------------------------------------------------
def _is_sympy(x) -> bool:
    return type(x).__module__.startswith("sympy")


def _solve_linear(rows, rhs):
    """Solve rows @ x = rhs for a full-column-rank system with at least as
    many equations as unknowns.  Returns (x, max_residual).  Exact for
    Fractions and sympy entries; least squares for floats."""
    flat = [v for row in rows for v in row] + list(rhs)
    if any(isinstance(v, float) for v in flat):
        import numpy as np
        A = np.array(rows, dtype=float)
        b = np.array(rhs, dtype=float)
        x = np.linalg.lstsq(A, b, rcond=None)[0]
        return list(x), float(np.max(np.abs(A @ x - b))) if len(b) else 0.0
    import sympy as sp

    def conv(v):
        if isinstance(v, Fraction):
            return sp.Rational(v.numerator, v.denominator)
        return sp.sympify(v)

    symbolic = any(_is_sympy(v) and getattr(v, "free_symbols", None) for v in flat)
    A = sp.Matrix([[conv(v) for v in row] for row in rows])
    b = sp.Matrix([conv(v) for v in rhs])
    n = A.shape[1]
    if symbolic:
        A = A.applyfunc(sp.cancel)
        b = b.applyfunc(sp.cancel)
    # pick n independent equations, solve, then test all of them
    chosen: list = []
    for i in range(A.shape[0]):
        trial = chosen + [i]
        if A.extract(trial, list(range(n))).rank() == len(trial):
            chosen = trial
        if len(chosen) == n:
            break
    if len(chosen) < n:
        raise ValueError("rank-deficient system: add longer words")
    x = A.extract(chosen, list(range(n))).LUsolve(b.extract(chosen, [0]))
    if symbolic:
        x = x.applyfunc(sp.cancel)
    res = A * x - b
    if symbolic:
        res = res.applyfunc(sp.cancel)
    out = []
    for v in x:
        if symbolic:
            out.append(v)
        else:
            out.append(Fraction(int(v.p), int(v.q)))
    worst = max((abs(v) for v in res), default=0) if not symbolic else \
        (0 if all(v == 0 for v in res) else sp.Integer(1))
    return out, worst


# ---------------------------------------------------------------------------
# the engine
# ---------------------------------------------------------------------------
def kappa(n: int) -> dict:
    """The pair coefficients of identity (1): kappa_ab = 4 (-1)^(n-b)
    C(n-1, b-1), slots 1 <= a < b <= n."""
    return {(a, b): 4 * (-1) ** (n - b) * comb(n - 1, b - 1)
            for a in range(1, n + 1) for b in range(a + 1, n + 1)}


def identity_coefficients(n: int) -> dict:
    """Every coefficient of the exact identity (1), keyed by the slot set
    T (a sorted tuple, |T| >= 2)."""
    from itertools import combinations
    return {T: (-2) ** len(T) * (-1) ** (n - max(T)) * comb(n - 1, max(T) - 1)
            for r in range(2, n + 1) for T in combinations(range(1, n + 1), r)}


def a_terms(n: int) -> list:
    """A_n as a list of (coefficient, i, k, j): coefficient x H^i Omega_k
    H^j, the N = infinity form (2)."""
    return [(c, a - 1, b - a - 1, n - b) for (a, b), c in sorted(kappa(n).items())]


class TowerChord:
    """Exact N = infinity words of the Omega_k blocks and the A_n channels
    at one (lambda, q).  lam and q may be floats, Fractions or sympy
    symbols; q = exp(-lam) is NOT imposed, because every identity here is
    a polynomial identity in the two independent variables."""

    def __init__(self, lam, q=None):
        if q is None:
            q = math.exp(-lam)
        self.lam, self.q = lam, q
        self.mu = lam / 2
        self._cache: dict = {}
        self._decomp: dict = {}
        self.norm: dict = {}

    # -- raw chord sums ----------------------------------------------------
    def _memo(self, key, fn):
        if key not in self._cache:
            self._cache[key] = fn()
        return self._cache[key]

    def h_moment(self, k: int):
        """<H^k>, chord units."""
        return self._memo(("h", k), lambda: sum(
            c * self.q ** nx for nx, c in h_table(k)))

    def matter_word(self, A: int, B: int, m):
        """G_m(A, B) = <H^A O_m H^B O_m> for a unit-normalized matter
        chord whose H-crossings weigh q^m."""
        def fn():
            qt = self.q ** m
            return sum(c * self.q ** nx * qt ** nm
                       for (nx, nm), c in matter_table(A, B))
        return self._memo(("g", A, B, m), fn)

    def omega_one(self, k: int, A: int):
        """<H^A Omega_k> under the chord rule (the O(p) conserved piece
        of the self-contracted term dropped)."""
        def fn():
            q, mu, two_lam = self.q, self.mu, 2 * self.lam
            return sum(c * two_lam ** r * q ** nx * (-mu if crossed else mu)
                       for (crossed, nx, r), c in omega_one_table(k, A))
        return self._memo(("o1", k, A), fn)

    def omega_word(self, k1: int, k2: int, A: int, B: int):
        """<H^A Omega_k1 H^B Omega_k2> under the chord rule."""
        def fn():
            q, mu, two_lam = self.q, self.mu, 2 * self.lam
            w = (mu + mu * mu, mu * mu - mu, mu * mu, -mu * mu)
            return sum(c * two_lam ** r * q ** nx * w[kind]
                       for (kind, nx, r), c in omega_word_table(k1, k2, A, B))
        return self._memo(("o2", k1, k2, A, B), fn)

    # -- q-Hermite polynomials (orthogonal for the chord measure) ----------
    def qint(self, n: int):
        return sum(self.q ** j for j in range(n)) if n else 0

    def hermite(self, n: int) -> dict:
        """h_n(x): h_{n+1} = x h_n - [n]_q h_{n-1}, as {power: coeff}."""
        def fn():
            H = [{0: 1}, {1: 1}]
            for m in range(1, n):
                nxt: dict = {}
                for k_, c in H[m].items():
                    nxt[k_ + 1] = nxt.get(k_ + 1, 0) + c
                for k_, c in H[m - 1].items():
                    nxt[k_] = nxt.get(k_, 0) - self.qint(m) * c
                H.append(nxt)
            return H[n]
        return self._memo(("herm", n), fn)

    def hermite_norm(self, n: int):
        """<h_n^2> = [n]_q!."""
        out = 1
        for m in range(1, n + 1):
            out = out * self.qint(m)
        return out

    def projection(self, k: int, j: int):
        """<h_j(H) Omega_k>: the conserved component of Omega_k along h_j."""
        return self._memo(("proj", k, j), lambda: sum(
            c * self.omega_one(k, a) for a, c in self.hermite(j).items()))

    def stationary_word(self, k1: int, k2: int, A: int, B: int):
        """The part of <H^A Omega_k1 H^B Omega_k2> carried by the
        projections of the two blocks onto polynomials in H."""
        def fn():
            tot = 0
            for j1 in range(k1 % 2, k1 + 3, 2):
                c1 = self.projection(k1, j1)
                if c1 == 0:
                    continue
                for j2 in range(k2 % 2, k2 + 3, 2):
                    c2 = self.projection(k2, j2)
                    if c2 == 0:
                        continue
                    hh = sum(ca * cb * self.h_moment(A + B + x + y)
                             for x, ca in self.hermite(j1).items()
                             for y, cb in self.hermite(j2).items())
                    tot = tot + c1 * c2 * hh / (self.hermite_norm(j1)
                                                * self.hermite_norm(j2))
            return tot
        return self._memo(("stat", k1, k2, A, B), fn)

    def dynamical_word(self, k1: int, k2: int, A: int, B: int):
        """<H^A Omega_k1 H^B Omega_k2> with both blocks projected off the
        polynomials in H."""
        return self.omega_word(k1, k2, A, B) - self.stationary_word(k1, k2, A, B)

    # -- the decomposition -------------------------------------------------
    def _sector_word(self, f1: dict, f2: dict, m: int, A: int, B: int):
        """sum_{u,v} f1_u f2_v G_m(A + i_u + i_v, B + j_u + j_v): the
        two-point word of two kernels acting on one unit matter chord."""
        tot = 0
        for (i1, j1), c1 in f1.items():
            for (i2, j2), c2 in f2.items():
                tot = tot + c1 * c2 * self.matter_word(A + i1 + i2,
                                                       B + j1 + j2, m)
        return tot

    def _model_word(self, k1: int, k2: int, A: int, B: int, upto: int):
        """The decomposition's prediction for the dynamical word, sectors
        2..upto."""
        tot = 0
        for m in range(2, upto + 1):
            f1 = self._decomp.get(k1, {}).get(m)
            f2 = self._decomp.get(k2, {}).get(m)
            if f1 and f2:
                tot = tot + self.norm[m] * self._sector_word(f1, f2, m, A, B)
        return tot

    @staticmethod
    def _words(parity: int, max_len: int) -> list:
        return [(a, s - a) for s in range(parity, max_len + 1, 2)
                for a in range(s + 1)]

    def decompose(self, kmax: int = 2, check_len: int = 4) -> dict:
        """Solve Omega_k = sum_m F_{k,m} . O^_m for k <= kmax.

        O^_m is the un-normalized matter chord of dimension m that first
        appears, with coefficient 1, in Omega_(m-2); self.norm[m] is its
        norm.  Sector m of Omega_k is found LINEARLY from the
        cross-correlator <Omega_(m-2) Omega_k>, and the top norm from
        <Omega_k Omega_k> at zero length.  Every other word up to
        `check_len` H's is then a check.

        Returns {"kernels": {k: {m: {(i, j): coeff}}}, "norm": {m: N_m},
        "max_residual": worst deviation over all checked words,
        "n_checked": their number}.  For exact inputs the residual is
        exactly 0.  Levels already decomposed on this instance are
        skipped: a repeat call reports the residual and count of the new
        levels only, while "kernels" covers every level."""
        worst = 0
        n_checked = 0

        def track(res):
            nonlocal worst, n_checked
            n_checked += 1
            if _is_sympy(res):
                import sympy as sp
                res = sp.cancel(res)
                if res != 0:
                    worst = sp.Integer(1)
            elif abs(res) > worst:
                worst = abs(res)

        for k in range(kmax + 1):
            if k in self._decomp:
                continue
            self._decomp[k] = {}
            for m in range(2, k + 2):
                k1 = m - 2
                deg = k + 2 - m
                classes = [(i, d - i) for d in range(deg, -1, -2)
                           for i in range(d, (d - 1) // 2, -1)]
                rows, rhs = [], []
                words = self._words((k1 + k) % 2, max(check_len, deg + 2))
                for (A, B) in words:
                    row = []
                    for (i, j) in classes:
                        g = self.matter_word(A + i, B + j, m)
                        if i != j:
                            g = g + self.matter_word(A + j, B + i, m)
                        row.append(self.norm[m] * g)
                    rows.append(row)
                    rhs.append(self.dynamical_word(k1, k, A, B)
                               - self._model_word(k1, k, A, B, m - 1))
                sol, res = _solve_linear(rows, rhs)
                n_checked += len(rows) - len(classes)
                if res > worst:         # float, Fraction, Rational or the 0/1 sentinel
                    worst = res
                kern: dict = {}
                for (i, j), c in zip(classes, sol):
                    kern[(i, j)] = c
                    kern[(j, i)] = c
                self._decomp[k][m] = kern
            top = k + 2
            n_top = (self.dynamical_word(k, k, 0, 0)
                     - self._model_word(k, k, 0, 0, top - 1))
            if _is_sympy(n_top):
                import sympy as sp
                n_top = sp.factor(sp.cancel(n_top))
            self.norm[top] = n_top
            self._decomp[k][top] = {(0, 0): 1}
            for (A, B) in self._words(0, check_len):
                track(self.dynamical_word(k, k, A, B)
                      - self._model_word(k, k, A, B, top))
            # the cross words with every lower block, all sectors
            for k1 in range(k):
                for (A, B) in self._words((k1 + k) % 2, check_len):
                    track(self.dynamical_word(k1, k, A, B)
                          - self._model_word(k1, k, A, B, top))
        return dict(kernels={k: {m: dict(f) for m, f in self._decomp[k].items()}
                             for k in range(kmax + 1)},
                    norm=dict(self.norm), max_residual=worst,
                    n_checked=n_checked)

    # -- the A_n channels --------------------------------------------------
    def a_kernels(self, n: int) -> dict:
        """F_{n,m} as {m: {(i, j): coeff}}: A_n = sum_m F_{n,m} . O^_m
        modulo polynomials in H.  Needs decompose(n - 2)."""
        self.decompose(n - 2)
        out: dict = {}
        for c, i, k, j in a_terms(n):
            for m, f in self._decomp[k].items():
                dst = out.setdefault(m, {})
                for (i2, j2), c2 in f.items():
                    key = (i + i2, j + j2)
                    dst[key] = dst.get(key, 0) + c * c2
        return out

    def a_word(self, n: int, alpha: int, beta: int, swapped: bool = False):
        """<H^alpha A_n H^beta A_n^dag> (or, swapped, <H^alpha A_n^dag
        H^beta A_n>), dynamical-content rule: every exactly conserved
        piece of A_n is dropped, so only combinations that annihilate
        polynomials in H -- the moments mu_k with k >= 1 -- are physical."""
        tot = 0
        terms = a_terms(n)
        for c1, i1, k1, j1 in terms:
            for c2, i2, k2, j2 in terms:
                if swapped:
                    w = self.omega_word(k1, k2, alpha + j1 + j2, i1 + beta + i2)
                else:
                    w = self.omega_word(k1, k2, alpha + i1 + i2, j1 + beta + j2)
                tot = tot + c1 * c2 * w
        return tot

    def a_moment(self, n: int, k: int):
        """mu_k = (1/2)[<ad^k(A_n) A_n^dag> + <ad^k(A_n^dag) A_n>] of the
        FULL A_n channel at N = infinity, k >= 1 (chord units)."""
        if k < 1:
            raise ValueError("a_moment needs k >= 1; mu_0 carries conserved "
                             "weight -- use sector_moment(n, 0)")
        tot = 0
        for s in range(k + 1):
            cf = (-1) ** s * comb(k, s)
            tot = tot + cf * (self.a_word(n, k - s, s)
                              + self.a_word(n, k - s, s, swapped=True))
        return tot / 2

    def a_cross_moment(self, n1: int, n2: int, k: int):
        """<ad^k(A_n1) A_n2^dag> at N = infinity, k >= 1 (one ordering, not
        symmetrized).  The cross-correlator of two levels: it isolates the
        sectors they share, and with n2 = 2 it is the cheapest object in
        which the Omega_2 remainder term appears at leading order."""
        if k < 1:
            raise ValueError("a_cross_moment needs k >= 1")
        tot = 0
        for s in range(k + 1):
            cf = (-1) ** s * comb(k, s)
            for c1, i1, k1, j1 in a_terms(n1):
                for c2, i2, k2, j2 in a_terms(n2):
                    tot = tot + cf * c1 * c2 * self.omega_word(
                        k1, k2, k - s + i1 + i2, j1 + s + j2)
        return tot

    def sector_moment(self, n: int, k: int, m: int = None):
        """mu_k of A_n from its matter-chord decomposition: the sum over
        sectors m (or one sector) of norm_m <|F_{n,m}|^2 (E - E')^k>.  For
        k >= 1 the sum over sectors equals a_moment(n, k); k = 0 gives the
        DYNAMICAL weight, with every conserved component removed."""
        kern = self.a_kernels(n)
        tot = 0
        for mm, f in kern.items():
            if m is not None and mm != m:
                continue
            part = 0
            for s in range(k + 1):
                cf = (-1) ** s * comb(k, s)
                w1 = 0
                w2 = 0
                for (i1, j1), c1 in f.items():
                    for (i2, j2), c2 in f.items():
                        w1 = w1 + c1 * c2 * self.matter_word(
                            k - s + i1 + i2, s + j1 + j2, mm)
                        w2 = w2 + c1 * c2 * self.matter_word(
                            k - s + j1 + j2, s + i1 + i2, mm)
                part = part + cf * (w1 + w2)
            tot = tot + self.norm[mm] * part / 2
        return tot


# ---------------------------------------------------------------------------
# closed forms
# ---------------------------------------------------------------------------
def _five(q):
    return 1 + q + q ** 2 + q ** 3 + q ** 4


def norm2(lam, q):
    """<O^_2 O^_2>: the dynamical weight of Omega_0 = A_2 / 4."""
    return lam * (1 - q * q + 2 * lam * q) / (2 * (1 + q))


def norm3(lam, q):
    return (lam * (1 + q) * (1 - q ** 6 + 2 * lam * q * q * (2 - q + 2 * q * q))
            / (2 * (1 + q * q) * (1 + q + q * q)))


def norm4(lam, q):
    return (lam * (1 - q) ** 2 * (1 + q) ** 2
            * ((1 - q ** 6) * (1 + q + 2 * q ** 2 + q ** 3 + q ** 4)
               + 2 * lam * q ** 3 * (2 + q + 2 * q ** 2 + q ** 3 + 2 * q ** 4))
            / (2 * (1 + q * q) * (1 - q + q * q) * _five(q)))


def omega_kernels_closed(lam, q) -> dict:
    """The kernels of Omega_0, Omega_1, Omega_2 in closed form, same
    layout as TowerChord.decompose()["kernels"]."""
    g2 = q / (1 + q * q)
    g3 = q / (1 - q + q * q)
    a = q * q * (1 + q + q * q) / ((1 + q * q) * _five(q))
    b = q * q * (1 + q) ** 2 / ((1 + q * q) * _five(q))
    c = (1 + q - q ** 2 - q ** 3 + 2 * lam * q * q * _five(q)) / _five(q)
    return {0: {2: {(0, 0): 1}},
            1: {2: {(1, 0): g2, (0, 1): g2}, 3: {(0, 0): 1}},
            2: {2: {(2, 0): a, (0, 2): a, (1, 1): b, (0, 0): c},
                3: {(1, 0): g3, (0, 1): g3}, 4: {(0, 0): 1}}}


def a2_split(lam, q=None) -> dict:
    """The proven split of the full A_2 channel at N = infinity."""
    if q is None:
        q = math.exp(-lam)
    stat = 4 * lam * lam * (1 - q) ** 2 / (1 + q)
    dyn = 16 * norm2(lam, q)
    mu2 = 16 * lam * (1 - q) * (1 - q * q + 2 * lam * q)
    return dict(lam=lam, q=q, h2_coefficient=2 * lam * (1 - q) / (1 + q),
                mu0_stationary=stat, mu0_dynamic=dyn, mu0=stat + dyn,
                stationary_fraction=stat / (stat + dyn), mu2=mu2,
                w2_raw=mu2 / (stat + dyn), w2_dynamic=2 * (1 - q * q))


def a3_mu2(lam, q):
    """mu_2 of the full A_3 channel at N = infinity."""
    return 8 * lam * (1 - q) * (
        25 + 5 * q - 38 * q ** 2 - 16 * q ** 3 + 13 * q ** 4 + 11 * q ** 5
        + lam * (46 * q + 14 * q ** 2 - 22 * q ** 3 - 26 * q ** 4))


def a3_mu4(lam, q):
    """mu_4 of the full A_3 channel at N = infinity."""
    return 8 * lam * (1 - q) ** 2 * (
        169 + 272 * q - 41 * q ** 2 - 384 * q ** 3 - 309 * q ** 4 - 6 * q ** 5
        + 143 * q ** 6 + 118 * q ** 7 + 38 * q ** 8
        + lam * (318 * q + 540 * q ** 2 + 276 * q ** 3 - 208 * q ** 4
                 - 350 * q ** 5 - 248 * q ** 6 - 88 * q ** 7))


def a4_mu2(lam, q):
    """mu_2 of the full A_4 channel at N = infinity."""
    return 16 * lam * (1 - q) * (
        154 + 78 * q - 276 * q ** 2 - 253 * q ** 3 + 76 * q ** 4 + 223 * q ** 5
        + 92 * q ** 6 - 31 * q ** 7 - 46 * q ** 8 - 17 * q ** 9
        + lam * (264 * q + 228 * q ** 2 - 184 * q ** 3 - 426 * q ** 4
                 - 128 * q ** 5 + 80 * q ** 6 + 112 * q ** 7 + 54 * q ** 8)
        + lam ** 2 * (80 * q ** 3 + 12 * q ** 4 - 48 * q ** 5 - 44 * q ** 6)
        + lam ** 3 * 8 * q ** 5)


# ---------------------------------------------------------------------------
# the exact spectral density of the A_n channel
# ---------------------------------------------------------------------------
def spectral_density(n: int, omegas, lam: float, n_theta: int = 1200) -> dict:
    """The exact N = infinity, T_B = infinity spectral density of the
    DYNAMICAL part of the full A_n channel, per matter-chord sector:

        S_m(omega) = norm_m int dE dE' rho(E) rho(E') V_m(E, E')
                     |F_{n,m}(E, E')|^2 delta(omega - (E - E')),

    symmetrized in omega (the two operator orderings).  rho and V_m are
    the DSSYK density of states and the Delta = m matter vertex (review
    2407.09396 Eqs. 2.21, 2.30, as coded in chord.py), and F_{n,m} is the
    kernel polynomial of a_kernels(n).  Chord units, <H^2> = 1.

    Returns {"omega": [...], "sectors": {m: [...]}, "total": [...]}.  The
    conserved part of A_n is a delta function at omega = 0 and is not
    included.  int S = sector_moment(n, 0) and int omega^2 S =
    a_moment(n, 2) are tested."""
    import numpy as np
    from chord import _vertex2_pairs, energy_theta, rho_theta
    q = math.exp(-lam)
    T = TowerChord(lam, q)
    kern = T.a_kernels(n)
    band = 2.0 / math.sqrt(1.0 - q)
    x, wq = np.polynomial.legendre.leggauss(n_theta)
    th1 = 0.5 * np.pi * (x + 1.0)
    w1 = 0.5 * np.pi * wq
    rho1 = rho_theta(th1, q)
    E1 = energy_theta(th1, q)
    omegas = np.atleast_1d(np.asarray(omegas, dtype=float))
    out = {m: np.zeros(len(omegas)) for m in kern}

    def kernel_sq(f, Ea, Eb):
        F = sum(c * Ea ** i * Eb ** j for (i, j), c in f.items())
        return F * F

    for i, om in enumerate(omegas):
        for sgn in (1.0, -1.0):
            E2 = E1 - sgn * om
            ok = np.abs(E2) < band * (1.0 - 1e-12)
            th2 = np.arccos(np.clip(E2[ok] / band, -1.0, 1.0))
            jac = band * np.sin(th2)
            base = w1[ok] * rho1[ok] * rho_theta(th2, q) / jac
            for m, f in kern.items():
                V = _vertex2_pairs(th1[ok], th2, q, float(m))
                out[m][i] += 0.5 * float(T.norm[m]) * float(np.sum(
                    base * V * kernel_sq(f, E1[ok], E2[ok])))
    total = sum(out.values())
    return dict(omega=omegas.tolist(),
                sectors={m: v.tolist() for m, v in out.items()},
                total=total.tolist())


def spectra_block(lam: float, levels=(2, 3, 4), n_omega: int = 201,
                  n_theta: int = 800) -> dict:
    """The exact dynamical spectral densities of A_2..A_4 at one lambda,
    on [0, 2 band], with their moment closures and the position of the
    maximum (omega = 0 for A_2: monotone)."""
    import numpy as np
    from chord import _round6
    q = math.exp(-lam)
    band = 2.0 / math.sqrt(1.0 - q)
    om = np.linspace(0.0, 2.0 * band, n_omega)
    T = TowerChord(lam, q)
    out = dict(lam=lam, band=band, omega=_round6(om),
               levels={})
    for n in levels:
        S = spectral_density(n, om, lam, n_theta=n_theta)
        tot = np.array(S["total"])
        out["levels"][str(n)] = dict(
            sectors={str(m): _round6(v) for m, v in S["sectors"].items()},
            total=_round6(tot),
            mu0_integral=float(2.0 * np.trapezoid(tot, om)),
            mu0_exact=float(T.sector_moment(n, 0)),
            mu2_integral=float(2.0 * np.trapezoid(tot * om ** 2, om)),
            mu2_exact=float(T.a_moment(n, 2)),
            omega_at_max=float(om[int(np.argmax(tot))]),
            min_value=float(tot.min()))
    return out


# ---------------------------------------------------------------------------
# finite-N gates
# ---------------------------------------------------------------------------
def _syk_path():
    path = str(HERE.parent / "syk-self-averaging")
    if path not in sys.path:
        sys.path.insert(0, path)


def finite_n_full(N: int, p: int, n: int, k: int) -> Fraction:
    """Exact finite-N mu_k of the FULL A_n channel, N mu^diag + N(N-1)
    mu^pair in chord units (exact_pair_moments engine, Fractions)."""
    _syk_path()
    import exact_pair_moments as epm
    pr = epm.channel_pair_moments(N, p, n, marks=('i', 'i', 'j', 'j'),
                                  ks=(k,), exact=True)
    dg = epm.channel_pair_moments(N, p, n, marks=('i', 'i', 'i', 'i'),
                                  ks=(k,), exact=True)
    epm._DP_CACHE.clear()
    epm._EXP_CACHE.clear()
    return N * dg[f"mu{k}"] + N * (N - 1) * pr[f"mu{k}"]


def finite_n_cross(N: int, p: int, n1: int, n2: int, k: int) -> Fraction:
    """Exact finite-N <ad^k(A_n1) A_n2^dag>, flavor-summed, chord units:
    the binomial ad-expansion of chord_moments._mu_ordered with two
    different levels, on the exact finite-N word tables."""
    _syk_path()
    import exact_pair_moments as epm

    def ordered(tab):
        tot = 0
        for s in range(k + 1):
            for r in range(n1 + 1):
                for rp in range(n2 + 1):
                    coeff = ((-1) ** (n2 + r + rp + s) * comb(k, s)
                             * comb(n1, r) * comb(n2, rp))
                    tot += coeff * tab.t4(k - s, n1 - r, r + s + n2 - rp, rp)
        return tot

    pr = ordered(epm.PairTab(N, p, ('i', 'i', 'j', 'j'), exact=True))
    dg = ordered(epm.PairTab(N, p, ('i', 'i', 'i', 'i'), exact=True))
    epm._DP_CACHE.clear()
    epm._EXP_CACHE.clear()
    return N * dg + N * (N - 1) * pr


def _gate_point(task):
    n, k, N, p = task
    t0 = time.time()
    if isinstance(n, tuple):
        mu = finite_n_cross(N, p, n[0], n[1], k)
    else:
        mu = finite_n_full(N, p, n, k)
    return dict(n=list(n) if isinstance(n, tuple) else n, k=k, N=N, p=p,
                mu=float(mu), mu_exact=str(mu),
                runtime_s=round(time.time() - t0, 1))


def h2_overlap(N: int, p: int):
    """Exact finite-N pieces of the H^2 projection of A_2, chord units:
    (<A_2 (H^2 - 1)>, <(H^2 - 1)^2>).  The intersection of two uniform
    p-subsets is hypergeometric; the two Wick pairings give
    4 sum_j P(j) j (1 + (-1)^j), and <H^4> - 1 = 1 + K/C (exact_wick)."""
    C = comb(N, p)
    P = [Fraction(comb(p, j) * comb(N - p, p - j), C) for j in range(p + 1)]
    ov = 4 * sum(P[j] * j * (1 + (-1) ** j) for j in range(p + 1))
    qN = sum(P[j] * (-1) ** j for j in range(p + 1))
    return ov, 1 + qN


def finite_n_a2_split(N: int, p: int, mu0, mu2) -> dict:
    """The exact finite-N full A_2 channel with its H^2 component removed:
    mu0 and mu2 are the exact full-channel moments (Fractions)."""
    ov, nrm = h2_overlap(N, p)
    stat = ov * ov / nrm
    dyn = mu0 - stat
    return dict(N=N, p=p, mu0_stationary=stat, mu0_dynamic=dyn,
                w2_dynamic=mu2 / dyn, w2_raw=mu2 / mu0)


def ed_commutant_split(N: int, p: int, n_inst: int, seed0: int = 1000) -> dict:
    """Dense ED: the full A_2 channel split per instance by removing the
    whole commutant of H (projected_channel.spectral_split), disorder-
    averaged, chord units.  This shares nothing with the chord rule."""
    import numpy as np
    _syk_path()
    from conserved_points import dense_majoranas, flavor_parts
    from pauli_strings import majorana_hamiltonian_fast
    from projected_channel import spectral_split
    from syk import coupling_variance
    psis = dense_majoranas(N)
    s2 = comb(N, p) * coupling_variance(N, p)            # <H^2>/dim
    rows = []
    for i in range(n_inst):
        H = majorana_hamiltonian_fast(N, p, np.random.default_rng(seed0 + i))
        A2 = 4 * sum(h @ h for h in flavor_parts(H, psis))
        r = spectral_split(H, A2)
        rows.append((r["mu0_traceless"] / s2 ** 2, r["mu0_stationary"] / s2 ** 2,
                     r["mu0_dynamic"] / s2 ** 2, r["mu2_dynamic"] / s2 ** 3))
    arr = np.array(rows)
    m = arr.mean(axis=0)
    e = arr.std(axis=0, ddof=1) / math.sqrt(n_inst)
    # ratio of ensemble means, with its error from the per-instance ratios
    ratios = arr[:, 3] / arr[:, 2]
    ov, nrm = h2_overlap(N, p)
    lam = 2.0 * p * p / N
    return dict(N=N, p=p, lam_chord=lam, n_inst=n_inst, seed0=seed0,
                mu0=float(m[0]), mu0_err=float(e[0]),
                mu0_stationary=float(m[1]), mu0_stationary_err=float(e[1]),
                mu0_dynamic=float(m[2]), mu0_dynamic_err=float(e[2]),
                mu2=float(m[3]), mu2_err=float(e[3]),
                w2_dynamic=float(m[3] / m[2]),
                w2_dynamic_err=float(ratios.std(ddof=1) / math.sqrt(n_inst)),
                w2_raw=float(m[3] / m[0]),
                h2_projection_exact=float(ov * ov / nrm),
                w2_dynamic_ninf=a2_split(lam)["w2_dynamic"],
                w2_raw_ninf=a2_split(lam)["w2_raw"])


def fixed_p_dynamic(p: int, Ns) -> list:
    """The exact finite-N full A_2 channel along FIXED p with its H^2
    component removed.  In the identity-subtracted normalization the
    width has a maximum (N = 24 at p = 4, full_channel_laws.py); that
    maximum is the conserved weight draining away, not dynamics.  The
    dynamical width falls monotonically from ~2 beyond the smallest
    sizes (N >= 12 at p = 4, N >= 16 at p = 6) and its physical value
    rises monotonically to 8p(p-2) J^2."""
    from math import factorial
    rows = []
    for N in Ns:
        mu0 = finite_n_full(N, p, 2, 0)
        mu2 = finite_n_full(N, p, 2, 2)
        sp_ = finite_n_a2_split(N, p, mu0, mu2)
        lam = 2.0 * p * p / N
        s2 = Fraction(comb(N, p) * factorial(p), N ** (p - 1))   # <H^2>/dim, J = 1
        rows.append(dict(
            p=p, N=N, lam_chord=lam, w2_raw=float(sp_["w2_raw"]),
            w2_dynamic=float(sp_["w2_dynamic"]),
            stationary_fraction=float(sp_["mu0_stationary"] / mu0),
            w2_dynamic_ninf=2.0 * (1.0 - math.exp(-2.0 * lam)),
            omega2_dynamic_physical=float(sp_["w2_dynamic"] * s2),
            omega2_limit=8 * p * (p - 2)))
    return rows


FIXED_P_GRID = {4: (10, 12, 14, 16, 20, 24, 32, 48, 64, 128, 256, 512, 1024),
                6: (14, 16, 20, 24, 36, 48, 72, 144, 288, 576)}

# the exact finite-N series behind the artifact: (n, k, lambda): p values
# n is a level, or a pair (n1, n2) for the cross moment <ad^k(A_n1) A_n2^dag>
GATE_GRID = {
    (2, 4, 2.0): (4, 6, 8, 10, 12, 14, 16, 18, 20),
    (3, 2, 2.0): (4, 6, 8, 10, 12, 14, 16, 18, 20),
    (3, 4, 2.0): (4, 6, 8, 10, 12),
    (4, 2, 2.0): (4, 6, 8, 10, 12),
    ((4, 2), 2, 2.0): (4, 6, 8, 10, 12, 14, 16, 18, 20),
    ((4, 2), 2, 1.0): (4, 6, 8, 10, 12, 14, 16, 18, 20),
}

# dense-ED commutant split: (N, p, number of instances)
ED_GRID = ((12, 4, 40), (14, 4, 30), (16, 4, 20), (18, 4, 10),
           (16, 6, 12), (18, 6, 8), (20, 6, 4))


def gate_key(n, k: int, lam: float) -> str:
    tag = f"cross{n[0]}{n[1]}" if isinstance(n, (tuple, list)) else f"n{n}"
    return f"{tag}_k{k}_lam{lam:g}"


def closed_form_target(n, k: int, lam: float) -> float:
    q = math.exp(-lam)
    if isinstance(n, (tuple, list)):
        return float(TowerChord(lam, q).a_cross_moment(n[0], n[1], k))
    if (n, k) == (2, 4):
        from full_channel_chord import mu4_closed
        return mu4_closed(lam, q)
    return {(3, 2): a3_mu2, (3, 4): a3_mu4, (4, 2): a4_mu2}[(n, k)](lam, q)


def remainder_control(lams=(1.0, 2.0)) -> dict:
    """<ad^2(A_4) A_2^dag> and mu_2(A_4) with the O(1) remainder of the
    self-contracted Omega_2 term switched off, next to the true values.

    The switch patches _omega_terms, which tables preloaded by
    build_word_tables would bypass, so it refuses to run with any."""
    global _omega_terms
    if _PRELOADED:
        raise RuntimeError("remainder_control needs freshly enumerated "
                           "tables; _PRELOADED is not empty")
    saved = _omega_terms
    out = {}
    try:
        _omega_terms = lambda sc: saved(sc) if sc[0] != sc[-1] else ()
        omega_word_table.cache_clear()
        omega_one_table.cache_clear()
        off = {lam: (float(TowerChord(lam).a_cross_moment(4, 2, 2)),
                     float(TowerChord(lam).a_moment(4, 2))) for lam in lams}
    finally:
        _omega_terms = saved
        omega_word_table.cache_clear()
        omega_one_table.cache_clear()
    for lam in lams:
        T = TowerChord(lam)
        out[f"lam{lam:g}"] = dict(
            cross42_with=float(T.a_cross_moment(4, 2, 2)),
            cross42_without=off[lam][0],
            mu2_a4_with=float(T.a_moment(4, 2)), mu2_a4_without=off[lam][1])
    return out


def series_limit(rows: list, p_min: int = 6) -> dict:
    """Polynomial-in-1/p least-squares limits of an exact series, every
    order from 2 up to (points - 2), capped at 5.  The spread over the
    orders is a sensitivity estimate, not a confidence interval.  With
    fewer than two orders there is nothing to compare, and the spread is
    None (not 0)."""
    import numpy as np
    pts = sorted((r["p"], r["mu"]) for r in rows if r["p"] >= p_min)
    x = np.array([1.0 / p for p, _ in pts])
    y = np.array([v for _, v in pts])
    fits = []
    for order in range(2, min(5, len(pts) - 2) + 1):
        A = np.vander(x, order + 1, increasing=True)
        c = np.linalg.lstsq(A, y, rcond=None)[0]
        fits.append(dict(order=order, limit=float(c[0])))
    lims = [f["limit"] for f in fits]
    return dict(p_used=[p for p, _ in pts], fits=fits,
                central=lims[-1] if lims else None,
                spread=(max(lims) - min(lims)) if len(lims) > 1 else None)


def _write_json(path: Path, payload: dict) -> None:
    import os
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)


def norm_scaling(lams=(Fraction(1, 5), Fraction(1, 10), Fraction(1, 20),
                       Fraction(1, 40), Fraction(1, 80)), kmax: int = 4) -> list:
    """The norms of the matter chords toward small lambda, in exact
    arithmetic at a rational q within 1e-9 of e^-lambda.  norm_m scales
    as lambda^(2 floor(m/2)): the ratios recorded are norm_m over that
    power."""
    rows = []
    for lam in lams:
        q = Fraction(math.exp(-float(lam))).limit_denominator(10 ** 9)
        T = TowerChord(lam, q)
        T.decompose(kmax, check_len=0)
        rows.append(dict(lam=float(lam), norm_over_power={
            str(m): float(v) / float(lam) ** (2 * (m // 2))
            for m, v in T.norm.items()}))
    return rows


def share_scaling(lams=(Fraction(1, 5), Fraction(1, 10), Fraction(1, 20),
                        Fraction(1, 40), Fraction(1, 80), Fraction(1, 160),
                        Fraction(1, 320), Fraction(1, 1280)),
                  kmax: int = 4) -> list:
    """Which chord carries the dynamical weight of A_3..A_(kmax+2) toward
    small lambda, in exact arithmetic at a rational q next to e^-lambda.
    The even chords fade linearly in lambda; the odd ones do not."""
    rows = []
    for lam in lams:
        q = Fraction(math.exp(-float(lam))).limit_denominator(10 ** 13)
        T = TowerChord(lam, q)
        T.decompose(kmax, check_len=0)
        row = dict(lam=float(lam), shares={})
        for n in range(3, kmax + 3):
            mu0 = {m: T.sector_moment(n, 0, m) for m in T.a_kernels(n)}
            tot = sum(mu0.values())
            row["shares"][str(n)] = {str(m): float(v / tot)
                                     for m, v in sorted(mu0.items())}
        rows.append(row)
    return rows


EXACT_POINTS = [(Fraction(2), Fraction(1, 7)), (Fraction(1, 3), Fraction(3, 5)),
                (Fraction(5, 2), Fraction(2, 9))]


def deep_level(kmax: int, workers: int = None) -> dict:
    """Levels n = kmax + 2 beyond singlet_tower.json's n <= 6: the same exact
    checks, the lambda = 2 shares, and the small-lambda norm and share
    scalings, with the word tables built in parallel first."""
    t0 = time.time()
    build_word_tables(needed_word_keys(kmax, 2), workers)
    tables_s = round(time.time() - t0, 1)
    exact = []
    for lam, q in EXACT_POINTS:
        T = TowerChord(lam, q)
        d = T.decompose(kmax, check_len=2)
        sub = {}
        for k in range(1, kmax + 1):
            qk = sum(q ** j for j in range(1, k + 1)) / (1 + q ** (k + 1))
            sub[str(k)] = d["kernels"][k][k + 1] == {(1, 0): qk, (0, 1): qk}
        exact.append(dict(
            lam=str(lam), q=str(q), max_residual=str(d["max_residual"]),
            n_checked=d["n_checked"],
            norms_positive=all(v > 0 for v in d["norm"].values()),
            subleading_closed_form=sub,
            kernels_symmetric=all(
                f[(j, i)] == c for k in range(kmax + 1)
                for f in d["kernels"][k].values() for (i, j), c in f.items())))
        print(f"  exact ({lam}, {q}): residual {d['max_residual']} on "
              f"{d['n_checked']} words", flush=True)
    return dict(kmax=kmax, level=kmax + 2, exact_checks=exact,
                scan_lam2=scan_row(2.0, kmax),
                norm_scaling=norm_scaling(kmax=kmax),
                share_scaling=share_scaling(kmax=kmax),
                tables_s=tables_s, runtime_s=round(time.time() - t0, 1))


def run_deep(kmaxes, out_path=None, workers: int = None) -> dict:
    """singlet_tower_deep.json: the levels past n = 6 (kmax = 5 is n = 7,
    kmax = 6 is n = 8).  Rows already in the file are kept."""
    out = Path(out_path) if out_path else HERE / "singlet_tower_deep.json"
    payload = {"levels": {}}
    if out.exists():
        with open(out) as fh:
            payload = json.load(fh)
    for kmax in kmaxes:
        print(f"level n = {kmax + 2}", flush=True)
        payload["levels"][str(kmax + 2)] = deep_level(kmax, workers)
        payload["provenance"] = (
            "singlet_tower.py --deep -- the matter-chord decomposition past "
            "n = 6: exact rational checks at three generic (lambda, q), the "
            "subleading coefficient q[k]_q/(1+q^(k+1)), sector shares at "
            "lambda = 2 and the small-lambda norm and share scalings")
        _write_json(out, payload)
    print(f"wrote {out}")
    return payload


def _exact_report(kmax: int, points) -> list:
    out = []
    for lam, q in points:
        T = TowerChord(lam, q)
        d = T.decompose(kmax, check_len=4 if kmax <= 2 else (3 if kmax == 3 else 2))
        closed = omega_kernels_closed(lam, q)
        out.append(dict(
            lam=str(lam), q=str(q), kmax=kmax,
            max_residual=str(d["max_residual"]), n_checked=d["n_checked"],
            closed_kernels_match=all(d["kernels"][k] == closed[k]
                                     for k in range(min(kmax, 2) + 1)),
            closed_norms_match=(d["norm"][2] == norm2(lam, q)
                                and d["norm"][3] == norm3(lam, q)
                                and d["norm"][4] == norm4(lam, q)),
            norm={str(m): str(v) for m, v in d["norm"].items()},
            moments_match=all(
                T.a_moment(n, k) == f(lam, q) == T.sector_moment(n, k)
                for n, k, f in ((3, 2, a3_mu2), (3, 4, a3_mu4), (4, 2, a4_mu2)))))
    return out


def scan_row(lam: float, kmax: int = 3) -> dict:
    """The tower at one lambda_chord: norms, the A_2 split, and the
    sector content of A_2..A_(kmax+2) (dynamical weight, width,
    kurtosis)."""
    q = math.exp(-lam)
    T = TowerChord(lam, q)
    T.decompose(kmax, check_len=1)
    row = dict(lam=lam, q=q, kmax=kmax,
               norm={str(m): float(v) for m, v in T.norm.items()},
               a2=a2_split(lam, q), channels={})
    for n in range(2, kmax + 3):
        kern = T.a_kernels(n)
        mu0 = {m: float(T.sector_moment(n, 0, m)) for m in kern}
        mu2 = float(T.sector_moment(n, 2))
        mu4 = float(T.sector_moment(n, 4))
        tot = sum(mu0.values())
        row["channels"][str(n)] = dict(
            mu0_dynamic=tot, mu2=mu2, mu4=mu4,
            sector_share={str(m): v / tot for m, v in mu0.items()},
            w2_dynamic=mu2 / tot, kurt_dynamic=mu4 * tot / mu2 ** 2)
    return row


def run(out_path=None, workers: int = None, resume: bool = False,
        grid: dict = None, ed_grid=None) -> dict:
    import multiprocessing as mp
    t0 = time.time()
    out = Path(out_path) if out_path else HERE / "singlet_tower.json"
    grid = GATE_GRID if grid is None else grid
    ed_grid = ED_GRID if ed_grid is None else ed_grid
    payload: dict = {}
    if resume and out.exists():
        with open(out) as fh:
            payload = json.load(fh)
    prior_runtime = float(payload.get("runtime_s", 0.0))

    print("exact decomposition checks", flush=True)
    pts = EXACT_POINTS
    payload["exact_checks"] = _exact_report(4, pts)
    # one point pushed to longer words for the blocks with closed forms
    T = TowerChord(*pts[0])
    deep = T.decompose(2, check_len=6)
    payload["exact_check_deep"] = dict(
        lam=str(pts[0][0]), q=str(pts[0][1]), kmax=2, check_len=6,
        max_residual=str(deep["max_residual"]), n_checked=deep["n_checked"])
    payload["scan"] = [scan_row(lam, 4) for lam in
                       (4.0, 3.0, 2.0, 1.5, 1.0, 0.5, 0.25, 0.1, 0.05)]
    payload["norm_scaling"] = norm_scaling()
    payload["share_scaling"] = share_scaling()

    # n = 2 against the committed exact finite-N full-channel rows
    with open(HERE / "full_channel_chord.json") as fh:
        fcc = json.load(fh)
    a2_gate = {}
    for key, g in fcc["gates"].items():
        lam = g["lam"]
        rows = []
        for r in g["rows"]:
            sp_ = finite_n_a2_split(r["N"], r["p"], Fraction(r["mu0_exact"]),
                                    Fraction(r["mu2_exact"]))
            rows.append(dict(p=r["p"], N=r["N"],
                             mu0_stationary=float(sp_["mu0_stationary"]),
                             mu0_dynamic=float(sp_["mu0_dynamic"]),
                             w2_dynamic=float(sp_["w2_dynamic"]),
                             w2_raw=float(sp_["w2_raw"])))
        lim = {name: series_limit([dict(p=r["p"], mu=r[name]) for r in rows],
                                  p_min=8)
               for name in ("mu0_stationary", "mu0_dynamic", "w2_dynamic")}
        cf = a2_split(lam)
        a2_gate[key] = dict(lam=lam, rows=rows, limits=lim,
                            closed_form={k: cf[k] for k in
                                         ("mu0_stationary", "mu0_dynamic",
                                          "w2_dynamic")})
    payload["a2_gate"] = a2_gate

    # exact finite-N series for the higher moments and levels
    gates = payload.setdefault("gates", {})

    def norm_n(n):
        return tuple(n) if isinstance(n, (tuple, list)) else n

    have = {(norm_n(r["n"]), r["k"], r["N"], r["p"])
            for g in gates.values() for r in g["rows"]}
    tasks = []
    for (n, k, lam), ps in grid.items():
        for p in ps:
            N = int(round(2 * p * p / lam))
            if (n, k, N, p) not in have:
                tasks.append((n, k, N, p))
    workers = workers or max(1, mp.cpu_count() - 1)
    print(f"{len(tasks)} exact finite-N points on {workers} workers", flush=True)

    def checkpoint():
        payload["runtime_s"] = round(prior_runtime + time.time() - t0, 1)
        _write_json(out, payload)

    if tasks:
        def words(t):          # H-nodes per word: the cost driver
            n = t[0]
            return t[1] + (sum(n) if isinstance(n, tuple) else 2 * n)

        tasks.sort(key=lambda t: -(t[3] ** (words(t) // 2) * t[2]))
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=min(workers, len(tasks))) as pool:
            for r in pool.imap_unordered(_gate_point, tasks):
                lam = 2.0 * r["p"] ** 2 / r["N"]
                key = gate_key(r["n"], r["k"], lam)
                g = gates.setdefault(key, dict(n=r["n"], k=r["k"], lam=lam,
                                               rows=[]))
                g["rows"].append(r)
                g["rows"].sort(key=lambda x: x["p"])
                print(f"  {key} p={r['p']:2d} N={r['N']:4d}: {r['mu']:.4f} "
                      f"({r['runtime_s']} s)", flush=True)
                checkpoint()
    for g in gates.values():
        g["closed_form"] = closed_form_target(g["n"], g["k"], g["lam"])
        g["limit"] = series_limit(g["rows"])
        if g["limit"]["central"] is not None:
            g["rel_dev"] = abs(g["limit"]["central"] / g["closed_form"] - 1.0)
    payload["gates_runtime_s"] = round(sum(
        r["runtime_s"] or 0.0 for g in gates.values() for r in g["rows"]), 1)
    # the Omega_2 remainder term, switched off: what the cross moment
    # would be without it (the gate that tells the two apart)
    payload["remainder_control"] = remainder_control()

    print("spectral densities", flush=True)
    payload["spectra"] = spectra_block(1.0)

    print("fixed-p dynamical width", flush=True)
    payload["fixed_p_dynamic"] = {str(p): fixed_p_dynamic(p, Ns)
                                  for p, Ns in FIXED_P_GRID.items()}

    print("ED commutant split", flush=True)
    payload["ed_commutant"] = [ed_commutant_split(N, p, n_inst)
                               for (N, p, n_inst) in ed_grid]
    payload["closed_forms"] = dict(
        a2_stationary="4 lam^2 (1-q)^2 / (1+q)",
        a2_dynamic="8 lam (1 - q^2 + 2 lam q) / (1+q)",
        a2_w2_dynamic="2 (1 - q^2)",
        norm2="lam (1 - q^2 + 2 lam q) / (2 (1+q))",
        norm3="lam (1+q) [1 - q^6 + 2 lam q^2 (2 - q + 2 q^2)] / "
              "(2 (1+q^2) (1+q+q^2))",
        norm4="lam (1-q)^2 (1+q)^2 [(1 - q^6)(1+q+2q^2+q^3+q^4) + 2 lam q^3 "
              "(2+q+2q^2+q^3+2q^4)] / (2 (1+q^2) (1-q+q^2) [5]_q)",
        omega1="[q/(1+q^2)] {H, O2} + O3",
        omega2="a (H^2 O2 + O2 H^2) + b H O2 H + c O2 + [q/(1-q+q^2)] {H, O3} "
               "+ O4;  a = q^2 [3]_q / ((1+q^2) [5]_q), b = q^2 (1+q)^2 / "
               "((1+q^2) [5]_q), c = (1 + q - q^2 - q^3 + 2 lam q^2 [5]_q)/[5]_q")
    provenance = (
        "singlet_tower.py -- the full A_n singlet channels at N = infinity "
        "as matter chords of dimension 2..n: exact decomposition checks in "
        "rational arithmetic, the lambda scan, the A_2 split against the "
        "exact finite-N full channel (rows of full_channel_chord.json, H^2 "
        "component removed) and against dense-ED commutant projection, and "
        "exact finite-N series for mu_4(A_2), mu_2 and mu_4 of A_3, "
        "mu_2(A_4) on double-scaling lines with their 1/p limits")
    # provenance first, as in the other artifacts (in place: checkpoint()
    # closes over this dict)
    rest = {k: v for k, v in payload.items() if k != "provenance"}
    payload.clear()
    payload["provenance"] = provenance
    payload.update(rest)
    checkpoint()
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--deep", type=int, nargs="+", metavar="KMAX",
                    help="levels past n = 6: KMAX = n - 2 (5 and 6)")
    args = ap.parse_args()
    if args.deep:
        run_deep(args.deep, workers=args.workers)
    elif args.run or args.resume:
        run(workers=args.workers, resume=args.resume)
    else:
        for lam in (2.0, 1.0, 0.5, 0.1):
            s = a2_split(lam)
            print(f"lambda_chord={lam:<4}: A_2 stationary fraction "
                  f"{s['stationary_fraction']:.4f}  w2 raw {s['w2_raw']:.6f}  "
                  f"w2 dynamic {s['w2_dynamic']:.6f}")
        row = scan_row(2.0)
        for n, ch in row["channels"].items():
            share = "  ".join(f"O{m}: {v:.3f}" for m, v in ch["sector_share"].items())
            print(f"A_{n} at lambda_chord=2: dynamical width {ch['w2_dynamic']:.4f}, "
                  f"kurtosis {ch['kurt_dynamic']:.3f}; weight by sector  {share}")
