#!/usr/bin/env python3
r"""
Exact finite-N pair moments of the A_n channels: the approach to the zero.
==========================================================================

The conserved-channel theorem (duality/chord_charges.py) proves that the
i != j pair moments of A_n = sum_i psi_i (ad_H)^n psi_i vanish at strict
N = infinity, every lambda, every matter weight.  This module computes
the SAME moments EXACTLY at finite N -- disorder-averaged by
Wick/Isserlis, subset sums collapsed to intersection-parity profiles.
No ED, no sampling noise: the measured 1/N approach to the proven zero
(the first bite of this task's tower extraction).

Method
------
E<Tr[H^a psi_i H^b psi_i H^c psi_j H^d psi_j]>/dim reduces by Isserlis to
a sum over perfect matchings of the K = a+b+c+d H-slots.  Each matched
pair carries one independent uniform p-subset.  The trace of any
resulting Majorana word obeys the chord-diagram sign law (collapse
argument; Gamma_c^2 = +1 in the i^{p/2} convention, so the all-disjoint
base trace is +1):

    trace = prod_{crossing H-chord pairs (al,be)} (-1)^{|c_al & c_be|}
          x prod_{H-chords al, marks e}           (-1)^{b_e[al] [e in c_al]}

where b_e[al] = parity of the number of e-marked fermion positions
strictly inside chord al's span.  The disorder average therefore needs
only E over independent uniform p-subsets of a parity product.  That is an
exact integer sum over Venn-cell occupancy profiles, memoized per
crossing/mark signature.  In chord units (divide a K-H word by
sigma_H^K = (C(N,p) sigma^2)^{K/2}) sigma cancels entirely:

    t4bar(a,b,c,d) = sum_matchings E[sign],    E exact in QQ.

The mu_k assembly (binomial ad-expansion, symmetrized, both orderings)
is IMPORTED from duality/chord_moments.py, so the finite-N and N =
infinity computations share it by construction.  t4bar -> chord t4 as
N -> infinity at fixed lambda, and every mu_k>=1 must fall to zero.

Scope: n <= 3 and moments mu_0, mu_2, mu_4 (words up to K = 10 H-nodes,
5 subsets).  Every signature goes through an exact element-wise DP, whose
state is the per-subset fill counts, with integers via CRT residue
tables.  Its cost is O(N (p+1)^m 2^m) per crossing set.  The closed
m <= 3 Venn-profile sum is retained as the in-tree reference the DP is
tested against.  The flavor sum gives exact finite-N
identities used as end-to-end gates: A_1 = -2pH per instance, so

    (N - 1) mu_k^pair(n=1) + mu_k^diag(n=1) = 0   (k = 2 and 4), every N,

with the diagonal channel (all four fermions on one flavor) computed by
the same machinery.

Validation (test_exact_pair_moments.py + an adversarial referee pass,
2026-08-19): literal Isserlis enumeration over ALL subset
tuples with Pauli-string traces at small (N, p) (exact match);
disorder-averaged dense ED (within errors); the n = 1 identity at
machine precision; and the chord limit (t4bar -> chord t4).  The referee
pass checked the sign law per configuration on ~5e5 adversarial
subset tuples, the factorization and canonicalization against naive
enumeration on all 311 possible m <= 3 signatures, and three independent
routes through p = 10.

Precision note, and the practical cap.  The per-word values are exact
rationals, but the mu_k assembly is float, so the alternating sums
inherit rounding noise that grows fast with the gross magnitude.  The
A_1 flavor-sum gate measures it directly.  On the committed scan:

    N          64      96     128     192     256
    mu_2 gate  2e-14   7e-13  1e-12   2e-12   1.7e-11
    mu_4 gate  2e-12   3e-12  9e-12   1.2e-10 9.1e-10

The mu_4 gate is the binding one: at N = 256 it sits at 9.1e-10 against
the 1e-9 threshold asserted in test_committed_campaign -- a 9% margin.
So N = 256 is the honest cap for the float assembly, NOT a round number
some way beyond it.  Extending the scan (or adding mu_6) requires
carrying _mu_ordered / _mu_ordered_swapped in Fractions; the word values
themselves are already exact and need no change.

Output: exact_pair_moments.json (fixed-p scan + the double-scaling line).
"""

from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "duality"))

from chord_moments import _mu_ordered, _mu_ordered_swapped, _trace_part  # noqa: E402


# ---------------------------------------------------------------------------
# diagrams: word layout, matchings, crossing/mark signatures
# ---------------------------------------------------------------------------
def _word_symbols(a: int, b: int, c: int, d: int, marks: tuple) -> list:
    """Circle word: H-blocks of sizes (a, b, c, d) separated by four marked
    fermion positions with distinguished-element labels `marks` (length 4;
    pair channel ('i','i','j','j'), diagonal ('i','i','i','i'))."""
    return (['H'] * a + [marks[0]] + ['H'] * b + [marks[1]]
            + ['H'] * c + [marks[2]] + ['H'] * d + [marks[3]])


def _matchings(items):
    if not items:
        yield []
        return
    first = items[0]
    for k in range(1, len(items)):
        rest = items[1:k] + items[k + 1:]
        for rem in _matchings(rest):
            yield [(first, items[k])] + rem


def _crosses(ch1, ch2) -> bool:
    (a, b), (c, d) = sorted(ch1), sorted(ch2)
    return (a < c < b) != (a < d < b)


def _signature(symbols, matching):
    """(m, S, marks_bits): chord count, frozenset of crossing chord-index
    pairs, and per-distinguished-element inside-parity bit vectors."""
    chords = [tuple(sorted(pair)) for pair in matching]
    m = len(chords)
    S = frozenset((al, be) for al in range(m) for be in range(al + 1, m)
                  if _crosses(chords[al], chords[be]))
    labels = sorted({s for s in symbols if s != 'H'})
    marks_bits = []
    for lab in labels:
        pos = [k for k, s in enumerate(symbols) if s == lab]
        bits = tuple(sum(1 for q in pos if lo < q < hi) % 2
                     for (lo, hi) in chords)
        marks_bits.append(bits)
    return m, S, tuple(marks_bits)


# ---------------------------------------------------------------------------
# the exact subset expectation per signature (integer profile sums)
# ---------------------------------------------------------------------------
_EXP_CACHE: dict = {}


def signature_expectation(N: int, p: int, m: int, S: frozenset,
                          marks_bits: tuple) -> Fraction:
    """E over m independent uniform p-subsets of [N] of

        prod_{(al,be) in S} (-1)^{|c_al & c_be|}
        x prod_{e, al} (-1)^{marks_bits[e][al] * [x_e in c_al]}

    with one distinguished element x_e per mark label.  Exact rational.
    The expectation factorizes over connected components of the crossing
    graph (the subsets are independent; a distinguished element whose
    bits vanish inside a component is just a generic element there).  So
    only genuinely 3-connected crossing components pay the full
    Venn-profile sum."""
    # union-find over the crossing graph
    parent = list(range(m))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (al, be) in S:
        ra, rb = find(al), find(be)
        if ra != rb:
            parent[ra] = rb
    comps: dict = {}
    for al in range(m):
        comps.setdefault(find(al), []).append(al)
    out = Fraction(1)
    for comp in comps.values():
        idx = sorted(comp)
        remap = {al: k for k, al in enumerate(idx)}
        S_c = frozenset(tuple(sorted((remap[a], remap[b])))
                        for (a, b) in S if a in remap and b in remap)
        bits_c = tuple(tuple(bits[al] for al in idx) for bits in marks_bits
                       if any(bits[al] for al in idx))
        out *= _expectation_core(N, p, len(idx), S_c, bits_c)
    return out


_CANON_CACHE: dict = {}


def _canonical(m: int, S: frozenset, marks_bits: tuple):
    """Canonical relabeling of the chords (m! permutations; m <= 5 in
    scope).  Pure in its arguments, and called once per word per signature,
    so it is memoized."""
    ck = (m, S, marks_bits)
    if ck in _CANON_CACHE:
        return _CANON_CACHE[ck]
    from itertools import permutations
    best = None
    for perm in permutations(range(m)):
        S2 = frozenset(tuple(sorted((perm[a], perm[b]))) for (a, b) in S)
        b2 = tuple(sorted(tuple(bits[perm.index(k)] for k in range(m))
                          for bits in marks_bits))
        key = (tuple(sorted(S2)), b2)
        if best is None or key < best:
            best = key
    _CANON_CACHE[ck] = best
    return best


# -- exact signed assignment counts by dynamic programming ------------------
# One element at a time: each generic element picks a membership pattern
# T subset of the m chords with sign (-1)^(number of S-pairs inside T).
# The DP state is the per-chord fill count.  Exact integers via CRT over
# 25-bit primes (numpy int64 residue arrays), so any m is feasible.  The
# direct Venn-profile sum is O((p+1)^(2^m - m - 1)) and stops at m = 3.

def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % a == 0:
            return n == a
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _first_primes_below(bound: int, count: int) -> list:
    out, n = [], bound - 1
    while len(out) < count:
        if _is_prime(n):
            out.append(n)
        n -= 2 if n % 2 else 1
    return out


_PRIMES = _first_primes_below(1 << 25, 16)
_DP_CACHE: dict = {}


def _dp_tables(n_elems: int, p: int, m: int, S_key: tuple):
    """Residue tables of the DP final state on the full (p+1)^m grid, one
    per prime.  It computes ALL targets at once (shared across membership
    patterns and across every word with the same crossing set S)."""
    import numpy as np
    key = (n_elems, p, m, S_key)
    if key in _DP_CACHE:
        return _DP_CACHE[key]
    # |signed count| <= prod_al C(n_elems, t_al) with every t_al <= p, and
    # C(n, t) peaks at t = n//2.  So the per-chord factor must be
    # C(n, min(p, n//2)), NOT C(n, p) (unsound for p > n//2)
    bound = 2 * (2 * comb(n_elems, min(p, n_elems // 2))) ** m
    signs = []
    for T in range(1 << m):
        par = sum(1 for (al, be) in S_key if (T >> al) & 1 and (T >> be) & 1)
        signs.append(-1 if par & 1 else 1)
    tables, modulus, used = [], 1, []
    for P in _PRIMES:
        if modulus > 2 * bound:
            break
        state = np.zeros((p + 1,) * m, dtype=np.int64)
        state[(0,) * m] = 1
        for _ in range(n_elems):
            new = state.copy()                       # T = empty set, sign +1
            for T in range(1, 1 << m):
                src = tuple(slice(0, p) if (T >> al) & 1 else slice(None)
                            for al in range(m))
                dst = tuple(slice(1, None) if (T >> al) & 1 else slice(None)
                            for al in range(m))
                if signs[T] > 0:
                    new[dst] += state[src]
                else:
                    new[dst] -= state[src]
            state = new % P
        tables.append(state)
        used.append(P)
        modulus *= P
    if modulus <= 2 * bound:
        raise ValueError("_dp_tables: CRT prime pool exhausted -- "
                         f"modulus {modulus} <= 2*bound {2 * bound}")
    # CRT coefficients Mi * (Mi^-1 mod P) do not depend on the target, so
    # fold them in here instead of recomputing a big-int modular inverse on
    # every _dp_count call.
    coeffs = [(modulus // P) * pow(modulus // P, -1, P) for P in used]
    out = (tables, coeffs, modulus)
    _DP_CACHE[key] = out
    return out


_S_CANON_CACHE: dict = {}


def _canonical_crossing_set(m: int, S_key: tuple):
    """(S_star, inv) with S_star the lexicographically smallest relabeling of
    the crossing set and inv the inverse permutation, so a target for S maps
    to target_star[b] = target[inv[b]].

    The DP table depends on the crossing set ONLY up to chord relabeling, so
    isomorphic signatures share one table instead of building ~5 redundant
    copies apiece."""
    if (m, S_key) in _S_CANON_CACHE:
        return _S_CANON_CACHE[(m, S_key)]
    from itertools import permutations
    best = None
    for perm in permutations(range(m)):
        cand = tuple(sorted(tuple(sorted((perm[a], perm[b])))
                            for (a, b) in S_key))
        if best is None or cand < best[0]:
            inv = [0] * m
            for a, b in enumerate(perm):
                inv[b] = a
            best = (cand, tuple(inv))
    _S_CANON_CACHE[(m, S_key)] = best
    return best


def _dp_count(n_elems: int, p: int, m: int, S_key: tuple,
              target: tuple) -> int:
    """Signed count over assignments of n_elems generic elements to
    membership patterns T (one per element), with exactly target[al]
    elements in chord al, weighted by prod (-1)^(#S-pairs inside T).
    Exact integer by CRT over the cached residue tables."""
    S_star, inv = _canonical_crossing_set(m, S_key)
    target = tuple(target[inv[b]] for b in range(m))
    tables, coeffs, modulus = _dp_tables(n_elems, p, m, S_star)
    total = 0
    for tab, c in zip(tables, coeffs):
        total += int(tab[target]) * c
    total %= modulus
    if total > modulus // 2:
        total -= modulus
    return total


def _expectation_core(N: int, p: int, m: int, S: frozenset,
                      marks_bits: tuple) -> Fraction:
    """One connected crossing component: membership patterns of the
    distinguished elements outside, the exact DP count inside."""
    key = (N, p, m, _canonical(m, S, marks_bits), "dp")
    if key in _EXP_CACHE:
        return _EXP_CACHE[key]
    n_dist = len(marks_bits)
    S_key = tuple(sorted(S))
    n_generic = N - n_dist
    total = 0
    for patterns in product(*[list(product((0, 1), repeat=m))
                              for _ in range(n_dist)]):
        base_parity = 0
        for (al, be) in S:
            for I in patterns:
                base_parity ^= I[al] & I[be]
        for e, bits in enumerate(marks_bits):
            I = patterns[e]
            for al in range(m):
                base_parity ^= bits[al] & I[al]
        target = tuple(p - sum(I[al] for I in patterns) for al in range(m))
        if any(t < 0 for t in target):
            continue
        cnt = _dp_count(n_generic, p, m, S_key, target)
        total += -cnt if base_parity else cnt
    out = Fraction(total, comb(N, p) ** m)
    _EXP_CACHE[key] = out
    return out


def _expectation_core_profile(N: int, p: int, m: int, S: frozenset,
                              marks_bits: tuple) -> Fraction:
    """The original Venn-profile sum (m <= 3 only) -- retained as the
    in-tree reference the DP path is tested against."""
    key = (N, p, m, _canonical(m, S, marks_bits), "profile")
    if key in _EXP_CACHE:
        return _EXP_CACHE[key]
    n_dist = len(marks_bits)
    free_cells = [frozenset(T) for r in range(2, m + 1)
                  for T in combinations(range(m), r)]
    total = 0
    n_generic = N - n_dist
    # membership patterns of each distinguished element in the m subsets
    for patterns in product(*[list(product((0, 1), repeat=m))
                              for _ in range(n_dist)]):
        base_parity = 0
        for (al, be) in S:
            for I in patterns:
                base_parity ^= I[al] & I[be]
        for e, bits in enumerate(marks_bits):
            I = patterns[e]
            for al in range(m):
                base_parity ^= bits[al] & I[al]
        dist_in = [sum(I[al] for I in patterns) for al in range(m)]
        # occupancy of the free (multi-subset) cells
        ranges = [range(p + 1)] * len(free_cells)
        for occ in product(*ranges):
            per_chord = [sum(o for o, T in zip(occ, free_cells) if al in T)
                         for al in range(m)]
            singles = [p - per_chord[al] - dist_in[al] for al in range(m)]
            if any(s < 0 for s in singles):
                continue
            used = sum(occ) + sum(singles)
            if used > n_generic:
                continue
            parity = base_parity
            for (al, be) in S:
                parity ^= (sum(o for o, T in zip(occ, free_cells)
                               if al in T and be in T) & 1)
            count, rem = 1, n_generic
            for s in list(occ) + singles:
                count *= comb(rem, s)
                rem -= s
            total += -count if parity else count
    out = Fraction(total, comb(N, p) ** m)
    _EXP_CACHE[key] = out
    return out


# ---------------------------------------------------------------------------
# normalized disorder-averaged trace words (chord units: /sigma_H^K)
# ---------------------------------------------------------------------------
def _psi_psi_sign(symbols) -> int:
    """The fermion-fermion crossing sign the H-chord law does not carry:
    for two distinct single-Majorana labels with two positions each,
    (-1) if their chords cross.  Supported mark patterns: one label with
    2 or 4 positions (sign +1), or two labels with 2 positions each.
    Anything else (odd multiplicities, > 2 labels) has no law here and
    raises."""
    pos = {}
    for k, s in enumerate(symbols):
        if s != 'H':
            pos.setdefault(s, []).append(k)
    counts = sorted(len(v) for v in pos.values())
    if counts in ([2], [4]):
        return 1
    if counts == [2, 2]:
        (x1, x2), (y1, y2) = (sorted(v) for v in pos.values())
        return -1 if (x1 < y1 < x2) != (x1 < y2 < x2) else 1
    raise ValueError(f"unsupported mark pattern {sorted(pos)} with "
                     f"multiplicities {counts}: each label needs an even "
                     "count, at most two labels")


def t4bar(N: int, p: int, a: int, b: int, c: int, d: int,
          marks=('i', 'i', 'j', 'j')) -> float:
    """E<Tr[H^a psi H^b psi H^c psi H^d psi]>/dim / sigma_H^(a+b+c+d),
    fermion flavors per `marks` (two labels of two fermions each, in any
    arrangement, or one label; other patterns raise).  Odd H-count words
    vanish (Isserlis)."""
    K = a + b + c + d
    if K % 2:
        return 0.0
    symbols = _word_symbols(a, b, c, d, marks)
    psi_sign = _psi_psi_sign(symbols)
    hpos = [k for k, s in enumerate(symbols) if s == 'H']
    tot = Fraction(0)
    for matching in _matchings(hpos):
        m, S, bits = _signature(symbols, matching)
        tot += signature_expectation(N, p, m, S, bits)
    return psi_sign * float(tot)


def t2bar(N: int, p: int, a: int, b: int) -> float:
    """E<Tr[H^a psi_i H^b psi_i]>/dim / sigma_H^(a+b)."""
    K = a + b
    if K % 2:
        return 0.0
    symbols = ['H'] * a + ['i'] + ['H'] * b + ['i']
    hpos = [k for k, s in enumerate(symbols) if s == 'H']
    tot = Fraction(0)
    for matching in _matchings(hpos):
        m, S, bits = _signature(symbols, matching)
        tot += signature_expectation(N, p, m, S, bits)
    return float(tot)


class PairTab:
    """Finite-N drop-in for chord_moments.MomentTable: normalized
    disorder-averaged words, i != j pair channel."""

    def __init__(self, N: int, p: int, marks=('i', 'i', 'j', 'j')):
        if p % 2:
            raise ValueError("p must be even")
        self.N, self.p, self.marks = N, p, marks
        self._t4: dict = {}
        self._t2: dict = {}

    def t4(self, a, b, c, d):
        key = (a, b, c, d)
        if key not in self._t4:
            self._t4[key] = t4bar(self.N, self.p, a, b, c, d, self.marks)
        return self._t4[key]

    def t2(self, a, b):
        if (a, b) not in self._t2:
            self._t2[(a, b)] = t2bar(self.N, self.p, a, b)
        return self._t2[(a, b)]


def channel_pair_moments(N: int, p: int, n: int,
                         marks=('i', 'i', 'j', 'j'),
                         ks=(0, 2)) -> dict:
    """Exact finite-N symmetrized traceless mu_k of the A_n channel (per
    flavor pair for the default marks).  Scope: K = k + 2n <= 10 H-nodes,
    i.e. at most 5 subsets per word (the DP count handles any m; the
    matching enumeration grows as (K-1)!!)."""
    if max(ks) + 2 * n > 10:
        raise ValueError("scope: k + 2n <= 10 (largest supported word)")
    tab = PairTab(N, p, marks)
    out = {}
    for k in ks:
        v1, g1 = _mu_ordered(tab, n, k)
        v2, g2 = _mu_ordered_swapped(tab, n, k)
        out[f"mu{k}"] = 0.5 * (v1 + v2)
        out[f"gross{k}"] = 0.5 * (g1 + g2)
    trX = _trace_part(tab, n)
    out["mu0"] -= trX * trX
    out["w2"] = (out["mu2"] / out["mu0"]) if out["mu0"] != 0 else float("nan")
    return out


# ---------------------------------------------------------------------------
# literal ground truth (small N, p): Isserlis over ALL subset tuples with
# Pauli-string traces -- shares nothing with the profile sums
# ---------------------------------------------------------------------------
def literal_t4bar(N: int, p: int, a: int, b: int, c: int, d: int,
                  marks=('i', 'i', 'j', 'j')) -> float:
    from pauli_strings import _I_POW
    from exact_wick import _trace_of_product
    K = a + b + c + d
    if K % 2:
        return 0.0
    elem = {'i': (0,), 'j': (1,)}
    symbols = _word_symbols(a, b, c, d, marks)
    hpos = [k for k, s in enumerate(symbols) if s == 'H']
    subsets = list(combinations(range(N), p))
    pref = _I_POW[((p // 2) * K) % 4]        # i^{p/2} per Gamma
    tot = 0.0
    for matching in _matchings(hpos):
        chord_of = {}
        for ci, (q1, q2) in enumerate(matching):
            chord_of[q1] = ci
            chord_of[q2] = ci
        for combo in product(subsets, repeat=len(matching)):
            word = [combo[chord_of[k]] if s == 'H' else elem[s]
                    for k, s in enumerate(symbols)]
            t = pref * _trace_of_product(N, word)
            tot += t.real
    return tot / comb(N, p) ** (K // 2)


# ---------------------------------------------------------------------------
# campaign
# ---------------------------------------------------------------------------
def run(out_path=None) -> dict:
    """Fixed p = 4, N rising (the lambda -> 0 corner + the 1/N tail), and
    the double-scaling line lambda_chord = 2 (p^2 proportional to N):
    exact mu_0, mu_2, mu_4 of the n = 1, 2, 3 pair channels, plus the
    n = 1 flavor-sum identity residuals (mu_2 and mu_4) as per-point
    gates."""
    t0 = time.time()
    points = []
    for (p, N) in [(4, 8), (4, 10), (4, 12), (4, 14), (4, 16), (4, 20),
                   (4, 24), (4, 32), (4, 48), (4, 64), (4, 96), (4, 128),
                   (4, 192), (4, 256)]:
        row = dict(p=p, N=N, lam_chord=2.0 * p * p / N, channels={})
        for n in (1, 2, 3):
            row["channels"][f"n={n}"] = channel_pair_moments(
                N, p, n, ks=(0, 2, 4))
        diag = channel_pair_moments(N, p, 1, marks=('i', 'i', 'i', 'i'),
                                    ks=(0, 2, 4))
        pair = row["channels"]["n=1"]
        for k in (2, 4):
            scale = max(abs(diag[f"mu{k}"]), 1e-300)
            row[f"a1_identity_residual_mu{k}"] = abs(
                (N - 1) * pair[f"mu{k}"] + diag[f"mu{k}"]) / scale
        row["a1_identity_residual"] = row["a1_identity_residual_mu2"]
        points.append(row)
        c2, c3 = row["channels"]["n=2"], row["channels"]["n=3"]
        print(f"p={p} N={N:4d} lam={row['lam_chord']:.3f}: "
              f"mu2(n=2)={c2['mu2']:+.6e}  mu4(n=2)={c2['mu4']:+.4e}  "
              f"mu2(n=3)={c3['mu2']:+.4e}  "
              f"A1={row['a1_identity_residual']:.0e}/"
              f"{row['a1_identity_residual_mu4']:.0e}", flush=True)
        # the caches key on n_elems = N - n_dist, so nothing is reused
        # across N points -- drop them to keep the sweep's memory flat
        _DP_CACHE.clear()
        _EXP_CACHE.clear()
    ds_points = []
    for (p, N) in [(4, 16), (6, 36), (8, 64), (10, 100)]:
        row = dict(p=p, N=N, lam_chord=2.0 * p * p / N,
                   n2=channel_pair_moments(N, p, 2, ks=(0, 2, 4)),
                   n3=channel_pair_moments(N, p, 3, ks=(0, 2)))
        ds_points.append(row)
        print(f"double-scaling lam=2: p={p} N={N}: "
              f"mu2(n=2)={row['n2']['mu2']:+.6e}  "
              f"mu2(n=3)={row['n3']['mu2']:+.6e}", flush=True)
        _DP_CACHE.clear()
        _EXP_CACHE.clear()
    payload = dict(
        provenance="exact_pair_moments.py -- exact finite-N Wick pair "
                   "moments of the A_n channels (CRT residue-table DP, "
                   "no ED): the measured approach to the conserved-channel "
                   "theorem's N = infinity zero (chord units, per flavor "
                   "pair; assembly imported from duality/chord_moments.py)",
        runtime_s=round(time.time() - t0, 1),
        fixed_p=points, double_scaling=ds_points)
    out = Path(out_path) if out_path else HERE / "exact_pair_moments.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


if __name__ == "__main__":
    run()
