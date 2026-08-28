#!/usr/bin/env python3
r"""
The flavor-DIAGONAL channel of the A_n bilinears at N = infinity.
=================================================================

The conserved-channel theorem (chord_charges.py) covers the i != j PAIR
amplitude of A_n = sum_i psi_i (ad_H)^n psi_i.  That is only half the
channel.  The full one is

    <A_n A_n> = N <X_i X_i> + N(N-1) <X_i X_j>,   X_i = psi_i ad_H^n psi_i,

and the diagonal piece <X_i X_i> is the same order.  It is the part the
ED bench says dominates at accessible sizes, and the part that does NOT
vanish for n >= 2.

This module computes the diagonal piece exactly at N = infinity.  That
closes the chord side of the gap between the theorem and the bench.

THE DIAGONAL CHORD LAW
----------------------
The diagonal word Tr[H^a psi_i H^b psi_i H^c psi_i H^d psi_i] has four
insertions of the SAME Majorana.  Start from the exact finite-N sign law
in syk-self-averaging/exact_pair_moments.py, which is validated against
literal Pauli-string traces:

    trace = prod_{crossing H-chords} (-1)^{|c_al & c_be|}
            x prod_{H-chords al} (-1)^{b[al] * [i in c_al]} ,

where b[al] is the PARITY of the number of the four mark positions lying
strictly inside chord al's span.

The subsets c_al are independent across chords.  So at N = infinity the
expectation factorizes.  You get the standard crossing weight q per H-H
crossing (q = e^{-lambda}, lambda = 2p^2/N), times, per H-chord,

    E[(-1)^{b[al] [i in c_al]}] = (1 - 2p/N)^{b[al]} -> q^{Delta_psi b[al]},

with Delta_psi = 1/p.

So the N = infinity diagonal amplitude is the PARITY-LAW diagram sum.  An
H-chord weighs q^{Delta_psi} if and only if it crosses the four-mark set
an ODD number of times.  Otherwise it weighs 1.

TWO HYPOTHESES THIS REPLACES.  Each is refuted by an exact anchor.

  * "Sum over the three Wick pairings of the four marks, with the
    standard q^{Delta} per matter-chord crossing."  At Delta -> 0 that
    sum gives 3 <Tr H^K>.  The true diagonal amplitude is <Tr H^K>
    exactly.  There is no pairing sum.  The parity weight is itself
    pairing-independent: for ANY pairing of the four marks into two
    chords, b[al] = (number of crossings with the two chords) mod 2.  So
    the three pairings collapse into one diagram family.
  * "One pairing, with weight q^{2 Delta} for chords crossing both matter
    chords."  The finite-N law depends on the marks only through the
    parity bits b[al], which are 0 or 1.  So a chord with TWO marks
    inside weighs exactly 1 at EVERY finite N.  For example
    Tr[H psi psi H psi psi] = Tr[H^2] gives t4(1,0,1,0) = 1 exactly at
    every (N, p), which is asserted.  The q^{2 Delta} weight predicts
    q^{2 Delta}, which is not 1.

psi^2 = 1 comes out automatically.  Two adjacent marks flip the flux
twice with no H-node between them, so a mark pair acts as the identity.
That is the same mechanism that collapses t4(a, 0, c, 0) to
<Tr H^{a+c}> exactly.

EVALUATION
----------
Two independent routes, asserted equal to machine precision.

  * `diag_word_amplitude`, an ordered-tag flux WALKER.  The state is the
    ordered tuple of open-chord tags, where a chord's tag is the mark
    parity, its "flux", at the moment it opened.  A mark flips the flux.
    Closing the chord at depth j costs q^{(number of chords opened after
    it and still open)} times q^{Delta (flux XOR tag_j)}.  This is exact,
    with no chord-lattice truncation, at a cost of about 2^{K/2} per
    word.
  * `bruteforce_parity_word`, a literal sum over the (K-1)!! matchings
    with the parity weight.  It mirrors chord.bruteforce_word.

The walker's 2-mark restriction is the ordinary one-matter-chord object.
It is asserted against the validated transfer-matrix engine
(chord_moments.MomentTable.t2) to machine precision.  That anchors the
crossing bookkeeping to the engine the whole module chain trusts.

The mu_k assembly is IMPORTED from chord_moments.py: the binomial
ad-expansion, symmetrized, both orderings, with mu_0 tracelessed by
<X>^2.  So the diagonal, the pair (chord_moments), and the finite-N
(exact_pair_moments) computations share it by construction.

GATES (test_diagonal_chord.py)
------------------------------
G1, the A_1 identity.  A_1 = -2pH per instance forces
(N-1) mu_k^pair(N) + mu_k^diag(N) = 0 exactly at every N, for k = 2 and
4.  At strict N = infinity every per-flavor piece on the double-scaling
line vanishes like 1/p.  So the identity's chord shadow is a CONVERGENCE
statement: diag_chord(lambda, 1/p) must equal -(N-1) mu_k^pair(N, p) up
to O(1/N), with the residual falling along the line.

(The exact finite-N identity itself is gated at machine precision in
test_exact_pair_moments.py.  A strict-N = infinity machine-precision
full-channel zero would need the pair channel's 1/N coefficient, which is
pair_moment_laws.py territory, not a chord-model object.)

G2, convergence.  Along lambda_chord = 2, the exact finite-N diagonal
moments at (p, N) = (4,16), (6,36), (8,64), (10,100) approach
channel_moments_diag(q = e^-2, Delta_psi = 1/p), with a deviation falling
monotonically.  It is O(1/N) with a large coefficient; the matter weight
itself contributes (1-2/p) against e^{-2/p}, which is O(1/N) on the line.

At FIXED p = 4, with the lattice lambda = 32/N matched exactly, the
deviation falls like 1/N^2.  Measured: 6.5e-2 -> 3.3e-4 over
N = 64 -> 1024.  The reason is that the O(1/N) mismatch of the weights is
then absorbed into q itself.

G3, the payoff.  The diagonal-channel width w2 = mu2/mu0 at
lambda_chord = 2.  w2_diag(Delta = 1/p) for p = 4..10 tracks the rising
finite-N FULL-channel series 1.044 -> 1.194, and its Delta -> 0 limit is
the number this module pins.

Scope, stated plainly.  The N = infinity FULL-channel width also needs
the pair channel's 1/N coefficient.  Both pieces are O(1/p) on the line,
and neither dominates in the strict limit.  This module supplies the
diagonal half exactly.

Output: diagonal_chord.json -- the campaign over the chord_moments lambda
grid, plus the three gate tables.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

from chord_moments import _mu_ordered, _mu_ordered_swapped, _trace_part

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# the parity-law word amplitudes
# ---------------------------------------------------------------------------
def word_diag(a: int, b: int, c: int, d: int) -> list:
    """Circle word for the diagonal 4-mark object: H-blocks (a, b, c, d)
    separated by four marks of the SAME Majorana ('P')."""
    return (['H'] * a + ['P'] + ['H'] * b + ['P']
            + ['H'] * c + ['P'] + ['H'] * d + ['P'])


def diag_word_amplitude(q: float, dpsi: float, word: list) -> float:
    """<Tr[word]> in chord units by the ordered-tag flux walker (exact).

    word: list over {'H', 'P'}; every 'P' flips the flux s in {0, 1} (the
    number of marks passed, mod 2 -- an even mark count is required so the
    trace closes).  A chord opened at flux t and closed at flux s has mark
    parity b = s XOR t and weighs q^{dpsi b}; its H-H crossing weight is
    q^{# chords opened after it that are still open at its close}, which
    reproduces the [n]_q rule chord-by-chord.  Cyclic invariance holds
    because matchings of the cut line biject with matchings of the circle
    and both the interleaving relation and the span parity are properties
    of the cyclic order (tested)."""
    n_marks = sum(1 for w in word if w == 'P')
    if n_marks % 2:
        raise ValueError("odd mark count cannot close the trace")
    if sum(1 for w in word if w == 'H') % 2:
        return 0.0
    states = {(): 1.0}
    s = 0
    for w in word:
        if w == 'P':
            s ^= 1
            continue
        new: dict = {}
        for tags, amp in states.items():
            opened = tags + (s,)
            new[opened] = new.get(opened, 0.0) + amp
            L = len(tags)
            for j in range(L):
                wgt = q ** (L - 1 - j) * q ** (dpsi * (s ^ tags[j]))
                closed = tags[:j] + tags[j + 1:]
                new[closed] = new.get(closed, 0.0) + amp * wgt
        states = new
    return states.get((), 0.0)


def bruteforce_parity_word(q: float, dpsi: float, word: list) -> float:
    """The same amplitude by literal (K-1)!! matching enumeration with the
    parity weight -- shares nothing with the walker but the microscopic
    law (mirrors chord.bruteforce_word)."""
    marks = [i for i, w in enumerate(word) if w == 'P']
    hpos = [i for i, w in enumerate(word) if w == 'H']
    if len(hpos) % 2:
        return 0.0

    def matchings(items):
        if not items:
            yield []
            return
        first = items[0]
        for k in range(1, len(items)):
            rest = items[1:k] + items[k + 1:]
            for rem in matchings(rest):
                yield [(first, items[k])] + rem

    total = 0.0
    for hm in matchings(hpos):
        w = 1.0
        spans = [tuple(sorted(ch)) for ch in hm]
        for i, (ai, bi) in enumerate(spans):
            for (aj, bj) in spans[i + 1:]:
                if (ai < aj < bi) != (ai < bj < bi):
                    w *= q
            if sum(1 for m in marks if ai < m < bi) % 2:
                w *= q ** dpsi
        total += w
    return total


# ---------------------------------------------------------------------------
# the moment table and channel assembly (assembly shared with chord_moments)
# ---------------------------------------------------------------------------
class DiagonalMomentTable:
    """Drop-in for chord_moments.MomentTable: normalized N = infinity
    words of the flavor-DIAGONAL channel (all four marks one Majorana),
    evaluated by the exact parity-law walker.  t2 (two marks) is the
    ordinary single-matter-chord object, anchored to the engine in the
    tests."""

    def __init__(self, q: float, delta_psi: float):
        self.q, self.d = q, delta_psi
        self._t4: dict = {}
        self._t2: dict = {}

    def t4(self, a: int, b: int, c: int, d: int) -> float:
        key = (a, b, c, d)
        if key not in self._t4:
            self._t4[key] = diag_word_amplitude(self.q, self.d,
                                                word_diag(a, b, c, d))
        return self._t4[key]

    def t2(self, a: int, b: int) -> float:
        if (a, b) not in self._t2:
            word = ['H'] * a + ['P'] + ['H'] * b + ['P']
            self._t2[(a, b)] = diag_word_amplitude(self.q, self.d, word)
        return self._t2[(a, b)]


def channel_moments_diag(q: float, delta_psi: float, n: int,
                         ks=(0, 2, 4, 6)) -> dict:
    """Exact N = infinity symmetrized traceless mu_k of the DIAGONAL A_n
    channel, per flavor, chord units -- the same assembly, subtraction,
    and invariants as chord_moments.channel_moments uses for the pair."""
    tab = DiagonalMomentTable(q, delta_psi)
    mu, gross = {}, {}
    for k in ks:
        v1, g1 = _mu_ordered(tab, n, k)
        v2, g2 = _mu_ordered_swapped(tab, n, k)
        mu[k] = 0.5 * (v1 + v2)
        gross[k] = 0.5 * (g1 + g2)
    trX = _trace_part(tab, n)
    mu[0] = mu[0] - trX * trX
    out = {f"mu{k}": mu[k] for k in ks}
    out.update({f"gross{k}": gross[k] for k in ks if k})
    out["w2"] = (mu[2] / mu[0]) if mu.get(0) else float("nan")
    if 4 in mu:
        out["kurt"] = (mu[4] * mu[0] / mu[2] ** 2) if mu.get(2) else float("nan")
    return out


# ---------------------------------------------------------------------------
# finite-N ground truth (exact Wick machinery from syk-self-averaging)
# ---------------------------------------------------------------------------
def _finite_n(N: int, p: int, n: int, marks, ks=(0, 2)) -> dict:
    sys.path.insert(0, str(HERE.parent / "syk-self-averaging"))
    import exact_pair_moments as epm
    out = epm.channel_pair_moments(N, p, n, marks=marks, ks=ks)
    epm._DP_CACHE.clear()
    epm._EXP_CACHE.clear()
    return out


def finite_n_diag(N: int, p: int, n: int, ks=(0, 2)) -> dict:
    return _finite_n(N, p, n, ('i', 'i', 'i', 'i'), ks)


def finite_n_pair(N: int, p: int, n: int, ks=(0, 2)) -> dict:
    return _finite_n(N, p, n, ('i', 'i', 'j', 'j'), ks)


# ---------------------------------------------------------------------------
# campaign
# ---------------------------------------------------------------------------
def w2_delta0(lam: float,
              deltas=(1.0 / 8, 1.0 / 12, 1.0 / 16, 1.0 / 24, 1.0 / 32,
                      1.0 / 48),
              n: int = 2) -> dict:
    """The Delta_psi -> 0 limit of the diagonal width w2 = mu2/mu0 at
    fixed lambda.  mu0 and mu2 are both O(Delta), so the ratio has a
    finite limit; w2(Delta) is fit linearly in Delta from the two
    smallest points (Richardson), with the fit spread over the point set
    recorded as the honest error."""
    q = math.exp(-lam)
    vals = [(d, channel_moments_diag(q, d, n, ks=(0, 2))["w2"])
            for d in sorted(deltas, reverse=True)]
    (d1, w1), (d0, w0) = vals[-2], vals[-1]
    extrap = w0 + (w0 - w1) * d0 / (d1 - d0)
    pairs = list(zip(vals[:-1], vals[1:]))
    spread = max(abs((wb + (wb - wa) * db / (da - db)) - extrap)
                 for (da, wa), (db, wb) in pairs)
    return dict(lam=lam, q=q, n=n, w2_of_delta=[[d, w] for d, w in vals],
                w2_delta0=extrap, extrap_spread=spread)


def run(p: int = 4,
        lams=(4.0, 3.2, 2.6667, 2.2857, 2.0, 1.7778, 1.6, 1.0, 0.5, 0.25,
              0.1),
        out_path=None) -> dict:
    """The campaign: the diagonal channel over the chord_moments lambda
    grid (n = 1..5, mu_0..mu_6), plus the three gate tables of the module
    docstring.  Runtime: minutes (the K = 16 walker words dominate)."""
    t0 = time.time()
    dpsi = 1.0 / p
    points = []
    for lam in lams:
        q = math.exp(-lam)
        row = dict(lam=lam, q=q, p=p, delta_psi=dpsi, channels={})
        for n in (1, 2, 3, 4, 5):
            row["channels"][f"n={n}"] = channel_moments_diag(q, dpsi, n)
        points.append(row)
        c2 = row["channels"]["n=2"]
        print(f"lam={lam}: diag n=2 mu2={c2['mu2']:+.6e} w2={c2['w2']:.4f}",
              flush=True)

    # G2: the lambda_chord = 2 line, chord vs exact finite-N diagonal
    line = []
    for (pp, N) in [(4, 16), (6, 36), (8, 64), (10, 100)]:
        ch = channel_moments_diag(math.exp(-2.0), 1.0 / pp, 2, ks=(0, 2))
        fN = finite_n_diag(N, pp, 2, ks=(0, 2))
        line.append(dict(p=pp, N=N,
                         chord_mu2=ch["mu2"], finiteN_mu2=fN["mu2"],
                         chord_w2=ch["w2"], finiteN_w2=fN["w2"],
                         rel_dev_mu2=abs(ch["mu2"] - fN["mu2"])
                         / abs(ch["mu2"])))
        print(f"line p={pp} N={N}: chord mu2={ch['mu2']:+.5f} "
              f"finiteN mu2={fN['mu2']:+.5f} "
              f"rel dev={line[-1]['rel_dev_mu2']:.3f}", flush=True)

    # G2 sharp: fixed p, lattice lambda matched exactly, at the WORD level
    # (dev ~ 1/N^2 there; the assembled mu_2 at fixed p is itself
    # O(lambda^2) ~ 1/N^2, the same order as its finite-N corrections, so
    # the sharp gate lives on the O(1) words, not the alternating sums)
    fixed_p = []
    sys.path.insert(0, str(HERE.parent / "syk-self-averaging"))
    import exact_pair_moments as epm
    for N in (64, 256, 1024):
        lamN = 2.0 * p * p / N
        row = dict(N=N, lam=lamN, words={})
        for (a, b, c, d) in [(1, 1, 1, 1), (2, 1, 0, 1)]:
            ch = diag_word_amplitude(math.exp(-lamN), dpsi,
                                     word_diag(a, b, c, d))
            fN = epm.t4bar(N, p, a, b, c, d, marks=('i', 'i', 'i', 'i'))
            row["words"][str((a, b, c, d))] = dict(
                chord=ch, finiteN=fN, dev=abs(ch - fN))
            epm._DP_CACHE.clear()
            epm._EXP_CACHE.clear()
        fixed_p.append(row)
        d1111 = row["words"]["(1, 1, 1, 1)"]["dev"]
        print(f"fixed p={p} N={N}: word (1,1,1,1) dev={d1111:.2e} "
              f"(N^2 dev = {N * N * d1111:.0f})", flush=True)

    # G1: the A_1 identity's chord shadow along the line
    a1 = []
    for (pp, N) in [(4, 16), (6, 36), (8, 64)]:
        ch = channel_moments_diag(math.exp(-2.0), 1.0 / pp, 1, ks=(0, 2, 4))
        pr = finite_n_pair(N, pp, 1, ks=(0, 2, 4))
        row = dict(p=pp, N=N)
        for k in (2, 4):
            row[f"resid_mu{k}"] = (abs(ch[f"mu{k}"] + (N - 1) * pr[f"mu{k}"])
                                   / abs(ch[f"mu{k}"]))
        a1.append(row)
        print(f"A1 gate p={pp} N={N}: resid mu2={row['resid_mu2']:.3f} "
              f"mu4={row['resid_mu4']:.3f}", flush=True)

    # G3: the diagonal width at lambda_chord = 2 and its Delta -> 0 limit,
    # next to the FULL-channel context.  The full width is a large
    # diagonal-vs-pair cancellation (N mu_k^diag + N(N-1) mu_k^pair with
    # opposite signs), so the diagonal alone does NOT approximate it; the
    # full-channel N = infinity value quoted here is a Richardson
    # extrapolation in 1/p of the exact finite-N series -- data, not a
    # chord computation.  The missing chord-side piece is the pair
    # channel's 1/N coefficient (pair_moment_laws.py territory).
    full = []
    for (pp, N) in [(4, 16), (6, 36), (8, 64), (10, 100)]:
        pr = finite_n_pair(N, pp, 2, ks=(0, 2))
        dg = finite_n_diag(N, pp, 2, ks=(0, 2))
        full.append([pp, (N * dg["mu2"] + N * (N - 1) * pr["mu2"])
                     / (N * dg["mu0"] + N * (N - 1) * pr["mu0"])])
    rich = [w2 + (w2 - w1) * (1.0 / pb) / (1.0 / pa - 1.0 / pb)
            for (pa, w1), (pb, w2) in zip(full, full[1:])]
    g3 = dict(w2_at_inverse_p=[[pp, channel_moments_diag(
        math.exp(-2.0), 1.0 / pp, 2, ks=(0, 2))["w2"]]
        for pp in (4, 6, 8, 10)],
        delta0=w2_delta0(2.0),
        full_channel_finiteN=full,
        full_channel_richardson=rich,
        full_channel_w2_inf_estimate=rich[-1],
        full_channel_estimate_spread=max(rich) - min(rich))

    payload = dict(
        provenance="diagonal_chord.py -- exact N=infinity flavor-DIAGONAL "
                   "A_n channel moments by the parity-law flux walker "
                   "(derived from the finite-N sign law of "
                   "exact_pair_moments.py; assembly shared with "
                   "chord_moments.py); gates: A_1 identity shadow, "
                   "finite-N convergence on and off the lambda_chord = 2 "
                   "line, and the Delta_psi -> 0 diagonal width",
        p=p, runtime_s=round(time.time() - t0, 1), points=points,
        line_lambda2=line, fixed_p=fixed_p, a1_gate=a1, g3=g3)
    out = Path(out_path) if out_path else HERE / "diagonal_chord.json"
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()
    if args.run:
        run()
    else:
        q = math.exp(-2.0)
        for n in (1, 2):
            print(n, channel_moments_diag(q, 0.25, n, ks=(0, 2, 4)))
        print("w2 Delta->0:", w2_delta0(2.0))
