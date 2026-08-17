#!/usr/bin/env python3
r"""
The exact N = infinity chord engine.
====================================

Double-scaled SYK at N = infinity (lambda = 2p^2/N fixed) is exactly
solvable by chord combinatorics (Berkooz-Isachenkov-Narovlansky-Torrents
1811.02584; conventions from the review Berkooz-Mamroud 2407.09396, both
equation-checked).  This module implements the solution as finite sparse
linear algebra.  No special functions.

  * Transfer matrix [review (2.15)-(2.16)]: T|n> = |n+1> + [n]_q |n-1>.
    Moments <Tr H^k>/<Tr H^2>^{k/2} = <0|T^k|0>; spectrum
    E(theta) = 2 cos(theta)/sqrt(1-q)  [review (2.21)].
  * Matter of dimension Delta [review (2.26)-(2.30)]: an H-chord crossing
    the M-chord weighs q^Delta.  The T = infinity 2-pt function is exactly
    <0| e^{+itT} q^{Delta n} e^{-itT} |0>  [review (2.29)]; the
    theta-integral (2.30) is coded separately as a cross-check.
  * One and two matter chords on the contour (stack construction: closing
    an object crosses everything opened after it that is still open).
    This reproduces the one-particle transfer matrices [review
    (5.18)-(5.20)] and evaluates ANY operator ordering on ANY complex
    contour -- crossed and Euclidean-folded included.  Matter-matter
    crossings weigh q^{Delta_1 Delta_2} times a fermionic sign (-1) for
    single-Majorana matter: the same minus that fold_bench.py derives.

Units: the chord Hamiltonian has <H^2> = 1.  Exact bridge to repo units
(psi^2 = 1, Var J = p!/N^{p-1}):  sigma_H^2 = C(N,p) p!/N^{p-1}, and
t_chord = t_repo * sigma_H.  In double scaling this gives
beta_chord = betaJ * e^{-lambda/8}/sqrt(lambda), with scriptJ = sqrt(2) p.

Truncation: finite chord lattice n <= nmax; q^n suppresses large n.
Campaign JSONs store paired truncations, and the tests assert stability.
Per-equation provenance lives in the function docstrings.  Validation:
test_chord.py (enumerator, closed forms, Krawtchouk microscopy, ED
anchors).
"""

from __future__ import annotations

import argparse
import json
import math
import time
from itertools import permutations
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent


def _round6(xs):
    """Arrays stored in campaign JSONs are rounded to 6 significant digits
    (the verdict scalars are kept at full precision) -- keeps the committed
    artifacts lean."""
    return [float(f"{x:.6g}") for x in xs]


# ---------------------------------------------------------------------------
# q-arithmetic
# ---------------------------------------------------------------------------
def qint(n: int, q: float) -> float:
    """[n]_q = (1-q^n)/(1-q); robust n at q = 1."""
    if q == 1.0:
        return float(n)
    return (1.0 - q ** n) / (1.0 - q)


def qints(nmax: int, q: float) -> np.ndarray:
    """[0]_q .. [nmax]_q as an array."""
    return np.array([qint(n, q) for n in range(nmax + 1)])


def qpoch(a: complex, q: float, n: int | None = None) -> complex:
    """(a; q)_n = prod_{k=0}^{n-1} (1 - a q^k); n = None means infinity
    (terminates when |a q^k| < 1e-18; requires |q| < 1)."""
    out = 1.0 + 0.0j
    k, term = 0, complex(a)
    while True:
        if n is not None and k >= n:
            break
        out *= (1.0 - term)
        term *= q
        k += 1
        if n is None and abs(term) < 1e-18:
            break
        if k > 100000:
            raise RuntimeError("qpoch did not terminate (q too close to 1)")
    return out


def log_qbinom(n: int, k: int, q: float) -> float:
    """log of the q-binomial [n; k]_q = [n]_q!/([k]_q![n-k]_q!) >= 0."""
    s = 0.0
    for j in range(1, n + 1):
        s += math.log(qint(j, q))
    for j in range(1, k + 1):
        s -= math.log(qint(j, q))
    for j in range(1, n - k + 1):
        s -= math.log(qint(j, q))
    return s


# ---------------------------------------------------------------------------
# units bridge (exact finite-N and double-scaled forms)
# ---------------------------------------------------------------------------
def sigma_H(N: int, p: int) -> float:
    """sqrt(<Tr H^2>/dim) in repo units (Var J = p!/N^{p-1}, psi^2 = 1):
    sigma_H^2 = C(N,p) p!/N^{p-1} = N prod_{k=1}^{p-1}(1 - k/N).  Exact."""
    prod = 1.0
    for k in range(1, p):
        prod *= (1.0 - k / N)
    return math.sqrt(N * prod)


def lambda_chord_of(N: int, p: int) -> float:
    """lambda_chord = 2 p^2 / N (the chord-literature convention; see
    dictionary.lambda_chord -- the factor-2 landmine is tested there)."""
    return 2.0 * p * p / N


def beta_chord_from_betaJ(betaJ: float, lam: float) -> float:
    """N = infinity bridge at fixed beta*scriptJ (scriptJ = sqrt(2) p):
    beta_chord = beta_phys sigma_H -> betaJ e^{-lam/8}/sqrt(lam)."""
    return betaJ * math.exp(-lam / 8.0) / math.sqrt(lam)


def q_of_lambda(lam: float) -> float:
    return math.exp(-lam)


# ---------------------------------------------------------------------------
# sector 0: no matter chords.  Symmetric tridiagonal in the scaled basis.
# ---------------------------------------------------------------------------
class Sector0:
    """Chord states |n>, n = 0..nmax, sqrt-scaled basis (T symmetric).

    T |n> = sqrt([n+1]) |n+1> + sqrt([n]) |n-1>   [review (2.15) scaled]
    """

    def __init__(self, nmax: int, q: float):
        self.nmax, self.q = nmax, q
        self.off = np.sqrt(qints(nmax, q)[1:])          # sqrt([1]..[nmax])
        self.E, self.V = eigh_tridiagonal(np.zeros(nmax + 1), self.off)
        self.c0 = self.V[0, :]                          # <n=0|theta_j>

    def vac(self) -> np.ndarray:
        v = np.zeros(self.nmax + 1, dtype=complex)
        v[0] = 1.0
        return v

    def evolve(self, v: np.ndarray, u: complex) -> np.ndarray:
        """e^{-u T} v for complex u (Re u > 0: Euclidean; Im: Lorentzian)."""
        return self.V @ (np.exp(-u * self.E) * (self.V.T @ v))

    def qn(self, delta: float) -> np.ndarray:
        """diagonal of q^{Delta n}."""
        return self.q ** (delta * np.arange(self.nmax + 1))

    def apply_T(self, v: np.ndarray) -> np.ndarray:
        out = np.zeros_like(v)
        out[1:] += self.off * v[:-1]
        out[:-1] += self.off * v[1:]
        return out


# ---------------------------------------------------------------------------
# sectors with matter chords (sqrt-scaled basis; T non-symmetric but bounded)
# ---------------------------------------------------------------------------
class SectorM:
    """One matter chord of dimension `delta` (or the top chord of two --
    the bottom-empty two-matter sector has the identical T with delta of the
    SECOND-opened chord; the first-opened chord is inert until it closes).

    States |n0, n1>: n0 H-chords below the matter chord (opened earlier),
    n1 above.  Scaled-basis action of one H node (derivation in module
    docstring; reduces to review (5.18)-(5.20) one-boundary case, and is
    validated exactly against the brute-force enumerator):

        open (top):            (n0, n1) -> (n0, n1+1)  coeff sqrt([n1+1])
        close from above:      (n0, n1) -> (n0, n1-1)  coeff sqrt([n1])
        close from below:      (n0, n1) -> (n0-1, n1)  coeff q^{Delta+n1} sqrt([n0])
    """

    def __init__(self, nmax: int, q: float, delta: float):
        self.nmax, self.q, self.delta = nmax, q, delta
        m = nmax + 1
        self.dim = m * m
        sq = np.sqrt(qints(nmax, q))
        rows, cols, vals = [], [], []
        for n0 in range(m):
            for n1 in range(m):
                i = n0 * m + n1
                if n1 + 1 <= nmax:                       # open on top
                    rows.append(n0 * m + n1 + 1)
                    cols.append(i)
                    vals.append(sq[n1 + 1])
                if n1 >= 1:                              # close from above
                    rows.append(n0 * m + n1 - 1)
                    cols.append(i)
                    vals.append(sq[n1])
                if n0 >= 1:                              # close from below
                    rows.append((n0 - 1) * m + n1)
                    cols.append(i)
                    vals.append(q ** (delta + n1) * sq[n0])
        self.T = sp.csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim))

    def evolve(self, v: np.ndarray, u: complex) -> np.ndarray:
        if u == 0:
            return v
        return expm_multiply(-u * self.T, v)

    def apply_T(self, v: np.ndarray) -> np.ndarray:
        return self.T @ v


# ---------------------------------------------------------------------------
# the contour walker
# ---------------------------------------------------------------------------
class Contour:
    """Evaluate <Tr[ ... ]>_J / <Tr 1> for an arbitrary ordering of matter
    operators (each label appearing exactly twice) and complex evolution
    segments, walking the trace contour with the chord stack.

    Usage: ops = [(label, u_after), ...] in TRACE order; the walker processes
    them RIGHT to LEFT (rightmost acts first on |0>), i.e. the cut sits just
    before the leftmost symbol.  `deltas[label]` gives the chord dimension;
    labels in `fermionic` are single-Majorana-type (odd), contributing the
    transposition sign (-1) per matter-matter crossing.

    Supported stack discipline (all 2-pt and 4-pt orderings, any fold):
    at most two matter chords open at once, and the second opens only while
    the region below the first is empty (guaranteed when the cut is at the
    first operator of a 4-pt word).
    """

    def __init__(self, q: float, nmax: int, deltas: dict,
                 fermionic: set | None = None):
        self.q, self.nmax = q, nmax
        self.deltas = deltas
        self.fermionic = fermionic if fermionic is not None else set()
        self.s0 = Sector0(nmax, q)
        self._sm = {}

    def sector(self, delta: float) -> SectorM:
        if delta not in self._sm:
            self._sm[delta] = SectorM(self.nmax, self.q, delta)
        return self._sm[delta]

    # -- maps between sectors (scaled-basis coefficients; module docstring) --
    def _open_first(self, v0: np.ndarray) -> np.ndarray:
        """0-matter |n> -> |n0=n, n1=0>, coeff 1."""
        m = self.nmax + 1
        v = np.zeros(m * m, dtype=complex)
        v[np.arange(m) * m] = v0
        return v

    def _close_single(self, v1: np.ndarray, delta: float) -> np.ndarray:
        """|n0, n1> -> q^{Delta n1} |n0+n1> (formal basis); in the scaled
        basis (v_s = sqrt([n0]![n1]!) v_f) the coefficient picks up
        sqrt([n0+n1]!/([n0]![n1]!)) = sqrt(qbinom(n0+n1, n1))."""
        m = self.nmax + 1
        out = np.zeros(m, dtype=complex)
        for n0 in range(m):
            for n1 in range(m):
                tot = n0 + n1
                if tot > self.nmax:
                    continue
                amp = v1[n0 * m + n1]
                if amp == 0:
                    continue
                c = self.q ** (delta * n1) * math.exp(
                    +0.5 * log_qbinom(tot, n1, self.q))
                out[tot] += c * amp
        return out

    def _open_second(self, v1: np.ndarray) -> np.ndarray:
        """1-matter |0, n1> -> 2-matter (m1=n1, m2=0).  Requires the
        below-region to be empty (asserted)."""
        m = self.nmax + 1
        big = np.abs(v1.reshape(m, m)[1:, :]).max() if m > 1 else 0.0
        norm = np.abs(v1).max()
        if norm > 0 and big > 1e-10 * norm:
            raise ValueError("second matter chord opened with occupied "
                             "below-region: unsupported stack shape")
        v = np.zeros(m * m, dtype=complex)
        v[np.arange(m) * m] = v1.reshape(m, m)[0, :]
        return v

    def _close_bottom_of_two(self, v2: np.ndarray, delta_bot: float,
                             delta_top: float, sign: float) -> np.ndarray:
        """2-matter (m1, m2) -> 1-matter (n0=m1, n1=m2), coeff
        sign * q^{D_bot D_top} * qt_bot^{m1+m2}   (crosses everything above:
        m1+m2 H-chords and the top matter chord).  Scaled coeff unchanged."""
        m = self.nmax + 1
        out = np.zeros(m * m, dtype=complex)
        qt = self.q ** delta_bot
        pref = sign * self.q ** (delta_bot * delta_top)
        for m1 in range(m):
            for m2 in range(m):
                amp = v2[m1 * m + m2]
                if amp == 0:
                    continue
                out[m1 * m + m2] += pref * qt ** (m1 + m2) * amp
        return out

    def _close_top_of_two(self, v2: np.ndarray, delta_top: float) -> np.ndarray:
        """2-matter (m1, m2) -> 1-matter (0, m1+m2), scaled coeff
        qt_top^{m2} * sqrt(qbinom(m1+m2, m2))  (regions merge above the
        remaining bottom chord)."""
        m = self.nmax + 1
        out = np.zeros(m * m, dtype=complex)
        for m1 in range(m):
            for m2 in range(m):
                tot = m1 + m2
                if tot > self.nmax:
                    continue
                amp = v2[m1 * m + m2]
                if amp == 0:
                    continue
                c = self.q ** (delta_top * m2) * math.exp(
                    +0.5 * log_qbinom(tot, m2, self.q))
                out[tot] += c * amp
        return out

    # -- the walk --
    def correlator(self, ops: list) -> complex:
        """ops = [(label, u_after), ...] in trace order; returns the
        normalized ensemble-averaged trace (Z-normalization is the caller's
        business -- divide by <0|e^{-beta T}|0> for thermal correlators)."""
        seq = list(reversed(ops))          # rightmost operator acts first
        state = self.s0.vac()
        open_stack = []                    # labels, bottom -> top
        sector = 0
        for k, (label, u_after) in enumerate(seq):
            # segment BEFORE this operator in ket-order is u of previous op;
            # we apply operator first, then its u_after (the segment that
            # sits to its LEFT in the trace word, i.e. later on the walk).
            if label is not None:
                if label not in open_stack:
                    # opening
                    if sector == 0:
                        state = self._open_first(state)
                        sector = 1
                    elif sector == 1:
                        state = self._open_second(state)
                        sector = 2
                    else:
                        raise ValueError(">2 open matter chords unsupported")
                    open_stack.append(label)
                else:
                    # closing
                    if sector == 1:
                        assert open_stack == [label]
                        state = self._close_single(state, self.deltas[label])
                        open_stack.remove(label)
                        sector = 0
                    else:
                        d_bot, d_top = (self.deltas[open_stack[0]],
                                        self.deltas[open_stack[1]])
                        if label == open_stack[0]:      # crossed closing
                            sgn = -1.0 if (
                                label in self.fermionic
                                and open_stack[1] in self.fermionic) else 1.0
                            state = self._close_bottom_of_two(
                                state, d_bot, d_top, sgn)
                            # remaining chord is the old top; its delta
                            # governs the 1-matter sector from here on
                        else:                            # nested closing
                            state = self._close_top_of_two(state, d_top)
                        open_stack.remove(label)
                        sector = 1
            if u_after != 0:
                if sector == 0:
                    state = self.s0.evolve(state, u_after)
                else:
                    delta_now = self.deltas[open_stack[-1]] if sector == 2 \
                        else self.deltas[open_stack[0]]
                    if sector == 1 and len(open_stack) == 1:
                        # general 1-matter evolution: T depends on the open
                        # chord's dimension
                        state = self.sector(delta_now).evolve(state, u_after)
                    else:
                        # two open: bottom chord inert, top chord's delta
                        state = self.sector(
                            self.deltas[open_stack[1]]).evolve(state, u_after)
        assert sector == 0 and not open_stack, "unbalanced operator labels"
        return complex(state[0])

    def apply_word_powers(self, word: list) -> complex:
        """<Tr[ word ]>: word is a list over {'H', labels}; exact integer
        powers of T (no exponentials) -- for moment validation.

        The trace is cyclic, and the walker's stack discipline requires the
        below-region of the first-opened matter chord to stay empty until it
        closes; rotating the word so it ENDS with a matter symbol guarantees
        this (the walk processes right-to-left, so that symbol opens on the
        true vacuum)."""
        word = list(word)
        if any(w != 'H' for w in word):
            while word[-1] == 'H':
                word = word[-1:] + word[:-1]
        seq = list(reversed(word))
        ops, i = [], 0
        # compress into (label, k_H_after) in reversed order, then rebuild
        # as correlator-style processing with T-powers
        state = self.s0.vac()
        open_stack, sector = [], 0
        for w in seq:
            if w == 'H':
                if sector == 0:
                    state = self.s0.apply_T(state)
                else:
                    delta_now = self.deltas[open_stack[1]] if sector == 2 \
                        else self.deltas[open_stack[0]]
                    state = self.sector(delta_now).apply_T(state)
            else:
                label = w
                if label not in open_stack:
                    if sector == 0:
                        state = self._open_first(state)
                        sector = 1
                    elif sector == 1:
                        state = self._open_second(state)
                        sector = 2
                    else:
                        raise ValueError(">2 open matter chords unsupported")
                    open_stack.append(label)
                else:
                    if sector == 1:
                        state = self._close_single(state, self.deltas[label])
                        open_stack.remove(label)
                        sector = 0
                    else:
                        d_bot, d_top = (self.deltas[open_stack[0]],
                                        self.deltas[open_stack[1]])
                        if label == open_stack[0]:
                            sgn = -1.0 if (
                                label in self.fermionic
                                and open_stack[1] in self.fermionic) else 1.0
                            state = self._close_bottom_of_two(
                                state, d_bot, d_top, sgn)
                        else:
                            state = self._close_top_of_two(state, d_top)
                        open_stack.remove(label)
                        sector = 1
        assert sector == 0 and not open_stack, "unbalanced word"
        return complex(state[0])


# ---------------------------------------------------------------------------
# brute-force chord-diagram enumerator (independent ground truth)
# ---------------------------------------------------------------------------
def bruteforce_word(q: float, word: list, deltas: dict,
                    fermionic: set | None = None) -> float:
    """<Tr[ word ]> by explicit summation over chord diagrams on the circle.

    word: list over {'H', matter labels}; every matter label appears exactly
    twice (those two nodes are contracted with each other -- one M-chord),
    and the H nodes are summed over all perfect matchings.  Weights:
    q per H-H crossing, q^{Delta_i} per H-M_i crossing,
    (sign_ij) q^{Delta_i Delta_j} per M_i-M_j crossing, where sign_ij = -1
    iff both labels are fermionic.  Crossing test: chords (a,b), (c,d) cross
    iff exactly one of c,d lies strictly between a and b along the circle.

    This shares NOTHING with the transfer-matrix engine except the
    microscopic Wick weights -- it is the ground truth the engine is
    validated against (test_chord.py)."""
    fermionic = fermionic if fermionic is not None else set()
    n = len(word)
    h_pos = [i for i, w in enumerate(word) if w == 'H']
    labels = sorted({w for w in word if w != 'H'})
    m_chords = []
    for lab in labels:
        pos = [i for i, w in enumerate(word) if w == lab]
        if len(pos) != 2:
            raise ValueError(f"label {lab} must appear exactly twice")
        m_chords.append((pos[0], pos[1], lab))
    if len(h_pos) % 2:
        return 0.0

    def crossings(c1, c2):
        (a, b), (c, d) = sorted(c1[:2]), sorted(c2[:2])
        in1 = (a < c < b)
        in2 = (a < d < b)
        return in1 != in2

    def matchings(items):
        if not items:
            yield []
            return
        a = items[0]
        for k in range(1, len(items)):
            b = items[k]
            rest = items[1:k] + items[k + 1:]
            for rem in matchings(rest):
                yield [(a, b)] + rem

    total = 0.0
    for hm in matchings(h_pos):
        w = 1.0
        # H-H crossings
        for i in range(len(hm)):
            for j in range(i + 1, len(hm)):
                if crossings(hm[i], hm[j]):
                    w *= q
        # H-M crossings
        for (a, b) in hm:
            for (c, d, lab) in m_chords:
                if crossings((a, b), (c, d)):
                    w *= q ** deltas[lab]
        # M-M crossings
        for i in range(len(m_chords)):
            for j in range(i + 1, len(m_chords)):
                ci, cj = m_chords[i], m_chords[j]
                if crossings(ci[:2], cj[:2]):
                    w *= q ** (deltas[ci[2]] * deltas[cj[2]])
                    if ci[2] in fermionic and cj[2] in fermionic:
                        w *= -1.0
        total += w
    return total


# ---------------------------------------------------------------------------
# closed forms (independent q-Pochhammer route)  [review (2.21), (2.30)]
# ---------------------------------------------------------------------------
def rho_theta(theta: np.ndarray, q: float) -> np.ndarray:
    """Density of states rho(theta) = (q; q)_inf |(e^{2 i theta}; q)_inf|^2
    / (2 pi) on [0, pi]  [review (2.21)]; vectorized via
    |(z; q)_inf|^2 = prod_k (1 - 2 q^k cos(2 theta) + q^{2k})."""
    theta = np.atleast_1d(np.asarray(theta, float))
    cq = abs(qpoch(q, q))
    c2 = np.cos(2.0 * theta)
    out = np.ones_like(theta)
    a = 1.0
    while a > 1e-18:
        out *= (1.0 - 2.0 * a * c2 + a * a)
        a *= q
    return cq * out / (2.0 * np.pi)


def energy_theta(theta: np.ndarray, q: float) -> np.ndarray:
    """E(theta) = 2 cos(theta)/sqrt(1-q)  [review (2.21)]."""
    return 2.0 * np.cos(theta) / math.sqrt(1.0 - q)


def vertex2(th1: float, th2: float, q: float, delta: float) -> float:
    """|<theta_1|M_Delta|theta_2>|^2 = (qt^2; q)_inf /
    prod_{s1,s2 = +-} (qt e^{i(s1 th1 + s2 th2)}; q)_inf,  qt = q^Delta
    [review (2.30)/(4.7)]."""
    qt = q ** delta
    num = qpoch(qt * qt, q).real
    den = 1.0 + 0.0j
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            den *= qpoch(qt * np.exp(1j * (s1 * th1 + s2 * th2)), q)
    return float(num / den.real)


def vertex2_matrix(th: np.ndarray, q: float, delta: float) -> np.ndarray:
    """|<th_i|M_Delta|th_j>|^2 on a theta grid, vectorized: the product over
    the four sign choices pairs into real quadratics,
    prod_k [1 - 2 qt q^k cos(th_i+th_j) + qt^2 q^{2k}]
           [1 - 2 qt q^k cos(th_i-th_j) + qt^2 q^{2k}]."""
    qt = q ** delta
    num = qpoch(qt * qt, q).real
    Cp = np.cos(np.add.outer(th, th))
    Cm = np.cos(np.subtract.outer(th, th))
    den = np.ones((len(th), len(th)))
    a = qt
    while a > 1e-18:
        den *= (1.0 - 2.0 * a * Cp + a * a) * (1.0 - 2.0 * a * Cm + a * a)
        a *= q
    return num / den


def g2_closed_form(beta1: complex, beta2: complex, q: float, delta: float,
                   n_theta: int = 200) -> complex:
    """The theta-integral (2.30): unnormalized
    G = int dth1 dth2 rho rho e^{-beta1 E1 - beta2 E2} |<th1|M|th2>|^2.
    Gauss-Legendre quadrature on [0, pi]; independent of the transfer-matrix
    route."""
    x, wq = np.polynomial.legendre.leggauss(n_theta)
    th = 0.5 * np.pi * (x + 1.0)
    w = 0.5 * np.pi * wq
    rho = rho_theta(th, q)
    E = energy_theta(th, q)
    V = vertex2_matrix(th, q, delta)
    f1 = w * rho * np.exp(-np.asarray(beta1) * E)
    f2 = w * rho * np.exp(-np.asarray(beta2) * E)
    return complex(f1 @ V @ f2)


# ---------------------------------------------------------------------------
# two-point function, spectral lines, and the Delta -> 0 (length) channel
# ---------------------------------------------------------------------------
def g2_engine(ts: np.ndarray, q: float, delta: float, nmax: int,
              beta: float = 0.0) -> np.ndarray:
    """G_Delta(t) at inverse temperature beta (chord units):
    <0| e^{-(beta/2) T} e^{+itT} q^{Delta n} e^{-itT} e^{-(beta/2) T} |0>
    / <0|e^{-beta T}|0>   (the symmetric thermal ordering; at beta = 0 this
    is exactly review (2.29) continued to real time)."""
    s = Sector0(nmax, q)
    qn = s.qn(delta)
    c = s.V.T @ s.vac().real
    if beta != 0.0:
        c = c * np.exp(-0.5 * beta * s.E)
    M = s.V.T @ (qn[:, None] * s.V)        # <theta_j| q^{Dn} |theta_k>
    Z = float(c @ c)
    out = np.empty(len(ts), dtype=complex)
    for i, t in enumerate(np.atleast_1d(ts)):
        ph = np.exp(-1j * t * s.E)
        out[i] = (c * ph.conj()) @ M @ (ph * c) / Z
    return out


def spectral_lines(q: float, delta: float, nmax: int, beta: float = 0.0):
    """Lehmann lines of S_Delta(omega): returns (omega_jk, w_jk) with
    om_jk = E_j - E_k and G(t) = sum w e^{+i om t} (the resummation used in
    the tests; for the symmetric weights here S(om) = S(-om), so the
    committed spectra are convention-free).  The truncated-T eigenbasis is
    the exact Gaussian quadrature of the theta-integrals."""
    s = Sector0(nmax, q)
    c = s.V.T @ s.vac().real
    if beta != 0.0:
        c = c * np.exp(-0.5 * beta * s.E)
    M = s.V.T @ (s.qn(delta)[:, None] * s.V)
    W = np.outer(c, c) * M
    Z = float(c @ c)
    om = np.subtract.outer(s.E, s.E)
    return om.ravel(), (W / Z).ravel()


def n_expectation(ts: np.ndarray, q: float, nmax: int) -> np.ndarray:
    """<n(t)> = <0| e^{+itT} n e^{-itT} |0>: the chord-number (wormhole
    length / size) growth curve.  The Delta -> 0 bilinear channel is
    EXACTLY  d G_Delta/d Delta|_0 = -lambda <n(t)>  (from (2.29):
    d q^{Delta n}/d Delta = -lambda n q^{Delta n})."""
    s = Sector0(nmax, q)
    c = s.V.T @ s.vac().real
    Nmat = s.V.T @ (np.arange(nmax + 1)[:, None] * s.V)
    out = np.empty(len(ts))
    for i, t in enumerate(np.atleast_1d(ts)):
        ph = np.exp(-1j * t * s.E)
        out[i] = float(np.real((c * ph.conj()) @ Nmat @ (ph * c)))
    return out


def length_channel_lines(q: float, nmax: int):
    """Lehmann lines of the Delta->0 channel: -lambda * <0|n(t)|0>-lines,
    i.e. S_L(omega) with  dG/dDelta|_0(t) = sum w e^{i(E_j-E_k)t}."""
    lam = -math.log(q)
    s = Sector0(nmax, q)
    c = s.V.T @ s.vac().real
    Nmat = s.V.T @ (np.arange(nmax + 1)[:, None] * s.V)
    W = -lam * np.outer(c, c) * Nmat
    om = np.subtract.outer(s.E, s.E)
    return om.ravel(), W.ravel()


# ---------------------------------------------------------------------------
# binned spectral function + structure diagnostics (rung-18 deliverable)
# ---------------------------------------------------------------------------
def binned_spectrum(om: np.ndarray, w: np.ndarray, bins: np.ndarray):
    """Histogram Lehmann lines into S(omega) on the given bin edges
    (weights summed, divided by bin width)."""
    idx = np.digitize(om, bins) - 1
    S = np.zeros(len(bins) - 1)
    ok = (idx >= 0) & (idx < len(S))
    np.add.at(S, idx[ok], np.real(w[ok]))
    return S / np.diff(bins)


def _vertex2_pairs(th1: np.ndarray, th2: np.ndarray, q: float,
                   delta: float) -> np.ndarray:
    """<th1|q^{Delta n}|th2> elementwise on paired theta arrays (the (2.30)
    kernel; real-factorized q-Pochhammer products)."""
    qt = q ** delta
    num = qpoch(qt * qt, q).real
    Cp = np.cos(th1 + th2)
    Cm = np.cos(th1 - th2)
    den = np.ones_like(Cp)
    a = qt
    while a > 1e-18:
        den *= (1.0 - 2.0 * a * Cp + a * a) * (1.0 - 2.0 * a * Cm + a * a)
        a *= q
    return num / den


def s_omega_smooth(omegas: np.ndarray, q: float, delta: float | None,
                   n_theta: int = 2000, eps_delta: float = 1e-5):
    """The EXACT continuous N = infinity spectral function at T_B = infinity:

        S(omega) = int dth1 rho(th1) rho(th2*) K(th1, th2*) / |dE/dth(th2*)|,
        th2* = arccos((E(th1) - omega)/band),   band = 2/sqrt(1-q),

    the reduced form of the double theta-integral with the energy delta
    function eaten analytically.  K is the (2.30) vertex for the O_Delta
    channel (delta = Delta), or its Delta-derivative at 0 for the fermion-
    bilinear/length channel (delta = None): dG/dDelta|_0 spectral density
    (finite for omega != 0; the omega = 0 conserved piece is the linear
    length growth, reported separately as <n(t)>).

    Normalization checks (tested): int S domega = 1 and
    int omega^2 S = 2(1-q^Delta) for the O_Delta channel."""
    band = 2.0 / math.sqrt(1.0 - q)
    x, wq = np.polynomial.legendre.leggauss(n_theta)
    th1 = 0.5 * np.pi * (x + 1.0)
    w1 = 0.5 * np.pi * wq
    rho1 = rho_theta(th1, q)
    E1 = energy_theta(th1, q)
    out = np.empty(len(omegas))
    for i, om in enumerate(np.atleast_1d(omegas)):
        E2 = E1 - om
        ok = np.abs(E2) < band * (1.0 - 1e-12)
        th2 = np.arccos(np.clip(E2[ok] / band, -1.0, 1.0))
        jac = band * np.sin(th2)
        if delta is None:
            # off the diagonal K_{Delta=0} = 0 (the identity channel is pure
            # theta_1 = theta_2 concentration), so dK/dDelta|_0 by two-point
            # Richardson is 2 K(eps)/eps - K(2 eps)/(2 eps).
            Kp = _vertex2_pairs(th1[ok], th2, q, eps_delta)
            Km = _vertex2_pairs(th1[ok], th2, q, 2 * eps_delta)
            K = 2.0 * Kp / eps_delta - Km / (2 * eps_delta)
        else:
            K = _vertex2_pairs(th1[ok], th2, q, delta)
        out[i] = float(np.sum(w1[ok] * rho1[ok] * rho_theta(th2, q) * K / jac))
    return out


def structure_report_smooth(x: np.ndarray, S: np.ndarray):
    """Peak/node diagnostics of a SMOOTH spectral curve on a grid: interior
    local maxima (features), the largest secondary-peak height relative to
    the global maximum, sign changes (nodes), and the location of the
    maximum.  On the exact continuum function these are genuine features,
    not binning artifacts (quadrature convergence asserted separately)."""
    S = np.asarray(S, float)
    imax = [i for i in range(1, len(S) - 1)
            if S[i] > S[i - 1] and S[i] >= S[i + 1]]
    peak0 = float(np.abs(S).max()) if len(S) else 0.0
    secondary = 0.0
    if len(imax) > 1 and peak0 > 0:
        heights = sorted((float(S[i]) for i in imax), reverse=True)
        secondary = heights[1] / peak0
    nodes = int(np.sum(np.abs(np.diff(np.sign(S[np.abs(S) > 0]))) > 1))
    return dict(n_local_maxima=len(imax),
                secondary_peak_fraction=secondary,
                n_sign_changes=nodes,
                omega_at_max=float(x[int(np.argmax(np.abs(S)))])
                if len(S) else 0.0)


# ---------------------------------------------------------------------------
# campaigns
# ---------------------------------------------------------------------------
def sech_spectral(omegas: np.ndarray, Jc: float, delta: float) -> np.ndarray:
    """Numeric Fourier transform of the semiclassical (large-p) form
    sech(Jc t)^{2 Delta}: the smooth continuum every channel approaches as
    lambda -> 0 IF no structure emerges.  S_sc(om) = (1/2pi) int dt
    e^{-i om t} sech^{2D}; even, positive, exponentially convergent."""
    T = 40.0 / Jc / max(2 * delta, 0.2)
    t = np.linspace(0.0, T, 40000)
    f = (1.0 / np.cosh(Jc * t)) ** (2 * delta)
    out = np.empty(len(omegas))
    for i, om in enumerate(np.atleast_1d(omegas)):
        out[i] = float(np.trapezoid(f * np.cos(om * t), t) / np.pi)
    return out


def run_spectra(lams=(2.0, 1.0, 0.5, 0.25, 0.1, 0.05),
                deltas=(1.0, 0.5, 0.25), out_path=None) -> dict:
    """Rung-18 first deliverable: the lambda -> 0 scan of the EXACT smooth
    N = infinity, T_B = infinity spectral functions (reduced theta-integral
    -- a continuous function, no binning artifacts), O_Delta channels and
    the Delta -> 0 fermion-bilinear/length channel.

    Two omega windows per point: a FINE window [0, 8 Jc] at the
    semiclassical scale (where a 't Hooft tower would sit under the flat
    dictionary: omega_chord = (omega/scriptJ) * Jc) and a COARSE window
    [0, 1.2 band] for the global shape.  Diagnostics: quadrature-doubling
    drift, moment closures (int S = 1, mu2 = 2(1-q^Delta), length-channel
    int omega^2 S_L = 2 lambda), structure counts on the fine window, and
    the deviation from the semiclassical sech^{2Delta} transform."""
    t0 = time.time()
    points = []
    for lam in lams:
        q = q_of_lambda(lam)
        band = 2.0 / math.sqrt(1.0 - q)
        Jc = math.sqrt(lam) * math.exp(lam / 8.0)
        n_th = int(min(4000, max(800, 60.0 / lam)))
        om_fine = np.linspace(0.0, min(8.0 * Jc, 2 * band), 320)
        # coarse grid covers the full support 2*band, with geometric small-
        # omega resolution (the length channel diverges ~1/om^2 there and
        # the moment integrands need it resolved)
        om_coarse = np.unique(np.concatenate([
            np.geomspace(1e-4 * Jc, Jc, 120),
            np.linspace(Jc, 2.05 * band, 320)]))
        row = dict(lam=lam, q=q, band=band, Jc_bridge=Jc, n_theta=n_th,
                   om_fine=[0.0, float(om_fine[-1]), len(om_fine)],
                   om_coarse=_round6(om_coarse),
                   channels={})
        for d in list(deltas) + [None]:
            key = "length" if d is None else f"delta={d}"
            om_f = om_fine if d is not None else om_fine[om_fine > 0]
            Sf = s_omega_smooth(om_f, q, d, n_th)
            Sf2 = s_omega_smooth(om_f, q, d, 2 * n_th)
            Sc = s_omega_smooth(om_coarse, q, d, n_th)
            drift = float(np.max(np.abs(Sf - Sf2)) /
                          max(np.max(np.abs(Sf2)), 1e-300))
            # moment closures on the symmetric extension (omega > 0 doubled)
            mu0 = float(2.0 * np.trapezoid(Sc, om_coarse))
            mu2 = float(2.0 * np.trapezoid(Sc * om_coarse ** 2, om_coarse))
            if d is None:
                mu2_closed = 2.0 * lam       # d mu2 / d Delta at 0
                mu0_closed = None            # omega->0 growth divergence
                sel = om_f > 0.5 * Jc        # exclude the 1/om^2 growth peak
                rep = structure_report_smooth(om_f[sel], Sf2[sel])
                ssc = np.array([])
                sech_dev = None
            else:
                mu2_closed = 2.0 * (1.0 - q ** d)
                mu0_closed = 1.0
                rep = structure_report_smooth(om_f, Sf2)
                ssc = sech_spectral(om_f, Jc, d)
                sech_dev = float(np.max(np.abs(Sf2 - ssc)) /
                                 max(np.max(ssc), 1e-300))
            row['channels'][key] = dict(
                om_fine=_round6(om_f),
                S_fine=_round6(Sf2), S_coarse=_round6(Sc),
                S_semiclassical_fine=_round6(ssc),
                quadrature_drift=drift, mu0=mu0, mu0_closed=mu0_closed,
                mu2=mu2, mu2_closed=mu2_closed, sech_dev=sech_dev, **rep)
        # length growth curve vs the semiclassical closed form (transfer
        # matrix route; truncation-clean window)
        nmax = int(max(80, min(3000, 30.0 / lam + 60)))
        ts = np.linspace(0.0, min(10.0 / Jc, 0.45 * nmax / band), 60)
        nt = n_expectation(ts, q, nmax)
        nt_sc = -(2.0 / lam) * np.log(1.0 / np.cosh(Jc * ts))
        row['length_growth'] = dict(ts=_round6(ts), nmax=nmax,
                                    n_of_t=_round6(nt),
                                    n_of_t_semiclassical=_round6(nt_sc))
        points.append(row)
        ch = row['channels']
        print(f"lam={lam}: n_theta={n_th} "
              f"maxima(D=0.5)={ch['delta=0.5']['n_local_maxima']} "
              f"maxima(length)={ch['length']['n_local_maxima']} "
              f"sech_dev(D=0.5)={ch['delta=0.5']['sech_dev']:.3f} "
              f"drift={ch['delta=0.5']['quadrature_drift']:.1e}", flush=True)
    payload = dict(
        provenance="chord.py run_spectra -- rung-18 first deliverable: "
                   "exact smooth N=infinity T_B=infinity spectral functions "
                   "(reduced theta-integral of review 2407.09396 (2.30); "
                   "length channel = dG/dDelta|_0), fine window at the "
                   "semiclassical scale 8*Jc where a 't Hooft tower would "
                   "sit, structure diagnostics on the doubled-quadrature "
                   "curve",
        runtime_s=round(time.time() - t0, 1), points=points)
    out = Path(out_path) if out_path else HERE / "chord_spectra.json"
    with open(out, "w") as fh:
        json.dump(payload, fh)
    print(f"wrote {out} ({payload['runtime_s']} s)")
    return payload


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--spectra", action="store_true",
                    help="run the rung-18 lambda->0 spectral scan")
    ap.add_argument("--demo", action="store_true",
                    help="quick sanity demo (moments + 2pt cross-checks)")
    args = ap.parse_args()
    if args.spectra:
        run_spectra()
    else:
        q = 0.5
        c = Contour(q, 24, {'A': 0.5})
        word = ['H', 'A', 'H', 'H', 'A', 'H']
        bf = bruteforce_word(q, word, {'A': 0.5})
        en = c.apply_word_powers(word)
        print(f"word {word}: brute force = {bf:.12f}, engine = {en.real:.12f}")
        g_eng = g2_engine(np.array([0.0]), q, 0.5, 40, beta=1.0)[0]
        Zc = Sector0(60, q)
        z = float(Zc.V[0] @ (np.exp(-1.0 * Zc.E) * Zc.V[0]))
        g_cf = g2_closed_form(0.5, 0.5, q, 0.5) / z
        print(f"thermal 2pt: engine = {g_eng.real:.10f}, "
              f"closed form (2.30) = {g_cf.real:.10f}")
