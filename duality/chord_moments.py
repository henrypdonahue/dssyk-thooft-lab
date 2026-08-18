#!/usr/bin/env python3
r"""
Exact N = infinity singlet-channel spectral moments (rung 15 proper).
=====================================================================

The DSSYK side of the moments-vs-sum-rules contact: EXACT N = infinity,
T_B = infinity spectral moments of the Majorana singlet bilinear channels

    A_n = sum_i psi_i (ad_H)^n psi_i ,        n = 1 .. 5,

computed from the chord expansion with two single-fermion matter chords
(Delta_psi = 1/p, weight qt = q^{1/p}).  Everything reduces to finite exact
combinations of mixed chord moments

    T4(a,b,c,d) = < Tr[ H^a psi_i H^b psi_i H^c psi_j H^d psi_j ] >,
    T2(a,b)     = < Tr[ H^a psi_i H^b psi_i ] >          (i != j, per pair)

evaluated by the validated transfer-matrix engine (chord.py) with EXACT
integer powers -- truncation plays no role (a word with K H-nodes never
occupies more than K/2 + 2 chords).

Assembly (module test = the operator identity A_1 = -2p H):

    ad^n psi = sum_r (-1)^r C(n,r) H^{n-r} psi H^r,
    X = psi ad^n psi ,   Y = X^dag = (-1)^n (ad^n psi) psi ,
    mu_k^sym = (1/2)[ Tr(ad^k(X) Y) + Tr(ad^k(Y) X) ]   (k even),

expanded into T4 words with binomial coefficients; mu_0 subtracts the
identity component <X><Y> (T2 words).  Shape invariants (normalization- and
flavor-count-free): w2 = mu_2/mu_0 (chord units have <H^2> = 1) and
kurt = mu_4 mu_0 / mu_2^2.

Error budget: the chord cancellations are EXACT algebraic identities of
the chord amplitudes at any finite matter weight Delta_psi -- there is no
finite-p matter-weight floor (asserted across a Delta_psi scan in
test_chord_moments.py), and they are now PROVEN: the chord representation
of every A_n insertion is a polynomial in the transfer matrix T, so every
moment mu_k (k >= 1) vanishes at every n (chord_charges.py -- the
dressed-propagator recursion; discovered numerically here first).  The
reported a1_residual is therefore the FLOAT-CANCELLATION yardstick
|mu_2|/gross (the residual of the alternating sum relative to its
pre-cancellation magnitude, ~1e-16), not a physical finite-p effect; it
calibrates how much of the measured 'zero' could be rounding noise.

The lambda -> 0 asymptote of the O_Delta channel kurtosis is the
semiclassical value 3 + 1/Delta (Fourier transform of sech^{2Delta}); the
exact closed forms

    mu_2(Delta) = 2 (1 - qt),
    mu_4(Delta) = 2(2+q) + 6 - 8(2+q) qt + 6(1+q) qt^2 ,   qt = q^Delta

are asserted against the engine in test_chord_moments.py.

Output: chord_moments.json (the comparison against the ED bench and the
Dirac moments.json trend lives in the tests + CHORD.md).
"""

from __future__ import annotations

import argparse
import json
import time
from functools import lru_cache
from math import comb
from pathlib import Path

from chord import Contour, q_of_lambda

HERE = Path(__file__).resolve().parent


class MomentTable:
    """Memoized mixed chord moments at fixed (q, delta_psi)."""

    def __init__(self, q: float, delta_psi: float, nmax: int = 14):
        self.q, self.d = q, delta_psi
        self.c = Contour(q, nmax, {'A': delta_psi, 'B': delta_psi},
                         fermionic={'A', 'B'})

    @lru_cache(maxsize=None)
    def t4(self, a: int, b: int, c: int, d: int) -> float:
        word = (['H'] * a + ['A'] + ['H'] * b + ['A']
                + ['H'] * c + ['B'] + ['H'] * d + ['B'])
        return self.c.apply_word_powers(word).real

    @lru_cache(maxsize=None)
    def t2(self, a: int, b: int) -> float:
        word = ['H'] * a + ['A'] + ['H'] * b + ['A']
        return self.c.apply_word_powers(word).real


def _mu_ordered(tab: MomentTable, n: int, k: int):
    """Tr(ad^k(X) Y)/pair with X = psi ad^n psi, Y = X^dag: word
    H^{k-s} A H^{n-r} A H^{r+s+n-r'} B H^{r'} B with coefficient
    (-1)^{n+r+r'+s} C(k,s) C(n,r) C(n,r').  Returns (value, gross) where
    gross = sum of |term| -- the cancellation yardstick (the amplitudes grow
    with [n]_q powers as q -> 1, so 'zero' must be judged against it)."""
    tot, gross = 0.0, 0.0
    for s in range(k + 1):
        for r in range(n + 1):
            for rp in range(n + 1):
                coeff = ((-1) ** (n + r + rp + s)
                         * comb(k, s) * comb(n, r) * comb(n, rp))
                term = coeff * tab.t4(k - s, n - r, r + s + n - rp, rp)
                tot += term
                gross += abs(term)
    return tot, gross


def _mu_ordered_swapped(tab: MomentTable, n: int, k: int):
    """Tr(ad^k(Y) X): word H^{k-s} A H^{r} A H^{n-r+s+r'} B H^{n-r'} B from
    Y = (-1)^n sum_r (-1)^r C(n,r) H^{n-r} psi H^r psi and
    X = sum_{r'} (-1)^{r'} C(n,r') psi H^{n-r'} psi H^{r'}: assembling
    Tr[H^{k-s} (H^{n-r} psi H^r psi) H^s (psi H^{n-r'} psi H^{r'})] and
    rotating the leading H-block around the trace.  Returns (value, gross)."""
    tot, gross = 0.0, 0.0
    for s in range(k + 1):
        for r in range(n + 1):
            for rp in range(n + 1):
                coeff = ((-1) ** (n + r + rp + s)
                         * comb(k, s) * comb(n, r) * comb(n, rp))
                # word: H^{k-s+n-r} A H^{r} A H^{s} B H^{n-r'} B H^{r'}
                # (trace-rotate the trailing H^{r'} to the front)
                term = coeff * tab.t4(rp + k - s + n - r, r, s, n - rp)
                tot += term
                gross += abs(term)
    return tot, gross


def _trace_part(tab: MomentTable, n: int) -> float:
    """<X> per flavor: sum_r (-1)^r C(n,r) T2(n-r, r) -- the identity
    component subtracted from mu_0."""
    return sum((-1) ** r * comb(n, r) * tab.t2(n - r, r)
               for r in range(n + 1))


def channel_moments(q: float, delta_psi: float, n: int) -> dict:
    """Symmetrized traceless mu_0, mu_2, mu_4, mu_6 of the A_n channel,
    per flavor pair, chord units."""
    tab = MomentTable(q, delta_psi)
    mu, gross = {}, {}
    for k in (0, 2, 4, 6):
        v1, g1 = _mu_ordered(tab, n, k)
        v2, g2 = _mu_ordered_swapped(tab, n, k)
        mu[k] = 0.5 * (v1 + v2)
        gross[k] = 0.5 * (g1 + g2)
    trX = _trace_part(tab, n)
    # <Y> = (-1)^n <reversed X-word> = <X> under the trace (cyclic + the
    # (-1)^n from dagger and the (-1)^r resummation give the same value)
    mu[0] = mu[0] - trX * trX
    # 'zero' judged against the pre-cancellation gross magnitude (float
    # noise scales with it as q -> 1)
    conserved = abs(mu[2]) < 1e-11 * max(gross[2], 1.0)
    return dict(mu0=mu[0], mu2=mu[2], mu4=mu[4], mu6=mu[6],
                gross2=gross[2], gross4=gross[4], gross6=gross[6],
                cancellation2=abs(mu[2]) / max(gross[2], 1e-300),
                cancellation6=abs(mu[6]) / max(gross[6], 1e-300),
                conserved=bool(conserved),
                w2=(mu[2] / mu[0]) if mu[0] != 0 else float('nan'),
                kurt=None if conserved
                else (mu[4] * mu[0] / mu[2] ** 2 if mu[2] != 0
                      else float('nan')))


def odelta_moments(q: float, delta: float) -> dict:
    """Closed forms for the random-operator O_Delta channel (checked
    against the engine in the tests), qt = q^Delta:

        mu_2 = 2 (1 - qt)
        mu_4 = 2(2+q) + 6 - 8(2+q) qt + 6(1+q) qt^2
        mu_6 = (70 + 42q + 6q^2 + 2q^3)
               - 4 (35 + 38q + 14q^2 + 3q^3) qt
               + 30 (1+q)(3 + 2q + q^2) qt^2
               - 20 (1+q)(1 + q + q^2) qt^3

    Structural pins: mu_k(qt = 1) = 0 identically (the identity operator),
    and mu_k(qt -> 0) = sum_s (-1)^s C(k,s) m_{k-s} m_s with m_k the
    Touchard-Riordan chord moments (the fully blocking matter chord
    factorizes the correlator through the vacuum)."""
    qt = q ** delta
    mu2 = 2.0 * (1.0 - qt)
    mu4 = 2.0 * (2.0 + q) + 6.0 - 8.0 * (2.0 + q) * qt + 6.0 * (1.0 + q) * qt * qt
    mu6 = ((70.0 + 42.0 * q + 6.0 * q * q + 2.0 * q ** 3)
           - 4.0 * (35.0 + 38.0 * q + 14.0 * q * q + 3.0 * q ** 3) * qt
           + 30.0 * (1.0 + q) * (3.0 + 2.0 * q + q * q) * qt * qt
           - 20.0 * (1.0 + q) * (1.0 + q + q * q) * qt ** 3)
    return dict(mu0=1.0, mu2=mu2, mu4=mu4, mu6=mu6, w2=mu2,
                kurt=(mu4 / mu2 ** 2) if mu2 != 0 else float('nan'),
                kurt_semiclassical=3.0 + 1.0 / delta)


def run(p: int = 4,
        lams=(4.0, 3.2, 2.6667, 2.2857, 2.0, 1.7778, 1.6, 1.0, 0.5, 0.25,
              0.1),
        out_path=None) -> dict:
    """The rung-15 campaign: exact N = infinity channel moments along the
    lambda grid (the first seven values are lambda_chord = 32/N for the ED
    bench N = 16..20 plus 8..14 anchors; the tail probes the lambda -> 0
    asymptote).  Runtime: seconds."""
    t0 = time.time()
    dpsi = 1.0 / p
    points = []
    for lam in lams:
        q = q_of_lambda(lam)
        row = dict(lam=lam, q=q, p=p, delta_psi=dpsi, channels={})
        for n in (1, 2, 3, 4, 5):
            row['channels'][f"n={n}"] = channel_moments(q, dpsi, n)
        row['channels']['O_delta=2/p'] = odelta_moments(q, 2.0 / p)
        row['channels']['O_delta=1'] = odelta_moments(q, 1.0)
        # the n=1 conservation residual (relative to the pre-cancellation
        # gross magnitude -- the float-noise yardstick)
        row['a1_residual'] = row['channels']['n=1']['cancellation2']
        points.append(row)
        c2, c5 = row['channels']['n=2'], row['channels']['n=5']
        print(f"lam={lam}: n=2 mu2/mu0={c2['w2']:.2e} "
              f"conserved={c2['conserved']}  n=5 canc6={c5['cancellation6']:.1e}"
              f"  a1_residual={row['a1_residual']:.2e}", flush=True)
    # the matter-weight scan: the theorem's 'any Delta_psi' clause, on disk
    dpsi_scan = []
    for dp in (0.1, 0.4, 0.7, 1.0):
        q = q_of_lambda(1.0)
        for n in (1, 2, 3, 4, 5):
            ch = channel_moments(q, dp, n)
            dpsi_scan.append(dict(
                lam=1.0, delta_psi=dp, n=n,
                cancellation2=ch['cancellation2'],
                cancellation4=abs(ch['mu4']) / max(ch['gross4'], 1e-300),
                cancellation6=ch['cancellation6'],
                conserved=ch['conserved']))
    payload = dict(
        provenance="chord_moments.py -- rung 15 proper: exact N=infinity "
                   "T_B=infinity singlet-channel spectral moments from the "
                   "chord expansion (Delta_psi = 1/p matter chords, exact "
                   "integer-power transfer-matrix words); the zeros are "
                   "proven in chord_charges.py (conserved-channel theorem); "
                   "ED comparison: syk-self-averaging/mn_majorana_bench.json",
        p=p, runtime_s=round(time.time() - t0, 1), points=points,
        dpsi_scan=dpsi_scan)
    out = Path(out_path) if out_path else HERE / "chord_moments.json"
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
        # demo: the A_1 identity residual and the n=2 invariants at one point
        q = q_of_lambda(2.0)
        for n in (1, 2):
            print(n, channel_moments(q, 0.25, n))
