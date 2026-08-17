"""Tests for the Majorana singlet-channel moments bench (rung 15 ED side).

The two exact identities that anchor the whole rung-15 story:
  * A_1 = -2p H identically (so the n = 1 channel is conserved at every N);
  * the binomial word-expansion used by duality/chord_moments.py is a
    per-instance OPERATOR identity (validated here against dense algebra,
    independent of any chord input).
Plus the committed-JSON assertions: the measured N-falloff of the n >= 2
channel moments toward the chord prediction w2 = 0 at N = infinity.
"""

import json
from math import comb
from pathlib import Path

import numpy as np
import pytest

from pauli_strings import majorana_hamiltonian_fast
from syk import majorana_operators

HERE = Path(__file__).resolve().parent


def _instance(N, p, seed):
    rng = np.random.default_rng(seed)
    H = majorana_hamiltonian_fast(N, p, rng)
    psis = [op.toarray().astype(complex) for op in majorana_operators(N)]
    return H, psis


def test_a1_identity():
    """A_1 = sum_i psi_i [H, psi_i] = -2p H exactly (any instance)."""
    for (N, p, seed) in [(8, 4, 0), (10, 4, 1), (12, 6, 2)]:
        H, psis = _instance(N, p, seed)
        A1 = sum(P @ (H @ P - P @ H) for P in psis)
        assert np.linalg.norm(A1 + 2 * p * H) < 1e-11 * np.linalg.norm(H)


@pytest.mark.parametrize("n,k", [(1, 2), (2, 2), (2, 4), (3, 2)])
def test_word_expansion_is_operator_identity(n, k):
    """Tr[ad_H^k(psi_i ad^n psi_i) (psi_j ad^n psi_j)^dag] equals the
    binomial-expanded sum of Tr[H^a psi_i H^b psi_i H^c psi_j H^d psi_j]
    words with coefficients (-1)^{n+r+r'+s} C(k,s) C(n,r) C(n,r') --
    the assembly used verbatim by duality/chord_moments.py."""
    N, p = 8, 4
    H, psis = _instance(N, p, 42)
    dim = H.shape[0]
    Hp = {m: np.linalg.matrix_power(H, m) for m in range(k + 2 * n + 1)}
    Pi, Pj = psis[0], psis[3]
    D = Pi.copy()
    for _ in range(n):
        D = H @ D - D @ H
    X = Pi @ D
    Dj = Pj.copy()
    for _ in range(n):
        Dj = H @ Dj - Dj @ H
    Y = (Pj @ Dj).conj().T
    adX = X.copy()
    for _ in range(k):
        adX = H @ adX - adX @ H
    direct = np.trace(adX @ Y) / dim
    tot = 0.0 + 0.0j
    for s in range(k + 1):
        for r in range(n + 1):
            for rp in range(n + 1):
                coeff = ((-1) ** (n + r + rp + s)
                         * comb(k, s) * comb(n, r) * comb(n, rp))
                w = (Hp[k - s] @ Pi @ Hp[n - r] @ Pi
                     @ Hp[r + s + n - rp] @ Pj @ Hp[rp] @ Pj)
                tot += coeff * np.trace(w) / dim
    assert abs(direct - tot) < 1e-10 * max(abs(direct), 1.0)


# ---------------------------------------------------------------------------
# committed-JSON assertions
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def bench():
    with open(HERE / "mn_majorana_bench.json") as fh:
        return json.load(fh)


def test_bench_n1_conserved(bench):
    """The n = 1 rows must show mu_2 = mu_4 = 0 at machine precision
    relative to mu_0 (the A_1 = -2pH conservation, per instance)."""
    for r in bench["points"]:
        if r["n"] == 1:
            assert abs(r["mu2"]) < 1e-10 * abs(r["mu0"])
            assert abs(r["mu4"]) < 1e-8 * abs(r["mu0"])


def test_bench_84_conservation_is_special_not_p_half_N(bench):
    """(N, p) = (8, 4): the n >= 2 channel moments vanish IDENTICALLY --
    an exact per-instance conservation [H, A_n] = 0 special to that point.
    Established by the adversarial review (2026-08-16): this is NOT a
    general p = N/2 phenomenon -- asserted live here: at (12, 6) the
    commutator is O(0.3), not zero.  Mechanism at (8, 4) unknown; the row
    is kept as the exhibition and excluded from trend statements."""
    for r in bench["points"]:
        if r["N"] == 8 and r["p"] == 4 and r["n"] >= 2:
            assert r["w2"] == 0.0
    # live: exact at (8,4), absent at the next p = N/2 point (12,6)
    for (N, p, expect_zero) in [(8, 4, True), (12, 6, False)]:
        H, psis = _instance(N, p, 11)
        A2 = sum(P @ (H @ (H @ P - P @ H) - (H @ P - P @ H) @ H)
                 for P in psis)
        r = np.linalg.norm(H @ A2 - A2 @ H) / max(np.linalg.norm(A2), 1e-300)
        if expect_zero:
            assert r < 1e-12
        else:
            assert r > 0.1


def test_bench_raw_w2_rises_toward_thooft_corner(bench):
    """Along fixed p (lambda = 2p^2/N falling), the RAW singlet-channel
    width w2 RISES: the tower channel turns on toward the fixed-p,
    lambda -> 0 corner where the 't Hooft match must live -- exactly where
    the strict double-scaled (chord) description says the per-pair
    amplitude is conserved.  (N = 2p excluded: the special (8,4)
    identity row.)"""
    for p in (4, 6):
        for n in (2, 3):
            rows = sorted((r for r in bench["points"]
                           if r["n"] == n and r["p"] == p and r["N"] > 2 * p),
                          key=lambda r: r["N"])
            w2 = [r["w2"] for r in rows]
            assert all(w2[i + 1] > w2[i] for i in range(len(w2) - 1)), \
                f"p={p} n={n}: {w2}"


def test_bench_matched_lambda_honest_negative(bench):
    """The honest fixed-lambda statement, post-review: at the near-matched
    pair (p=4, N=10, lam=3.20) <-> (p=6, N=22, lam=3.27) the RAW w2
    RISES with N (0.30 -> 0.51 for n = 2).  The approach to the chord
    N = infinity zero is NOT yet visible at dense-ED sizes -- expected,
    since the chord statement is about the i != j PAIR amplitude while
    the measured channel is flavor-diagonal-dominated at these N (the
    pair part of mu_0 is a few percent and opposite in sign; see
    mn_majorana_bench.py docstring).  A previous (1-2p/N)^2 'artifact
    correction' that manufactured a falling trend was removed by the
    adversarial review: its motivating 'p = N/2 zero' does not exist for
    the p = 6 series and its exponent was underived."""
    for n in (2, 3):
        r4 = next(r for r in bench["points"]
                  if r["p"] == 4 and r["N"] == 10 and r["n"] == n)
        r6 = next(r for r in bench["points"]
                  if r["p"] == 6 and r["N"] == 22 and r["n"] == n)
        assert abs(r4["lambda_chord"] - r6["lambda_chord"]) < 0.1
        # the honest negative: raw w2 larger at larger N (matched lambda)
        assert r6["w2"] > r4["w2"]
