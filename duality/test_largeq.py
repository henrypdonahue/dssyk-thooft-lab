#!/usr/bin/env python3
"""Tests for the large-q anchor (rung 17): conventions bridge, closed forms
vs the ED bench, and the structure of the would-be tower."""

import json
from math import comb, factorial
from pathlib import Path

import numpy as np
import pytest

from largeq_anchor import (v_of_betaJ, G_infT_repo, quantization_sets,
                           gr_tower)

_BENCH = Path(__file__).resolve().parent.parent / "syk-self-averaging" / "largeq_bench.json"


def test_v_solver_limits():
    """v -> betaJ/pi at weak coupling, v -> 1 in the conformal limit."""
    assert v_of_betaJ(0.01) == pytest.approx(0.01 / np.pi, rel=1e-2)
    assert v_of_betaJ(1000.0) > 0.99
    # the defining relation itself
    for bj in (0.5, 3.0, 30.0):
        v = v_of_betaJ(bj)
        assert np.pi * v / np.cos(np.pi * v / 2) == pytest.approx(bj, rel=1e-10)


@pytest.mark.skipif(not _BENCH.exists(), reason="run largeq_bench.py first")
def test_conventions_bridge_H2():
    """The exact repo <-> MS bridge implies <H^2>/dim = C(N,p) p!/N^{p-1};
    the bench must confirm it within errors (locks the scriptJ = sqrt(2) p
    identification)."""
    for r in json.loads(_BENCH.read_text()):
        N, p = r["N"], r["p"]
        exact = comb(N, p) * factorial(p) / N ** (p - 1)
        assert r["m2_mean"] == pytest.approx(exact, rel=0.05)


@pytest.mark.skipif(not _BENCH.exists(), reason="run largeq_bench.py first")
def test_G_prediction_at_1_over_p_accuracy():
    """G(t) = sech(sqrt(2) p t)^{2/p} must match ED with NO free parameters
    at the O(1/p) level.  (The p = 4 deviations do not shrink with N: the
    formula's intrinsic 1/p error survives the lambda -> 0 limit, and at
    lambda ~ 1 the two correction families partially cancel -- measured
    sequence ~0.07-0.10 vs 1/p = 0.25.  The assertion is the magnitude, not a
    spurious monotonicity.)"""
    devs = {}
    for r in json.loads(_BENCH.read_text()):
        if r["p"] != 4:
            continue
        ts = np.array(r["ts"])
        G = np.array(r["G_re"])
        mask = (ts > 0) & (ts <= 1.0)
        pred = G_infT_repo(ts[mask], r["p"])
        devs[r["N"]] = float(np.max(np.abs(G[mask] - pred)))
    assert len(devs) >= 3
    assert all(d < 0.15 for d in devs.values())    # O(1/p), no free params


def test_quantization_sets_alternate_and_degenerate():
    """CMS (3.3): M^e and M^o interleave (the CP-alternation slot), and at
    T = inf (v -> 0) degenerate to odd/even integers respectively."""
    Me, Mo = quantization_sets(0.05, m_max=8.0)
    for k, m in enumerate(Me[:4]):
        assert m == pytest.approx(2 * k + 1, abs=0.02)   # odd integers
    for k, m in enumerate(Mo[:3]):
        assert m == pytest.approx(2 * k + 2, abs=0.02)   # even integers
    # interleaving at intermediate coupling
    Me, Mo = quantization_sets(0.5, m_max=8.0)
    merged = sorted([(m, "e") for m in Me] + [(m, "o") for m in Mo])
    kinds = [k for _, k in merged]
    assert all(a != b for a, b in zip(kinds, kinds[1:]))


def test_gr_tower_structure():
    """h_n -> 2n+1 + O(1/q) and c_n^2 ~ 1/q^2: the tower decouples at
    leading large q -- the structural finding of rung 17."""
    for q in (8.0, 16.0):
        rows = gr_tower(q, 4)
        for r in rows:
            # 2 eps_n <= 2*(4/2)/q = 4/q, largest at n = 1
            assert r["h"] == pytest.approx(2 * r["n"] + 1, abs=5.0 / q)
    r8 = gr_tower(8.0, 2)[0]["c2"]
    r16 = gr_tower(16.0, 2)[0]["c2"]
    assert r8 / r16 == pytest.approx(4.0, rel=0.01)      # ~ 1/q^2


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
