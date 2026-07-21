#!/usr/bin/env python3
"""Correctness gates for the charge-sector-blocked moment pipeline
(moments_fast.py): every structural claim in its derivation is asserted
against the dense path (moments_pipeline / dirac), sharing no blocked code.

Run from this directory:  pytest -q test_moments_fast.py -m 'not slow'
(the Nc = 10 end-to-end gate is slow: dense reference takes ~2 min).
"""

import time

import numpy as np
import pytest

from dirac import (draw_couplings, dirac_hamiltonian, annihilators, mn_tower,
                   particle_hole, cp_conjugate)
from moments_pipeline import instance_moments
from moments_fast import (sector_indices, blocked_hamiltonian,
                          blocked_annihilators, blocked_particle_hole,
                          cp_conjugate_blocked, blocked_mn_tower,
                          instance_moments_blocked, ensemble_moments_blocked,
                          rows_from_instances)

NC, P, SEED = 6, 4, 11


def _setup(Nc=NC, p=P, seed=SEED):
    rng = np.random.default_rng(seed)
    couplings = draw_couplings(Nc, p, rng, "c_symmetric")
    idx, pos, sec = sector_indices(Nc)
    return couplings, idx, pos, sec


def _rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def test_blocked_hamiltonian_is_dense_slice():
    """H conserves Q: dense H vanishes off the popcount blocks, and the
    blocked assembly reproduces each diagonal block exactly."""
    couplings, idx, pos, sec = _setup()
    H = dirac_hamiltonian(NC, P, couplings=couplings)
    Hb = blocked_hamiltonian(NC, P, couplings, idx, pos, sec)
    off = H.copy()
    for q in range(NC + 1):
        block = H[np.ix_(idx[q], idx[q])]
        assert np.abs(Hb[q] - block).max() < 1e-12
        off[np.ix_(idx[q], idx[q])] = 0.0
    assert np.abs(off).max() == 0.0          # exactly zero off the blocks


def test_derivative_chain_preserves_corner_and_tower_blocks():
    """c^(n) = i[H, c^(n-1)] stays in the (q-1, q) corner (asserted on the
    dense chain), and the blocked M_n tower equals the dense tower's
    diagonal blocks."""
    couplings, idx, pos, sec = _setup()
    H = dirac_hamiltonian(NC, P, couplings=couplings)
    cs = annihilators(NC)
    # corner-preservation of the dense chain for one mode
    Cn = cs[2]
    for n in range(1, 4):
        Cn = 1j * (H @ Cn - Cn @ H)
        leak = Cn.copy()
        for q in range(1, NC + 1):
            leak[np.ix_(idx[q - 1], idx[q])] = 0.0
        assert np.abs(leak).max() < 1e-10 * max(np.abs(Cn).max(), 1e-300)
    # blocked tower vs dense tower blocks
    towers_dense = mn_tower(H, cs, 3)
    Hb = blocked_hamiltonian(NC, P, couplings, idx, pos, sec)
    cmaps = blocked_annihilators(NC, idx, pos)
    towers = blocked_mn_tower(NC, Hb, cmaps, 3)
    for n in range(4):
        scale = max(np.abs(towers_dense[n]).max(), 1e-300)
        for q in range(NC + 1):
            block = towers_dense[n][np.ix_(idx[q], idx[q])]
            assert np.abs(towers[n][q] - block).max() < 1e-10 * scale
    # sanity: M_0 = Q blockwise (popcount on the diagonal)
    for q in range(NC + 1):
        assert np.abs(towers[0][q] - q * np.eye(len(idx[q]))).max() < 1e-12


def test_blocked_particle_hole_matches_dense():
    """C is a signed complement permutation (sector q -> Nc - q mirror), and
    the blocked CP conjugation of a block-diagonal operator equals the dense
    U X U^dag on every block."""
    couplings, idx, pos, sec = _setup()
    U = particle_hole(NC)
    ph = blocked_particle_hole(NC, idx, pos)
    # U supported ONLY on the mirror blocks, and matching signs/permutation
    for q in range(NC + 1):
        r, sg = ph[q]
        B = U[np.ix_(idx[NC - q], idx[q])]
        Bref = np.zeros_like(B)
        Bref[r, np.arange(len(r))] = sg
        assert np.abs(B - Bref).max() == 0.0
    mask = np.zeros_like(U)
    for q in range(NC + 1):
        mask[np.ix_(idx[NC - q], idx[q])] = 1.0
    assert np.abs(U * (1 - mask)).max() == 0.0
    # blocked conjugation of a genuinely non-Hermitian block-diagonal op (M_3)
    H = dirac_hamiltonian(NC, P, couplings=couplings)
    M3 = mn_tower(H, annihilators(NC), 3)[3]
    M3c_dense = cp_conjugate(U, M3)
    Hb = blocked_hamiltonian(NC, P, couplings, idx, pos, sec)
    cmaps = blocked_annihilators(NC, idx, pos)
    M3b = blocked_mn_tower(NC, Hb, cmaps, 3, keep={3})[3]
    M3c = cp_conjugate_blocked(NC, ph, M3b)
    scale = max(np.abs(M3).max(), 1e-300)
    for q in range(NC + 1):
        block = M3c_dense[np.ix_(idx[q], idx[q])]
        assert np.abs(M3c[q] - block).max() < 1e-10 * scale


def _assert_instances_match(Nc, p, seed, tol=1e-9):
    t0 = time.perf_counter()
    dense = instance_moments(Nc, p, seed)
    t_dense = time.perf_counter() - t0
    t0 = time.perf_counter()
    blocked = instance_moments_blocked(Nc, p, seed)
    t_blocked = time.perf_counter() - t0
    print(f"\n  Nc={Nc} seed={seed}: dense {t_dense:.2f} s, "
          f"blocked {t_blocked:.2f} s (x{t_dense / t_blocked:.1f})")
    assert set(dense) == set(blocked)
    for key in dense:
        d, b = dense[key], blocked[key]
        for name in ("mu0", "mu2", "mu4", "escale"):
            assert _rel(d[name], b[name]) < tol, (key, name, d[name], b[name])
        # mu1 is an exact zero of the symmetrized weights: compare both on
        # the pipeline's normalized odd-moment scale
        scale = np.sqrt(d["mu0"] * d["mu2"])
        assert abs(d["mu1"]) / scale < 1e-10
        assert abs(b["mu1"]) / scale < 1e-10


def test_instance_moments_identical_nc6():
    """End-to-end gate: identical mu_k and escale vs the dense pipeline at
    Nc = 6 (both channels, n = 2, 3), to better than 1e-9 relative."""
    _assert_instances_match(6, 4, 100)


def test_instance_moments_identical_nc8():
    """Same end-to-end gate at Nc = 8."""
    _assert_instances_match(8, 4, 100)


@pytest.mark.slow
def test_instance_moments_identical_nc10():
    """Same end-to-end gate at Nc = 10 (dense reference ~2 min)."""
    _assert_instances_match(10, 4, 100)


def test_ensemble_rows_match_pipeline_format():
    """ensemble_moments_blocked emits rows with exactly the pipeline's keys
    and sane invariants (kurt >= 1 is Cauchy-Schwarz, exact)."""
    from moments_pipeline import ensemble_moments
    rows_d = ensemble_moments(Nc=6, p=4, n_inst=3, seed0=7)
    rows_b = ensemble_moments_blocked(Nc=6, p=4, n_inst=3, seed0=7)
    assert len(rows_b) == len(rows_d) == 4
    for rd, rb in zip(rows_d, rows_b):
        assert set(rd) == set(rb)
        for k in ("Nc", "p", "n", "channel", "n_inst"):
            assert rd[k] == rb[k]
        for k in ("mu0", "mu0_err", "w2", "w2_err", "kurt", "kurt_err"):
            assert _rel(rd[k], rb[k]) < 1e-8, (k, rd[k], rb[k])
        assert rb["odd_moment_check"] < 1e-10
        assert rb["kurt"] >= 1.0
        assert rb["w2"] > 0


def test_rows_from_instances_sem_ddof1():
    """SEM is std(ddof=1)/sqrt(n) of the per-instance invariants (house
    convention), checked against a direct recomputation."""
    insts = [instance_moments_blocked(6, 4, s) for s in (3, 4, 5)]
    rows = rows_from_instances(6, 4, insts)
    r = next(r for r in rows if r["n"] == 2 and r["channel"] == "right")
    w2s = [i[(2, "right")]["mu2"] /
           (i[(2, "right")]["mu0"] * i[(2, "right")]["escale"] ** 2)
           for i in insts]
    assert _rel(r["w2"], np.mean(w2s)) < 1e-12
    assert _rel(r["w2_err"], np.std(w2s, ddof=1) / np.sqrt(3)) < 1e-12
