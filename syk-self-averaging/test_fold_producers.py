"""The fold producers are executed, not just read back from JSON.

`fold_bench.py` and `fold_edge_probe.py` write the artifacts that carry the
fold verdicts, but every other test in the repo reads only the
committed JSON.  That leaves the producing code itself unexercised: a sign
or scale corruption inside `bilocal4_crossed`, `hz_fold` or the connected
ratio would leave the whole suite green and only show up the next time
someone regenerated the data.

These tests close that loop.  Both producers are deterministic given their
seed, so the committed rows can be reproduced bit-for-bit rather than
merely smoke-tested.
"""

import json
from pathlib import Path

import numpy as np
import pytest

import fold_bench
import fold_edge_probe

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def scan_row():
    path = HERE / "fold_scan.json"
    if not path.exists():
        pytest.skip("run fold_bench.py --scan first")
    pts = json.loads(path.read_text())["points"]
    return min((p for p in pts if p["p"] == 4), key=lambda q: q["N"])


def test_fold_bench_reproduces_its_committed_row(scan_row):
    """The cheapest committed scan point (p = 4, N = 8, 192 instances, ~3 s)
    recomputed from scratch must match the artifact exactly -- same seed,
    same ensemble, same numbers."""
    out = fold_bench.run_point(
        N=scan_row["N"], p=scan_row["p"], betaJ=scan_row["betaJ"],
        n_inst=scan_row["n_inst"], seed=scan_row["seed"],
        fold_n=scan_row["fold_n"])
    for key in ("toc", "otoc", "fold"):
        assert out[key]["re"] == pytest.approx(scan_row[key]["re"],
                                               rel=1e-12, abs=1e-12)
    # and the qualitative structure the verdicts rest on
    assert out["toc"]["re"] > 0                  # time-ordered
    assert out["otoc"]["re"] < 0                 # crossed: scrambling sign
    assert out["fold"]["re"] < 0
    assert np.isfinite(out["fold"]["re"])        # bounded spectrum (HZ)


def test_fold_bench_growth_is_monotone_at_small_size():
    """Re script-F/F_d decreases with Lorentzian time -- finite-N SYK
    scrambles -- at a size small enough to run in the fast suite."""
    out = fold_bench.run_point(N=8, p=4, betaJ=2.0, n_inst=8, seed=11)
    seq = [out["otoc"]["re"]] + [out[f"growth_thL={t}"]["re"]
                                 for t in (0.4, 0.8, 1.2)]
    assert all(b < a for a, b in zip(seq, seq[1:])), seq


def test_hz_fold_matches_the_committed_angles(scan_row):
    """The fold map the producer applies is the one the artifacts record."""
    folded = fold_bench.hz_fold(fold_bench.CONFIG_CROSSED, scan_row["fold_n"])
    assert list(np.real(folded)) == pytest.approx(scan_row["fold_thetas"])


def test_fold_edge_probe_runs_and_finds_the_edge_correlation():
    """The fold's finite-N mechanism: the per-instance folded connected
    fraction tracks the sample's UPPER spectral edge, while the unfolded
    control does not.  At N = 8 with 24 instances (~0.05 s) the sign and
    the ordering are already clear; the committed campaign quotes
    -0.86 vs +0.15 at N = 12 with 48."""
    out = fold_edge_probe.run(N=8, p=4, n_inst=24, seed=3)
    assert out["corr_Emax_folded"] < -0.3
    assert abs(out["corr_Emax_folded"]) > abs(out["corr_Emax_unfolded"])
    assert out["folded_frac_std"] > out["unfolded_frac_std"]


def test_committed_edge_probe_verdict():
    """The committed artifact says the same thing, more sharply."""
    path = HERE / "fold_edge_probe.json"
    if not path.exists():
        pytest.skip("run fold_edge_probe.py first")
    d = json.loads(path.read_text())
    assert d["corr_Emax_folded"] < -0.8
    assert abs(d["corr_Emax_unfolded"]) < 0.3
    assert abs(d["corr_negEmin_folded"]) < abs(d["corr_Emax_folded"])
    assert d["folded_frac_std"] > 10 * d["unfolded_frac_std"]
