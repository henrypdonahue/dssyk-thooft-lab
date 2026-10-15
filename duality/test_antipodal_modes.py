"""Tests for the antipodal mode separation (antipodal_modes.json).

Pinned here: the sign structure, which is a flipped positive component
plus a smaller unflipped negative one wherever the flip operates; the
the harmonic-tower test on the physical half-circle, where kappa_2/kappa_1 is far
from 2 and the harmonic weights are comparable; and the recorded fit
degeneracies.
"""

import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def data():
    with open(HERE / "antipodal_modes.json") as fh:
        return json.load(fh)


def test_sign_structure_where_flip_operates(data):
    """For lambda in [0.8, 1.5] the constrained two-mode fit finds a
    positive flipped amplitude and a smaller negative unflipped one in
    the antipodal frame."""
    for pt in data["points"]:
        if 0.8 <= pt["lam"] <= 1.5:
            con = pt["antipodal"]["constrained"]
            assert con["A1"] > 0 > con["A2"]
            assert abs(con["A2"]) < abs(con["A1"])
            # and the fit actually describes the curve
            assert pt["antipodal"]["resid_over_range_con"] < 0.05


def test_lambda06_fit_degeneracy_recorded(data):
    """The lambda = 0.6 two-mode fit collapses (amplitudes blow up with
    swapped signs) -- kept as the recorded instrument failure, not
    silently dropped."""
    pt = next(p for p in data["points"] if p["lam"] == 0.6)
    con = pt["antipodal"]["constrained"]
    assert not (con["A1"] > 0 > con["A2"]) or abs(con["A1"]) > 10


def test_second_harmonic_subdominant_not_a_clean_tower(data):
    """The even-harmonic mechanism predicts kappa_2 = 2 kappa_1 with
    |A_2| << |A_1|.  On the physical half-circle the second harmonic is
    subdominant (|A_2/A_1| ~ 0.3-0.6, model-dependent) and grows faster
    than the first, kappa_2/kappa_1 ~ 1.2-1.5: between the tower's 2 and
    the ~1 that the withdrawn full-circle DFT reported.  Neither refuted
    nor established; the model dependence (two vs three harmonics) and
    the 5-10% residuals are recorded, and the fit is re-run live at one
    lambda."""
    import numpy as np
    from antipodal_modes import harmonics
    for h in data["harmonics"]:
        assert h["deltas"][-1] < h["pole_delta"]          # physical half only
        for nh, m in h["models"].items():
            assert 0.2 < min(m["A2_over_A1"]) and max(m["A2_over_A1"]) < 0.7
            assert 1.1 < m["kappa2_over_kappa1"] < 1.7, (h["lam"], nh, m["kappa2_over_kappa1"])
            assert max(m["resid_over_range"]) < 0.12
            assert m["kappa_1"] > 0 and m["kappa_2"] > m["kappa_1"]
        lo, hi = h["kappa2_over_kappa1_range"]
        assert lo < 2.0 and hi < 2.0                         # not a clean tower
    row = next(h for h in data["harmonics"] if h["lam"] == 1.0)
    live = harmonics(1.0)
    assert np.allclose(live["R"], row["R"], rtol=1e-9, atol=1e-12)
    for nh in ("2", "3"):
        assert live["models"][nh]["kappa2_over_kappa1"] == pytest.approx(
            row["models"][nh]["kappa2_over_kappa1"], rel=1e-9)


def test_censored_fits_are_flagged_and_stay_flagged(data):
    """`kappa_at_grid_max` marks a fit whose kappa saturated the search
    grid ceiling -- an instrument failure, not a rate.  Pinning the
    censored set means a regeneration that silently changes which fits
    saturate shows up as a test failure instead of looking identical to a
    clean run."""
    flags = {}
    for pt in data["points"]:
        for frame in ("antipodal", "unshifted"):
            for fit in ("constrained", "free"):
                flags[(pt["lam"], frame, fit)] = \
                    pt[frame][fit]["kappa_at_grid_max"]
    assert sum(flags.values()) == 9, sorted(k for k, v in flags.items() if v)
    # every free-fit row saturates except the two smallest lambda
    assert flags[(0.8, "antipodal", "constrained")] is True
    # the rows test_sign_structure_where_flip_operates relies on for their
    # AMPLITUDES are fine; only lambda = 0.8's kappa is censored
    for lam in (1.5, 1.2, 1.0):
        assert flags[(lam, "antipodal", "constrained")] is False


def test_harmonics_reproduces_its_committed_row(data):
    """Execute the producer, do not just read its output.

    Every other test here reads antipodal_modes.json, so a regression in
    the half-circle harmonic fit (the primary experiment behind the
    harmonic-tower verdict) would keep the suite green until someone
    regenerated the data.  The engine is deterministic; one committed
    harmonics row reproduces (~2.5 s).

    TOLERANCE.  The chain is a chord transfer-matrix build, a complex
    contour evaluation, then a least-squares fit.  It is deterministic on one machine,
    but the last two digits move between platforms (different libm and
    BLAS): the measured macOS-arm64 against linux-x86-64 spread on this
    row is 1.2e-12 relative.  So rel=1e-9, three orders above that spread
    and nine below anything a real regression would do -- a wrong sign,
    weight or contour moves these by O(1).  A 1e-8 perturbation of the
    producer is still caught."""
    import antipodal_modes

    want = max(data["harmonics"], key=lambda h: h["lam"])
    got = antipodal_modes.harmonics(want["lam"], betaJ=want["betaJ"],
                                    n_delta=want["n_delta"],
                                    thLs=tuple(want["thLs"]))
    for gr, wr in zip(got["R"], want["R"]):
        assert gr == pytest.approx(wr, rel=1e-9, abs=1e-12)
    for nh in want["models"]:
        for key in ("abs_A0", "abs_A1", "abs_A2", "A2_over_A1",
                    "kappa_1", "kappa_2", "kappa2_over_kappa1"):
            assert got["models"][nh][key] == pytest.approx(
                want["models"][nh][key], rel=1e-9, abs=1e-12), (nh, key)
