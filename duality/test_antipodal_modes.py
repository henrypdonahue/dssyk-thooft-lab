"""Tests for the antipodal mode separation (antipodal_modes.json).

Pinned here: the sign structure, which is a flipped positive component
plus a smaller unflipped negative one wherever the flip operates; the
REFUTATION of the harmonic-tower mechanism, where kappa_2/kappa_1 is far
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


def test_harmonic_tower_refuted(data):
    """The even-harmonic mechanism predicts TWO things: a rate
    kappa_2 = 2 kappa_1, and a second harmonic that is subdominant.  The
    weights refute it outright.

    The refutation leads with the weight comparison because that is the
    robust half: |c_2| / |c_1| ~ 1 at both ends of every window, so there
    is no harmonic hierarchy at all.  A tower would need |c_2| << |c_1|.
    """
    for h in data["harmonics"]:
        for key in ("c2_over_c1_first", "c2_over_c1_last"):
            assert 0.5 < h[key] < 2.0, (h["lam"], key, h[key])


def test_harmonic_rate_ratio_is_contaminated_not_evidence(data):
    """The companion rate comparison kappa_2/kappa_1 is NOT a growth-rate
    measurement, and this test records why rather than quoting it as one.

    The kappa_1 > 0.05 gate does not guarantee a measurable growth window:
    on every row that passes it, |c_1| is non-monotone across the window
    (lambda = 1.0: 16.7, 59.9, 80.6, 65.3; lambda = 0.6: 965, 1570, 501,
    1230), so the fitted kappa are window averages over a rising-then-
    falling curve.  The ratios do come out near 1 rather than 2, which is
    consistent with the refutation, but the weight comparison above is
    what carries it."""
    gated = [h for h in data["harmonics"] if h["kappa_1"] > 0.05]
    assert gated, "no row passes the rate gate"
    for h in gated:
        c1 = h["abs_c1"]
        rising = all(b > a for a, b in zip(c1, c1[1:]))
        assert not rising, (h["lam"], c1)      # the contamination, pinned
        ratio = h["kappa_2"] / h["kappa_1"]
        assert 0.5 < ratio < 1.6, (h["lam"], ratio)   # near 1, not 2


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
    the shift-circle DFT (the primary, fitting-free experiment behind the
    harmonic-tower refutation) would keep the suite green until someone
    regenerated the data.  The engine is deterministic; one committed
    harmonics row reproduces (~2.5 s).

    TOLERANCE.  The chain is a chord transfer-matrix build, a complex
    contour evaluation, then a DFT.  It is deterministic on one machine,
    but the last two digits move between platforms (different libm and
    BLAS): the measured macOS-arm64 against linux-x86-64 spread on this
    row is 1.2e-12 relative.  So rel=1e-9, three orders above that spread
    and nine below anything a real regression would do -- a wrong sign,
    weight or contour moves these by O(1).  A 1e-8 perturbation of the
    producer is still caught."""
    import antipodal_modes

    want = max(data["harmonics"], key=lambda h: h["lam"])
    got = antipodal_modes.harmonics(want["lam"], betaJ=want["betaJ"],
                                    M=want["M"], thLs=tuple(want["thLs"]))
    for key in ("abs_c0", "abs_c1", "abs_c2", "abs_c3"):
        assert got[key] == pytest.approx(want[key], rel=1e-9, abs=1e-12)
    for key in ("kappa_1", "kappa_2", "kappa_3",
                "c2_over_c1_last", "c2_over_c1_first"):
        assert got[key] == pytest.approx(want[key], rel=1e-9, abs=1e-12)
