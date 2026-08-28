"""Tests for the chord-side Euclidean fold at N = infinity.

The committed chord_fold.json is the N = infinity answer to the question
fold_scan.json, decider #1, posed at ED sizes.  Three verdicts are
encoded here.

* The UNFOLDED crossed OTOC ratio is lambda-stable and real.  The control
  behaves at N = infinity exactly as it did in the ED scan.

* The FOLDED connected fraction stays O(1) as lambda -> 0.  There is no
  lambda ~ 1/N suppression.  So the N-scaled folded ratio diverges away
  from the continued Streicher closed form, like 1/lambda.  Harlow-Zhao's
  assumption that continuation commutes with the semiclassical limit
  fails IN THE DOUBLE-SCALED THEORY ITSELF.  The ED N-obstruction was
  real, and at N = infinity it is a lambda-obstruction.

* On the fold, the 1/N convergence itself degrades.  At matched lambda
  the ED folded values sit 1.5-1.6x above the exact N = infinity ones,
  while unfolded values are within 10%.  The fold amplifies finite-N
  corrections, which is consistent with a soft-spectral-edge mechanism.
"""

import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def data():
    with open(HERE / "chord_fold.json") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def ed_scan():
    with open(HERE.parent / "syk-self-averaging" / "fold_scan.json") as fh:
        return json.load(fh)


def test_truncation_clean(data):
    for pt in data["points"]:
        assert pt["truncation_drift_fold"] < 1e-4


def test_unfolded_control_stable_and_real(data):
    """The unfolded crossed OTOC ratio: lambda-stable within a narrow band,
    scrambling-side sign.  (The im == 0 check is a SANITY check
    only: for real-theta configurations every chord segment and matrix is
    real, so reality is manifest by construction -- the nontrivial reality
    statement is the Lorentzian growth branch, asserted in
    test_fold_breaks_reality_like_the_continuation.)"""
    res = [pt["otoc"]["re"] for pt in data["points"]]
    ims = [abs(pt["otoc"]["im"]) for pt in data["points"]]
    assert all(-1.75 < r < -1.60 for r in res)
    assert max(res) - min(res) < 0.1
    assert all(i < 1e-8 for i in ims)


def test_folded_connected_fraction_not_suppressed(data):
    """THE decider: as lambda -> 0 the unfolded connected fraction dies
    ~ lambda (the 1/N hierarchy), the folded one does NOT -- it stays O(1)
    and grows."""
    pts = sorted(data["points"], key=lambda p: -p["lam"])
    for pt in pts:
        lam = pt["lam"]
        otoc_frac = abs(pt["otoc_connected_fraction"]["re"])
        fold_frac = abs(pt["fold_connected_fraction"]["re"])
        assert fold_frac > 0.15
        # unfolded fraction = |otoc ratio| * lam/(2 p^2), tiny at small lam
        assert otoc_frac < 0.25 * lam
    small = [pt for pt in pts if pt["lam"] <= 0.8]
    for pt in small:
        assert (abs(pt["fold_connected_fraction"]["re"])
                > 5 * abs(pt["otoc_connected_fraction"]["re"]))
    # and the folded fraction GROWS toward the semiclassical corner
    fr = [abs(pt["fold_connected_fraction"]["re"]) for pt in pts]
    assert fr[-1] > fr[0]


def test_folded_ratio_departs_from_continuation(data):
    """The N-scaled folded ratio diverges from the continued closed form as
    lambda falls (|chord/closed| exceeds 10 by lambda = 0.5), while the
    unfolded ratio/closed stays O(1)."""
    pts = sorted(data["points"], key=lambda p: -p["lam"])
    ratios = [abs(pt["fold"]["re"] / pt["fold"]["closed_re"]) for pt in pts]
    assert all(ratios[i + 1] > ratios[i] - 1e-9 for i in range(len(ratios) - 1))
    assert ratios[-1] > 10.0
    for pt in pts:
        assert abs(pt["otoc"]["re"] / pt["otoc"]["closed_re"]) < 2.0


def test_fold_breaks_reality_like_the_continuation(data):
    """Lorentzian branch: the unfolded growth curve stays exactly real; the
    folded one develops Im < 0 -- the same symmetry-breaking pattern the ED
    lab measured, amplified in magnitude.

    The continued closed form's Im is asserted only in MAGNITUDE (growing
    with thL).  The folded angles lie far outside Streicher's principal
    window (theta_12^- = -3 pi), so its branch -- and hence its overall
    sign -- is a choice, not a prediction.  See
    test_antiscrambling.py::test_fold_angles_leave_principal_window."""
    for pt in data["points"]:
        for key, v in pt.items():
            if isinstance(v, dict) and key.startswith("growth_thL"):
                assert abs(v["im"]) < 1e-8
        ims = []
        for thL in (0.8, 1.2):
            v = pt[f"fold_thL={thL}"]
            assert v["im"] < -1e-3              # chord side: unambiguous
            assert v["closed_im"] != 0.0
            ims.append(abs(v["closed_im"]))
        assert ims[1] > ims[0]


def test_ed_bridge_unfolded_close_folded_far(data, ed_scan):
    """At matched lambda with the exact finite-N bridge: the chord and ED
    UNFOLDED ratios agree to ~10%; the FOLDED ones disagree by ~1.5-1.6x
    (growing N=8 -> 18): the fold amplifies finite-N corrections -- the
    quantitative statement that the fold destroys the 1/N hierarchy, from
    both ends."""
    ed = {(r["p"], r["N"]): r for r in ed_scan["points"]}
    for row in data["ed_bridge_points"]:
        p, N = row["p"], row["N_ed"]
        ed_otoc = ed[(p, N)]["otoc"]["re"]
        ed_fold = ed[(p, N)]["fold"]["re"]
        assert abs(row["otoc"]["re"] / ed_otoc - 1.0) < 0.12
        fold_gap = ed_fold / row["fold"]["re"]
        assert 1.3 < fold_gap < 1.9
    gap8 = ed[(4, 8)]["fold"]["re"] / next(
        r for r in data["ed_bridge_points"] if r["N_ed"] == 8)["fold"]["re"]
    gap18 = ed[(6, 18)]["fold"]["re"] / next(
        r for r in data["ed_bridge_points"] if r["N_ed"] == 18)["fold"]["re"]
    assert gap18 > gap8


def test_configs_match_ed_lab(data, ed_scan):
    """Same theta configurations and fold map as the ED laboratory."""
    assert data["config_crossed"] == ed_scan["points"][0]["config_crossed"]
    assert data["fold_thetas"] == ed_scan["points"][0]["fold_thetas"]


def test_soft_edge_mechanism(ed_scan):
    """The measured mechanism behind the folded finite-N gap
    (syk-self-averaging/fold_edge_probe.py): the per-instance folded
    connected fraction is dominated by the instance's UPPER spectral edge
    (|corr| ~ 0.86), the unfolded control is not, and the folded scatter
    is ~35x larger -- folded correlators at ED sizes measure edge
    fluctuations, which the sharp-edged chord theory lacks."""
    path = HERE.parent / "syk-self-averaging" / "fold_edge_probe.json"
    with open(path) as fh:
        d = json.load(fh)
    assert d["corr_Emax_folded"] < -0.7
    assert abs(d["corr_Emax_unfolded"]) < 0.4
    assert abs(d["corr_Emax_folded"]) > 2 * abs(d["corr_negEmin_folded"])
    assert (d["folded_frac_std"] / d["unfolded_frac_std"]) > 10


def test_chord_fold_reproduces_its_committed_row(data):
    """Execute the producer, do not just read its output.

    Every other test here reads the committed JSON, which leaves
    chord_fold.py itself unexercised -- a sign or scale corruption in the
    contour assembly would keep the suite green until someone regenerated
    the data.  The engine is deterministic, so the cheapest committed row
    (largest lambda, smallest chord lattice, ~0.1 s) reproduces exactly.
    """
    import chord_fold

    row = max(data["points"], key=lambda p: p["lam"])
    out = chord_fold.run_point(row["lam"], row["p"], row["betaJ"],
                               row["nmax"])
    for key in ("toc", "otoc", "fold"):
        assert out[key]["re"] == pytest.approx(row[key]["re"],
                                               rel=1e-9, abs=1e-12)
    for key in ("fold_connected_fraction", "otoc_connected_fraction"):
        assert out[key]["re"] == pytest.approx(row[key]["re"],
                                               rel=1e-9, abs=1e-12)
    # the structure the verdict rests on, recomputed from scratch
    assert out["toc"]["re"] > 0
    assert out["otoc"]["re"] < 0
    assert out["otoc"]["im"] == pytest.approx(0.0, abs=1e-12)
