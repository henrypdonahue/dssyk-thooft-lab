"""Tests for the exact N = infinity singlet-channel moments (rung 15).

The load-bearing verdict encoded here: at strict N = infinity, fixed
lambda, the H-derived singlet bilinear channels A_n = sum_i psi_i ad^n
psi_i carry ZERO spectral weight away from omega = 0 -- mu_2 = mu_4 = 0 to
machine precision for n = 1, 2, 3 at every lambda tested, while the
random-operator O_Delta channels propagate (mu_2 = 2(1-q^Delta) != 0).
The word-expansion assembly behind this is validated as a per-instance
operator identity in syk-self-averaging/test_mn_majorana.py; the engine
amplitudes behind it are validated against the brute-force enumerator in
test_chord.py.  Together with the ED bench falloff, this is the exact
chord-level confirmation of the rung-17 finding (the tower cancels at
leading order): the M_n tower is a subleading-in-1/N effect at every
coupling, not just at large q.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from chord import q_of_lambda, spectral_lines
from chord_moments import channel_moments, odelta_moments

HERE = Path(__file__).resolve().parent


@pytest.mark.parametrize("lam", [2.0, 1.0, 0.5])
@pytest.mark.parametrize("n", [1, 2, 3])
def test_an_channels_conserved_at_leading_order(lam, n):
    """mu_2 and mu_4 of every A_n channel vanish identically in the chord
    amplitudes at finite lambda and finite Delta_psi = 1/p -- not a small
    number: an algebraic cancellation (checked at machine precision
    relative to the pre-cancellation gross magnitude)."""
    q = q_of_lambda(lam)
    ch = channel_moments(q, 0.25, n)
    assert ch["mu0"] > 0
    assert ch["conserved"]
    assert ch["cancellation2"] < 1e-12       # machine-precision cancellation
    assert abs(ch["mu4"]) < 1e-11 * max(ch["gross4"], 1.0)


@pytest.mark.parametrize("dpsi", [0.1, 0.25, 0.5, 1.0])
def test_an_conservation_any_matter_weight(dpsi):
    """The cancellation holds at ANY matter weight (added post-review: the
    'any weight' clause of CHORD.md was previously asserted only at
    Delta_psi = 1/4)."""
    q = q_of_lambda(1.0)
    for n in (1, 2):
        ch = channel_moments(q, dpsi, n)
        assert ch["conserved"]
        assert ch["cancellation2"] < 1e-12


@pytest.mark.parametrize("q,delta", [(0.5, 0.25), (0.8, 0.5), (0.3, 1.0)])
def test_odelta_closed_forms(q, delta):
    """The propagating O_Delta channel: engine Lehmann moments equal the
    closed forms mu_2 = 2(1-qt), mu_4 = 2(2+q)+6-8(2+q)qt+6(1+q)qt^2."""
    om, w = spectral_lines(q, delta, 40)
    mu2 = float(np.sum(np.real(w) * om ** 2))
    mu4 = float(np.sum(np.real(w) * om ** 4))
    cf = odelta_moments(q, delta)
    assert mu2 == pytest.approx(cf["mu2"], rel=1e-10)
    assert mu4 == pytest.approx(cf["mu4"], rel=1e-10)


def test_odelta_kurtosis_semiclassical_limit():
    """kurt(lambda -> 0) -> 3 + 1/Delta (Fourier transform of
    sech^{2Delta}), approached linearly in lambda."""
    for delta in (0.25, 0.5, 1.0):
        target = 3.0 + 1.0 / delta
        k1 = odelta_moments(q_of_lambda(0.2), delta)["kurt"]
        k2 = odelta_moments(q_of_lambda(0.1), delta)["kurt"]
        k3 = odelta_moments(q_of_lambda(0.05), delta)["kurt"]
        assert abs(k3 - target) < abs(k2 - target) < abs(k1 - target)
        assert abs(k3 - target) < 0.2 * target


def test_committed_moments_json():
    """The committed campaign JSON: A_1 residual at machine precision on
    every row (the identity holds in the chord algebra, not only at
    lambda -> 0), and the O_Delta w2 closed form on every row."""
    with open(HERE / "chord_moments.json") as fh:
        data = json.load(fh)
    assert len(data["points"]) >= 8
    for row in data["points"]:
        assert row["a1_residual"] < 1e-12
        q = row["q"]
        cf = 2.0 * (1.0 - q ** (2.0 / row["p"]))
        assert row["channels"]["O_delta=2/p"]["w2"] == pytest.approx(cf,
                                                                     rel=1e-12)
        for n in (2, 3):
            ch = row["channels"][f"n={n}"]
            assert ch["conserved"]
            assert ch["cancellation2"] < 1e-12


def test_ed_bench_scope():
    """Cross-module scope statement, post-review: the chord computation is
    the i != j flavor-PAIR amplitude; the ED bench measures the full
    physical channel, which at dense-ED sizes is flavor-diagonal-dominated
    (the pair part of mu_0 is a few percent and opposite in sign at
    N = 12-14, p = 4 -- measured during the adversarial review).  So the
    bench does NOT yet probe the chord cancellation: at the near-matched-
    lambda pair the RAW w2 rises with N.  What the bench establishes is
    the raw fixed-p turn-on toward the 't Hooft corner (asserted in
    syk-self-averaging/test_mn_majorana.py); here the matched pair's
    existence and the honest negative are pinned."""
    path = HERE.parent / "syk-self-averaging" / "mn_majorana_bench.json"
    with open(path) as fh:
        bench = json.load(fh)
    r4 = next(r for r in bench["points"]
              if r["p"] == 4 and r["N"] == 10 and r["n"] == 2)
    r6 = next(r for r in bench["points"]
              if r["p"] == 6 and r["N"] == 22 and r["n"] == 2)
    assert abs(r4["lambda_chord"] - r6["lambda_chord"]) < 0.1
    assert r6["w2"] > r4["w2"]      # the approach to the chord zero is
    #                                 NOT visible at dense-ED sizes
