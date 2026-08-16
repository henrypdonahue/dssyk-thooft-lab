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
import math
from pathlib import Path

import numpy as np
import pytest

from chord import q_of_lambda, spectral_lines
from chord_moments import channel_moments, odelta_moments, MomentTable

HERE = Path(__file__).resolve().parent


@pytest.mark.parametrize("lam", [2.0, 1.0, 0.5])
@pytest.mark.parametrize("n", [1, 2, 3])
def test_an_channels_conserved_at_leading_order(lam, n):
    """mu_2 and mu_4 of every A_n channel vanish identically in the chord
    amplitudes at finite lambda and finite Delta_psi = 1/p -- not a small
    number: an algebraic cancellation (checked at machine precision
    relative to mu_0)."""
    q = q_of_lambda(lam)
    ch = channel_moments(q, 0.25, n)
    assert ch["mu0"] > 0
    assert ch["conserved"]
    assert ch["cancellation2"] < 1e-12       # machine-precision cancellation
    assert abs(ch["mu4"]) < 1e-11 * max(ch["gross4"], 1.0)


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


def test_ed_bench_consistent_with_chord_zero():
    """Cross-module verdict (stated with both halves): the RAW ED w2 rises
    along fixed p as lambda falls (the tower channel turns on toward the
    't Hooft corner -- where the chord description says it is conserved),
    while the p = N/2-artifact-corrected w2 falls with N within each
    fixed-p series and across the near-matched-lambda pair -- the finite-N
    trend consistent with the chord cancellation at fixed lambda.  The
    full assertions live in syk-self-averaging/test_mn_majorana.py; here
    the cross-module anchor pair is pinned."""
    path = HERE.parent / "syk-self-averaging" / "mn_majorana_bench.json"
    with open(path) as fh:
        bench = json.load(fh)

    def corrected(r):
        return r["w2"] / (1.0 - 2.0 * r["p"] / r["N"]) ** 2

    r4 = next(r for r in bench["points"]
              if r["p"] == 4 and r["N"] == 10 and r["n"] == 2)
    r6 = next(r for r in bench["points"]
              if r["p"] == 6 and r["N"] == 22 and r["n"] == 2)
    assert abs(r4["lambda_chord"] - r6["lambda_chord"]) < 0.1
    assert corrected(r6) < 0.6 * corrected(r4)
