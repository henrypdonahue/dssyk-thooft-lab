"""Validation battery for the chord engine (rung-18 toolchain).

Three independent legs hold it up:
  1. the brute-force chord-diagram enumerator (shares only the microscopic
     Wick weights with the engine, none of the transfer-matrix bookkeeping);
  2. the published closed forms: theta-integral 2-pt (review 2407.09396
     (2.30)), mu_2 = 2(1-q^Delta), Touchard-Riordan moments;
  3. the microscopic Krawtchouk check that the chord weight q = e^{-2p^2/N}
     is the right N -> infinity limit of the exact finite-N crossing average
     (no chords involved at all).
Plus the repo bridges: the large-q sech^{2Delta} anchor and the finite-N ED
anchor from fold_bench.json.
"""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from chord import (Contour, Sector0, bruteforce_word, g2_closed_form,
                   g2_engine, spectral_lines, length_channel_lines,
                   n_expectation, rho_theta, energy_theta,
                   sigma_H, lambda_chord_of)

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# 1. engine == brute-force enumerator, exactly, for every supported topology
# ---------------------------------------------------------------------------
WORDS = [
    # 2-pt words
    (['A', 'A'], {'A': 0.7}, set()),
    (['A', 'H', 'H', 'A'], {'A': 0.7}, set()),
    (['H', 'A', 'H', 'H', 'A', 'H'], {'A': 0.3}, set()),
    (['A', 'H', 'H', 'H', 'H', 'A'], {'A': 1.2}, set()),
    # uncrossed (adjacent) 4-pt
    (['A', 'A', 'B', 'B'], {'A': 0.5, 'B': 0.25}, set()),
    (['A', 'H', 'A', 'H', 'B', 'H', 'B', 'H'], {'A': 0.5, 'B': 0.25}, set()),
    # nested 4-pt
    (['A', 'B', 'B', 'A'], {'A': 0.5, 'B': 0.25}, set()),
    (['A', 'H', 'B', 'H', 'H', 'B', 'H', 'A'], {'A': 0.5, 'B': 1.0}, set()),
    # crossed 4-pt (the OTOC topology), bosonic and fermionic
    (['A', 'B', 'A', 'B'], {'A': 0.5, 'B': 0.25}, set()),
    (['A', 'H', 'B', 'H', 'A', 'H', 'B', 'H'], {'A': 0.5, 'B': 0.25}, set()),
    (['A', 'H', 'B', 'H', 'A', 'H', 'B', 'H'], {'A': 0.25, 'B': 0.25},
     {'A', 'B'}),
    (['H', 'A', 'B', 'H', 'H', 'A', 'B', 'H'], {'A': 0.4, 'B': 0.9},
     {'A', 'B'}),
    # same-label pairs of a repeated operator: two A-chords, crossed pattern
    # is not expressible with one label; keep to distinct labels here.
]


@pytest.mark.parametrize("q", [0.0, 0.3, 0.7, 0.95])
@pytest.mark.parametrize("word,deltas,ferm", WORDS)
def test_engine_matches_bruteforce(q, word, deltas, ferm):
    bf = bruteforce_word(q, word, deltas, ferm)
    en = Contour(q, 16, deltas, ferm).apply_word_powers(word).real
    assert en == pytest.approx(bf, rel=1e-12, abs=1e-12)


def test_fermionic_crossed_is_minus_bosonic():
    """In the crossed topology the two matter chords cross exactly once, so
    the fermionic amplitude is exactly minus the bosonic one."""
    q, word, deltas = 0.6, ['A', 'H', 'B', 'H', 'A', 'H', 'H', 'B'], \
        {'A': 0.5, 'B': 0.5}
    b = Contour(q, 16, deltas, set()).apply_word_powers(word).real
    f = Contour(q, 16, deltas, {'A', 'B'}).apply_word_powers(word).real
    assert f == pytest.approx(-b, rel=1e-13)


def test_delta_zero_matter_is_transparent():
    """A Delta = 0 matter chord weighs nothing: any word equals the word
    with that pair deleted."""
    q = 0.55
    with_m = Contour(q, 16, {'A': 0.0, 'B': 0.7}).apply_word_powers(
        ['A', 'H', 'B', 'H', 'A', 'H', 'B', 'H']).real
    without = Contour(q, 16, {'B': 0.7}).apply_word_powers(
        ['H', 'B', 'H', 'H', 'B', 'H']).real
    assert with_m == pytest.approx(without, rel=1e-13)


# ---------------------------------------------------------------------------
# 2. pure-H moments: Touchard-Riordan values and the theta-measure
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q", [0.0, 0.3, 0.7, 0.95])
def test_moments_touchard(q):
    """m_2 = 1, m_4 = 2 + q, m_6 = 5 + 6q + 3q^2 + q^3 (chord counting),
    engine == brute force == closed values."""
    c = Contour(q, 12, {})
    m2 = c.apply_word_powers(['H'] * 2).real
    m4 = c.apply_word_powers(['H'] * 4).real
    m6 = c.apply_word_powers(['H'] * 6).real
    assert m2 == pytest.approx(1.0, abs=1e-14)
    assert m4 == pytest.approx(2.0 + q, rel=1e-14)
    assert m6 == pytest.approx(5.0 + 6 * q + 3 * q * q + q ** 3, rel=1e-14)
    assert m6 == pytest.approx(bruteforce_word(q, ['H'] * 6, {}), rel=1e-12)


@pytest.mark.parametrize("q", [0.3, 0.7, 0.9])
def test_theta_measure_matches_lattice(q):
    """The q-Pochhammer density rho(theta) [review (2.21)] reproduces the
    lattice moments: int rho = 1, int rho E^2 = 1, int rho E^4 = 2 + q."""
    th = (np.arange(2000) + 0.5) * np.pi / 2000
    w = np.pi / 2000
    rho = rho_theta(th, q)
    E = energy_theta(th, q)
    assert w * rho.sum() == pytest.approx(1.0, rel=1e-6)
    assert w * (rho * E ** 2).sum() == pytest.approx(1.0, rel=1e-6)
    assert w * (rho * E ** 4).sum() == pytest.approx(2.0 + q, rel=1e-6)


def test_band_edges():
    """Truncated-T spectrum lies inside [-2/sqrt(1-q), 2/sqrt(1-q)]."""
    for q in (0.2, 0.6, 0.9):
        s = Sector0(300, q)
        assert np.max(np.abs(s.E)) <= 2.0 / math.sqrt(1.0 - q) + 1e-9


# ---------------------------------------------------------------------------
# 3. the microscopic weight: q = e^{-2p^2/N} is the exact crossing average
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [4, 6])
def test_qeff_krawtchouk(p):
    """Two independent uniformly-random p-subsets I, J of {1..N}: the exact
    Majorana-string crossing average is E[(-1)^{|I cap J|}] =
    sum_m (-1)^m C(p,m) C(N-p,p-m) / C(N,p).  The chord weight q =
    exp(-2p^2/N) must be its N -> infinity limit at fixed lambda -- this
    validates the engine's microscopic input with no chords involved."""
    def qeff(N):
        return sum((-1) ** m * math.comb(p, m) * math.comb(N - p, p - m)
                   for m in range(p + 1)) / math.comb(N, p)
    errs = []
    for N in (160, 320, 640, 1280):
        q_chord = math.exp(-2.0 * p * p / N)
        errs.append(abs(qeff(N) - q_chord))
    # deviation is O(p^4/N^2): each N-doubling shrinks it ~4x (asymptotic
    # from below) -- this measured law is itself the honest statement of
    # HOW GOOD the chord weight is at finite N
    assert errs[-1] < 1e-3
    assert errs[-2] / errs[-1] > 3.5
    assert errs[-3] / errs[-2] > 3.3


# ---------------------------------------------------------------------------
# 3b. the CONTOUR walker (Contour.correlator) -- the path that produces
# chord_fold.json -- pinned three ways (added post-review: previously only
# apply_word_powers was enumerator-pinned)
# ---------------------------------------------------------------------------
def test_correlator_thermal_2pt_vs_closed_form():
    """Contour.correlator with asymmetric Euclidean segments vs the
    independent q-Pochhammer theta-integral (2.30)."""
    from chord import Contour, g2_closed_form
    q, delta, nmax = 0.6, 0.5, 60
    s = Sector0(nmax, q)
    for b1, b2 in [(0.7, 0.3), (1.4, 0.2)]:
        c = Contour(q, nmax, {'A': delta})
        amp = c.correlator([('A', b1), ('A', b2)])
        z = float(s.V[0] @ (np.exp(-(b1 + b2) * s.E) * s.V[0]))
        cf = g2_closed_form(b1, b2, q, delta, n_theta=240).real
        assert amp.real / z == pytest.approx(cf / z, rel=1e-8)
        assert abs(amp.imag) < 1e-12


def test_correlator_lorentzian_vs_eig():
    """Contour.correlator at a complex (Lorentzian) 2-pt contour vs the
    eigenbasis route g2_engine."""
    from chord import Contour, g2_engine
    q, delta, nmax = 0.5, 0.25, 80
    t = 1.3
    c = Contour(q, nmax, {'A': delta})
    # G(t) = <0| e^{+itT} q^{Dn} e^{-itT} |0>: segments u = -it (left of
    # first op along the trace) and u = +it
    amp = c.correlator([('A', -1j * t), ('A', +1j * t)])
    ref = g2_engine(np.array([t]), q, delta, nmax)[0]
    assert amp == pytest.approx(ref, rel=1e-9)


@pytest.mark.parametrize("labels,path_word", [
    (['A', 'B', 'A', 'B'], None),          # crossed
    (['A', 'B', 'B', 'A'], None),          # nested
    (['A', 'A', 'B', 'B'], None),          # adjacent
])
def test_correlator_4pt_vs_taylor_resummation(labels, path_word):
    """Contour.correlator on 4-pt orderings with four DISTINCT complex
    segments vs a Taylor-in-segments resummation of the enumerator-pinned
    apply_word_powers: sum over (k1..k4) of prod (-u_j)^{k_j}/k_j! times
    the word amplitude with H^{k_j} placed before operator j.  Pins the
    segment-placement convention, the sector-Delta selection, and the
    complex-u evolution path in one shot."""
    from math import factorial
    from chord import Contour
    q = 0.55
    deltas = {'A': 0.5, 'B': 0.25}
    ferm = {'A', 'B'}
    us = [0.22 + 0.11j, 0.17 - 0.09j, 0.25 + 0.05j, 0.14 + 0.13j]
    c = Contour(q, 20, deltas, ferm)
    amp = c.correlator(list(zip(labels, us)))
    # Taylor resummation over total H-order <= K
    K = 15
    total = 0.0 + 0.0j
    cw = Contour(q, 20, deltas, ferm)
    from itertools import product
    for ks in product(range(K + 1), repeat=4):
        if sum(ks) > K:
            continue
        coeff = 1.0 + 0.0j
        for u, k in zip(us, ks):
            coeff *= (-u) ** k / factorial(k)
        word = []
        for lab, k in zip(labels, ks):
            word += ['H'] * k + [lab]
        total += coeff * cw.apply_word_powers(word)
    assert amp == pytest.approx(total, rel=1e-9)


# ---------------------------------------------------------------------------
# 4. the two-point function: engine == theta-integral closed form (2.30)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q,delta", [(0.3, 0.5), (0.7, 0.25), (0.7, 1.0),
                                     (0.9, 0.5)])
def test_g2_engine_vs_closed_form(q, delta):
    """Thermal Euclidean 2-pt: transfer-matrix route vs the independent
    q-Pochhammer theta-integral [review (2.30)], several (beta1, beta2)."""
    s = Sector0(200, q)
    for b1, b2 in [(0.5, 0.5), (1.5, 0.5), (2.0, 1.0)]:
        beta = b1 + b2
        z = float(s.V[0] @ (np.exp(-beta * s.E) * s.V[0]))
        cf = g2_closed_form(b1, b2, q, delta, n_theta=240).real / z
        # engine via the eigenbasis (asymmetric ordering done directly)
        c = s.V.T @ np.eye(len(s.E))[:, 0]
        M = s.V.T @ ((q ** (delta * np.arange(len(s.E))))[:, None] * s.V)
        en = float((c * np.exp(-b1 * s.E)) @ M @ (np.exp(-b2 * s.E) * c)) / z
        assert en == pytest.approx(cf, rel=2e-5)


def test_g2_basics():
    """G(0) = 1 exactly; Delta = 0 is the identity channel: G = 1 at all t;
    G is real and |G| <= 1 at T = infinity (q^{Delta n} is a contraction)."""
    ts = np.linspace(0.0, 8.0, 30)
    g0 = g2_engine(ts, 0.6, 0.0, 80)
    assert np.allclose(g0, 1.0, atol=1e-12)
    g = g2_engine(ts, 0.6, 0.5, 80)
    assert g[0] == pytest.approx(1.0, abs=1e-12)
    assert np.max(np.abs(g.imag)) < 1e-10
    assert np.max(np.abs(g)) <= 1.0 + 1e-10


@pytest.mark.parametrize("q,delta", [(0.4, 0.5), (0.8, 0.25), (0.8, 2.0)])
def test_mu2_closed_form(q, delta):
    """Second spectral moment mu_2 = -G''(0) = 2(1 - q^Delta) exactly
    (hand-derived from (2.29); the engine's Lehmann lines must reproduce it
    to machine precision -- truncation is exact for moments)."""
    om, w = spectral_lines(q, delta, 40)
    mu2 = float(np.sum(np.real(w) * om ** 2))
    assert mu2 == pytest.approx(2.0 * (1.0 - q ** delta), rel=1e-10)
    mu0 = float(np.sum(np.real(w)))
    assert mu0 == pytest.approx(1.0, rel=1e-12)


def test_length_channel_is_dDelta():
    """dG/dDelta|_0 = -lambda <n(t)> exactly: numerical Delta-derivative of
    the engine 2-pt vs the direct n-expectation route."""
    q, nmax = 0.6, 120
    lam = -math.log(q)
    ts = np.linspace(0.0, 5.0, 12)
    eps = 1e-5
    gp = g2_engine(ts, q, eps, nmax).real
    gm = g2_engine(ts, q, 2 * eps, nmax).real
    # Richardson: dG/dD|_0 = (4(G(eps)-1)/... use two-point extrapolation
    d1 = (gp - 1.0) / eps
    d2 = (gm - 1.0) / (2 * eps)
    dG = 2 * d1 - d2
    nt = n_expectation(ts, q, nmax)
    assert np.allclose(dG, -lam * nt, rtol=2e-4, atol=1e-8)
    assert nt[0] == pytest.approx(0.0, abs=1e-12)
    assert np.all(np.diff(nt) > -1e-9)      # monotone growth


def test_length_lines_match_curve():
    """Lehmann lines of the length channel resum to -lambda <n(t)>."""
    q, nmax = 0.5, 100
    lam = -math.log(q)
    om, w = length_channel_lines(q, nmax)
    ts = np.linspace(0.0, 4.0, 9)
    resum = np.array([np.sum(w * np.exp(1j * om * t)) for t in ts]).real
    nt = n_expectation(ts, q, nmax)
    assert np.allclose(resum, -lam * nt, rtol=1e-9, atol=1e-10)


# ---------------------------------------------------------------------------
# 5. the large-q (lambda -> 0) sech anchor: parameter-free bridge
# ---------------------------------------------------------------------------
def test_sech_anchor_smalllambda():
    """At small lambda the exact chord 2-pt must approach the large-p SYK
    closed form G = sech(Jc t)^{2 Delta} with the DERIVED bridge
    Jc = sqrt(lambda) e^{lambda/8} in chord units (module docstring; the
    e^{lambda/8} is the double-scaled sigma_H correction).  Tolerance is the
    expected O(lambda, Delta/p-corrections) scale, measured tight."""
    lam, delta = 0.05, 0.5
    q = math.exp(-lam)
    Jc = math.sqrt(lam) * math.exp(lam / 8.0)
    ts = np.linspace(0.0, 2.0 / Jc, 15)
    g = g2_engine(ts, q, delta, 700).real
    sech = (1.0 / np.cosh(Jc * ts)) ** (2 * delta)
    dev = np.max(np.abs(g - sech))
    assert dev < 3 * lam                    # O(lambda) agreement
    # and the deviation halves when lambda halves (first order in lambda)
    lam2 = lam / 2
    q2, Jc2 = math.exp(-lam2), math.sqrt(lam2) * math.exp(lam2 / 8.0)
    ts2 = np.linspace(0.0, 2.0 / Jc2, 15)
    g2v = g2_engine(ts2, q2, delta, 1000).real
    sech2 = (1.0 / np.cosh(Jc2 * ts2)) ** (2 * delta)
    dev2 = np.max(np.abs(g2v - sech2))
    assert dev2 < 0.65 * dev


# ---------------------------------------------------------------------------
# 6. the finite-N ED anchor (fold_bench.json): chord vs exact diagonalization
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 7. the committed rung-18 scan: the discrete-vs-continuum verdict
# ---------------------------------------------------------------------------
def test_structure_detector_catches_planted_features():
    """The detector itself is falsifiable: structure_report_smooth must
    flag a planted interior peak and a planted node (guards the verdict
    test below against a broken-detector false-green, a hole demonstrated
    by the adversarial review)."""
    from chord import structure_report_smooth
    x = np.linspace(0.0, 10.0, 500)
    smooth = np.exp(-x)
    assert structure_report_smooth(x, smooth)["n_local_maxima"] == 0
    peaked = smooth.copy()
    peaked[250] += 5.0
    assert structure_report_smooth(x, peaked)["n_local_maxima"] >= 1
    noded = smooth - 0.1
    assert structure_report_smooth(x, noded)["n_sign_changes"] >= 1


def test_committed_spectra_verdict():
    """chord_spectra.json, the rung-18 first deliverable: at N = infinity,
    T_B = infinity, every channel (O_Delta and the Delta->0 fermion-
    bilinear/length channel) is a STRUCTURELESS CONTINUUM at every lambda
    down to 0.05, converging to the smooth semiclassical sech^{2Delta}
    transform linearly in lambda.

    Post-review hardening: the structure diagnostics are RECOMPUTED here
    from the committed S_fine/om_fine and S_coarse/om_coarse arrays (fine
    window at the semiclassical scale AND the full-support coarse curve),
    and the sech deviation is recomputed from the stored arrays -- the
    stored verdict scalars are cross-checked against them, so a corrupted
    or wrongly-diagnosed committed curve fails the suite."""
    from chord import structure_report_smooth
    with open(HERE / "chord_spectra.json") as fh:
        data = json.load(fh)
    pts = sorted(data["points"], key=lambda p: -p["lam"])
    assert pts[-1]["lam"] <= 0.05
    sech_devs = []
    for pt in pts:
        om_coarse = np.array(pt["om_coarse"])
        for key, ch in pt["channels"].items():
            om_f = np.array(ch["om_fine"])
            S_f = np.array(ch["S_fine"])
            rep = structure_report_smooth(om_f, S_f)
            assert rep["n_local_maxima"] == 0, (pt["lam"], key)
            assert rep["n_local_maxima"] == ch["n_local_maxima"]
            assert rep["n_sign_changes"] == ch["n_sign_changes"] == 0
            # full-support coarse curve: no interior structure hiding
            # beyond the fine window (length channel: skip its om->0
            # growth region, as in the campaign)
            S_c = np.array(ch["S_coarse"])
            sel = om_coarse > (0.5 * pt["Jc_bridge"] if key == "length"
                               else 0.0)
            # ignore float-rounding wiggles in the deep tail
            floor = 1e-8 * np.abs(S_c).max()
            sel &= np.abs(S_c) > floor
            rep_c = structure_report_smooth(om_coarse[sel], S_c[sel])
            assert rep_c["n_local_maxima"] == 0, (pt["lam"], key, "coarse")
            assert ch["quadrature_drift"] < 1e-3
            if ch["mu0_closed"] is not None:
                # 5e-3: omega-grid trapezoid floor on the narrowest
                # (lam = 0.05) spectral functions, measured 3.1e-3
                assert abs(ch["mu0"] - ch["mu0_closed"]) < 5e-3
            assert abs(ch["mu2"] - ch["mu2_closed"]) < 2e-3 * max(
                ch["mu2_closed"], 1.0)
            if ch["sech_dev"] is not None:
                ssc = np.array(ch["S_semiclassical_fine"])
                dev = float(np.max(np.abs(S_f - ssc)) / max(ssc.max(), 1e-300))
                assert dev == pytest.approx(ch["sech_dev"], rel=1e-3)
        sech_devs.append(pt["channels"]["delta=0.5"]["sech_dev"])
    # linear-in-lambda convergence to the semiclassical continuum
    # (the quoted numbers are the Delta = 1/2 channel; the worst channel,
    # Delta = 1, obeys the same law at ~2x the coefficient)
    assert all(sech_devs[i + 1] < sech_devs[i] for i in range(len(pts) - 1))
    assert sech_devs[-1] < 0.01
    for pt, dev in zip(pts, sech_devs):
        assert dev < 0.25 * pt["lam"]
    for pt in pts:
        d1 = pt["channels"]["delta=1.0"]["sech_dev"]
        assert d1 < 0.45 * pt["lam"]


def test_committed_length_growth():
    """The chord-number (length) growth curve <n(t)> approaches
    -(2/lambda) log sech(Jc t) linearly in lambda -- the de Sitter length
    growth, exact in lambda, converging to its semiclassical form."""
    with open(HERE / "chord_spectra.json") as fh:
        data = json.load(fh)
    devs = []
    for pt in sorted(data["points"], key=lambda p: -p["lam"]):
        g = pt["length_growth"]
        nt = np.array(g["n_of_t"])
        sc = np.array(g["n_of_t_semiclassical"])
        sel = sc > 0.5
        devs.append((pt["lam"],
                     float(np.max(np.abs(nt[sel] - sc[sel]) / sc[sel]))))
    for lam, dev in devs:
        assert dev < 0.2 * lam + 0.005
    assert all(devs[i + 1][1] < devs[i][1] for i in range(len(devs) - 1))


def test_ed_anchor_thermal_g2():
    """Chord thermal 2-pt at (p=6, N=14, betaJ=2) vs the disorder-averaged
    ED values committed in syk-self-averaging/fold_bench.json, through the
    EXACT finite-N units bridge t_chord = t_repo * sigma_H(N,p).  The chord
    answer is the strict N = infinity, fixed-lambda limit, so the gap is the
    genuine finite-N correction -- measured here at the few-percent level
    (asserted < 5%), far inside the ~30-55% floor of the large-q closed
    forms at the same sizes (ANTISCRAMBLING.md)."""
    path = HERE.parent / "syk-self-averaging" / "fold_bench.json"
    with open(path) as fh:
        data = json.load(fh)
    pt = next(p for p in data["points"]
              if p["N"] == 14 and p["p"] == 6 and p["betaJ"] == 2.0)
    N, p, betaJ = pt["N"], pt["p"], pt["betaJ"]
    lam = lambda_chord_of(N, p)
    q = math.exp(-lam)
    sH = sigma_H(N, p)
    beta_phys = betaJ / (math.sqrt(2.0) * p)
    beta_c = beta_phys * sH
    s = Sector0(80, q)
    c = s.V[0, :]
    z = float(c @ (np.exp(-beta_c * s.E) * c))
    delta = 1.0 / p
    M = s.V.T @ ((q ** (delta * np.arange(81)))[:, None] * s.V)
    for th, ed, err in zip(pt["th_g2"], pt["g2_mean"], pt["g2_err"]):
        tau = beta_c * th / (2.0 * np.pi)
        ch = float((c * np.exp(-(beta_c - tau) * s.E)) @ M
                   @ (np.exp(-tau * s.E) * c)) / z
        assert abs(ch - ed) < 0.05 * abs(ed) + 3 * err, \
            f"theta={th}: chord {ch:.4f} vs ED {ed:.4f} +- {err:.4f}"
