"""Gates for the Dirac/complex-SYK chord layer (dirac_chord.py).

Four legs, mirroring the structure of test_chord.py.

1. THE CONVENTIONS BRIDGE.  q_Dirac = e^{-p^2/Nc} against dictionary.py's
   converters, which is e^{-lambda_paper}, and is the Majorana-literature
   q at N = 2 Nc.  Plus the m2 closed forms against seeded ED means, and
   against exact_wick_dirac's independent route.

2. THE MICROSCOPIC LEG, the Dirac Krawtchouk analog.  The exact finite-Nc
   crossed/uncrossed monomial trace average tends to e^{-p^2/Nc}.  Seeded
   compressed-register sampling; no chords involved.

3. GATE G1.  Transfer-matrix moments <0|T^{2k}|0> against
   disorder-averaged Dirac ED Tr H^{2k}/dim ratios at matched q(Nc),
   converging as Nc grows.  Measured with c_symmetric, seeds 1000+, 6
   draws, where deviation = ED - chord:

       r4:  Nc=8 +0.198,  Nc=10 +0.050,  Nc=12 +0.004 (4 draws, slow)
       r6:  Nc=8 +1.273,  Nc=10 +0.339,  Nc=12 +0.041

   The tolerances below are these deviations with about 50% headroom.
   They bound the finite-Nc convergence claim.  They are not machine
   precision.

4. GATE G2.  The M_1 channel is conserved, because M_1 = -(p/2) i H
   exactly.  Its correlator is time-independent, and its Lehmann weight
   off omega = 0 is zero.

GATE G3 IS NOT RUN.  That would be the CP-resolved M_n tower against the
moments_pipeline data.  The tower's chord rules are unpublished; see the
dirac_chord docstring.  The scope test pins that claim so it cannot
silently rot.
"""

import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from chord import bruteforce_word, g2_closed_form
from dictionary import lambda_paper, lambda_chord, q_parameter
from dirac_chord import (MN_TOWER_CHORD_RULES_PUBLISHED,
                         charge_sector_asymmetry, delta_charged,
                         dirac_contour, dirac_sector0, g2_charged,
                         hamiltonian_moment_ratios, lambda_dirac,
                         m1_channel_correlator, m1_offdiagonal_weight,
                         m2_mean_dirac, neutral_bilinear_moments,
                         neutral_bilinear_spectrum, probe_is_fermionic,
                         q_dirac, q_dirac_charge, sigma_H_dirac)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "syk-self-averaging"))

from dirac import (apply_monomial, dirac_hamiltonian,      # noqa: E402
                   draw_couplings)

P = 4  # the repo's fixed-p trajectory


# ---------------------------------------------------------------------------
# 1. conventions bridge
# ---------------------------------------------------------------------------
def test_q_dirac_vs_dictionary():
    """q_Dirac = e^{-lambda_paper(Nc,p)} exactly, and equals the chord-
    literature q at N = 2 Nc Majoranas: the Dirac side has NO factor-2
    landmine (dictionary.py correction #4 is Majorana-specific)."""
    for Nc, p in [(8, 4), (10, 4), (12, 4), (20, 6)]:
        q = q_dirac(Nc, p)
        assert q == pytest.approx(math.exp(-lambda_paper(Nc, p)), rel=1e-14)
        assert q == pytest.approx(q_parameter(lambda_chord(2 * Nc, p)),
                                  rel=1e-14)
        assert lambda_dirac(Nc, p) == pytest.approx(p * p / Nc, rel=1e-14)


def test_q_charge_sector_limits():
    """q(Qs) [BNR (5.17)]: Qs = 0 gives q_dirac; |Qs| -> 1/2 freezes the
    sector (q -> 0); symmetric in Qs."""
    Nc, p = 12, 4
    assert q_dirac_charge(Nc, p, 0.0) == pytest.approx(q_dirac(Nc, p))
    assert q_dirac_charge(Nc, p, 0.3) == pytest.approx(
        q_dirac_charge(Nc, p, -0.3), rel=1e-14)
    assert q_dirac_charge(Nc, p, 0.49) < 1e-10
    with pytest.raises(ValueError):
        q_dirac_charge(Nc, p, 0.5)


def test_m2_closed_form_vs_exact_wick_dirac():
    """m2_mean_dirac('c_symmetric') equals exact_wick_dirac.m2_mean_exact
    (numeric trace route over the pair operators; shares no algebra)."""
    from exact_wick_dirac import m2_mean_exact
    for Nc in (6, 8):
        assert m2_mean_dirac(Nc, P, "c_symmetric") == pytest.approx(
            m2_mean_exact(Nc, P), rel=1e-10)


@lru_cache(maxsize=None)
def _ed_moments(Nc: int, ensemble: str, n_seeds: int):
    """Disorder-mean (m2, m4, m6) of Tr H^{2k}/dim and the per-seed m4/m2^2
    list, seeded 1000+j.  Cached across tests."""
    rows = []
    for j in range(n_seeds):
        rng = np.random.default_rng(1000 + j)
        H = dirac_hamiltonian(Nc, P,
                              couplings=draw_couplings(Nc, P, rng, ensemble))
        dim = 1 << Nc
        m2 = float(np.vdot(H, H).real) / dim
        H2 = H @ H
        m4 = float(np.vdot(H2, H2).real) / dim
        H3 = H @ H2
        m6 = float(np.vdot(H3, H3).real) / dim
        rows.append((m2, m4, m6))
    arr = np.array(rows)
    mean = arr.mean(axis=0)
    r4_per_seed = arr[:, 1] / arr[:, 0] ** 2
    return tuple(mean), tuple(r4_per_seed)


@pytest.mark.parametrize("ensemble,rel_tol", [("generic", 0.08),
                                              ("c_symmetric", 0.20)])
def test_m2_closed_forms_vs_ed(ensemble, rel_tol):
    """E[Tr H^2/dim] closed forms vs 6-draw seeded ED means at Nc = 8.
    The tolerance covers disorder fluctuation of the mean (per-instance
    rms(m2)/m2 ~ 10% for c_symmetric at this size); a z-score guard with
    the measured SEM backs it up."""
    (m2, _, _), _ = _ed_moments(8, ensemble, 6)
    closed = m2_mean_dirac(8, P, ensemble)
    assert abs(m2 - closed) / closed < rel_tol
    assert sigma_H_dirac(8, P, ensemble) == pytest.approx(math.sqrt(closed))


# ---------------------------------------------------------------------------
# 2. microscopic leg: the Dirac crossing average (no chords involved)
# ---------------------------------------------------------------------------
def _word_trace(ops, modes):
    """tr(prod ops)/dim on the compressed register of involved modes.
    Sound because every mode in a trace-nonvanishing word hosts an even
    number of operators, so spectator-mode Jordan-Wigner strings cancel."""
    relab = {m: i for i, m in enumerate(modes)}
    nc = len(modes)
    dim = 1 << nc
    states0 = np.arange(dim, dtype=np.int64)
    cur = states0.copy()
    sign = np.ones(dim)
    valid = np.ones(dim, dtype=bool)
    for (A, B) in reversed(ops):
        Ar = tuple(sorted(relab[a] for a in A))
        Br = tuple(sorted(relab[b] for b in B))
        v, rows, s = apply_monomial(nc, Ar, Br, cur)
        valid &= v
        sign = sign * s
        cur = np.where(v, rows, cur)
    good = valid & (cur == states0)
    return float(sign[good].sum()) / dim


def _crossing_average(Nc, k, n_samp, seed):
    """E[tr(O_a O_b O_abar O_bbar)] / E[tr(O_a O_abar O_b O_bbar)] over
    independent uniform k-subsets (A1,B1,A2,B2): crossed vs uncrossed
    H-chord 4-words, the finite-Nc object whose limit BNR (4.3) says is
    q = e^{-4k^2/Nc}.  Disjoint-support samples contribute identically
    to both orders (even operators commute), so the ratio's deviation
    from 1 is driven purely by index-coincidence events."""
    rng = np.random.default_rng(seed)
    num = den = 0.0
    for _ in range(n_samp):
        A1, B1, A2, B2 = (tuple(sorted(rng.choice(Nc, k, replace=False)))
                          for _ in range(4))
        modes = sorted(set(A1) | set(B1) | set(A2) | set(B2))
        num += _word_trace([(A1, B1), (A2, B2), (B1, A1), (B2, A2)], modes)
        den += _word_trace([(A1, B1), (B1, A1), (A2, B2), (B2, A2)], modes)
    return num / den


def test_microscopic_crossing_average():
    """Measured (seed 7, 2500 samples, k = 2): the crossing average vs
    e^{-4k^2/Nc} deviates by a few 1e-2 and shrinks with Nc.  Tolerances
    = expected O(k^4/Nc^2) systematic + sampling noise, with headroom."""
    k = P // 2
    devs = {}
    for Nc, tol in [(16, 0.05), (32, 0.03), (64, 0.025)]:
        q = math.exp(-4.0 * k * k / Nc)
        r = _crossing_average(Nc, k, 2500, seed=7)
        devs[Nc] = abs(r - q)
        assert devs[Nc] < tol, f"Nc={Nc}: crossing {r} vs q {q}"
    assert devs[64] < devs[16]


# ---------------------------------------------------------------------------
# 3. gate G1: transfer-matrix H-moments vs ED at matched q
# ---------------------------------------------------------------------------
def test_transfer_matrix_moment_polynomials():
    """<0|T^{2k}|0> from exact T-powers equals the chord polynomials
    2 + q and 5 + 6q + 3q^2 + q^3, and equals the brute-force diagram
    enumerator (the object BNR (4.3) equates the Dirac moments to)."""
    for q in (0.1353, 0.2636, 0.7):
        r = hamiltonian_moment_ratios(q)
        assert r[2] == pytest.approx(1.0, abs=1e-12)
        assert r[4] == pytest.approx(2.0 + q, rel=1e-12)
        assert r[6] == pytest.approx(5 + 6 * q + 3 * q ** 2 + q ** 3,
                                     rel=1e-12)
        assert r[6] == pytest.approx(bruteforce_word(q, ['H'] * 6, {}),
                                     rel=1e-12)


def _g1_deviations(Nc, ensemble, n_seeds):
    (m2, m4, m6), r4_seeds = _ed_moments(Nc, ensemble, n_seeds)
    q = q_dirac(Nc, P)
    r = hamiltonian_moment_ratios(q)
    dev4 = m4 / m2 ** 2 - r[4]
    dev6 = m6 / m2 ** 3 - r[6]
    sem4 = float(np.std(r4_seeds, ddof=1) / math.sqrt(len(r4_seeds)))
    return dev4, dev6, sem4


def test_g1_moments_converge_csym():
    """G1 (fast half): c_symmetric ED ratios at Nc = 8, 10 vs chord values
    at q(Nc).  Asserts the measured finite-Nc deviations (docstring table)
    within ~50% headroom AND that they shrink with Nc."""
    dev4_8, dev6_8, _ = _g1_deviations(8, "c_symmetric", 6)
    dev4_10, dev6_10, _ = _g1_deviations(10, "c_symmetric", 6)
    assert abs(dev4_10) < 0.08
    assert abs(dev6_10) < 0.50
    assert abs(dev4_10) < abs(dev4_8)
    assert abs(dev6_10) < abs(dev6_8)


def test_g1_generic_shares_limit_but_slower():
    """The generic ensemble has the same N = infinity chord limit but
    converges more slowly (index coincidences inside a term are frequent,
    O(k^2/Nc), at ED sizes).  Measured at Nc = 10: |dev4| = 0.23 generic
    vs 0.05 c_symmetric."""
    dev4_gen, _, _ = _g1_deviations(10, "generic", 6)
    dev4_sym, _, _ = _g1_deviations(10, "c_symmetric", 6)
    assert abs(dev4_gen) < 0.35            # still heading to the limit
    assert abs(dev4_sym) < abs(dev4_gen)   # the CP ensemble converges faster


@pytest.mark.slow
def test_g1_moments_nc12_slow():
    """G1 (slow half): Nc = 12, 4 draws.  Measured dev4 = +0.004 (SEM
    0.006), dev6 = +0.041.  Asserts the near-closure and the z-score of
    the r4 deviation."""
    dev4_10, dev6_10, _ = _g1_deviations(10, "c_symmetric", 6)
    dev4, dev6, sem4 = _g1_deviations(12, "c_symmetric", 4)
    assert abs(dev4) < 0.03
    assert abs(dev6) < 0.15
    assert abs(dev4) < abs(dev4_10)
    assert abs(dev6) < abs(dev6_10)
    assert abs(dev4) < 4.0 * sem4          # consistent with closure at SEM


# ---------------------------------------------------------------------------
# 4. gate G2: the conserved M_1 channel
# ---------------------------------------------------------------------------
def test_g2_m1_channel_conserved():
    """The chord image of M_1 = -(p/2) i H is a transfer-matrix insertion.
    Its correlator is time-independent to machine precision (computed by
    explicit evolution, commutativity never invoked) at T_B = infinity
    and at finite beta, and its Lehmann weight off omega = 0 vanishes."""
    Nc = 12
    ts = np.linspace(0.0, 6.0, 13)
    for beta in (0.0, 1.5):
        g = m1_channel_correlator(ts, Nc, P, nmax=48, beta=beta)
        assert np.abs(g - g[0]).max() < 1e-10
        assert np.abs(g.imag).max() < 1e-10
        # G(0) = <T^2>_beta > 0
        assert g[0].real > 0.0
    assert m1_offdiagonal_weight(Nc, P, nmax=48) < 1e-10
    assert m1_offdiagonal_weight(Nc, P, nmax=48, beta=1.5) < 1e-10


def test_g2_ed_anchor_operator_identity():
    """The SYK-side anchor of G2: M_1 + (p/2) i H = 0 exactly on a
    c_symmetric instance (so conservation of the M_1 channel is an
    operator identity, not an approximation).  One cheap instance."""
    from dirac import annihilators, mn_tower
    rng = np.random.default_rng(5)
    H = dirac_hamiltonian(6, P,
                          couplings=draw_couplings(6, P, rng, "c_symmetric"))
    M1 = mn_tower(H, annihilators(6), 1)[1]
    assert np.abs(M1 + (P / 2) * 1j * H).max() < 1e-10


# ---------------------------------------------------------------------------
# charged-probe machinery: reductions and cross-checks
# ---------------------------------------------------------------------------
def test_delta_and_fermion_parity():
    assert delta_charged(1, 1, 4) == pytest.approx(0.5)
    assert delta_charged(2, 2, 4) == pytest.approx(1.0)
    assert delta_charged(1, 0, 4) == pytest.approx(0.25)
    assert not probe_is_fermionic(1, 1)
    assert probe_is_fermionic(1, 0)
    with pytest.raises(ValueError):
        delta_charged(0, 0, 4)


def test_charge_asymmetry_reductions():
    """B1 = B2 = 0 at Qs = 0 for any (a, b); at a = b the two sides remain
    unequal only through Qs (and swap under Qs -> -Qs with a <-> b:
    charge conjugation)."""
    Nc = 12
    assert charge_sector_asymmetry(Nc, P, 2, 1, 0.0) == (1.0, 1.0)
    s1, s2 = charge_sector_asymmetry(Nc, P, 2, 1, 0.2)
    t1, t2 = charge_sector_asymmetry(Nc, P, 1, 2, -0.2)
    assert s1 == pytest.approx(t1, rel=1e-12)
    assert s2 == pytest.approx(t2, rel=1e-12)
    assert s1 != pytest.approx(1.0)


def test_g2_charged_reduces_to_majorana_form():
    """At Qs = 0 the fixed-charge 2-pt amplitude [BNR (7.23)] equals the
    Majorana theta-integral [review (2.30)] at (q_dirac, (a+b)/p)."""
    Nc = 12
    for (a, b), (b1, b2) in [((1, 1), (0.7, 0.3)), ((2, 1), (0.5, 0.5))]:
        got = g2_charged(b1, b2, Nc, P, a, b)
        want = g2_closed_form(b1, b2, q_dirac(Nc, P),
                              delta_charged(a, b, P))
        assert got.real == pytest.approx(want.real, rel=1e-10)


def test_dirac_contour_matches_theta_integral():
    """The transfer-matrix contour walker at the Dirac parameters (via
    dirac_contour) agrees with the independent q-Pochhammer route for a
    thermal 2-pt of the neutral bilinear -- the Dirac-layer version of
    test_chord.test_correlator_thermal_2pt_vs_closed_form."""
    Nc, nmax = 12, 60
    c = dirac_contour(nmax, Nc, P, {'A': (1, 1)})
    s = dirac_sector0(nmax, Nc, P)
    b1, b2 = 0.6, 0.4
    amp = c.correlator([('A', b1), ('A', b2)])
    z = float(s.V[0] @ (np.exp(-(b1 + b2) * s.E) * s.V[0]))
    cf = g2_charged(b1, b2, Nc, P, 1, 1, n_theta=240)
    assert amp.real / z == pytest.approx(cf.real / z, rel=1e-8)
    assert abs(amp.imag) < 1e-12


def test_neutral_bilinear_moment_closures():
    """The independent neutral bilinear channel: int S = 1 and
    int omega^2 S = 2(1 - q^{2/p}) by quadrature of the exact smooth
    spectral function at the Dirac q (the review's tested closures,
    rewired through the Dirac layer)."""
    Nc = 12
    q = q_dirac(Nc, P)
    band = 2.0 / math.sqrt(1.0 - q)
    om = np.linspace(1e-4, 2.05 * band, 400)
    S = neutral_bilinear_spectrum(om, Nc, P, n_theta=800)
    mu0 = 2.0 * np.trapezoid(S, om)
    mu2 = 2.0 * np.trapezoid(S * om ** 2, om)
    ref = neutral_bilinear_moments(Nc, P)
    assert mu0 == pytest.approx(ref["mu0"], rel=2e-2)
    assert mu2 == pytest.approx(ref["mu2"], rel=2e-2)


# ---------------------------------------------------------------------------
# gate G3 scope: pinned, not fabricated
# ---------------------------------------------------------------------------
def test_g3_scope_is_pinned():
    """The M_n tower (n >= 2) chord rules are unpublished -- the module
    must keep saying so (and naming its sources) until someone derives
    them.  If that constant ever flips, this test demands the derivation
    exist (and a real G3 gate replace it)."""
    assert MN_TOWER_CHORD_RULES_PUBLISHED is False
    doc = __import__("dirac_chord").__doc__
    for needle in ("2006.13983", "2407.09396", "2302.11882", "2512.07715",
                   "NOT established"):
        assert needle in doc
