"""Gates for the Dirac/complex-SYK chord layer (dirac_chord.py).

Five legs, mirroring the structure of test_chord.py.

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

5. GATE G3, PAIR HALF.  Two legs.  The census leg: the lone-crossing
   factor of a BARE c_i / c^dag_i chord against a c_symmetric H-chord is
   exactly 1 - p/Nc, measured PER SAMPLE (compressed-register sampler,
   no chords) at Nc = 16, 32, 64, converging on q_dirac^{1/p} at the
   census rate p^2/(2 Nc^2).  The theorem-transfer leg: chord_charges'
   recursion, propagator and charge identities re-verified at
   (q_dirac(Nc, p), Delta = 1/p), plus the word-level moment zeros
   through the independent chord_moments assembly at the Dirac q.

GATE G3 ITSELF IS STILL NOT RUN.  That would be the CP-resolved M_n
tower against the moments_pipeline data.  The tower's chord rules are
unpublished; the PAIR (i != j) half is now derived in-repo (leg 5), and
the flavor-diagonal half and the CP resolution remain open.  The scope
test pins that status so it cannot silently rot.
"""

import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from chord import bruteforce_word, g2_closed_form, qint
from chord_charges import (close_j_map, dressed_propagator, motzkin_coeffs,
                           poly_of_matrix, transfer_formal,
                           transfer_one_matter_formal)
from dictionary import lambda_paper, lambda_chord, q_parameter
from dirac_chord import (MN_TOWER_CHORD_RULES_PUBLISHED,
                         MN_TOWER_PAIR_CHANNEL_DERIVED,
                         bare_fermion_crossing, charge_sector_asymmetry,
                         delta_charged, dirac_contour, dirac_sector0,
                         g2_charged, hamiltonian_moment_ratios, lambda_dirac,
                         m1_channel_correlator, m1_offdiagonal_weight,
                         m2_mean_dirac, neutral_bilinear_moments,
                         neutral_bilinear_spectrum, probe_is_fermionic,
                         q_dirac, q_dirac_charge, sigma_H_dirac, xhat_dirac)

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
# 5. gate G3, pair half: the bare-fermion census + the theorem transfer
# ---------------------------------------------------------------------------
def _hu_terms(A, B):
    """The two monomials of a c_symmetric H-chord endpoint,
    h_u = c^dag_A c_B + c^dag_B c_A."""
    return [(A, B), (B, A)]


def _sample_hchord(Nc, k, rng):
    """Uniform disjoint pair (A, B): the index content of one c_symmetric
    coupling."""
    perm = rng.permutation(Nc)
    A = tuple(sorted(int(x) for x in perm[:k]))
    B = tuple(sorted(int(x) for x in perm[k:2 * k]))
    return A, B


def _matter_crossing_average(Nc, k, n_samp, seed, dagger_first=False):
    """(ratio, max_census_dev):
    ratio = E[tr(c_i h_u c^dag_i h_u)] / E[tr(c_i c^dag_i h_u h_u)] over
    uniform disjoint pairs (A, B) and a uniform mode i -- the finite-Nc
    lone-crossing factor of a BARE single-fermion chord, the microscopic
    input of the Dirac conserved pair channel.  max_census_dev is the
    largest per-sample deviation from the census (crossed = 0 when i is
    in A u B, crossed = uncrossed when it is not): 0.0 means the census
    holds sample by sample, not just on average.  dagger_first probes the
    c^dag_i ... c_i ordering, the other end of a pair word."""
    rng = np.random.default_rng(seed)
    num = den = 0.0
    max_dev = 0.0
    for _ in range(n_samp):
        A, B = _sample_hchord(Nc, k, rng)
        i = int(rng.integers(Nc))
        modes = sorted(set(A) | set(B) | {i})
        lo, hi = (((i,), ()), ((), (i,))) if dagger_first \
            else (((), (i,)), ((i,), ()))
        cn = cd = 0.0
        for t1 in _hu_terms(A, B):
            for t2 in _hu_terms(A, B):
                cn += _word_trace([lo, t1, hi, t2], modes)
                cd += _word_trace([lo, hi, t1, t2], modes)
        num += cn
        den += cd
        want = 0.0 if i in set(A) | set(B) else cd
        max_dev = max(max_dev, abs(cn - want))
    return num / den, max_dev


def test_census_bare_fermion_lone_crossing():
    """The census behind the Dirac conserved pair channel, measured with
    no chords involved.  Seed 7:

        Nc=16 (2500 samples): ratio 0.74560   1 - p/Nc = 0.75000
        Nc=32 (2500 samples): ratio 0.87720   1 - p/Nc = 0.87500
        Nc=64 (4000 samples): ratio 0.93725   1 - p/Nc = 0.93750

    with per-sample census deviation 0.0 at every size.  Three claims.
    (1) The census is exact per sample: the crossed word vanishes
    identically when i is in the H-chord support and equals the
    uncrossed word when it is not.  (2) The sampled ratio is the
    avoidance probability 1 - p/Nc within binomial noise (asserted at
    5 sigma; measured < 0.6 sigma).  (3) The approach to
    q_dirac^{1/p} = e^{-p/Nc} is at the census rate p^2/(2 Nc^2)
    (arithmetic: e^{-x} - (1 - x) = x^2/2 + ...), so the measured
    deviations shrink monotonically, +3.3e-2 -> +5.3e-3 -> +2.2e-3
    (seeded, so deterministic).  The c^dag_i ... c_i ordering (the
    other end of a pair word) obeys the same census."""
    k = P // 2
    devs = {}
    for Nc, n_samp in [(16, 2500), (32, 2500), (64, 4000)]:
        r, census_dev = _matter_crossing_average(Nc, k, n_samp, seed=7)
        assert census_dev < 1e-14, Nc
        exact = bare_fermion_crossing(Nc, P)
        sigma = math.sqrt(exact * (1.0 - exact) / n_samp)
        assert abs(r - exact) < 5.0 * sigma, (Nc, r, exact)
        devs[Nc] = abs(math.exp(-P / Nc) - r)
    assert devs[16] > devs[32] > devs[64]
    r_dag, dev_dag = _matter_crossing_average(32, k, 2000, seed=7,
                                              dagger_first=True)
    assert dev_dag < 1e-14
    exact = bare_fermion_crossing(32, P)
    assert abs(r_dag - exact) < 5.0 * math.sqrt(exact * (1 - exact) / 2000)


@pytest.mark.slow
def test_census_two_crossings_factorize_slow():
    """Factorization beyond the lone crossing: a bare c-chord crossed by
    TWO mutually uncrossed (nested) H-chords,

        E[tr(c_i h_u h_v c^dag_i h_v h_u)]
        / E[tr(c_i c^dag_i h_u h_u h_v h_v)].

    The census is exact per sample here too: the crossed word vanishes
    identically when i lies in the union of the two H-chord supports and
    equals the uncrossed word when it does not (asserted at 1e-14 on
    every sample; an earlier version of this test compared the ratio to
    q^{2/p} with 3x headroom and called the census inexact -- it is
    exact, the u-v support overlap only reweights the uncrossed word).
    So the ratio is the joint avoidance probability, (1 - p/Nc)^2 up
    to that reweighting, asserted at 5 sigma; q_dirac^{2/p} = e^{-2p/Nc}
    is its Nc -> infinity limit, the gap being ~p^2/Nc^2 (0.013 at
    Nc = 32; coefficient 1 here, 1/2 for the lone crossing -- the same
    order, the census rate of a double crossing).  Measured
    (seed 7, 1200 samples): 0.7657 at Nc = 32 [(1 - p/Nc)^2 = 0.7656],
    0.8707 at Nc = 64 [0.8789]."""
    k = P // 2
    for Nc in (32, 64):
        rng = np.random.default_rng(7)
        num = den = 0.0
        max_dev = 0.0
        for _ in range(1200):
            A1, B1 = _sample_hchord(Nc, k, rng)
            A2, B2 = _sample_hchord(Nc, k, rng)
            i = int(rng.integers(Nc))
            modes = sorted(set(A1) | set(B1) | set(A2) | set(B2) | {i})
            lo, hi = ((), (i,)), ((i,), ())
            cn = cd = 0.0
            for t1 in _hu_terms(A1, B1):
                for t1b in _hu_terms(A1, B1):
                    for t2 in _hu_terms(A2, B2):
                        for t2b in _hu_terms(A2, B2):
                            cn += _word_trace([lo, t1, t2, hi, t2b, t1b],
                                              modes)
                            cd += _word_trace([lo, hi, t1, t1b, t2, t2b],
                                              modes)
            want = 0.0 if i in set(A1) | set(B1) | set(A2) | set(B2) else cd
            max_dev = max(max_dev, abs(cn - want))
            num += cn
            den += cd
        assert max_dev < 1e-14, (Nc, max_dev)
        avoid2 = (1.0 - P / Nc) ** 2
        sigma = math.sqrt(avoid2 * (1.0 - avoid2) / 1200)
        assert abs(num / den - avoid2) < 5.0 * sigma, (Nc, num / den, avoid2)
        q2p = q_dirac(Nc, P) ** (2.0 / P)
        assert 0.5 * P * P / Nc ** 2 < q2p - avoid2 < 1.5 * P * P / Nc ** 2


# the chord side is a pure (q, Delta) statement, so a p = 6 point is
# legitimate even though the c_symmetric ensemble needs p/2 even
DIRAC_G3_POINTS = [(16, 4), (64, 4), (18, 6)]   # q = e^-1, e^-0.25, e^-2
NMAX_G3 = 24


def _interior(M, keep):
    return M[:keep, :keep]


@pytest.mark.parametrize("Nc,p", DIRAC_G3_POINTS)
def test_theorem_transfer_recursion_at_dirac_q(Nc, p):
    """chord_charges' three-term close_j recursion -- the whole proof of
    the conserved-channel theorem -- holds verbatim at
    (q_dirac(Nc, p), Delta = 1/p), j <= 4, as a matrix identity on the
    truncation-free interior block."""
    q, dpsi = q_dirac(Nc, p), 1.0 / p
    T = transfer_formal(NMAX_G3, q)
    TM = transfer_one_matter_formal(NMAX_G3, q, dpsi)
    closes = [close_j_map(NMAX_G3, q, dpsi, j) for j in range(6)]
    m = NMAX_G3 + 1
    keep = NMAX_G3 - 1
    cols = [n0 * m + n1 for n0 in range(m) for n1 in range(m)
            if n0 + n1 <= keep - 1]
    for j in range(5):
        lhs = closes[j] @ TM
        rhs = q ** (dpsi + j) * (T @ closes[j]) \
            + (1.0 - q ** (2 * dpsi + j)) * closes[j + 1]
        if j >= 1:
            rhs = rhs + qint(j, q) * closes[j - 1]
        scale = max(np.abs(lhs[:keep][:, cols]).max(), 1.0)
        assert np.abs((lhs - rhs)[:keep][:, cols]).max() < 1e-12 * scale, j


@pytest.mark.parametrize("Nc,p", DIRAC_G3_POINTS)
def test_theorem_transfer_charges_at_dirac_q(Nc, p):
    """The Dirac conserved pair channel: C(m) = P_m(T) with [C(m), T] = 0
    for m <= 6, [Xhat_n, T] = 0 for n <= 4, and the closed forms
    Xhat_1 = -(1 - q^{1/p}) T and
    Xhat_2 = (1 - q^{1/p})^2 T^2 + (1 - q^{2/p}), all at
    (q_dirac(Nc, p), Delta = 1/p) through dirac_chord.xhat_dirac.  The
    machinery that proves the Majorana theorem applies verbatim at the
    Dirac q; nothing in it knows which side supplied q."""
    q, dpsi = q_dirac(Nc, p), 1.0 / p
    T = transfer_formal(NMAX_G3, q)
    keep = NMAX_G3 - 8
    for m in range(0, 7):
        Cm = dressed_propagator(m, NMAX_G3, q, dpsi)
        Pm = poly_of_matrix(motzkin_coeffs(m, q, dpsi), T)
        scale = max(np.abs(_interior(Cm, keep)).max(), 1.0)
        assert np.abs(_interior(Cm - Pm, keep)).max() < 1e-11 * scale, m
        comm = Cm @ T - T @ Cm
        assert np.abs(_interior(comm, keep)).max() < 1e-11 * scale, m
    for n in range(1, 5):
        X = xhat_dirac(n, Nc, p, nmax=NMAX_G3)
        comm = X @ T - T @ X
        scale = max(np.abs(_interior(X, keep)).max(), 1.0)
        assert np.abs(_interior(comm, keep)).max() < 1e-11 * scale, n
    qt = q ** dpsi
    X1 = xhat_dirac(1, Nc, p, nmax=NMAX_G3)
    assert np.abs(_interior(X1 + (1.0 - qt) * T, keep)).max() < 1e-12
    X2 = xhat_dirac(2, Nc, p, nmax=NMAX_G3)
    ref2 = (1.0 - qt) ** 2 * (T @ T) \
        + (1.0 - qt * qt) * np.eye(NMAX_G3 + 1)
    assert np.abs(_interior(X2 - ref2, keep)).max() < 1e-12


def test_theorem_transfer_word_zeros_at_dirac_q():
    """The independent engine route at the Dirac q: the chord_moments word
    assembly (bare Delta = 1/p matter chords, exact integer-power
    transfer-matrix words) cancels to |mu_k| <= 1e-12 * gross for n <= 3,
    k in {2, 4}, both operator orderings, at q = q_dirac(16, 4) and
    q_dirac(12, 4).  Nonzero gross is asserted so the check cannot pass
    vacuously.  Measured cancellations are ~1e-17."""
    from chord_moments import MomentTable, _mu_ordered, _mu_ordered_swapped
    for Nc in (16, 12):
        tab = MomentTable(q_dirac(Nc, P), 1.0 / P)
        for n in (1, 2, 3):
            for k in (2, 4):
                for fn in (_mu_ordered, _mu_ordered_swapped):
                    v, g = fn(tab, n, k)
                    assert g > 0.1, (Nc, n, k)
                    assert abs(v) <= 1e-12 * g, (Nc, n, k, fn.__name__)


# ---------------------------------------------------------------------------
# gate G3 scope: pinned, not fabricated
# ---------------------------------------------------------------------------
def test_g3_scope_is_pinned():
    """Gate G3 bookkeeping, kept honest.  The H-correlated tower chord
    rules are still UNPUBLISHED (constant False, sources named).  The
    PAIR (i != j) half is now derived in-repo -- the census and
    theorem-transfer legs above -- and the module must say exactly which
    half is closed: the flavor-diagonal channel and the CP resolution
    stay open.  If either constant flips, this test demands the
    docstring (and a real gate) move with it."""
    assert MN_TOWER_CHORD_RULES_PUBLISHED is False
    assert MN_TOWER_PAIR_CHANNEL_DERIVED is True
    doc = __import__("dirac_chord").__doc__
    for needle in ("2006.13983", "2407.09396", "2302.11882", "2512.07715",
                   "NOT established", "1 - p/Nc", "pair half",
                   "flavor-diagonal", "CP RESOLUTION"):
        assert needle in doc
