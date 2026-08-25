# DSSYK self-averaging, the Mₙ tower, and the QED sector

The falsifiable corner of the paper: the §3.5 self-averaging claims, plus
the operator content its U(N)/Dirac "standard model" story needs.
Model (Eq. 3.2): H = i^(p/2) Σ J ψψ…ψ with Var(J) = p!/N^(p−1),
λ = p²/N, ⟨W⟩ = Tr(W)/dim.

| claim | eq. | verdict |
|---|---|---|
| singlet Var⟨W⟩ ~ p!/Nᵖ | 3.27 | confirmed; rate exactly 2/C(N,p) |
| adjoint violations ≤ exp(−a√N) | 3.28 | verified exactly in double scaling (`exact_wick.py`) |

## Self-averaging

**Eq. 3.27.** The singlet rate is exactly 2/C(N,p) — one over the number
of couplings — which equals the paper's p!/Nᵖ times an O(1) factor. The
m₂ identity holds regardless of the paper's diagrammatics, so the
non-tautological content is m₄ and m₆ self-averaging at the same rate;
ED confirms that at N ≤ 18, within ~10% statistics.

**Eq. 3.28, settled beyond ED.** B_jk = Tr(ψⱼHψₖH)/dim is quadratic in
the couplings, so Wick contraction gives closed forms at any (N, p):

    Var(B_jk) = 4σ⁴ C(N−2, p−1),    E[B_jj] = σ² C(N,p) (1 − 2p/N)

verified against literal contraction enumeration (exact) and ED (~1%).
Along p = √N this gives ln rms = −¼ √N ln N, which beats exp(−a√N) for
every a, scanned to N = 1600.

At fixed p — the regime where the spectrum match lives — emergence is
only **polynomial**: the violation amplitude is N^(−(p−1)/2) and the
per-element violation N^(−(p+1)/2), so reaching 10⁻⁶ needs N ≈ 826 at
p = 4. The emergent symmetry is weakest exactly where the spectra are to
be compared. E[B_jj] = 0 at p = N/2, which is what produces the ratio
artifact noted below.

## The Mₙ tower and CP

Exact and asserted (`dirac.py`, `mn_spectroscopy.py`): [H, Q] = 0;
M₀ = Q; **M₁ = −(p/2)iH** (the paper's M₁ = H fails at O(1));
**M₂ = −Σ ċ†ċ** (the paper's + sign is wrong); Mₙ† = (−1)ⁿMₙ only for
n ≤ 2. CP must be the **unitary** particle-hole map, since M₁ ∝ iH rules
out an antiunitary one. CP of the traceless M₀ (= Q − Nc/2, with
C Q C = Nc − Q) is −1 and CP(M₁) = +1, both exact; raw Mₙ for n ≥ 2 are
CP-mixtures with 15–65% wrong-channel content, so the paper's
CP = (−1)^(n+1) lives in CP-*resolved* correlators.
`mn_spectroscopy.py` computes them: conserved n = 0, 1 exact at ω = 0,
broad continua at λ = 1.6.

## Campaigns

**U(N)/QED** (`qed_campaign.py`): U(N)-class self-averaging made exact
(the Dirac analog of the Wick closed forms, verified by enumeration and
ED); disorder-averaged charging curves with full covariance and GLS fits;
charged-vs-singlet spectral scales. The capacitive-vs-confining question
is stated with ΔAIC rather than asserted.

**Baryons** (`baryons.py`): ε-tensor baryons at Nc = 4–10 — exact 2-point
collapse onto charge sectors, mass extensivity, and the edge-curvature
hierarchy. First data on the paper's underived inter-baryon-force claim,
with hierarchy claims kept inside their error bars.

**Freeness** (`freeness.py`): Cui–Kolchmeyer argue early/late observer
algebras combine as a *free product* past the scrambling time. For the
pair E = iψ₀ψ₁, F(t) = iψ₂(t)ψ₃(t), freeness means
spec(E + F(t)) = the arcsine law. Measured at N = 12–24: the spectrum
lands on the arcsine law, the W₁ plateau falls like dim^(−0.9) (a
rigid-spectrum floor, far below sampling noise), and all alternating
words → 0. This is the first exhibition of the free-product transition in
a microscopic dS-candidate model. Caveat: the crossover sits at the
*dissipation* time (t𝒥 ≈ 1.2, nearly N-independent), so ED cannot
resolve the ln N scrambling-time separation where the claim properly
lives, and word plateaus at N ≥ 20 are statistics-limited. The crossover
is not quite N-independent either: t𝒥 runs 1.191, 1.172, 1.233, 1.273 over
N = 12–24, and the N = 16 → 24 rise is ~4σ against a measured resolution
of ~0.011. A drift is resolved; four points over a factor 2 in N simply
cannot say whether it is ln N or finite-size.

**Shape-invariant moments** (`moments_pipeline.py`, `moments_fast.py`):
CP-resolved spectral moment invariants (w₂, kurtosis) along fixed p = 4 —
map-robust targets for the exact N = ∞ computation.
**Large-q ED bench** (`largeq_bench.py`): ground truth for
`duality/largeq_anchor.py`.

## The approach to the conserved-channel zero, exact

The chord theorem (`duality/chord_charges.py`) makes every i ≠ j *pair*
moment of Aₙ vanish at N = ∞ (the flavor-diagonal part is the same order
and does not vanish for n ≥ 2 — see `duality/CHORD.md`). `exact_pair_moments.py` computes the same moments
**exactly at finite N** — Wick over H-slot matchings, with the signed
subset sums evaluated by a CRT residue-table dynamic program. No ED, no
sampling noise. Scope: n ≤ 3 and μ₀/μ₂/μ₄, words up to K = 10 H-nodes.

Along fixed p = 4 every |μ_k| with k ≥ 1 falls monotonically to the zero,
with tail slopes 2.96, 3.90, 4.82 for n = 1, 2, 3 — that is,
μ₂ ~ N^−(n+2), so the n = 2 channel is the N⁻⁴ one. (μ₀ of the n = 2
channel is not monotone: it changes sign near N = 12, the exact version
of an effect first seen in ED.) Along the double-scaling line
λ_chord = 2, |μ₂| ≈ c/(pN) with c = 28.9 → 29.8 for n = 2 over
(p, N) = (4,16)…(10,100), drifting slowly upward, and c ≈ 200–205 for
n = 3, flat to 2% and non-monotone. These are observed regularities; the
analytic 1/N chord correction behind them is the open theory half.

Gates: A₁ = −2pH forces (N−1) μ_k^pair + μ_k^diag = 0 at k = 2 and 4,
held at machine precision. Validated against the closed m ≤ 3
Venn-profile sums (retained in-tree as the reference), literal Isserlis
enumeration (exact match), dense ED, and the chord limit at the
Krawtchouk N⁻² rate. Figure: `pair_moments.png`.

**Majorana singlet-moment bench** (`mn_majorana_bench.py`) is the ED side
of `duality/chord_moments.py`. The raw width w₂ rises along fixed p
toward the λ → 0 corner — the tower channel turning on. Scope: the chord
result is an i ≠ j pair amplitude, and at dense-ED sizes the measured
channel is flavor-diagonal-dominated, so the bench does not yet probe the
N = ∞ zero; at matched λ the raw w₂ rises with N, which is the honest
negative. The (8, 4) point shows an exact per-instance conservation
absent at (12, 6) and (16, 8); the mechanism is unknown.

Numbers live in `results.txt` / `results.json`; figures and data JSONs sit
next to their scripts.

## Why the numbers are trustworthy

Instance-level ED caps near N = 20 and cannot reach double scaling, so
every double-scaling statement here is exact combinatorics rather than
extrapolation, and ED-window fits carry error bars and ΔAIC verdicts. The
construction was validated before the physics: Clifford algebra and
parity hold to machine precision, and the N mod 8 RMT pattern
(GOE/GUE/GSE) is reproduced with parity-sector resolution. Fast bitwise
builders are equivalence-tested against independent sparse constructions.
Read `rms_off`, not the p = N/2-artifacted ratio.

## Files

| file | purpose |
|------|---------|
| `syk.py` | sparse Jordan–Wigner operators (slow reference builders) |
| `pauli_strings.py` | bitwise Hamiltonian assembly (10–100×, equivalence-tested) |
| `dirac.py` | complex SYK: [H,Q] = 0, the Mₙ tower, unitary CP machinery |
| `exact_wick.py` | closed-form Eq. 3.28 statistics beyond ED; plus exact E[m4], Var(m4) with relVar(m4) = 4 relVar(m2)(1 + o(1)) |
| `exact_wick_dirac.py` | the Dirac counterpart: exact Var(B_ij) for U(N) emergence (rms ~ Nc^(−(p−1)/2)), and ED-free singlet-moment anchors |
| `exact_pair_moments.py` | exact finite-N pair moments of Aₙ (Wick + CRT residue-table DP, no ED; n ≤ 3, μ₀/μ₂/μ₄) |
| `mn_spectroscopy.py` | CP-resolved Mₙ spectral functions |
| `moments_pipeline.py` | CP-resolved moment invariants with statistics |
| `moments_fast.py` | charge-sector-blocked reimplementation (~140× at Nc = 12; dense equivalence measured at 5e-15 against a 1e-9 gate) |
| `qed_campaign.py` | U(N) self-averaging + charging curves |
| `baryons.py` | baryon sector: collapse, extensivity, hierarchy |
| `largeq_bench.py` | ED ground truth for the large-q anchor |
| `fold_bench.py` | folded/complex-time 4-pt ED lab; `--scan` is the fixed-λ decider |
| `fold_edge_probe.py` | the fold's spectral-edge mechanism, measured |
| `freeness.py` | the classical → free-product transition |
| `mn_majorana_bench.py` | Majorana singlet-channel moments (ED side of the chord computation) |
| `validate_syk.py`, `self_averaging.py` | RMT validation; the §3.5 measurement |
| `plot_self_averaging.py`, `plot_freeness.py`, `plot_pair_moments.py` | figures (read the JSONs) |
| `test_*.py` | asserting suites (179 tests, 168 fast) |

## Usage

```bash
pip install -r requirements.txt
python3 validate_syk.py         # RMT validation (~25 s)
python3 self_averaging.py       # the 3.27 measurement (~2 min)
python3 exact_wick.py           # exact 3.28 + double-scaling scan (~1 s)
python3 exact_pair_moments.py   # the approach to the chord zero (~6 min)
python3 dirac.py                # Dirac/QED identities + CP report (~5 s)
python3 mn_spectroscopy.py      # Mn spectral functions (seconds)
python3 qed_campaign.py         # QED campaign
python3 baryons.py              # baryon campaign
python3 freeness.py             # freeness campaign (hours)
python3 mn_majorana_bench.py    # singlet-moment bench (~1 h)
pytest -q -m 'not slow'         # 168 fast asserts (~30 s)
```

## Scope

This settles §3.5 — both equations, one of them exactly — and the
operator-level prerequisites of the U(N) story. The headline spectrum
match itself lives at fixed p, λ → 0 (see `duality/CHORD.md` for why),
with the CP-resolved spectroscopy as its pipeline.
