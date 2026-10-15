# Self-averaging, the Mₙ tower, and the QED sector

This is the falsifiable corner of the paper. The §3.5 self-averaging
claims, plus the operator content the U(N)/Dirac "standard model" story
needs.

The model is the paper's Eq. 3.2: H = i^(p/2) Σ J ψψ…ψ, with
Var(J) = p!/N^(p−1), λ = p²/N, and ⟨W⟩ = Tr(W)/dim.

| claim | eq. | verdict |
|---|---|---|
| singlet Var⟨W⟩ ~ p!/Nᵖ | 3.27 | holds; the rate is exactly 2/C(N,p) |
| adjoint violations ≤ exp(−a√N) | 3.28 | exact in double scaling (`exact_wick.py`) |

## Self-averaging

The paper argues from two diagrams. We can do better than that.

**Eq. 3.27.** The singlet rate is exactly 2/C(N,p). That is one over the
number of couplings, and it equals the paper's p!/Nᵖ times an O(1)
factor. The m₂ identity holds whatever the diagrammatics say, so the
content that is not tautological is m₄ and m₆ self-averaging at the same
rate. ED gives that at N ≤ 18, within about 10% statistics.

**Eq. 3.28, settled beyond ED.** B_jk = Tr(ψⱼHψₖH)/dim is quadratic in
the couplings. So Wick contraction gives closed forms at any (N, p):

    Var(B_jk) = 4σ⁴ C(N−2, p−1),    E[B_jj] = σ² C(N,p) (1 − 2p/N)

Both are checked against literal contraction enumeration, which is exact,
and against ED, to about 1%.

Along p = √N this gives ln rms = −¼ √N ln N. That beats exp(−a√N) for
every a, scanned to N = 1600. The paper's bound is right, and
conservative.

**But at fixed p, emergence is only polynomial.** The violation amplitude
rms(Bⱼₖ) falls as N^(−(p−1)/2), and the per-element violation
rms(Bⱼₖ)/E[Bⱼⱼ] as N^(−(p+1)/2). A per-element violation of 10⁻⁶ needs
N ≈ 826 at p = 4 (exactly 1.01 × 10⁻⁶ there). The paper places the
spectrum match at fixed p, so the emergent symmetry is weakest exactly
where the spectra are to be compared.

E[B_jj] = 0 at p = N/2. That is what produces the ratio artifact noted
below.

## The Mₙ tower and CP

Exact, and asserted in `dirac.py` and `mn_spectroscopy.py`:

- [H, Q] = 0, and M₀ = Q.
- **M₁ = −(p/2)iH.** The paper's M₁ = H fails at O(1).
- **M₂ = −Σ ċ†ċ.** The paper's + sign is wrong.
- Mₙ† = (−1)ⁿMₙ only for n ≤ 2.

CP has to be the **unitary** particle-hole map, because M₁ ∝ iH rules out
an antiunitary one. CP of the traceless M₀ (= Q − Nc/2, with C Q C =
Nc − Q) is −1, and CP(M₁) = +1. Both are exact.

Raw Mₙ for n ≥ 2 are CP mixtures, with 15–65% wrong-channel content. So
the paper's CP = (−1)^(n+1) lives in CP-*resolved* correlators, not in
the raw ones. `mn_spectroscopy.py` computes those: n = 0 and n = 1 are
conserved and exact at ω = 0, and there are broad continua at λ = 1.6.

## Campaigns

**U(N)/QED** (`qed_campaign.py`). U(N)-class self-averaging made exact,
which is the Dirac analog of the Wick closed forms, checked by
enumeration and by ED. Disorder-averaged charging curves with full
covariance and GLS fits. Charged-against-singlet spectral scales. The
capacitive-against-confining question is stated with ΔAIC, not asserted.

**Baryons** (`baryons.py`). ε-tensor baryons at Nc = 4–10. Exact 2-point
collapse onto charge sectors, mass extensivity, and the edge-curvature
hierarchy. First data on the paper's underived inter-baryon-force claim,
with the hierarchy claims kept inside their error bars.

**Freeness** (`freeness.py`). Cui–Kolchmeyer argue that early and late
observer algebras combine as a *free product* past the scrambling time.
Take the pair E = iψ₀ψ₁ and F(t) = iψ₂(t)ψ₃(t). Freeness means
spec(E + F(t)) is the arcsine law.

Measured at N = 12–24, the spectrum lands on the arcsine law, the W₁
plateau falls like dim^(−0.9), and every alternating word goes to zero.
The plateau is a rigid-spectrum floor, far below sampling noise. To our
knowledge this is the first exhibition of the free-product transition in
a microscopic dS-candidate model.

The caveat is real. The crossover sits at the *dissipation* time
(t𝒥 ≈ 1.2, nearly N-independent), so ED cannot resolve the ln N
scrambling-time separation where the claim properly lives. Word plateaus
at N ≥ 20 are statistics-limited. The crossover is not quite
N-independent either: t𝒥 runs 1.191, 1.172, 1.233, 1.273 over N = 12–24,
and the N = 16 → 24 rise is about 4σ against a rerun-measured resolution
of 0.011, and about 7σ on the artifact's own propagated per-point errors
(0.005–0.013, `freeness.verdict_scalars`). A drift is resolved. Four
points over a factor of 2 in N cannot say whether it is ln N or finite
size.

**Shape-invariant moments** (`moments_pipeline.py`, `moments_fast.py`).
CP-resolved spectral moment invariants (w₂, kurtosis) along fixed p = 4.
These are map-robust targets for the exact N = ∞ computation.

**Large-q ED bench** (`largeq_bench.py`). Ground truth for
`duality/largeq_anchor.py`.

## The approach to the pair-channel zero, exact

The chord theorem (`duality/chord_charges.py`) makes every i ≠ j *pair*
moment of Aₙ vanish at N = ∞. The diagonal part is the same order and
does not vanish for n ≥ 2. Read `duality/CHORD.md` before quoting that.

`exact_pair_moments.py` computes the same moments **exactly at finite
N**. Wick over the H-slot matchings, with the signed subset sums
evaluated by a CRT residue-table dynamic program. No ED, no sampling
noise. Scope: n ≤ 3, moments μ₀/μ₂/μ₄, words up to K = 10 H-nodes.

Along fixed p = 4 every |μ_k| with k ≥ 1 falls monotonically to zero. The
fitted tail slopes are 2.96, 3.90, 4.82 for n = 1, 2, 3.

(μ₀ of the n = 2 channel is not monotone. It changes sign near N = 12,
which is the exact counterpart of an effect first seen in ED. The word is
"counterpart" because the μ₀ subtraction convention differs:
ensemble-level here, per-instance in the ED benches. That gap is of order
Cov(⟨X⟩,⟨Y⟩) and matters only at small N. Both conventions show the sign
change, and μ_k≥1 needs no subtraction at all.)

Along the double-scaling line λ_chord = 2, |μ₂| ≈ c/(pN). For n = 2,
c = 28.9 → 29.8 over (p, N) = (4,16)…(10,100), drifting slowly upward.
For n = 3, c ≈ 200–205, flat to 2% and not monotone.

**Those float regularities are now exact laws** (`pair_moment_laws.py`).
P(N) = 2μ_k C(N,p)^m is a polynomial in N, so interpolation over
exact-Fraction samples gives the closed form of μ_k(N) at fixed p,
certified out of sample. Results:

- The deficits are exactly n + 2 through n = 4. The fitted slopes above
  were their shadows.
- The leading coefficients have certified closed forms in p:
  A(n=1) = −16p²(p−1) and A(n=2) = −64p³(p−1)(4p−5).
- μ₄ dies one power faster, and at n = 1 shares its leading coefficient
  with μ₂(n=2).
- The exact double-scaling constants c₂(p) rise over p = 4..16; the
  ladder below carries them to p = 48 and identifies the limit.

**The ladder** (`pair_moment_ladder.py`) carries this to the K = 12 word
frontier, with K = k + 2n the word length:

- deficit = K/2 + 1 at every one of the 55 laws, k = 2..10 — so
  μ₂ ~ N^−(n+2) through n = 5, μ₄ ~ N^−(n+3) through n = 4, and
  μ₆ ~ N^−(n+4) through n = 3;
- A(n=1, k) = A(n=2, k−2) exactly at every rung and every p, and every
  deeper step differs (K = 10, p = 2: 21 : 21 : 14 : 9 × (−2¹⁷) over
  n = 1..4) — the coincidence above is the first step of a ladder and
  nothing more;
- A(n=3) = −256p⁴(p−1)(34p² − 108p + 91), certified at p = 18, so the
  k = 2 family is −2^(2n+2) p^(n+1)(p−1) Q_n(p) with deg Q_n = n − 1
  through n = 3 (three data points; not a proven pattern);
  the other K = 8 family is A(1,6) = A(2,4) = −256p⁴(p−1)(34p² − 93p + 64),
  same prefactor, same leading 34p², certified at p = 18;
- c₂(p) exact to p = 48 (31.654), extrapolating to 32.3341891(1);
  c₃(p) exact to p = 32 (`double_scaling_law.json`), extrapolating to
  223.63(3).

**Those two limits are closed forms** (`double_scaling_law.py`). On the
line λ_chord = λ, with q = e^(−λ),

    c₂(λ) = 4λ(1−q)²(5+3q),    c₃(λ) = 4λ(1−q)³(35 + 56q + 35q² + 10q³),

checked on three lines each (λ = 1, 2, 4; six digits at K = 6). They
are the Δψ → 0 slopes λP_K(q) of the N = ∞ *diagonal* chord channel
(`duality/diagonal_chord.py`), which depend on the word length K only:
the O(Δ) parts of the pair and diagonal halves cancel on every
double-scaling line, the n ≥ 2 shadow of the exact A₁ identity, and the
finite-N table in the artifact shows the cancellation directly. The
λ-expansion of the closed form fixes the large-p growth of every 1/N
coefficient; the N⁻⁵ one is a closed form too,
B(2,2;p) = 64p³(p−1)(11p³ − 36p² + 29p − 3), fit on p = 4..12 and
reproducing the committed p = 2 value out of sample.

What is left for the chord side is to *derive* these targets, the
cancellation above all.

**The full channel at fixed p** (`full_channel_laws.py`). Summing the
two halves exactly, the full A₂ width along fixed p peaks at
λ_chord ≈ 4/3 (N = 24 at p = 4) and then dies as w₂ = 8p(p−2)/N in
bandwidth units, with the coefficient 8p(p−2) fit on p = 2, 4, 6 and
confirmed out of sample at p = 8, 10, 12 (the two leading
coefficients are 64p²(p−1)²/N² and 512p³(p−1)²(p−2)/N³). In physical
units (J = 1) that is an N-independent mean-square frequency 8p(p−2):
the channel keeps a finite O(J) width into the λ → 0 corner. At fixed
λ the N = ∞ full channel is exact too (`duality/full_channel_chord.py`).

Read the peak with its denominator. μ₀ is identity-subtracted, and A₂
overlaps the conserved H². At N = ∞ the split is a theorem
(`duality/singlet_tower.py`): A₂ is a conserved H² piece plus a matter
chord of dimension 2, with dynamical width 2(1 − q²). At finite N the H²
component can be removed in exact arithmetic, and then the width has no
maximum beyond the smallest sizes: it sits near 2 at the small-N end
and falls monotonically to the same 8p(p−2)/N. The peak at N = 24 is
conserved weight leaving the denominator (85% of μ₀ at N = 10, 19% at
N = 24 for p = 4).

**Gates.** A₁ = −2pH forces (N−1)μ_k^pair + μ_k^diag = 0 at k = 2 and 4,
and that holds at machine precision. The computation is checked against
the closed m ≤ 3 Venn-profile sums (kept in-tree as the reference),
literal Isserlis enumeration (exact match), dense ED, and the chord limit
at the Krawtchouk N⁻² rate. Figure: `pair_moments.png`.

**Majorana singlet-moment bench** (`mn_majorana_bench.py`) is the ED side
of `duality/chord_moments.py`. The raw width w₂ rises along fixed p
inside the bench's window, N ≤ 20. The exact engine shows that window
ends before the peak (N = 24 at p = 4), after which the width dies as
8p(p−2)/N. The rise inside the window is the near side of that peak,
not the tower turning on (ERRATA E4).

Scope, stated plainly. The chord result is an i ≠ j pair amplitude. At
dense-ED sizes the measured channel is diagonal-dominated, so the bench
does not probe the N = ∞ pair zero. At matched λ the raw w₂ rises with
N, which is the honest negative for the raw ratio. The full channel is
another matter: with the commutant of H projected out per instance
(`projected_channel.py`), ED does sit on the N = ∞ full-channel law —
dynamical width 2.03, 2.03, 2.01 at p = 6, N = 16–20 against
2(1 − q²) = 2.00, where the raw width is 0.29–0.41. The (8, 4) point shows an exact
per-instance conservation that is absent at (12, 6) and (16, 8). The
mechanism is now proven (`conserved_points.py`): the spin-trace map
Σᵢψᵢ X ψᵢ weighs Clifford grade g by (−1)^g(N − 2g), every power of H
at N = 8 has grades {0, 4, 8} only, so every Aₙ is a polynomial in H
with central coefficients. The same algebra gives two exact identities
at every (N, p), A₁ = −2ΣᵢHᵢ and A₂ = 4ΣᵢHᵢ², with Hᵢ the part of H
containing flavor i, and they are the first cases of a general one:
every Aₙ is a fixed combination of the words Hⁿ with Hᵢ on two or more
slots, summed over i (n = 3: 4Σᵢ(HᵢHHᵢ + HHᵢ² − 2Hᵢ²H) − 8ΣᵢHᵢ³;
asserted for n ≤ 5). That identity is what the N = ∞ tower law is built
on.

Numbers live in `results.txt` and `results.json`. Figures and data JSONs
sit next to their scripts.

## Why the numbers are trustworthy

Instance-level ED caps near N = 20 and cannot reach double scaling. So
every double-scaling statement here is exact combinatorics, not
extrapolation. ED-window fits carry error bars and ΔAIC verdicts.

The construction was checked before the physics. Clifford algebra and
parity hold to machine precision. The N mod 8 RMT pattern (GOE/GUE/GSE)
comes out right, with parity-sector resolution. Fast bitwise builders are
equivalence-tested against independent sparse constructions.

Read `rms_off`, not the ratio, which is artifacted at p = N/2.

## Files

| file | purpose |
|------|---------|
| `syk.py` | sparse Jordan–Wigner operators (slow reference builders) |
| `pauli_strings.py` | bitwise Hamiltonian assembly (10–100×, equivalence-tested) |
| `dirac.py` | complex SYK: [H,Q] = 0, the Mₙ tower, unitary CP machinery |
| `exact_wick.py` | closed-form Eq. 3.28 statistics beyond ED, plus exact E[m4] and Var(m4) with relVar(m4) = 4 relVar(m2)(1 + o(1)) |
| `exact_wick_dirac.py` | the Dirac counterpart: exact Var(B_ij) for U(N) emergence (rms ~ Nc^(−(p−1)/2)), and ED-free singlet-moment anchors |
| `exact_pair_moments.py` | exact finite-N pair moments of Aₙ (Wick + CRT residue-table DP, no ED; n ≤ 3, μ₀/μ₂/μ₄; exact-Fractions assembly to K = 12) |
| `pair_moment_laws.py` | the exact 1/N laws by interpolation: deficits n+2 through n = 4, closed forms A(n=1,2)(p) certified at p = 14, exact double-scaling c(p) |
| `pair_moment_ladder.py` | the μ₆ / K = 12 extension and the word-length ladder: deficit K/2+1 at all 55 laws, the exact n = 1 → 2 coefficient step, the K = 12 frontier, A(n=3)(p) certified at p = 18, c₂(p) exact to p = 48 (multiprocessing; `--refit`, `--ds-extend [--n 3]`) |
| `double_scaling_law.py` | the closed forms of the double-scaling constants, c_{n,2}(λ) = λP_{2n+2}(e^(−λ)), the exact sequences on three lines, the finite-N diagonal–pair cancellation, the B(2,2;p) closed form (multiprocessing; `--resume`) |
| `full_channel_laws.py` | the FULL A_n channel at fixed p, exact: the peak at λ_chord ≈ 4/3 and the certified law N w₂ → 8p(p−2) |
| `conserved_points.py` | the (8, 4) per-instance conservation explained (spin-trace map, Cl₈ grades), the identities A₁ = −2ΣHᵢ, A₂ = 4ΣHᵢ², and the general mixed-word identity for every Aₙ |
| `projected_channel.py` | finite-N ED split of a channel into stationary and dynamical weight (the whole commutant of H removed, degeneracies included) |
| `mn_spectroscopy.py` | CP-resolved Mₙ spectral functions |
| `moments_pipeline.py` | CP-resolved moment invariants with statistics |
| `moments_fast.py` | charge-sector-blocked reimplementation (~140× at Nc = 12; dense equivalence measured at 5e-15 against a 1e-9 gate) |
| `qed_campaign.py` | U(N) self-averaging and charging curves |
| `baryons.py` | the baryon sector: collapse, extensivity, hierarchy |
| `largeq_bench.py` | ED ground truth for the large-q anchor |
| `fold_bench.py` | folded and complex-time 4-pt ED lab; `--scan` is the fixed-λ decider |
| `fold_edge_probe.py` | the fold's spectral-edge mechanism, measured |
| `freeness.py` | the classical → free-product transition |
| `mn_majorana_bench.py` | Majorana singlet-channel moments (the ED side of the chord computation) |
| `validate_syk.py`, `self_averaging.py` | RMT validation; the §3.5 measurement |
| `plot_self_averaging.py`, `plot_freeness.py`, `plot_pair_moments.py` | figures (they read the JSONs) |
| `test_*.py` | the asserting suites (265 tests, 254 fast) |

## Usage

```bash
pip install -e ".[dev]"      # from the repo root, once
python3 validate_syk.py         # RMT validation                  (~25 s)
python3 self_averaging.py       # the 3.27 measurement            (~2 min)
python3 exact_wick.py           # exact 3.28 + double-scaling scan (~1 s)
python3 exact_pair_moments.py   # the approach to the chord zero  (~6 min)
python3 pair_moment_laws.py --deep  # the exact 1/N laws          (~27 min)
python3 pair_moment_ladder.py       # the K = 12 ladder campaign   (~7.5 h, 7 workers; --resume)
python3 pair_moment_ladder.py --ds-extend 32,36,40,44,48   # c2(p) to p = 48
python3 dirac.py                # Dirac/QED identities + CP report (~5 s)
python3 mn_spectroscopy.py      # Mn spectral functions           (seconds)
python3 qed_campaign.py         # the QED campaign
python3 baryons.py              # the baryon campaign
python3 freeness.py             # the freeness campaign           (hours)
python3 mn_majorana_bench.py    # the singlet-moment bench        (~1 h)
pytest -q -m 'not slow'         # 254 fast tests                  (~45 s)
```

## Scope

This settles §3.5, both equations, one of them exactly. It also settles
the operator-level prerequisites of the U(N) story.

The paper places the spectrum match at fixed p, λ → 0. See
`exact_wick.py` for why the emergent gauge symmetry is only polynomially
accurate there, and `duality/CHORD.md` for where the N = ∞ weight
actually sits. The CP-resolved spectroscopy is the pipeline for it.
