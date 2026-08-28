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
is N^(−(p−1)/2) and the per-element violation is N^(−(p+1)/2). Reaching
10⁻⁶ needs N ≈ 826 at p = 4. Fixed p is where the spectrum match lives,
so the emergent symmetry is weakest exactly where the spectra have to be
compared.

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
The plateau is a rigid-spectrum floor, far below sampling noise. This is
the first exhibition of the free-product transition in a microscopic
dS-candidate model.

The caveat is real. The crossover sits at the *dissipation* time
(t𝒥 ≈ 1.2, nearly N-independent), so ED cannot resolve the ln N
scrambling-time separation where the claim properly lives. Word plateaus
at N ≥ 20 are statistics-limited. The crossover is not quite
N-independent either: t𝒥 runs 1.191, 1.172, 1.233, 1.273 over N = 12–24,
and the N = 16 → 24 rise is about 4σ against a measured resolution of
0.011. A drift is resolved. Four points over a factor of 2 in N cannot
say whether it is ln N or finite size.

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
- μ₄ dies one power faster, and shares its leading coefficient with
  μ₂(n+1).
- The exact double-scaling constants run c₂ = 28.92 → 30.53 over
  p = 4..16, still rising.

What is left for the chord side is to *derive* these targets.

**Gates.** A₁ = −2pH forces (N−1)μ_k^pair + μ_k^diag = 0 at k = 2 and 4,
and that holds at machine precision. The computation is checked against
the closed m ≤ 3 Venn-profile sums (kept in-tree as the reference),
literal Isserlis enumeration (exact match), dense ED, and the chord limit
at the Krawtchouk N⁻² rate. Figure: `pair_moments.png`.

**Majorana singlet-moment bench** (`mn_majorana_bench.py`) is the ED side
of `duality/chord_moments.py`. The raw width w₂ rises along fixed p
toward the λ → 0 corner. That is the tower channel turning on.

Scope, stated plainly. The chord result is an i ≠ j pair amplitude. At
dense-ED sizes the measured channel is diagonal-dominated, so the bench
does not probe the N = ∞ zero. At matched λ the raw w₂ rises with N,
which is the honest negative. The (8, 4) point shows an exact
per-instance conservation that is absent at (12, 6) and (16, 8). The
mechanism is unknown.

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
| `test_*.py` | the asserting suites (189 tests, 178 fast) |

## Usage

```bash
pip install -r requirements.txt
python3 validate_syk.py         # RMT validation                  (~25 s)
python3 self_averaging.py       # the 3.27 measurement            (~2 min)
python3 exact_wick.py           # exact 3.28 + double-scaling scan (~1 s)
python3 exact_pair_moments.py   # the approach to the chord zero  (~6 min)
python3 pair_moment_laws.py --deep  # the exact 1/N laws          (~27 min)
python3 dirac.py                # Dirac/QED identities + CP report (~5 s)
python3 mn_spectroscopy.py      # Mn spectral functions           (seconds)
python3 qed_campaign.py         # the QED campaign
python3 baryons.py              # the baryon campaign
python3 freeness.py             # the freeness campaign           (hours)
python3 mn_majorana_bench.py    # the singlet-moment bench        (~1 h)
pytest -q -m 'not slow'         # 178 fast asserts                (~30 s)
```

## Scope

This settles §3.5, both equations, one of them exactly. It also settles
the operator-level prerequisites of the U(N) story.

The headline spectrum match itself lives at fixed p, λ → 0. See
`exact_wick.py` for why the emergent gauge symmetry is only polynomially
accurate there, and `duality/CHORD.md` for where the N = ∞ weight
actually sits. The CP-resolved spectroscopy is the pipeline for it.
