# The exact N = ∞ chord engine, and what it says about the tower

*Code: `chord.py` (engine + spectral scan), `chord_moments.py` (singlet
moments), `chord_charges.py` (the conserved-channel theorem),
`chord_fold.py` (the fold at N = ∞), and the ED side
`syk-self-averaging/mn_majorana_bench.py`. Every verdict below is asserted
from the committed JSONs in `test_chord*.py`.*

## The engine

Double-scaled SYK at N = ∞ (λ = 2p²/N fixed) is exactly solvable by chord
combinatorics (Berkooz et al. 1811.02584; conventions from the review
2407.09396, both equation-checked). `chord.py` implements the full
solution as finite sparse linear algebra: the transfer matrix
T|n⟩ = |n+1⟩ + [n]q|n−1⟩, matter chords with weight q^Δ, and one- and
two-matter sectors that evaluate **any operator ordering on any
complex-time contour** — crossed and Euclidean-folded included. No special
functions.

Three independent legs hold it up (`test_chord.py`):

1. A brute-force chord-diagram enumerator. It shares nothing with the
   transfer-matrix bookkeeping and agrees to machine precision on every
   topology, fermionic signs included.
2. The published closed forms: the θ-integral 2-pt function (review
   Eq. 2.30), μ₂ = 2(1−q^Δ), and the Touchard–Riordan moments.
3. A chord-free check: the exact finite-N crossing average (a Krawtchouk
   sum) approaches q = exp(−2p²/N) at the measured O(p⁴/N²) rate.

Units bridge, exact and tested: σ_H² = C(N,p) p!/N^(p−1), giving
β_chord = βJ · exp(−λ/8)/√λ at N = ∞. Against finite-N ED the chord
thermal 2-pt agrees to **0.6–1.6%** (the large-q closed forms sit at a
30–55% floor at the same sizes).

Stated caps: at most 2 open matter chords (covers all 2-pt and 4-pt
orderings); folded contours below λ ≈ 0.5 exceed float headroom (the
λ = 0.3 point fails truncation-doubling; drifts stored per row in
`chord_fold.json`).

## Discrete tower or continuum? — Continuum. Sharper: the channels are conserved.

`chord.py --spectra` computes the exact smooth spectral functions S_Δ(ω)
at T_B = ∞ for Δ ∈ {1, ½, ¼} and the Δ → 0 length channel, at
λ = 2.0 → 0.05, with a fine window at the semiclassical scale where a
't Hooft tower would sit.

- **Zero interior maxima, zero nodes — every channel, every λ.** The exact
  spectral function converges to the smooth sech^2Δ transform linearly in
  λ (Δ = ½: deviation 0.375 → 0.006, about 0.13 λ; worst channel Δ = 1:
  about 0.25 λ). The mass-encoding *nodes* of BOOST.md are absent too.
  The tests recompute these diagnostics from the committed curves.
- The chord-number growth ⟨n(t)⟩ — the de Sitter length — approaches
  −(2/λ) ln sech(Jc t) at the same linear rate. New exact-in-λ data.

**The moments make the verdict algebraic — and it is now a theorem**
(`chord_moments.py`, proof in `chord_charges.py`). For the paper's own
channels — the singlet bilinears A_n = Σᵢ ψᵢ (ad_H)ⁿ ψᵢ, assembled by a
binomial expansion verified as a per-instance operator identity against
dense ED (`test_mn_majorana.py`):

> **Measured: μ₂ = μ₄ = μ₆ = 0 to machine precision (relative to the
> pre-cancellation gross sum), for n = 1 … 5, at every λ in [0.1, 4] and
> across the matter-weight scan Δψ ∈ [0.1, 1] (committed in
> `chord_moments.json`). Proven: every spectral moment μ_k, k ≥ 1, of
> every A_n pair channel vanishes at N = ∞.** The mechanism: in a
> sequential pair word the fermion pair factors through the zero-chord
> sector as a dressed propagator C(m) = close∘T_Mᵐ∘open, and a
> three-term recursion with scalar coefficients (close_j T_M =
> q^(Δ+j) T close_j + [j] close_{j−1} + (1−q^(2Δ+j)) close_{j+1}) forces
> C(m) to be a **polynomial in the transfer matrix T** — the Motzkin-path
> polynomial with flat weight q^(Δ+j)T, down [j], up 1−q^(2Δ+j). The
> chord representation of A_n is then itself a polynomial in T:
> conserved. Explicitly, X̂₁ = −(1−q^Δ)T (the N = ∞ image of the exact
> law A₁ = −2pH) and X̂₂ = (1−q^Δ)²T² + (1−q^(2Δ)). Random O_Δ
> operators *do* propagate (μ₂ = 2(1−q^Δ), kurtosis → 3 + 1/Δ; closed
> forms, tested) — their insertion q^(Δn̂) does **not** commute with T.
> The H-correlated bilinears the duality needs cannot propagate; the
> random ones can. Every step is asserted in `test_chord_charges.py`.

This is the exact-in-λ version of the known large-q cancellation
(LARGEQ.md), now with its algebraic origin.

**Scope, sharpened by adversarial review.** The chord moments are i ≠ j
flavor-*pair* amplitudes. The ED bench measures the full channel, which at
dense-ED sizes is flavor-diagonal-dominated (the pair part of μ₀ is a few
percent and opposite in sign at N = 12–14). So the bench does not yet
probe the chord cancellation. What it does establish:

- Along fixed p with λ falling, the raw width w₂ **rises** (n = 2:
  0.30 → 1.20 over N = 10 → 20 at p = 4). The tower channel turns on
  exactly where the match must be made.
- The honest negative: at the near-matched-λ pair (p=4, N=10) ↔
  (p=6, N=22) the raw w₂ also rises with N. The approach to the N = ∞
  zero is not yet visible at dense-ED sizes. A former (1−2p/N)²
  "artifact correction" was removed by the review: the (8,4) all-moments
  zero it leaned on is an exact per-instance conservation special to
  (8,4) — mechanism unknown, absent at (12,6) and (16,8). Firm fixed-λ
  asymptotics need N ≫ 2p, beyond dense ED.

**Consequence.** At N = ∞ the singlet channels are neither discrete nor
continuous — they are conserved, with all weight at ω = 0. At finite N
they are smooth continua whose weight grows toward the fixed-p, λ → 0
corner. The meson tower, if the duality produces one, is a **1/N effect,
invisible at every order of the double-scaled expansion**. It must be
extracted at fixed p, λ → 0 — where 't Hooft mesons are themselves the
leading large-N states. Nothing here falsifies the duality; it locates it.

Remaining caveats: the chord evaluation uses the standard leading-chord
matter weight q^(1/p); the conserved-channel zeros hold at any weight, so
this is not the source. The Dirac CP-resolved data (`moments.json`) awaits
the complex-SYK chord machinery — the Dirac ↔ Majorana channel map is not
one-to-one.
