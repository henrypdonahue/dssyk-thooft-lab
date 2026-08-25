# The exact N = ∞ chord engine, and what it says about the tower

*Code: `chord.py` (engine + spectral scan), `chord_moments.py` (singlet
moments), `chord_charges.py` (the conserved-channel theorem),
`chord_fold.py` (the fold at N = ∞), and the ED side
`syk-self-averaging/mn_majorana_bench.py`. Every verdict below is
asserted from the committed JSONs in `test_chord*.py`.*

## The engine

Double-scaled SYK at N = ∞ (λ = 2p²/N fixed) is exactly solvable by chord
combinatorics (Berkooz et al. 1811.02584); conventions follow the review
2407.09396, and both are equation-checked. `chord.py` implements the
solution as finite sparse linear algebra: the transfer matrix
T|n⟩ = |n+1⟩ + [n]_q|n−1⟩, matter chords of weight q^Δ, and one- and
two-matter sectors. Those sectors evaluate **any operator ordering on any
complex-time contour** — crossed and Euclidean-folded included. No
special functions.

Three independent legs hold it up (`test_chord.py`):

1. A brute-force chord-diagram enumerator, sharing nothing with the
   transfer-matrix bookkeeping. It agrees to machine precision on every
   topology, fermionic signs included.
2. The published closed forms: the θ-integral 2-pt function (review
   Eq. 2.30), μ₂ = 2(1−q^Δ), and the Touchard–Riordan moments.
3. A chord-free check: the exact finite-N crossing average (a Krawtchouk
   sum) approaches q = exp(−2p²/N) at the measured O(p⁴/N²) rate.

The units bridge is exact and tested: σ_H² = C(N,p) p!/N^(p−1), giving
β_chord = βJ · exp(−λ/8)/√λ at N = ∞. Against finite-N ED the chord
thermal 2-pt agrees to **0.6–1.6%**, where the large-q closed forms sit
at a 30–55% floor at the same sizes.

Stated caps: at most two open matter chords, which covers all 2-pt and
4-pt orderings; folded contours below λ ≈ 0.5 exceed float headroom, and
the λ = 0.3 point fails truncation-doubling. Per-row drifts are stored in
`chord_fold.json`.

## Discrete tower or continuum? Continuum — and the pair channel is conserved

`chord.py --spectra` computes the exact smooth spectral functions S_Δ(ω)
at T_B = ∞, covering Δ ∈ {1, ½, ¼} and the Δ → 0 length channel over
λ = 2.0 → 0.05, with a fine window at the semiclassical scale where a
't Hooft tower would sit.

- **Zero interior maxima, zero nodes — every channel, every λ.** The
  exact spectral function converges to the smooth sech^2Δ transform
  linearly in λ: for Δ = ½ the deviation falls 0.375 → 0.006, about
  0.13 λ; the worst channel, Δ = 1, is about 0.25 λ. The mass-encoding
  *nodes* of [BOOST.md](BOOST.md) are absent too.
- The chord-number growth ⟨n(t)⟩ is the de Sitter length. It approaches
  −(2/λ) ln sech(J_c t) at the same linear rate — new exact-in-λ data.

The moments make the verdict algebraic, and it is a **theorem**
(`chord_moments.py` measures it, `chord_charges.py` proves it). The
paper's own channels are the singlet bilinears Aₙ = Σᵢ ψᵢ (ad_H)ⁿ ψᵢ,
assembled by a binomial expansion that `test_mn_majorana.py` verifies as
a per-instance operator identity against dense ED. For these channels:

> **Every spectral moment μ_k, k ≥ 1, of the i ≠ j PAIR amplitude of
> every Aₙ vanishes at N = ∞ — every λ, every matter weight, both
> orderings, odd k included.** (Scope matters here: the full channel also
> has a flavor-diagonal part of the same order, which does *not* vanish
> for n ≥ 2. See "What it means" below.) Measured: μ₂ = μ₄ = μ₆ = 0 to machine precision (1e-17 to
> 1e-18 relative to the pre-cancellation gross sum) for n = 1…5, on the
> λ grid [0.1, 4] at Δψ = ¼ and across the matter-weight scan
> Δψ ∈ [0.1, 1] at λ = 1.

The mechanism: in a sequential pair word the fermion pair factors through
the zero-chord sector as a dressed propagator C(m) = close∘T_Mᵐ∘open. The
three-term recursion

    close_j T_M = q^(Δ+j) T close_j + [j] close_{j−1} + (1−q^(2Δ+j)) close_{j+1}

has *scalar* coefficients — both matching conditions reduce to q-number
identities independent of the chord occupation — which forces C(m) to be
the Motzkin-path polynomial in T (flat weight q^(Δ+j)T, down [j], up
1−q^(2Δ+j)). The chord representation of Aₙ is then itself a polynomial
in T, hence conserved. Explicitly Â₁ = −(1−q^Δ)T, the N = ∞ image of the
exact law A₁ = −2pH, and Â₂ = (1−q^Δ)²T² + (1−q^(2Δ)).

The contrast that makes it physical: a random operator O_Δ inserts
q^(Δn̂), which does **not** commute with T, so random matter propagates
(μ₂ = 2(1−q^Δ), kurtosis → 3 + 1/Δ) while the H-correlated bilinears the
duality needs cannot. Those closed forms are tested;
`test_chord_charges.py` asserts every step of the proof.

This is the exact-in-λ version of the known large-q cancellation
([LARGEQ.md](LARGEQ.md)), now with its algebraic origin.

## The finite-N approach to the zero

`syk-self-averaging/exact_pair_moments.py` computes the same moments
exactly at finite N — Wick over H-slot matchings with the chord-diagram
sign law, the signed subset sums evaluated by a CRT residue-table dynamic
program. No ED, no sampling noise; the whole 18-point scan takes 332 s,
with the largest row (N = 256) about 2 s per channel. Scope:
n ≤ 3 and μ₀/μ₂/μ₄, words up to K = 10 H-nodes.

Along **fixed p = 4** every μ_k with k ≥ 1 falls monotonically to the
zero, with tail slopes 2.96, 3.90, 4.82 for n = 1, 2, 3 — i.e.
μ₂ ~ N^−(n+2), so the n = 2 channel is the N⁻⁴ one. (μ₀ of the n = 2
channel is not monotone: it changes sign near N = 12, the exact version
of an effect first seen in ED.)

Along the **double-scaling line λ_chord = 2**, |μ₂| ≈ c/(pN) with
c = 28.9, 28.9, 29.4, 29.8 for n = 2 over (p, N) = (4,16)…(10,100), and
c = 203, 200, 202, 205 for n = 3 — the n = 2 constant drifting slowly
upward, the n = 3 one flat to 2% and non-monotone. These are observed
regularities; the analytic 1/N chord correction behind them is open, and
these are the numbers it must reproduce.

Validation: literal Isserlis enumeration (exact match), dense ED, the
chord limit at the Krawtchouk N⁻² rate, and the per-row gate
(N−1) μ_k^pair + μ_k^diag = 0 forced by A₁ = −2pH, held at 1e-15 to
1e-11 for k = 2 and 1e-15 to 9e-10 for k = 4 (the k = 4 gate is what sets
the N = 256 float-assembly cap; see `exact_pair_moments.py`). Figure: `../syk-self-averaging/pair_moments.png`.

## Scope

The chord moments are i ≠ j flavor-*pair* amplitudes at strict N = ∞. The
ED bench measures the full channel, which at dense-ED sizes is
flavor-diagonal-dominated but not overwhelmingly so: the pair share of μ₀
is −0.6% at N = 12, −17.2% at N = 14 and −34% at N = 16 — opposite in
sign, and *growing fast* with N. So the bench does **not** probe the
chord cancellation, and the growth of that share is exactly why the pair
theorem cannot be lifted to the full channel. What the bench does establish is that along fixed p with λ
falling, the raw width w₂ *rises*: 0.30 → 1.20 over N = 10 → 20 at
p = 4, n = 2. The tower channel turns on exactly where the match must be
made. At the near-matched-λ pair (p=4, N=10) ↔ (p=6, N=22) the raw w₂
also rises with N; the approach to the N = ∞ zero is not visible at
dense-ED sizes, and firm fixed-λ asymptotics need N ≫ 2p.

One loose end worth naming: the (N, p) = (8, 4) point shows an exact
per-instance conservation of all moments, absent at (12, 6) and (16, 8).
The mechanism is unknown and it is not a p = N/2 law.

Remaining caveats: the chord evaluation uses the standard leading-chord
matter weight q^(1/p), but the conserved-channel zeros hold at any
weight, so that is not the source. The Dirac CP-resolved data
(`../syk-self-averaging/moments.json`) awaits complex-SYK chord
machinery, and the Dirac ↔ Majorana channel map is not one-to-one.

## What it means — and what it does not

The theorem is about the **pair** amplitude. The full physical channel is

    ⟨AₙAₙ⟩ = N·⟨XᵢXᵢ⟩ + N(N−1)·⟨XᵢXⱼ⟩,   Xᵢ = ψᵢ(ad_H)ⁿψᵢ,

and the two pieces are the *same order*: the per-pair amplitude is smaller
than the diagonal one by exactly the 1/N that the N(N−1) multiplicity
restores. So conserving the pair part does **not** conserve the channel.

Evaluated exactly at finite N (`exact_pair_moments.py`, both flavor
patterns, no ED), along λ_chord = 2:

| (p, N) | (4,16) | (6,36) | (8,64) | (10,100) |
|---|---|---|---|---|
| w₂ = μ₂/μ₀ of the **full** A₂ channel | 1.044 | 1.118 | 1.164 | 1.194 |

It rises toward ≈ 1.2 — it does not fall to zero. **The full singlet
channel propagates at N = ∞ for n ≥ 2, and its N = ∞ weight is carried
entirely by the flavor-diagonal part.** The control is exact: at n = 1 the
two pieces cancel to machine precision (w₂ = 3e-16), which is A₁ = −2pH.

So the sharp statement this module establishes is a *flavor* statement,
not an invisibility statement:

> At N = ∞ the genuinely two-flavor (i ≠ j) part of every singlet
> bilinear channel is conserved and carries no spectral weight. Whatever
> the tower has at N = ∞ is diagonal.

An earlier reading of this result — that the tower is a 1/N effect
invisible at every order of the double-scaled expansion, so the spectrum
match must be made at fixed p, λ → 0 — **does not follow** and is
withdrawn. The separate fixed-p argument still stands on its own footing:
`syk-self-averaging/exact_wick.py` shows the emergent gauge symmetry is
only polynomially accurate there, whatever the tower does.
