# The N = ∞ chord engine, and what it says about the tower

*Code: `chord.py` (the engine and the spectral scan), `chord_moments.py`
(singlet moments), `chord_charges.py` (the conserved-channel theorem),
`diagonal_chord.py` (the diagonal law), `chord_fold.py` (the fold at
N = ∞). The ED side is `syk-self-averaging/mn_majorana_bench.py`. Every
verdict below is asserted from a committed JSON in `test_chord*.py`.*

## What the bulk is

The exact N = ∞ solution says the bulk is a length.

The state |n⟩ means the bulk is n chords long. The chord Hamiltonian T
does one thing. It lengthens the bulk or it shortens it:

    T|n⟩ = |n+1⟩ + [n]_q |n−1⟩,    [n]_q = (1 − qⁿ)/(1 − q)

Open a new chord, or close one of the n open ones. Closing the k-th one
means crossing k−1 others, so the closing weight sums to [n]_q.

A particle in that bulk is a chord stretched across it. The bulk grows
and shrinks, the particle gets dragged along, and that is what puts
weight at nonzero frequency in a correlator. That is what "propagates"
means here.

Keep this picture. Everything below asks which probes are particles, and
which ones are the bulk itself.

## The engine

Double-scaled SYK at N = ∞ (λ = 2p²/N fixed) is solved by chord
combinatorics (Berkooz et al.
[1811.02584](https://arxiv.org/abs/1811.02584)). Conventions follow the
review [2407.09396](https://arxiv.org/abs/2407.09396). Both are
equation-checked.

`chord.py` writes that solution as finite sparse linear algebra. It has
the transfer matrix, matter chords of weight q^Δ, and the one- and
two-matter sectors. Those sectors evaluate any operator ordering on any
complex-time contour. Crossed and Euclidean-folded contours are included.
There are no special functions.

Three independent legs hold the engine up (`test_chord.py`):

1. A brute-force chord-diagram enumerator. It shares nothing with the
   transfer-matrix bookkeeping. It agrees to machine precision on every
   topology, fermionic signs included.
2. The published closed forms. The θ-integral 2-point function (review
   Eq. 2.30), μ₂ = 2(1−q^Δ), and the Touchard–Riordan moments.
3. A chord-free check. The exact finite-N crossing average is a
   Krawtchouk sum. It approaches q = exp(−2p²/N) at the measured
   O(p⁴/N²) rate.

The units bridge is exact and tested: σ_H² = C(N,p) p!/N^(p−1), which
gives β_chord = βJ · exp(−λ/8)/√λ at N = ∞. Against finite-N ED the
chord thermal 2-point function agrees to 0.6–1.6%. The large-q closed
forms sit at a 30–55% floor at the same sizes.

Stated caps. At most two matter chords are open at once, which covers
every 2-point and 4-point ordering. Folded contours below λ ≈ 0.5 exceed
float headroom. The λ = 0.3 point fails truncation-doubling. Per-row
drifts live in `chord_fold.json`.

## Is the spectrum a tower or a continuum?

It is a continuum.

`chord.py --spectra` computes the exact spectral functions S_Δ(ω) at
T_B = ∞. It covers Δ ∈ {1, ½, ¼} and the Δ → 0 length channel over
λ = 2.0 → 0.05. There is a fine window at the semiclassical scale, where
a 't Hooft tower would sit.

- **No interior maxima. No nodes. Every channel, every λ.** The exact
  spectral function converges to the smooth sech^2Δ transform linearly in
  λ. For Δ = ½ the deviation falls 0.375 → 0.006, which is about 0.13 λ.
  The worst channel is Δ = 1, at about 0.25 λ. The mass-encoding nodes of
  [BOOST.md](BOOST.md) are absent too.
- The chord-number growth ⟨n(t)⟩ is the de Sitter length. It approaches
  −(2/λ) ln sech(J_c t) at the same linear rate. That is new exact-in-λ
  data.

## The pair channel is frozen

Now take the paper's own channels, the singlet bilinears

    Aₙ = Σᵢ ψᵢ (ad_H)ⁿ ψᵢ.

You build Aₙ out of the Hamiltonian. You hit a fermion with H, n times.
`test_mn_majorana.py` verifies the binomial assembly as a per-instance
operator identity against dense ED.

`chord_moments.py` measures the moments. `chord_charges.py` proves the
result:

> **Every spectral moment μ_k, k ≥ 1, of the i ≠ j pair amplitude of
> every Aₙ vanishes at N = ∞.** Every λ, every matter weight, both
> orderings, odd k included.

Measured: μ₂ = μ₄ = μ₆ = 0 to machine precision, 1e-17 to 1e-18 relative
to the pre-cancellation gross sum, for n = 1…5, on the λ grid [0.1, 4] at
Δψ = ¼, and across the matter-weight scan Δψ ∈ [0.1, 1] at λ = 1.

The scope is load-bearing. Read "What it means" below before quoting it.

### Why it happens

In a sequential pair word the fermion pair factors through the zero-chord
sector. What is left is a dressed propagator,

    C(m) = close ∘ T_Mᵐ ∘ open.

C(m) is the whole block ψ Hᵐ ψ, seen from the chord line. A three-term
recursion controls it:

    close_j T_M = q^(Δ+j) T close_j + [j] close_{j−1} + (1−q^(2Δ+j)) close_{j+1}

The coefficients are **scalars**. That is the whole result. Both matching
conditions collapse to q-number identities that do not depend on the
chord occupation. So C(m) is the Motzkin-path polynomial in T, with flat
weight q^(Δ+j)T, down weight [j], and up weight 1−q^(2Δ+j).

If those coefficients had kept any occupation dependence, the family
would not close. C(m) would be a genuine operator on the chord line, and
the channel would propagate like anything else.

The chord image of Aₙ is therefore a polynomial in T. A polynomial in the
Hamiltonian commutes with the Hamiltonian. It does not tick. Explicitly:

    Â₁ = −(1−q^Δ) T,        Â₂ = (1−q^Δ)²T² + (1−q^(2Δ))

Â₁ ∝ T is the N = ∞ image of the exact finite-N law A₁ = −2pH. That is a
real check. The n = 1 case lands on something known independently.

### What does propagate

A random operator O_Δ inserts q^(Δn̂). That does **not** commute with T,
because T changes n. So random matter propagates: μ₂ = 2(1−q^Δ), and the
kurtosis goes to 3 + 1/Δ.

The bulk holds particles perfectly well. The discriminator is not the
size of the operator. It is whether the operator's couplings are
correlated with H's. Probes built out of the metric are not particles.

Those closed forms are tested. `test_chord_charges.py` asserts every step
of the proof. This is the exact-in-λ version of the known large-q
cancellation ([LARGEQ.md](LARGEQ.md)), now with its algebraic origin.

## The finite-N approach to the zero

`syk-self-averaging/exact_pair_moments.py` computes the same moments
exactly at finite N. It runs Wick over the H-slot matchings with the
chord-diagram sign law, and evaluates the signed subset sums by a CRT
residue-table dynamic program. There is no ED and no sampling noise. The
18-point scan takes 332 s. The largest row, N = 256, takes about 2 s per
channel. Scope: n ≤ 3, moments μ₀/μ₂/μ₄, words up to K = 10 H-nodes.

Along **fixed p = 4** every μ_k with k ≥ 1 falls monotonically to zero.
The fitted tail slopes are 2.96, 3.90 and 4.82 for n = 1, 2, 3. Those
slopes still drift upward at N = 256. Read them as consistent with
μ₂ ~ N^−(n+2), not as a measurement of the exponent. The exponent itself
is proven separately, below. (μ₀ of the n = 2 channel is not monotone. It
changes sign near N = 12, which is the exact version of an effect first
seen in ED.)

Along the **double-scaling line λ_chord = 2**, |μ₂| ≈ c/(pN). For n = 2,
c = 28.9, 28.9, 29.4, 29.8 over (p, N) = (4,16)…(10,100). For n = 3,
c = 203, 200, 202, 205. The n = 2 constant drifts slowly upward. The
n = 3 one is flat to 2% and not monotone.

`pair_moment_laws.py` then replaces the fits with certified closed forms.
μ₂ ~ N^−(n+2) holds exactly through n = 4. The coefficients are
A(n=1) = −16p²(p−1) and A(n=2) = −64p³(p−1)(4p−5). The double-scaling
constant runs c₂ = 28.92 → 30.53 over p = 4..16.

Validation of the finite-N engine: literal Isserlis enumeration (exact
match), dense ED, the chord limit at the Krawtchouk N⁻² rate, and the
per-row gate (N−1)μ_k^pair + μ_k^diag = 0 forced by A₁ = −2pH. That gate
holds at 1e-15 to 1e-11 for k = 2, and 1e-15 to 9e-10 for k = 4. The
k = 4 gate is what sets the N = 256 float-assembly cap. Figure:
`../syk-self-averaging/pair_moments.png`.

## What it means, and what it does not

The theorem is about the **pair** channel. The full physical channel has
two parts:

    ⟨AₙAₙ⟩ = N·⟨XᵢXᵢ⟩ + N(N−1)·⟨XᵢXⱼ⟩,    Xᵢ = ψᵢ(ad_H)ⁿψᵢ

The pair channel is the coherent one. Two different constituents, which
is what a meson looks like. The diagonal channel is an incoherent pile of
N single-flavor pieces.

The two parts are the **same order**. The per-pair amplitude is smaller
than the diagonal one by exactly the 1/N that the N(N−1) multiplicity
restores. So freezing the pair part does not freeze the channel.

Evaluate both exactly at finite N, along λ_chord = 2:

| (p, N) | (4,16) | (6,36) | (8,64) | (10,100) |
|---|---|---|---|---|
| w₂ = μ₂/μ₀, **full** A₂ channel | 1.044 | 1.118 | 1.164 | 1.194 |

It rises toward ≈ 1.2. It does not fall to zero. All four values are
recomputed on demand and pinned in
`syk-self-averaging/test_exact_pair_moments.py`. The committed
`exact_pair_moments.json` stores the pair channel only.

The control is exact. At n = 1 the two parts cancel to machine precision,
w₂ = 3e-16. That is A₁ = −2pH.

The diagonal channel now has its own N = ∞ chord law
(`diagonal_chord.py`). An H-chord weighs q^Δψ if and only if it crosses
the four-mark set an odd number of times. The law is derived from the
finite-N sign law, whose lone-chord mark factor is exactly (1 − 2p/N).
There is no pairing sum and no matter-matter chord. The diagonal A₂ width
at λ_chord = 2 is w₂ = 4.67(9) at Δψ → 0. The full-channel width is
≈ 1.31(5), by Richardson extrapolation of the exact finite-N series
above. The gap between 4.67 and 1.31 is a large diagonal-against-pair
cancellation. Its remaining chord-side piece is the pair channel's 1/N
coefficient, supplied exactly at fixed p by `pair_moment_laws.py`.

So the sharp statement is about flavor:

> At N = ∞ the coherent i ≠ j part of every singlet bilinear channel is
> frozen and carries no weight. Whatever the tower has at N = ∞ is
> diagonal.

The tower propagates. It is not invisible. All of its N = ∞ weight sits
in the channel with the least meson-like structure. In 't Hooft's model a
meson is a color singlet, and the whole content of "singlet" is coherence
across constituents. That is the question worth putting to the authors.

**A withdrawn reading.** An earlier version of this note concluded that
the tower is a 1/N effect, invisible at every order of the double-scaled
expansion, and that the spectrum match must therefore be made at fixed p,
λ → 0. That does not follow from this theorem, and it is withdrawn. The
fixed-p argument still stands on its own footing:
`syk-self-averaging/exact_wick.py` finds the emergent gauge symmetry only
polynomially accurate there, whatever the tower does.

## Scope and loose ends

The chord moments are i ≠ j pair amplitudes at strict N = ∞.

The ED bench measures the full channel. At dense-ED sizes that channel is
diagonal-dominated, but not overwhelmingly. The pair share of μ₀ is
−0.6% at N = 12, −17.2% at N = 14, and −34% at N = 16. It is opposite in
sign, and it grows fast with N. So the bench does not probe the chord
cancellation. That growth is also why the pair theorem cannot be lifted
to the full channel.

What the bench does find: along fixed p with λ falling, the raw width w₂
rises from 0.30 to 1.20 over N = 10 → 20 at p = 4, n = 2. The tower
channel turns on where the match has to be made. At the near-matched-λ
pair (p=4, N=10) ↔ (p=6, N=22) the raw w₂ also rises with N. The approach
to the N = ∞ zero is not visible at dense-ED sizes. Firm fixed-λ
asymptotics need N ≫ 2p.

One loose end. The (N, p) = (8, 4) point shows exact per-instance
conservation of all moments. That is absent at (12, 6) and (16, 8). The
mechanism is unknown, and it is not a p = N/2 law.

Two more caveats. The chord evaluation uses the standard leading-chord
matter weight q^(1/p), but the pair-channel zeros hold at any weight, so
the weight is not the source. The Dirac CP-resolved data
(`../syk-self-averaging/moments.json`) waits on the complex-SYK chord
layer, and the Dirac ↔ Majorana channel map is not one-to-one.
