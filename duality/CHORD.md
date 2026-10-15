# The N = ∞ chord engine, and what it says about the tower

*Code: `chord.py` (the engine and the spectral scan), `chord_moments.py`
(singlet moments), `chord_charges.py` (the conserved-channel theorem),
`diagonal_chord.py` (the diagonal law), `full_channel_chord.py` (the
full A₂ moments), `singlet_tower.py` (the full channel resolved into
matter chords, every n), `chord_fold.py` (the fold at N = ∞). The ED side is `syk-self-averaging/mn_majorana_bench.py`. Every
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

Four independent checks hold the engine up (`test_chord.py`): three
legs that share nothing with the transfer-matrix bookkeeping, and the
finite-N ED comparison below them.

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
of the proof. It parallels the known large-q cancellation of the bilinear
tower at leading order ([LARGEQ.md](LARGEQ.md)). The parallel is
qualitative: the theorem holds at fixed λ for the pair amplitude, and how
the λ → 0 limit relates to the large-p expansion at fixed p is not
established.

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

`pair_moment_laws.py` then replaces the fits with certified closed forms,
and `pair_moment_ladder.py` carries them to the K = 12 word frontier
(k + 2n ≤ 12, so μ₆ through n = 3 and μ₂ through n = 5). Every law has
deficit K/2 + 1: μ₂ ~ N^−(n+2) exactly through n = 5. The leading
coefficient is shared across the n = 1 → 2 step at every rung,
A(n=1, k) = A(n=2, k−2), and across no deeper step. The k = 2
coefficients are A(n=1) = −16p²(p−1), A(n=2) = −64p³(p−1)(4p−5), and
A(n=3) = −256p⁴(p−1)(34p² − 108p + 91); the other K = 8 family is
A(1,6) = A(2,4) = −256p⁴(p−1)(34p² − 93p + 64). The double-scaling
constant runs c₂ = 28.92 → 31.65 over p = 4..48, exact, and its limit
is the closed form c₂(λ) = 4λ(1−q)²(5+3q) = 32.3341891 at λ_chord = 2
(`../syk-self-averaging/double_scaling_law.py`): the Δψ → 0 slope of
the diagonal channel below, on every double-scaling line, one rung up
too (c₃(λ) = 4λ(1−q)³(35+56q+35q²+10q³)).
(The float scan above is n ≤ 3 with μ₀/μ₂/μ₄; the exact-Fractions
path behind the laws has no such cap short of K = 12.)

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

| (p, N) | (4,16) | (6,36) | (8,64) | (10,100) | (20,400) | (40,1600) |
|---|---|---|---|---|---|---|
| w₂ = μ₂/μ₀, **full** A₂ channel | 1.044 | 1.118 | 1.164 | 1.194 | 1.255 | 1.286 |

It rises toward 1.3169. It does not fall to zero. The first four values
are recomputed on demand and pinned in
`syk-self-averaging/test_exact_pair_moments.py`; the series to p = 40
is in `full_channel_chord.json`. The committed `exact_pair_moments.json`
stores the pair channel only.

The control is exact. At n = 1 the two parts cancel to machine precision,
w₂ = 3e-16. That is A₁ = −2pH.

The diagonal channel has its own N = ∞ chord law
(`diagonal_chord.py`). An H-chord weighs q^Δψ if and only if it crosses
the four-mark set an odd number of times. The law is derived from the
finite-N sign law, whose lone-chord mark factor is exactly (1 − 2p/N).
There is no pairing sum and no matter-matter chord. Its Δψ → 0 slopes
are closed forms, λP_K(q) with K the word length, and the diagonal A₂
width at λ_chord = 2, Δψ → 0, is P₆/P₄ = (1−q)(5+3q) = 4.674.

The full channel has its own N = ∞ law too (`full_channel_chord.py`),
and it is the physical one. Sum the flavor first: the exact identity
A₂ = 4ΣᵢHᵢ² (`../syk-self-averaging/conserved_points.py`) makes A₂ the
operator H² with each pair of H-chords weighted by the number of indices
they share, |c ∩ c′|, and in the double-scaled limit those intersections
are independent Poisson(λ/2) variables. Every full-channel word is a
finite matching sum with the usual q per crossing and the Poisson
intersection moments on the two weighted pairs. All of it is polynomial:

    μ₀ = 4λ[λ(1+q) + 2(1−q)],   μ₂ = 16λ(1−q)[1 − q² + 2λq],

μ₄ likewise, so w₂ = 1.316906 at λ_chord = 2 (μ₀ = 32 exactly), which
the exact finite-N series above reaches to 1e-5. At small λ,
w₂ = 4λ − 4λ² + … and the kurtosis → 7/2: a dimension-2 matter chord,
as ψ∂²ψ should be.

That closed form is what the 4.674 and the 1.3169 were both shadows of.
Per flavor, N μ^diag and N(N−1) μ^pair are each O(√N) on a
double-scaling line and cancel identically at that order — the pair
channel's double-scaling constant is the diagonal channel's Δ-slope —
and the O(1) remainder belongs to neither half.

### The full channel, resolved

Those moments still carry something that does not move. Subtracting the
mean takes the identity out of A₂. It does not take out H², and A₂
overlaps H²:

    ⟨A₂ (H² − 1)⟩ = 2λ(1 − q),      ⟨(H² − 1)²⟩ = 1 + q.

Remove that component and what is left is, word for word, a matter chord
of dimension 2 (`singlet_tower.py`; this one is a theorem):

    A₂ − ⟨A₂⟩ = [2λ(1−q)/(1+q)] (H² − 1)  +  ν O₂,
    ν² = 8λ(1 − q² + 2λq)/(1+q).

O₂ is an ordinary matter chord. An H-chord that crosses it weighs q².
So every moment of the part that moves is the O_Δ closed form at Δ = 2,
and in particular

    w₂^dyn = μ₂/ν² = 2(1 − q²),    exactly, at every λ.

At λ_chord = 2 that is 1.963369. The 1.316906 above is the same number
diluted: 32.9% of μ₀ there is the conserved H² piece. The "Δ = 2 at small
λ" remark above is the small-λ face of an identity that holds everywhere.

The proof is short. In a two-insertion word, sort the chord diagrams by
whether both chords of one A₂ end on the other A₂. If they do, the pair
is a bundle that every separating H-chord crosses twice: a Δ = 2 matter
chord. If they do not, swapping which of the two adjacent slots carries
which chord flips their mutual crossing and nothing else, and that
involution turns the remaining diagrams into those of a plain H² − 1,
times a number.

Three independent checks, none of which use the proof:

- The exact finite-N full channel with its H² component removed (exact
  rationals) extrapolates in 1/p to w₂^dyn = 1.96336 on λ_chord = 2
  and 1.7293 on λ_chord = 1, against 1.963369 and 1.729329.
- Dense ED, per instance, with the whole commutant of H projected out.
  At p = 6, N = 16–20 the dynamical width is 2.03, 2.03, 2.01 against the
  N = ∞ value 2.00. The raw width there is 0.29–0.41. So the bench does
  reach the N = ∞ law, once the conserved weight is out of the
  denominator.
- The spectral density built from the decomposition and the closed-form
  matter vertex integrates back to the moments.

Two things in this note were the dilution, not dynamics. The raw width's
maximum near λ_chord ≈ 1 at N = ∞: the dynamical width is monotone and
saturates at 2. And the fixed-p maximum at N = 24 (λ_chord = 4/3): with
the H² component removed the width falls monotonically along fixed p,
and in units of J the mean-square frequency rises monotonically to
8p(p−2).

### The tower above n = 2

The same construction runs for every n. An exact operator identity
(`../syk-self-averaging/conserved_points.py`, any N, per instance)
writes Aₙ as mixed words in H and the flavor parts Hᵢ. In the
double-scaled limit only the pairs survive, and Aₙ is assembled from the
blocks Ωₖ = Σᵢ Hᵢ Hᵏ Hᵢ. Each block is solved (`singlet_tower.py`):

    Ω₀ = O₂
    Ω₁ = [q/(1+q²)] {H, O₂} + O₃
    Ω₂ = a(H²O₂ + O₂H²) + b H O₂ H + c O₂ + [q/(1−q+q²)] {H, O₃} + O₄
    Ω₃ = … + O₅,    Ω₄ = … + O₆

modulo polynomials in H, with Oₘ an un-normalized matter chord of
dimension m. So

> at N = ∞ every Aₙ computed (n ≤ 6) is a polynomial in H plus matter
> chords of dimension 2 through n, each dressed by a polynomial in H
> acting from the left and the right. One new chord, of dimension n,
> enters at level n.

Read the weights before reading that as a tower of independent
particles. The new chords are faint. At λ_chord = 2 the dynamical weight
of A₃ is 85% O₂ and 15% O₃; A₄ is 83%, 15% and 1.6% O₄; A₅ carries 0.2%
in O₅ and A₆ 0.02% in O₆. Across λ_chord from 0.05 to 4 the weight of
A₄, A₅, A₆ outside O₂ and O₃ stays between 0.9% and 5.4%. The norm of Oₘ
goes as λ^(2⌊m/2⌋): λ², 2λ², (52/5)λ⁴, about 12.5λ⁴, then λ⁶.

Norms are not shares, though, because the dressings carry powers of the
frequency. The dressing of O₂ inside A₃ is (3/2)[H, A₂], a total
derivative, plus an anticommutator admixture −(1−q)²/(3(1+q²)) of it.
A derivative costs a factor of the width, which is O(λ). So toward small
λ the even chords fade out of every level above the second, linearly in
λ, and the odd ones stay: A₃ and A₄ tend to O₃ alone, A₅ to 99.5% O₃
plus 0.53% O₅, A₆ to 99.2% O₃ plus 0.80% O₅ (exact arithmetic down to
λ = 1/1280). A₂ is O₂ at every λ. The dimensions that survive the
semiclassical limit are 2, 3, 5: the pattern of the large-q SYK bilinear
spectrum h = 2, 3, 5, 7, … That is an observation on four levels, not a
derivation.

What these continua look like in a Rindler frame is in
[HORIZON.md](HORIZON.md): as λ → 0 each chord of dimension Δ becomes the
quasinormal ladder |Γ(Δ + iω/2𝒥)|² of a scale-invariant operator, with no
massive component.

n = 2 is proven. n = 3..6 are exact identities verified in computer
algebra: each kernel is solved linearly from cross-correlators and every
remaining two-point word then has zero residual, symbolically for
k ≤ 2 and in exact rational arithmetic for k = 3, 4. The closed forms
for μ₂ of A₃ and A₄ are gated by the exact finite-N engine, to 1.5e-5
and 2e-3; the cross moment ⟨ad²(A₄)A₂†⟩ to 7e-6. A proof for general n
is open. Figure: `singlet_tower.png`.

So the sharp statement is:

> At N = ∞ the coherent i ≠ j part of every singlet bilinear channel is
> frozen. The full channel is not: it is a conserved polynomial in H
> plus matter chords of dimension 2..n. The split of its weight into
> "diagonal" and "pair" is bookkeeping that cancels at leading order.

The tower propagates. It is not invisible, and its weight is not
"all diagonal" (ERRATA E1, E2). What stands is the spec. A dictionary has to
produce, from the bulk side, a conserved H² admixture of the stated
size, a dimension-2 continuum with the stated norm, and above it chords
of dimension 3, 4, 5 with norms that are closed forms in λ. None of that
is a free parameter, and none of it is a discrete line: each sector is
the smooth matter-chord continuum of the spectral scan above.

The fixed-p question stands on its own footing, separate from this
theorem: `syk-self-averaging/exact_wick.py` finds the emergent gauge
symmetry only polynomially accurate there, whatever the tower does.

## Scope and loose ends

The chord moments are i ≠ j pair amplitudes at strict N = ∞.

The ED bench measures the full channel. At dense-ED sizes that channel is
diagonal-dominated, but not overwhelmingly. The pair share of μ₀ is
−0.6% at N = 12, −17.2% at N = 14, and −34% at N = 16. It is opposite in
sign, and it grows fast with N. So the bench does not probe the chord
cancellation. That growth is also why the pair theorem cannot be lifted
to the full channel.

What the bench does find: along fixed p with λ falling, the raw width w₂
rises from 0.30 to 1.20 over N = 10 → 20 at p = 4, n = 2. That window
ends before the peak. The exact engine continues the same channel
(`../syk-self-averaging/full_channel_laws.py`): the full width tops out
at N = 24 (λ_chord = 4/3) and then dies as 8p(p−2)/N (the coefficient
fit on p = 2, 4, 6 and confirmed out of sample at p = 8, 10, 12). In physical units that is an N-independent
mean-square frequency 8p(p−2)·J², a finite O(J) width all the way into
the corner. The bench's rising width was the near side of that peak
(ERRATA E4). The peak itself is the
conserved H² weight draining out of the denominator (85% of μ₀ at
N = 10, 19% at N = 24, 0.01% at N = 1024): with that component removed
the width has no maximum. It falls from about 2 once N ≥ 12, tracking
the N = ∞ law 2(1 − q²) at the small-N end, and its physical value
rises monotonically to the same 8p(p−2). At the near-matched-λ
pair (p=4, N=10) ↔ (p=6, N=22) the raw w₂ also rises with N. The
approach to the N = ∞ pair zero is not visible at dense-ED sizes. Firm
fixed-λ asymptotics need N ≫ 2p.

The (N, p) = (8, 4) point shows exact per-instance conservation of all
moments, absent at (12, 6) and (16, 8). It is not a p = N/2 law, and the
mechanism is now proven (`../syk-self-averaging/conserved_points.py`):
the spin-trace map Σᵢψᵢ X ψᵢ weighs Clifford grade g by (−1)^g(N − 2g),
every power of H at N = 8 has grades {0, 4, 8} only, so every Aₙ is a
polynomial in H with central coefficients. Two exact identities fall
out for every (N, p): A₁ = −2ΣHᵢ and A₂ = 4ΣHᵢ², Hᵢ the part of H
containing flavor i — the second is the operator behind the full-channel
chord law above.

Two more caveats. The chord evaluation uses the standard leading-chord
matter weight q^(1/p), but the pair-channel zeros hold at any weight, so
the weight is not the source. On the Dirac side the pair half now
transfers (`dirac_chord.py`): the bare-fermion census gives a lone
crossing of exactly 1 − p/Nc, the Majorana sign census at N = 2Nc, so
the theorem holds verbatim at (q_dirac, Δ = 1/p). The Dirac CP-resolved
data (`../syk-self-averaging/moments.json`) still waits on the diagonal
half and the CP resolution, and the Dirac ↔ Majorana channel map is not
one-to-one.
