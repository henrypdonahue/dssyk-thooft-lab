# Road map — stress-testing DSSYK∞ ↔ 't Hooft ([arXiv:2607.05678](https://arxiv.org/abs/2607.05678))

What the paper claims. What this repo has settled. What is left. The
literature snapshot is current as of 2026-09-09; the arXiv versions and
the citation count were re-checked on 2026-10-14. Withdrawn and
corrected claims are collected in [ERRATA.md](ERRATA.md).

---

## Part I — The paper, honestly assessed

*Holograms and Standard Models* (Miyashita–Sekino–Susskind, July 2026) is
a short talk-style manifesto. It is not a technical paper, and it does
not pretend to be. Here is its chain of claims, with the paper's own
confidence levels attached.

1. **Assumed.** Their words: "unproved but we will assume it." DSSYK∞ =
   JT-de Sitter static-patch holography.
2. **Quasi-derived** (§3.5, from two example diagrams, with the
   systematics deferred to [arXiv:2501.09423](https://arxiv.org/abs/2501.09423)). Perfect self-averaging,
   Var⟨W⟩ ≈ p!/Nᵖ < e^(−a√N). Hence *emergent* O(N)/U(N) symmetry of
   single disorder instances in the flat-space limit.
3. **Argued by counting.** The entropy match S = N log 2 gives
   dim matter = N. So the bulk dual is the 't Hooft model: U(N) QCD₂ plus
   QED, one fundamental flavor, at renormalized quark mass m_q² → 0, with
   the dictionary N_ym = N_syk, ᾱ = ḡ²N = p², ḡ² = λ = p²/N.
4. **Qualitative only.** Their words: "We don't have a quantitative
   calculation of the DSSYK spectrum." The singlet tower
   Mₙ = ψ†(dⁿ/dtⁿ)ψ maps to 't Hooft mesons, with CP = (−1)ⁿ⁺¹, and with
   n = 0 (photon = Q) and n = 1 (graviton = H) non-dynamical. The paper
   names the quantitative spectrum match as "the best test of the
   duality." **No published work has run it** (as of 2026-10-14).
5. **A logical inversion, not a computation** (§5, "Holographic
   Tuning"). At T_B = ∞ the entropy is exactly log dim H, so the vacuum
   energy gets no corrections, and the coupling ḡ²(N) becomes a *derived*
   quantity. Its finite-N form is only a sketched figure.

### Two structural weaknesses the paper does not address

**The boost-against-mass gap.** The SYK Hamiltonian is the Rindler
*boost* generator. A fixed-mass particle has a *continuous*
boost-frequency spectrum. How discrete 't Hooft masses are supposed to
appear as features of Mₙ correlators at T_B = ∞ is an unstated piece of
dictionary, and every proposed test presupposes it, including the
authors' own.

`duality/BOOST.md` calibrates that on the exactly solvable case. A single
mass produces a continuum with mass-encoding *nodes*, not a peak. And
`duality/boost_inverse.py` prices the gap: on that kernel, moment
invariants carry regulator leverage ≈ 1, so with the state map unfixed
a single invariant carries no mass digits at all, while node positions
are exactly regulator-independent (under any multiplicative regulator)
and still transfer half the data digits. So the
moment route — the one the Mₙ correlators naturally supply — is blocked
outright until the authors supply that map, and the test survives only
through node-level data. (Free scalar, 1+1, single particle, sharp
worldline, e^(−εω) regulator: the calibration bounds what shape data
can show; it does not fix the map.)

**Does the toy even have the disease?** QCD₂ is super-renormalizable. Its
only mass renormalization, m_q² = m_bare² − g²/π, is finite. So the 4d
fine-tuning drama of §5 is transplanted into a model where it barely
exists. And the vacuum-energy/entropy relation ε = 3/(8S) is only ever
written in D = 4, never derived for the JT-dS where the paper works.

`thooft-target/condensate.py` measures the 2d counterpart. It is a single
finite number: ε(0) = −0.0407 Nc g².

### The regime problem

The regime where the spectrum match lives is **fixed p, λ → 0**. That is
*not* the double-scaling regime where self-averaging is
super-exponential.

At fixed p the suppression is only the power law N⁻ᵖ. So gauge-symmetry
emergence is quantitatively weakest exactly where the match has to be
made. This repo makes that a printed, tested table rather than a remark.

---

## Part II — What is settled

### The yardstick is real

An independent re-derivation of the whole Galerkin assembly. A
re-verification of the hypersingular Chebyshev-U identity by regularized
quadrature. And a check of the anchor constants against FLZ's published
tables ([arXiv:0905.2280](https://arxiv.org/abs/0905.2280)), digit for
digit.

The worst deviation from FLZ over n ≤ 29 is 1.0e-12, and that is FLZ's
own last printed digit, not this solver: there is now a 40-level
reference table good to 17–21 digits
(`reference_spectrum_polished.json`: exact-rational Rayleigh quotients
of double-precision Ritz vectors, with the measured K⁻⁶ truncation law
extrapolated away), and against it the double-precision solver is at
9e-14 for every level while FLZ's printed values are off by up to one
unit in their twelfth decimal. The older 24-level mpmath table agrees
with it to at least its own claimed digits at every row. All sixteen of
FLZ's Table 3 sum rules (s = 1 through 8, both parities) come out at the
table's last digit, about 2e-12, once the transcribed Eq. 1.6 tail is
used; 7ζ(3) and 2 among them.

Off the duality point the Chebyshev basis hardwires the α = 0 endpoint
exponent ½, while the true exponent solves πβ cot(πβ) = −α. So
convergence there is only algebraic. `jacobi_solver.py` fixes that with a
matched-exponent basis [x(1−x)]^β: the ground state is below 1e-11 by
K = 32, and the LM sum rules come out to about 3e-8 at s = 2 and 4e-12
at s ≥ 4.

Eigenfunctions are checked (Parseval 3e-6, weak-form Rayleigh quotient
~1e-12 on tanh-sinh nodes, cross-solver overlap > 0.9999). FLZ's spectral-determinant constants d±
are consumed too (`flz_dconstants.py`): their Eq. 4.21 products
converge to d₊ = 0.9631784563984(1), d₋ = 0.43358263983325(1), inside
the 5e-13 rounding of both 12-digit values and one to two digits past
them, with the exact ratio √2/π met to 2e-15. Pinning their definition corrected a
mischaracterization — d± are spectrum functionals, not wavefunction
data — so the external eigenfunction-anchor slot is genuinely empty.
The decay constants are new data. No published wavefunction tables
exist.

### Self-averaging, exactly

`exact_wick.py` gives Var(B_jk) = 4σ⁴C(N−2,p−1) in closed form at any
(N, p). It is checked by literal contraction enumeration, by ED, and by
Monte Carlo.

In double scaling ln rms = −¼√(λN) ln N, which beats e^(−a√N) for every
a, scanned to N = 1600. The paper's Eq. 3.28 is not just right. It is
conservative. The singlet rate of Eq. 3.27 is exactly 2/C(N,p).

A second anchor, exact E[m₄] and Var(m₄), gives
relVar(m₄) = 4·relVar(m₂)(1+o(1)). The Dirac/U(N) counterpart is in
`exact_wick_dirac.py`, with rms ~ Nc^(−(p−1)/2).

**Worth telling the authors.** The "systematic derivation" that §3.5
defers to [arXiv:2501.09423](https://arxiv.org/abs/2501.09423) **is not
there**. That paper never mentions self-averaging. Its counting rule
(Eq. 4.26) is magnitude-only and part-empirical (their §6.2), and its n
≪ p approximation fails precisely in the diagram-B regime, with error
~e^p. That is harmless for the bound and fatal for prefactors.
`exact_wick.py` appears to be the only complete prefactor-exact
derivation. It agrees with both papers at leading order, and it adds the
exact pinning factor 4p(N−p)/N(N−1).

### The dictionary needs four corrections

Each one is forced, and each one is tested (`duality/dictionary.py`,
`syk-self-averaging/dirac.py`).

**M₁ = −(p/2)iH, not H.** It follows from anti-Hermiticity. Q̇ = 0 forces
M₁† = −M₁, so M₁ cannot equal the Hermitian H. Verified to 1.3e-15.

**M₂ = −Σċ†ċ.** The paper's + sign is wrong.

**CP is the unitary particle-hole map.** M₁ ∝ iH rules out an antiunitary
one.

**λ = p²/N against λ_chord = 2p²/N.** This factor of two is internal to
the Susskind corpus. [2501.09423](https://arxiv.org/abs/2501.09423) Eq.
3.6 uses 2p²/N.

Raw Mₙ for n ≥ 2 are CP mixtures with 15–65% wrong-channel content. So
CP = (−1)^(n+1) only lives in CP-resolved correlators.

### The singlet channels at N = ∞, and where the tower lives

This is the headline result. It has two halves, a conserved pair
amplitude and a propagating full channel, and the scope of each matters.

**The pair theorem.** The chord image of the i ≠ j *pair* part of every
Aₙ = Σψ(ad_H)ⁿψ is a polynomial in the chord Hamiltonian T. So every
spectral moment of every Aₙ pair channel vanishes at N = ∞: every λ,
every matter weight, both orderings, odd k included.

The proof is a three-term recursion whose coefficients are scalars. Both
matching conditions reduce to q-number identities that do not depend on
the chord occupation. The explicit charges are Â₁ = −(1−q^Δ)T, which is
the N = ∞ image of the exact law A₁ = −2pH, and
Â₂ = (1−q^Δ)²T² + (1−q^(2Δ)).

Random O_Δ operators insert q^(Δn̂), which does *not* commute with T, so
they propagate. That contrast is what makes the result physical rather
than a normalization artifact.

**The scope, which is load-bearing.** This is the i ≠ j *pair* amplitude.
The full channel also has a flavor-diagonal part of the same order. The
per-pair amplitude is smaller by exactly the 1/N that the N(N−1)
multiplicity restores. And that diagonal part does **not** vanish for
n ≥ 2.

Evaluated exactly at finite N, the full A₂ width rises along λ_chord = 2
(1.044 → 1.286 over (p,N) = (4,16)…(40,1600)). At n = 1 the two parts
cancel to 3e-16, which reproduces A₁ = −2pH.

**The full channel at N = ∞ is now an exact chord law of its own**
(`duality/full_channel_chord.py`). Sum the flavor first: the exact
identity A₂ = 4Σᵢ Hᵢ² makes A₂ the operator H² with each pair of
H-chords weighted by the number of indices they share, and in the
double-scaled limit those intersections are independent Poisson(λ/2)
variables. Every full-channel word is then a finite matching sum, and

    μ₀ = 4λ[λ(1+q) + 2(1−q)],    μ₂ = 16λ(1−q)[1 − q² + 2λq],

with μ₄ a polynomial too. At λ_chord = 2 that is μ₀ = 32 exactly and
w₂ = 1.316906; the exact finite-N series extrapolates to 1.31690(1). At
small λ it is w₂ = 4λ − 4λ² + …, kurtosis → 7/2: as λ → 0 the full A₂
channel is a dimension-2 matter chord, as ψ∂²ψ should be.

**And that channel is now resolved, at every λ**
(`duality/singlet_tower.py`, `notes/singlet-tower-note.md`). Mean
subtraction leaves the conserved H² inside A₂. Take it out and the rest
is exactly a matter chord of dimension 2 — a theorem, proved by a
sign-reversing involution on the chord diagrams:

    A₂ − ⟨A₂⟩ = [2λ(1−q)/(1+q)](H² − 1) + ν O₂,
    ν² = 8λ(1 − q² + 2λq)/(1+q),      w₂^dyn = 2(1 − q²).

At λ_chord = 2 the conserved share of μ₀ is 32.9% and the dynamical
width is 1.963369; the 1.316906 is that width diluted. Three checks that
do not use the proof agree: the exact finite-N channel with its H²
component removed (1.96336), dense ED with the whole commutant of H
projected out (2.03, 2.03, 2.01 at p = 6, N = 16–20, against 2.00), and
the spectral density rebuilt from the closed-form matter vertex.

The levels above follow the same pattern. An exact operator identity
writes every Aₙ through the blocks Ωₖ = ΣᵢHᵢHᵏHᵢ, and each Aₙ computed
(n ≤ 6) is a polynomial in H plus matter chords of dimension 2..n with
polynomial dressings — zero residual in computer algebra, a proof only
at n = 2. One new chord enters per level, and it is faint: 14.9%, 1.6%,
0.18%, 0.02% of the dynamical weight at λ_chord = 2 for n = 3..6, with
norm ∝ λ^(2⌊n/2⌋). Toward small λ the even chords fade out of A₃..A₆
linearly in λ and the odd ones stay (A₅ → 99.5% O₃ + 0.53% O₅): the
dimensions that survive are 2, 3, 5, the pattern of the large-q SYK
bilinear spectrum. That last is an observation, not a derivation.

So the tower does propagate at N = ∞. The per-flavor split of its weight
into "diagonal" and "pair" is bookkeeping: on a double-scaling line the
two halves are each O(√N) in the flavor-summed normalization and cancel
identically at that order (their O(1/p) leading pieces are the same
function of λ, see below), and the O(1) remainder, the physical channel,
belongs to neither (ERRATA E2). The theorem does not make the tower a
1/N effect, and it does not by itself place the match at fixed p
(ERRATA E1). The fixed-p question stands separately, on the
self-averaging side.

**The finite-N laws are now exact, out to the K = 12 word frontier.**
`exact_pair_moments.py` measures the approach to the pair zero by Wick
plus a CRT residue-table dynamic program, with no ED.
`pair_moment_laws.py` turns those float regularities into certified
closed forms, and `pair_moment_ladder.py` extends them to K = 12 words.
The organizing structure is the word length K = k + 2n. Every computed
law — 55 of them, k = 2..10, K = 4..12 — has deficit exactly K/2 + 1,
so μ₂ ~ N^−(n+2) is certified through n = 5, μ₄ ~ N^−(n+3) through
n = 4, and μ₆ ~ N^−(n+4) through n = 3. The leading coefficient obeys
A(n=1, k) = A(n=2, k−2) *exactly*, at every rung and every p, while
every deeper step differs (at K = 10, p = 2 the coefficients run
21 : 21 : 14 : 9 × (−2¹⁷) over n = 1..4). So ad-shifting is exact at
leading 1/N order for the first step of the ladder only. The k = 2
family is A(n) = −2^(2n+2) p^(n+1) (p−1) Q_n(p) with Q₁ = 1,
Q₂ = 4p − 5, and the newly identified Q₃ = 34p² − 108p + 91, certified
at p = 18. The other K = 8 family is
A(1,6) = A(2,4) = −256p⁴(p−1)(34p² − 93p + 64), certified at p = 18 too:
the same prefactor and the same leading 34p² as A(3,2), differing only
below. The double-scaling constant is exact to p = 48,
c₂ = 28.92 → 31.65, and its limit is now a closed form
(`syk-self-averaging/double_scaling_law.py`):

    c₂(λ) = lim |μ₂^pair| pN = 4λ(1−q)²(5+3q),   q = e^(−λ),

32.3341891 at λ_chord = 2, which the exact sequence's 1/p extrapolation
reproduces to 1e-8 relative (and 9.7555 at λ = 1, 77.944 at λ = 4, both
checked on their own lines). It is not 97/3. The same holds one rung
up: c₃(λ) = 4λ(1−q)³(35 + 56q + 35q² + 10q³), 223.648 at λ = 2, which
the n = 3 sequence to p = 32 reproduces to 1e-4. These constants are the
Δψ → 0 slopes of the N = ∞ *diagonal* chord channel, λP_K(q) with K the
word length: the O(Δ) parts of the diagonal and pair halves cancel on
every double-scaling line, the n ≥ 2 shadow of the exact A₁ identity.
The λ-expansion of the closed form fixes the large-p growth of every
exact rational 1/N coefficient, and the N⁻⁵ coefficient turned out to be
a closed form too: B(2,2;p) = 64p³(p−1)(11p³ − 36p² + 29p − 3), fit on
p = 4..12 and reproducing the committed p = 2 value out of sample. What
the chord-side derivation now has to reproduce is that cancellation, and
the fixed-p polynomials Q_n(p).

### The tower at the stretched horizon

[`duality/HORIZON.md`](duality/HORIZON.md). The flat-space limit of
2506.18054 is λ fixed, energies fixed in string units, so the N = ∞
results above describe the singlets in exactly that limit. In those units
(J_MSS = pJ) the chords' ladder spacing is a = 2𝒥 = 2√2 J_MSS: string
scale, the scale of the meson masses and of the stretched horizon of
their §6. As λ → 0 each chord of dimension Δ is the quasinormal ladder
|Γ(Δ + iω/a)|², reached with deviation (Δ/4)λ (measured Δ = 2..6), and
the exact decay rates plateau at Δa. At a Rindler horizon a field of
mass M gives the entire density K_{iω/a}(M/a)², and a ladder is exactly
the superposition of those with mass density ∝ M^(2Δ−1): gapless and
scale-invariant. With the exact decomposition (and the invertibility of
the mass transform), the channels' Rindler mass spectrum is gapless power
laws with Δ = 2, 3, 5, with no gap and no discrete masses. The gapped,
discrete 't Hooft tower gives an entire density instead. That is
consistent with their §6 picture of a deconfined stretched horizon, and a
sign that the meson tower is not visible in the N = ∞ two-point functions
of the Sₙ in their Majorana form (the Dirac full channel is still open).

At fixed p, N → ∞, the paper's calculational regime
(`syk-self-averaging/thooft_limit.py`), the exact leading 1/N
coefficients give the channels' physical moments. A₂ has
⟨ω²⟩ = 8p(p − 2), ⟨ω⁴⟩ = 32p²(p − 2)(7p − 15) and kurtosis
7/2 − 1/(2(p − 2)), approaching the Δ = 2 ladder as p grows. That is the
expected large-p ↔ small-λ correspondence. What is new is the size of the
finite-p corrections, ⟨ω²⟩ at (1 − 2/p) of the ladder value, which is
large at the p where comparisons are made.

### The dS axis

DSSYK∞ **scrambles at every temperature**. The growth coefficient is
a(v) = −1/cos(πv/2), real and negative, tending to −1 at the hologram
point, with |a| growing with β. A dS observer needs the opposite sign, so
the dictionary owes a sign-flipping map.

The Harlow–Zhao Euclidean fold is not that map. At ED sizes the folded
connected fraction grows from 0.365 to 0.786 over N = 8 → 20, while λ
shrinks. At N = ∞ it stays O(1) and grows as λ → 0, while the unfolded
one dies like λ. Neither statement needs a closed form. The mechanism
behind the bad finite-N convergence is measured: folded correlators at ED
sizes track the sample's upper spectral edge at correlation −0.86.

One candidate map does work, semiclassically. The antipodal
half-fake-period shift Δ = π/v flips the sign exactly in the closed form,
and for λ ≤ 1.7 in the exact theory. The same shift makes the Wightman
sector exactly thermal at β_eff = β + β_fake. One map, both of
Cui–Kolchmeyer's dS-observer signatures. Geometrically it is NV's ±iβ/4
dressing, and `duality/NV-DRESSING.md` pins that identification at the
contour level with equation-level citations.

Read the closed-form flip for what it is. Δ enters the growing mode's
phase, so half a period negates it, and the algebra is short. The content
is in the exact-λ boundary and in the deep campaign.

That campaign finds a transient anti-scrambling epoch growing with
ln(1/λ), not eternal dS. It locates the chaos-bound evasion: the phase
winds at constant magnitude, so there is no strip-KMS structure. And it
measures the frame's 2-pt: it is the T = ∞ correlator continued to the
complex time τ_s + it, the large-q |1/cos(𝒥z)|^(2Δ), matched ever
better toward semiclassics, with the antipodal point itself sitting on
the fake-period singularity that the bounded chord spectrum cuts off.
The observer is thermal by the exact detailed-balance identity, not as a
Gibbs state.

### Other results

`freeness.py` exhibits the free-product transition. spec(E + F(t)) lands
on the arcsine law, and the W₁ plateau falls like dim^(−0.9). To our
knowledge that is the first time in a microscopic dS-candidate model.

`qed_campaign.py` closes the gap between Majorana numerics and Dirac
claims. `baryons.py` gives first data on the wholly underived
inter-baryon-force claim. `condensate.py` transcribes the exact
condensate and gives first α → −1 GMOR data.

---

## Part III — What is open, and what each would take

### Next work items, in priority order

1. **Where are the mesons?** Done at N = ∞ for the Majorana channels
   ([duality/HORIZON.md](duality/HORIZON.md)): the orthogonalized tower is
   a pure chord, Õₙ = 4Ôₙ, and as λ → 0 every singlet channel is a sum of
   quasinormal ladders |Γ(Δ + iω/a)|², Δ = 2, 3, 5, whose Rindler mass
   spectrum is gapless (density ∝ M^(2Δ−1)), with no gap and no discrete
   masses. So the paper's small-λ prediction cannot be met by the N = ∞
   two-point functions of the Majorana Sₙ. Next, theory-gated: whether
   1/N corrections open a gap or produce discrete masses. Moments cannot
   decide that (they do not fix analytic structure); it needs the 1/N
   correction to the chord density itself.
2. **The 't Hooft limit, by shape.** Done for A₂, with the A₃ ratio
   (`syk-self-averaging/thooft_limit.py`): at fixed p the channels'
   moments approach the λ → 0 ladders as p grows, with exact closed-form
   finite-p corrections that are large at small p. Open: A₄ and beyond
   (the K = 12 word frontier limits the moments available).
3. **Test the 2, 3, 5 coincidence.** At T = ∞ each chord is already the
   ladder of a dimension-Δ scale-invariant operator (item 1); at small λ
   and low temperature it should give conformal τ^(−2Δ) correlators.
   Check whether O₃ inside A₃ behaves as the large-p h = 3 bilinear, and
   whether the O₅/O₃ weight in A₅ (0.53% / 99.5%) matches the ratio of
   large-p OPE coefficients (Gross–Rosenhaus). Either outcome turns an
   observation into a result.
4. **The Dirac full channel and its CP resolution** (below). This is what
   the paper's own CP assignment waits on.
5. **A proof of the decomposition for general n**, the λ^(2⌊n/2⌋) norm law,
   and which dimensions survive λ → 0. At every level computed, the first
   sub-leading coefficient is (q + … + q^k)/(1 + q^(k+1)), which suggests
   extending the n = 2 involution bundle by bundle.
6. **A chord-side derivation of the finite-N pair laws** (below).
7. The O(1/q) residues, the NV constrained-sector step and the bulk-side
   items, as below.

### The tower extraction

Part II has what is settled: the pair theorem, the exact finite-N pair
laws and double-scaling constants, the full A₂ channel at N = ∞ with its
split into a conserved piece and a Δ = 2 chord, and the n = 3..6
decomposition. What remains is below.

**A derivation of the finite-N pair laws.** At fixed p the 1/N laws are
exact by interpolation: 55 laws through the K = 12 word frontier, deficit
K/2 + 1 at every one, the leading coefficient shared across the
n = 1 → 2 step and no other, and the k = 2 family
−2^(2n+2) p^(n+1)(p−1) Q_n(p) with Q₁ = 1, Q₂ = 4p − 5,
Q₃ = 34p² − 108p + 91. The double-scaling constants are closed forms,
c₂(λ) = 4λ(1−q)²(5+3q) and c₃(λ) = 4λ(1−q)³(35 + 56q + 35q² + 10q³),
equal to the Δψ → 0 slopes of the diagonal chord channel. None of this is
derived on the chord side. The sharpest target is the O(Δ) cancellation
between the diagonal and pair halves on every double-scaling line, which
makes the pair channel's constant the diagonal channel's slope; it is
measured to 1e-8 at K = 6, not derived. The p = 4 frontier laws behind
`pair_moment_ladder.py --deep` are compute-gated (the (5,2) law alone is
past ten hours on seven workers); they would add two checks to a rule
that already holds at 55 of 55.

**A proof of the decomposition for general n.** The n = 2 split is a
theorem. n = 3..6 are exact identities verified in computer algebra on
finite word sets. Open: a proof for general n, the reason the new
chord's norm scales as λ^(2⌊n/2⌋), and a derivation of which dimensions
survive λ → 0.

**Moments → spectrum.** The forward direction (`duality/boost_forward.py`):
boost-kernel moments diverge, S → 1/ω, the regulated invariants sweep the
whole SYK band as the regulator turns, and the coupling-weighted tower
sits within 5% of a single mass at the invariant level. The ε*
consistency table in `boost_forward.json` is the recorded spec. The
inverse direction (`duality/boost_inverse.py`) answers "how many
moments" on the solvable kernel: at least as many as there are levels
(the double-precision level count is bounded by min(K, L)), with L
moments sufficing only at aε ≥ 0.4, and each level costs 0.8–2.0 more
digits of moment precision than the last
(log₁₀(σ₁/σⱼ) ≈ a(j−1) + b(j−1)², a = 1.18 → 0.44 over aε = 0.1 → 1.6,
b ≈ 0.05–0.07). Double precision buys 9–12 of 12 levels; 8 digits buy
6–10. Resolving a δ-split pair takes a third moment known to ~0.05 δ².
Moment invariants carry regulator leverage ≈ 1.0, so with the state map
unfixed a single invariant carries no mass digits, while node positions
carry exactly zero leverage under any multiplicative regulator. The
node verdict is quantitative, conditional on the e^(−εω) regulator and
the free kernel.

What the SYK side supplies at N = ∞ is no longer open: the full Aₙ
channels are a conserved polynomial in H plus matter chords of dimension
2..n whose norms and dressings are closed forms in λ, so their spectral
densities are known — smooth continua, no lines, no nodes. At fixed p the
full width is an exact law, w₂ = 8p(p−2)/N + O(1/N²) in bandwidth units,
an N-independent mean-square frequency 8p(p−2) in units of J²
(`syk-self-averaging/full_channel_laws.py`). Which finite-N object
carries the tower, and how a dictionary maps this structure to masses,
remain open.

**The Dirac full channel.** The complex-SYK chord machinery has its
H-sector built (`duality/dirac_chord.py`, on Berkooz–Narovlansky–Raj
[2006.13983](https://arxiv.org/abs/2006.13983)). The chord sum is the
Majorana one at q = e^(−p²/Nc) = e^(−λ_paper), so the factor-2 landmine
is Majorana-only. It is verified microscopically and against ED (m₄/m₂²
deviation +0.198 → +0.004 over Nc = 8 → 12, c_symmetric).
M₁ = −(p/2)iH is a transfer-matrix insertion, so it is an exactly
conserved channel. Charged random probes carry Δ = (a+b)/p, with
fixed-charge q(𝒬) and BNR asymmetry factors.

The pair (i ≠ j) half of the Dirac question is closed. The tower's pair
words carry bare single-fermion chords, the lone-crossing census is
exactly 1 − p/Nc → q_dirac^(1/p), and the conserved-channel theorem
applies verbatim at (q_dirac, Δ = 1/p): every Dirac pair moment vanishes
at Nc = ∞, with X̂₁ = −(1−q^(1/p))T the image of M₁ = −(p/2)iH. No
published source provides H-correlated chord rules, and the pair half
did not need them.

Still open, and what the CP-resolved data waits on: the full
(flavor-summed) channel at the Dirac q, and the CP resolution of the
Dirac Mₙ. The Majorana route (sum the flavor first, so the fermion marks
disappear and Aₙ becomes an intersection-weighted Hⁿ) does not transfer
as it stands: cᵢcᵢ† = 1 − nᵢ replaces ψᵢ² = 1, so rewriting the Dirac
M₂ = −Σᵢ(∂H/∂cᵢ†)†(∂H/∂cᵢ†) as an intersection-weighted H² leaves an
occupation factor 1 − nᵢ between its two halves. That factor is the
thing to handle, and it is where the charge dependence will enter.
Comparison to `thooft-target` then goes through the BOOST.md kernel and
`dictionary.py`.

### The O(1/q) residues

The large-q singlet 4-point function has Sommerfeld–Watson poles in only
two places: m = 0, which is Hamiltonian exchange and the entire
time-ordered part, and m = ±v, the scramblon.

The would-be tower cancels identically at leading order, and its
couplings are O(1/q²). It first appears at O(1/q), organized on the
exact-in-temperature skeleton 𝓜ᵉ/𝓜ᵒ, whose finite-temperature residues
are unpublished. Computing them is the open task. And how the 1/q
expansion at fixed λ relates to the 1/N one deserves a clean answer.

The skeleton itself is now a derived object rather than a transcription
(`duality/largeq_kernel.py`): the Maldacena–Stanford finite-T kernel
eigenproblem at k = 1, discretized by parity-split sine-Galerkin,
reproduces the CMS quantization sets to 3e-12 across v = 0.05..0.95, and
locates m = 0 and m = ±v as the exact non-Dirichlet solutions whose
eigenfunctions degenerate. That is the perturbable substrate the O(1/q)
computation needs; the residues on it are still the open task.

### The tomperature map

The antipodal shift passes both necessary conditions semiclassically. But
a sign flip plus effective thermality under a contour shift is not a
dictionary.

What would upgrade it is the algebra story: why the observer's operators
live at the antipodal point. That is exactly NV's construction to
complete.

### The NV discriminating computation

The OTOC sign of Narovlansky–Verlinde's dressed operators. They engineer
KMS at the dS temperature, which is the no-go's hypothesis, so they must
break another hypothesis of the bound. This is the sharpest next
computation on the dS axis, and it needs the two-sided constrained chord
sector, which is adjacent to the existing engine.

`duality/NV-DRESSING.md` now pins the identification at the
contour level: their Eqs. 33–34 and 61–62 against the fake-period shifts,
with the composition algebra and the 2-point thermality tested. What
remains is the constrained-sector step, which is whether NV's BRST
projection renormalizes the crossed growth coefficient. Conditional on it
acting trivially, the answer is already measured: anti-scrambling,
semiclassical, transient.

### The bulk side

**JT-dS ε = 3/(8S) — settled as far as it settles**
(`duality/jt_ds_vacuum.py`). The general-D law is
ε_D ≡ ρG^(D/(D−2)) = c_D S^(−2/(D−2)) with c₄ = 3/8 exactly,
sympy-certified; the paper's relation is its D = 4 instance, and the
1/S power is dimension-specific. The direct D → 2 continuation
degenerates fourfold (no dimensionless ρG^a, G_μν ≡ 0, Λ and c_D vanish,
the ℓ-elimination dies), machine-checked. What JT-dS actually supplies
is the exact first law δS = −2πℓ²ρ = −δE/T_dS (the dilaton shifts, ℓ is
rigid; the sign of SVVV [2203.00700](https://arxiv.org/abs/2203.00700)
Eq. 9), plus the reduction-inherited family — S⁻² from a dS₃ parent,
machine-checked end to end; S⁻¹ from near-Nariai dS₄ (MTY
[1904.01911](https://arxiv.org/abs/1904.01911)) — whose exponent the
intrinsic JT data (φ₀, φ_r, ℓ) cannot choose. That parent choice is a
named dictionary input, alongside the tomperature map. The measured 2d
condensate now has a first-law image, δS = +2π(0.0407)N_c(gℓ)²: finite,
with nothing for Holographic Tuning to cancel, on both sides of the
first law.

**The finite-N coupling curve ḡ²(N).** Their Fig. 17 is a sketch with no
formula, yet the whole set-it-and-forget-it inversion rests on it. Match
a bulk observable at a sequence of finite N and extrapolate. Gated on the
tower extraction.

**The order-of-limits problem.** T_B → ∞, N → ∞, λ fixed against λ → 0,
and the boost-spectrum → invariant-mass reorganization. The
chaos-to-integrability transition at small λ
([2403.01950](https://arxiv.org/abs/2403.01950),
[2601.09801](https://arxiv.org/abs/2601.09801)) is the likely mechanism
zone. This is where the duality most plausibly *fails*, and a sharp
statement of where it holds would be a real contribution. The first
exact probe of that corner is in (`syk-self-averaging/full_channel_laws.py`):
along fixed p the identity-subtracted A₂ width in bandwidth units peaks
at λ_chord ≈ 4/3 (N = 24 at p = 4) and then dies as 8p(p−2)/N, while the
N = ∞ ratio at fixed λ (`duality/full_channel_chord.py`) peaks near
λ_chord ≈ 1 and falls to 4λ as λ → 0. The ED bench's N ≤ 20 window sat on
the rising side of that peak (ERRATA E4). Both peaks are the conserved H² share
of μ₀ changing (`duality/singlet_tower.py`): with it removed the
dynamical width is 2(1 − q²) at N = ∞, monotone in λ, and along fixed p
it falls monotonically from about 2 to the same 8p(p−2)/N. In physical
units the channel keeps a finite O(J) width, √(8p(p−2)) ≈ 2√2·p, all the
way into the corner, approached monotonically from below.

**Freeness sharpening.** The ln N onset-time question needs sizes beyond
ED.

### Out of scope

Proving DSSYK∞ = JT-dS. The authors assume it themselves.

The 1+1 → 3+1 uplift. Contrary to how it is sometimes described, *no*
uplift construction exists in the source. §§5–6 are explicit that it is a
hope.

Instance-level ED in the true double-scaling regime (p ~ √N, N ≳ 100,
dim 2⁵⁰). That is the honest wall. The exact-Wick and chord computations
are the ways around it.

---

## Context

The QCD₂ side is effectively solved beyond FLZ. Ambrosino–Komatsu recast
't Hooft's equation as a TQ-Baxter integrable system
([2312.15598](https://arxiv.org/abs/2312.15598),
[2406.11078](https://arxiv.org/abs/2406.11078)), and
Litvinov–Meshcheriakov ([2409.11324](https://arxiv.org/abs/2409.11324))
give independent high-precision spectra. So numerical effort should go
almost entirely to the DSSYK side.

The DSSYK toolkit is mature. Chord techniques
([1811.02584](https://arxiv.org/abs/1811.02584); review Berkooz–Mamroud
[2407.09396](https://arxiv.org/abs/2407.09396); path integral
[2403.05980](https://arxiv.org/abs/2403.05980)). And the
chaos-to-integrability transition at small λ is directly reusable for
the λ → 0 regime the duality needs.

The duality sits in contested territory. Narovlansky–Verlinde
([2310.16994](https://arxiv.org/abs/2310.16994)) and the sine-dilaton
program of Blommaert–Mertens–Papalini
([2404.03535](https://arxiv.org/abs/2404.03535)) propose *different*
bulk duals. Okuyama ([2505.08116](https://arxiv.org/abs/2505.08116))
reaches dS JT gravity from DSSYK by a further route. And Rahman–Susskind
([2407.12988](https://arxiv.org/abs/2407.12988)) argues the
identifications conflict. A clean meson-tower match, or mismatch, would
be the first quantitative discriminator among all three programs.
`duality/DISCRIMINATOR.md` lays out that axis.

**Convention landmine.** This paper uses λ = p²/N with no q. The chord
literature uses λ = 2p²/N, q = e^(−λ). Any reuse of chord results must
pass through a tested conversion.

**Priority.** As of 2026-10-14 (INSPIRE, arXiv) the target paper is
still v1 and has exactly one citation, Cui–Kolchmeyer
[2607.13665](https://arxiv.org/abs/2607.13665). That paper argues that
dS observer OTOCs *anti-scramble*, which poses a structural challenge to
observer-centric static-patch holography. Its
same-week companion Harlow–Zhao
[2607.14215](https://arxiv.org/abs/2607.14215) proposes the
Euclidean-fold mechanism, and poses the SYK test that
`duality/ANTISCRAMBLING.md` answers.

The authors' promised meson-spectrum computation has not appeared. The
first-mover race on the headline test is still open.
