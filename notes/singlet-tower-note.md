# The singlet tower of double-scaled SYK at N = ∞: a conserved piece plus matter chords

Henry Donahue — working technical note, 4 October 2026; revised 14 and 15 October and 10 November 2026.

**Purpose.** State, with its proof and its checks, what the Majorana
singlet channels Aₙ = Σᵢ ψᵢ (ad_H)ⁿ ψᵢ are at strict double scaling. This
is the Majorana form of the operators Mₙ of Miyashita–Sekino–Susskind
[2607.05678v1, §4.3](https://arxiv.org/abs/2607.05678v1). The n = 2
theorem and its proof also appear in the LaTeX manuscript
[conserved-channel-note](conserved-channel-note/main.tex) §4. This note
does not extract a meson mass spectrum, and it does not cover the Dirac
tower. Novelty
has not been established by an exhaustive literature review.

Code: `duality/singlet_tower.py` (and `duality/conformal_tower.py` for
§4a). Tests: `duality/test_singlet_tower.py`. Data:
`duality/singlet_tower.json`, `duality/singlet_tower_deep.json`.

## 1. Conventions

{ψᵢ, ψⱼ} = 2δᵢⱼ, p even, H = i^(p/2) Σ_c J_c ψ_c over p-subsets c, with
independent Gaussian J_c. The trace is normalized. λ = λ_chord = 2p²/N,
q = e^(−λ), and energies are in chord units, E[Tr H²]/dim = 1. Hᵢ is the
part of H whose monomials contain flavor i, so ψᵢHψᵢ = H − 2Hᵢ and
ΣᵢHᵢ = pH.

A *matter chord of dimension Δ* is the standard double-scaled object: an
operator whose two-point words are chord sums in which every H-chord
crossing the matter chord weighs q^Δ. Its moments are the closed forms
μ₂ = 2(1 − q^Δ), … of the chord literature
([1811.02584](https://arxiv.org/abs/1811.02584),
[2407.09396](https://arxiv.org/abs/2407.09396)).

## 2. An exact operator identity

For a set T of slots let W_T = Σᵢ Πₜ (Hᵢ if t ∈ T, else H), a word of
length n. For every n ≥ 2, every N and every even p, per realization,

```text
Aₙ = Σ_{|T|≥2} (−2)^|T| (−1)^(n − max T) C(n−1, max T − 1) W_T.
```

n = 2 is A₂ = 4ΣᵢHᵢ². n = 3 is
A₃ = 4Σᵢ(HᵢHHᵢ + HHᵢ² − 2Hᵢ²H) − 8ΣᵢHᵢ³.

*Derivation.* S(X) = ΣᵢψᵢXψᵢ multiplies a Clifford monomial of grade g
by N − 2g when p is even. A product of m Hamiltonian monomials has grade
|c₁ △ … △ c_m| = Σ_T (−2)^(|T|−1) |∩_T c_t|. Expanding the nested
commutator gives Aₙ = Σⱼ(−1)ʲC(n,j) S(Hⁿ⁻ʲ)Hʲ, the |T| ≤ 1 terms cancel
in the binomial sum, and |∩_T c_t| = Σᵢ Π_{t∈T}[i ∈ c_t] turns H into Hᵢ
on the slots of T. Checked against dense matrices for n ≤ 5
(`syk-self-averaging/test_conserved_points.py`).

A corollary used below: A₃ − A₃† = 3[H, A₂], exactly.

## 3. The n = 2 theorem

**Statement.** At strict double scaling, for every λ > 0,

```text
A₂ − ⟨A₂⟩ = [2λ(1−q)/(1+q)] (H² − 1)  +  ν O₂,
ν² = 8λ(1 − q² + 2λq)/(1+q),
```

where O₂ is a unit-normalized matter chord of dimension 2, uncorrelated
with every function of H. Equivalently, for all A, B ≥ 0,

```text
⟨H^A A₂ H^B A₂⟩_c = [2λ(1−q)/(1+q)]² ⟨H^(A+B)(H²−1)²⟩ + ν² G₂(A, B),
```

with G₂ the two-point word of a Δ = 2 matter chord.

**Consequences.**

```text
μ₀ = 4λ²(1−q)²/(1+q)  [conserved]  +  ν²  [dynamical]
w₂^dyn = μ₂/ν² = 2(1 − q²),          every λ.
```

At λ = 2: μ₀ = 32, of which 10.5364 (32.9%) is conserved; the dynamical
width is 1.963369. The identity-subtracted ratio μ₂/μ₀ = 1.316906 quoted
in earlier notes is that width diluted by the conserved weight. The
dynamical kurtosis and every higher moment are the Δ = 2 matter-chord
values. The complete spectral density of the dynamical part is the
Δ = 2 matter-chord density, a smooth monotone continuum.

**Proof.** A₂ = 4Ω₀ with Ω₀ = ΣᵢHᵢ² = Σ_{c,c′}|c ∩ c′| J_cJ_c′ψ_cψ_c′:
the word H² with its two monomials weighted by the number X of indices
they share. At double scaling the intersections of distinct chords are
independent Poisson(μ), μ = λ/2, and a crossing of two chords carries
(−1)^X. So the weighted pair carries E[X(−1)^(sX)] = μ or −μq according
as its chords cross (s = 1) or not, and E[X²(−1)^(sX)] = μ + μ² or
q(μ² − μ) when both insertions weight the same pair. A contraction of
the two slots of one insertion with each other is the c-number ⟨A₂⟩.

In ⟨H^A Ω₀ H^B Ω₀⟩ sort the chord diagrams into two classes.

*Class I: both chords of one insertion end on the other.* The two chords
form a bundle. Every H-chord separating the insertions crosses both
strands and weighs q². Nested plus crossed, the bundle weighs
(μ + μ²) + q(μ² − μ). The class sums to [μ(1−q) + μ²(1+q)] G₂(A, B).

*Class II: they do not.* At each insertion, exchanging which of its two
adjacent slots carries which chord is an involution on this class. No
node lies between adjacent slots, so the exchange flips the mutual
crossing of the two chords and changes no other crossing. It maps the
weight μ to −μq. On the plain word H² the same exchange maps 1 to q.
The two involutions (one per insertion) commute. Hence class II equals
[μ(1−q)/(1+q)]² times class II of ⟨H^A(H²−1)H^B(H²−1)⟩.

Class I of that plain word is (1+q)G₂. Adding and subtracting it,

```text
⟨H^AΩ₀H^BΩ₀⟩ = [μ(1−q)/(1+q)]² ⟨H^(A+B)(H²−1)²⟩
              + [μ(1−q) + μ²(1+q) − μ²(1−q)²/(1+q)] G₂(A, B),
```

and the bracket is λ(1 − q² + 2λq)/(2(1+q)) = ν²/16. ∎

The only inputs are the identity of §2 and the Poisson limit. The proof
uses no property of the q-Hermite polynomials; H² − 1 appears because it
is what a two-slot insertion without self-contraction is.

## 4. The levels above

In the double-scaled limit the triple intersections in §2 are O(1/p), so

```text
Aₙ ≃ Σ_{a<b} κ_ab H^(a−1) Ω_(b−a−1) H^(n−b),   Ωₖ = Σᵢ Hᵢ Hᵏ Hᵢ,
κ_ab = 4(−1)^(n−b) C(n−1, b−1),
```

up to exactly conserved operators. One subtlety is not optional. When
the first and last slot of Ωₖ contract with each other, X = p, and the
enclosing chord's crossing average depends on the grade of what it
encloses: p(1 − 2p/N)^grade = p qᵏ[1 + (2λ/p)Σ_pairs X + …]. The O(p)
term multiplies a polynomial in H. The O(1) remainder is 2λ times an
X-weight on each enclosed pair, and it enters at k = 2.

With Ôₘ the un-normalized matter chord of dimension m (norm Nₘ), and
modulo polynomials in H:

```text
Ω₀ = Ô₂
Ω₁ = [q/(1+q²)] {H, Ô₂} + Ô₃
Ω₂ = a(H²Ô₂ + Ô₂H²) + b HÔ₂H + c Ô₂ + [q/(1−q+q²)] {H, Ô₃} + Ô₄
Ω₃ .. Ω₆: kernels computed; new chords Ô₅ .. Ô₈

a = q²(1+q+q²)/((1+q²)[5]),  b = q²(1+q)²/((1+q²)[5]),
c = (1 + q − q² − q³ + 2λq²[5])/[5],   [5] = 1+q+q²+q³+q⁴

N₂ = λ(1 − q² + 2λq) / (2(1+q))
N₃ = λ(1+q)[1 − q⁶ + 2λq²(2 − q + 2q²)] / (2(1+q²)(1+q+q²))
N₄ = λ(1−q)²(1+q)²[(1−q⁶)(1+q+2q²+q³+q⁴) + 2λq³(2+q+2q²+q³+2q⁴)]
     / (2(1+q²)(1−q+q²)[5])
```

The first sub-leading coefficient of every block computed is
(q + … + qᵏ)/(1 + q^(k+1)): a plain chord end inside a bundle of k + 1
strands crosses only the strands on its way out.

**Status.** These are exact identities verified in computer algebra, not
proofs. Each sector's kernel is solved linearly from the cross-correlator
with a lower block; every remaining two-point word is then a check. The
residual is identically zero: symbolically in (λ, q) for k ≤ 2 (words
through four H's; through six at one rational point), and in exact
rational arithmetic at three generic points for k = 3..6 (k = 5, 6, the
levels n = 7, 8, are in `duality/singlet_tower_deep.json`). λ and q are
treated as independent, because the identities are polynomial in both.
The independent checks per point are 37 at n ≤ 6, 55 at n = 7 and 78
at n = 8; the stored word counts are about three times larger, because
a word and its mirror give the same equation and some words are already
used in the solves (ERRATA E26). What a zero residual tests is the
structure, a polynomial in H plus dressed matter chords, not the values:
it stays exactly zero under several wrong chord rules (the remainder
off or rescaled, other pair weights), because the strands running
between the two insertions always form a matter chord. The values (the
norms, the dressings and so the shares) rest on the derivation of the
rule, the n = 2 theorem, and the finite-N gates at n = 3 and 4, which do
discriminate (⟨ad²A₄A₂†⟩ at λ = 1: 155.03 with the remainder against
150.65 without). The subleading
coefficient q[k]_q/(1 + q^(k+1)) holds exactly through k = 6. A proof for
general k is open.

**What the levels contain.** Every Aₙ computed (n ≤ 8) is a polynomial in
H plus matter chords of dimension 2 through n, each dressed by a
polynomial of degree n − m in H acting from the left and right. One new
chord, of dimension n, enters at level n. Its weight is small:

| channel | O₂ | O₃ | O₄ | O₅ | O₆ | O₇ | O₈ |
|---|---|---|---|---|---|---|---|
| A₂ | 100% | | | | | | |
| A₃ | 85.1% | 14.9% | | | | | |
| A₄ | 83.0% | 15.3% | 1.6% | | | | |
| A₅ | 76.9% | 20.3% | 2.5% | 0.18% | | | |
| A₆ | 72.5% | 22.7% | 4.4% | 0.39% | 0.02% | | |
| A₇ | 68.5% | 24.6% | 6.0% | 0.84% | 0.05% | 0.002% | |
| A₈ | 65.1% | 25.8% | 7.5% | 1.4% | 0.14% | 0.006% | 0.0002% |

(shares of the dynamical weight at λ = 2). Over λ from 0.05 to 4 the
weight of A₄, A₅, A₆ outside O₂ and O₃ stays between 0.9% and 5.4%
(A₇ and A₈ were computed at λ = 2 and small λ only; there it is 6.9% and
9.1%).

The norms are graded: Nₘ scales as λ^(2⌊m/2⌋), with N₂ → λ², N₃ → 2λ²,
N₄ → (52/5)λ⁴ from the closed forms, N₅ → (88/7)λ⁴, N₆ → (4120/7)λ⁶ and
N₇ → 720λ⁶ as exact limits (§4a), and N₈ ∝ λ⁸ measured (N₇ and N₈ were
the law's predictions before they were computed).

Shares do not follow norms, because the dressings carry powers of the
frequency. The O₂ dressing in A₃ is (3/2)[H, A₂], a total derivative,
plus an anticommutator admixture −(1−q)²/(3(1+q²)) of it, and a
derivative costs a factor of the width, O(λ). Measured in exact
arithmetic down to λ = 1/1280:

| λ = 1/1280 | O₂ | O₃ | O₄ | O₅ | O₆ | O₇ | O₈ |
|---|---|---|---|---|---|---|---|
| A₃ | ∝ λ | 99.65% | | | | | |
| A₄ | ∝ λ | 99.61% | ∝ λ | | | | |
| A₅ | ∝ λ | 99.19% | ∝ λ | 0.53% | | | |
| A₆ | ∝ λ | 98.96% | ∝ λ | 0.80% | ∝ λ | | |
| A₇ | ∝ λ | 98.40% | ∝ λ | 1.38% | ∝ λ | 0.006% | |
| A₈ | ∝ λ | 97.80% | ∝ λ | 1.99% | ∝ λ | 0.014% | ∝ λ |

(∝ λ marks the fading even chords, which carry the rest, e.g. 0.26% O₂
and 0.02% O₄ in A₅. A₃ and A₄ tend to O₃ alone; the exact limits of
A₅..A₇ are in §4a, e.g. A₅ → 2047/2058 O₃ + 11/2058 O₅.) A₂ is O₂ at every λ. So the
dimensions that survive the semiclassical limit are 2, 3, 5 and 7. That
is the pattern of the bilinear spectrum of large-p SYK (large q in the
usual notation), h = 2, 3, 5, 7, …
([1604.07818](https://arxiv.org/abs/1604.07818)). Levels 7 and 8 were
computed to test it after it was noticed on levels 3–6: the pattern
predicted that O₇ survives and O₈ fades, and both held. §4a reduces the
pattern to the norm law and shows that, as λ → 0 and with the chord
rule's grade remainder off, the chords are conformal primaries of weight
m in the operator expansion of the bilocal; the full rule adds one
Liouville rung.

Closed forms for the moments follow, for example

```text
μ₂(A₃) = 8λ(1−q)[25 + 5q − 38q² − 16q³ + 13q⁴ + 11q⁵
                  + λ(46q + 14q² − 22q³ − 26q⁴)].
```

## 4a. The semiclassical limit: a conformal OPE plus one Liouville rung

Code: `duality/conformal_tower.py`. Tests: `duality/test_conformal_tower.py`.
Data: `duality/conformal_tower.json`.

Write σ = (E + E′)/2 and ω = E − E′ for the left and right energies of a
dressing F_{n,m}(E, E′), so that ad_H acts as ω. As λ → 0 a chord of
dimension m is the ladder sech^(2m)(𝒥t), with 𝒥 = √λ in chord units
(`duality/HORIZON.md`, whose ladder spacing is 2𝒥). Write a ≡ √λ in this
section, so ω ~ a while σ ~ 1. Scale ω = a w. Then, at every level taken to the
limit (n ≤ 7),

```text
F_{n,m}(E, E′) = a^(n−m) P_{n,m}(w) (1 + O(a)),
```

with P_{n,m} a polynomial with rational coefficients and no σ dependence.
The norms have rational limits too: Nₘ/λ^(2⌊m/2⌋) → νₘ = 1, 2, 52/5,
88/7, 4120/7, 720 for m = 2..7. "Exact" means exact rational chord sums
at λ = 1/(3·2^j), j = 10..16, polynomial extrapolation to λ = 0 and
rational reconstruction (denominators ≤ 1000; the largest met is 35).
Each scaled coefficient is a rational function of q = e^(−λ) times an
integer power of λ, so it is analytic in λ and the extrapolation has no
logarithms or half-integer powers to miss. The residuals are below
1e-23, no limit moves by more than 2e-14 when an end point is dropped
(a hidden divergence would move it at O(1)), and the terms that must
vanish do so to below 1e-26. It is not a symbolic proof.

**Without the remainder, Aₙ is an SL(2) operator expansion.** The chord
rule of §4 has one ingredient that is not pairwise: the remainder, an
O(1) term that appears when the first and last slots of Ωₖ contract with
each other. Switch it off. Then every P_{n,m} through n = 7 (28
coefficients below the diagonal) is exactly the prediction of a
one-dimensional conformal OPE of the bilocal,

```text
Σᵢ ψᵢ(0)ψᵢ(t) = Σₘ Cₘ (δu)^m ₁F₁(m; 2m; δu ∂_u) Oₘ(u),
u = tanh(at),  δu = u(t) − u(0),  Oₘ^(t) = (du/dt)^m Oₘ^(u),
```

with each Oₘ a primary of weight h = m and Δ_ψ → 0. The leading
coefficient is the plain block coefficient 4C(n,m)(m)ₙ₋ₘ/(2m)ₙ₋ₘ (for
example 6, 36/5, 80/7, 60/7 for (n,m) = (3,2), (4,2), (5,3), (6,2)). Each
lower coefficient is the correction for the frame:
A₄ ⊃ (36/5 w² + 16/5) O₂ and A₆ ⊃ (60/7 w⁴ + 96/7 w² + 64/7) O₂. The
frame tanh(at) is the one in which the chord two-point function takes
conformal form: sech^(2m)(a(t₁ − t₂)) = (u₁′u₂′)^m (1 − u₁u₂)^(−2m), the
two-point function with one point inverted (u₂ → 1/u₂, orientation
reversing, so up to a sign (−1)^m). e^(2at) is
Möbius-equivalent and gives the same numbers; tan(at) multiplies the
w^(n−m−2j) coefficient by (−1)^j (A₆ ⊃ (60/7 w⁴ − 96/7 w² + 64/7) O₂) and
fails. So, with the remainder off and through n = 7, the dressings are
exactly those of a weight-m primary in the bilocal's SL(2) OPE. Some
SL(2) structure at λ → 0 is expected: Lin–Stanford
([2307.15725](https://arxiv.org/abs/2307.15725)) find that a subalgebra
of the chord algebra simplifies to sl₂ in the semiclassical limit, and
the one-particle chord space carries the positive discrete series of
U_q(su(1,1)) (van der Heijden–Verlinde–Xu
[2511.08743](https://arxiv.org/abs/2511.08743)), the natural route to a
general-n proof. What is computed here is the frame, the OPE data and
the remainder's rung. The semiclassical T = ∞ four-point function (Xu
[2411.04251](https://arxiv.org/abs/2411.04251)) is a possible route to
deriving the coefficients. The
full rule adds a non-conformal rung, so the full channels are not
exactly conformal.

**The remainder is one Liouville rung.** With the remainder on, only the
lower coefficients change: A₄ ⊃ (36/5 w² + 56/5) O₂,
A₅ ⊃ (80/7 w² + 136/7) O₃. At leading order the change is exactly one
rung of the linearized large-p Liouville equation
∂₁∂₂δg = −2𝒥²e^g δg (e^g = sech², 𝒥 = a), acting once on the OPE
bilocal (the first Born term, not the resummed ladder):

```text
B(0, t) = B_OPE(0, t) − 2a² ∫_{0<t₁<t₂<t} sech²(a(t₂ − t₁)) B_OPE(t₁, t₂),
```

with B(t₁, t₂) = Σⱼψⱼ(t₁)ψⱼ(t₂). Nothing is fitted: the coupling −2 and
the shape sech² are Liouville's. This reproduces all 13 below-diagonal
coefficients the remainder changes through n = 7. Ten see only the
rung's value at τ = 0, and each fixes the coupling at −2 on its own; the
other three, (n, m) = (6, 2), (7, 2), (7, 3), see its τ² term, and each
needs exactly sech²'s −1. That the rung acts once fits the chord rule,
in which each Ωₖ insertion carries the remainder at most once.
Iterating the rung agrees through n − m = 3 and differs first, at
n − m = 4, in exactly those three coefficients. (An earlier version
iterated it and read the difference as sector-dependent τ² weights,
−5/6 and −9/10; ERRATA E24.) The next level is a prediction, recorded in
`conformal_tower.json`: for example A₈ ⊃ (28/3 w⁶ + 9892/105 w⁴ +
93488/105 w² + 76352/35) O₂, where only the w⁰ term sees sech²'s τ⁴
coefficient. The remainder is O(1) at λ → 0: it lowers O₅'s share of A₅
from 11/1911 to 11/2058, by about 7%.

**Exact semiclassical shares; why odd dimensions survive.** The weight
of the dressed Oₘ in Aₙ is Nₘ⟨|P_{n,m}|²⟩ ∝ λ^(n−m+2⌊m/2⌋). That is
λ^(n−1) for every odd m and λ^n for every even m. So the survival of
2, 3, 5, 7 follows from the norm law Nₘ ∝ λ^(2⌊m/2⌋) and the scaling
F_{n,m} ~ a^(n−m) with a nonzero limit (A₂ is O₂ alone). Both are
established here through n = 7 and neither is proved; given them, the
survival is a consequence rather than a separate observation. The
odd-only pattern is itself what a generalized free Majorana fermion
gives (its bilinear primaries have h = 2Δψ + 2n + 1, odd as Δψ → 0,
with h = 2 the energy); what the chords add is the norm law and the
coefficients below. The
surviving shares are exact rationals:

| λ → 0 | O₃ | O₅ | O₇ |
|---|---|---|---|
| A₅ | 2047/2058 (99.4655%) | 11/2058 (0.5345%) | |
| A₆ | 20443/20608 (99.1993%) | 165/20608 (0.8007%) | |
| A₇ | 540648/548261 (98.6114%) | 7578/548261 (1.3822%) | 5/78323 (0.0064%) |

The λ = 1/1280 values in the table of §4 approach these.

**The T = ∞ OPE coefficients, against the IR ones.** In the
normalization of the conformal block expansion, with νₘ the norm limits
above, Cₘ² = 16νₘ/m!² × λ^(2⌊m/2⌋−m). The even ones are O(1): C₂² = 4, C₄² = 13/45,
C₆² = 103/5670. The odd ones are O(1/λ): λC₃² = 8/9, λC₅² = 22/1575,
λC₇² = 1/2205. Compare Gross–Rosenhaus
([1702.08016](https://arxiv.org/abs/1702.08016) (2.15)): the large-p
couplings c_h² of the IR conformal SYK tower h = 3, 5, 7 (their q is
our p). Exactly, for all three,

```text
Cₕ² = N (1 − k_c(h))² cₕ²,     k_c(h) = 2/(h(h−1)),   N = 2p²/λ,
```

where k_c is the p = ∞ conformal ladder-kernel eigenvalue
(Maldacena–Stanford [1604.07818](https://arxiv.org/abs/1604.07818)
(3.76)), taken as the regular function 2/(h(h−1)). (At large p the
kernel's poles at h = 1 + 2/p + 2n land on the odd integers, so the
pointwise limit there is singular.) The form (1 − k_c)² was read off
from the ratios at h = 5 and 7, and it is one of several equivalent
factorizations: it is also GR's formula with their ε_h² replaced by
(1 + 2/(h(h−1)))². The form that needs no GR is

```text
Cₕ² = N c²_GFF(h) (1 + k_c(h))/(1 − k_c(h)),
c²_GFF(h) = 4Δψ² k_c(h) Γ(h)²/Γ(2h − 1),   Δψ = 1/p,
```

with c²_GFF the Δψ → 0 coefficient of the h-th bilinear primary of a
generalized free fermion (the block expansion of
F₀ = −χ^(2Δ) + (χ/(1−χ))^(2Δ), checked in the tests): the free
coefficient times the ladder factor (1 + k_c)/(1 − k_c) =
1 + 2k_c + 2k_c² + …. It is the same relation rewritten, not further
evidence. Since the norms, and so the Cₕ, are identical with the
remainder on and off, that ladder factor is not the rung above. The
absolute normalization at h = 3, with exactly the factor N, is the one
out-of-sample check, given ψ² = 1 and GR's four-point normalization of
c_h². That is a single number, and the form was chosen after seeing
the h = 5 and 7 data, so the relation is suggestive, not established,
and it is not derived. It ties the infinite-temperature chord tower to
the IR bilinear tower operator by operator, with no continuous free
parameter. It predicts that the
dimension-9 norm tends to (21482496/143)λ⁸ ≈ 1.5·10⁵λ⁸, beyond the
levels computed. The even chords have no counterpart in the large-p
tower (whose only even member is h = 2), and their T = ∞ couplings are
a power of λ smaller.

## 5. Checks that do not use the proof

1. *Exact finite N, n = 2.* The exact rational full-channel moments on
   two double-scaling lines, with the H² component removed in exact
   arithmetic (the overlap is a hypergeometric sum), extrapolate in 1/p
   to the theorem: at λ = 2, w₂^dyn → 1.96336 against 1.963369; at
   λ = 1, 1.7293 against 1.729329.
2. *Dense ED, n = 2.* Per realization, with the entire commutant of H
   projected out. The stationary weight equals the H² projection within
   errors. At p = 6, N = 16, 18, 20 the dynamical width is 2.03, 2.03,
   2.01 against 2(1 − q²) = 2.00; the identity-subtracted width there is
   0.29–0.41.
3. *Exact finite N, higher levels.* μ₄(A₂), μ₂(A₃), μ₄(A₃), μ₂(A₄) and
   the cross moment ⟨ad²(A₄)A₂†⟩ approach the closed forms on
   double-scaling lines. The cross moment isolates the remainder rule of
   §4: 155.03 with it and 150.65 without at λ = 1; the exact series to
   p = 20 extrapolates to 155.027.
4. *Spectral density.* The decomposition pushed through the closed-form
   matter vertex integrates back to μ₀^dyn and μ₂.

The 1/p limits are least-squares extrapolations; their spreads are
sensitivity estimates, not confidence intervals.

## 6. What changes in the earlier notes

(The full list of corrections is in [ERRATA.md](../ERRATA.md); E2–E4
are the ones this note settles.)

- "μ₂/μ₀ = 1.316906 at λ = 2" stands as a ratio of identity-subtracted
  moments. As a width it is superseded by 2(1 − q²) = 1.963369.
- The maximum of that ratio near λ ≈ 1 at N = ∞, and its maximum at
  N = 24 along fixed p = 4, are the conserved share of μ₀ changing. With
  the H² component removed the width is monotone in both cases. The
  fixed-p asymptote 8p(p−2)/N and the physical mean-square frequency
  8p(p−2)J² are unaffected.
- "A complete projected large-N spectral density is still open" is
  closed for the Majorana A₂ channel and, with the verification status
  of §4, for A₃ through A₈.

## 7. What is not claimed

No meson masses. Each sector is a continuum in boost frequency, and the
relation of boost-frequency data to invariant masses needs the
state/time/operator dictionary (see `duality/BOOST.md`). No statement
about the Dirac Mₙ: there cᵢcᵢ† = 1 − nᵢ replaces ψᵢ² = 1, and the flavor
sum leaves an occupation factor between the two halves of M₂ instead of
removing the fermion marks. No statement at fixed p beyond
the exact laws already in `syk-self-averaging/full_channel_laws.py`; the
N = ∞ results here hold at fixed λ.

For a dictionary the content is a specification: a conserved H²
admixture of stated size, a dimension-2 continuum of stated norm, and
chords of dimension 3, 4, … whose norms are closed forms for m ≤ 4 and
exact rationals at any point for m ≤ 8, and scale
as λ^(2⌊m/2⌋), with the even ones dropping out of the higher levels as
λ → 0. In the paper's dictionary ḡ² = p²/N, so ḡ² = λ/2 in this note's
Majorana chord convention. The paper's complex model has its own chord
convention (q = e^(−4p²/N) with p creators per term, 2006.13983), so the
factor is convention-dependent.

## 8. What this implies for labelling the tower

The paper pairs the n-th operator Mₙ with the n-th meson. The small-λ
shares of §4 bear on that labelling, in the Majorana channels and at
fixed λ.

As λ → 0, A₃ through A₈ each put at least 97% of their dynamical
weight in one sector, the dimension-3 chord O₃, with dressings that
differ from level to level. So in that limit the raw Aₙ with n ≥ 3 do
not each open a new propagating sector. A one-to-one map from raw
operators to new states cannot hold there.

The decomposition also says what can be labelled. It is triangular: Oₙ
first appears in Aₙ, with constant dressing and coefficient exactly 4. So
removing from Aₙ its overlap with the commutant of H and with every
dressed lower chord f(H) Oₘ g(H), m < n, leaves Õₙ = 4Ôₙ at N = ∞, of
weight 16 Nₙ: a matter chord of dimension n whose spectral density is the
standard Δ = n continuum (`duality/horizon_tower.py`, checked
n = 2..5). The natural candidate for "the n-th member of the tower" is
therefore the dimension-n chord, not the raw Aₙ. Both start at 2: the
paper's dynamical tower begins at n = 2, and so do the matter chords.

Whether a dimension-n chord carries a meson's quantum numbers, and at
what mass, is the dictionary question. It is not answered here. What a
chord looks like in a Rindler frame is: as λ → 0 its density is the
quasinormal ladder |Γ(n + iω/2𝒥)|² of a scale-invariant operator of
dimension n, whose Rindler mass spectrum is gapless (density ∝ M^(2n−1)),
with no gap and no discrete masses (`duality/HORIZON.md`). These
statements are about the Majorana Aₙ; the Dirac Mₙ are not covered.

## Reproduce

```bash
cd duality
python singlet_tower.py                  # headline numbers, seconds
python -m pytest -q test_singlet_tower.py
python singlet_tower.py --run            # artifact and finite-N gates
python singlet_tower.py --deep 5 6       # levels n = 7, 8 (~25 min on 6 cores)
python conformal_tower.py                # §4a, the λ → 0 limit through n = 5, seconds
python conformal_tower.py --run          # its artifact, n ≤ 7 (~4 min, one core)
python -m pytest -q test_conformal_tower.py
```
