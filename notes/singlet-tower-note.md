# The singlet tower of double-scaled SYK at N = ∞: a conserved piece plus matter chords

Henry Donahue — working technical note, 4 October 2026; revised 14 October 2026.

**Purpose.** State, with its proof and its checks, what the Majorana
singlet channels Aₙ = Σᵢ ψᵢ (ad_H)ⁿ ψᵢ are at strict double scaling. This
is the Majorana form of the operators Mₙ of Miyashita–Sekino–Susskind
[2607.05678v1, §4.3](https://arxiv.org/abs/2607.05678v1). The n = 2
theorem and its proof also appear in the LaTeX manuscript
[conserved-channel-note](conserved-channel-note/main.tex) §4. This note
does not extract a mass spectrum, and it does not cover the Dirac
tower. Novelty
has not been established by an exhaustive literature review.

Code: `duality/singlet_tower.py`. Tests: `duality/test_singlet_tower.py`.
Data: `duality/singlet_tower.json`.

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
Ω₃, Ω₄: kernels computed; new chords Ô₅, Ô₆

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
rational arithmetic at three generic points for k = 3 and 4. λ and q are
treated as independent, because the identities are polynomial in both. A
proof for general k is open.

**What the levels contain.** Every Aₙ computed (n ≤ 6) is a polynomial in
H plus matter chords of dimension 2 through n, each dressed by a
polynomial of degree n − m in H acting from the left and right. One new
chord, of dimension n, enters at level n. Its weight is small:

| channel | O₂ | O₃ | O₄ | O₅ | O₆ |
|---|---|---|---|---|---|
| A₂ | 100% | | | | |
| A₃ | 85.1% | 14.9% | | | |
| A₄ | 83.0% | 15.3% | 1.6% | | |
| A₅ | 76.9% | 20.3% | 2.5% | 0.18% | |
| A₆ | 72.5% | 22.7% | 4.4% | 0.39% | 0.02% |

(shares of the dynamical weight at λ = 2). Over λ from 0.05 to 4 the
weight of A₄, A₅, A₆ outside O₂ and O₃ stays between 0.9% and 5.4%.

The norms are graded: Nₘ scales as λ^(2⌊m/2⌋), with N₂ → λ², N₃ → 2λ²,
N₄ → (52/5)λ⁴ from the closed forms and N₅ ≈ 12.5λ⁴, N₆ ∝ λ⁶ measured.

Shares do not follow norms, because the dressings carry powers of the
frequency. The O₂ dressing in A₃ is (3/2)[H, A₂], a total derivative,
plus an anticommutator admixture −(1−q)²/(3(1+q²)) of it, and a
derivative costs a factor of the width, O(λ). Measured in exact
arithmetic down to λ = 1/1280:

| λ → 0 | O₂ | O₃ | O₄ | O₅ | O₆ |
|---|---|---|---|---|---|
| A₃ | ∝ λ | → 100% | | | |
| A₄ | ∝ λ | → 100% | ∝ λ | | |
| A₅ | ∝ λ | → 99.5% | ∝ λ | → 0.53% | |
| A₆ | ∝ λ | → 99.2% | ∝ λ | → 0.80% | → 0 |

A₂ is O₂ at every λ. So the dimensions that survive the semiclassical
limit are 2, 3 and 5. That is the pattern of the bilinear spectrum of
large-q SYK, h = 2, 3, 5, 7, …
([1604.07818](https://arxiv.org/abs/1604.07818)). It is recorded here as
an observation on four levels. It is not derived, and the map between a
matter-chord dimension at infinite temperature and a conformal dimension
has not been checked beyond this coincidence.

Closed forms for the moments follow, for example

```text
μ₂(A₃) = 8λ(1−q)[25 + 5q − 38q² − 16q³ + 13q⁴ + 11q⁵
                  + λ(46q + 14q² − 22q³ − 26q⁴)].
```

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
  of §4, for A₃ through A₆.

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
chords of dimension 3, 4, … whose norms are closed forms in λ and scale
as λ^(2⌊m/2⌋), with the even ones dropping out of the higher levels as
λ → 0. In the paper's dictionary ḡ² = p²/N, so ḡ² = λ/2 in this note's
chord convention.

## 8. What this implies for labelling the tower

The paper pairs the n-th operator Mₙ with the n-th meson. The small-λ
shares of §4 bear on that labelling, in the Majorana channels and at
fixed λ.

As λ → 0, A₃, A₄, A₅ and A₆ each put at least 99% of their dynamical
weight in one sector, the dimension-3 chord O₃, with dressings that
differ from level to level. So in that limit the raw Aₙ with n ≥ 3 do
not each open a new propagating sector. A one-to-one map from raw
operators to new states cannot hold there.

The decomposition also says what can be labelled. It is triangular: Oₙ
first appears in Aₙ, and its dressing there has degree n − n = 0. So
removing from Aₙ its overlap with the commutant of H and with every
dressed lower chord f(H) Oₘ g(H), m < n, leaves a multiple of Oₙ at
N = ∞: a matter chord of dimension n whose spectral density is the
standard Δ = n continuum. The natural candidate for "the n-th member of the
tower" is therefore the dimension-n chord, not the raw Aₙ. Both start
at 2: the paper's dynamical tower begins at n = 2, and so do the matter
chords.

Whether a dimension-n chord carries a meson's quantum numbers, and at
what mass, is the dictionary question. It is not answered here. These
statements are about the Majorana Aₙ; the Dirac Mₙ are not covered.

## Reproduce

```bash
cd duality
python singlet_tower.py                  # headline numbers, seconds
python -m pytest -q test_singlet_tower.py
python singlet_tower.py --run            # artifact and finite-N gates
```
