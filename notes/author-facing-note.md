# Full-channel moments and normalization checks for the DSSYK–'t Hooft proposal

Henry Donahue — working technical note, 16 September 2026; §2 and the
priorities revised 4 October 2026.

**Purpose.** Provide specific calculations relevant to the singlet-spectrum
question in Miyashita–Sekino–Susskind, [2607.05678v1, §4.3](https://arxiv.org/abs/2607.05678v1),
with the distinction between fixed-p and double-scaled limits described in
[2506.18054v3, §5](https://arxiv.org/abs/2506.18054v3).
No meson mass match or no-go theorem is claimed. Novelty has not been
established by an exhaustive literature review.

## 1. Definitions and exact operator checks

In the Majorana calculation, {ψᵢ,ψⱼ} = 2δᵢⱼ, p is even,
H = i^(p/2)Σ_c J_c ψ_c, and Aₙ = Σᵢψᵢ(ad_H)ⁿψᵢ.
The trace is normalized by Hilbert-space dimension. Let Hᵢ be the part
of H whose monomials contain index i. Conjugation by ψᵢ gives
ψᵢHψᵢ = H−2Hᵢ, and ΣᵢHᵢ = pH. Consequently

```text
A₁ = −2pH,                    A₂ = 4ΣᵢHᵢ².
```

For A₂, expanding the double commutator gives
Σᵢ(ψᵢH²ψᵢ−2ψᵢHψᵢH+H²). Substituting H−2Hᵢ and summing cancels
the mixed commutator term because Σᵢ[Hᵢ,H] = [pH,H] = 0.
These are per-realization identities, independently checked against dense
matrices in `conserved_points.py` and `test_full_channel_chord.py`.
They are the first two cases of a general identity,
Aₙ = Σ_{|T|≥2} (−2)^|T| (−1)^(n−max T) C(n−1, max T−1) W_T, with W_T the
word Hⁿ carrying Hᵢ on the slots of T and summed over i (checked for
n ≤ 5; derivation in [singlet-tower-note.md](singlet-tower-note.md) §2).

The paper's U(N) discussion instead uses complex fermions and
Mₙ = Σᵢcᵢ†(d/dt)ⁿcᵢ. For the repository's normal-ordered,
charge-conserving p-body Hamiltonian and d/dt = i[H,·], monomial counting
gives M₁ = −(p/2)iH. Thus conservation agrees with the paper, while
its literal M₁ = H needs a normalization/phase convention. Differentiating
M₁ gives M₂ = −Σᵢċᵢ†ċᵢ. The Majorana A₂ result below has not been
identified with the complete Dirac CP-resolved meson channel.

## 2. Leading pair cancellation does not eliminate the full channel

In the standard sequential matter-chord representation,
C(m) = close T_M^m open is a polynomial in the chord Hamiltonian T.
The three-term recursion in `chord_charges.py` proves this for every m;
it is not just an n = 1,2,3 numerical observation. Hence the resulting
leading two-flavor insertion has vanishing commutator moments at every
positive chord coupling and admissible matter weight.

This is not a proof that the flavor-summed Aₙ is conserved. For Xᵢ =
ψᵢ(ad_H)ⁿψᵢ, the full correlator has both N diagonal terms and N(N−1)
pair terms. Subleading per-flavor terms can survive their multiplicities.
The pair cross-correlator is not a positive spectral measure on its own.
The claim that the physical tower is invisible to all orders of double
scaling is therefore withdrawn.

A direct calculation avoids this order-of-limits issue by summing flavors
first. In A₂, each pair of Hamiltonian monomials is weighted by its
intersection size. At strict double scaling, pairwise intersections of
finitely many distinct p-subsets converge to independent Poisson(λ/2)
variables, with λ = 2p²/N and q = e^−λ. Higher shared intersections vanish
in this limit. For X with that Poisson law,

```text
E[X] = λ/2;   E[X(−1)^X] = −(λ/2)q;
E[X²] = λ/2 + λ²/4;   E[X²(−1)^X] = q(λ²/4−λ/2).
```

In units E[Tr H²]/dim = 1, the matching sum gives the full,
ensemble-mean-subtracted Majorana A₂ moments

```text
μ₀ = 4λ[λ(1+q) + 2(1−q)]
μ₂ = 16λ(1−q)[1−q² + 2λq].
```

At λ = 2 these give μ₀ = 32, μ₂ ≈ 42.141005 and μ₂/μ₀ ≈ 1.316906.
The positive μ₂ establishes dynamical weight in this channel. The
symbolic matching sum and closed forms agree; exact finite-N rational
Wick sums on two double-scaling lines provide a separate convergence
check. Extrapolation spreads are numerical sensitivity estimates, not
confidence intervals. Code: `duality/full_channel_chord.py`.

**The conserved component, resolved at N = ∞.** Mean subtraction
removes the identity, not all stationary weight: H² overlaps A₂, with
⟨A₂(H²−1)⟩ = 2λ(1−q) and ⟨(H²−1)²⟩ = 1+q. At strict double scaling the
remainder is exactly a matter chord of dimension 2:

```text
A₂ − ⟨A₂⟩ = [2λ(1−q)/(1+q)] (H² − 1) + ν O₂,   ν² = 8λ(1 − q² + 2λq)/(1+q)
μ₀ = 4λ²(1−q)²/(1+q)  [conserved]  +  ν²  [dynamical],   μ₂/ν² = 2(1 − q²).
```

The proof is a sign-reversing involution on the chord diagrams
([singlet-tower-note.md](singlet-tower-note.md) §3). At λ = 2 the
conserved share is 32.9% and the dynamical width is 1.963369; the
μ₂/μ₀ ≈ 1.316906 above is that width diluted. So the projected large-N
spectral density of this channel is the known Δ = 2 matter-chord
continuum. Independent checks: exact finite-N rational moments with the
H² component removed extrapolate to 1.96336; dense ED with the full
commutant of H projected out per realization (`projected_channel.py`)
gives 2.03, 2.03, 2.01 at p = 6, N = 16–20 against 2.00.

For n = 3..6 the same construction gives Aₙ as a polynomial in H plus
matter chords of dimension 2..n with polynomial dressings. Those are
exact identities verified in computer algebra on finite sets of words,
not proved. The chord new at level n is small (14.9%, 1.6%, 0.18%,
0.02% of the dynamical weight at λ = 2) and its norm scales as
λ^(2⌊n/2⌋) at small λ. Shares behave differently from norms, because the
dressings carry powers of the frequency: as λ → 0 the even chords fade
out of A₃..A₆ linearly in λ while the odd ones stay (A₃, A₄ → O₃;
A₅ → 99.5% O₃ + 0.53% O₅). The surviving dimensions, 2, 3, 5, match the
pattern of the large-q SYK bilinear spectrum; that is an observation on
four levels, not a derivation.

The same engine computes fixed-p finite-N moments without ED. For the
tested even p = 2,…,12, the coefficients support
μ₂/μ₀ = 8p(p−2)/N + O(N⁻²) in bandwidth units. Multiplying by
σ_H² = C(N,p)p!J²/N^(p−1) gives a limiting physical mean-square frequency
8p(p−2)J². This is an asymptotic statement, not a finite-N equality or a
proof of the general-p coefficient. The maximum of the identity-
subtracted ratio along fixed p (N = 24 at p = 4), and the N ≤ 20 ED rise
before it, reflect the conserved share of μ₀ falling from 85% to 19% to
zero: with the H² component removed the width has no maximum beyond the
smallest sizes and its physical value rises monotonically to 8p(p−2)J².

## 3. A separate exact self-averaging check

For j≠k and Bⱼₖ = Tr(ψⱼHψₖH)/dim, Gaussian Wick contraction gives

```text
Var(Bⱼₖ) = 4σ⁴ C(N−2,p−1),    σ² = p!J²/N^(p−1).
```

A nonzero trace requires the two coupling subsets to differ by exchanging
j and k; the two cross-Wick contractions supply the factor four.
At fixed p the absolute RMS is asymptotically
2p!J²/sqrt((p−1)!) × N^−((p−1)/2). Relative to E[Bⱼⱼ], the exponent
is (p+1)/2. The singlet statistic is different:
relVar(Tr H²/dim) = 2/C(N,p). These formulas establish concentration for
specified observables, not all aspects of emergent gauge symmetry.
Code: `syk-self-averaging/exact_wick.py`.

## 4. Interpretation and next steps

The most direct question for the proposal is how these full, appropriately
projected channels map to the Dirac CP-resolved bulk meson observables.
A boost generator's spectral density is not an invariant-mass spectrum;
the repository's free-particle boost-kernel examples illustrate that
obstruction but do not determine the interacting dictionary.

The contour/OTOC program is a separate investigation. The current fold
data test selected prescriptions and a converged range λ_chord≥0.5;
they do not prove failure of every Euclidean fold. Harlow–Zhao
[2607.14215v2, §6](https://arxiv.org/abs/2607.14215v2) explicitly makes
commutation of analytic continuation with the semiclassical limit an
assumption. Two-point detailed balance under an imaginary-time shift is
not a demonstration of a complete observer state or all-correlator KMS.

Priorities before a stronger claim:

1. Release the full test and reproduction manifests, and reconcile generated
   artifacts against the exact source commit.
2. Prove the matter-chord decomposition for general n (done for n = 2;
   verified, not proved, for n = 3..6), explain the λ^(2⌊n/2⌋) scaling of
   the new chord's norm, and derive the small-λ survival of the odd
   dimensions.
3. Compute the complete Dirac CP-resolved counterpart, with its ensemble
   specified. The disjoint-coupling CP-symmetric subensemble differs from
   the generic Dirac ensemble and that difference need not vanish at fixed λ.
   The Majorana route does not transfer directly: cᵢcᵢ† = 1 − nᵢ replaces
   ψᵢ² = 1, so the flavor sum leaves an occupation factor between the two
   halves of M₂ instead of removing the fermion marks.
4. Supply analytic degree bounds for interpolation claims across p, or keep
   them explicitly limited to the tested p values.
5. Obtain the state/time/operator dictionary needed to compare boost data
   with invariant masses. No amount of numerical precision alone supplies it.

## Reproduction and literature scope

See `REPRODUCIBILITY.md` and `scripts/validate.py`; the runner records the
working-tree hashes and executes the Makefile in a preserved snapshot.
The 't Hooft target is checked against FLZ [0905.2280](https://arxiv.org/abs/0905.2280).
The chord framework follows [1811.02584](https://arxiv.org/abs/1811.02584)
and [2407.09396](https://arxiv.org/abs/2407.09396). For large interaction
order, Gross–Rosenhaus [1702.08016, Eqs. 2.14–2.15](https://arxiv.org/abs/1702.08016)
distinguishes OPE coefficients cₙ = O(1/p) from their squared spectral
weights cₙ² = O(1/p²); these powers must not be interchanged or identified
with a finite-temperature full-channel expansion without a derivation.
