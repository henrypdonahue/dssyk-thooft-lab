# Numerical tests of DSSYK∞ and the 't Hooft proposal

This repository supplies reproducible calculations for the proposed
DSSYK∞–'t Hooft correspondence of Miyashita, Sekino and Susskind
([2506.18054v3](https://arxiv.org/abs/2506.18054v3),
[2607.05678v1](https://arxiv.org/abs/2607.05678v1)). It does **not**
establish or refute the duality, and it does not extract a matched meson
mass spectrum.

## Where to start

- [notes/author-facing-note.md](notes/author-facing-note.md): the
  technical note. Definitions, results, limitations, priorities.
- [notes/singlet-tower-note.md](notes/singlet-tower-note.md): the N = ∞
  structure of the singlet channels, with the proof and the checks.
- [notes/conserved-channel-note/](notes/conserved-channel-note): the
  LaTeX manuscript of the two theorems below.
- [road_map.md](road_map.md): what the paper claims, what is settled
  here, and the next work items.
- [ERRATA.md](ERRATA.md): every claim this repository has withdrawn or
  corrected, with what replaced it.

## The main result: the singlet channels at N = ∞

The paper's meson tower is to be read from the singlet operators
Mₙ = Σᵢ cᵢ†(d/dt)ⁿcᵢ. Their Majorana form is Aₙ = Σᵢ ψᵢ (ad_H)ⁿ ψᵢ. At
strict double scaling (N = ∞ with λ = λ_chord = 2p²/N fixed, q = e^(−λ),
normalized trace and E[Tr H²]/dim = 1) the repository establishes the
following.

**1. The pair amplitude is conserved (theorem).** The i ≠ j part of every
Aₙ channel has a chord image that is a polynomial in the chord
Hamiltonian T. A three-term recursion with scalar coefficients forces
this. So every spectral moment μ_k, k ≥ 1, of that part vanishes, at every
λ and every matter weight ([proof](notes/conserved-channel-note/main.tex),
[code](duality/chord_charges.py)).

**2. The full A₂ channel is a conserved piece plus a dimension-2 chord
(theorem).** The pair part is not the whole channel. Summing the flavor
first, the exact identity A₂ = 4ΣᵢHᵢ² makes A₂ the word H² with its two
monomials weighted by the indices they share, which are Poisson(λ/2) at
N = ∞. Then

```text
A₂ − ⟨A₂⟩ = [2λ(1−q)/(1+q)] (H² − 1) + ν O₂,   ν² = 8λ(1 − q² + 2λq)/(1+q)
μ₀ = 4λ²(1−q)²/(1+q)  [conserved]  +  ν²  [dynamical]
μ₂/ν² = 2(1 − q²)
```

O₂ is a matter chord of dimension 2, uncorrelated with every function of
H. At λ = 2, 32.9% of μ₀ is conserved and the dynamical width is
1.963369. The identity-subtracted ratio μ₂/μ₀ = 1.316906 is that width
diluted. The proof is a sign-reversing involution on chord diagrams
([note](notes/singlet-tower-note.md) §3, [code](duality/singlet_tower.py)).
Two checks do not use it. Exact finite-N rational moments with the H²
component removed extrapolate to 1.96336. Dense ED with the full
commutant of H projected out per realization gives 2.03, 2.03, 2.01 at
p = 6, N = 16–20, against 2.00.

**3. The levels above are matter chords too (verified, not proved).** For
n = 3..6, each Aₙ is a polynomial in H plus matter chords of dimension
2..n with polynomial dressings: zero residual in exact computer algebra
on every checked word. The chord that is new at level n is faint: 14.9%,
1.6%, 0.18% and 0.02% of the dynamical weight for n = 3..6 at λ = 2, with
norm ∝ λ^(2⌊n/2⌋) at small λ. As λ → 0 the even chords fade out of
A₃..A₆ and the odd ones remain (A₃, A₄ → O₃; A₅ → 99.5% O₃ + 0.53% O₅).
The surviving dimensions, 2, 3, 5, follow the pattern of the large-p
SYK bilinear spectrum h = 2, 3, 5, …; that is an observation, not a
derivation.

**What follows.** At N = ∞ every Aₙ channel is a conserved part plus
smooth matter-chord continua, with no discrete lines. Those continua are
in the frequency conjugate to H, which in the proposal generates boosts.
A particle of fixed mass already has a continuous boost spectrum
([duality/BOOST.md](duality/BOOST.md)), so turning these channels into
meson masses needs the state, time and operator dictionary. Small λ adds
one more point: the raw A₃..A₆ there all probe essentially the same Δ = 3
sector, so the natural "n-th member" of the tower is the dimension-n
chord, not the raw Aₙ ([note](notes/singlet-tower-note.md) §8).

**4. In the frame where the mesons should be, the mass spectrum is gapless
(Majorana form, N = ∞, small λ).** The paper's flat-space limit is λ
fixed with energies fixed in string units, J_MSS = pJ. Each chord of
dimension Δ becomes, as λ → 0, the quasinormal ladder |Γ(Δ + iω/a)|² at a
string-scale spacing (a = 2√2 J_MSS in this repo's Majorana translation).
At a Rindler horizon a field of mass M gives the entire density
K_{iω/a}(M/a)², and the ladder is exactly the superposition of those with
mass density ∝ M^(2Δ−1): a gapless, scale-invariant spectrum reaching
M = 0. With the exact decomposition, that fixes the channels' Rindler
mass spectrum as gapless power laws with Δ = 2, 3, 5, with no gap and no
discrete masses, unlike the gapped, discrete 't Hooft tower. That is
consistent with the paper's own picture of a deconfined, conformal
stretched horizon, and it means the meson tower is not visible in these
N = ∞ two-point functions ([duality/HORIZON.md](duality/HORIZON.md)).
At fixed p, N → ∞, the paper's calculational regime, the exact moments
approach the same ladders as p grows, with large closed-form finite-p
corrections ([syk-self-averaging/thooft_limit.py](syk-self-averaging/thooft_limit.py)).

## Results and their scope

| Calculation | Supported conclusion | Limitation |
|---|---|---|
| Exact Majorana singlet identities | A₁ = −2pH and A₂ = 4ΣᵢHᵢ² per realization, any N; a general mixed-word identity for every Aₙ (checked n ≤ 5) | The Dirac Mₙ need their own channel and CP treatment |
| Pair (i ≠ j) chord amplitude | Conserved at N = ∞ for every n, λ and matter weight (theorem) | Not the full flavor-summed channel |
| Full A₂ channel at N = ∞ | Conserved H² part plus a Δ = 2 matter chord (theorem), checked against exact finite N and projected ED | Majorana only; a boost-frequency continuum is not a mass |
| Aₙ, n = 3..6, at N = ∞ | Polynomial in H plus matter chords of dimension 2..n | Verified in computer algebra on finite word sets; no general-n proof |
| Singlet tower in a Rindler frame | As λ → 0 each chord is the quasinormal ladder \|Γ(Δ + iω/a)\|², the response of a gapless mass density ∝ M^(2Δ−1); the channels' Rindler mass spectrum has no gap and no discrete masses | Majorana only; free-field calibration at fixed distance from the horizon; a = 2𝒥 is the natural, not a supplied, identification |
| 't Hooft limit (fixed p, N → ∞) | Exact moments of A₂ (and the A₃ ratio) approach the λ → 0 ladders as p grows, with large closed-form finite-p corrections | Laws fit on low p, checked on higher p; the expected large-p ↔ small-λ consistency, and moments do not decide whether the finite-p spectrum is gapped |
| Fixed p, finite N | Exact rational moment laws; μ₂/μ₀ = 8p(p−2)/N + O(N⁻²) at each tested even p ≤ 12 | Not proved for all p; not an exact finite-N equality |
| Self-averaging | Exact finite-N variance of the adjoint Bⱼₖ and the singlet m₂, with an m₄ extension | Specific observables, not a theorem about full emergent gauge symmetry |
| 't Hooft target | 40 levels to 17–21 digits; all of FLZ's tables and sum rules reproduced | An accurate target does not supply the SYK/bulk dictionary |
| Complex-time OTOCs | DSSYK∞ scrambles at every temperature; the Euclidean fold fails at ED sizes and at N = ∞; a half-fake-period shift flips the sign semiclassically | Checked windows only; no universal fold no-go, observer dictionary or all-correlator KMS theorem |

### Fixed p

The proposal places the 't Hooft limit at fixed p, N → ∞, outside strict
double scaling. There the repository computes exact rational moments
without ED. For the tested even p, μ₂/μ₀ = 8p(p−2)/N + O(N⁻²) in units
E[Tr H²]/dim = 1, so the physical mean-square frequency tends to
8p(p−2)J² in the repository's coupling normalization. The identity-
subtracted ratio has a maximum along fixed p (N = 24 at p = 4). That
maximum is the conserved share of μ₀ changing. With the H² component
removed in exact arithmetic, the width falls monotonically from N = 12 at
p = 4 and from N = 16 at p = 6, and its physical value rises
monotonically toward 8p(p−2)J²
([code](syk-self-averaging/full_channel_laws.py)).
`syk-self-averaging/projected_channel.py` does the same split per
realization in ED, with degeneracies and tolerance sensitivity reported.

### Exact self-averaging statistics

For `Bⱼₖ = Tr(ψⱼHψₖH)/dim`, j ≠ k, with independent Gaussian couplings
of variance σ² = p! J²/N^(p−1),

```text
Var(Bⱼₖ) = 4σ⁴ C(N−2,p−1).
```

So the absolute off-diagonal RMS scales as
`2p! J² / sqrt((p−1)!) × N^−((p−1)/2)` at fixed p. Relative to the
mean diagonal observable its exponent is `(p+1)/2`. Separately,
`relVar(Tr H²/dim) = 2/C(N,p)`. These are distinct statistics. Along
the double-scaling line the violation beats the paper's e^(−a√N) for
every a; at fixed p it is only a power law
([code](syk-self-averaging/exact_wick.py),
[manuscript](notes/exact-wick-note/main.tex)).

### The operator dictionary

For Dirac fermions, with `d/dt = i[H,·]` and p/2 creators per Hamiltonian
monomial, `M₁ = −(p/2)iH` and `M₂ = −Σċ†ċ`. Anti-Hermiticity fixes the
phase of M₁; the coefficient follows from monomial counting. CP must be
the unitary particle-hole map, and for n ≥ 2 the raw Mₙ mix both CP
channels. The paper's λ = p²/N is half the chord literature's
λ_chord = 2p²/N ([dictionary](duality/dictionary.py),
[Dirac implementation](syk-self-averaging/dirac.py)).

## Reproduce

Use Python 3.9–3.12 in a virtual environment:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[dev]'
make test-all PY=python
make reproduce PY=python
make reproduce-slow PY=python
```

Scientific dependencies and direct test/plot dependencies are pinned in
[pyproject.toml](pyproject.toml). This is not a complete transitive lock
or an identical BLAS guarantee. CI covers Python 3.9 and 3.12 on Linux.

Run tests separately in each directory (the three modules and
`scripts/`), as the Makefile does; the code uses flat script imports.
`test` excludes slow tests (687 fast tests, 14 of them for the
runner); `test-all` includes them. `reproduce` runs the fast
deterministic producers and plots in a few minutes; `reproduce-slow`
runs every campaign and reference generator and takes many hours.
Producers overwrite their output files, so use a scratch copy or the
logged snapshot runner, which also compares regenerated artifacts with
the committed ones ([REPRODUCIBILITY.md](REPRODUCIBILITY.md)).

The ~1% ED agreement sometimes quoted refers to selected thermal
two-point benchmarks (0.6–1.6%; a test allows 5% plus statistical
error). It does not certify every chord observable to 1%, especially
folded four-point functions.

## Repository guide

- [thooft-target/](thooft-target): the reference spectrum (40 levels to
  17–21 digits), independent solvers, external anchors, eigenfunctions
  and the condensate.
- [syk-self-averaging/](syk-self-averaging): sparse and bitwise ED, exact
  Wick statistics, exact finite-N moment laws, finite-size and disorder
  campaigns.
- [duality/](duality): the chord engines, the singlet-tower
  decomposition, the dictionary, the boost-kernel calibration,
  large-interaction-order calculations and the contour/OTOC diagnostics.
- [notes/](notes): the technical note, the singlet-tower note and two
  LaTeX manuscripts.
- [scripts/validate.py](scripts/validate.py): the logged snapshot runner.
- [papers/README.md](papers/README.md): source versions and retrieval.

## Open problems

The main theoretical gaps are the state/time/operator dictionary that
relates boost-frequency data to invariant masses, a proof of the
matter-chord decomposition for general n, and the full Dirac
CP-resolved channel. [road_map.md](road_map.md) Part III lists them with
the next computational steps.

Code and original documentation are [MIT licensed](LICENSE). Third-party
papers are excluded from version control and keep their own licenses.
