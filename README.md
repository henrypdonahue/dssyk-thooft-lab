# Numerical tests of DSSYK∞ and the 't Hooft proposal

This repository supplies reproducible calculations relevant to the proposed
DSSYK∞–'t Hooft correspondence of Miyashita, Sekino and Susskind
([2506.18054v3](https://arxiv.org/abs/2506.18054v3),
[2607.05678v1](https://arxiv.org/abs/2607.05678v1)). It does **not** establish
or refute the duality, and does not extract a matched meson mass spectrum.

Start with the [technical note](notes/author-facing-note.md) for definitions,
results and limitations, and the [singlet-tower note](notes/singlet-tower-note.md)
for the N = ∞ structure of the singlet channels with its proof; see
[road_map.md](road_map.md) for the larger program.

## Results and their scope

| Calculation | Supported conclusion | Limitation |
|---|---|---|
| Exact Majorana singlet identities | A₁ = −2pH and A₂ = 4ΣᵢHᵢ², with ψᵢ² = 1 | The Dirac Mₙ tower needs its own channel and CP treatment |
| Sequential two-flavor chord insertion | A three-term recursion makes its chord image a polynomial in T for every n | This leading pair amplitude is not the full flavor-summed correlator |
| Full, mean-subtracted Majorana A₂ channel | Nonzero strict-double-scaling moments from intersection-weighted chords, cross-checked with exact finite-N Wick sums. At N = ∞ the channel is exactly a conserved H² component plus a matter chord of dimension 2 (proved) | Majorana only; a matter-chord continuum in boost frequency is not a meson mass |
| Majorana singlet tower Aₙ, n ≤ 6 | An exact finite-N operator identity; at N = ∞ each Aₙ is a polynomial in H plus matter chords of dimension 2..n with polynomial dressings, the new chord at level n having norm ∝ λ^(2⌊n/2⌋) | n ≥ 3 is verified in exact computer algebra on finite word sets, not proved; no general-n proof; not the Dirac Mₙ |
| Self-averaging | Exact finite-N variance for the specified Bⱼₖ observable and singlet m₂; an m₄ extension | Not a theorem for all observables or an individual realization's full emergent gauge symmetry |
| 't Hooft solver | Literature-anchored spectrum, sum rules and convergence checks | An accurate target spectrum does not supply the missing SYK/bulk dictionary |
| Complex-time OTOCs | Contour-dependent sign and fold diagnostics in checked parameter windows | No universal fold no-go, full observer dictionary or all-correlator KMS theorem |

### The full A₂ channel is not frozen

Use λ = λ_chord = 2p²/N, q = exp(−λ), normalized trace and
E[Tr H²]/dim = 1. With the ensemble mean removed, the strict-double-scaling
intersection-weighted chord calculation gives

```text
μ₀ = 4λ[λ(1+q) + 2(1−q)]
μ₂ = 16λ(1−q)[1−q² + 2λq].
```

At λ = 2, μ₀ = 32 and μ₂/μ₀ ≈ 1.316906. The positive second moment
rules out zero propagating weight for this full channel. It does not
establish a Regge tower or the Dirac M₂ result.
[Code and checks](duality/full_channel_chord.py),
[flavor-split discussion](duality/CHORD.md).

The mean-subtracted μ₀ contains conserved weight. At strict double
scaling the split is exact
([proof and checks](notes/singlet-tower-note.md),
[code](duality/singlet_tower.py)):

```text
A₂ − ⟨A₂⟩ = [2λ(1−q)/(1+q)] (H² − 1) + ν O₂,   ν² = 8λ(1 − q² + 2λq)/(1+q)
μ₀ = 4λ²(1−q)²/(1+q)  [conserved]  +  ν²  [dynamical]
μ₂/ν² = 2(1 − q²).
```

O₂ is a matter chord of dimension 2, uncorrelated with every function of
H. At λ = 2 the conserved share of μ₀ is 32.9% and the dynamical width is
1.963369; the ratio 1.316906 is that width diluted. The statement is a
theorem within the double-scaled chord rules. It is checked without the
proof against exact finite-N rational moments with the H² component
removed (1.96336 at λ = 2) and against dense ED with the full
commutant of H projected out per realization (2.03, 2.03, 2.01 at p = 6,
N = 16–20, against 2.00). The dynamical spectral density is therefore the
known Δ = 2 matter-chord continuum: smooth, with no discrete lines.

For n = 3..6 the same construction gives each Aₙ as a polynomial in H
plus matter chords of dimension 2..n with polynomial dressings. Those
decompositions are exact identities verified in computer algebra on
finite sets of words (zero residual), not proved. The chord that is new
at level n is small: 14.9%, 1.6%, 0.18% and 0.02% of the dynamical
weight for n = 3, 4, 5, 6 at λ = 2, with norm ∝ λ^(2⌊n/2⌋) at small λ.
As λ → 0 the even chords' shares of A₃..A₆ vanish linearly and the odd
ones remain (A₅ → 99.5% O₃ + 0.53% O₅), so the surviving dimensions are
2, 3, 5; the resemblance to the large-q SYK bilinear spectrum
h = 2, 3, 5, … is an observation, not a derivation.

**Do not quote the earlier claim that the meson tower is invisible at every
order of double scaling.** Taking a leading per-flavor limit before the
flavor sum loses terms that survive the sum. The pair contribution is a
cross-correlator, not a separately positive spectral measure. Its negative
moments are compatible with positivity of the full autocorrelator.

At fixed p, exact rational calculations for the tested even p values
support the A₂ bandwidth-width asymptotic
`μ₂/μ₀ = 8p(p−2)/N + O(N⁻²)`. This is not an exact finite-N equality,
and checks at finitely many p do not prove the coefficient for all p.
The corresponding physical mean-square frequency approaches
`8p(p−2) J²` in the repository's coupling normalization. The maximum of
the identity-subtracted ratio along fixed p (N = 24 at p = 4) is the
conserved share of μ₀ changing: with the H² component removed in exact
arithmetic the width has no maximum beyond the smallest sizes (it falls
monotonically from N = 12 at p = 4 and from N = 16 at p = 6), and its
physical value rises monotonically toward 8p(p−2) J². The asymptote is
unaffected.

`syk-self-averaging/projected_channel.py` removes all stationary energy
blocks at finite N, including degeneracies. It reports stationary
weight, dynamical moments and tolerance sensitivity. It is an ED
diagnostic; the asymptotic statement is the N = ∞ split above, and the
two agree where they overlap.

### Exact self-averaging statistics

For `Bⱼₖ = Tr(ψⱼHψₖH)/dim`, j ≠ k, with independent Gaussian couplings
of variance σ² = p! J²/N^(p−1),

```text
Var(Bⱼₖ) = 4σ⁴ C(N−2,p−1).
```

Consequently the **absolute off-diagonal RMS** scales as
`2p! J² / sqrt((p−1)!) × N^−((p−1)/2)` at fixed p. Relative to the
mean diagonal observable its exponent is `(p+1)/2`. Separately,
`relVar(Tr H²/dim) = 2/C(N,p)`. These are distinct statistics.
[Derivation](syk-self-averaging/exact_wick.py).

For Dirac fermions, with `d/dt = i[H,·]` and p/2 creators per Hamiltonian
monomial, `M₁ = −(p/2)iH` and `M₂ = −Σċ†ċ`. Anti-Hermiticity fixes the
phase issue in M₁; the coefficient requires the monomial-counting identity.
[Conventions](duality/dictionary.py), [Dirac implementation](syk-self-averaging/dirac.py).

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
[pyproject.toml](pyproject.toml). This is not a complete transitive lock or
an identical BLAS guarantee. CI covers Python 3.9 and 3.12 on Linux.

Run tests separately in the three module directories, as the Makefile does;
the code uses flat script imports. `test` excludes slow tests; `test-all`
includes them. `reproduce` runs the fast deterministic producers and plots;
`reproduce-slow` runs all listed campaigns and reference generators and
can take many hours. Producer commands overwrite their output files: use a
scratch copy or the logged snapshot runner described in
[REPRODUCIBILITY.md](REPRODUCIBILITY.md). Regenerated JSON runtime fields
and floating-point last digits need not be byte-identical.

The ~1% ED agreement sometimes quoted refers to selected thermal two-point
benchmarks (a test allows 5% plus statistical error). It does not certify
all chord observables to 1%, especially folded four-point functions.

## Repository guide

- [thooft-target/](thooft-target): reference spectrum (40 levels to 17+
  digits), alternative solvers, external anchors, eigenfunctions and
  condensate calculations.
- [syk-self-averaging/](syk-self-averaging): sparse/bitwise ED, exact Wick
  statistics, flavor-resolved moments, finite-size and disorder campaigns.
- [duality/](duality): chord engines, dictionary, boost-kernel calibration,
  large-interaction-order calculations and contour diagnostics.
- [notes/](notes): derivations, the concise author-facing note and the
  singlet-tower note.
- [papers/README.md](papers/README.md): source versions and retrieval.

The main theoretical gaps are the observer/state/time dictionary, the
relation of boost-frequency data to invariant masses, a proof of the
matter-chord decomposition for general n, and the complete Dirac
CP-resolved full-channel calculation. Numerical work remains useful
before those gaps close: independent implementations, convergence and
uncertainty studies and reproducible release artifacts. See the
prioritized review checklist in the technical note.

Code and original documentation are [MIT licensed](LICENSE). Third-party
papers are excluded from version control and retain their own licenses.
