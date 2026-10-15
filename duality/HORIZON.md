# The singlet tower at the stretched horizon

*Code: `horizon_tower.py`. Data: `horizon_tower.json`. Tests:
`test_horizon_tower.py`. Builds on `singlet_tower.py` and
[BOOST.md](BOOST.md).*

## The question

Miyashita–Sekino–Susskind predict that "the spectrum of the DSSYK∞
singlets for small λ should match that of the 't Hooft model"
([2506.18054v3](https://arxiv.org/abs/2506.18054) §5.3), in a flat-space
limit that holds λ fixed as N → ∞ with energies fixed in string units
(their §2). `singlet_tower.py` gives the content of the singlet channels
in exactly that limit: a conserved polynomial in H plus matter chords of
dimension 2..n. This note asks what that content looks like in the frame
where the mesons are supposed to be, and what mass spectrum it
corresponds to there. Everything here is for the Majorana form Aₙ of the
paper's Sₙ.

## Results

**1. The orthogonalized tower is a pure chord.** The chord that is new at
level n enters Aₙ with a constant dressing and coefficient exactly 4 (the
T = {1, n} term of the Aₙ identity). So removing from Aₙ its overlap
with functions of H and with every dressed lower chord leaves Õₙ = 4Ôₙ, a
matter chord of dimension n with weight 16 Nₙ (checked against the sector
moments, n = 2..5). Its share of Aₙ at λ = 2 is 100%, 14.9%, 1.6%, 0.18%
for n = 2..5.

**2. Each chord is a quasinormal ladder.** As λ → 0 the T = ∞ two-point
function of a chord of dimension Δ tends to sech^(2Δ)(𝒥t). Its spectral
density is

    ρ_Δ(ω) ∝ |Γ(Δ + iω/a)|²,     a = 2𝒥,

with poles at ω = −ia(Δ + k): a ladder with spacing a and lowest decay
rate Δa. At the fixed scale 𝒥 the exact N = ∞ density deviates from it
by (Δ/4)λ + O(λ²) for Δ = 2..6. Almost all of that is a scale shift: the
best-fit scale is 𝒥(1 − (Δ/4)λ), and there the shape deviates by about
0.045λ for every Δ. So one common spacing a holds for all channels at
leading order only. In time, the exact decay rate −d ln G/dt plateaus at
Δa: 2.00 and 3.00 (units of a) for Δ = 2, 3 at λ = 0.03.

**3. Full channels are sums of ladders.** A dressing is a polynomial in ω
at small λ, so a channel's density is a sum of polynomial multiples of
its sectors' ladders. Fitting the exact A₃ and A₄ densities to that form
leaves a residual that falls about tenfold per halving of λ (1.7e-2 →
1.0e-5 for A₃ over λ = 0.4 → 0.05). This is a consistency check only.
Since |Γ(m + 1 + ix)|² = (m² + x²)|Γ(m + ix)|², the basis is just an even
polynomial times ρ₂, and a massive-field density also fits it on a finite
window, at a floor of about 5e-4 at best.

**4. The ladder is string-scale.** In the Majorana translation of the
paper's units (J_MSS = pJ with this repo's Var(J) = p!J²/N^(p−1)),
𝒥 = √2 J_MSS exactly, so a = 2√2 J_MSS. The O(1) factor depends on that
translation; the paper's complex model changes it. Either way the ladder
sits at the string scale: the scale of the meson masses, and of the
stretched horizon of their §6, which sits a string length from the
Rindler horizon with a string-scale local Unruh temperature.

**5. At finite λ the ladder ends in a power law.** The spectrum is
bounded, so the exponential decay cannot last. For Δ = 2 the exact decay
rate peaks at aτ ≈ 4.7–5.6, at a height that approaches Δa as λ falls
(1.81a, 1.87a, 1.94a over λ = 1, 0.5, 0.3); for Δ = 3 the peaks sit at
aτ ≈ 4.3–4.9 and their height is not monotone in λ. Then a power-law tail
takes over. For λ ≥ 0.5 its exponent is −3.0 (−3.01 and −3.05 for
Δ = 2, −3.02 and −3.03 for Δ = 3), with drift up to 0.13 under doubling
the truncation. The tail is monotone: the density of states has
square-root band edges, and pairs of states near the same edge make the
spectral density non-analytic at ω = 0 (ω² ln|ω|). At λ = 0.3 the fitted
exponent moves by 0.6–0.9 under doubling the truncation, and below that
the crossover sits under 1e-16: not resolved.

## The calibration: mass spectra at a Rindler horizon

The T = ∞ correlator Tr(O(t)O)/dim is real, even and regular at t = 0.
That is the structure of a two-sided correlator: sech^(2Δ)(𝒥t) first
becomes singular at t = iβ_R/2, β_R = 2π/a. For a probe at Rindler
radius 1/a:

| bulk content | two-sided correlator | spectral density | poles |
|---|---|---|---|
| free field of mass M | K₀((2M/a) cosh(aτ/2)) | K_{iω/a}(M/a)² / (πa) | none (entire in ω) |
| mass density ρ(M) | ∫ρ(M) K₀((2M/a) cosh(aτ/2)) dM | ∫ρ(M) K_{iω/a}(M/a)² dM / (πa) | none if ρ is gapped |
| ρ(M) ∝ M^(2Δ−1), gapless | ∝ cosh^(−2Δ)(aτ/2) | ∝ \|Γ(Δ + iω/a)\|² | ladder −ia(Δ + k) |

The first row is DLMF 10.32.17. Its density has real zeros (the nodes of
BOOST.md), and its decay rate grows without bound,
(M/a) sinh(aτ/2) + ¼ tanh(aτ/2) in units of a. The last row is the
Mellin transform of K_ν²:
∫₀^∞ z^(2Δ−1) K_{ix}(z)² dz = 2^(2Δ−3) Γ(Δ)² |Γ(Δ + ix)|² / Γ(2Δ),
checked numerically. So a ladder is the response of a gapless,
scale-invariant mass spectrum reaching down to M = 0. Its poles come from
the M → 0 end, while any gapped spectrum gives an entire density.

**Reading the tower against it.** With a = 2𝒥, the exact decomposition
says each Aₙ is a polynomial in H plus matter chords, nothing else. Each
chord's density converges to its ladder. The map ρ(M) ↦ S(ω) is
invertible: the cosine representation turns S into the K₀-transform of
ρ, which is injective. So, as λ → 0, the channels' Rindler mass spectrum
is fixed: a power-law density M^(2m−1)(c₀ + c₁M² + …) with Δ = 2, 3, 5,
with no gap and no discrete masses. The 't Hooft model at m_q² = 0 is
gapped and discrete (FLZ). Its tower, built with the repo's masses and
decay constants, gives an entire density, and a decay rate that grows
like its lightest meson's (13.9a at aτ = 8 for M₀ = a/2).

The decay-rate plateau illustrates this but does not prove it. The
plateau at Δa is the slowest, M → 0 part, and a gapped admixture would
decay faster and hide underneath it. The exclusion comes from the
decomposition and the invertibility, not from the plateau.

In de Sitter language: real Δ > 1 means M²ℓ² = Δ(1 − Δ) < 0 in dS₂
(−2, −6, −20 for Δ = 2, 3, 5), or m² = 4Δ(1 − Δ) < 0 in the dS₃
convention of Narovlansky–Verlinde
[2310.16994](https://arxiv.org/abs/2310.16994) Eq. 3. A massive particle
with Mℓ > ½ has complex Δ = ½ + i√(M²ℓ² − ¼), with quasinormal modes
whose real part is ≈ ±M. The tower has none.

## In the paper's calculational regime

MSS do the 't Hooft comparison at fixed p, N → ∞, not at fixed λ.
`syk-self-averaging/thooft_limit.py` takes the exact fixed-p laws to that
limit. The channels' physical moments approach the ladders as p grows,
which is the expected large-p ↔ small-λ correspondence and a consistency
check of it. The new content is the exact size of the finite-p
corrections, and they are large:
- A₂'s ⟨ω²⟩ is (1 − 2/p) times its ladder value, half of it at p = 4.
- A₂'s kurtosis is 7/2 − 1/(2(p − 2)).
- A₃'s ⟨ω⁴⟩/⟨ω²⟩ grows toward the Δ = 3 value 40p², but at the p
  computed it is still between the Δ = 2 and Δ = 3 values.

So at finite ᾱ = p² the 't Hooft-limit channels are the ladders only
roughly. The moments do not decide whether the finite-p spectrum is gapped.

## What it means for the proposal

At N = ∞ and small λ, the Majorana singlet channels of DSSYK∞ have a
gapless, scale-invariant Rindler mass spectrum, dimensions 2, 3, 5 at the
string-scale ladder spacing, with no gap and no discrete masses. That
fits the picture of [2506.18054](https://arxiv.org/abs/2506.18054) §6,
in which the stretched horizon is the hot, deconfined, effectively
conformal region and the mesons are what escapes it. It also means that,
if the duality holds, the discrete meson tower is not visible in these
channels at N = ∞. The prediction "the singlet spectrum at small λ
matches the 't Hooft spectrum" then has to refer to something other than
the N = ∞ two-point functions of the Sₙ in this Majorana form: finite-N
corrections, the Dirac channels (not yet computed), operators built to
reach past the stretched horizon, or a dictionary that is not a local
probe at fixed distance from the horizon. Naming which is the
dictionary's job.

## Scope

- N = ∞ at fixed λ, then λ → 0; at finite λ the spectrum is bounded and
  correlators end in a t⁻³ tail (result 5).
- Majorana Aₙ, n ≤ 6 (the decomposition is a theorem at n = 2 and
  verified in computer algebra for n = 3..6).
- Free-field calibration with a probe at a fixed distance 1/a from the
  horizon. A probe smeared down to the horizon, radius ρ ≪ 1/M, sees a
  conformal window up to τ ≈ (2/a) ln(1/(Mρ)) before the mass shows.
- a = 2𝒥 (ladder spacing equals acceleration) is the natural
  identification, not one the proposal supplies; the spacing depends on
  Δ at O(λ) (result 2).

## Reproduce

```bash
cd duality
python horizon_tower.py                # all tables, a few minutes
python -m pytest -q test_horizon_tower.py
```
