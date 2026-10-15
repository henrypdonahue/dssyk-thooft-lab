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
where the mesons are supposed to be, and whether a massive meson could be
hiding in it.

## Results

**1. The orthogonalized tower is a pure chord.** The chord that is new at
level n enters Aₙ with a constant dressing and coefficient exactly 4. So
removing from Aₙ its overlap with functions of H and with every dressed
lower chord leaves Õₙ = 4Ôₙ, a matter chord of dimension n with weight
16 Nₙ (checked against the sector moments, n = 2..5). Its share of Aₙ at
λ = 2 is 100%, 14.9%, 1.6%, 0.18% for n = 2..5.

**2. Each chord is a quasinormal ladder.** As λ → 0 the T = ∞ two-point
function of a chord of dimension Δ tends to sech^(2Δ)(𝒥t). Its spectral
density is

    ρ_Δ(ω) ∝ |Γ(Δ + iω/a)|²,     a = 2𝒥,

with poles at ω = −ia(Δ + k): a ladder with spacing a and lowest decay
rate Δa. Against the exact N = ∞ density, the maximum deviation is
(Δ/4)λ + O(λ²) for every Δ = 2..6 (dev/λ = 0.504, 0.757, 1.011, 1.266,
1.522 at λ = 0.05). In time, the exact correlators' decay rate
−d ln G/dt plateaus at Δa: 2.00 and 3.00 (units of a) for Δ = 2, 3 at
λ = 0.03.

**3. Full channels inherit the ladders.** The dressing of a sector is a
polynomial in ω at small λ, so a full channel's density is a sum of
polynomial multiples of the ladders of its sectors. Fitting the exact A₃
and A₄ densities to that form leaves a residual that falls about tenfold
per halving of λ: 1.7e-2, 1.9e-3, 1.6e-4, 1.0e-5 for A₃ over
λ = 0.4 → 0.05.

**4. The ladder is string-scale.** MSS measure energies in units of
J_MSS = pJ, with this repo's Var(J) = p!J²/N^(p−1) (their §2), and call
J_MSS the inverse string length. Then 𝒥 = √2 J_MSS exactly, and the
ladder spacing is a = 2√2 J_MSS. That is the scale of the meson masses
and of the stretched horizon of their §6, which sits a string length
from the Rindler horizon with a string-scale local Unruh temperature.
(The fundamental fermion, Δ = 1/p, decays at a/p, the cosmic scale of
Rahman–Susskind [2401.08555](https://arxiv.org/abs/2401.08555).)

## The calibration: what a mass looks like at a Rindler horizon

The T = ∞ correlator Tr(O(t)O)/dim is real, even, and regular at t = 0.
So it compares with two-sided correlators. For a probe at Rindler radius
1/a:

| bulk content | two-sided correlator | spectral density | poles |
|---|---|---|---|
| scale-invariant operator, dimension Δ | cosh^(−2Δ)(aτ/2) | \|Γ(Δ + iω/a)\|² | ladder −ia(Δ + k) |
| free field of mass M | K₀((2M/a) cosh(aτ/2)) | K_{iω/a}(M/a)² / (πa) | none (entire in ω) |

The second row is DLMF 10.32.17 (checked to 14 digits; the density is
also computed independently by quadrature). The massive density has real
zeros, the nodes of BOOST.md, and its decay rate grows without bound,
(M/a) sinh(aτ/2) + 1/4 in units of a. A gapped tower coupled with
positive weights is a positive sum of such terms. The 't Hooft model at
m_q² = 0 is gapped (FLZ), and with the repo's own masses and decay
constants its rate tracks its lightest meson: 13.9a at aτ = 8 for
M₀ = a/2. The one-sided responses differ only by the entire factor
e^(πω/a), so the dichotomy does not depend on that choice.

**Reading the tower against it.** With a = 2𝒥, the N = ∞ singlet
channels are, as λ → 0, exactly the scale-invariant responses with
Δ = 2, 3, 5. The exact decomposition says Aₙ contains nothing but a
polynomial in H and matter chords. Each chord's density converges to its
ladder. So the channels have no entire, K²-type component for any M > 0.
The decay rates show it directly: they plateau at Δa, where any massive
field's rate keeps growing. A least-squares fit in frequency is a
consistency check only: a massive density also fits the ladder basis to a
few 1e-3 on a finite window, but at a fixed floor, while the singlet
channels go below it (A₃ at λ = 0.05 is 200 times below the best massive
floor).

In de Sitter language the same statement reads: real Δ > 1 means
M²ℓ² = Δ(1 − Δ) < 0 in dS₂ (−2, −6, −20 for Δ = 2, 3, 5), or
m² = 4Δ(1 − Δ) < 0 in the dS₃ convention of Narovlansky–Verlinde
[2310.16994](https://arxiv.org/abs/2310.16994) Eq. 3. A massive particle
with Mℓ > ½ has complex Δ = ½ + i√(M²ℓ² − ¼), and quasinormal modes
with real part ≈ ±M. The tower has none.

## What it means for the proposal

At N = ∞ and small λ, the singlet channels of DSSYK∞ carry the analytic
structure of scale-invariant operators of dimension 2, 3 and 5 at the
string-scale ladder spacing, and no massive component. That fits the
picture of [2506.18054](https://arxiv.org/abs/2506.18054) §6, in which
the stretched horizon is the hot, deconfined, effectively conformal region
and the mesons are what escapes it. It also means that, if the duality
holds, the mesons are not visible in these channels at N = ∞. The
prediction "the singlet spectrum at small λ matches the 't Hooft
spectrum" then has to refer to something other than the N = ∞ two-point
functions of the Sₙ: finite-N corrections, operators built to reach past
the stretched horizon, or a dictionary that is not a local probe at fixed
distance from the horizon. Naming which is the dictionary's job.

## Scope

- N = ∞ at fixed λ, then λ → 0. At finite λ the spectrum is bounded,
  each correlator eventually decays as a power law, and the decay rate
  turns over (at λ = 0.5 it peaks near 1.9a for Δ = 2).
- Majorana Aₙ, n ≤ 6 (the decomposition is a theorem at n = 2 and
  verified in computer algebra for n = 3..6).
- The calibration is for free fields probed at a fixed distance 1/a from
  the horizon. A probe smeared down to the horizon, radius ρ ≪ 1/M, sees
  a conformal window up to τ ≈ (2/a) ln(1/(Mρ)) before the mass shows.
- a = 2𝒥 (ladder spacing equals acceleration) is the natural
  identification, not one the proposal supplies.

## Reproduce

```bash
cd duality
python horizon_tower.py                # all tables, about 30 s
python -m pytest -q test_horizon_tower.py
```
