# The NV dressing is the antipodal shift: the identification note

*Companion to `DISCRIMINATOR.md` and `ANTISCRAMBLING.md`. Source read
for this note: Narovlansky–Verlinde, "Double-scaled SYK and de Sitter
Holography," [arXiv:2310.16994](https://arxiv.org/abs/2310.16994), cited
below by their equation numbers. `test_nv_dressing.py` asserts every
checkable statement. The conditional ones are labeled.*

`DISCRIMINATOR.md` names the discriminating computation on the dS axis:
the OTOC sign of NV's dressed operators. `road_map.md` notes that the
antipodal shift *is* their dressing geometrically, so the answer may
already sit in `signflip_probe.json`.

This note settles what can be settled at the contour level. It states
precisely what remains, and it draws the conditional consequence.

## 1. What NV's construction says

NV work in a doubled, constrained DSSYK with λ = 2p²/N and q = e^(−λ).
That is this repo's chord convention, so there is no factor-of-two
landmine on this axis. The pieces used here:

**Physical operators** (§2.3). Gravitationally dressed by a time integral
over the two sides, schematically ∫dt O^(1−D)_L(t) O^D_R(τ−t), their
Eq. 32.

**The ±iβ_dS/4 dressing** (their Eq. 33). Two variants of each physical
operator:

    O^+(τ) = e^(+β_dS H/4) O_phys(τ) e^(−β_dS H/4)
    O^−(τ) = e^(−β_dS H/4) O_phys(τ) e^(+β_dS H/4)

They describe it as "a 90° rotation along the thermal circle," with the
two types "sitting on opposite sides of the half thermal circle."

**Composition** (their Eq. 34). O^−(τ) = O^+(τ + iβ_dS/2). Two
quarter-turns make one half-turn.

**β_dS is fixed by KMS** (§2.3, verbatim): "β_dS is a free parameter for
now, that will be adjusted later on to match the KMS property of the
correlation functions in the de Sitter vacuum state |Ψ_dS⟩." In §3.3 they
set β_dS = 2π, so T_dS = 𝕁/2π in their 𝕁 = 1 units (their Eq. 10).

**The antipodal map** (their Eqs. 61–62). The second field species is the
antipodal involution of the first, φ_−(x) = φ_+(x^Ω). Restricted to the
pole, φ_−(τ) = φ_+(τ + iβ_dS/2). The half-period shift *is* their
antipodal map at the pole.

**Thermality** (their Eq. 69). The Unruh detector on dressed correlators
obeys detailed balance, Ṗ(E_i→E_j)/Ṗ(E_j→E_i) = e^(−β_dS ΔE).

## 2. The map to this repo's contours

This repo's probe (`signflip_probe.py`) displaces the crossed pair by Δ
in Euclidean angle on the thermal circle, with θ-period 2π against chord
time β_chord. The fake period is 2π/v in θ, so β_fake = β/v in matched
time units.

The identification is a dictionary of displacements. It is exact at the
level of where insertions sit on the contour, given one identification of
scales.

| NV ([2310.16994](https://arxiv.org/abs/2310.16994)) | this repo |
|---|---|
| β_dS (Eq. 10, fixed to 2π by KMS matching, §2.3/§3.3) | β_fake = β/v (chord units: `beta_fake_c` = β_c/v) |
| single dressing e^(±β_dS H/4), Eq. 33, "a 90° rotation" | quarter-fake-period shift Δ = π/(2v), giving a(Δ) = −i·a(0), a literal 90° rotation of the growth coefficient |
| O^+ against O^−, relative displacement iβ_dS/2, Eq. 34 | half-fake-period shift Δ = π/v, the sign-flipping antipodal shift |
| antipodal involution at the pole, Eqs. 61–62 | the same Δ = π/v displacement, the antipodal shift on the fake disk |
| detailed balance at β_dS, Eq. 69 | per-line detailed balance at β_eff = β + β_fake, which goes to β_fake at the hologram point β → 0 |

Two things are asserted in `test_nv_dressing.py`, at machine precision.

**The composition algebra.** Dress the two probe pairs by ±π/(2v), which
is NV's O^+/O^− assignment. That equals the one-sided half-period shift
identically, by translation invariance of the corrected `otoc_ratio`. And
it equals −a(0) by the rotation law a(Δ) = a(0)e^(−ivΔ). A single quarter
dressing gives −i·a(0), which is NV's 90° rotation realized on the growth
coefficient. Residuals are about 1e-15.

**Thermality at NV's temperature, exact in λ.** `detailed_balance`
re-derives β_fake_c = β_c/v and β_eff = β_c(1 + 1/v) from the bridge. The
per-line law holds at float precision, and the independent
contour-engine check holds at ≤ 2e-13 on every stored row (≤ 8e-15 on
the λ = 2 rows). At the hologram point
(βJ = 0.02) the relative gap (β_eff − β_fake)/β_fake is exactly
v ≈ 0.006, the absolute gap being β_c itself. The
dressed frame is thermal at the tomperature. That is NV's Eq. 69 with
β_dS ↔ β_fake.

## 3. What closes, and what does not

**What closes, and is tested.**

*Contour kinematics, exact.* NV's mixed (+,−) insertions displace the
probe pair by half the fake period. That displacement is this repo's
antipodal shift. This is exact bookkeeping, not asymptotics.

*The 2-point sector, exact in λ.* The shifted frame's detailed balance at
β_eff = β + β_fake holds at every λ. NV's KMS-matching step, the reason
their β_dS is 2π, is reproduced by the shift with β_dS ↔ β_fake.

*The 4-point sign, closed form only.* The flip a(0) → −a(0) under the
half-period displacement is exact in the large-q closed form at every
temperature. It survives in the exact theory semiclassically, λ ≤ 1.7 in
the tested window (`signflip_probe.json`).

**What does not close.**

NV's operators are dressed in the *doubled, constrained* theory. Their
Eq. 33 acts on the physical operators of Eq. 32, which carry the
two-sided time-integral dressing and live in the constraint-projected
sector. This repo's probe inserts ordinary one-sided matter chords on a
shifted contour.

So the identification assumes the constrained-sector projection acts
trivially on the crossed 4-point function's growth coefficient. To settle
that, you need the two-sided constrained chord sector, which is the
computation `DISCRIMINATOR.md` names as adjacent to the engine. Until it
is done, the statements above are about a contour realization of NV's
displacement algebra. They are not about NV's dressed operators
themselves.

There is also a unit-matching gap. β_dS ↔ β_fake is the natural
identification, because both are the KMS temperature the construction
engineers. But the quantitative bridge between NV's 𝕁 = 1 units and this
repo's β_chord units is untested. That is the same caveat
`DISCRIMINATOR.md` records for sine-dilaton's sin θ ↔ v.

## 4. The consequence, stated conditionally

Suppose the identification holds at the exact level. That is, the
constraint projection does not renormalize the crossed channel's growth
coefficient. Then everything in `signflip_probe.json` and
`signflip_deep.json` transfers to NV's dressed operators verbatim.

- **Their dressed OTOC anti-scrambles semiclassically.** A transient
  epoch with the dS-observer sign, its window growing monotonically with
  ln(1/λ) (t·Jc ≈ 0.52 → 1.10 over λ = 1.0 → 0.3, sub-logarithmic at the
  small-λ end), and scrambling coming back afterward. Not eternal dS.
- **Their KMS engineering evades the chaos bound by the measured
  mechanism.** The coefficient's phase winds around the fake circle at
  constant magnitude. So there is no operator placement with the
  definite-sign strip structure the MSS proof requires (`mss_strip` in
  `signflip_deep.json`). The frame's thermality is detailed balance, not
  strip-analytic KMS. The evasion is consistent, because the bound's
  hypothesis fails. The frame is not a Gibbs-like state: its 2-pt is the
  T = ∞ correlator continued to z = τ_s + it, the large-q
  |1/cos(𝒥z)|^(2Δ), which the exact theory matches ever better toward
  semiclassics, while the naive Gibbs form fails (`thermal_shape` in
  `signflip_deep.json`).
- **The flip is semiclassical only**, λ ≤ 1.7 in the tested window. If
  NV's proposal needs c < 0 at strong coupling too, the identification
  predicts it fails there. That is falsifiable.

Suppose instead the constrained sector changes the answer. Then the
difference is itself the discriminator. It would isolate exactly what the
projection contributes at 4-point level, with the unconstrained baseline
already measured here.

## 5. What this does to the discriminator axis

`DISCRIMINATOR.md`'s three-program table sharpens.

**NV (program 2).** The still-unpublished discriminating OTOC now has a
measured candidate value: anti-scrambling, semiclassical, transient. It
is conditional on one derivation, the constrained-sector step, rather
than on a full computation from scratch. The burden shifts. NV's
construction must either endorse the contour realization, and inherit the
numbers including the λ ≤ 1.7 boundary and the transient window, or
exhibit what the constraint changes.

**MSS flat-space (program 1).** Unchanged. The tomperature-map debt
stands. This note makes the candidate map's algebra pedigree explicit
rather than folkloric.

**Sine-dilaton (program 3).** Unchanged, but note the parallel. Both
quantitative unit bridges, their sin θ ↔ v and NV's β_dS ↔ β_fake, are
the same kind of open item.

The spec framing survives. These are the conditions a dS dictionary must
meet, and the antipodal contour shift is the one candidate that meets the
necessary ones semiclassically. The constrained-sector computation is
still the sharpest next step on this axis. This note reduces it from
"compute the NV OTOC" to "check one projection property against a
measured baseline."
