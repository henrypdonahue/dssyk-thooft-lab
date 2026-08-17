# The KMS/OTOC axis: discriminating the three competing duals

*Companion to `ANTISCRAMBLING.md` and `CHORD.md`. Sources read for this
note: Narovlansky–Verlinde 2310.16994, Blommaert–Mertens–Papalini
2404.03535, Cui–Kolchmeyer 2607.13665, Harlow–Zhao 2607.14215. Repo-side
numbers are asserted in `test_discriminator.py`.*

## The axis

Cui–Kolchmeyer: a dS observer's OTOC must **anti-scramble** (their c < 0).
The chaos bound forces c > 0 for any unitary system in a KMS state at the
dS temperature. So a dS hologram must break one of: unitarity, dS-KMS, or
another hypothesis of the bound. This repo's facts on that axis:

- Under the flat dictionary, DSSYK∞ **scrambles at every temperature**:
  Re a = −1 exactly, rate λ_L = 2πv/β.
- The Harlow–Zhao fold — the one proposed repair — **fails**: the folded
  connected part is not 1/N-suppressed at ED sizes, and at N = ∞ it grows
  toward the semiclassical corner (`chord_fold.json`).

## Where each program sits

**1. MSS flat-space (this repo's target, 2607.05678).** The hologram is
the T = ∞ maximal-entropy state — not KMS at the dS temperature, so the
no-go is dodged in letter. But the dictionary then owes the map that
produces a dS observer's thermal physics and its c < 0: the "tomperature"
map. Our data sharpens the debt: the map must flip the OTOC sign, and the
fold cannot be that map. **Status: open pressure. This is the proposal's
weakest axis.**

**2. Narovlansky–Verlinde doubled DSSYK (2310.16994).** Their operators
are dressed and shifted by ±i β_dS/4, with β_dS "adjusted... to match the
KMS property" (their §2.3). That engineers KMS at the dS temperature —
exactly the no-go's hypothesis. So NV must show their constrained dressing
breaks another hypothesis of the bound. **The discriminating computation,
unpublished by anyone: the OTOC sign of their dressed operators.** It
needs the two-sided constrained sector — adjacent to this repo's chord
engine, and the sharpest next computation on this axis.

**3. Sine-dilaton gravity (2404.03535).** Their bulk is a black hole:
"a smooth black hole with Hawking temperature equal to the fake
temperature" (their §1; the fake temperature is their Eq. (1.7),
β_BH = 2π/sin θ), with "all matter correlators... thermal at the fake
temperature," and DSSYK as a sub-maximal scrambler. That is qualitatively
what this repo measures: sub-maximal scrambling with Re a = −1 and rate
λ_L = 2π/β_fake. Honesty note: with β_fake defined as β/v, that identity
is definitional. The quantitative tie sin θ ↔ v in matched units is
untested here — a naive chord-saddle check does not close under the
β_chord bridge. Open item. CK pressure attaches only to this program's
further dS3 reinterpretation, not to the black-hole reading.

## A candidate map that works — semiclassically (`signflip_probe.py`)

One shift does flip the sign: displace the probe pair by **half the fake
period** (Δ = π/v in Euclidean angle) — an antipodal shift on the fake
disk, and geometrically what NV's ±i β/4 dressing does. Measured
(`signflip_probe.json`, seven assertions in `test_signflip.py`):

- **Closed form: the flip is exact.** The growth coefficient obeys
  a(Δ) = a(0)·exp(−ivΔ), so Re a goes −1 → +1 at Δ = π/v, at every
  temperature tested, and at a second crossed configuration to machine
  precision. The shifted 2-pt stays finite and positive (margin π/2 from
  the propagator pole), and the time-ordered 4-pt is unchanged in
  magnitude under the shift — only the crossed correlator flips.
- **The antipodal frame is thermal at the tomperature — exactly.** The
  shifted Wightman weights obey per-line detailed balance at
  β_eff = β + β_fake, at any λ (an exact spectral identity, verified to
  1e-12). At the T = ∞ hologram point this is β_eff = β_fake: the same
  shift that flips the OTOC sign makes the 2-pt sector exactly thermal
  at the fake temperature. One map, both of CK's dS-observer
  requirements.
- **Exact in λ: the flip is semiclassical.** The antipodal trend turns
  anti-scrambling for λ ≤ 1.7 and not for λ ≥ 1.8 (this configuration
  and window), while the unshifted control scrambles at every λ.
- **The phase period is set by v, not by the growth rate.** The measured
  magnitude exponent κ ≈ 0.83–0.88 exceeds v ≈ 0.47, yet the flip works
  with Δ = π/v and fails with Δ = π/κ everywhere: phase and magnitude
  rates of the growing mode decouple at finite λ. This is why no rate
  recalibration rescues strong coupling.

Honesty: a sign flip plus effective thermality under a contour shift are
*necessary* conditions for this candidate map, not a dictionary. The
algebra story — why the observer's operators live at the antipodal
point — is exactly NV's construction to complete. But this is the first
map tested here that produces both dS-observer signatures in the exact
theory anywhere in parameter space.

## Does it pan out? The deep campaign (`--deep`, `signflip_deep.json`)

Four attack tests, run at the actual hologram point (β = 0):

- **A semiclassical dS epoch, not eternal dS.** The antipodal frame
  anti-scrambles in a transient window that grows as λ falls
  (t·Jc ≈ 0.57 → 1.14 over λ = 1.0 → 0.3, consistent with ~½ ln(1/λ));
  scrambling reasserts afterward. The unshifted control scrambles from
  the start at every λ.
- **The MSS evasion is located.** The growth coefficient's phase winds
  around the fake circle while its magnitude stays constant: there is no
  Euclidean operator placement with the definite-sign strip structure
  the chaos-bound proof requires. The frame's thermality is detailed
  balance, not strip-analytic KMS — that is how c < 0 coexists with the
  theorem.
- **The costs, measured.** The shifted contour is *fold-adjacent*: its
  correlator shapes are spectral-edge-dominated (|G(0)| > 1, growing
  along the shift), match neither the continued-sech nor a Gibbs form
  (deviations grow toward semiclassics, the fold pattern), and show no
  fake-period recurrence. What survives edge domination is exactly the
  ratio identities: detailed balance (exact) and the trend-level flip.
- **Rate extraction is contaminated, and we tested why**
  (`antipodal_modes.py`). The natural mechanism — the unflippable *even*
  scramblon harmonic (rate 2κ) ending the epoch — is **refuted**:
  shift-circle Fourier analysis gives κ₂/κ₁ ≈ 0.9–1.2, not 2, with
  comparable harmonic weights. The exact shift-response has continuous
  band content, not v-harmonics (the harmonic tower is only the λ → 0
  skeleton — the same phase/magnitude decoupling seen above). What
  survives: a flipped positive component plus a smaller unflipped
  negative one wherever the flip operates. The dS-rate question stays
  open, now with its obstruction understood.

Net: the mechanism is real and semiclassical — an anti-scrambling epoch
with exact effective detailed balance, at Susskind's own tomperature —
but the observer frame is not a Gibbs-like state, and its correlators
inherit the edge sensitivity of folded contours. Whether a full dS
dictionary can live on such a frame is precisely the NV-algebra
question.

## Verdict

Measured, DSSYK∞ is a fake-temperature scrambler — sine-dilaton's reading,
qualitatively. Making it a *de Sitter* hologram requires a sign-flipping
state/time map. The fold is excluded in-model; the antipodal
half-fake-period shift now works semiclassically (above) and connects
directly to NV's dressing — whose dressed OTOC sign remains the
discriminating computation. This axis and the meson-tower axis are
independent: a proposal can pass one and fail the other.

*Caveat: NV and sine-dilaton statements are quoted from their papers, not
re-derived. The repo-side tested content is Re a = −1 and the fold
verdicts.*
