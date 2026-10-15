# The KMS/OTOC axis: telling the three duals apart

*Companion to `ANTISCRAMBLING.md` and `CHORD.md`. Sources read for this
note: Narovlansky–Verlinde
[2310.16994](https://arxiv.org/abs/2310.16994),
Blommaert–Mertens–Papalini
[2404.03535](https://arxiv.org/abs/2404.03535), Cui–Kolchmeyer
[2607.13665](https://arxiv.org/abs/2607.13665), Harlow–Zhao
[2607.14215](https://arxiv.org/abs/2607.14215). `test_discriminator.py`
asserts the repo-side numbers.*

## The axis

Cui–Kolchmeyer set up a clean fork. A de Sitter observer's OTOC must
anti-scramble, c < 0. The chaos bound forces c > 0 for any unitary system
in a KMS state at the dS temperature.

So a dS hologram has to break one of three things. Unitarity. dS-KMS. Or
another hypothesis of the bound.

Two facts from this repo sit on that axis:

- Under the flat dictionary, DSSYK∞ **scrambles at every temperature**.
  The growth coefficient is a(v) = −1/cos(πv/2), real and negative,
  tending to −1 at the T = ∞ hologram point. The rate is λ_L = 2πv/β.
  |a| grows with β, so the sign is not marginal anywhere.
- The Harlow–Zhao fold is the one proposed repair, and it **fails**. The
  folded connected part is not 1/N-suppressed at ED sizes, and at N = ∞
  it grows toward the semiclassical corner (`chord_fold.json`). Neither
  statement needs a closed form.

## Where each program sits

### 1. MSS flat-space (this repo's target, [2607.05678](https://arxiv.org/abs/2607.05678))

The hologram is the T = ∞ maximal-entropy state. It is not KMS at the dS
temperature. So the no-go is dodged, in letter.

But then the dictionary owes the map that produces a dS observer's
thermal physics and its c < 0. That is the "tomperature" map. This
repo's data sharpens the debt: the map must flip the OTOC sign, and the
fold cannot be that map.

**Status: open pressure. This is the proposal's weakest axis.**

### 2. Narovlansky–Verlinde doubled DSSYK ([2310.16994](https://arxiv.org/abs/2310.16994))

They dress their operators and shift them by ±iβ_dS/4, with β_dS
"adjusted... to match the KMS property" (their §2.3).

That engineers KMS at the dS temperature, which is exactly the no-go's
hypothesis. So NV must break another hypothesis of the bound. The
discriminating computation is the OTOC sign of their dressed operators,
and it is still unpublished. It needs the two-sided constrained sector,
which is adjacent to this repo's chord engine.

`NV-DRESSING.md` now pins the identification at the contour level, with
equation-level citations (their Eqs. 33–34 and 61–62 against the
fake-period shifts, tested). Conditional on the constrained-sector
projection acting trivially on the crossed channel, the answer is
already measured here: anti-scrambling, semiclassical, transient.

### 3. Sine-dilaton gravity ([2404.03535](https://arxiv.org/abs/2404.03535))

Their bulk is a black hole. Their words: "a smooth black hole with
Hawking temperature equal to the fake temperature" (their §1). The fake
temperature is their Eq. (1.7), β_BH = 2π/sin θ. They report "all matter
correlators... thermal at the fake temperature," and they read DSSYK as
a sub-maximal scrambler.

That is qualitatively what this repo measures: sub-maximal scrambling,
a(v) < 0, rate λ_L = 2π/β_fake.

One honesty note. With β_fake defined as β/v, that identity is
definitional. The quantitative tie sin θ ↔ v in matched units is untested
here, and a naive chord-saddle check does not close under the β_chord
bridge. Open item. CK pressure attaches only to this program's further
dS3 reinterpretation, not to the black-hole reading.

## A candidate map that works, semiclassically

*Code: `signflip_probe.py`. Data: `signflip_probe.json`. Nine tests in
`test_signflip.py`.*

One shift does flip the sign. Displace the probe pair by **half the fake
period**, Δ = π/v in Euclidean angle. That is an antipodal shift on the
fake disk, and it is geometrically what NV's ±iβ/4 dressing does.

**The flip is exact in the closed form.** The growth coefficient obeys
a(Δ) = a(0)·exp(−ivΔ). With a(0) = −1/cos(πv/2) real, the antipodal shift
is a clean negation. Re a goes from a(0) to −a(0): −1.353 → +1.353 at
βJ = 2, and −1.030 → +1.030 at βJ = 0.5, tending to −1 → +1 at the T = ∞
hologram point. This holds at every temperature tested, to machine
precision.

Read that flip for what it is. Δ enters the growing mode's phase, so
Re a(Δ) = a(0)cos(vΔ), and half a period negates it. The algebra is
short. The content is in the next three items.

The shifted 2-point function stays finite and positive up to a *measured*
propagator pole at θ = π + π/v, with margin π/2 from the endpoint. The
time-ordered 4-point function does not change magnitude under the shift.
Only the crossed correlator flips.

**The antipodal frame is thermal at the tomperature, exactly.** The
shifted Wightman weights obey per-line detailed balance at
β_eff = β + β_fake, at any λ. The per-line law is exact by construction,
so it proves nothing on its own. The check with teeth is independent: the
shifted correlator from the contour engine, a code path with no Lehmann
weights, matches the per-line-shifted sum to ≤ 2e-13 on every stored row.
At the T = ∞ hologram point this is β_eff = β_fake.

One map, both of Cui–Kolchmeyer's dS-observer signatures.

**Exact in λ, the flip is semiclassical.** The antipodal trend turns
anti-scrambling for λ ≤ 1.7, and not for λ ≥ 1.8, in this configuration
and window. The unshifted control scrambles at every λ.

**The phase period is set by v, not by an effective growth rate.** A
three-point log-difference on the θ_L ≤ 2 window gives an effective
exponent κ_eff ≈ 0.83–0.88 for the exact chord OTOC. The same estimator
on the large-q closed form itself gives 0.94 on that window, and reaches
v ≈ 0.47 only for θ_L ≳ 6, because the cosh(vθ_L) growing/decaying pair
contaminates early differences. So κ_eff > v is an early-window artifact,
not a rate decoupling (the exact values sit below the closed form's own
early-window value, ratio 0.88–0.93), and an earlier version of this
note that read it as "phase and magnitude rates decouple" is corrected.
What stands: the flip works with Δ = π/v and fails with Δ = π/κ_eff
everywhere, and π/κ_eff ≈ 3.8 lands near π/(2v) = 3.3, where the
rotation law puts Re a = 0. No rate recalibration rescues strong
coupling.

**Honesty.** A sign flip plus effective thermality under a contour shift
are *necessary* conditions for this map. They are not a dictionary. The
algebra story, why the observer's operators live at the antipodal point,
is NV's construction to complete. But this is the first map tested here
that produces both dS-observer signatures in the exact theory anywhere in
parameter space.

## Does it hold up? The deep campaign

*Run with `--deep`. Data: `signflip_deep.json`. The window and shape
probes run at the actual hologram point, β = 0.*

**A semiclassical dS epoch, not eternal dS.** The antipodal frame
anti-scrambles in a transient window. The window grows as λ falls:
t·Jc ≈ 0.52 → 1.10 over λ = 1.0 → 0.3 (the maximum of the antipodal
curve, located by interpolation on a grid that now runs past every
window; the earlier 0.57 → 1.14 were grid points, one of them censored).
The growth is monotone in ln(1/λ) but not a fixed ½ ln(1/λ): the
pairwise slopes are 0.66, 0.38, 0.34, falling as λ falls, against a
quantization of ±0.07–0.14. Scrambling comes back afterward. The
unshifted control scrambles from the start at every λ.

**The MSS evasion is located.** The growth coefficient's phase winds
around the fake circle while its magnitude stays constant. There is no
Euclidean operator placement with the definite-sign strip structure that
the chaos-bound proof requires. The frame's thermality is detailed
balance, not strip-analytic KMS. That is how c < 0 coexists with the
theorem.

**The shapes, measured.** The shifted 2-pt is the T = ∞ correlator
continued to the complex time z = τ_s + it, so its large-q form is
|1/cos(𝒥z)|^(2Δ): sech in real time, the fake-period growth 1/cos(𝒥τ)
in Euclidean time, singular at the antipodal point itself. For t > 0 the
exact chord 2-pt *is* that form, with a deviation that shrinks toward
semiclassics (0.46 → 0.18 over λ = 2 → 0.25 at t ≥ 0.25/𝒥), and the
Euclidean 2-pt tracks the fake-period growth the same way (0.49 → 0.26).
At t = 0 the exact value is finite, |G(0)| > 1: the bounded chord
spectrum cuts off the sec(π/2) singularity, an edge effect, not a fold
pathology. The naive Gibbs guess, the observer thermal at β𝒥 = π in the
ordinary finite-temperature form, fails and fails worse toward
semiclassics (0.42 → 1.03). So the tomperature thermality is the exact
detailed-balance identity of the shifted sector, not a Gibbs state. An
earlier version of this note continued the real-time form with t and
τ_s swapped and read the resulting mismatch as "shapes fail worse toward
semiclassics"; that instrument was wrong, and the negative is withdrawn.

**The end of the epoch, and `antipodal_modes.py`.** The natural
mechanism is the unflippable *even* scramblon harmonic, at rate 2κ,
ending the epoch. On the physical half of the fake circle,
Δ ∈ [0, π/v], a least-squares harmonic fit finds a second harmonic that
is subdominant, |A₂/A₁| ≈ 0.3–0.6, and grows faster than the first,
κ₂/κ₁ ≈ 1.2–1.5 at λ = 1.5, 1.0, 0.6. The numbers depend on the model
(two against three harmonics) at the 10% level and the fits leave 5–10%
residuals: the response is not a clean sum of a few v-harmonics. So the
mechanism is neither refuted nor established; a faster-growing
unflippable component is there, at a rate between 1 and 2 times the
first harmonic's. An earlier version of this note reported the
mechanism refuted (κ₂/κ₁ ≈ 1, comparable weights) from a DFT over the
full fake circle; the samples past the 2-pt pole at Δ = π/2 + π/v were
two orders of magnitude larger than the physical ones and dominated
every DFT bin, so that instrument measured a spike, and the refutation
is withdrawn. What survives the two-mode fit is a flipped positive
component plus a smaller unflipped negative one, wherever the flip
operates. The dS-rate question stays open.

**Net.** The mechanism is real and semiclassical. It is an
anti-scrambling epoch with exact effective detailed balance, at
Susskind's own tomperature. But the observer frame is not a Gibbs-like
state, and its correlators inherit the edge sensitivity of folded
contours. Whether a full dS dictionary can rest on such a frame is
precisely the NV-algebra question.

## Verdict

Measured, DSSYK∞ is a fake-temperature scrambler. That is sine-dilaton's
reading, qualitatively.

Making it a *de Sitter* hologram needs a sign-flipping state and time
map. The fold is excluded in-model. The antipodal half-fake-period shift
works semiclassically, and it connects directly to NV's dressing. That
dressing's OTOC sign is still the discriminating computation.

This axis and the meson-tower axis are independent. A proposal can pass
one and fail the other.

*Caveat: the NV and sine-dilaton statements are quoted from their papers,
not re-derived. The repo-side tested content is the scrambling sign,
a(v) = −1/cos(πv/2) < 0, and the fold verdicts.*
