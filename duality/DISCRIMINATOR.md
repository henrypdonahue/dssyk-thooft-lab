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

## Verdict

Measured, DSSYK∞ is a fake-temperature scrambler — sine-dilaton's reading,
qualitatively. Making it a *de Sitter* hologram requires either a
sign-flipping state/time map that no current proposal supplies (the fold
is excluded in-model), or a demonstration that constrained dressing evades
the chaos bound (NV; computable, unpublished). This axis and the
meson-tower axis are independent: a proposal can pass one and fail the
other.

*Caveat: NV and sine-dilaton statements are quoted from their papers, not
re-derived. The repo-side tested content is Re a = −1 and the fold
verdicts.*
