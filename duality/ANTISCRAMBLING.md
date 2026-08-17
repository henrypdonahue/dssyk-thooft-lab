# Anti-scrambling vs DSSYK∞: the sign test and the Euclidean fold

*Code: `antiscrambling.py` (verified closed-form transcription),
`syk-self-averaging/fold_bench.py` (ED lab), `chord_fold.py` (N = ∞).
Verdicts asserted in `test_antiscrambling.py` and `test_chord_fold.py`.*

## Why this exists

Two July 2026 papers made the OTOC sector a sharp external constraint:

- **Cui–Kolchmeyer** ([2607.13665](https://arxiv.org/abs/2607.13665)):
  dS shockwaves produce a time *advance*, so a dS observer's OTOC must
  have coefficient c < 0 ("anti-scrambling"). AdS/black holes have c > 0.
  The chaos bound forces c > 0 for any unitary system in a KMS state at
  the dS temperature. The target paper's T = ∞ hologram is not such a KMS
  state, so it dodges the letter of the no-go — but then its dictionary
  must produce c < 0 some other way.
- **Harlow–Zhao** ([2607.14215](https://arxiv.org/abs/2607.14215)):
  a proposed mechanism. Systems with spectra bounded above and below allow
  **Euclidean-folded** correlators (their Eq. 6.9), and the fold flips the
  relevant sign — if the folded correlator equals the analytic
  continuation of the semiclassical answer. They pose the test: *"study
  the validity of this assumption in a concrete model such as the SYK
  model."* This repo ran it.

## Result 1: DSSYK∞ scrambles — the flat dictionary cannot anti-scramble

Streicher's exact large-q closed forms (all temperatures, transcription
verified; bridge 𝒥 = √2 p, ED-locked) give, at the symmetric crossed
configuration:

> F/F_d grows as a(v)·exp(v θ_L), with a(v) = −exp(−iπv/2)/cos(πv/2),
> so **Re a(v) = −1 exactly, at every temperature**,
> rate λ_L = 2πv/β → 2𝒥 at T = ∞.

That is the scrambling sign (c > 0), unsuppressed at the hologram point.

What ED adds, and what it does not: ED confirms the *growth* at every
measured (N, βJ), beyond 100σ. The *sign* is analytic
(Streicher-transcribed) and rides on a crossed-operator ordering fixed
empirically — ED does not measure the sign independently.

**Consequence.** With boundary time read as observer time, DSSYK∞ has the
wrong OTOC sign to be a dS observer, at every temperature. The
"tomperature" gap sharpens from "unstated" to "must flip the sign."

## Result 2: the fold works mechanically; the large-N question opened

`fold_bench.py` evaluates the HZ fold (6.9, n = 1). At finite N the
ordered-trace correlator is entire in the four angles, so the fold is
unambiguous: freeze the operator order, move the arguments. Findings at
p = 6, N = 12–16:

1. **Mechanics work.** Every folded correlator is finite (the bounded
   spectrum — HZ's requirement).
2. **The fold response tracks the continuation in form.** The unfolded ED
   OTOC is exactly real; the folded one develops a negative, monotone
   imaginary part — matching the continued closed form to ~10% at the
   best point, and in sign and trend everywhere.
3. **But the folded Euclidean 4-pt drifts away from the continued form as
   N grows** (ratio 1.05 → 1.21 → 1.42 for N = 12 → 16), while unfolded
   observables converge. λ-artifact, or a real obstruction?

## Decider #1: a genuine N-obstruction (`fold_bench.py --scan`)

p = 4, N = 8–20 at βJ = 2, plus a matched-λ anchor (`fold_scan.json`):

- The unfolded control is N-stable (ED/closed 1.80 → 1.84). The folded
  ratio explodes (1.52 → 8.19), same points, same instances.
- The drift grows as λ shrinks — opposite to a semiclassical correction —
  and at matched λ = 2.0 the larger-N point drifts further
  (1.700 ± 0.004 vs 1.522 ± 0.027).
- The folded connected fraction reaches O(1): −0.37 → −0.79 over
  N = 8 → 20. The 1/N hierarchy itself collapses under the fold.

Stated cap: dense ED stops at N = 20, so a large-N crossover could not be
excluded from ED alone. That question went to the chord engine.

## Decider #2: the fold fails at N = ∞ — the obstruction is structural

`chord_fold.py` evaluates the fold exactly in the double-scaled theory,
mirroring `fold_bench.py` point for point (same angles, same fold; the
transposition minus is *derived* here as the fermionic chord-crossing
sign). Scan λ = 4.0 → 0.5, truncation-doubled, drifts ≤ 2e−5. Verdict:

1. **The control behaves.** The unfolded ratio is λ-stable
   (−1.649 … −1.677), exactly real (manifest for Euclidean contours in the
   chord representation; the growth branch stays real too), and ~10% from
   matched-λ ED — the expected finite-N gap.
2. **The folded connected part is not suppressed at N = ∞.** As λ → 0 the
   unfolded connected fraction dies like λ; the folded one stays O(1) and
   grows (−0.21 → −0.33 over λ = 2.7 → 0.5). The N-scaled folded ratio
   runs away from the continued form like 1/λ: −1.86 → −21.3 against a
   constant −1.92. The folded imaginary branch keeps the predicted sign
   but its magnitude runs away the same way (about 1× at λ = 4, 3.5× at
   λ = 2, 40–50× at λ = 0.5). **The continuation assumption fails in the
   double-scaled theory itself, on both branches.**
3. **Why finite N converges so badly on folds — measured.** At matched λ,
   ED folded values sit 1.5–1.6× above the exact N = ∞ ones while
   unfolded gaps stay ≤ 10%. `fold_edge_probe.py`: the per-instance
   folded fraction tracks the sample's *upper spectral edge* at
   correlation −0.86 (unfolded control +0.15; lower edge −0.23), with
   35× the scatter. Folded correlators at ED sizes measure edge
   fluctuations — a finite-N effect the sharp-edged chord theory lacks.

Stated cap: folded contours below λ ≈ 0.5 exceed float headroom (the
λ = 0.3 point fails truncation-doubling and is excluded; the committed
fold drifts grow 2e-15 → 2e-5 over λ = 4 → 0.5). The λ-trend over
[0.5, 4] is unambiguous.

Consequence for HZ: in the one bounded-spectrum, exactly solvable
candidate model, folded correlators are finite (their requirement holds)
but do not approach the continued semiclassical answer. The fold does not
turn DSSYK∞ into an anti-scrambler.

Also run: the negative-temperature variant (HZ 6.12) is the H → −H
reflection. The mean-zero coupling measure is J → −J symmetric, which
flips H for every even p, so the disorder-averaged (6.9) and (6.12)
prescriptions coincide here.

## Standing caveats

Large-q closed forms at ED sizes carry a measured O(1/p, λ) floor (2-pt
deviation ratios ~0.3–0.55), so all finite-N 4-pt comparisons are
sign/trend-level, read against same-size unfolded calibration. The
connected ratio includes the disorder-replica piece (that is what the
ladder resums). The crossed-path sign convention is fixed empirically
(the alternative is off by 25×), and the fold configuration in the JSON
is asserted equal to `fold_map` of the crossed configuration.
