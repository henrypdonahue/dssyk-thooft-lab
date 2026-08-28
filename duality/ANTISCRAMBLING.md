# Anti-scrambling and DSSYK∞: the sign test, and the fold

*Code: `antiscrambling.py` (closed forms),
`syk-self-averaging/fold_bench.py` (the ED lab), `chord_fold.py`
(N = ∞). Verdicts are asserted in `test_antiscrambling.py` and
`test_chord_fold.py`.*

## The constraint

Two papers from July 2026 turned the OTOC sector into a sharp external
test.

**Cui–Kolchmeyer** ([2607.13665](https://arxiv.org/abs/2607.13665)).
A de Sitter shockwave produces a time *advance*, not a delay. So a dS
observer's OTOC must have coefficient c < 0. They call it
anti-scrambling. AdS and black holes have c > 0.

The chaos bound then forces c > 0 for any unitary system in a KMS state
at the dS temperature. The paper's T = ∞ hologram is not such a state, so
it avoids the letter of the no-go. But then its dictionary owes c < 0
some other way.

**Harlow–Zhao** ([2607.14215](https://arxiv.org/abs/2607.14215)).
A proposed mechanism. A system with a spectrum bounded above and below
admits Euclidean-folded correlators (their Eq. 6.9). The fold flips the
relevant sign, *if* the folded correlator equals the analytic
continuation of the semiclassical answer.

They pose the test directly: "study the validity of this assumption in a
concrete model such as the SYK model." This repo ran it.

## DSSYK∞ scrambles

Streicher's exact large-q closed forms cover all temperatures. The bridge
𝒥 = √2 p is locked against ED. At the symmetric crossed configuration
ℱ/F_d grows like a(v)·e^(v θ_L), with

> **a(v) = −1/cos(πv/2)**, real and negative at every temperature.

|a| = 1.030 at βJ = 0.5, 1.353 at βJ = 2, and 3.194 at βJ = 8. It tends
to a = −1 at the T_B = ∞ hologram point, where the growth is not
suppressed at all. The rate is λ_L = 2πv/β, which goes to 2𝒥 at T = ∞.

That is the scrambling sign, c > 0. It holds at the hologram point, and
its magnitude *grows* with β. The closed form is checked against a
late-time fit of the correlator itself, to 1e-15.

ED measures the growth at every (N, βJ), beyond 100σ. The sign
is analytic. It is transcribed from Streicher, and it inherits a
crossed-operator ordering fixed empirically. ED does not measure the sign
independently. That is the one soft spot on this axis, and it is worth
knowing.

**What follows.** Read boundary time as observer time, and DSSYK∞ has the
wrong OTOC sign to be a dS observer, at every temperature. The
"tomperature" gap sharpens. It was "unstated." It is now "must flip the
sign."

### A note on transcribing Streicher

Streicher's crossed-OTOC display (his Eq. 3.5) is written *after* the
time-translation freedom is fixed. His Appendix A says it: "we use
time-translation invariance to set θ₃₄⁺ = 0."

So the φ₁₂⁺ in that display is the squeeze of θ₁₂⁺ − θ₃₄⁺, not of θ₁₂⁺.
The two agree in that gauge and nowhere else. Using the gauge-fixed
display at general angles silently breaks time-translation invariance.
Anyone reusing these formulas should take the gauge-invariant
combination.

Three checks pin it, and all three are asserted:

1. **Translation invariance.** ℱ/F_d is a ratio of thermal correlators,
   so it depends only on time differences. The chord engine and ED are
   manifestly like that. The gauge-fixed form drifts by a factor of about
   2 under a rigid 0.5 shift. The invariant one does not move.
2. **A structural identity.** At the coincidence θ₂ = θ₃ the crossed and
   uncrossed branches must differ by exactly −2. The invariant form gives
   −2.000000 at every configuration and every v.
3. **The exact N = ∞ chord engine.** As λ → 0 the double-scaled theory
   *is* the large-q limit, so chord/closed must approach 1. Measured:
   0.964 → 0.979 over λ = 2 → 0.5, the committed grid, tracking the TOC
   control 0.683 → 0.945.

## How good is the closed form?

Quantitative, not merely sign-level.

| comparison | ratio |
|---|---|
| ED / closed, unfolded crossed, βJ = 2, N = 12 → 16 | 1.084 → 1.048 |
| ED / closed across the p = 4 scan, N = 8 → 20 | 1.063 → 1.083 (spread < 0.03) |
| N = ∞ chord / closed, λ = 4 → 0.5 | 0.973 → 0.979 |

So crossed-OTOC comparisons here are good to 5–8% against ED, and 2–4%
against the exact chord engine.

## The fold works mechanically. The connected part does not shrink.

`fold_bench.py` evaluates the HZ fold (6.9, n = 1). At finite N the
ordered-trace correlator is entire in the four angles, so the fold is
unambiguous. Freeze the operator order, move the arguments. At p = 6,
N = 12–16:

1. **The mechanics work.** Every folded correlator is finite. That is
   HZ's requirement, and it holds.
2. **The response has the right shape.** The unfolded ED OTOC is exactly
   real. The folded one develops a monotone imaginary part, and the
   continued closed form tracks its magnitude to within 7% at the largest
   N (|ED/closed| = 1.072, 1.052, 1.019 at N = 16).
3. **But the folded Euclidean 4-point function drifts away from the
   continued form as N grows**, while the unfolded observables converge.

### Why the continued form is a magnitude only

The fold sends the configuration far outside the window where Streicher's
expression is the principal branch. θ₁₂⁻ and θ₃₄⁻ go to −3π, and
(θ₁₂⁺ − θ₃₄⁺) goes to −5π, against a principal window of (0, 2π).

So the continued value at the fold is a **branch choice**. Its overall
sign is not a prediction. Both verdicts below therefore lead with the
closed-form-free observable, the connected fraction (F − F_d)/F_d. The
ratio to the continued form is a secondary, magnitude-only diagnostic.

### At ED sizes: an N problem, not a λ artifact

p = 4, N = 8–20 at βJ = 2, plus a matched-λ anchor (`fold_scan.json`):

- **No closed form used.** The folded connected fraction |(F − F_d)/F_d|
  grows monotonically from 0.365 to 0.786 over N = 8 → 20, *while*
  λ = 16/N shrinks. The 1/N hierarchy collapses under the fold, in the
  direction opposite to a semiclassical correction.
- The unfolded control is N-stable (ED/closed 1.063 → 1.083, spread
  < 0.03). So the fold is the culprit, not the closed form.
- Secondary. The |folded ratio| against the continued form explodes from
  4.28 to 23.05 over the same range. At matched λ the larger-N point
  drifts further: |−4.79| at (p=6, N=18) against |−4.28| at (p=4, N=8).
  That anchor is λ_chord = 4.0; `fold_scan.json` quotes it as
  λ = p²/N = 2.0.

Dense ED stops at N = 20, so ED alone could not exclude a large-N
crossover. That question went to the chord engine.

### At N = ∞: the problem is structural

`chord_fold.py` evaluates the fold exactly in the double-scaled theory.
It mirrors `fold_bench.py` point for point, with the same angles and the
same fold. The transposition minus sign is *derived*, as the fermionic
chord-crossing sign. Scan λ = 4.0 → 0.5, truncation-doubled, drifts
≤ 2e−5.

1. **The control behaves.** The unfolded ratio is λ-stable
   (−1.649 … −1.677) and exactly real, which is manifest for Euclidean
   contours in the chord representation. It sits 2–4% from the closed
   form and about 10% from matched-λ ED.
2. **The folded connected part is not suppressed, and you need no closed
   form to see it.** As λ → 0 the *unfolded* connected fraction dies like
   λ, from −0.208 to −0.026 over λ = 4 → 0.5. The *folded* one stays O(1)
   and grows, from −0.233 to −0.333. There is no 1/N suppression to be
   had at N = ∞. The continuation assumption fails in the double-scaled
   theory itself. Secondary, magnitude only: the N-scaled folded ratio
   against the continued form runs |2.73| → |31.26|, which goes like 1/λ.
3. **Why finite N converges so badly on folds.** At matched λ, ED folded
   values are 1.52× (p=4, N=8) and 1.62× (p=6, N=18) the exact N = ∞
   ones. The unfolded gaps are 1.10× and 1.03×. In
   `syk-self-averaging/fold_edge_probe.py` the per-instance folded
   fraction tracks the sample's *upper spectral edge* at correlation
   −0.86, with 35× the scatter. The unfolded control is +0.15, and the
   lower edge is −0.23. So folded correlators at ED sizes measure edge
   fluctuations. That is a finite-N effect, and the sharp-edged chord
   theory does not have it.

Cap. Folded contours below λ ≈ 0.5 exceed float headroom. The λ = 0.3
point fails truncation-doubling and is excluded. Committed fold drifts
grow from 2e-15 to 2e-5 over λ = 4 → 0.5. The λ-trend over [0.5, 4] is
unambiguous.

**For Harlow–Zhao.** In the one bounded-spectrum, exactly solvable
candidate model, folded correlators are finite. Their requirement holds.
But the connected part they need to be 1/N-suppressed is not suppressed,
at ED sizes or at N = ∞. The fold does not turn DSSYK∞ into an
anti-scrambler.

The negative-temperature variant (HZ 6.12) is the H → −H reflection. The
mean-zero coupling measure is symmetric under J → −J, which flips H for
every even p. So the disorder-averaged (6.9) and (6.12) prescriptions
coincide here.

## Standing caveats

The crossed-path sign convention is fixed empirically. The alternative is
off by 25×.

The connected ratio includes the disorder-replica piece. That is what the
ladder resums.

The fold configuration in the JSON is asserted equal to `fold_map` of the
crossed configuration.

At the folded angles the continued closed form is a branch choice. Only
its magnitude is read as a prediction. The verdicts rest on the
closed-form-free connected fractions.
