# Errata

Every claim this repository has published and then withdrawn or
corrected is listed here, once. Each entry gives the earlier statement,
what replaced it and the commit. The other documents state only the
current result and cite the entry number (E1, E2, …) where a reader may
have seen the older version.

Dates are commit dates. Work committed on 2026-10-14 was done between
2026-09-01 and 2026-10-04.

## The singlet tower

**E1. "The meson tower is a 1/N effect, invisible at every order of the
double-scaled expansion, so the match must be made at fixed p, λ → 0."**
Withdrawn 2026-08-25 (`0525e64`). The conserved-channel theorem covers
the i ≠ j *pair* amplitude only. ⟨AₙAₙ⟩ = N⟨XᵢXᵢ⟩ + N(N−1)⟨XᵢXⱼ⟩, and
the two terms are the same order, so removing the pair part leaves a
propagating channel. The fixed-p argument stands separately, on the
self-averaging side (`exact_wick.py`).

**E2. "Whatever weight the tower carries at N = ∞ is entirely
flavor-diagonal."** Withdrawn when the full channel was solved
(committed 2026-10-14, `aa66ab8`). On a double-scaling line the diagonal
and pair halves are each O(√N) in the flavor-summed normalization and
cancel at that order. The physical channel is the O(1) remainder, which
belongs to neither half (`duality/full_channel_chord.py`).

**E3. "The full A₂ width at N = ∞ is ≈ 1.19 (or a Richardson estimate of
the finite-N series)."** Superseded 2026-10-14. The exact value is the
closed form μ₂/μ₀ = 1.316906 at λ_chord = 2. As a *width* of the moving
part, that ratio is itself superseded by 2(1 − q²) = 1.963369: about a
third of μ₀ is the conserved H² admixture (`duality/singlet_tower.py`).

**E4. "The tower channel turns on toward the λ → 0 corner"** (the raw
ED width rising along fixed p over N ≤ 20). Withdrawn 2026-10-14. That
window is the rising side of a peak at N = 24 (p = 4), and the peak is
conserved H² weight leaving the denominator. With the H² component
removed the width falls monotonically along fixed p, to 8p(p−2)/N
(`syk-self-averaging/full_channel_laws.py`).

**E5. "The ad-shift coincidence A(n=1, k) = A(n=2, k−2) holds at every
step of the ladder."** Refuted 2026-08-28 (`501bd34`). It holds exactly
at the n = 1 → 2 step and at no deeper one: at K = 10, p = 2 the
coefficients run 21 : 21 : 14 : 9 over n = 1..4.

**E6. "The double-scaling constant c₂ has no identified closed form"
and "A(n=3, p) has degree above 5 in p."** Superseded 2026-10-14 and
2026-08-28. c₂(λ) = 4λ(1−q)²(5+3q), and
A(n=3, p) = −256p⁴(p−1)(34p² − 108p + 91), certified at p = 18.

## The dS / OTOC axis

**E7. "Re a = −1 exactly, at every temperature," and "ED/closed-form gap
1.80–1.84."** Corrected 2026-08-25 (`0525e64`). Streicher's crossed-OTOC
display is written in the gauge θ₃₄⁺ = 0, and using it at general angles
broke time-translation invariance. The gauge-invariant coefficient is
a(v) = −1/cos(πv/2), real and negative, tending to −1 only at the
T = ∞ point. The ED gap is 1.05–1.08. The scrambling verdict survived
and strengthened.

**E8. "ED confirms the scrambling sign."** Downgraded 2026-08-06
(`e38dfcf`). ED confirms the growth. The sign is analytic, transcribed
from Streicher with an empirically fixed crossed ordering.

**E9. "The even scramblon harmonic is refuted as the mechanism that ends
the antipodal anti-scrambling epoch."** Withdrawn 2026-10-14 (`127b50f`).
The full-circle DFT behind it was dominated by samples past the 2-pt
pole at Δ = π/2 + π/v. Least-squares fits on the physical half circle
find a subdominant second harmonic (|A₂/A₁| ≈ 0.3–0.6) growing faster
than the first (κ₂/κ₁ ≈ 1.2–1.5). The mechanism is neither refuted nor
established.

**E10. "Phase and magnitude growth rates decouple at finite λ (κ > v)."**
Withdrawn 2026-10-14 (`127b50f`). The same three-point estimator on the
large-q closed form gives 0.94 on that window, so κ_eff > v is an
early-window artifact.

**E11. "The antipodal observer's 2-pt shapes fail, worse toward
semiclassics."** Withdrawn 2026-10-14 (`127b50f`). The continued form had
t and τ_s swapped. For t > 0 the antipodal 2-pt is the continued
large-q form |1/cos(𝒥z)|^(2Δ), with a deviation shrinking toward
semiclassics (0.46 → 0.18 over λ = 2 → 0.25).

**E12. "The anti-scrambling window grows 0.57 → 1.14 (t·J_c)."**
Corrected 2026-10-14. Those were grid points, one of them censored. The
interpolated maxima are 0.52 → 1.10 over λ = 1.0 → 0.3.

## The boost-to-mass calibration

**E13. "Without the state map the authors' best test cannot be run."**
Restated 2026-08-28 (`1e74201`). Only the moment route is blocked. Node
positions carry zero regulator leverage, so the test survives through
node-level data.

**E14. "Singular values floor near 1e-17 at aε = 0.1."** Corrected
2026-08-31 (`8d5fca0`). That was double-precision rounding of certified
moments. Assembled at certified precision, the conditioning law is clean
(a = 1.18, b = 0.070).

## The 't Hooft target

**E15. "FLZ's d± are eigenfunction data."** Corrected 2026-08-27
(`7f873cf`). They are spectral-determinant constants (FLZ Eq. 4.21), so
no external wavefunction anchor exists.

**E16. "The d± reproduction sits on a K = 300 double-precision floor
that worsens with K."** Superseded 2026-10-14 (`582a111`). The floor was
the direct solve. The reciprocal solve converges in K:
d₊ = 0.9631784563984(1), d₋ = 0.43358263983325(1).

## Campaigns

**E17. "The freeness crossover time is N-independent; any ln N drift is
below resolution."** Corrected 2026-08-25 (`8628159`). A drift is
resolved (t𝒥 = 1.191, 1.172, 1.233, 1.273 over N = 12–24, about 4σ).
Four points cannot separate ln N from finite size.

**E18. "The generic Dirac ensemble's particle-hole asymmetry comes from
the diagonal couplings J_AA."** Corrected 2026-08-25 (`8628159`). The
cause is normal-ordering of overlapping monomials. No number changed.

**E19. "Best fit exp(−aN) at p = 6" in the self-averaging verdicts, and a
baryon-mass ΔAIC verdict printed past its own one-decade guard.**
Retracted 2026-08-06 (`e38dfcf`) and 2026-08-25 (`8628159`). Both now
report "cannot discriminate."
