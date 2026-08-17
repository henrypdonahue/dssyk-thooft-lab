# What a mass looks like in boost frequency

The DSSYK Hamiltonian generates **boosts**, and every spectrum-match plan
implicitly assumes "peaks in Mn correlators = meson masses." This note
calibrates that assumption on the exactly solvable case — one free
particle of mass m in 1+1d, probed along a uniformly accelerated worldline
(acceleration a) — and the assumption fails in an instructive way.

Exact result (`boost_mass.py`, verified against an independent Fourier
transform of the Wightman function to ~1e-12):

    S(ω) = (1/πa) · exp(πω/a) · K_{iω/a}(m/a)²

## The four calibration facts

1. **One discrete mass gives a continuum.** S(ω) has support at all ω —
   no pole, no peak at ω = m. Broad continua in `mn_spectroscopy.py` do
   **not** falsify a discrete meson tower. Peak-hunting is the wrong
   observable.
2. **The mass lives in the node pattern.** The Bessel kernel oscillates
   in ω with zeros whose positions shift with m (first node ω ≈ 3.0 at
   m/a = 1 vs ≈ 4.4 at m/a = 2). Masses are extractable from the node
   structure or from moments via the kernel — a linear inverse problem
   for a tower.
3. **Heavy is faint.** S(0) ~ exp(−2m/a): higher tower members are
   exponentially suppressed. Finite-precision spectroscopy sees only the
   first few masses.
4. **A sharp open tension.** The bulk response is KMS-asymmetric at the
   Unruh temperature. The T = ∞ hologram's symmetrized correlator is
   exactly ω-symmetric. An asymmetric response can only emerge through a
   nontrivial state/frequency map — the "tomperature" identification
   (Rahman–Susskind 2401.08555). Until that map is fixed, comparing DSSYK
   spectral shapes to bulk masses is convention-ambiguous. This is the
   gating question of the headline test, stated precisely.

## Implication for the headline test

"Discrete tower or continuum?" must be read carefully: boost-frequency
data looks continuous either way. The discriminating observables are (a)
the node structure and its m/a scaling, (b) moment ratios against the
kernel-transformed 't Hooft tower, and (c) the KMS structure once the
state map is fixed. `boost_mass.py` provides the kernel;
`dictionary.py` the tower.

*Scope: free scalar, 1+1d, single particle, sharp worldline. Composite
operators smear the kernel; interactions dress it. This calibration
bounds what shape data can show — it does not fix the tomperature map.*
