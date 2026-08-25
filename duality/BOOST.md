# What a mass looks like in boost frequency

The DSSYK Hamiltonian generates **boosts**. Every spectrum-match plan
implicitly assumes "peaks in Mn correlators = meson masses." This note
calibrates that assumption on the exactly solvable case. That case is one
free particle of mass m in 1+1d, probed along a uniformly accelerated
worldline (acceleration a). The assumption fails in an instructive way.

Exact result (`boost_mass.py`, verified against an independent Fourier
transform of the Wightman function to ~1e-12):

    S(ω) = (1/πa) · exp(πω/a) · K_{iω/a}(m/a)²

## The four calibration facts

1. **One discrete mass gives a continuum.** S(ω) has support at all ω.
   It has no pole and no peak at ω = m. Broad continua in
   `mn_spectroscopy.py` do **not** falsify a discrete meson tower.
   Peak-hunting is the wrong observable.
2. **The mass lives in the node pattern.** The Bessel kernel oscillates
   in ω with zeros whose positions shift with m. The first node is at
   ω ≈ 3.0 for m/a = 1, and at ω ≈ 4.4 for m/a = 2. Masses are
   extractable from the node structure or from moments via the kernel.
   This is a linear inverse problem for a tower.
3. **Heavy is faint.** S(0) ~ exp(−2m/a). Higher tower members are
   exponentially suppressed. Finite-precision spectroscopy sees only the
   first few masses.
4. **A sharp open tension.** The bulk response is KMS-asymmetric at the
   Unruh temperature. The T = ∞ hologram's symmetrized correlator is
   exactly ω-symmetric. An asymmetric response can only emerge through a
   nontrivial state/frequency map. That map is the "tomperature"
   identification (Rahman–Susskind 2401.08555). Until that map is fixed,
   comparing DSSYK spectral shapes to bulk masses is convention-ambiguous.
   This is the gating question of the headline test, stated precisely.

## Implication for the headline test

Read "Discrete tower or continuum?" carefully. Boost-frequency data looks
continuous either way. The discriminating observables are (a) the node
structure and its m/a scaling, (b) moment ratios against the
kernel-transformed 't Hooft tower, and (c) the KMS structure once the
state map is fixed. `boost_mass.py` provides the kernel.
`dictionary.py` provides the tower.

*Scope: free scalar, 1+1d, single particle, sharp worldline. Composite
operators smear the kernel. Interactions dress it. This calibration
bounds what shape data can show. It does not fix the tomperature map.*
