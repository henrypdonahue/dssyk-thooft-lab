# What a mass looks like in boost frequency

The DSSYK Hamiltonian generates **boosts**, not translations.

Every spectrum-match plan quietly assumes something: that peaks in Mₙ
correlators are meson masses. This note calibrates that assumption on the
one case you can solve exactly. The assumption fails, and it fails in a
useful way.

The exactly solvable case is one free particle of mass m in 1+1
dimensions, probed along a uniformly accelerated worldline of
acceleration a. `boost_mass.py` gives

    S(ω) = (1/πa) · exp(πω/a) · K_{iω/a}(m/a)²

verified against an independent Fourier transform of the Wightman
function to about 1e-12.

## Four calibration facts

1. **One discrete mass gives a continuum.** S(ω) has support at every ω.
   There is no pole and no peak at ω = m. So broad continua in
   `mn_spectroscopy.py` do not falsify a discrete meson tower.
   Peak-hunting is the wrong observable.
2. **The mass lives in the node pattern.** The Bessel kernel oscillates
   in ω, and the zeros move with m. The first node sits at ω ≈ 3.0 for
   m/a = 1, and at ω ≈ 4.4 for m/a = 2. You can extract masses from the
   node structure, or from moments through the kernel. For a tower that
   is a linear inverse problem.
3. **Heavy is faint.** S(0) ~ exp(−2m/a). Higher tower members are
   exponentially suppressed. Finite-precision spectroscopy sees only the
   first few masses.
4. **A sharp open tension.** The bulk response is KMS-asymmetric at the
   Unruh temperature. The T = ∞ hologram's symmetrized correlator is
   exactly ω-symmetric. An asymmetric response can only appear through a
   nontrivial state and frequency map. That map is the "tomperature"
   identification (Rahman–Susskind [2401.08555](https://arxiv.org/abs/2401.08555)). Until it is fixed,
   comparing DSSYK spectral shapes to bulk masses is
   convention-ambiguous. This is the gating question of the headline
   test, stated precisely.

## What this means for the headline test

Read "Is the spectrum a tower or a continuum?" in
[CHORD.md](CHORD.md) with fact 1 in hand. Boost-frequency data looks
continuous either way.

The observables that can tell the difference are:

- the node structure and how it scales with m/a,
- moment ratios against the kernel-transformed 't Hooft tower,
- the KMS structure, once the state map is fixed.

`boost_mass.py` supplies the kernel. `dictionary.py` supplies the tower.

## The forward check has now been run

`boost_forward.py` pushes the 't Hooft masses and decay constants forward
through the boost kernel. The result is structural.

Raw moments of the boost kernel diverge, because S(ω) → 1/ω. So moment
invariants are regulator-dominated: kurt(ε) sweeps 1.8 to 25 across the
whole SYK band as the regulator turns. And at the invariant level the
coupling-weighted tower differs from a single ground-state mass by no
more than 5%.

So moments cannot see the tower. The spec the state map must supply is
the ε* consistency table in `boost_forward.json`: one smearing scale has
to fit every channel, not one scale per channel. The tower discriminant
is the nodes.

*Scope: free scalar, 1+1 dimensions, single particle, sharp worldline.
Composite operators smear the kernel. Interactions dress it. This
calibration bounds what shape data can show. It does not fix the
tomperature map.*
