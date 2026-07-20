# What a mass looks like in boost frequency (rung 13 calibration)

The DSSYK Hamiltonian generates **boosts**, and every spectrum-match plan implicitly
assumes "peaks in $M_n$ correlators = meson masses." This note calibrates that
assumption on the exactly-solvable case — one free particle of mass $m$ in 1+1d,
probed along a uniformly accelerated worldline (acceleration $a$) — and it fails in an
instructive way.

**Exact result** (`boost_mass.py`, derived from the plane-wave expansion and verified
against an independent Fourier transform of the $K_0$ Wightman function to ~1e-12):

$$S(\omega) \;=\; \frac{1}{\pi a}\, e^{\pi\omega/a}\, K_{i\omega/a}(m/a)^2 .$$

## The four calibration facts

1. **One discrete mass ⇒ a continuum.** $S(\omega)$ is supported on all $\omega$ — no
   pole, no peak at $\omega = m$. So broad continua in `mn_spectroscopy.py` do **not**
   falsify a discrete meson tower, and peak-hunting at fixed setup is the wrong
   observable.
2. **The mass lives in the node pattern.** $K_{i\omega/a}(m/a)$ oscillates in $\omega$
   with zeros; $S$ has near-nodes whose positions shift with $m$ (first node
   $\omega\approx3.0$ at $m/a=1$ vs $\approx4.4$ at $m/a=2$). Masses are extractable
   from the node/oscillation structure or moments via the $K$-kernel — a linear
   inverse problem, $S_{\rm tot}(\omega)=\sum_n c_n\,S(\omega; M_n, a)$ for a tower.
3. **Heavy = faint.** $S(0)\sim e^{-2m/a}$: at fixed $a$, higher tower members are
   exponentially Boltzmann-suppressed. Any finite-precision spectroscopy sees only the
   first few masses.
4. **A sharp open tension.** The bulk response is KMS-asymmetric at the Unruh
   temperature ($S(\omega)=e^{2\pi\omega/a}S(-\omega)$, exact above). At
   $T_B=\infty$ the hologram obeys $S_{AA^\dagger}(\omega)=S_{A^\dagger A}(-\omega)$
   — the KMS factor degenerates to 1, and the symmetrized correlator (the
   convention-stable object used in `mn_spectroscopy`/`moments_pipeline`) is exactly
   $\omega$-symmetric. An asymmetric Unruh response can only emerge through a
   nontrivial state/frequency map — the "tomperature" identification
   (Rahman–Susskind 2401.08555). Until that map is fixed, comparing DSSYK spectral
   shapes with bulk masses is convention-ambiguous. **This is the gating question of
   rung 13/18**, now stated precisely.

## Implication for the headline test (rung 18)

The rung-18 deliverable "discrete tower or continuum?" must be reformulated:
the boost-frequency data will look continuous either way. The discriminating
observables are (a) the node/oscillation structure and its $m/a$ scaling, (b)
moment ratios against the $K$-kernel-transformed 't Hooft tower, and (c) the
KMS/symmetry structure once the state map is fixed. `boost_mass.py` provides the
kernel; `duality/dictionary.py` the tower.

*Scope caveats:* free scalar, 1+1d, single particle, sharp worldline. Composite
(bilinear) operators smear the kernel; interactions dress it. The calibration
bounds what can be concluded from shape data — it does not yet fix the
tomperature map.
