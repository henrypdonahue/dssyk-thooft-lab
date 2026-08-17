# Stress-testing the DSSYK∞ ↔ 't Hooft duality

Computational legwork for Miyashita–Sekino–Susskind, *"Holograms and Standard
Models"* ([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)): DSSYK at
infinite temperature, claimed dual to the 't Hooft model (large-$N$ 2d QCD).
No new theory — target numbers and stress-tests, each cross-checked against an
independent method or analytic anchor, every verdict asserted by a test.

Full assessment, task ladder, and progress log: [`road_map.md`](road_map.md).
The paper: [`2607.05678v1.pdf`](2607.05678v1.pdf).

## Modules

| module | role | tests |
|---|---|---|
| [`thooft-target/`](thooft-target) | the yardstick: 't Hooft meson spectrum (FLZ-anchored to 2e-12), condensate/GMOR corner, vacuum energy | 42 |
| [`syk-self-averaging/`](syk-self-averaging) | the falsifiable SYK-side claims: self-averaging, $M_n$/CP sector, QED + baryon campaigns, freeness onset, singlet-moment bench | 103 |
| [`duality/`](duality) | the dictionary + fingerprints, boost-to-mass calibration, large-q anchor, **the exact $N=\infty$ chord engine**, OTOC/fold laboratory, KMS-axis discriminator | 151 |

`pip install -e ".[dev]"`, then `make test` (~3 min) reruns every assertion.

## Headline results (one line each; details behind the links)

- **The $M_n$ tower is invisible at strict double scaling** — the singlet
  channels carry zero propagating weight at $N=\infty$ ($A_1=-2pH$ exact;
  $n=2,3$ a machine-precision cancellation across the λ- and weight-scans,
  tested); the tower is a $1/N$ effect living at fixed $p$, λ→0, where
  symmetry emergence is only polynomial. The paper's "best test" is now
  precisely located. → [`duality/CHORD.md`](duality/CHORD.md)
- **Self-averaging verified exactly**, sharper than claimed
  ($\ln\mathrm{rms}=-\tfrac14\sqrt N\ln N$); the emergent-symmetry counting
  confirmed at leading order. → `syk-self-averaging/exact_wick.py`
- **Dictionary corrections** (each forced and tested): $M_1=-(p/2)\,iH$, the
  $M_2$ sign, CP = the unitary particle-hole map,
  $n_{\rm paper}=n_{\rm FLZ}+2$, the $\lambda$ factor-2 landmine.
  → [`duality/`](duality)
- **What a de Sitter dictionary must accomplish** (constraints, not
  obituaries): flip the OTOC sign ($\mathrm{Re}\,a=-1$ at every $T$); not via
  the Euclidean fold (a λ-obstruction at $N=\infty$; at ED sizes the fold's
  slow $1/N$ convergence is measured to be spectral-edge domination); evade
  the chaos bound that engineered dS-KMS walks into.
  → [`duality/ANTISCRAMBLING.md`](duality/ANTISCRAMBLING.md),
  [`duality/DISCRIMINATOR.md`](duality/DISCRIMINATOR.md)
- **The 2d vacuum energy is one finite number**,
  $\varepsilon(0)=-0.0407\,N_cg^2$ — nothing to tune in this bulk; plus first
  chiral-corner solver data (GMOR slope + measured NLO).
  → `thooft-target/condensate.py`
- **'t Hooft yardstick ready**: spectrum to 12+ digits, decay constants (new
  data), α≠0 spectrally accurate. → [`thooft-target/`](thooft-target)

## Status

Open, and why (details in `road_map.md`): the tower extraction itself (now a
subleading-in-$1/N$ computation — the live race with the authors), the
$O(1/q)$ residues, the tomperature map (sign-gated), the JT-dS
$\varepsilon=3/(8S)$ derivation. Out of scope: proving DSSYK∞ = JT-dS
(assumed by the authors), the $1{+}1\to3{+}1$ uplift.
