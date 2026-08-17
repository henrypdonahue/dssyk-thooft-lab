# Stress-testing the DSSYK∞ ↔ 't Hooft duality

Miyashita, Sekino, and Susskind propose that DSSYK at infinite temperature
is dual to the 't Hooft model — large-N QCD in two dimensions
([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)). This repo tests
that proposal with numbers. No new theory. Every claim is cross-checked
against an independent method, and every verdict is asserted by a test.

Full assessment, task ladder, and progress log: [road_map.md](road_map.md).
The paper: [2607.05678v1.pdf](2607.05678v1.pdf).

## Modules

| module | role | tests |
|---|---|---|
| [thooft-target/](thooft-target) | the yardstick: 't Hooft meson spectrum (12+ digits), condensate, vacuum energy | 42 |
| [syk-self-averaging/](syk-self-averaging) | the SYK-side claims: self-averaging, the Mn/CP sector, QED and baryon campaigns, freeness, singlet moments | 103 |
| [duality/](duality) | the dictionary and fingerprints, the exact N = ∞ chord engine, the OTOC/fold laboratory, the KMS discriminator | 151 |

Install once: `pip install -e ".[dev]"`. Then `make test` (about 3 minutes)
reruns every assertion.

## Main results

- **The meson tower is invisible at strict double scaling.** The singlet
  channels carry zero propagating weight at N = ∞: exact for n = 1
  (A₁ = −2pH), machine-precision cancellation for n = 2, 3, at every λ and
  matter weight tested. The tower is a 1/N effect. It lives at fixed p,
  λ → 0 — exactly where symmetry emergence is weakest. The paper's "best
  test" is now located. → [duality/CHORD.md](duality/CHORD.md)
- **Self-averaging holds exactly**, sharper than the paper claims
  (ln rms = −¼ √N ln N). The emergent-symmetry counting checks out.
  → `syk-self-averaging/exact_wick.py`
- **The operator dictionary needs four corrections**, each forced and
  tested: M₁ = −(p/2)iH, the sign of M₂, CP as the unitary particle-hole
  map, and the λ-convention factor of two. → [duality/](duality)
- **A de Sitter dictionary must meet a measured spec.** It must flip the
  OTOC sign (Re a = −1 at every temperature). The Euclidean fold cannot do
  it — it fails even at N = ∞. Engineered dS-KMS meets the chaos bound.
  → [duality/ANTISCRAMBLING.md](duality/ANTISCRAMBLING.md),
  [duality/DISCRIMINATOR.md](duality/DISCRIMINATOR.md)
- **The 2d vacuum energy is one finite number**: ε(0) = −0.0407 Nc g².
  This bulk has nothing to tune. Plus first data in the chiral corner
  (GMOR slope, with its next-order correction measured).
  → `thooft-target/condensate.py`
- **The 't Hooft yardstick is ready**: spectrum to 12+ digits, decay
  constants (new data), spectral accuracy at nonzero quark mass.
  → [thooft-target/](thooft-target)

## Open

The tower extraction itself (a 1/N computation — the live race with the
authors), the O(1/q) residues, the tomperature map, and the JT-dS entropy
relation. Out of scope: proving DSSYK∞ = JT-dS (the authors assume it) and
the 1+1 → 3+1 uplift. Details in [road_map.md](road_map.md).
