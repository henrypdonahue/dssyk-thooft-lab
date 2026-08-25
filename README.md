# Stress-testing the DSSYK∞ ↔ 't Hooft duality

Miyashita, Sekino, and Susskind propose that DSSYK at infinite temperature
is dual to the 't Hooft model — large-N QCD in two dimensions
([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)). This repo tests
that proposal with numbers. No new theory. Every claim is cross-checked
against an independent method, and a test asserts every verdict.

Assessment of the paper, the open problems, and why they are hard:
[road_map.md](road_map.md).

## Modules

| module | role | tests |
|---|---|---|
| [thooft-target/](thooft-target) | the yardstick: 't Hooft meson spectrum (12+ digits), condensate, vacuum energy | 45 |
| [syk-self-averaging/](syk-self-averaging) | the SYK-side claims: self-averaging, the Mₙ/CP sector, QED and baryon campaigns, freeness, singlet moments | 179 |
| [duality/](duality) | the dictionary, the exact N = ∞ chord engine, the conserved-channel theorem, the OTOC/fold laboratory, the KMS discriminator | 195 |

Install once: `pip install -e ".[dev]"`. Then `make test` runs the 408
fast assertions in about 40 s; `make test-all` adds the slow ED sweeps,
419 in total.

## Main results

**The two-flavor part of the singlet tower is conserved at N = ∞ —
proven.** For every Aₙ = Σψ(ad_H)ⁿψ, the chord representation of the
i ≠ j *pair* amplitude is a polynomial in the chord Hamiltonian, so all
its spectral moments vanish at every λ and every matter weight. A
three-term recursion proves it; A₁ = −2pH is the n = 1 case. Measured to
n = 5 and μ₆ across the λ and matter-weight grids, and the finite-N
approach is computable exactly (Wick, no ED): |μ₂| ~ N^−(n+2) at fixed p,
and ≈ c/(pN) with c ≈ 29 (n = 2) along λ_chord = 2.
**Scope, and it matters:** the full channel also carries a
flavor-diagonal part of the same order, which does *not* vanish for
n ≥ 2 — the full A₂ width rises to w₂ ≈ 1.19 along λ_chord = 2. So the
tower propagates at N = ∞, with all of its weight diagonal. What is
established is a flavor statement, not an invisibility statement.
→ [duality/CHORD.md](duality/CHORD.md)

**Self-averaging holds exactly**, and more sharply than the paper claims:
ln rms = −¼ √(λN) ln N in double scaling. The emergent-symmetry counting
checks out. At fixed p, where the spectrum match must be made, the
suppression is only polynomial.
→ `syk-self-averaging/exact_wick.py`

**The operator dictionary needs four corrections.** Each is forced and
tested: M₁ = −(p/2)iH, the sign of M₂, CP as the unitary particle-hole
map, and the λ-convention factor of two.
→ [duality/](duality)

**A de Sitter dictionary must meet a measured spec, and one candidate map
passes it semiclassically.** The map must flip the OTOC sign, because
DSSYK∞ scrambles: the growth coefficient is a(v) = −1/cos(πv/2) < 0 at
every temperature. The Euclidean fold fails to supply that flip, even at
N = ∞. The antipodal half-fake-period shift supplies it, and makes the
observer thermal at the tomperature — exactly in the closed form, and in
the exact theory for λ ≲ 1.7.
→ [duality/ANTISCRAMBLING.md](duality/ANTISCRAMBLING.md),
[duality/DISCRIMINATOR.md](duality/DISCRIMINATOR.md)

**The 2d vacuum energy is one finite number**: ε(0) = −0.0407 Nc g².
This bulk has nothing to tune. First data in the chiral corner too: the
GMOR slope with its next-order correction measured.
→ `thooft-target/condensate.py`

**The 't Hooft yardstick is ready**: spectrum to 12+ digits, decay
constants (new data), spectral accuracy at nonzero quark mass.
→ [thooft-target/](thooft-target)

## Open

Four problems stay open: extracting the tower itself (a 1/N computation,
and a live race with the authors), the O(1/q) residues, the tomperature
map, and the JT-dS entropy relation. Two are out of scope: proving
DSSYK∞ = JT-dS, which the authors assume, and the 1+1 → 3+1 uplift.
[road_map.md](road_map.md) says what each would take.
