# Stress-testing the DSSYK∞ ↔ 't Hooft duality

Start with the picture the proposal is built on.

A JT static patch has a Schwarzian boundary somewhere between you and the
horizon. Inside it there is ordinary AdS bulk. Between it and the horizon
there is a thin layer. Susskind calls that layer the "sliver."

Now heat the hologram. The Schwarzian surface moves toward the horizon.
The sliver grows. The bulk shrinks. At infinite Boltzmann temperature the
sliver fills the whole static patch, and it has de Sitter geometry.

Is there life in the sliver?

Miyashita, Sekino and Susskind answer yes
([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)). The life is
't Hooft mesons. Their claim: DSSYK at infinite temperature is dual to
the 't Hooft model, which is large-N QCD in two dimensions. They name one
test as the best one. Match the two mass spectra.

Nobody has run that test. This repo builds the parts you need to run it,
and measures the pieces that are already reachable.

No new theory here. Every number is cross-checked against an independent
method. Every verdict has a test that asserts it. The assessment of the
paper and of the open problems is in [road_map.md](road_map.md).

## Modules

| module | role | tests |
|---|---|---|
| [thooft-target/](thooft-target) | the yardstick: the 't Hooft meson spectrum to 12+ digits, the condensate, the vacuum energy | 45 |
| [syk-self-averaging/](syk-self-averaging) | the SYK-side claims: self-averaging, the Mₙ/CP sector, QED and baryon campaigns, freeness, the exact 1/N laws | 189 |
| [duality/](duality) | the dictionary, the N = ∞ chord engines, the conserved-channel theorem, the OTOC laboratory, the KMS discriminator | 236 |

Install once with `pip install -e ".[dev]"`. Then `make test` runs the 458
fast tests in about 90 seconds. `make test-all` adds the slow ED sweeps,
470 in total.

## Main results

### The pair channel is frozen at N = ∞

This one is proven, not measured.

The paper's singlet channels are Aₙ = Σψ(ad_H)ⁿψ. You build them out of
the Hamiltonian. You hit a fermion with H, n times.

At N = ∞ the chord picture says what such a probe is. It is a polynomial
in the chord Hamiltonian, and nothing else. A polynomial in the
Hamiltonian is conserved. It does not tick. All of its weight sits at
zero frequency.

A three-term recursion proves it, for every n, every λ and every matter
weight. The n = 1 case is the exact operator identity A₁ = −2pH. The
engine measures the same zeros to n = 5 and μ₆.

The finite-N approach to that zero is now an exact law, not a fit
(`syk-self-averaging/pair_moment_laws.py`). μ₂ ~ N^−(n+2) holds exactly
through n = 4. The coefficients are certified: A(n=1) = −16p²(p−1) and
A(n=2) = −64p³(p−1)(4p−5). On the double-scaling line the constant runs
c₂ = 28.92 → 30.53 over p = 4..16.

**The scope matters, and it is the whole story.** ⟨AₙAₙ⟩ has two parts.
The pair channel (i ≠ j) is the coherent one. The diagonal channel
(i = j) is an incoherent pile of N single-flavor pieces. They are the
same order in N. The theorem kills the first. The second survives, and
for n ≥ 2 it grows: the full A₂ width rises to w₂ ≈ 1.19 along
λ_chord = 2. The diagonal channel now has its own N = ∞ chord law
(`duality/diagonal_chord.py`), w₂ = 4.67(9) at λ_chord = 2.

So the tower does propagate at N = ∞. All of its weight is diagonal.
This is a statement about flavor, not a statement that the tower is
invisible.
→ [duality/CHORD.md](duality/CHORD.md)

### Self-averaging holds exactly

The paper argues from two diagrams. `exact_wick.py` gives the closed form
at any (N, p), with the prefactor.

In double scaling, ln rms = −¼ √(λN) ln N. That beats the paper's own
bound. The emergent-symmetry counting checks out.

At fixed p the suppression is only polynomial. That is where the spectrum
match has to be made, so it is worth knowing.
→ `syk-self-averaging/exact_wick.py`

### The dictionary needs four corrections

Each one is forced, and each one has a test. M₁ = −(p/2)iH. The sign of
M₂. CP as the unitary particle-hole map. The λ-convention factor of two.

A silent factor of two would wreck any match. The dictionary was built
before the DSSYK side existed, so there is exactly one place for it to
hide.
→ [duality/](duality)

### A dS dictionary must flip the OTOC sign

DSSYK∞ scrambles. The growth coefficient is a(v) = −1/cos(πv/2), which is
negative at every temperature. A de Sitter observer needs the opposite
sign. So the dictionary owes a map that flips it.

The Harlow–Zhao Euclidean fold does not supply that flip. It fails at ED
sizes and it fails at N = ∞.

One map does work, semiclassically. Shift the probe pair by half the fake
period. The flip is exact in the closed form, and it holds in the exact
theory for λ ≲ 1.7. The same shift makes the observer thermal at the
tomperature.

That shift is not a new invention. It is NV's ±iβ_dS/4 dressing, at the
contour level, equation-cited and tested
(`duality/NV-DRESSING.md`). Conditional on one constrained-sector step,
their discriminating OTOC is already measured here.
→ [duality/ANTISCRAMBLING.md](duality/ANTISCRAMBLING.md),
[duality/DISCRIMINATOR.md](duality/DISCRIMINATOR.md)

### The 2d vacuum energy is one number

ε(0) = −0.0407 Nc g². This bulk has nothing to tune. There is first data
in the chiral corner too: the GMOR slope with its next-order correction.
→ `thooft-target/condensate.py`

### The yardstick is ready

The spectrum runs to 12+ digits. The decay constants are new data. The
solver stays accurate at nonzero quark mass.
→ [thooft-target/](thooft-target)

## What is still open

**Extract the tower.** The exact 1/N laws and the diagonal chord law are
in hand. The chord-side derivation is not. Neither is the step from
moments to a spectrum, and the boost-kernel forward check says that step
needs node-level data, or the state map fixed first. This is a live race
with the authors.

**The O(1/q) residues.** The tower cancels at leading order in large q.
It first appears at O(1/q), on a skeleton whose finite-temperature
residues are unpublished.

**The tomperature map.** Now reduced to one constrained-sector step in
NV's construction.

**The JT-dS entropy relation.** The paper writes ε = 3/(8S) only in D = 4.

Two problems are out of scope. Proving DSSYK∞ = JT-dS, which the authors
assume. And the 1+1 → 3+1 uplift, which no source constructs.

[road_map.md](road_map.md) says what each one would take.
