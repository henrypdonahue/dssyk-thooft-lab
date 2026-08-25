# Road map — stress-testing DSSYK∞ ↔ 't Hooft (arXiv:2607.05678)

What the paper claims, what this repo has settled, and what is left. The
literature snapshot is current as of 2026-08-24.

---

## Part I — The paper, honestly assessed

*Holograms and Standard Models* (Miyashita–Sekino–Susskind, July 2026) is
a short, talk-style manifesto rather than a technical paper. Its chain of
claims, with the paper's own confidence levels:

1. **Assumed** — their words, "unproved but we will assume it":
   DSSYK∞ = JT-de Sitter static-patch holography.
2. **Quasi-derived** (§3.5, two example diagrams; systematics deferred to
   arXiv:2501.09423): perfect self-averaging, Var⟨W⟩ ≈ p!/Nᵖ < e^(−a√N),
   hence *emergent* O(N)/U(N) symmetry of single disorder instances in
   the flat-space limit.
3. **Argued by counting** (entropy match S = N log 2 ⇒ dim matter = N):
   the bulk dual is the 't Hooft model — U(N) QCD₂ + QED, one fundamental
   flavor, at renormalized quark mass m_q² → 0, with dictionary
   N_ym = N_syk, ᾱ = ḡ²N = p², ḡ² = λ = p²/N.
4. **Qualitative only** — their words, "We don't have a quantitative
   calculation of the DSSYK spectrum": the singlet tower
   Mₙ = ψ†(dⁿ/dtⁿ)ψ ↔ 't Hooft mesons, CP = (−1)ⁿ⁺¹, with n = 0
   (photon = Q) and n = 1 (graviton = H) non-dynamical. The paper names
   the quantitative spectrum match as "the best test of the duality."
   **It has never been run, by anyone.**
5. **A logical inversion, not a computation** (§5, "Holographic Tuning"):
   at T_B = ∞ the entropy is exactly log dim H, so the vacuum energy gets
   no corrections, and the coupling ḡ²(N) becomes a *derived* quantity —
   its finite-N form is only a sketched figure.

### Two structural weaknesses the paper does not address

**The boost-vs-mass gap.** The SYK Hamiltonian is the Rindler *boost*
generator, and a fixed-mass particle has a *continuous* boost-frequency
spectrum. How discrete 't Hooft masses are supposed to appear as features
of Mₙ correlators at T_B = ∞ is an unstated piece of dictionary that
every proposed test — including the authors' — presupposes.
`duality/BOOST.md` calibrates it on the exactly solvable case and finds
that a single mass produces a continuum with mass-encoding *nodes*, not a
peak.

**Does the toy even have the disease?** QCD₂ is super-renormalizable; its
only mass renormalization, m_q² = m_bare² − g²/π, is finite. The 4d
fine-tuning drama of §5 is transplanted into a model where it barely
exists, and the vacuum-energy/entropy relation ε = 3/(8S) is only ever
written in D = 4, never derived for the JT-dS where the paper works.
`thooft-target/condensate.py` measures the 2d counterpart: it is the
single finite number ε(0) = −0.0407 Nc g².

### The regime problem

The regime where the spectrum match lives — **fixed p, λ → 0** — is *not*
the double-scaling regime where self-averaging is super-exponential. At
fixed p the suppression is only the power law N⁻ᵖ, so gauge-symmetry
emergence is quantitatively weakest exactly where the match is to be
made. This repo makes that a printed, tested table rather than a remark.

---

## Part II — What is settled

### The 't Hooft yardstick is real

An independent re-derivation of the whole Galerkin assembly, a
re-verification of the hypersingular Chebyshev-U identity by regularized
quadrature, and a check of the anchor constants against FLZ's published
tables (arXiv:0905.2280) digit-for-digit. Worst deviation over n ≤ 29 is
4.5e-12, at the level of FLZ's own last printed digit at the top of the
range; against their Padé the ground state agrees to 1.1e-15. Sum rules
7ζ(3) and 2 reproduce to ~1e-10.

Off the duality point the Chebyshev basis hardwires the α = 0 endpoint
exponent ½ while the true exponent solves πβ cot(πβ) = −α, so
convergence there is only algebraic. `jacobi_solver.py` fixes this with
a matched-exponent basis [x(1−x)]^β: ground state below 1e-11 by K = 32,
LM sum rules to ~1e-12. Eigenfunctions are validated (Parseval 3e-6,
weak-form Rayleigh ~1e-6, cross-solver overlap > 0.9999), and the decay
constants are new data — no published wavefunction tables exist.

### Self-averaging, exactly

`exact_wick.py` gives Var(B_jk) = 4σ⁴C(N−2,p−1) in closed form at any
(N, p), verified by literal contraction enumeration, ED, and Monte Carlo.
In double scaling ln rms = −¼√(λN) ln N, which beats e^(−a√N) for every
a, scanned to N = 1600 — the paper's Eq. 3.28 is not just right but
conservative. The singlet rate of Eq. 3.27 is exactly 2/C(N,p). A second
anchor, exact E[m₄] and Var(m₄), gives relVar(m₄) = 4·relVar(m₂)(1+o(1)).
The Dirac/U(N) counterpart is in `exact_wick_dirac.py`, with rms ~
Nc^(−(p−1)/2).

Worth telling the authors: the "systematic derivation" §3.5 defers to
arXiv:2501.09423 **is not there**. That paper never mentions
self-averaging; its counting rule (Eq. 4.26) is magnitude-only and
part-empirical (their §6.2), and its n ≪ p approximation fails precisely
in the diagram-B regime (error ~e^p — harmless for the bound, fatal for
prefactors). `exact_wick.py` appears to be the only complete
prefactor-exact derivation, and it *confirms* both papers at leading
order while adding the exact pinning factor 4p(N−p)/N(N−1).

### The dictionary needs four corrections

Each forced and tested (`duality/dictionary.py`,
`syk-self-averaging/dirac.py`): M₁ = −(p/2)iH, not H — it follows from
anti-Hermiticity, since Q̇ = 0 forces M₁† = −M₁, so M₁ cannot equal the
Hermitian H (verified to 1.3e-15); M₂ = −Σċ†ċ, the paper's + sign being
wrong; CP as the *unitary* particle-hole map, since M₁ ∝ iH rules out an
antiunitary one; and the λ = p²/N vs λ_chord = 2p²/N factor of two, which
is internal to the Susskind corpus (2501.09423 Eq. 3.6 uses 2p²/N).

Raw Mₙ for n ≥ 2 are CP-mixtures with 15–65% wrong-channel content, so
CP = (−1)^(n+1) only lives in CP-resolved correlators.

### The conserved-channel theorem, and where the tower actually lives

The headline result. The chord representation of every
Aₙ = Σψ(ad_H)ⁿψ is a polynomial in the chord Hamiltonian T, so **every
spectral moment of every Aₙ pair channel vanishes at N = ∞** — every λ,
every matter weight, both orderings, odd k included. The proof is a
three-term recursion whose coefficients are scalars; the two matching
conditions reduce to q-number identities independent of the chord
occupation. Explicit charges: Â₁ = −(1−q^Δ)T (the N = ∞ image of the
exact law A₁ = −2pH) and Â₂ = (1−q^Δ)²T² + (1−q^(2Δ)). Random O_Δ
operators, whose insertion q^(Δn̂) does *not* commute with T, do
propagate — that contrast is what makes it physical.

Scope, which is load-bearing: this is the i ≠ j *pair* amplitude. The
full channel also has a flavor-diagonal part of the same order — the
per-pair amplitude is smaller by exactly the 1/N that the N(N−1)
multiplicity restores — and it does **not** vanish for n ≥ 2. Evaluated
exactly at finite N, the full A₂ width rises to w₂ ≈ 1.19 along
λ_chord = 2 (1.044 → 1.194 over (p,N) = (4,16)…(10,100)), while at n = 1
the two parts cancel to 3e-16, reproducing A₁ = −2pH. So the tower does
propagate at N = ∞ and its weight is entirely diagonal. The conclusion
that the tower is a 1/N effect invisible at every order — and that the
match must therefore be made at fixed p, λ → 0 — does not follow from
this theorem and is withdrawn; the fixed-p argument survives separately
on the self-averaging side.

`exact_pair_moments.py` measures the finite-N approach to the pair zero
exactly (Wick + a CRT residue-table DP, no ED): fixed-p tail slopes
2.96, 3.90, 4.82 for n = 1, 2, 3 — i.e. μ₂ ~ N^−(n+2) — and along
λ_chord = 2, |μ₂| ≈ c/(pN) with c ≈ 29 (n = 2) and c ≈ 200–205 (n = 3).
Those are the numbers the analytic 1/N chord correction must reproduce.

Scope, stated plainly: this is the i ≠ j flavor-*pair* amplitude at
strict N = ∞. The flavor-diagonal channel and the Dirac case are open.

### The dS axis

DSSYK∞ **scrambles at every temperature**: a(v) = −1/cos(πv/2), real and
negative, tending to −1 at the hologram point, with |a| growing with β.
A dS observer needs the opposite sign, so the dictionary owes a
sign-flipping map.

The Harlow–Zhao Euclidean fold is not that map. At ED sizes the folded
connected fraction grows 0.365 → 0.786 over N = 8 → 20 while λ shrinks,
and at N = ∞ it stays O(1) and grows as λ → 0 while the unfolded one dies
like λ. Both statements are closed-form-free. The mechanism behind the
bad finite-N convergence is measured: folded correlators at ED sizes
track the sample's upper spectral edge at correlation −0.86.

One candidate map does work, semiclassically: the antipodal
half-fake-period shift Δ = π/v flips the sign exactly in the closed form
(a(Δ) = a(0)e^(−ivΔ), so the flip is a clean negation) and for λ ≤ 1.7 in
the exact theory, *and* makes the Wightman sector exactly thermal at
β_eff = β + β_fake. One map, both of Cui–Kolchmeyer's dS-observer
signatures. It is NV's ±iβ/4 dressing, geometrically. The deep campaign
finds a transient anti-scrambling epoch growing like ~½ln(1/λ) rather
than eternal dS, locates the chaos-bound evasion (the phase winds at
constant magnitude, so there is no strip-KMS structure), and records the
costs honestly: the shifted contour is fold-adjacent and edge-dominated,
and rate extraction is mode-contaminated.

### Other results

`freeness.py` exhibits the free-product transition — spec(E + F(t))
landing on the arcsine law, W₁ plateau falling like dim^(−0.9) — the
first time in a microscopic dS-candidate model. `qed_campaign.py` closes
the "Majorana numerics vs Dirac claims" gap. `baryons.py` gives first
data on the wholly-underived inter-baryon-force claim. `condensate.py`
transcribes the exact condensate and gives first α → −1 GMOR data.

---

## Part III — What is open, and what each would take

### The tower extraction (the headline; a live race)

Now located at fixed p, λ → 0, as a 1/N computation from the
double-scaled viewpoint. The exact finite-N pair moments are its first
bite and supply the target laws. Routes:

- The **analytic 1/N chord correction** that reproduces μ₂ ~ N^−(n+2) at
  fixed p and ≈ c/(pN) on the double-scaling line. This is the theory
  half, and the most valuable single result available.
- Extending the exact computation to n ≥ 4 or μ₆ needs K ≥ 12 H-node
  words, whose matching count grows as (K−1)!!. The DP handles any number
  of subsets; the matching enumeration is the wall.
- The **flavor-diagonal channel** at N = ∞. The theorem covers only i ≠ j
  pairs, while the ED bench says the diagonal part dominates at
  accessible sizes. A chord treatment would close the gap between the
  theorem and the bench.
- **Moments → spectrum** is an inverse problem, and with μ_k → 0 at
  N = ∞ the question of which finite-N object carries the tower, and how
  many moments a Regge trajectory needs, is unanswered.
- The **Dirac/complex-SYK chord machinery**, which the CP-resolved data
  is waiting on.

Comparison to `thooft-target` then goes through the BOOST.md kernel — a
mass is a Bessel-kernel continuum with nodes, not a peak — and
`dictionary.py`.

### The O(1/q) residues

The large-q singlet 4-point function has Sommerfeld–Watson poles only at
m = 0 (Hamiltonian exchange, the entire time-ordered part) and m = ±v
(the scramblon). The would-be tower cancels identically at leading order
and its couplings are O(1/q²). The tower first appears at O(1/q),
organized on the exact-in-temperature skeleton 𝓜ᵉ/𝓜ᵒ, whose
finite-temperature residues are unpublished. Computing them is the open
task, and the conceptual question of how the 1/q expansion at fixed λ
relates to the 1/N one deserves a clean answer.

### The tomperature map

The antipodal shift passes both necessary conditions semiclassically, but
a sign flip plus effective thermality under a contour shift is not a
dictionary. What would upgrade it is the algebra story — why the
observer's operators live at the antipodal point — which is exactly NV's
construction to complete.

### The NV discriminating computation

The OTOC sign of Narovlansky–Verlinde's dressed operators. They engineer
KMS at the dS temperature, which is the no-go's hypothesis, so they must
break another hypothesis of the bound. This is the sharpest next
computation on the dS axis and needs the two-sided constrained chord
sector — adjacent to the existing engine. Since the antipodal shift *is*
their dressing geometrically, the answer may already be predictable from
`signflip_probe.py`, which would itself be worth establishing.

### The bulk side

- **JT-dS ε = 3/(8S)**: the paper writes it only in D = 4. Deriving the
  JT-dS analog is the still-open theory half of the vacuum-energy story.
- **The finite-N coupling curve ḡ²(N)**: Fig. 17 is a sketch with no
  formula, yet the whole set-it-and-forget-it inversion rests on it.
  Match a bulk observable at a sequence of finite N and extrapolate —
  gated on the tower extraction.
- **The order-of-limits problem**: T_B → ∞, N → ∞, λ fixed vs → 0, and
  the boost-spectrum → invariant-mass reorganization. The
  chaos-to-integrability transition at small λ (2403.01950, 2601.09801)
  is the likely mechanism zone, and this is where the duality most
  plausibly *fails*. A sharp statement of where it holds would be a real
  contribution.
- **Freeness sharpening**: the ln N onset-time question needs sizes
  beyond ED.

### Out of scope

- Proving DSSYK∞ = JT-dS — assumed by the authors themselves.
- The 1+1 → 3+1 uplift. Contrary to how it is sometimes described, *no*
  uplift construction exists in the source; §§5–6 are explicit that this
  is a hope.
- Instance-level ED in the true double-scaling regime (p ~ √N, N ≳ 100,
  dim 2⁵⁰) — the honest wall. The exact-Wick and chord computations are
  the ways around it.

---

## Context

The QCD₂ side is effectively solved beyond FLZ: Ambrosino–Komatsu recast
't Hooft's equation as a TQ-Baxter integrable system (2312.15598,
2406.11078), and Litvinov–Meshcheriakov (2409.11324) give independent
high-precision spectra. Numerical effort should go almost entirely to the
DSSYK side.

The DSSYK toolkit is mature: chord techniques (1811.02584; review
Berkooz–Mamroud 2407.09396; path integral 2403.05980), and the
chaos-to-integrability transition at small λ is directly reusable for the
λ → 0 regime the duality needs.

The duality sits in contested territory. Narovlansky–Verlinde
(2310.16994) and the sine-dilaton program (Blommaert et al. 2404.03535,
2505.08116) propose *different* bulk duals, and Rahman–Susskind 2407.12988
argues the identifications conflict. A clean meson-tower match or mismatch
would be the first quantitative discriminator among all three programs.
`duality/DISCRIMINATOR.md` lays out that axis.

**Convention landmine:** this paper uses λ = p²/N with no q; the chord
literature uses λ = 2p²/N, q = e^(−λ). Any reuse of chord results must
pass through a tested conversion.

**Priority.** As of 2026-08-24 the target paper has exactly one citation
— Cui–Kolchmeyer 2607.13665, which argues dS observer OTOCs
*anti-scramble* and poses a structural challenge to observer-centric
static-patch holography. Its same-week companion Harlow–Zhao 2607.14215
proposes the Euclidean-fold mechanism and poses the SYK test that
`duality/ANTISCRAMBLING.md` answers. The authors' promised meson-spectrum
computation has not appeared, so the first-mover race on the headline
test is still open.
