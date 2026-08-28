# Road map — stress-testing DSSYK∞ ↔ 't Hooft ([arXiv:2607.05678](https://arxiv.org/abs/2607.05678))

What the paper claims. What this repo has settled. What is left. The
literature snapshot is current as of 2026-08-24.

---

## Part I — The paper, honestly assessed

*Holograms and Standard Models* (Miyashita–Sekino–Susskind, July 2026) is
a short talk-style manifesto. It is not a technical paper, and it does
not pretend to be. Here is its chain of claims, with the paper's own
confidence levels attached.

1. **Assumed.** Their words: "unproved but we will assume it." DSSYK∞ =
   JT-de Sitter static-patch holography.
2. **Quasi-derived** (§3.5, from two example diagrams, with the
   systematics deferred to [arXiv:2501.09423](https://arxiv.org/abs/2501.09423)). Perfect self-averaging,
   Var⟨W⟩ ≈ p!/Nᵖ < e^(−a√N). Hence *emergent* O(N)/U(N) symmetry of
   single disorder instances in the flat-space limit.
3. **Argued by counting.** The entropy match S = N log 2 gives
   dim matter = N. So the bulk dual is the 't Hooft model: U(N) QCD₂ plus
   QED, one fundamental flavor, at renormalized quark mass m_q² → 0, with
   the dictionary N_ym = N_syk, ᾱ = ḡ²N = p², ḡ² = λ = p²/N.
4. **Qualitative only.** Their words: "We don't have a quantitative
   calculation of the DSSYK spectrum." The singlet tower
   Mₙ = ψ†(dⁿ/dtⁿ)ψ maps to 't Hooft mesons, with CP = (−1)ⁿ⁺¹, and with
   n = 0 (photon = Q) and n = 1 (graviton = H) non-dynamical. The paper
   names the quantitative spectrum match as "the best test of the
   duality." **Nobody has run it.**
5. **A logical inversion, not a computation** (§5, "Holographic
   Tuning"). At T_B = ∞ the entropy is exactly log dim H, so the vacuum
   energy gets no corrections, and the coupling ḡ²(N) becomes a *derived*
   quantity. Its finite-N form is only a sketched figure.

### Two structural weaknesses the paper does not address

**The boost-against-mass gap.** The SYK Hamiltonian is the Rindler
*boost* generator. A fixed-mass particle has a *continuous*
boost-frequency spectrum. How discrete 't Hooft masses are supposed to
appear as features of Mₙ correlators at T_B = ∞ is an unstated piece of
dictionary, and every proposed test presupposes it, including the
authors' own.

`duality/BOOST.md` calibrates that on the exactly solvable case. A single
mass produces a continuum with mass-encoding *nodes*, not a peak.

**Does the toy even have the disease?** QCD₂ is super-renormalizable. Its
only mass renormalization, m_q² = m_bare² − g²/π, is finite. So the 4d
fine-tuning drama of §5 is transplanted into a model where it barely
exists. And the vacuum-energy/entropy relation ε = 3/(8S) is only ever
written in D = 4, never derived for the JT-dS where the paper works.

`thooft-target/condensate.py` measures the 2d counterpart. It is a single
finite number: ε(0) = −0.0407 Nc g².

### The regime problem

The regime where the spectrum match lives is **fixed p, λ → 0**. That is
*not* the double-scaling regime where self-averaging is
super-exponential.

At fixed p the suppression is only the power law N⁻ᵖ. So gauge-symmetry
emergence is quantitatively weakest exactly where the match has to be
made. This repo makes that a printed, tested table rather than a remark.

---

## Part II — What is settled

### The yardstick is real

An independent re-derivation of the whole Galerkin assembly. A
re-verification of the hypersingular Chebyshev-U identity by regularized
quadrature. And a check of the anchor constants against FLZ's published
tables ([arXiv:0905.2280](https://arxiv.org/abs/0905.2280)), digit for
digit.

The worst deviation over n ≤ 29 is 4.5e-12, which is the level of FLZ's
own last printed digit at the top of the range. Against their Padé the
ground state agrees to 1.1e-15. The sum rules 7ζ(3) and 2 come out to
about 1e-10.

Off the duality point the Chebyshev basis hardwires the α = 0 endpoint
exponent ½, while the true exponent solves πβ cot(πβ) = −α. So
convergence there is only algebraic. `jacobi_solver.py` fixes that with a
matched-exponent basis [x(1−x)]^β: the ground state is below 1e-11 by
K = 32, and the LM sum rules come out to about 1e-12.

Eigenfunctions are checked (Parseval 3e-6, weak-form Rayleigh ~1e-6,
cross-solver overlap > 0.9999). FLZ's spectral-determinant constants d±
are consumed too (`flz_dconstants.py`): their Eq. 4.21 products
reproduce both 12-digit values at the anchor's own rounding floor, and
the exact ratio √2/π to 4e-14. Pinning their definition corrected a
mischaracterization — d± are spectrum functionals, not wavefunction
data — so the external eigenfunction-anchor slot is genuinely empty.
The decay constants are new data. No published wavefunction tables
exist.

### Self-averaging, exactly

`exact_wick.py` gives Var(B_jk) = 4σ⁴C(N−2,p−1) in closed form at any
(N, p). It is checked by literal contraction enumeration, by ED, and by
Monte Carlo.

In double scaling ln rms = −¼√(λN) ln N, which beats e^(−a√N) for every
a, scanned to N = 1600. The paper's Eq. 3.28 is not just right. It is
conservative. The singlet rate of Eq. 3.27 is exactly 2/C(N,p).

A second anchor, exact E[m₄] and Var(m₄), gives
relVar(m₄) = 4·relVar(m₂)(1+o(1)). The Dirac/U(N) counterpart is in
`exact_wick_dirac.py`, with rms ~ Nc^(−(p−1)/2).

**Worth telling the authors.** The "systematic derivation" that §3.5
defers to [arXiv:2501.09423](https://arxiv.org/abs/2501.09423) **is not
there**. That paper never mentions self-averaging. Its counting rule
(Eq. 4.26) is magnitude-only and part-empirical (their §6.2), and its n
≪ p approximation fails precisely in the diagram-B regime, with error
~e^p. That is harmless for the bound and fatal for prefactors.
`exact_wick.py` appears to be the only complete prefactor-exact
derivation. It agrees with both papers at leading order, and it adds the
exact pinning factor 4p(N−p)/N(N−1).

### The dictionary needs four corrections

Each one is forced, and each one is tested (`duality/dictionary.py`,
`syk-self-averaging/dirac.py`).

**M₁ = −(p/2)iH, not H.** It follows from anti-Hermiticity. Q̇ = 0 forces
M₁† = −M₁, so M₁ cannot equal the Hermitian H. Verified to 1.3e-15.

**M₂ = −Σċ†ċ.** The paper's + sign is wrong.

**CP is the unitary particle-hole map.** M₁ ∝ iH rules out an antiunitary
one.

**λ = p²/N against λ_chord = 2p²/N.** This factor of two is internal to
the Susskind corpus. [2501.09423](https://arxiv.org/abs/2501.09423) Eq.
3.6 uses 2p²/N.

Raw Mₙ for n ≥ 2 are CP mixtures with 15–65% wrong-channel content. So
CP = (−1)^(n+1) only lives in CP-resolved correlators.

### The conserved-channel theorem, and where the tower lives

This is the headline result, and its scope is the whole story. Read both
halves.

**The theorem.** The chord representation of every Aₙ = Σψ(ad_H)ⁿψ is a
polynomial in the chord Hamiltonian T. So every spectral moment of every
Aₙ **pair** channel vanishes at N = ∞. Every λ, every matter weight, both
orderings, odd k included.

The proof is a three-term recursion whose coefficients are scalars. Both
matching conditions reduce to q-number identities that do not depend on
the chord occupation. The explicit charges are Â₁ = −(1−q^Δ)T, which is
the N = ∞ image of the exact law A₁ = −2pH, and
Â₂ = (1−q^Δ)²T² + (1−q^(2Δ)).

Random O_Δ operators insert q^(Δn̂), which does *not* commute with T, so
they propagate. That contrast is what makes the result physical rather
than a normalization artifact.

**The scope, which is load-bearing.** This is the i ≠ j *pair* amplitude.
The full channel also has a flavor-diagonal part of the same order. The
per-pair amplitude is smaller by exactly the 1/N that the N(N−1)
multiplicity restores. And that diagonal part does **not** vanish for
n ≥ 2.

Evaluated exactly at finite N, the full A₂ width rises to w₂ ≈ 1.19 along
λ_chord = 2 (1.044 → 1.194 over (p,N) = (4,16)…(10,100)). At n = 1 the
two parts cancel to 3e-16, which reproduces A₁ = −2pH.

So the tower does propagate at N = ∞, and all of its weight is diagonal.

**Withdrawn.** The conclusion that the tower is a 1/N effect, invisible
at every order, and that the match must therefore be made at fixed p,
λ → 0, does not follow from this theorem. It is withdrawn. The fixed-p
argument survives separately, on the self-averaging side.

**The finite-N laws are now exact.** `exact_pair_moments.py` measures the
approach to the pair zero by Wick plus a CRT residue-table dynamic
program, with no ED. `pair_moment_laws.py` then turns those float
regularities into certified closed forms: μ₂ ~ N^−(n+2) exactly through
n = 4, with A(n=1) = −16p²(p−1) and A(n=2) = −64p³(p−1)(4p−5), and the
exact double-scaling constants c₂ = 28.92 → 30.53 over p = 4..16. Those
are the numbers the chord-side derivation has to reproduce.

### The dS axis

DSSYK∞ **scrambles at every temperature**. The growth coefficient is
a(v) = −1/cos(πv/2), real and negative, tending to −1 at the hologram
point, with |a| growing with β. A dS observer needs the opposite sign, so
the dictionary owes a sign-flipping map.

The Harlow–Zhao Euclidean fold is not that map. At ED sizes the folded
connected fraction grows from 0.365 to 0.786 over N = 8 → 20, while λ
shrinks. At N = ∞ it stays O(1) and grows as λ → 0, while the unfolded
one dies like λ. Neither statement needs a closed form. The mechanism
behind the bad finite-N convergence is measured: folded correlators at ED
sizes track the sample's upper spectral edge at correlation −0.86.

One candidate map does work, semiclassically. The antipodal
half-fake-period shift Δ = π/v flips the sign exactly in the closed form,
and for λ ≤ 1.7 in the exact theory. The same shift makes the Wightman
sector exactly thermal at β_eff = β + β_fake. One map, both of
Cui–Kolchmeyer's dS-observer signatures. Geometrically it is NV's ±iβ/4
dressing, and `duality/NV-DRESSING.md` pins that identification at the
contour level with equation-level citations.

Read the closed-form flip for what it is. Δ enters the growing mode's
phase, so half a period negates it, and the algebra is short. The content
is in the exact-λ boundary and in the deep campaign.

That campaign finds a transient anti-scrambling epoch growing like
~½ln(1/λ), not eternal dS. It locates the chaos-bound evasion: the phase
winds at constant magnitude, so there is no strip-KMS structure. And it
records the costs. The shifted contour is fold-adjacent and
edge-dominated, and rate extraction is mode-contaminated.

### Other results

`freeness.py` exhibits the free-product transition. spec(E + F(t)) lands
on the arcsine law, and the W₁ plateau falls like dim^(−0.9). That is the
first time in a microscopic dS-candidate model.

`qed_campaign.py` closes the gap between Majorana numerics and Dirac
claims. `baryons.py` gives first data on the wholly underived
inter-baryon-force claim. `condensate.py` transcribes the exact
condensate and gives first α → −1 GMOR data.

---

## Part III — What is open, and what each would take

### The tower extraction (the headline, and a live race)

Two halves, and they are different computations.

The **pair channel** is a 1/N computation from the double-scaled
viewpoint. The **diagonal channel**, which carries the tower's N = ∞
weight, is a strict-N = ∞ chord computation.

**The 1/N law at fixed p is now exact**, and needs no chord derivation
(`pair_moment_laws.py`). P(N) = 2μ_k C(N,p)^m is a polynomial in N, so
interpolation gives the closed forms, certified out of sample.
μ₂ ~ N^−(n+2) exactly through n = 4, with A(n=1) = −16p²(p−1) and
A(n=2) = −64p³(p−1)(4p−5) certified at p = 14. μ₄ dies one power faster,
deficit n+3, sharing its leading coefficient with μ₂(n+1).

What remains for theory: the chord-side derivation that *reproduces*
these exact targets, and the large-p limit of the double-scaling
constant. c₂(p) = 28.92 → 30.53 over p = 4..16 is exact, and still
rising.

**Extending to μ₆** needs K = 12 H-node words. That is unlocked now — the
exact-Fractions assembly removed the float N = 256 cap and raised the
word guard — but the (K−1)!! matching enumeration (10395 at K = 12, with
m = 6 chord canonicalization) makes it hours per point in pure Python.
n = 4 (K = 10) is done.

**The flavor-diagonal channel at N = ∞ is done**
(`duality/diagonal_chord.py`). It is the parity-law diagram sum, derived
from the finite-N sign law and triple-gated: walker against brute force,
a transfer-matrix anchor, and finite-N convergence at the O(1/p) line
rate. The diagonal A₂ width at λ_chord = 2 is w₂ = 4.67(9) at Δψ → 0. The
full-channel limit of ≈ 1.31(5) is a diagonal-against-pair cancellation
whose missing chord-side piece is the pair channel's 1/N coefficient.

**Moments → spectrum.** The forward direction was run first
(`duality/boost_forward.py`): boost-kernel moments diverge, S → 1/ω, the
regulated invariants sweep the whole SYK band as the regulator dial
turns, and the coupling-weighted tower sits within 5% of single-mass at
the invariant level. The ε* consistency table in `boost_forward.json` is
the recorded spec.

**"How many moments" is now answered on the solvable kernel**
(`duality/boost_inverse.py`). Exactly as many as there are levels — the
singular spectrum of the weights→moments map saturates at K = L — and
each level costs 0.8–2.4 digits of moment precision more than the last
(log₁₀(σ₁/σⱼ) ≈ a(j−1) + b(j−1)², a = 1.54 → 0.44 over aε = 0.1 → 1.6,
b ≈ 0.05). Double precision buys 9–12 of 12 levels; 8 digits buy 6–10.
Resolving a δ-split pair takes a third moment known to ~0.05 δ².
Decisively: moment invariants carry regulator leverage ≈ 1.0, so with
the state map unfixed they carry zero mass digits, while node positions
carry exactly zero leverage. The "needs node-level data" verdict is now
quantitative, conditional on the e^(−εω) regulator and the free kernel.
Which finite-N object carries the tower remains open.

**The Dirac/complex-SYK chord machinery** has its H-sector built
(`duality/dirac_chord.py`, on Berkooz–Narovlansky–Raj
[2006.13983](https://arxiv.org/abs/2006.13983)). The chord sum is the
Majorana one at q = e^(−p²/Nc) = e^(−λ_paper), so the factor-2 landmine
is Majorana-only. It is verified microscopically and against ED (m₄/m₂²
deviation +0.198 → +0.004 over Nc = 8 → 12, c_symmetric). M₁ = −(p/2)iH
is a transfer-matrix insertion, so it is an exactly conserved channel.
Charged random probes carry Δ = (a+b)/p, with fixed-charge q(𝒬) and BNR
asymmetry factors.

The pair (i ≠ j) half of the M_{n≥2} chord question is now closed
in-repo. The tower's pair words carry bare single-fermion chords, the
lone-crossing census is exactly 1 − p/Nc → q_dirac^(1/p), and the
conserved-channel theorem applies verbatim at (q_dirac, Δ = 1/p): every
Dirac pair moment vanishes at Nc = ∞, with X̂₁ = −(1−q^(1/p))T the image
of M₁ = −(p/2)iH. No published source provides H-correlated chord
rules, and the pair half did not need them.

Still open, and what the CP-resolved data is actually waiting on: the
flavor-diagonal channel at the Dirac q, and the CP resolution of the
Dirac Mₙ.

Comparison to `thooft-target` then goes through the BOOST.md kernel — a
mass is a Bessel-kernel continuum with nodes, not a peak — and through
`dictionary.py`.

### The O(1/q) residues

The large-q singlet 4-point function has Sommerfeld–Watson poles in only
two places: m = 0, which is Hamiltonian exchange and the entire
time-ordered part, and m = ±v, the scramblon.

The would-be tower cancels identically at leading order, and its
couplings are O(1/q²). It first appears at O(1/q), organized on the
exact-in-temperature skeleton 𝓜ᵉ/𝓜ᵒ, whose finite-temperature residues
are unpublished. Computing them is the open task. And how the 1/q
expansion at fixed λ relates to the 1/N one deserves a clean answer.

The skeleton itself is now a derived object rather than a transcription
(`duality/largeq_kernel.py`): the Maldacena–Stanford finite-T kernel
eigenproblem at k = 1, discretized by parity-split sine-Galerkin,
reproduces the CMS quantization sets to 3e-12 across v = 0.05..0.95, and
locates m = 0 and m = ±v as the exact non-Dirichlet solutions whose
eigenfunctions degenerate. That is the perturbable substrate the O(1/q)
computation needs; the residues on it are still the open task.

### The tomperature map

The antipodal shift passes both necessary conditions semiclassically. But
a sign flip plus effective thermality under a contour shift is not a
dictionary.

What would upgrade it is the algebra story: why the observer's operators
live at the antipodal point. That is exactly NV's construction to
complete.

### The NV discriminating computation

The OTOC sign of Narovlansky–Verlinde's dressed operators. They engineer
KMS at the dS temperature, which is the no-go's hypothesis, so they must
break another hypothesis of the bound. This is the sharpest next
computation on the dS axis, and it needs the two-sided constrained chord
sector, which is adjacent to the existing engine.

`duality/NV-DRESSING.md` now pins the identification at the
contour level: their Eqs. 33–34 and 61–62 against the fake-period shifts,
with the composition algebra and the 2-point thermality tested. What
remains is the constrained-sector step, which is whether NV's BRST
projection renormalizes the crossed growth coefficient. Conditional on it
acting trivially, the answer is already measured: anti-scrambling,
semiclassical, transient.

### The bulk side

**JT-dS ε = 3/(8S) — settled as far as it settles**
(`duality/jt_ds_vacuum.py`). The general-D law is
ε_D ≡ ρG^(D/(D−2)) = c_D S^(−2/(D−2)) with c₄ = 3/8 exactly,
sympy-certified; the paper's relation is its D = 4 instance, and the
1/S power is dimension-specific. The direct D → 2 continuation
degenerates fourfold (no dimensionless ρG^a, G_μν ≡ 0, Λ and c_D vanish,
the ℓ-elimination dies), machine-checked. What JT-dS actually supplies
is the exact first law δS = −2πℓ²ρ = −δE/T_dS (the dilaton shifts, ℓ is
rigid; the sign of SVVV [2203.00700](https://arxiv.org/abs/2203.00700)
Eq. 9), plus the reduction-inherited family — S⁻² from a dS₃ parent,
machine-checked end to end; S⁻¹ from near-Nariai dS₄ (MTY
[1904.01911](https://arxiv.org/abs/1904.01911)) — whose exponent the
intrinsic JT data (φ₀, φ_r, ℓ) cannot choose. That parent choice is a
named dictionary input, alongside the tomperature map. The measured 2d
condensate now has a first-law image, δS = +2π(0.0407)N_c(gℓ)²: finite,
with nothing for Holographic Tuning to cancel, on both sides of the
first law.

**The finite-N coupling curve ḡ²(N).** Their Fig. 17 is a sketch with no
formula, yet the whole set-it-and-forget-it inversion rests on it. Match
a bulk observable at a sequence of finite N and extrapolate. Gated on the
tower extraction.

**The order-of-limits problem.** T_B → ∞, N → ∞, λ fixed against λ → 0,
and the boost-spectrum → invariant-mass reorganization. The
chaos-to-integrability transition at small λ
([2403.01950](https://arxiv.org/abs/2403.01950),
[2601.09801](https://arxiv.org/abs/2601.09801)) is the likely mechanism
zone. This is where the duality most plausibly *fails*, and a sharp
statement of where it holds would be a real contribution.

**Freeness sharpening.** The ln N onset-time question needs sizes beyond
ED.

### Out of scope

Proving DSSYK∞ = JT-dS. The authors assume it themselves.

The 1+1 → 3+1 uplift. Contrary to how it is sometimes described, *no*
uplift construction exists in the source. §§5–6 are explicit that it is a
hope.

Instance-level ED in the true double-scaling regime (p ~ √N, N ≳ 100,
dim 2⁵⁰). That is the honest wall. The exact-Wick and chord computations
are the ways around it.

---

## Context

The QCD₂ side is effectively solved beyond FLZ. Ambrosino–Komatsu recast
't Hooft's equation as a TQ-Baxter integrable system
([2312.15598](https://arxiv.org/abs/2312.15598),
[2406.11078](https://arxiv.org/abs/2406.11078)), and
Litvinov–Meshcheriakov ([2409.11324](https://arxiv.org/abs/2409.11324))
give independent high-precision spectra. So numerical effort should go
almost entirely to the DSSYK side.

The DSSYK toolkit is mature. Chord techniques
([1811.02584](https://arxiv.org/abs/1811.02584); review Berkooz–Mamroud
[2407.09396](https://arxiv.org/abs/2407.09396); path integral
[2403.05980](https://arxiv.org/abs/2403.05980)). And the
chaos-to-integrability transition at small λ is directly reusable for
the λ → 0 regime the duality needs.

The duality sits in contested territory. Narovlansky–Verlinde
([2310.16994](https://arxiv.org/abs/2310.16994)) and the sine-dilaton
program of Blommaert–Mertens–Papalini
([2404.03535](https://arxiv.org/abs/2404.03535)) propose *different*
bulk duals. Okuyama ([2505.08116](https://arxiv.org/abs/2505.08116))
reaches dS JT gravity from DSSYK by a further route. And Rahman–Susskind
([2407.12988](https://arxiv.org/abs/2407.12988)) argues the
identifications conflict. A clean meson-tower match, or mismatch, would
be the first quantitative discriminator among all three programs.
`duality/DISCRIMINATOR.md` lays out that axis.

**Convention landmine.** This paper uses λ = p²/N with no q. The chord
literature uses λ = 2p²/N, q = e^(−λ). Any reuse of chord results must
pass through a tested conversion.

**Priority.** As of 2026-08-24 the target paper has exactly one
citation: Cui–Kolchmeyer [2607.13665](https://arxiv.org/abs/2607.13665),
which argues that dS observer OTOCs *anti-scramble*, and poses a
structural challenge to observer-centric static-patch holography. Its
same-week companion Harlow–Zhao
[2607.14215](https://arxiv.org/abs/2607.14215) proposes the
Euclidean-fold mechanism, and poses the SYK test that
`duality/ANTISCRAMBLING.md` answers.

The authors' promised meson-spectrum computation has not appeared. The
first-mover race on the headline test is still open.
