# Road map — stress-testing DSSYK∞ ↔ 't Hooft (arXiv:2607.05678)

*Written 2026-07-17, after a full multi-agent re-review: both papers' claims re-read from
source, both test suites re-run, every top-level README claim audited against paper + code,
and a literature sweep (INSPIRE, arXiv) as of this date.*

---

## Part I — Where the project stands

### The paper, honestly assessed

*Holograms and Standard Models* (Miyashita–Sekino–Susskind, July 2026) is a short,
talk-style manifesto, not a technical paper. Its chain of claims:

1. **Assumed** (their words: "unproved but we will assume it"): DSSYK∞ = JT-de Sitter
   static-patch holography.
2. **Quasi-derived** (§3.5, two example diagrams; systematics deferred to arXiv:2501.09423):
   perfect self-averaging, Var⟨W⟩ ≈ p!/Nᵖ < e^(−a√N), hence *emergent* O(N)/U(N) symmetry
   of single disorder instances in the flat-space limit.
3. **Argued by counting** (entropy match S = N log 2 ⇒ dim matter = N): the bulk dual is
   the 't Hooft model — U(N) QCD₂ + QED, one fundamental flavor, at renormalized quark
   mass m_q² → 0, with dictionary N_ym = N_syk, ᾱ = ḡ²N = p², ḡ² = λ = p²/N.
4. **Qualitative only** (their words: "We don't have a quantitative calculation of the
   DSSYK spectrum"): the singlet tower Mₙ = ψ†(dⁿ/dtⁿ)ψ ↔ 't Hooft mesons, CP = (−1)ⁿ⁺¹,
   with n=0 (photon = Q) and n=1 (graviton = H) non-dynamical. The paper names the
   quantitative spectrum match as "the best test of the duality." **It has never been run,
   by anyone.**
5. **A logical inversion, not a computation** (§5): "Holographic Tuning" — at T_B = ∞ the
   entropy is exactly log dim H, so the vacuum energy gets no corrections; the coupling
   ḡ²(N) is a *derived* quantity (its finite-N form is only a sketched figure).

Two structural weaknesses the review surfaced that the paper does not address:

- **The boost-vs-mass gap.** The SYK Hamiltonian is the Rindler *boost* generator, and a
  fixed-mass particle has a *continuous* boost-frequency spectrum. How discrete 't Hooft
  masses are supposed to appear as features of Mₙ correlators at T_B = ∞ is an unstated
  dictionary that every proposed test (including the authors') presupposes.
- **Does the toy even have the disease?** QCD₂ is super-renormalizable — its only mass
  renormalization (m_q² = m_bare² − g²/π) is finite. The 4d fine-tuning drama of §5 is
  transplanted into a model where it barely exists, and the vacuum-energy/entropy relation
  ε = 3/(8S) is only ever written in D=4, never derived for JT-dS where the paper works.

Also relevant: the regime where the spectrum match lives (**fixed p, λ → 0**) is *not* the
double-scaling regime where self-averaging is super-exponential — at fixed p the
suppression is only the power law N⁻ᵖ. Gauge-symmetry emergence is quantitatively weaker
exactly where the match is to be made. Nobody (paper or repo) has drawn this consequence.

### The code, verdicts confirmed by independent re-derivation

**`thooft-target` — solid. The yardstick is real.**
An agent re-derived the whole Galerkin assembly by hand, re-verified the hypersingular
Chebyshev-U identity by independent regularized quadrature, and — critically — fetched
FLZ (arXiv:0905.2280) itself and extracted the PDF text: the anchor constants in the tests
match FLZ's published tables digit-for-digit (they are not self-generated numbers labeled
"FLZ"). Fresh recomputation: worst FLZ-table deviation 2.0e-12 (at the *rounding floor of
FLZ's own tables*, so "12–13 digits" is a floor set by the reference, not the solver;
ground state actually agrees to ~8.5e-15), sum rules 7ζ(3) and 2 reproduced to ~9e-11.
13/13 tests pass in ~4 s. Real weaknesses: for α ≠ 0 the basis hardwires the α = 0 endpoint
exponent ½ while the true exponent solves πβ cot(πβ) = −α, so off the duality point
convergence is only algebraic (the README's "should be equally accurate" is not right);
eigenvectors are computed but never validated; `thooft_highprec.py` is needlessly O(K³)
slow; a latent divide-by-zero sits in `_f()` at k = ±1.

**`syk-self-averaging` — honest, correct core, weaker statistics than advertised.**
The Majorana/Jordan–Wigner construction, parity handling, and N mod 8 RMT validation are
correct (GSE Kramers doubling verified explicitly; parity-sector subtlety handled right).
The exact identity relVar(m₂) = 2/C(N,p) and its relation to the paper's p!/Nᵖ (factor
Π(1−k/N)⁻¹ ≈ e^(+p²/2N)) check to machine precision. 7 tests pass (+1 slow deselected).
Real weaknesses: **no error bars anywhere** — with n = 200 realizations the sampling noise
is ~10%, and actual ED deviations from 2/C(N,p) are 2.4–9.4% (the README's "~1%" is wrong);
the fit-quality model comparison is statistically meaningless (5 points, <1 decade; at
p = 6 `results.txt` even reports best fit e^(−aN), contradicting the target form, without
comment); `plot_self_averaging.py` plots the ratio metric the module *itself* declares
untrustworthy; Hamiltonian assembly by sparse matmuls is 10–100× slower than a bitwise
Pauli-string method (this caps N and realization counts); and the entire Dirac/complex
sector is dead code — no complex H is ever built, `test_charge_conserved_complex` never
tests [H,Q] = 0. Consequence: **every executed numeric is Majorana/O(N), while every
headline claim is Dirac/U(N) QCD+QED** — the U(1)/electromagnetic half has zero coverage.

### The top-level README — audit results

All four numbered "load-bearing corrections" are **CONFIRMED** (the M₁ = −(p/2)iH identity
was verified numerically to 1.3e-15 using the repo's own `dirac_operators` — it follows
from anti-Hermiticity: Q̇ = 0 forces M₁† = −M₁, so M₁ cannot equal Hermitian H). But the
audit found real defects to fix:

- **REFUTED as worded:** "The infinite-T kernel is not in print." The exact DSSYK 4-point
  function at infinite temperature *is* published (Berkooz–Isachenkov–Narovlansky–Torrents,
  arXiv:1811.02584, chord transfer-matrix), and large-q 4-point functions exist in closed
  form at all temperatures (Streicher 1911.10171; Choi–Mezei–Sárosi 1912.00004).
  "The SYK ladder is always solved in the low-T conformal regime" is also false
  (Maldacena–Stanford treat general βJ). What survives, and matters: **no published work
  extracts the bilinear-tower masses from it, and discrete-vs-continuum is genuinely open.**
  This is good news — the headline computation starts from a published object, not a new
  kernel solve.
- Misattributed factor: the paper's Eq. (3.27) reads p!/Nᵖ with *no* factor 2; the 2 in
  "the paper's 2p!/Nᵖ" is the repo's exact prefactor.
- Fabricated quote: 2506.18054 §5.3 actually says "At the moment we don't know the masses
  of the states created by the operators in (5.20)" — the README's quoted phrase "are not
  currently known" appears in neither version. (Substance correct, quote wrong.)
- Notation: "ᾱ = g²N = p²" conflates the paper's dimensionful α = g²N with dimensionless
  ᾱ = ḡ²N; correct chain is ᾱ = ḡ²N = p².
- Overcompressed: "ED … down to ~10⁻⁴⁶ at N=640" — that number is pure combinatorics
  (2/C(640,26)); ED caps at N ≈ 18–20. "All four GOE/GUE/GSE classes" is three classes
  across four N mod 8 values. And "the validated Var⟨W⟩ < e^(−a√N) applies…" contradicts the
  module's own honesty note — what is validated is the 2/C(N,p) rate; e^(−a√N) is exactly
  what ED could *not* reach.

### Context (literature as of 2026-07-17)

- The paper has **0 citations** and no blog/forum coverage; the only public discussion is
  Susskind's IAS talk (May 2026). The authors promise the meson-mass computation in
  "a future publication" — **there is a first-mover race on the headline test.**
- The QCD₂ side is effectively *solved* beyond FLZ: Ambrosino–Komatsu recast 't Hooft's
  equation as a TQ-Baxter integrable system (2312.15598, 2406.11078); Litvinov–Meshcheriakov
  (2409.11324) give independent high-precision spectra. Numerical effort should go almost
  entirely to the DSSYK side.
- The DSSYK toolkit is mature: chord techniques (1811.02584; review Berkooz–Mamroud
  2407.09396; path integral 2403.05980), and the chaos-to-integrability transition at small
  λ (2403.01950, 2601.09801) is directly reusable for the λ → 0 regime the duality needs.
- The duality sits in **contested territory**: Narovlansky–Verlinde (2310.16994) and the
  sine-dilaton program (Blommaert et al. 2404.03535, 2505.08116) propose *different* bulk
  duals, and Rahman–Susskind 2407.12988 argues the identifications conflict. A clean
  meson-tower match (or mismatch) would be the first quantitative discriminator among all
  three programs — high payoff either way.
- ⚠ Convention landmine: this paper uses λ = p²/N with no q; the chord literature uses
  λ = 2p²/N, q = e^(−λ). Any reuse of chord results must pass through a tested conversion.
- **Update 2026-07-23:** the paper is no longer uncited. Cui–Kolchmeyer 2607.13665
  ("A de Sitter Anti-Scrambling Algebra", MIT/IAS) cites it as [12] while arguing dS
  observer OTOCs *anti-scramble* (negative eikonal phase, KMS failure at the dS
  temperature, Hartle–Hawking not a trace on the crossed product) — a structural
  challenge to observer-centric static patch holography. The same-week companion
  Harlow–Zhao 2607.14215 proposes the QM mechanism: bounded-spectrum systems with
  Euclidean-folded correlators — DSSYK's defining features — and poses the fold test
  "in a concrete model such as the SYK model" as open. Both papers on arXiv
  (local copies: papers/); this
  repo's response is rungs 25–26 (`duality/ANTISCRAMBLING.md`). Time-ordered physics
  (the entire spectrum-match program, rungs 13–18) is explicitly untouched by these
  papers; the OTOC sector is a NEW discriminator axis for rung 24.

---

## Part II — The ladder

Ordered from an afternoon's work to open research problems. Each rung is (mostly)
independent; the tiers gate the ones that aren't.

### Progress log (2026-07-20; all reviewed by an adversarial 41-agent pass, findings applied)

- ✅ **1–4** (Tier 0): git history; README defects fixed; deps pinned; `pyproject.toml`;
  solver hardening (guards, interleave assertion, ddof=1, JSON-driven plots).
- ✅ **5–6**: `dirac.py` — complex SYK, [H,Q]=0, M₀=Q, M₁=−(p/2)iH; *new*:
  M₂=−Σċ†ċ (paper's + sign wrong), Mₙ†=(−1)ⁿMₙ only n≤2; CP is the *unitary*
  particle-hole map (law: C cᵢC† = (−1)^(Nc−1)cᵢ†); CP(M₀)=−1, CP(M₁)=+1 exact;
  raw n≥2 are CP-mixtures ⇒ the sharp claim lives in CP-resolved correlators.
- ✅ **7–8**: error bars + ΔAIC everywhere; bitwise Pauli-string assembly (10–100×,
  suite 2 min → 4 s).
- ✅ **9**: eigenfunctions validated (Parseval 3e-6; weak-form Rayleigh ~1e-6;
  cross-solver overlap >0.9999; decay constants — all fₙ>0 under φ(0⁺)>0, new data).
  External anchors: FLZ s=3,4 to 3e-12; LM α≠0 tables at their rounding floor; LM
  exact α-dependent sum rules to 1e-5–2e-4. **α≠0 gap closed.**
- ✅ **10**: `duality/` dictionary + the two parameter-free fingerprints, tested.
- ✅ **11** (also answers 16 analytically): `exact_wick.py` —
  Var(B_jk)=4σ⁴C(N−2,p−1) exact ⇒ **Eq. 3.28 verified in double scaling**
  (ln rms = −¼√N lnN, scanned to N=1600); fixed-p asymptote N^(−(p−1)/2);
  E[B_jj]∝(1−2p/N) derives the p=N/2 artifact.
- ✅ **14** (first version): `mn_spectroscopy.py` — first CP-resolved Mₙ spectral
  functions (single realization, λ=1.6): conserved n=0,1 machine-exact at ω=0;
  broad continua for n≥2. Deliberately *not* the λ→0 match.
- ✅ **12** (audited 2026-07-20): the "systematic derivation" behind §3.5 is
  published **nowhere** — 2501.09423 never mentions self-averaging, its counting
  rule (Eq. 4.26) is magnitude-only, part-empirical (their §6.2), and its n≪p
  approximation fails exactly in the diagram-B regime (error ~e^p; harmless for
  the bound, fatal for prefactors). `exact_wick.py` is the only complete
  prefactor-exact derivation we're aware of; it *confirms* both papers at leading
  order (their diagram-B = σ⁴C(N,p); ours adds the exact pinning factor
  4p(N−p)/N(N−1)) — a small publishable note. Also: the λ factor-2 clash is
  internal to the Susskind corpus (2501.09423 Eq. 3.6: λ=2p²/N). And
  2511.10907's entropy correction (T_c = Λ√N per McLerran–Sen; stretched
  horizon at Planck, not string, depth) leaves §5's Holographic Tuning
  untouched — 2607.05678 already cites and incorporates it; but 2501.09423
  §6.1–6.2 is superseded, and rungs 19/21 must use T_c = Λ√N.
- ✅ **13** (calibration level): `duality/boost_mass.py` + `BOOST.md` — exact
  boost spectral function of one mass, S = (1/πa)e^(πω/a)K_{iω/a}(m/a)²,
  verified to 1e-12 by an independent route. One discrete mass ⇒ a *continuum*
  with mass-encoding *nodes* (peak-hunting is the wrong observable); heavy
  states Boltzmann-faint; the ω-symmetric T_B=∞ hologram vs KMS-asymmetric
  bulk response is the precisely-stated open map (tomperature). Rung 18's
  "discrete vs continuum" question reformulated accordingly.
- ✅ **15** (groundwork): `moments_pipeline.py` — disorder-averaged CP-resolved
  spectral-moment shape invariants (w₂, kurtosis) with errors along fixed
  p = 4, Nc = 6→12 (Nc = 12 via the charge-sector-blocked moments_fast.py,
  landed 2026-08-17 from ws/moments-fast after review; dense equivalence
  asserted at 1e-9, measured 5e-15); the map-robust targets for the exact N=∞ chord/Wick
  computation. (Also fixed en route: for n ≥ 3, Tr(M(t)M†) is genuinely
  ω-asymmetric since M_n is not ±Hermitian — the symmetrized correlator is the
  convention-stable object.)
- ✅ **17** (structure settled; residues = the remaining work):
  `duality/largeq_anchor.py` + `LARGEQ.md` — equation-verified large-q closed
  forms with the exact conventions bridge (𝒥 = √2·p; G(t) = sech(√2pt)^(2/p)
  matches ED with no free parameters at O(1/p)). **Finding: the λ→0 anchor
  cannot come from the leading large-q 4-pt function** — its SW plane has
  poles only at m=0 (Hamiltonian exchange = the entire TOC) and m=±v
  (scramblon); the would-be tower (quantization sets 𝓜^e/𝓜^o, the
  all-temperature k(h)=1) cancels at leading order and its couplings are
  O(1/q²) (GR 2.15). The tower's even/odd alternation = the CP slot; at T=∞
  the sets → integers (consistent with BOOST.md). Rung 17 proper = computing
  the unpublished finite-T O(1/q) residues on that skeleton. En route: the
  Majorana M₁ analog Σψ[H,ψ] = −2pH (asserted); ED bench `largeq_bench.json`.
- ✅ **19** (2026-07-21, adversarially reviewed): `qed_campaign.py` — the
  Dirac/QED campaign. U(N)-class self-averaging made exact (Dirac analog of
  the Wick closed forms, verified by enumeration + ED); disorder-averaged
  charging curves E0(q) with full mean-vector covariance and GLS fits
  (capacitive vs confining stated via ΔAIC, not asserted); charged-vs-singlet
  spectral scales. Closes the "Majorana numerics vs Dirac claims" gap.
- ✅ **23** (2026-07-21, adversarially reviewed): `baryons.py` — baryon-sector
  first contact. ε-tensor baryons at Nc=4–8: exact 2-point collapse onto
  charge sectors, M_B extensivity, edge-curvature hierarchy; hierarchy claims
  kept honest (spectral widths, rejected-fit error inflation). First data on
  the wholly-underived inter-baryon-force claim.
- ✅ thooft polish (2026-07-21, adversarially reviewed): `jacobi_solver.py` —
  matched-exponent Jacobi basis [x(1-x)]^β, β from πβcotπβ=−α: spectral
  accuracy at α≠0 (ground state <1e-11 by K=32 vs the Chebyshev solver's
  algebraic K^(−4β)); β=½ anchor vs closed form ~5e-13; LM sum rules s≥4 to
  ~1e-12. Removes the α≠0 precision limitation flagged in Part I.
- ✅ **25** (2026-07-23): `duality/antiscrambling.py` + `ANTISCRAMBLING.md` —
  the OTOC sign test forced by Cui–Kolchmeyer 2607.13665. Closed form
  (Streicher transcription, LaTeX-source level): the growing OTOC term has
  Re a(v) = −1 EXACTLY at every temperature, λ_L → 2𝒥 at T_B=∞ — the
  scrambling/AdS sign (CK's c>0). ED independently confirms the OTOC
  magnitude *growth* at every measured (N, βJ); the scrambling *sign* itself
  is analytic (Streicher-transcribed), inheriting an empirically-fixed
  crossed-operator ordering, not an independent ED sign measurement. A dS
  observer needs c<0: **under the flat dictionary DSSYK∞ cannot
  anti-scramble** — the tomperature map (13 remainder) is now gated by
  "must flip the OTOC sign," a sharper statement.
- ✅ **26** (first data, 2026-07-23): `syk-self-averaging/fold_bench.py` —
  the Harlow–Zhao Euclidean-fold test (2607.14215 §6) they explicitly pose
  for SYK. Folded correlators are finite (bounded spectrum, as HZ require)
  and the fold's Lorentzian response tracks the analytic continuation
  (Im 𝓕/F_d to ~10% at the best point, sign/trend everywhere); but the
  folded Euclidean 4-pt drifts AWAY from the continued closed form as N
  grows at fixed p (ratio 1.05→1.21→1.42 for N=12→14→16) while unfolded
  observables converge — λ-corrections vs a genuine Stokes obstruction is
  not separable at ED sizes. Deciders: fixed-λ scan; chord-side folded
  4-pt at N=∞ (adjacent to rung 18 tooling).
- ✅ **26 decider #1** (2026-08-12): `fold_bench.py --scan` — the fixed-λ
  scan, p=4 N=8..20 (λ=16/N ∈ [0.8,2]) + matched-λ anchor
  (p=6,N=18)↔(p=4,N=8). Verdict: **a genuine N-obstruction, not a
  λ-artifact.** Control: the unfolded OTOC ratio is N-stable (1.80→1.84)
  while the folded ratio explodes (1.52→8.19); the drift GROWS as λ shrinks
  (opposite to a semiclassical correction); at matched λ=2.0 the larger-N
  point drifts further (1.700±0.004 vs 1.522±0.027) despite smaller 1/p
  error; and the folded connected fraction (F−F_d)/F_d reaches O(1)
  (−0.37→−0.79) — the 1/N hierarchy itself collapses under the fold. At ED
  sizes HZ's "continuation commutes" assumption fails, increasingly badly
  with system size; honest cap N ≤ 20 stated (crossed 4-pt is O(N²dim³)).
  Final word at N=∞ stays with the chord-side fold (decider #2, rung-18
  toolchain).
- ✅ **16** (numerical link, 2026-08-06): `exact_wick.py` `fixed_p_scan` —
  the fixed-p, λ→0 map made explicit (the analytic answer was rung 11's):
  local slopes converge onto N^(−(p−1)/2) with prefactor 2p!/√((p−1)!);
  per-element violation rms(B_jk)/E[B_jj] ~ N^(−(p+1)/2) (10⁻⁶ needs
  N ≈ 826 at p = 4). The Part-I consequence ("emergence is weakest where
  the match is made") is now a printed, tested table, not a remark.
- ✅ **27** (first data, 2026-08-12): `syk-self-averaging/freeness.py` —
  the freeness onset CK pose (their §4), measured. For the bosonic pair
  E = iψ₀ψ₁, F(t) = iψ₂(t)ψ₃(t), freeness ⇔ spec(E+F) = arcsine law
  (Bernoulli ⊞ Bernoulli); classical (t=0, exact) = binomial. At
  N = 12..24: the spectrum LANDS on the arcsine law — W1 plateau falls
  strictly, ~dim^(−0.9) (rigid-spectrum 1/dim floor, not sampling
  dim^(−1/2)); all alternating words → 0. Nobody had exhibited the
  free-product transition in a microscopic dS-candidate model. Honest
  caveats: the distributional crossover sits at the DISSIPATION time
  (t𝒥 ≈ 1.2, nearly N-independent) — ED cannot resolve the ln N
  scrambling-time separation where CK's statement properly lives; word
  plateaus at N ≥ 20 are SEM-limited.
- ✅ **18 tooling + first deliverable** (2026-08-16): `duality/chord.py` —
  the exact N=∞ chord engine (transfer matrix + matter-chord contour
  evaluator, no Askey–Wilson functions), validated three independent ways
  (brute-force chord-diagram enumerator; published closed forms (2.30)/
  Touchard–Riordan; chord-free Krawtchouk microscopy of q = e^(−2p²/N),
  deviation O(p⁴/N²) measured) + ED anchors (thermal 2-pt to 0.6–1.6% at
  N=12–16 — vs the ~30–55% large-q floor). λ→0 scan (`chord_spectra.json`):
  every channel a STRUCTURELESS CONTINUUM at every λ down to 0.05 — zero
  maxima, zero nodes at the semiclassical scale — converging to the
  sech^{2Δ} transform linearly in λ; ⟨n(t)⟩ length growth exact-in-λ.
- ✅ **15 proper** (2026-08-16): `duality/chord_moments.py` +
  `syk-self-averaging/mn_majorana_bench.py` — exact N=∞ moments of the
  H-derived singlet channels A_n = Σψᵢ(ad_H)ⁿψᵢ via chord words (assembly
  verified as a per-instance operator identity): **μ₂ = μ₄ = 0 to machine
  precision, n = 1,2,3, every λ ∈ [0.1, 4], any matter weight** — the
  singlet bilinear tower carries NO propagating weight at strict double
  scaling (A₁ = −2pH is the exact anchor; the O_Δ random channels DO
  propagate: μ₂ = 2(1−q^Δ), kurt → 3+1/Δ). ED bench (p = 4, 6; N ≤ 22):
  raw w₂ RISES along fixed p toward the λ→0 corner (the tower turns on
  exactly where the match lives). THE TOWER IS A 1/N EFFECT INVISIBLE AT
  EVERY ORDER OF THE DOUBLE-SCALED EXPANSION — the exact-in-λ
  generalization of rung 17's cancellation; the match lives strictly
  outside strict double scaling. Adversarial-review corrections
  (2026-08-16): the chord moments are i≠j PAIR amplitudes while the
  dense-ED channel is flavor-diagonal-dominated (pair μ₀ a few percent,
  opposite sign at N=12–14), so the bench does NOT yet probe the chord
  zero — at the matched-λ pair the raw w₂ rises with N (honest negative,
  asserted as such); a former (1−2p/N)² "artifact correction" was removed
  (its motivating zero is an exact per-instance conservation SPECIAL to
  (N,p)=(8,4), mechanism unknown — not a p=N/2 law); firm fixed-λ
  asymptotics need N ≫ 2p, beyond dense ED.
- ✅ **26 decider #2** (2026-08-16, the fold question CLOSED):
  `duality/chord_fold.py` — the HZ fold exactly at N=∞ (fold_bench
  conventions mirrored; the crossed transposition minus DERIVED as the
  fermionic matter-chord crossing sign). Unfolded control: λ-stable
  (−1.65…−1.67), exactly real, ~10% from matched-λ ED. Folded: the
  connected fraction stays O(1) and GROWS as λ→0 (−0.21→−0.33 over
  λ = 2.7→0.5) — no 1/N suppression at N=∞ — so the N-scaled folded ratio
  runs from the continued closed form ∝ 1/λ (−1.86→−21.3 vs −1.92
  constant): HZ's "continuation commutes" fails IN THE DOUBLE-SCALED
  THEORY ITSELF (both branches: the folded Lorentzian Im amplification
  grows ~1×→~45× over λ = 4→0.5, ∝ 1/λ). ED sits ×1.5–1.6 above the
  exact N=∞ folded values at matched λ (unfolded ≤ 10%): the fold
  amplifies finite-N corrections — mechanism MEASURED (`fold_edge_probe`:
  per-instance folded fraction tracks the upper spectral edge at
  corr −0.86; ×35 scatter). Cap: λ ≥ 0.5 (float headroom, measured).
- ✅ **21 bulk half** (2026-08-16): `thooft-target/condensate.py` — LM
  2409.11324 (6.1)–(6.3) transcribed (exact condensate at any mass;
  Zhitnitsky chiral value −Nc√(g²/12π) reproduced through the I(α)
  near-singularity); GMOR corner vs the Jacobi solver — FIRST α→−1 data:
  2λ₀ → (2/π)√(a/3) with measured NLO ≃ √a (solver 10-digit converged);
  vacuum energy by Feynman–Hellmann: the ENTIRE 2d counterpart of §5's
  fine-tuning drama is the single finite number ε(0) = −0.0407 Nc g² —
  super-renormalizability confronted: nothing to tune. (JT-dS ε = 3/(8S)
  derivation = the still-open theory half.)
- ✅ **24 note** (2026-08-16): `duality/DISCRIMINATOR.md` — the KMS/OTOC
  axis across the three programs, quote-anchored (NV §2.3 engineers
  dS-KMS → faces the CK no-go head-on, discriminating computation = the
  OTOC sign of their dressed operators, unpublished; sine-dilaton §1:
  correlators thermal at the FAKE temperature, black-hole bulk,
  sub-maximal scrambling — the one reading consistent with this repo's
  data as-is (β_fake·λ_L = 2π is definitional given β_fake ≡ β/v; the
  quantitative θ↔v tie to their (1.7) is untested — open); MSS
  flat-space owes the
  sign-flipping map, and the fold is now excluded in-model).
- ✅ **13/24 probe** (2026-08-17): `duality/signflip_probe.py` — the first
  candidate tomperature map that works: the antipodal half-fake-period
  shift (Δ = π/v; NV's ±iβ/4 dressing, geometrically) flips the OTOC sign
  EXACTLY in the closed form (a(Δ) = a(0)e^{−ivΔ}, Re a: −1 → +1, every
  temperature, two configs; shifted 2-pt finite/positive; TOC unchanged)
  and, exact-in-λ, for λ ≤ 1.7 (not λ ≥ 1.8; unshifted control scrambles
  at every λ; truncation-clean). PLUS: the same shift makes the Wightman
  sector EXACTLY thermal at β_eff = β + β_fake (per-line detailed
  balance, any λ) → the tomperature at the T=∞ hologram point — one map,
  both CK dS-observer signatures. Phase period set by v, not the
  measured growth κ ≈ 0.83–0.88 > v: no recalibration rescues strong
  coupling. Necessary, not sufficient — the sign-gate of the tomperature
  map is OPEN semiclassically; the NV dressed-OTOC remains the closer.
  DEEP campaign (--deep): at the hologram point (β=0) the frame
  anti-scrambles in a TRANSIENT window growing ~½ln(1/λ) (t·Jc 0.57→1.14
  over λ 1.0→0.3) — a semiclassical dS EPOCH, scrambling reasserts after;
  MSS evasion located (phase winds at constant |a|: no strip-KMS);
  measured costs: the shifted contour is FOLD-ADJACENT (edge-dominated
  shapes, no Gibbs/sech form, no fake-period recurrence — ratio
  identities survive, shapes don't); rate extraction mode-contaminated.
- ⬜ Open: 18 remainder (the tower extraction itself, now located at fixed
  p, λ→0 — a 1/N computation from the double-scaled viewpoint; with the
  conserved-channel theorem in hand, the sharpest tractable first bite is
  the EXACT finite-N disorder-averaged pair moments of A_n by Wick — a
  3-pairing Isserlis sum over three coupling pairs, exact_wick-style, no
  ED noise — which would exhibit the measured approach to the proven
  N = ∞ zero), 17 proper
  (the O(1/q) residues — theory), 13 remainder (tomperature map: the
  antipodal shift passes the sign gate semiclassically, see above; the
  full dictionary is open), 20, 21 theory half (JT-dS ε = 3/(8S)), 22,
  27 sharpening (ln N onset, beyond ED); NV-sector OTOC sign (the
  rung-24 discriminating computation; needs the two-sided constrained
  chord extension).

- ✅ **The conserved-channel THEOREM** (2026-08-17, `duality/chord_charges.py`;
  upgrades rung 15 proper from measurement to proof): in a sequential pair
  word the fermion pair factors through the zero-chord sector as the
  dressed propagator C(m) = close∘T_Mᵐ∘open, and the three-term recursion
  close_j T_M = q^(Δ+j) T close_j + [j] close_{j−1} + (1−q^(2Δ+j)) close_{j+1}
  (scalar coefficients — two n₁-independence miracles) forces C(m) = P_m(T),
  the Motzkin-path polynomial (flat q^(Δ+j)T, down [j], up 1−q^(2Δ+j)).
  Hence the chord representation of EVERY A_n is a polynomial in T:
  **every spectral moment μ_k (k ≥ 1, both orderings, odd k included) of
  every A_n pair channel vanishes at N = ∞, every λ, every matter
  weight.** Explicit charges: X̂₁ = −(1−q^Δ)T (the A₁ = −2pH image),
  X̂₂ = (1−q^Δ)²T² + (1−q^(2Δ)). The O_Δ contrast is now algebraic:
  q^(Δn̂) does not commute with T, so random matter propagates while the
  H-derived tower cannot. Verified: recursion as a matrix identity
  (j ≤ 5), C(m) ≡ P_m(T) and [C(m),T] = 0 (m ≤ 8), closed forms, and the
  word-level zeros extended to n ≤ 5, k ≤ 6 (|μ_k|/gross ≤ 3e-17) across
  the λ grid and the committed Δψ scan (18 new tests + extended
  chord_moments.json). Scope unchanged: the i≠j pair amplitude at strict
  N = ∞; finite-N and the flavor-diagonal channel stay with the ED bench.
- ✅ **Full-repo review pass** (2026-08-17, three parallel adversarial
  agents + independent verification; all physics claims survived —
  pauli strings bitwise-exact vs a from-scratch construction, Wick
  closed forms 1.000000000000 vs fresh enumeration, chord engine
  unbroken at 6e-16 on fresh words, thooft numerics reproduced to the
  digit). Fixes applied: the antipodal detailed-balance check made
  NON-CIRCULAR (independent contour-engine cross-check, ≤1e-13, stored
  as `engine_check_dev` and asserted; figure panel (d) now plots
  residuals, not an identity); SEM convention unified to
  std(ddof=1)/sqrt(n) (three files had /sqrt(n-1) — committed JSONs
  predate the fix, ~10–15% conservative); sine_solver header identity
  corrected (it stated <u,FP v> = −B, the true identity is
  −B − <u,Vv>; code was always right); jacobi n_ts aliasing guard +
  interleave bounds check added; chord.py walker deduplicated and dead
  code removed; ChordCorrelators now caches its Contour (signflip
  probe 4x faster); README precision claims aligned with committed
  artifacts (thooft result-table digits, LM s=2 vs s≥4 rates, Nc range
  4–10, dense-check 5e-15 with its 1e-9 gate); stale test counts and
  timings corrected everywhere.
- ✅ **Branch landings** (2026-08-17, each agent-reviewed before merge):
  ws/moments-fast (charge-sector-blocked Dirac moments pipeline, 8×+
  measured, dense-checked 5e-15, the Nc = 12 production point),
  ws/wick-m4 (exact E[m4]/Var(m4): relVar(m4) = 4·relVar(m2)(1+o(1)) —
  the second self-averaging anchor; verified by independent re-derivation
  + exact enumeration to 1e-13 + MC), and ws/wick-dirac (exact Dirac-side
  U(N)-emergence rate Var(B_ij), rms ~ Nc^{-(p-1)/2}, MC-verified ≤1.3σ;
  ED-free singlet-moment anchors). All three prior-session worktree
  branches are now merged and deleted.

### Tier 0 — Hygiene (hours; do immediately)

1. **`git init` and commit.** The repo has a `.gitignore` but is not a git repository.
   The README's own "earlier prose had this sign wrong" claim is unverifiable precisely
   because there is no history. Everything below deserves provenance.
2. **Fix the audited README defects** (the five items above: the "not in print" sentence,
   the factor-2 attribution, the misquote, ᾱ = ḡ²N, the ED/combinatorics and "four classes"
   phrasings, and the "validated e^(−a√N)" self-contradiction). The repo's credibility is
   its honesty; these are the only places it currently overclaims.
3. **Actually pin dependencies** (`==`, not `>=` — both READMEs say "pinned"; neither is),
   one `pyproject.toml` so modules can import each other, register the `slow` marker in
   both modules.
4. **Small hardening:** guard `_f()` at k = ±1; assert the even/odd interleaving inside
   `spectrum()` instead of assuming it (it can fail silently toward the chiral limit);
   `ddof=1` in variance estimates; build plot data from `results.txt` programmatically
   instead of hardcoding.

### Tier 1 — Low-hanging fruit (laptop-minutes to a day each; real scientific value)

5. **Encode the Mₙ operator/CP audit as tests** — the cheapest real physics on the table.
   M₀ = Q, M₁ = −(p/2)iH (two-line check, already verified ad hoc during this review, not
   yet in the harness), **build the complex SYK Hamiltonian and finally test [H,Q] = 0**
   (the operators already exist in `syk.py`; today the "QED seed" is untested dead code).
   Pins correction #3 permanently and opens the Dirac sector.
6. **Construct the anti-unitary CP/T operator in ED and *measure* CP(Mₙ) = (−1)ⁿ⁺¹.**
   In 0+1d "CP" is currently a naming convention; a concrete operator (charge conjugation ×
   time reversal in the Dirac model) turns half the spectrum-match claim into data.
7. **Error bars + honest model selection in `syk-self-averaging`:** bootstrap or analytic
   √(2/(n−1)) bars on every relVar and rms_off; AIC over a candidate set that *includes*
   the true double-scaling form e^(−a√N·lnN); state plainly when the data cannot
   discriminate (currently it cannot, for the symmetry-violation fits). Fix the figure to
   plot `rms_off` (the signal the module itself trusts).
8. **Bitwise Pauli-string Hamiltonian assembly.** Each Majorana monomial is one signed
   Pauli string computable by XOR/phase bookkeeping — no sparse matmuls. 10–100× speedup;
   unlocks n_inst ≈ 2000 (3% error bars), sector-resolved ED at N = 20–22, and Krylov
   methods to N ≈ 28–34 later.
9. **External validation of `thooft-target` at α ≠ 0** (its own admitted gap): digitize
   Hanson–Peccei–Prasad or Brower–Spence–Weis massive-quark eigenvalues and assert against
   them. Add more FLZ sum rules (Appendix A lists higher-s certificates) and an
   eigenfunction-level check (wavefunctions/decay constants — the duality comparison will
   need matrix elements, not just eigenvalues).
10. **The conversion-dictionary layer, built stubbed now:** one small module encoding
    n_paper = n_FLZ + 2, the Mₙ normalization, CP conventions, and the λ = p²/N vs
    λ_chord = 2p²/N, q = e^(−λ) factor-2 landmine — with tests, *before* any DSSYK number
    exists to compare.

### Tier 2 — Medium difficulty (days to weeks; each is a self-contained result)

11. **Exact-Wick adjoint variance — test e^(−a√N) where ED can't go.** B_jk is quadratic in
    J, so Var(B_jk) reduces by Wick to a finite combinatorial sum over coupling pairs
    weighted by ±1/0 traces of Majorana strings. Evaluable at arbitrary (N,p) including
    p ~ √N at N in the hundreds — a direct test of the symmetry-violation claim (Eq. 3.28)
    in the regime that motivated it. Same trick gives exact relVar(m₄) as a second
    non-tautological anchor.
12. **Verify the load-bearing combinatorics.** Symbolically enumerate index contractions at
    small (N,p) to check diagram A = N², diagram B = p!/N^(p−2) — the entire emergent-gauge
    story rests on this counting — and audit arXiv:2501.09423 (the companion carrying the
    systematic derivation; never inspected by anyone in this review). Also check
    arXiv:2511.10907 (Susskind's "Correction to a Wrong Claim" about *where* the DSSYK-dS
    entropy sits) against §5's entropy placement — if the horizon bookkeeping moved, the
    Holographic Tuning argument may need re-examination.
13. **Boost-to-mass dictionary calibration.** Solve the free massive Dirac fermion in
    Rindler space at T = ∞ and work out exactly how an invariant mass gap imprints on
    boost-generator spectral functions (which are continuous per mass). Define the DSSYK
    observable that actually encodes a 't Hooft mass *before* hunting for peaks. Every
    spectrum-match plan silently assumes "discrete boost-frequency poles = masses"; this
    rung is what makes rung 16 interpretable.
14. **Single-realization Mₙ spectroscopy by ED** — the cheapest version of "the best test."
    Compute ⟨Mₙ(t)Mₘ(0)⟩ at T = ∞ in single Dirac-SYK instances (self-averaging legitimizes
    single-shot at the measured 1/C(N,p) rate — cite it as such, not as e^(−a√N)) for
    N ≤ ~22 with rung 8's fast assembly. Nobody — paper authors included — has run even
    this. Whatever structure appears (or doesn't) is publishable guidance.
15. **N = ∞ singlet spectral moments vs 't Hooft sum rules.** Disorder covariances and
    moments of the singlet correlator have exact finite-(N,p) chord/Wick expressions —
    ED-free, evaluable at N ~ 10²–10³ along p = √(λN). Compare against `thooft-target`'s
    sum rules (7ζ(3), 2) through the rung-10 dictionary. Robust to the discrete-vs-continuum
    unknown, and the first genuine number-for-number N = ∞ contact between the two sides.
16. **Map symmetry violation in the regime where the match lives:** fixed p, λ → 0 —
    quantify the *power-law* N⁻ᵖ suppression there (not the double-scaling e^(−a√N)),
    connecting the syk module's story to the thooft module's regime. No one has linked the
    two numerically, and the linkage bounds how "exact" the emergent gauge symmetry is
    where the spectra are to be compared.

### Tier 3 — Hard (weeks to months; the real next results)

17. **Large-q anchor + Regge continuation.** The closed-form large-q 4-point function
    (Streicher; Choi–Mezei–Sárosi) at T = ∞, Sommerfeld–Watson/Regge-continued, as the
    analytic λ → 0 anchor — the one limit where the match should become exact. Validates
    any numerical pipeline before it's trusted, and may already reveal the tower.
18. **The headline: extract the singlet/bilinear spectrum from the *published* exact DSSYK
    4-point function** (Berkooz et al. 1811.02584, chord transfer-matrix) at T_B = ∞,
    small λ. First deliverable is binary and publishable either way: **discrete tower or
    continuum?** If discrete: masses vs the FLZ/`thooft-target` yardstick through the
    dictionary. Falsifiable fingerprints, in order of sharpness: (a) a *single* Regge-like
    tower (any daughter trajectory falsifies), (b) alternating CP, (c) absence of dynamical
    n = 0,1 states, (d) the interleave splittings already sitting in
    `spectrum_double_precision.csv` (2λ adjacent even/odd gaps 1.0167, 1.0029, 1.0012, … —
    a parameter-free fingerprint that also forces a decision on what the paper's "two
    degenerate trajectories" wording means: the CSV shows an interleaved, *non*-degenerate
    tower), (e) the Regge slope τ = p²m_q². Note the race: the authors have promised this
    computation.
19. ✅ **Dirac-SYK QED campaign:** charged-sector correlators for the linear-Coulomb /
    U(1)-confinement signature (the 2012.12326 story), plus U(N)-class self-averaging —
    closing the Majorana-numerics/Dirac-claims gap wholesale. *(Done — progress log.)*

### Tier 4 — Hardest (open-ended; each could be a paper)

20. **The order-of-limits problem in full.** T_B → ∞, N → ∞, λ fixed-vs-→0, and the
    boost-spectrum → invariant-mass reorganization (rung 13's dictionary made rigorous).
    The chaos-to-integrability transition at small λ (2403.01950, 2601.09801) is likely
    the mechanism zone. This is where the duality most plausibly *fails*, and a sharp
    statement of where it holds would be a real contribution.
21. **Holographic Tuning made falsifiable in D=2.** Derive the JT-dS analog of ε = 3/(8S);
    compute the 't Hooft-model vacuum energy under the dictionary (the chiral condensate is
    known exactly at large N) and test whether interactions shift the dimensionless vacuum
    energy. Confront honestly the fact that QCD₂ is super-renormalizable — whether the toy
    can even exhibit the disease it claims to cure. First bulk-side test of §5.
22. **The finite-N coupling curve ḡ²(N).** Fig. 17 is a sketch with no formula, yet the
    whole set-it-and-forget-it inversion rests on it. Compute it: match a bulk observable
    (e.g. a mass ratio from rung 18's pipeline) at a sequence of finite N and extrapolate.
23. ✅ **Baryon-sector first contact.** ε-tensor baryon operators in small-N Dirac SYK
    (N = 4–8): mass scaling ~N, baryon-baryon correlations — the only route to data on the
    wholly-underived "inter-baryon forces are gravitational / O(√N) mesons" claim.
    *(Done at first-contact level — progress log.)*
24. **Discriminate the competing duals.** Whatever rung 18 produces, run the same
    observable against the Narovlansky–Verlinde and sine-dilaton predictions where they
    differ. The meson tower is currently the only proposed observable sharp enough to
    separate the three programs. *(2026-07-23: the CK anti-scrambling/KMS constraint is
    a second discriminator axis — check which programs formulate the hologram as a KMS
    state at the dS temperature.)*

### 2026-07-23 additions (the CK/HZ response ladder)

25. ✅ **The OTOC sign test.** Does DSSYK∞ scramble or anti-scramble under the flat
    dictionary? *(Done — scrambles, Re a = −1 exactly at all temperatures;
    `duality/ANTISCRAMBLING.md`.)*
26. **The Euclidean-fold laboratory.** HZ's fold prescription evaluated exactly in
    finite-N SYK vs the continued closed forms. *(First data done — mechanics work,
    Lorentzian response tracks. Decider #1 done 2026-08-12: the Euclidean drift is a
    genuine N-obstruction, not a λ-artifact — the folded connected part is not
    1/N-suppressed; progress log.)* Remaining decider: the folded 4-pt in the exact
    chord theory at N = ∞ (the fold is a contour statement; chord correlators are
    analytic in the chord angles) — adjacent to the rung-18 toolchain and would answer
    HZ exactly.
27. ✅ **Freeness onset.** CK show the early/late-operator algebra becomes a free product
    beyond the scrambling time. Freeness is a concrete statement about mixed moments of
    a(0), b(T) — directly measurable with the existing Pauli-string ED at T_B = ∞.
    Nobody has exhibited the free-product transition in a microscopic dS-candidate
    model; self-contained, publishable either way. *(Done at first-data level —
    progress log: arcsine-law convergence with a ~1/dim floor. Remaining
    sharpening: the ln N onset-time question needs sizes beyond ED.)*

### Blocked / not worth attempting from here

- Proving DSSYK∞ = JT-dS (assumed by the authors themselves).
- The 1+1 → 3+1 uplift — contrary to the paper's framing in secondary discussion, *no*
  uplift construction exists in the source; §§5–6 are explicit that this is a hope.
- Anything requiring instance-level ED in the true double-scaling regime (p ~ √N, N ≳ 100,
  dim 2⁵⁰) — the honest wall. Rungs 11 and 15 are the ways around it.
