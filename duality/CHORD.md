# The exact N = ∞ chord engine, and what it says about the tower

*Rungs 18 (tooling + first deliverable), 15 proper, and 26 decider #2
(2026-08-16). Code: `chord.py` (engine + spectral scan),
`chord_moments.py` (singlet-channel moments), `chord_fold.py` (the HZ fold
at N = ∞), ED side `syk-self-averaging/mn_majorana_bench.py`. Verdicts
asserted from the committed JSONs in `test_chord*.py`.*

## The engine (rung-18 toolchain, now in place)

Double-scaled SYK at N = ∞, λ_chord = 2p²/N fixed, is exactly solvable by
chord combinatorics (Berkooz–Isachenkov–Narovlansky–Torrents 1811.02584;
conventions from the review 2407.09396, both fetched and equation-checked).
`chord.py` implements the whole solution as finite sparse linear algebra —
transfer matrix `T|n⟩ = |n+1⟩ + [n]_q|n−1⟩`, matter chords with weight
q^Δ, one- and two-matter stack sectors evaluating **any operator ordering
on any complex-time contour** (crossed/OTO and Euclidean-folded included),
with no Askey–Wilson special functions anywhere.

Held up by three independent legs, all in `test_chord.py`:

1. a brute-force chord-diagram enumerator (explicit sum over pairings with
   crossing weights — shares nothing with the transfer-matrix bookkeeping):
   engine == enumerator to machine precision on every supported topology,
   fermionic signs included;
2. the published closed forms: the θ-integral 2-pt (review Eq. 2.30) via an
   independent q-Pochhammer route, μ₂ = 2(1−q^Δ), Touchard–Riordan moments,
   the q-Hermite measure's moment closures;
3. microscopy with no chords at all: the exact finite-N crossing average
   (Krawtchouk sum) approaches q = e^{−2p²/N} at the measured O(p⁴/N²)
   rate — the chord weight's finite-N quality quantified.

Units bridge (exact, tested): σ_H(N,p)² = C(N,p)p!/N^{p−1}, so
β_chord = βJ·e^{−λ/8}/√λ at N = ∞. Against finite-N ED
(`fold_bench.json`, p = 6, N = 12–16, βJ = 2) the chord thermal 2-pt agrees
to **0.6–1.6%** — compare the ~30–55% floor of the large-q closed forms at
the same sizes. The engine is the repo's best finite-N anchor by an order
of magnitude, and it is exact at N = ∞.

Honest caps, stated: the contour evaluator supports ≤ 2 simultaneous open
matter chords (all 2-pt/4-pt orderings); folded contours below λ ≈ 0.5
exceed float cancellation headroom (measured, see `chord_fold.py`).

## Rung 18 first deliverable: discrete tower or continuum? — **Continuum. And, sharper: the tower channels are conserved.**

`chord.py --spectra` (`chord_spectra.json`): the exact smooth spectral
functions S_Δ(ω) (reduced θ-integral, a continuous function — no binning)
for Δ ∈ {1, ½, ¼} and the Δ→0 fermion-bilinear/length channel
∂G/∂Δ|₀ = −λ⟨n(t)⟩, at λ = 2.0 → 0.05, each with a fine ω-window at the
semiclassical scale 8·Jc where a 't Hooft tower would sit under the flat
dictionary (ω_chord = (ω/𝒥)·Jc, Jc = √λ e^{λ/8}).

- **Zero interior maxima, zero nodes, every channel, every λ.** No
  resonance structure emerges as λ → 0; the exact spectral function
  converges to the smooth semiclassical sech^{2Δ} transform **linearly in
  λ** (max deviation 0.375 → 0.006 over λ = 2 → 0.05, ≈ 0.13·λ). The
  rung-13 alternative (mass-encoding *nodes*) is absent too.
- The chord-number growth ⟨n(t)⟩ — the de Sitter length — approaches
  −(2/λ)ln sech(Jc t) at the same linear-in-λ rate (dev 0.31 → 0.012):
  the exact-in-λ length-growth curve, new data.

**Rung 15 proper** (`chord_moments.py`, `chord_moments.json`) makes the
verdict algebraic. The exact N = ∞ spectral moments of the *actual paper
channels* — the H-derived singlet bilinears A_n = Σᵢψᵢ(ad_H)ⁿψᵢ, assembled
from mixed chord moments by a binomial expansion that is verified as a
per-instance **operator identity** against dense ED algebra
(`test_mn_majorana.py`) — give

> **μ₂ = μ₄ = 0 to machine precision, for n = 1, 2, 3, at every λ from 4.0
> to 0.1 and any matter weight.** The n = 1 case is the exact conservation
> law A₁ = −2pH; the n ≥ 2 cancellations are the same phenomenon: at
> strict N = ∞, fixed λ, **the singlet bilinear tower carries no
> propagating spectral weight at all.** Random O_Δ operators *do* propagate
> (μ₂ = 2(1−q^Δ), kurtosis → 3 + 1/Δ semiclassically — closed forms,
> tested); the H-correlated bilinears the duality actually needs do not.

This is the exact-in-λ, all-orders-in-the-chord-limit version of the
rung-17 finding (the tower cancels at leading order in the large-q 4-pt).
The Majorana ED bench (`mn_majorana_bench.json`, p = 4 and 6, N up to 22)
shows the finite-N side of the same coin:

- along **fixed p with λ falling** — the 't Hooft-corner trajectory — the
  raw channel width w₂ *rises* (n = 2: 0.30 → 1.20 over N = 10 → 20 at
  p = 4): the tower channel **turns on** exactly where the match must be
  made, as a finite-N effect;
- at **fixed λ with N growing**, after stripping the known p = N/2
  proximity factor (1−2p/N)² (the exact_wick artifact; the (p=4, N=8) row
  sits exactly on its zero and is kept as the exhibition), the corrected
  w₂ falls monotonically with N in both p-series and across the
  near-matched-λ anchor (p=4, N=10, λ=3.20) ↔ (p=6, N=22, λ=3.27) —
  consistent with the chord zero. Firm asymptotics need N ≫ 2p, beyond
  dense ED: stated cap.

**Consequence for the headline program.** At T_B = ∞ the binary question
"discrete tower or continuum?" has answer: *at N = ∞ the singlet channels
are neither — they are conserved densities with all spectral weight at
ω = 0; at any finite N they are smooth continua whose weight grows toward
the fixed-p, λ→0 corner.* The meson tower, if the duality produces one, is
a **1/N effect invisible at every order of the strict double-scaled
expansion** — it must be extracted along fixed p, λ → 0 (where 't Hooft
mesons are themselves the leading large-N states), not from the N = ∞
chord theory. This sharpens Part I's "emergence is weakest where the match
is made": the match lives strictly outside the double-scaled description.
The λ→0 anchor (rung 17's O(1/q) residues) and rung 22's ḡ²(N) curve
inherit this as their working regime; nothing here falsifies the duality —
it locates it.

Caveats, in one place: the chord evaluation keeps the matter weight at the
finite-p value q^{1/p} (the standard leading-chord treatment) — the
conserved-channel zeros hold for *any* weight, so they are not an artifact
of that choice. The Dirac CP-resolved groundwork (`moments.json`) is kept
for the eventual charged-chord comparison; the channel identification
Dirac ↔ Majorana is not 1:1 and awaits the complex-SYK chord machinery.
