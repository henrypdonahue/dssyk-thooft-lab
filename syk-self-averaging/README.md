# DSSYK self-averaging, the Mn tower, and the QED sector

Tests of the falsifiable corner of the paper: the §3.5 self-averaging
claims, plus the operator content its U(N)/Dirac "standard model" story
needs. Model (Eq. 3.2): H = i^(p/2) Σ J ψψ...ψ with Var(J) = p!/N^(p−1),
λ = p²/N, ⟨W⟩ = Tr(W)/dim.

| claim | eq. | verdict |
|---|---|---|
| singlet Var⟨W⟩ ~ p!/N^p | 3.27 | confirmed; rate exactly 2/C(N,p) |
| adjoint violations ≤ exp(−a√N) | 3.28 | verified exactly in double scaling (`exact_wick.py`) |

## Results

**Eq. 3.27.** The singlet rate is exactly 2/C(N,p) — one over the number
of couplings. It equals the paper's p!/N^p times an O(1) factor. Caveats:
the m₂ identity holds regardless of the paper's diagrammatics; the
non-tautological content is m₄ and m₆ self-averaging at the same rate
(confirmed by ED at N ≤ 18, within ~10% statistics).

**Eq. 3.28 — settled beyond ED.** B_jk = Tr(ψⱼHψₖH)/dim is quadratic in
the couplings, so Wick gives closed forms at any (N, p):

    Var(B_jk) = 4σ⁴ C(N−2, p−1),    E[B_jj] = σ² C(N,p) (1 − 2p/N),

verified against literal contraction enumeration (exact) and ED (~1%).
Along p = √N: ln rms = −¼ √N ln N — beats exp(−a√N) for every a, scanned
to N = 1600. At fixed p (the match regime, λ → 0) emergence is only
**polynomial**: violation amplitude N^(−(p−1)/2), per-element violation
N^(−(p+1)/2). Reaching 10⁻⁶ needs N ≈ 826 at p = 4: the emergent symmetry
is weakest exactly where the spectrum match is to be made. E[B_jj] = 0 at
p = N/2 derives the ratio artifact flagged below.

**Mn tower and CP, measured** (`dirac.py`, `mn_spectroscopy.py`). Exact
and asserted: [H, Q] = 0; M₀ = Q; M₁ = −(p/2)iH (the paper's M₁ = H fails
at O(1)); M₂ = −Σ ċ†ċ (the paper's + sign is wrong); Mn† = (−1)ⁿMn only
for n ≤ 2. CP must be the **unitary** particle-hole map (M₁ ∝ iH rules
out antiunitary). CP(M₀) = −1 and CP(M₁) = +1 are exact; raw Mn for n ≥ 2
are CP-mixtures (15–65% wrong channel), so CP = (−1)^(n+1) lives in
CP-resolved correlators. `mn_spectroscopy.py` computes them — first data:
conserved n = 0, 1 exact at ω = 0; broad continua at λ = 1.6.

**U(N)/QED campaign** (`qed_campaign.py`): U(N)-class self-averaging made
exact (Dirac analog of the Wick closed forms, verified by enumeration and
ED); disorder-averaged charging curves with full covariance and GLS fits;
charged-vs-singlet spectral scales. The capacitive-vs-confining question
is stated with ΔAIC, not asserted.

**Baryon sector, first contact** (`baryons.py`): epsilon-tensor baryons at
Nc = 4–8 — exact 2-point collapse onto charge sectors, mass extensivity,
and the edge-curvature hierarchy. First data on the paper's underived
inter-baryon-force claim; hierarchy claims kept honest.

**Shape-invariant moments** (`moments_pipeline.py`): CP-resolved spectral
moment invariants (w₂, kurtosis) along fixed p = 4 — map-robust targets.
**Large-q ED bench** (`largeq_bench.py`): ground truth for
`duality/largeq_anchor.py`.

**Freeness onset, first data** (`freeness.py`): Cui–Kolchmeyer argue
early/late observer algebras combine as a *free product* past the
scrambling time. For the pair E = iψ₀ψ₁, F(t) = iψ₂(t)ψ₃(t), freeness
means spec(E + F(t)) = the arcsine law. Measured at N = 12–24: the
spectrum lands on the arcsine law, with the W₁ plateau falling like
dim^(−0.9) (a rigid-spectrum floor, far below sampling noise), and all
alternating words → 0. Caveats: the crossover sits at the *dissipation*
time (t𝒥 ≈ 1.2, nearly N-independent) — ED cannot resolve the ln N
scrambling-time separation where the claim properly lives; word plateaus
at N ≥ 20 are statistics-limited. First exhibition of the free-product
transition in a microscopic dS-candidate model.

**Majorana singlet-moment bench** (`mn_majorana_bench.py`): the ED side
of the exact chord computation in `duality/chord_moments.py`. The raw
width w₂ rises along fixed p toward the λ → 0 corner — the tower channel
turning on. Scope, per adversarial review: the chord result is an i ≠ j
pair amplitude, and at dense-ED sizes the measured channel is
flavor-diagonal-dominated, so the bench does not yet probe the chord
N = ∞ zero (at matched λ the raw w₂ rises with N — the honest negative).
The (8,4) point is an exact per-instance conservation special to (8,4);
mechanism unknown.

Numbers: `results.txt` / `results.json`. Figures and data JSONs live next
to their scripts.

## Why the numbers are trustworthy

- Instance-level ED caps near N = 20 and cannot reach double scaling —
  every double-scaling statement here is exact combinatorics, not
  extrapolation. ED-window fits carry error bars and ΔAIC verdicts.
- Construction validated before physics: Clifford algebra and parity to
  machine precision; the N mod 8 RMT pattern (GUE/GSE/GOE) reproduced
  with parity-sector resolution.
- Fast bitwise builders are equivalence-tested against independent sparse
  constructions. Read `rms_off`, not the p = N/2-artifacted ratio.

## Files

| file | purpose |
|------|---------|
| `syk.py` | sparse Jordan–Wigner operators (slow reference builders) |
| `pauli_strings.py` | bitwise Hamiltonian assembly (10–100×, equivalence-tested) |
| `dirac.py` | complex SYK: [H,Q] = 0, the Mn tower, unitary CP machinery |
| `exact_wick.py` | closed-form Eq. 3.28 statistics, beyond ED |
| `mn_spectroscopy.py` | CP-resolved Mn spectral functions |
| `moments_pipeline.py` | CP-resolved moment invariants with statistics |
| `moments_fast.py` | charge-sector-blocked reimplementation (8x+ faster, dense-checked to 1e-15): the Nc = 12 production point and the fixed-p trend figure |
| `qed_campaign.py` | U(N) self-averaging + charging curves |
| `baryons.py` | baryon sector: collapse, extensivity, hierarchy |
| `largeq_bench.py` | ED ground truth for the large-q anchor |
| `fold_bench.py` | folded/complex-time 4-pt ED lab; `--scan` = the fixed-λ decider |
| `fold_edge_probe.py` | the fold's spectral-edge mechanism, measured |
| `freeness.py` | the classical → free-product transition |
| `mn_majorana_bench.py` | Majorana singlet-channel moments (ED side of the chord computation) |
| `validate_syk.py`, `self_averaging.py` | RMT validation; the §3.5 measurement |
| `plot_self_averaging.py`, `plot_freeness.py` | figures (read the JSONs) |
| `test_*.py` | asserting suites |

## Usage

```bash
pip install -r requirements.txt
python3 validate_syk.py         # RMT validation (~7 s)
python3 self_averaging.py       # the 3.27 measurement (~2 min)
python3 exact_wick.py           # exact 3.28 + double-scaling scan (~1 s)
python3 dirac.py                # Dirac/QED identities + CP report (~5 s)
python3 mn_spectroscopy.py      # Mn spectral functions (seconds)
python3 qed_campaign.py         # QED campaign
python3 baryons.py              # baryon campaign
python3 freeness.py             # freeness campaign (hours)
python3 mn_majorana_bench.py    # singlet-moment bench (~1 h)
pytest -q -m 'not slow'         # fast asserts (~1 min)
```

## Scope

This settles §3.5 (both equations, one exactly) and the operator-level
prerequisites of the U(N) story. The headline spectrum match itself lives
at fixed p, λ → 0 (see `duality/CHORD.md` for why), with the CP-resolved
spectroscopy as its pipeline.
