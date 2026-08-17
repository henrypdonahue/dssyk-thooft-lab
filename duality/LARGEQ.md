# The large-q anchor: what the λ → 0 corner actually contains

`largeq_anchor.py` encodes the equation-verified large-q SYK closed forms
(Maldacena–Stanford 1604.07818; Streicher 1911.10171; Choi–Mezei–Sárosi
1912.00004; Gross–Rosenhaus 1702.08016) with the exact conventions bridge
to this repo (𝒥 = √2 p), validated against the ED bench
(`syk-self-averaging/largeq_bench.json`).

## The structural finding

The large-q singlet 4-point function is known in closed form at all
temperatures. Its Sommerfeld–Watson structure has poles only at m = 0
(Hamiltonian exchange — the entire time-ordered part) and m = ±v (the
scramblon, OTOC only). The would-be bilinear tower — the quantization
sets 𝓜ᵉ/𝓜ᵒ (CMS 3.3) — **cancels identically at leading order**, and its
couplings are O(1/q²) (GR 2.15).

**Consequence: the λ → 0 anchor of the 't Hooft match cannot come from
the leading large-q 4-point function.** The tower first appears at
O(1/q), organized on the exact-in-temperature skeleton 𝓜ᵉ/𝓜ᵒ, whose
finite-temperature residues are unpublished. Computing them is the open
task — and it is where the quantitative content of the λ → 0 corner
lives.

## Encouraging structure

- 𝓜ᵉ and 𝓜ᵒ **interleave** — even/odd alternation, exactly the
  CP-alternation slot of the 't Hooft tower (`dictionary.py`).
- At T = ∞ the sets degenerate to integers in 2π/β units — consistent
  with `BOOST.md`: masses are not boost-frequency peak positions.
- The ED bench confirms the transcription end to end: ⟨H²⟩ matches the
  bridge exactly, and G(t) = sech(√2 p t)^(2/p) matches ED with no free
  parameters at the O(1/p) level (measured deviations 0.07–0.10 at p = 4
  against 1/p = 0.25; the deviation plateaus rather than vanishing as
  λ → 0 — the intrinsic 1/p error survives).

## En-route finding

The "obvious" Majorana singlet (i/N) Σ ψ ψ̇ is exactly conserved:
Σ ψ[H, ψ] = −2pH — the Majorana counterpart of M₁ = −(p/2)iH, asserted in
the tests. The first dynamical Majorana member is Σ ψ̇ψ̇, used by the
bench.
