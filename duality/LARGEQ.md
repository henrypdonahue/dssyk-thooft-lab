# The large-q anchor: what is really in the λ → 0 corner

`largeq_anchor.py` encodes the large-q SYK closed forms,
equation-checked against their sources (Maldacena–Stanford
[1604.07818](https://arxiv.org/abs/1604.07818); Streicher
[1911.10171](https://arxiv.org/abs/1911.10171); Choi–Mezei–Sárosi
[1912.00004](https://arxiv.org/abs/1912.00004); Gross–Rosenhaus
[1702.08016](https://arxiv.org/abs/1702.08016)). It applies this repo's
conventions bridge, 𝒥 = √2 p. The forms are checked against the ED bench
in `syk-self-averaging/largeq_bench.json`.

## The structural finding

The large-q singlet 4-point function is known in closed form at all
temperatures. Look at its Sommerfeld–Watson structure and you find poles
in only two places. At m = 0, which is Hamiltonian exchange and accounts
for the entire time-ordered part. And at m = ±v, which is the scramblon,
and appears in the OTOC only.

The would-be bilinear tower is the quantization sets 𝓜ᵉ/𝓜ᵒ (CMS 3.3).
It **cancels identically at leading order**. Its couplings are O(1/q²)
(GR 2.15).

**What follows.** The λ → 0 anchor of the 't Hooft match cannot come from
the leading large-q 4-point function. The tower first appears at O(1/q),
organized on the exact-in-temperature skeleton 𝓜ᵉ/𝓜ᵒ, whose
finite-temperature residues are unpublished. Computing them is the open
task. That is where the quantitative content of the λ → 0 corner lives.

## Encouraging structure

𝓜ᵉ and 𝓜ᵒ **interleave**. That even/odd alternation is exactly the
CP-alternation slot of the 't Hooft tower (`dictionary.py`).

At T = ∞ the sets degenerate to integers in 2π/β units. That is
consistent with [BOOST.md](BOOST.md): masses are not boost-frequency peak
positions.

The ED bench checks the transcription end to end. ⟨H²⟩ matches the
bridge exactly. G(t) = sech(√2 p t)^(2/p) matches ED with no free
parameters at the O(1/p) level. The measured deviations are 0.065–0.095
at p = 4, against 1/p = 0.25. They do not shrink as λ → 0. The intrinsic
1/p error survives.

## One finding along the way

The obvious Majorana singlet (i/N) Σ ψ ψ̇ is exactly conserved. The
identity is Σ ψ[H, ψ] = −2pH. This is the Majorana counterpart of
M₁ = −(p/2)iH, and it is asserted in the tests.

So the first dynamical Majorana member is Σ ψ̇ψ̇, which is what the bench
uses.
