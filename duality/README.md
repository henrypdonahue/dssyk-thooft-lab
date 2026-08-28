# The dictionary, and the N = ∞ engines

Before you can compare a DSSYK number against a 't Hooft number, you have
to fix conventions. An index offset, a factor of two, a normalization.
Get one wrong and the match is meaningless.

So every convention lives here, once, with a test. The dictionary was
written **before** the DSSYK side existed. That way a mistake has exactly
one place to hide.

| encoded fact | source |
|---|---|
| n_paper = n_FLZ + 2, the photon/graviton relabeling | correction #1 |
| CP = (−1)^(n+1), consistent with FLZ parity through the offset | paper §4.2–4.3 |
| M₁ = −(p/2) iH (`m1_coefficient`) | correction #3, tested in `syk-self-averaging` |
| paper λ = p²/N, chord λ = 2p²/N, q = exp(−λ_chord) | the factor-2 landmine |
| ᾱ = ḡ²N = p², ḡ² = λ, τ = p² mq² | paper §4.1–4.2 |

## What a DSSYK tower would have to reproduce

Two fingerprints. Both are parameter-free, and both come from
`thooft-target/reference_spectrum.json`.

1. **Mass-squared ratios** Mn²/Mη²: 1, 2.3794, 3.7285, 5.0892, …
2. **Interleave splittings.** One alternating-CP trajectory, with gaps
   near unity: 1.0167, 0.9944, 1.0029, … → 1.

Degenerate CP doublets would falsify the identification. This also
settles the paper's "two degenerate Regge trajectories" wording. The
slopes are equal and the trajectories interleave. Adjacent levels are
not degenerate.

## The rest of the module

**[CHORD.md](CHORD.md)** and `chord.py`, `chord_moments.py`,
`chord_charges.py`, `chord_fold.py`. The N = ∞ chord engine: a transfer
matrix plus matter-chord contours, held up four independent ways
(brute-force enumerator, published closed forms, Krawtchouk microscopy,
finite-N ED to 0.6–1.6%). The spectra have no structure at any λ. The
i ≠ j pair part of every singlet channel is frozen at N = ∞, and that is
a theorem: its chord image is a polynomial in the chord Hamiltonian. The
diagonal part is the same order and does not vanish for n ≥ 2, so the
full channel propagates. CHORD.md carries the measured widths. The
Euclidean fold fails in the double-scaled theory itself. Figure:
`plot_chord.py`.

`diagonal_chord.py`. The diagonal channel at N = ∞. An H-chord weighs
q^Δψ if and only if it crosses the four-mark set an odd number of times.
That law is derived from the finite-N sign law, and it is triple-gated.
The diagonal A₂ width is w₂ = 4.67(9) at λ_chord = 2, Δψ → 0.

**[ANTISCRAMBLING.md](ANTISCRAMBLING.md)** and `antiscrambling.py`. The
OTOC sign test, and the fold. DSSYK∞ scrambles at every temperature,
a(v) = −1/cos(πv/2) < 0. So the flat dictionary cannot give a de Sitter
observer's anti-scrambling. The Harlow–Zhao fold does not repair it. It
fails at ED sizes (`../syk-self-averaging/fold_bench.py --scan`) and it
fails at N = ∞ (`chord_fold.py`).

**[DISCRIMINATOR.md](DISCRIMINATOR.md)**. Where the three competing duals
sit: MSS flat-space, Narovlansky–Verlinde, and sine-dilaton. Measured
against the chaos-bound no-go and this repo's data. Sine-dilaton's
black-hole reading agrees with the data qualitatively. Its quantitative
θ ↔ v tie is untested here.

`signflip_probe.py`. The first candidate dS map that works. Shift the
probe pair by half the fake period. That flips the OTOC sign, exactly in
the closed form and in the exact theory for λ ≤ 1.7, and it makes the
2-point sector thermal at the tomperature. Necessary, not sufficient.
Mode separation is in `antipodal_modes.py`. Figure: `plot_signflip.py`.

**[NV-DRESSING.md](NV-DRESSING.md)** and `test_nv_dressing.py`. The
identification note. NV's ±iβ_dS/4 dressing
([2310.16994](https://arxiv.org/abs/2310.16994) Eqs. 33–34, 61–62) is
the fake-period contour shift, exactly, at the contour and 2-point
levels. The composition algebra and the thermality are tested. What
remains is the constrained-sector step.

**[BOOST.md](BOOST.md)** and `boost_mass.py`. Boost-to-mass calibration.
One mass gives a boost-frequency continuum with nodes, not a peak. So
peak-hunting is the wrong observable. The tomperature map is the stated
open question.

`boost_forward.py`. The forward tower check: 't Hooft masses and decay
constants pushed through the boost kernel. The raw moments diverge
(S → 1/ω). The regulated invariants sweep the whole SYK band as the
regulator turns. The coupling-weighted tower sits within 5% of a single
mass at the invariant level. So the tower discriminant is the nodes, not
the moments.

`boost_inverse.py`. The inverse direction, measured at certified
precision. The weights-from-moments map saturates at one moment per
level and is super-geometrically ill-conditioned (each level costs
~1–2 more digits). Moment invariants carry regulator leverage ≈ 1, so
they carry no mass digits until the state map fixes ε; node positions
carry leverage exactly zero. The node verdict, with its price list.

**[LARGEQ.md](LARGEQ.md)** and `largeq_anchor.py`. The large-q, λ → 0
anchor, with the conventions bridge 𝒥 = √2 p locked against ED. The
meson tower cancels at leading order. The anchor lives at O(1/q).
`largeq_kernel.py` derives the 𝓜ᵉ/𝓜ᵒ skeleton from the
Maldacena–Stanford kernel eigenproblem and reproduces the transcribed
sets to 3e-12 — the perturbable substrate for the O(1/q) task.

`jt_ds_vacuum.py`. The JT-dS analog of ε = 3/(8S), sympy-certified:
the general-D identity with c₄ = 3/8, the fourfold D → 2 degeneration,
the exact first law δS = −2πℓ²ρ, and the reduction-inherited family
whose exponent is a dictionary input.

`dirac_chord.py`. The complex-SYK chord layer (Berkooz–Narovlansky–Raj
[2006.13983](https://arxiv.org/abs/2006.13983)). The chord sum is the
Majorana one at q = e^(−p²/Nc) = e^(−λ_paper), so the factor-2 landmine
is Majorana-only. It is ED-verified. M₁ = −(p/2)iH is an exactly
conserved transfer-matrix insertion. Charged probes carry Δ = (a+b)/p at
fixed charge q(𝒬). Chord rules for the H-correlated M_{n≥2} tower are
unpublished; the pair half never needed them — the bare-fermion census
(lone crossing exactly 1 − p/Nc, the Majorana sign census at N = 2Nc)
transfers the conserved-channel theorem verbatim to (q_dirac, Δ = 1/p),
freezing every Mₙ pair channel at Nc = ∞. The diagonal half and the CP
resolution are pinned as the remaining gap.

```bash
python3 dictionary.py            # the yardstick in both conventions
python3 antiscrambling.py        # the OTOC sign table and the fold map
python3 chord.py --spectra       # the lambda -> 0 spectral scan  (~4 min)
python3 chord_moments.py --run   # exact N = inf moments          (~5 s)
python3 chord_charges.py         # the conserved-channel theorem  (~1 s)
python3 chord_fold.py            # the fold at N = inf            (~10 s)
python3 signflip_probe.py        # the antipodal map              (~15 s)
python3 signflip_probe.py --deep # its deep campaign              (~15 s)
python3 jt_ds_vacuum.py          # the JT-dS vacuum-energy note    (~5 s)
python3 largeq_kernel.py         # the derived large-q skeleton    (~5 min)
python3 boost_inverse.py         # moments -> spectrum conditioning (~3 min)
pytest -q                        # 284 tests
```
