# Draft note: the conserved-channel theorem

`main.tex` is a ~4-page standalone proof that the **two-flavor (i ≠ j)**
part of the singlet bilinear channels A_n = Σψ(ad_H)ⁿψ cannot propagate at
N = ∞ in double-scaled SYK: its chord representation is a polynomial in the
chord Hamiltonian, via the dressed-propagator recursion of
`duality/chord_charges.py`.

Scope is the point of the note, so it is stated in the abstract: the full
channel also carries a flavor-diagonal part of the same order, which does
*not* vanish for n ≥ 2 (the full A₂ width rises to ≈ 1.19 along
λ_chord = 2, while the n = 1 parts cancel to 3e-16 — that is A₁ = −2pH).
The result is therefore a statement about flavor structure: whatever weight
the tower of arXiv:2607.05678 carries at N = ∞ is entirely diagonal.

**Status: draft for the author's review — not for distribution.** Resolve
the two `TODO(author)` markers (affiliation, acknowledgments) before
circulating.

This machine has no TeX toolchain. Compile on Overleaf or with
`pdflatex main.tex` (twice, for references). Every displayed equation is
machine-verified in `duality/test_chord_charges.py`.
