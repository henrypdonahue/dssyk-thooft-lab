# Draft note: the conserved-channel theorem

`main.tex` is a standalone proof, about four pages. It proves that the
**pair (i ≠ j)** part of the singlet bilinear channels Aₙ = Σψ(ad_H)ⁿψ
cannot propagate at N = ∞ in double-scaled SYK. Its chord representation
is a polynomial in the chord Hamiltonian, through the dressed-propagator
recursion of `duality/chord_charges.py`.

Scope is the point of the note, so it sits in the abstract. The full
channel also carries a flavor-diagonal part of the same order, and that
part does *not* vanish for n ≥ 2. The full A₂ width rises to ≈ 1.19 along
λ_chord = 2, while the n = 1 parts cancel to 3e-16, which is A₁ = −2pH.

So the result is a statement about flavor structure. Whatever weight the
tower of arXiv:2607.05678 carries at N = ∞ is entirely diagonal.

**Status: a draft for the author's review. Not for distribution.**
Resolve the two `TODO(author)` markers, affiliation and
acknowledgments, before circulating.

This machine has no TeX toolchain. Compile on Overleaf, or with
`pdflatex main.tex` run twice for the references. Every displayed
equation is machine-verified in `duality/test_chord_charges.py`.
