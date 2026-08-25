# Draft note: the conserved-channel theorem

`main.tex` is a ~4-page standalone proof that the singlet bilinear channels
A_n = Σψ(ad_H)ⁿψ cannot propagate at N = ∞ in double-scaled SYK. Their
chord representation is a polynomial in the chord Hamiltonian (the
dressed-propagator recursion of `duality/chord_charges.py`). Consequence:
the meson tower of arXiv:2607.05678 is a 1/N effect. The spectrum match
lives at fixed p, λ → 0.

**Status: draft for the author's review — not for distribution.** Resolve
the two `TODO(author)` markers (affiliation, acknowledgments) before
circulating.

This machine has no TeX toolchain. Compile on Overleaf or with
`pdflatex main.tex` (twice, for references). Every displayed equation is
machine-verified in `duality/test_chord_charges.py`.
