# Draft note: the conserved-channel theorem

`main.tex` — a ~4-page standalone proof that the singlet bilinear channels
A_n = Σψ(ad_H)ⁿψ cannot propagate at N = ∞ in double-scaled SYK: their
chord representation is a polynomial in the chord Hamiltonian (the
dressed-propagator recursion of `duality/chord_charges.py`). Consequence:
the meson tower of arXiv:2607.05678 is a 1/N effect; the spectrum match
lives at fixed p, λ → 0.

**Status: draft for the author's review — not for distribution.** Two
`TODO(author)` markers (affiliation, acknowledgments) must be resolved
before circulating.

No TeX toolchain on this machine: compile on Overleaf or with
`pdflatex main.tex` (twice, for references). Every displayed equation is
machine-verified in `duality/test_chord_charges.py`.
