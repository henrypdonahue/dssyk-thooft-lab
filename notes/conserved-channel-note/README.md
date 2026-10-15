# Draft note: the singlet bilinear channels at N = ∞

`main.tex` is a standalone manuscript, about seven pages, on the channels
Aₙ = Σᵢ ψᵢ (ad_H)ⁿ ψᵢ of double-scaled SYK at strict N = ∞. It proves two
theorems.

1. **The pair (i ≠ j) amplitude is conserved.** Its chord image is a
   polynomial in the chord Hamiltonian, by a three-term recursion with
   scalar coefficients (`duality/chord_charges.py`). So every moment
   μ_k, k ≥ 1, vanishes, at every λ and every matter weight.
2. **The full A₂ channel is a conserved H² piece plus a dimension-2 matter
   chord.** A₂ − ⟨A₂⟩ = [2λ(1−q)/(1+q)](H² − 1) + ν O₂, so the
   propagating part has width 2(1 − q²) at every λ (1.963369 at
   λ_chord = 2). The proof is a sign-reversing involution on chord
   diagrams (`duality/singlet_tower.py`).

The two fit together. The flavor-diagonal and pair halves are each
O(√N) on a double-scaling line and cancel at that order, so the
physical channel is the O(1) remainder that Theorem 2 describes. The
note also records, as verified but not proved, the n = 3..6 levels: a
polynomial in H plus matter chords of dimension 2..n.

`notes/singlet-tower-note.md` is the longer companion, with the block
kernels, the norm table and the small-λ shares.

**Status: a draft for the author's review. Not for distribution.**
Resolve the two `TODO(author)` markers, affiliation and
acknowledgments, before circulating.

Build with `tectonic main.tex`, or `pdflatex main.tex` run twice, or on
Overleaf. The displayed equations are machine-checked in
`duality/test_chord_charges.py`, `duality/test_full_channel_chord.py`
and `duality/test_singlet_tower.py`.
