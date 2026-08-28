# The 't Hooft meson spectrum at the massless point

This is the **yardstick**. Meson masses and CP assignments, to high
precision, on the side of the duality that is already solved.

The numbers are checked against Fateev–Lukyanov–Zamolodchikov (FLZ,
[arXiv:0905.2280](https://arxiv.org/abs/0905.2280)). At nonzero quark
mass they are also checked against Litvinov–Meshcheriakov
([arXiv:2409.11324](https://arxiv.org/abs/2409.11324)).

## The problem, and the mass point

't Hooft's equation for the light-cone wavefunction φ(x) on [0, 1] is
FLZ Eq. 1.1:

    2π²λ φ(x) = α (1/x + 1/(1−x)) φ(x) − FP∫ φ(y)/(y−x)² dy

The integral is a Hadamard finite part. The boundary conditions are
φ(0) = φ(1) = 0. Per quark, α = π m_bare²/g² − 1. Masses are
Mn² = 2π g² λn, and we report 2λn = Mn²/(π g²).

The duality needs **renormalized m_q² = 0**. Since m_q² = m_bare² − g²/π,
that means α = 0.

Be careful here. α = 0 is FLZ's exactly studied case, and it is a regular
gapped theory. It is **not** the chiral limit, which is α = −1. The
endpoint behavior is φ ~ √(x(1−x)), from πβ cot(πβ) = −α.

## Method

Rayleigh–Ritz in a basis that already has the exact endpoint behavior:
Chebyshev-U times √(1−ξ²). That basis diagonalizes the hypersingular
operator exactly.

What comes out is a symmetric-definite eigenproblem with **exact
rational** matrix elements. So convergence is variational, and monotone
from above.

The two parity sectors interleave into one alternating-CP trajectory,
CP = (−1)^(n+1).

For α ≠ 0 there is a second solver, `jacobi_solver.py`. It uses the
matched-exponent Jacobi basis [x(1−x)]^β, with β solving
πβ cot(πβ) = −α exactly. That restores spectral accuracy off the duality
point, where the Chebyshev basis is only algebraic.

## What is checked (45 tests)

| check | result |
|---|---|
| FLZ Tables 1–2 (n ≤ 29) | worst deviation ~2e-12 at the tested K = 300, which is the level of FLZ's own last printed digit. It is platform-dependent at the 12th digit, ≲1e-11 across configurations; `validate.py` prints the live figure. Ground state against their Padé: 8.5e-15 at K = 300, and 3.7e-16 from the high-precision reference |
| exact sum rules s = 2, 3, 4, both parities | 7ζ(3) and 2 to ~1e-10; the s = 3, 4 closed forms to ~3e-12 |
| α ≠ 0, external (`external_anchors.json`) | LM tables matched at their 6-digit floor; LM exact sum rules to 1e-5–2e-4 |
| eigenfunctions (`eigenfunctions.py`) | Parseval to 3e-6; weak-form Rayleigh quotient to 1e-6–6e-6 for n ≤ 2; cross-solver overlap > 0.9999 |
| independent solver (`sine_solver.py`) | a different basis, no shared code; it checks the quark-mass machinery |
| Jacobi solver at α ≠ 0 | β = ½ anchor to ~5e-13; α = ±0.5 ground state below 1e-11 by K = 32; LM sum rules s = 2 to ~3e-8, s ≥ 4 to ~4e-12 |
| asymptotics | 2λn → n + 3/4; matrix elements against quadrature to machine precision |

**The decay constants are new data.** No published wavefunction tables
exist. In the φ(0⁺) > 0 convention every fn is positive and decreasing,
with f₀ = +0.9428.

## Precision notes

Matrix elements are exact rationals, so the accuracy floor is basis
truncation. `reference_spectrum.json` carries 12–15 stable digits by
self-convergence. The external certification is the FLZ tables and the
sum rules. Trust about 11 digits of the CSV.

The Chebyshev solver's 12-digit guarantee holds at α = 0 only. Off that
point it has the wrong endpoint exponent. Use the Jacobi solver there.

## Files and usage

| file | purpose |
|------|---------|
| `thooft_spectrum.py` | the main solver (`--csv` regenerates the CSV) |
| `jacobi_solver.py` | the matched-exponent solver for α ≠ 0 |
| `thooft_highprec.py`, `generate_reference.py` | mpmath solver, reference builder |
| `sine_solver.py` | the independent cross-check solver |
| `eigenfunctions.py` | wavefunction checks and decay constants |
| `condensate.py` | the exact condensate (LM Eq. 6.1), the GMOR chiral corner (first data at α → −1: 2λ₀ → (2/π)√(a/3), next order ≈ √a measured), and the vacuum-energy curve ε(0) = −0.0407 Nc g² |
| `flz_reference.py` | FLZ Tables 1–2 and the Padé anchor, one provenance-checked copy |
| `external_anchors.json` | FLZ sum rules and LM tables, with provenance inside |
| `test_thooft.py`, `test_jacobi.py`, `test_condensate.py`, `validate.py` | the asserting suites, plus the human-readable gate |
| `jacobi_convergence.json`, `condensate.json` | measured convergence and residual records; they set the test tolerances |
| `reference_spectrum.json`, `spectrum_double_precision.csv`, `plot_spectrum.py`, `regge_trajectory.png` | outputs |

```bash
pip install -r requirements.txt
python3 thooft_spectrum.py   # print the spectrum
python3 validate.py          # the validation report (exit code = pass/fail)
python3 eigenfunctions.py    # wavefunction checks and decay constants
python3 condensate.py        # condensate, GMOR corner, vacuum energy
pytest -q                    # 45 tests (~3 s)
```

## The result

2λn = Mn²/(π g²) at the duality point. The full table is in
`reference_spectrum.json`.

| n | 2λn | parity | CP |
|--:|--------------:|:------:|:--:|
| 0 | 0.737061746292690 | sym | − |
| 1 | 1.753731336917499 | anti | + |
| 2 | 2.748160912370602 | sym | − |
| 3 | 3.751057581705472 | anti | + |
| 4 | 4.749295381037520 | sym | − |
| 5 | 5.750492623648781 | anti | + |
