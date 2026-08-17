# The 't Hooft meson spectrum at the renormalized-massless point

The **yardstick side** of the duality test: meson masses and CP
assignments to high precision, validated against Fateev–Lukyanov–
Zamolodchikov (FLZ, [arXiv:0905.2280](https://arxiv.org/abs/0905.2280))
and, at nonzero quark mass, Litvinov–Meshcheriakov
([arXiv:2409.11324](https://arxiv.org/abs/2409.11324)).

## Problem and mass point

't Hooft's equation (FLZ Eq. 1.1) for the light-cone wavefunction φ(x) on
[0, 1]:

    2π²λ φ(x) = α (1/x + 1/(1−x)) φ(x) − FP∫ φ(y)/(y−x)² dy

with a Hadamard finite-part integral, φ(0) = φ(1) = 0, and
α = π mq²/g² − 1 per quark. Masses are Mn² = 2π g² λn; we report
2λn = Mn²/(π g²).

The duality needs **renormalized mq² = 0**. With mq² = m_bare² − g²/π this
is α = 0 — FLZ's exactly studied case, a regular gapped theory, **not**
the chiral limit (α = −1). Endpoint behavior: φ ~ √(x(1−x)), from
πβ cot(πβ) = −α.

## Method

Rayleigh–Ritz in a basis with the exact endpoint behavior: Chebyshev-U
times √(1−ξ²), which diagonalizes the hypersingular operator exactly. The
result is a symmetric-definite eigenproblem with **exact rational** matrix
elements, so convergence is variational (monotone from above). The two
parity sectors interleave into one alternating-CP trajectory,
CP = (−1)^(n+1).

For α ≠ 0 a second solver (`jacobi_solver.py`) uses the matched-exponent
Jacobi basis [x(1−x)]^β with β solving πβ cot(πβ) = −α exactly. That
restores spectral accuracy off the duality point, where the Chebyshev
basis is only algebraic.

## Validation (42 tests)

| check | result |
|---|---|
| FLZ Tables 1–2 (n ≤ 29) | worst deviation 2e-12 — FLZ's own rounding floor |
| Exact sum rules s = 2, 3, 4, both parities | 7ζ(3) and 2 to ~1e-10; s = 3, 4 closed forms to ~3e-12 |
| α ≠ 0, external (`external_anchors.json`) | LM tables matched at their 6-digit floor; LM exact sum rules to 1e-5–2e-4 |
| Eigenfunctions (`eigenfunctions.py`) | Parseval to 3e-6; weak-form Rayleigh quotient to ~1e-6; cross-solver overlap > 0.9999 |
| Independent solver (`sine_solver.py`) | different basis, no shared code; validates the quark-mass machinery |
| Jacobi solver at α ≠ 0 | β = ½ anchor to ~5e-13; α = ±0.5 ground state < 1e-11 by K = 32; LM sum rules to ~1e-12 |
| Asymptotics | 2λn → n + 3/4; matrix elements vs quadrature to machine precision |

**Decay constants** (new data — no published wavefunction tables exist):
in the φ(0⁺) > 0 convention, all fn are positive and decreasing,
f₀ = +0.9428.

## Precision notes

- Matrix elements are exact rationals, so the accuracy floor is basis
  truncation. `reference_spectrum.json` carries 12–15 stable digits by
  self-convergence; external certification is the FLZ tables and sum
  rules. Trust ~11 digits of the CSV.
- The Chebyshev solver's 12-digit guarantee holds at α = 0 only. Off that
  point it has the wrong endpoint exponent; use the Jacobi solver there.

## Files and usage

| file | purpose |
|------|---------|
| `thooft_spectrum.py` | main solver (`--csv` regenerates the CSV) |
| `jacobi_solver.py` | matched-exponent solver for α ≠ 0 |
| `thooft_highprec.py`, `generate_reference.py` | mpmath solver, reference builder |
| `sine_solver.py` | independent cross-check solver |
| `eigenfunctions.py` | wavefunction checks + decay constants |
| `condensate.py` | exact condensate (LM Eq. 6.1), the GMOR chiral corner (first data at α → −1: 2λ₀ → (2/π)√(a/3), next order ≈ √a measured), and the vacuum-energy curve: ε(0) = −0.0407 Nc g² |
| `external_anchors.json` | FLZ sum rules + LM tables (provenance inside) |
| `test_thooft.py`, `test_jacobi.py`, `test_condensate.py`, `validate.py` | asserting suites; human-readable gate |
| `reference_spectrum.json`, `spectrum_double_precision.csv`, `plot_spectrum.py` | outputs |

```bash
pip install -r requirements.txt
python3 thooft_spectrum.py   # print the spectrum
python3 validate.py          # validation report (exit code = pass/fail)
python3 eigenfunctions.py    # wavefunction checks + decay constants
python3 condensate.py        # condensate + GMOR corner + vacuum energy
pytest -q                    # 42 tests (~1 min)
```

## Result

2λn = Mn²/(π g²) at the duality point (full table in
`reference_spectrum.json`):

| n | 2λn | parity | CP |
|--:|--------------:|:------:|:--:|
| 0 | 0.737061746292690 | sym | − |
| 1 | 1.753731336917500 | anti | + |
| 2 | 2.748160912370600 | sym | − |
| 3 | 3.751057581705400 | anti | + |
| 4 | 4.749295381037500 | sym | − |
| 5 | 5.750492623648700 | anti | + |
