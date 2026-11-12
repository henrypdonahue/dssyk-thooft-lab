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

## What is checked (61 tests)

| check | result |
|---|---|
| FLZ Tables 1–2 (n ≤ 29) | worst deviation 1.1e-12 at K = 300 (`validate.py` prints the live figure). That is one unit of FLZ's own last printed digit: against the 17-digit polished table below, FLZ's printed values are off by up to 1.0e-12 (n = 19), and this solver by 9e-14 at K = 800 (4e-13 at K = 300). Ground state against their Padé: 3e-16 |
| exact sum rules, both parities | all sixteen FLZ Table 3 values, s = 1 (regularized) through s = 8, to within one unit of the table's last printed digit with the transcribed Eq. 1.6 tail (≤ 6e-13 for every G₋ and for G₊ at s ≤ 4; 5e-12 to 7e-11 for the s = 5–8 G₊ entries, which are printed to 1e-11–1e-10); 7ζ(3) and 2 to 2e-12 (1e-10 with the leading tail only); all six s = 1, 3, 4 closed forms against the table to 1e-12 |
| α ≠ 0, external (`external_anchors.json`), Chebyshev solver | LM tables matched at their 6-digit floor; LM exact sum rules to 1e-5–2e-4 |
| eigenfunctions (`eigenfunctions.py`) | Parseval to 3e-6 (a completeness check); weak-form Rayleigh quotient on tanh-sinh nodes to ~1e-12 for n ≤ 5 in both sectors (1e-6 to 6e-6 on Gauss–Legendre, the endpoint-limited original); cross-solver overlap > 0.9999 |
| independent solver (`sine_solver.py`) | a different basis, no shared code; it checks the quark-mass machinery |
| Jacobi solver at α ≠ 0 | β = ½ anchor to ~5e-13; α = ±0.5 ground state below 1e-11 by K = 32; LM sum rules s = 2 to ~3e-8, s ≥ 4 to ~4e-12 |
| asymptotics | 2λn → n + 3/4; matrix elements against quadrature to machine precision |
| FLZ d± constants (`flz_dconstants.py`, their Eq. 4.21 products) | converged in K: d₊ = 0.9631784563984(1), d₋ = 0.43358263983325(1) at K = 800 with the head n ≤ 160, inside the ±5e-13 rounding of FLZ's 12 printed digits and one to two digits past them; d₋/d₊ = √2/π to 2e-15; Jacobi cross-solver to its 1.5e-8 quadrature floor |
| polished reference (`reference_spectrum_polished.json`) | 40 levels to 17–21 digits (2λ₀ = 0.737061746292689589581); the measured truncation rate is K⁻⁶ at every level; the older 24-level table agrees to at least its own claimed digits at every row |

**The decay constants are new data.** No published wavefunction tables
exist. In the φ(0⁺) > 0 convention every fn is positive and decreasing,
with f₀ = +0.9428101504585 (thirteen stable digits), and the tower obeys
fn (n + 3/4) → 2/π from above, with a +0.064/(n + 3/4) correction
(`eigenfunctions.decay_constant_asymptotics`).

FLZ's d± are spectral-determinant constants, not wavefunction data. They
test the eigenvalues only. No external wavefunction anchor exists, so the
wavefunctions are validated internally only.

## Precision notes

Matrix elements are exact rationals, so the accuracy floor is basis
truncation, and truncation falls like K⁻⁶ (measured, every level).

There are two reference tables.

`reference_spectrum_polished.json` is the sharp one: 40 levels, 17–21
digits. Each level is the exact-rational Rayleigh quotient of a
double-precision Ritz vector, evaluated in mpmath. The quotient is
stationary at an eigenvector, so a vector good to 1e-13 gives the
eigenvalue of the truncated problem to 1e-26, at O(K²) per level. A
three-rung basis ladder (K = 600, 900, 1350) then extrapolates the
truncation away; the size of that extrapolation is the stored error bar.
About seven minutes (`generate_reference.py --polish`).

`reference_spectrum.json` is the original: a full mpmath eigensolve at
K = 160 and 220, 24 levels, 12–15 stable digits by self-convergence. It
shares no eigenvector with the polished table, and the two agree to at
least the older table's claimed digits at every row. It is the file
`duality/` reads, and it is unchanged.

The double-precision solver solves the reciprocal problem at α = 0, so
its low levels carry no roundoff growth with K. The CSV is that solve at
K = 800, written to 13 decimals (measured within 1e-13 of the polished
table over all 40 rows). A 14th decimal would be BLAS-dependent noise.

The Chebyshev solver's guarantee holds at α = 0 only. Off that point it
has the wrong endpoint exponent. Use the Jacobi solver there.

## Files and usage

| file | purpose |
|------|---------|
| `thooft_spectrum.py` | the main solver (`--csv` regenerates the CSV) |
| `jacobi_solver.py` | the matched-exponent solver for α ≠ 0 |
| `thooft_highprec.py`, `generate_reference.py` | mpmath solver, the Rayleigh-quotient polish, both reference builders (`--polish` for the 40-level table) |
| `sine_solver.py` | the independent cross-check solver |
| `eigenfunctions.py` | wavefunction checks and decay constants |
| `condensate.py` | the exact condensate (LM Eq. 6.1), the GMOR chiral corner (this repo's solvers at α → −1: 2λ₀ → (2/π)√(a/3), next order ≈ √a measured), and the vacuum-energy curve ε(0) = −0.0407 Nc g² |
| `flz_reference.py` | FLZ Tables 1–2 and the Padé anchor, one provenance-checked copy |
| `external_anchors.json` | FLZ Table 3 sum rules and their closed forms, the d± constants, LM tables and LM sum rules, with provenance inside; all consumed except the LM s = 1 sums (stored, not consumed; the JSON says why) |
| `test_thooft.py`, `test_jacobi.py`, `test_condensate.py`, `validate.py` | the asserting suites, plus the human-readable gate |
| `jacobi_convergence.json`, `condensate.json` | measured convergence and residual records; they set the test tolerances |
| `reference_spectrum_polished.json`, `reference_spectrum.json`, `spectrum_double_precision.csv`, `plot_spectrum.py`, `regge_trajectory.png` | outputs |

Two label conventions coexist, and both are pinned by the tests. The
reference JSON tables label the sectors `sym`/`anti`. `spectrum()` and
the CSV label the same sectors `even`/`odd` (parity under x → 1−x).
`reference_spectrum.json` is a bare list of rows with no provenance block,
because `duality/` reads that format. The other JSON outputs carry
`provenance` and `runtime_s`.

```bash
pip install -e ".[dev]"   # from the repo root, once
python3 thooft_spectrum.py   # print the spectrum
python3 validate.py          # the validation report (exit code = pass/fail)
python3 eigenfunctions.py    # wavefunction checks and decay constants
python3 condensate.py        # condensate, GMOR corner, vacuum energy
python3 flz_dconstants.py    # the FLZ d± anchor, consumed
pytest -q                    # 61 tests (~10 s)
```

## The result

2λn = Mn²/(π g²) at the duality point. The full table is in
`reference_spectrum_polished.json`.

| n | 2λn | parity | CP |
|--:|--------------:|:------:|:--:|
| 0 | 0.737061746292690 | sym | − |
| 1 | 1.753731336917498 | anti | + |
| 2 | 2.748160912370600 | sym | − |
| 3 | 3.751057581705467 | anti | + |
| 4 | 4.749295381037510 | sym | − |
| 5 | 5.750492623648764 | anti | + |
