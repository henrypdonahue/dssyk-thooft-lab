# The 't Hooft meson spectrum at the renormalized-massless point

The **yardstick side** of the DSSYK$_\infty$ ↔ 't Hooft duality test
([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)): meson masses $M_n^2$ and CP
assignments to high precision, validated against Fateev–Lukyanov–Zamolodchikov
(FLZ, [arXiv:0905.2280](https://arxiv.org/abs/0905.2280)) and, at $\alpha\neq0$,
Litvinov–Meshcheriakov ([arXiv:2409.11324](https://arxiv.org/abs/2409.11324)).

## Problem and mass point

't Hooft's equation (FLZ Eq. 1.1) for the light-cone wavefunction $\phi(x)$, $x\in[0,1]$:

$$2\pi^2\lambda\,\phi(x)=\left(\tfrac{\alpha_1}{x}+\tfrac{\alpha_2}{1-x}\right)\phi(x)
-⨍_0^1 dy\,\tfrac{\phi(y)}{(y-x)^2},\qquad \alpha_i=\tfrac{\pi m_i^2}{g^2}-1,$$

with a Hadamard finite-part integral, $\phi(0)=\phi(1)=0$, $M_n^2=2\pi g^2\lambda_n$.
We report $2\lambda_n = M_n^2/(\pi g^2)$.

The duality needs **renormalized $m_q^2=0$**: with $m_q^2=m_{\rm bare}^2-g^2/\pi$ this is
$\alpha\equiv\pi m_q^2/g^2=0$ — FLZ's exactly-studied case (bare mass $g/\sqrt\pi$), a
regular gapped theory, **not** the chiral limit ($\alpha=-1$). Endpoint behavior
$\phi\sim\sqrt{x(1-x)}$ (from $\pi\beta\cot\pi\beta=-\alpha$).

## Method

Rayleigh–Ritz in the basis with the exact endpoint behavior: $x=(1+\xi)/2$,
$\phi=\sqrt{1-\xi^2}\sum a_mU_m(\xi)$ (Chebyshev-U), which diagonalizes the
hypersingular operator exactly, $⨍\sqrt{1-\eta^2}\,U_m/(\xi-\eta)^2 = -\pi(m+1)U_m$.
Result: a generalized symmetric-definite eigenproblem $D\,a=\Lambda S\,a$
($\Lambda=2\pi^2\lambda$) with **exact rational** matrix elements
($D=\pi^2(m{+}1)\delta$; $S$ the SPD Gram matrix; quark mass enters as $4\alpha I_{nm}$,
closed form in `thooft_spectrum.py`). Variational ⇒ monotone convergence from above.
$D,S$ block-diagonalize under $x\to1-x$: the two CP sectors, interleaving into one
alternating-CP trajectory, CP $=(-1)^{n+1}$.

For $\alpha\neq0$ a second production solver (`jacobi_solver.py`) uses the
*matched-exponent* Jacobi basis $[x(1-x)]^\beta P_k$, with $\beta$ solving
$\pi\beta\cot\pi\beta=-\alpha$ exactly — spectral accuracy off the duality
point, where the Chebyshev basis is only algebraic.

## Validation (asserted in `test_thooft.py` + `test_jacobi.py`, 36 tests)

| check | result |
|---|---|
| FLZ Tables 1–2 ($n\le29$) | worst dev 2e-12 (FLZ's rounding floor; ground state ~15 digits vs their Padé) |
| Exact sum rules $s=2,3,4$, both parities | $7\zeta(3)$, $2$ to ~1e-10; $s=3,4$ ζ-closed-forms to ~3e-12 — six whole-spectrum certificates |
| **α≠0, external** (`external_anchors.json`) | LM tables matched at their 6-digit floor (≲2e-5, four α values); LM *exact* α-dependent sum rules to 1e-5–2e-4 via $G^{(2)}(\alpha)-G^{(2)}(0)$ |
| Eigenfunctions (`eigenfunctions.py`) | Parseval $\sum f_n^2=1$ to 3e-6, $\sum(\int x\phi)^2=\tfrac13$ to 2e-6; weak-form Rayleigh quotient recovers eigenvalues from eigenvectors to ~1e-6; cross-solver overlap >0.9999 |
| Independent solver (`sine_solver.py`) | different basis + Gagliardo weak form, no shared code; validates the quark-mass machinery (low-precision by design) |
| **Jacobi solver at α≠0** (`jacobi_solver.py`) | β=½ anchor vs exact closed form ~5e-13; α=0 spectrum vs Chebyshev ~1e-12; α=±0.5 ground state <1e-11 by K=32 (old solver: algebraic ~$K^{-4\beta}$, still 1e-7..1e-5 at K=128); LM sum rules s≥4 to ~1e-12 relative |
| Asymptotics, matrix elements | $2\lambda_n\to n+\tfrac34$; closed forms vs quadrature to machine precision |

**Decay constants** (new data — no published wavefunction tables exist; checked FLZ/AK/LM):
in the physical $\phi_n(0^+)>0$ convention, $f_n=\int\phi_n$ are all positive and
monotonically decreasing, $f_0=+0.9428$.

## Precision notes

- Matrix elements are exact rationals: the accuracy floor is basis truncation, not
  arithmetic. `reference_spectrum.json` (mpmath, $K=160$ vs $220$) carries 12–15 stable
  digits — a self-convergence estimate (both variational, same bias), mildly optimistic;
  external certification is the FLZ tables + sum rules. The CSV prints 14 digits; trust ~11.
- The Chebyshev solver's 12-digit guarantee is an $\alpha=0$ statement: off the
  duality point it hardwires the wrong endpoint exponent and converges only
  algebraically (measured 1e-5–2e-4 in the LM sum-rule residuals). The
  matched-exponent Jacobi solver removes that limitation (quadrature floor
  ~1e-12; honest caveat on subleading endpoint structure in its docstring).

## Files & usage

| file | purpose |
|------|---------|
| `thooft_spectrum.py` | main double-precision solver (takes `alpha`; `--csv` regenerates the CSV) |
| `jacobi_solver.py` | matched-exponent Jacobi solver for α≠0 (`__main__` regenerates `jacobi_convergence.json`) |
| `thooft_highprec.py` / `generate_reference.py` | mpmath solver / certified reference builder |
| `sine_solver.py` | independent cross-check solver |
| `eigenfunctions.py` | wavefunction-level validation + decay constants |
| `external_anchors.json` | FLZ sum rules s=1..8 + LM α≠0 tables and exact sum rules (provenance inside) |
| `test_thooft.py`, `test_jacobi.py` / `validate.py` | asserting suites / human-readable gate (exit code) |
| `reference_spectrum.json`, `spectrum_double_precision.csv`, `jacobi_convergence.json`, `plot_spectrum.py` | outputs |

```bash
pip install -r requirements.txt
python3 thooft_spectrum.py     # print the spectrum
python3 validate.py            # validation report (exit code = pass/fail)
python3 eigenfunctions.py      # wavefunction-level checks + decay constants
pytest -q                      # 36 asserting tests (~1 min; Jacobi sweeps dominate)
```

## Result

$2\lambda_n = M_n^2/(\pi g^2)$ at the duality point (full table in
`reference_spectrum.json`):

| $n$ | $2\lambda_n$ | parity | CP |
|----:|--------------:|:------:|:--:|
| 0 | 0.737061746292690 | sym | $-$ |
| 1 | 1.753731336917500 | anti | $+$ |
| 2 | 2.748160912370600 | sym | $-$ |
| 3 | 3.751057581705400 | anti | $+$ |
| 4 | 4.749295381037500 | sym | $-$ |
| 5 | 5.750492623648700 | anti | $+$ |
