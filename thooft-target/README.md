# The 't Hooft meson spectrum at the renormalized-massless point

A clean, reproducible numerical solution of **'t Hooft's equation** for the
meson mass spectrum of large-$N$ QCD in $1{+}1$ dimensions, at the special
quark-mass point required by the DSSYK$_\infty$ ↔ 't Hooft duality of
Miyashita–Sekino–Susskind ([arXiv:2607.05678](https://arxiv.org/abs/2607.05678),
*"Holograms and Standard Models"*).

This is the **target ("yardstick") side** of the duality test: the meson masses
$M_n^2$ and their CP assignments, to high precision, validated against the exact
analytic results of Fateev–Lukyanov–Zamolodchikov (FLZ,
[arXiv:0905.2280](https://arxiv.org/abs/0905.2280)).

## The problem

't Hooft's equation (FLZ Eq. 1.1) for the meson light-cone wavefunction
$\phi(x)$, $x\in[0,1]$:

$$
2\pi^2\lambda\,\phi(x)=\left(\frac{\alpha_1}{x}+\frac{\alpha_2}{1-x}\right)\phi(x)
-\;⨍_0^1 dy\,\frac{\phi(y)}{(y-x)^2},
\qquad \alpha_i=\frac{\pi m_i^2}{g^2}-1,
$$

with a **Hadamard finite-part** (double-pole) integral, boundary conditions
$\phi(0)=\phi(1)=0$, and physical meson masses

$$
M_n^2 = 2\pi g^2\,\lambda_n .
$$

We report the dimensionless eigenvalue $2\lambda_n = M_n^2/(\pi g^2)$.

### Which mass point?

The duality needs the **renormalized quark mass $m_q^2 = 0$**. In the Susskind
et al. conventions the renormalized and bare masses differ by the self-energy,
$m_q^2 = m_{\text{bare}}^2 - g^2/\pi$, so

$$
\alpha \equiv \frac{\pi m_q^2}{g^2} = 0 .
$$

This is **exactly** FLZ's exactly-solved case $\alpha_1=\alpha_2=0$ (bare mass
$m=g/\sqrt\pi$). It is a regular theory with a mass gap and is **not** the chiral
limit (that is $\alpha=-1$, where the massless "pion" $\phi(x)=1$, $M^2=0$
appears). At $\alpha=0$ the endpoint behavior is $\phi(x)\sim\sqrt{x(1-x)}$
(exponent $s=\tfrac12$ from $\pi s\cot\pi s=-\alpha$).

## Method

A **Rayleigh–Ritz / Galerkin spectral method** in the basis that builds in the
correct endpoint behavior (the standard Chebyshev / Multhopp approach used by
't Hooft, Hanson–Peccei–Prasad, and FLZ). Writing $x=(1+\xi)/2$ and

$$
\phi(x)=\sqrt{1-\xi^2}\sum_m a_m\,U_m(\xi),\qquad U_m=\text{Chebyshev 2nd kind},
$$

uses the exact hypersingular identity

$$
⨍_{-1}^{1} d\eta\,\frac{\sqrt{1-\eta^2}\,U_m(\eta)}{(\xi-\eta)^2}=-\pi(m+1)\,U_m(\xi).
$$

The eigenproblem becomes a **generalized symmetric-definite** one,
$D\,a=\Lambda\,S\,a$ with $\Lambda=2\pi^2\lambda$,

$$
D_{nm}=\pi^2(m+1)\,\delta_{nm},\qquad
S_{nm}=\tfrac12\!\left[f(n-m)-f(n+m+2)\right],\quad
f(k)=\begin{cases}\dfrac{2}{1-k^2}&k\text{ even}\\[4pt]0&k\text{ odd}\end{cases}
$$

Here $D$ is diagonal and $S$ is the (SPD) Gram matrix of the basis in
$L^2[0,1]$. Both matrix elements are **exact rationals** (times $\pi^2$), so the
only error sources are basis truncation and arithmetic precision — both
controllable. Because it is variational, eigenvalues converge **monotonically
from above** as the basis grows.

**Parity / CP.** $S$ and $D$ block-diagonalize under $\xi\to-\xi$, i.e.
$x\to 1-x$: even-degree $m$ → $\phi$ symmetric, odd-degree $m$ → antisymmetric.
These are the two CP trajectories. The levels interleave into one Regge
trajectory with alternating CP, $\text{CP}=(-1)^{n+1}$ (Susskind et al.).

### The quark-mass (potential) term, for $\alpha\neq0$

The duality point is $\alpha=0$, where the potential drops out entirely. So the
solver can *only* be validated against FLZ at that one point — unless we also
implement general $\alpha$ and check the potential machinery independently.
Both solvers therefore take an `alpha` argument (default `0.0`) for the
equal-mass case $\alpha_1=\alpha_2=\alpha$. Projecting
$\alpha\!\left(\tfrac1x+\tfrac1{1-x}\right)\phi$ onto the basis adds a matrix
$P=4\alpha\,I$ to $D$, with

$$
I_{nm}=\int_{-1}^{1}U_n(\xi)U_m(\xi)\,d\xi
      =2\big[H_Q-H_R\big],\quad Q=\tfrac{n+m}{2}+1,\ R=\tfrac{|n-m|}{2},\quad
      H_p=\sum_{k=1}^{p}\tfrac1{2k-1}
$$

for $n+m$ even (zero otherwise, so the CP blocks stay decoupled). The factor $4$
is exact: $(1-\xi^2)\!\left(\tfrac1x+\tfrac1{1-x}\right)=4$ under $x=(1+\xi)/2$.
At $\alpha=0$, $P=0$ and everything reduces to the block above.

## Files

| file | purpose |
|------|---------|
| `thooft_spectrum.py`  | double-precision solver (numpy/scipy) — the main entry point; takes `alpha` |
| `thooft_highprec.py`  | arbitrary-precision solver (mpmath); takes `alpha` |
| `sine_solver.py`      | **independent** sine-basis cross-check solver (different basis, different singularity handling) |
| `test_thooft.py`      | asserting `pytest` regression suite (FLZ tables, sum rules, cross-check, …) |
| `validate.py`         | human-readable validation report; exits non-zero if any check regresses |
| `generate_reference.py` | builds the certified reference by comparing two basis sizes (~5 min) |
| `plot_spectrum.py`    | Regge-trajectory plot colored by CP |
| `reference_spectrum.json` | certified high-precision reference table (12–15 digits/level) |
| `spectrum_double_precision.csv` | 40-level spectrum table (double precision; printed to 14 digits, accurate to ~11 — see the precision note) |
| `regge_trajectory.png` | the plot |
| `requirements.txt` | pinned dependencies |

Dependencies: `numpy`, `scipy` (core); `mpmath` (arbitrary precision);
`matplotlib` (plot); `pytest` (tests). `pip install -r requirements.txt`.

## Usage

```bash
python3 thooft_spectrum.py     # print the spectrum (double precision)
python3 validate.py            # human-readable validation report (exit code = pass/fail)
pytest -q                      # asserting regression suite (~5 s)
python3 sine_solver.py         # independent cross-check vs the main solver
python3 thooft_highprec.py     # arbitrary-precision spectrum
```

The solvers take a quark-mass argument for the equal-mass case, e.g.
`spectrum(alpha=0.0)` (the duality point, default) or `alpha=2.0` (heavier
quarks); `solve_sector(K, parity, alpha)` likewise.

## Validation

Against FLZ ([arXiv:0905.2280](https://arxiv.org/abs/0905.2280)), all asserted
in `test_thooft.py`:

- **Eigenvalue tables (Tables 1 & 2):** reproduced for $n\le 29$ to
  $\sim 12$–$13$ digits in double precision; worst deviation $\sim 3\times10^{-12}$.
- **Ground state:** $2\lambda_0 = 0.737061746292690\ldots$, matching FLZ's Padé
  value (Eq. 4.37) to published precision.
- **Exact sum rules (FLZ Eq. 1.8):** independent whole-spectrum check,
  $$\sum_m \lambda_{2m}^{-2}=7\zeta(3),\qquad \sum_m \lambda_{2m+1}^{-2}=2,$$
  reproduced to $\sim10^{-10}$ — plus the $s=3,4$ closed forms
  ($-\tfrac43\pi^2+28\zeta(3)$, …, FLZ Eq. 1.8/Table 3) to $\sim3\times10^{-12}$:
  six independent whole-spectrum certificates in total.
- **Eigenfunctions** (`eigenfunctions.py`): Parseval completeness of the
  decay-constant overlaps ($\sum f_n^2 = 1$ to $3\times10^{-6}$,
  $\sum(\int x\phi_n)^2=\tfrac13$ to $2\times10^{-6}$), independent weak-form
  Rayleigh-quotient recovery of every checked eigenvalue from its eigenvector
  ($\sim10^{-6}$), cross-solver ground-state overlap $>0.9999$. First
  decay-constant table ($f_0=+0.9428$; in the physical $\phi_n(0^+)>0$ sign
  convention all $f_n$ are positive and monotonically decreasing) — no
  published wavefunction tables exist to compare against (checked: FLZ, AK, LM).
- **Asymptotics:** $2\lambda_n\to n+\tfrac34$ ('t Hooft/WKB leading behavior).
- **Independent cross-check** (`sine_solver.py`): a second solver in the
  $\sin(k\pi x)$ basis, sharing *no* code, basis, or closed-form matrix elements
  with the main one — it evaluates the finite-part operator through its exact
  weak (Gagliardo-seminorm) form by Gauss–Legendre quadrature. It reproduces the
  Chebyshev spectrum at $\alpha=0,\,0.5,\,2.0$, converging as $O(1/N)$ from above
  (it lacks the $\sqrt{x(1-x)}$ endpoint, so it is a low-precision check by
  design). This is what independently validates the **quark-mass term**, which
  is invisible at the duality point $\alpha=0$.
- **Matrix elements:** the closed-form potential matrix $I_{nm}$ is checked
  against direct numerical quadrature to machine precision.

### A note on precision

The matrix elements are **exact rationals**, so arithmetic is never the
bottleneck — the accuracy floor is set entirely by the **basis truncation** $K$.
Consequently the mpmath solver buys nothing over double precision at small $K$;
its advantage appears only at large $K$ (the `reference_spectrum.json` table is
built at $K=160$ vs $K=220$ and agrees to $\sim15$ digits). The `stable_digits`
column there measures agreement *between those two basis sizes*, i.e. the
converged-digit count. Because both are variational upper bounds biased the same
way, this is a self-convergence estimate, not a two-sided bound — so it is mildly
optimistic. The external accuracy is instead certified by the FLZ tables and the
sum rules above.

The double-precision `spectrum_double_precision.csv` prints 14 digits per level,
but only about **11** of them are trustworthy — double-precision arithmetic and
basis truncation set the floor. Treat the trailing digits as padding, not
precision; the certified high-precision table is `reference_spectrum.json`.

**Two honesty notes.** (1) The headline value below, $0.737061746292690$, is
FLZ's published Padé estimate; the solver's *own* computed value is
$0.7370617462926896\ldots$ (`reference_spectrum.json`), agreeing with FLZ to
~15 digits. (2) The high-precision ($\sim$12-digit) guarantee is an $\alpha=0$
statement: it rests on the FLZ tables and sum rules, which exist only at the
duality point. Away from $\alpha=0$ the method is *not* equally
accurate: the basis hardwires the $\alpha=0$ endpoint exponent $\tfrac12$, while the
true exponent $\beta(\alpha)$ solves $\pi\beta\cot\pi\beta=-\alpha$, so for
$\alpha\neq0$ convergence degrades from spectral to algebraic. The $\alpha\neq0$
machinery **is** now externally anchored (`external_anchors.json`, asserted in the
tests): the Litvinov–Meshcheriakov tables (arXiv:2409.11324) are reproduced at their
6-digit rounding floor ($\lesssim2\times10^{-5}$ on $2\lambda$, $|\alpha|\le0.5$,
$n<10$), and their *exact* $\alpha$-dependent sum rules are matched to
$10^{-5}$–$2\times10^{-4}$ via the absolutely-convergent difference
$G^{(2)}(\alpha)-G^{(2)}(0)$ — that residual *is* the measured algebraic-convergence
limitation. (A Jacobi-type $[x(1-x)]^{\beta}$ basis would restore spectral accuracy
off the duality point.)

## Result (summary)

The lowest meson mass-squared eigenvalues at the massless point, $2\lambda_n =
M_n^2/(\pi g^2)$, with parity under $x\to1-x$ and CP $=(-1)^{n+1}$:

| $n$ | $2\lambda_n = M_n^2/(\pi g^2)$ | parity | CP |
|----:|-------------------------------:|:------:|:--:|
| 0 | 0.737061746292690 | symmetric | $-$ |
| 1 | 1.753731336917500 | antisymmetric | $+$ |
| 2 | 2.748160912370600 | symmetric | $-$ |
| 3 | 3.751057581705400 | antisymmetric | $+$ |
| 4 | 4.749295381037500 | symmetric | $-$ |
| 5 | 5.750492623648700 | antisymmetric | $+$ |

(Full, higher-precision table in `reference_spectrum.json`.)
