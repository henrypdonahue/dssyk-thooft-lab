# Stress-testing the DSSYK∞ ↔ 't Hooft duality

Computational legwork for the duality of Miyashita–Sekino–Susskind,
*"Holograms and Standard Models"* ([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)),
which claims **double-scaled SYK at infinite temperature** (DSSYK∞) is dual to the
**'t Hooft model** (large-$N$ 2D QCD). No new theory — just the target numbers and
numerical stress-tests, each cross-checked against an independent method or analytic anchor.

Paper source and PDF: [`arXiv-2607.05678v1/`](arXiv-2607.05678v1), [`2607.05678v1.pdf`](2607.05678v1.pdf).

## Modules

**[`thooft-target/`](thooft-target) — the yardstick ✅**
Solves 't Hooft's singular integral equation for the meson spectrum at the
renormalized-massless point ($\alpha=0$) the duality needs. Chebyshev-U spectral solver;
reproduces Fateev–Lukyanov–Zamolodchikov to ~12–13 digits, sum rules to $10^{-10}$;
independently confirmed by a sine-basis solver. Every comparison is measured against this.

**[`syk-self-averaging/`](syk-self-averaging) — the falsifiable claim ✅**
Brute-force ED test of §3.5's variance suppression. Validated against $N\bmod 8$ RMT Bott
periodicity (all four GOE/GUE/GSE classes). Singlet self-averaging confirmed at its root:
relative variance is exactly $2/\binom{N}{p}$ (the paper's $p!/N^p$), super-exponential
along the double-scaling line down to $\sim10^{-46}$ at $N=640$. The sharp $e^{-a\sqrt N}$
symmetry-violation form lives beyond ED's reach — an honest wall.

Each folder is self-contained (own README, deps, solver, cross-check, `pytest`).
Run `pip install -r requirements.txt` then `pytest -q` inside either
(`-m 'not slow'` skips the long ED sweeps).

## Status & corrections

| Module | Verdict | Tests | Anchor |
|---|---|---|---|
| thooft-target | No known bugs; independently re-derived | 13/13 | FLZ to 2e-12; sum rules ~9e-11 |
| syk-self-averaging | Honest; claims survive scrutiny | 8/8 | RMT Bott periodicity, all 4 classes |

Load-bearing corrections for any DSSYK ↔ 't Hooft comparison:

1. **Off-by-two.** The paper's dynamical tower starts at $n=2$ ($n=0,1$ are the
   non-dynamical photon/graviton = charge/Hamiltonian). FLZ's ground state is the paper's
   $n=2$: apply $n_{\text{paper}} = n_{\text{FLZ}} + 2$.
2. **Where the match lives.** At **small $\lambda$, fixed $p$, $N\to\infty$**
   ('t Hooft's $\bar\alpha=g^2N=p^2$ fixed) — *not* the $\lambda$-fixed double-scaling line
   where self-averaging sits. At finite $\lambda$ the mesons interact, so an *exact* match
   is only expected as $\lambda\to0$ ($q\to1$).
3. **Normalization.** $M_0=Q$ exactly, but $M_1=-(p/2)\,iH$ (not just $H$) — compare peaks
   in the same units.
4. **Self-averaging.** Exact singlet relative variance is $2/\binom{N}{p}$; equals the
   paper's $2\,p!/N^p$ times $\approx e^{+p^2/2N}$ (earlier prose had this sign wrong; no
   numeric output was affected).

## The reframe that unblocks the headline test

The headline test — does DSSYK reproduce the 't Hooft spectrum? — is *not* "invent a new
formalism." The DSSYK singlet/meson spectrum **is** the SYK fermion-bilinear four-point
spectrum, governed by the standard melonic ladder eigenvalue equation $k(h)=1$; the whole
$M_n=\psi^\dagger(d^n/dt^n)\psi$ tower is frequency-weighting of *one* object, the bilinear
4-point function. So it's "connect known SYK four-point machinery and do one
infinite-temperature solve."

The masses remain genuinely uncomputed (arXiv:2506.18054 §5.3 confirms they "are not
currently known"), and one dragon remains: the SYK ladder is always solved in the low-$T$
conformal regime, but the flat-space limit is the opposite corner ($T_B\to\infty$,
$\beta J\to0$) where that ansatz fails. The infinite-$T$ kernel is not in print, and whether
its spectrum is a **discrete tower or a continuum is unconfirmed** — the real open risk.

## Directions (ranked)

**Do first (near-free).** Self-averaging ⇒ single-realization spectroscopy: the meson
spectral function is a $U(N)$ singlet, so the validated $\mathrm{Var}\langle W\rangle<
e^{-a\sqrt N}$ applies to it — single-shot ED is legitimate, no ensemble averaging needed.
Alongside it, an $M_n$ operator-identity + CP audit ($M_0=Q$, $M_1=-(p/2)iH$,
$CP=(-1)^{n+1}$) on the existing harness — laptop-minutes, pins correction #3.

**The real next result.** $N=\infty$ singlet spectral **moments** vs 't Hooft **sum rules** —
the first genuine number-for-number $N=\infty$ comparison. ED-free (chord/Wick
combinatorics), robust to the discrete-vs-continuum unknown, de-risks the full solve.
`thooft-target` already computes the sum rules ($7\zeta(3)$, $2$). Encode the corrections
above once, in a small comparison-dictionary layer, built stubbed *before* the DSSYK number
exists.

**Anchor, then destination.** Large-$q$ closed-form 4-point + Sommerfeld–Watson/Regge
continuation as the analytic $\lambda\to0$ anchor (the one limit where the match should be
exact). Then the full infinite-$T$ bilinear-ladder solve ($\{h:k(h)=1\}$) — gated behind the
moments/anchor, whose first job is discrete-vs-continuum. A "no discrete tower" result would
itself be publishable, about *where* the duality holds.

**Infra minimum.** One `pyproject.toml` so the modules can import each other (today
`self_averaging.py`'s `from syk import …` only works from inside its folder) and pinned deps
(numpy 2.0.2, scipy 1.13.1, mpmath 1.3.0) — `eigh(D,S)` is version-sensitive at the 12th
digit. Everything else (CI, `make reproduce`, dashboards) is release-hardening, not a gate.

## Blocked / out of scope

The DSSYK∞ = JT-dS identification, the $1{+}1\to3{+}1$ uplift, and the $\bar g^2(N)$
fine-tuning story are asserted in the paper with no functional form to compute against.
