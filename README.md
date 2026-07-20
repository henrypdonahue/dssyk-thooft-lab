# Stress-testing the DSSYK∞ ↔ 't Hooft duality

Computational legwork for the duality of Miyashita–Sekino–Susskind,
*"Holograms and Standard Models"* ([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)),
which claims **double-scaled SYK at infinite temperature** (DSSYK∞) is dual to the
**'t Hooft model** (large-$N$ 2D QCD). No new theory — just the target numbers and
numerical stress-tests, each cross-checked against an independent method or analytic anchor.

Paper source and PDF: [`arXiv-2607.05678v1/`](arXiv-2607.05678v1), [`2607.05678v1.pdf`](2607.05678v1.pdf).
Full project assessment and the tiered task ladder: [`road_map.md`](road_map.md).

## Modules

**[`thooft-target/`](thooft-target) — the yardstick ✅**
Solves 't Hooft's singular integral equation for the meson spectrum at the
renormalized-massless point ($\alpha=0$) the duality needs. Chebyshev-U spectral solver;
reproduces Fateev–Lukyanov–Zamolodchikov to ~12–13 digits (the floor is FLZ's published
tables, not the solver — the ground state agrees to ~15), sum rules to $10^{-10}$;
cross-checked by an independent sine-basis solver (low-precision by design, 2–4 digits —
it validates the quark-mass machinery, not the digit count). Every comparison is measured
against this.

**[`syk-self-averaging/`](syk-self-averaging) — the falsifiable claim ✅**
ED test of §3.5's variance suppression, validated against $N\bmod 8$ RMT Bott
periodicity (all three GOE/GUE/GSE classes, across four $N\bmod 8$ values). Singlet
self-averaging rate identified exactly: relative variance $2/\binom{N}{p}$, ED-confirmed
at $N\le18$ (to within the $\sim$10% sampling noise of 200 realizations); followed as
pure combinatorics along the double-scaling line it is super-exponential,
$\sim10^{-46}$ by $N=640$. The sharp $e^{-a\sqrt N}$ symmetry-violation form lives
beyond instance-level ED's reach — an honest wall (see `exact_wick.py` for the way
around it).

**[`duality/`](duality) — the comparison dictionary ✅**
Every convention a DSSYK ↔ 't Hooft comparison needs ($n$-offset, CP, the
$\lambda$-vs-$\lambda_{\rm chord}$ factor 2, $M_n$ normalization), encoded once and
unit-tested *before* any DSSYK number exists, plus the two parameter-free falsifiable
fingerprints: mass-squared ratios and the alternating-CP interleave splittings.

Each folder is self-contained (own README, deps, solver, cross-check, `pytest`).
Run `pip install -r requirements.txt` then `pytest -q` inside any of them
(`-m 'not slow'` skips the long ED sweeps).

## Status & corrections

| Module | Verdict | Tests | Anchor |
|---|---|---|---|
| thooft-target | No known bugs; independently re-derived; eigenfunctions validated | 17/17 | FLZ to 2e-12; sum rules ~9e-11; Parseval to 3e-6 |
| syk-self-averaging | Honest; Eq. 3.28 now verified *exactly* in double scaling | 24 | RMT Bott periodicity; Wick-vs-enumeration exact |
| duality | Conventions + fingerprints encoded before any DSSYK number exists | 8/8 | FLZ reference table |

Load-bearing corrections for any DSSYK ↔ 't Hooft comparison:

1. **Off-by-two.** The paper's dynamical tower starts at $n=2$ ($n=0,1$ are the
   non-dynamical photon/graviton = charge/Hamiltonian). FLZ's ground state is the paper's
   $n=2$: apply $n_{\text{paper}} = n_{\text{FLZ}} + 2$.
2. **Where the match lives.** At **small $\lambda$, fixed $p$, $N\to\infty$**
   ('t Hooft's $\bar\alpha=\bar g^2N=p^2$ fixed) — *not* the $\lambda$-fixed double-scaling
   line where self-averaging sits. At finite $\lambda$ the mesons interact (2506.18054's
   own conclusions: the fixed-$\bar g$ limit is "a theory of interacting meson-strings"),
   so an *exact* match is only expected as $\lambda\to0$ ($q\to1$). Note this regime also
   weakens the symmetry-emergence bound: at fixed $p$ the singlet variance is the power
   law $N^{-p}$ and the symmetry-violation amplitude is
   $\mathrm{rms}(B_{jk})\sim N^{-(p-1)/2}$ (exact: `exact_wick.py`), not $e^{-a\sqrt N}$.
3. **Normalization.** $M_0=Q$ exactly, but $M_1=-(p/2)\,iH$ (not the paper's $M_1=H$;
   forced by symmetry — $\dot Q=0$ makes $M_1$ anti-Hermitian) — compare peaks in the
   same units. Encoded as a test in `syk-self-averaging/test_syk.py`.
4. **Self-averaging.** Exact singlet relative variance is $2/\binom{N}{p}$; this is the
   paper's $p!/N^p$ (Eq. 3.27 carries no factor 2) times the exact factor
   $2\prod_{k=1}^{p-1}(1-k/N)^{-1}\approx 2\,e^{+p^2/2N}$ — an $O(1)$ correction on the
   double-scaling line.

## The reframe that unblocks the headline test

The headline test — does DSSYK reproduce the 't Hooft spectrum? — is *not* "invent a new
formalism." The DSSYK singlet/meson spectrum **is** the SYK fermion-bilinear four-point
spectrum, governed by the standard melonic ladder eigenvalue equation $k(h)=1$; the whole
$M_n=\psi^\dagger(d^n/dt^n)\psi$ tower is frequency-weighting of *one* object, the bilinear
4-point function. So it's "connect known SYK four-point machinery and do one
infinite-temperature solve."

The masses remain genuinely uncomputed — arXiv:2506.18054 §5.3: "At the moment we don't
know the masses of the states created by the operators" (their $S_n$ tower = our $M_n$),
promised for "a future publication". The starting object, however, **is** in print: the
exact all-energy DSSYK 4-point function at infinite temperature exists in chord
transfer-matrix form (Berkooz–Isachenkov–Narovlansky–Torrents, arXiv:1811.02584), and
large-$q$ 4-point functions are known in closed form at all temperatures (Streicher
1911.10171; Choi–Mezei–Sárosi 1912.00004). What is *not* in print is the extraction of
the bilinear-tower masses from it, and whether that spectrum is a **discrete tower or a
continuum is unconfirmed** — the real open risk. (⚠ convention: the chord literature
uses $\lambda_{\rm chord}=2p^2/N$, $q=e^{-\lambda_{\rm chord}}$; this paper's
$\lambda=p^2/N$. See `duality/`.)

## Directions (ranked)

**Do first (near-free).** Self-averaging ⇒ single-realization spectroscopy: the meson
spectral function is a $U(N)$ singlet, so the *measured* $1/\binom{N}{p}$ singlet rate
applies to it (the sharper $e^{-a\sqrt N}$ form is a double-scaling statement beyond
instance-level ED) — single-shot ED is legitimate, no ensemble averaging needed.
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

**Infra minimum — done.** Deps pinned (numpy 2.0.2, scipy 1.13.1, mpmath 1.3.0 —
`eigh(D,S)` is version-sensitive at the 12th digit), `pyproject.toml` at the root, git
history since 2026-07-17. Scripts run from inside each module directory; cross-module
access is file-based (`duality/` reads the certified reference table). Remaining
(CI, `make reproduce`, dashboards) is release-hardening, not a gate.

## Blocked / out of scope

The DSSYK∞ = JT-dS identification, the $1{+}1\to3{+}1$ uplift, and the $\bar g^2(N)$
fine-tuning story are asserted in the paper with no functional form to compute against.
