# Stress-testing the DSSYK∞ ↔ 't Hooft duality

Computational legwork for Miyashita–Sekino–Susskind, *"Holograms and Standard Models"*
([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)): DSSYK at infinite temperature is
claimed dual to the 't Hooft model (large-$N$ 2D QCD). No new theory — target numbers and
numerical stress-tests, each cross-checked against an independent method or analytic anchor.

- Paper: [`arXiv-2607.05678v1/`](arXiv-2607.05678v1), [`2607.05678v1.pdf`](2607.05678v1.pdf)
- Assessment + task ladder + **progress log**: [`road_map.md`](road_map.md)

## Modules

| module | role | tests | headline anchor |
|---|---|---|---|
| [`thooft-target/`](thooft-target) | the yardstick: 't Hooft meson spectrum at the duality point $\alpha=0$ | 20 | FLZ to 2e-12; six exact sum rules to ~1e-10; α≠0 externally anchored (LM tables + exact sum rules) |
| [`syk-self-averaging/`](syk-self-averaging) | the falsifiable claims: §3.5 self-averaging + the $M_n$/CP/QED sector | 24 | RMT Bott periodicity; exact-Wick closed forms vs enumeration and ED |
| [`duality/`](duality) | the comparison dictionary + falsifiable fingerprints, built before any DSSYK number exists | 8 | FLZ reference table |

Each folder is self-contained (README, pinned deps, solver, cross-check, `pytest`).
`pip install -r requirements.txt`, then `pytest -q` inside any of them
(`-m 'not slow'` skips long ED sweeps).

## Findings so far (short version; details in `road_map.md`)

**Survives scrutiny — often sharper than the paper stated:**
- Singlet self-averaging (Eq. 3.27): rate exactly $2/\binom{N}{p}$, super-exponential in
  double scaling ($\sim10^{-46}$ at $N=640$).
- Symmetry violation (Eq. 3.28): **verified exactly** in the double-scaling regime via
  closed-form Wick combinatorics — $\mathrm{Var}(B_{jk})=4\sigma^4\binom{N-2}{p-1}$, so
  $\ln\mathrm{rms}=-\tfrac14\sqrt N\ln N$, beating $e^{-a\sqrt N}$ for every $a$
  (`exact_wick.py`; previously an ED-unreachable "honest wall").

**Corrected (bookkeeping, not structure):** the four numbered corrections below, plus
$M_2=-\sum\dot\psi^\dagger\dot\psi$ (paper's $+$ sign is wrong), CP must be the *unitary*
particle-hole map, and raw $M_{n\ge2}$ are CP-mixtures — CP lives in resolved correlators.

**New data:** first CP-resolved $M_n$ spectral functions (`mn_spectroscopy.py`); first
't Hooft decay-constant table (all $f_n>0$ in the $\phi(0^+)>0$ convention).

**Still untested:** the headline spectrum match itself, and whether the DSSYK singlet
spectrum is a discrete tower or a continuum — the real open risk.

## Load-bearing corrections for any comparison

1. **Off-by-two.** The paper's dynamical tower starts at $n=2$ ($n=0,1$ = charge/Hamiltonian,
   non-dynamical): $n_{\text{paper}} = n_{\text{FLZ}} + 2$.
2. **Where the match lives.** Fixed $p$, $\lambda\to0$ ($\bar\alpha=\bar g^2N=p^2$ fixed) —
   *not* the double-scaling line. There the mesons interact (per 2506.18054 itself), and
   symmetry emergence is only polynomial: singlet variance $N^{-p}$, violation amplitude
   $N^{-(p-1)/2}$ (exact).
3. **Normalization.** $M_0=Q$; $M_1=-(p/2)\,iH$, not the paper's $M_1=H$ (forced:
   $\dot Q=0$ makes $M_1$ anti-Hermitian). Tested in `syk-self-averaging`.
4. **Conventions.** Exact rate $2/\binom{N}{p}$ = paper's $p!/N^p$ times
   $2\prod(1-k/N)^{-1}$; and the chord literature uses $\lambda_{\rm chord}=2p^2/N$,
   $q=e^{-\lambda_{\rm chord}}$ vs the paper's $\lambda=p^2/N$. Encoded in `duality/`.

## The headline path

The DSSYK meson spectrum **is** the fermion-bilinear 4-point spectrum; the $M_n$ tower is
frequency-weighting of one object. The masses are unpublished (2506.18054 §5.3: "we don't
know the masses…", promised for a future paper — a live race), but the starting object is
in print: the exact $T=\infty$ chord 4-point function (arXiv:1811.02584) and closed-form
large-$q$ 4-point functions. Remaining ladder (rungs in `road_map.md`): boost-to-mass
dictionary (13) → $N=\infty$ moments vs sum rules (15) → large-$q$ Regge anchor (17) →
extract the tower from the chord 4-point function (18). Discrete-vs-continuum is the first
deliverable and is publishable either way.

## Blocked / out of scope

DSSYK∞ = JT-dS itself (assumed by the authors), the $1{+}1\to3{+}1$ uplift, and the
$\bar g^2(N)$ fine-tuning curve — asserted with no functional form to compute against.
