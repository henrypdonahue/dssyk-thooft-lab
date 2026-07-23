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
| [`thooft-target/`](thooft-target) | the yardstick: 't Hooft meson spectrum, $\alpha=0$ (Chebyshev) and $\alpha\neq0$ (matched-exponent Jacobi) | 36 | FLZ to 2e-12; six exact sum rules to ~1e-10; α≠0 spectrally accurate, externally anchored (LM tables + exact sum rules) |
| [`syk-self-averaging/`](syk-self-averaging) | the falsifiable claims: §3.5 self-averaging + the $M_n$/CP sector + the U(N)/QED campaign (rung 19) + the baryon sector (rung 23) | 86 | RMT Bott periodicity; exact-Wick closed forms vs enumeration and ED |
| [`duality/`](duality) | the comparison dictionary + fingerprints (built before any DSSYK number exists), boost-to-mass calibration, large-q anchor, anti-scrambling sign test (rungs 25–26) | 38 | FLZ reference table; exact Rindler spectral function; ED-locked large-q bridge; folded-correlator ED lab |

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
't Hooft decay-constant table (all $f_n>0$ in the $\phi(0^+)>0$ convention);
CP-resolved spectral-moment invariants along fixed $p$ (`moments_pipeline.py`);
U(N)-class self-averaging exact + charging curves with GLS statistics
(`qed_campaign.py`, rung 19); baryon-sector first contact — exact 2-point collapse,
$M_B$ extensivity, edge-curvature hierarchy (`baryons.py`, rung 23).

**New calibration (rung 13, `duality/BOOST.md`):** one discrete mass appears in
boost frequency as a *continuum with mass-encoding nodes* — peak-hunting is the
wrong observable, and the ω-symmetric $T_B=\infty$ hologram vs KMS-asymmetric bulk
response ("tomperature" map) is the precisely-stated gating unknown.

**Audited (rung 12):** the systematics behind §3.5 were never published anywhere —
`exact_wick.py` is the only prefactor-exact derivation in existence and *confirms*
the papers at leading order; the λ factor-2 clash is internal to the Susskind corpus;
the 2511.10907 entropy correction leaves §5 untouched.

**New external constraint (2026-07-23, `duality/ANTISCRAMBLING.md`):**
Cui–Kolchmeyer 2607.13665 (cites the target paper) argue dS observer OTOCs
must *anti-scramble*; Harlow–Zhao 2607.14215 propose bounded-spectrum QM with
Euclidean-folded correlators as the mechanism and pose the SYK test as open.
This repo ran both: DSSYK∞ **scrambles** at every temperature
($\mathrm{Re}\,a=-1$ exactly; ED-confirmed) — the flat dictionary cannot be a
dS observer at OTOC level, sharpening the tomperature gap to "must flip the
sign" — and the first SYK fold data (mechanics work; the Lorentzian response
tracks the continuation; the Euclidean $N$-trend is the open question).
Time-ordered physics — the entire spectrum-match program — is untouched.

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
