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
| [`thooft-target/`](thooft-target) | the yardstick: 't Hooft meson spectrum, $\alpha=0$ (Chebyshev) and $\alpha\neq0$ (matched-exponent Jacobi); condensate/GMOR corner + vacuum energy (rung 21) | 42 | FLZ to 2e-12; six exact sum rules to ~1e-10; α≠0 spectrally accurate, externally anchored (LM tables + exact sum rules); LM (6.1) condensate + chiral-corner GMOR data |
| [`syk-self-averaging/`](syk-self-averaging) | the falsifiable claims: §3.5 self-averaging + the $M_n$/CP sector + the U(N)/QED campaign (rung 19) + the baryon sector (rung 23) + freeness onset (rung 27) + the Majorana singlet-moment bench (rung 15) | 103 | RMT Bott periodicity; exact-Wick closed forms vs enumeration and ED |
| [`duality/`](duality) | the comparison dictionary + fingerprints, boost-to-mass calibration, large-q anchor, anti-scrambling sign test (rungs 25–26), **the exact $N=\infty$ chord engine** (rungs 15/18/26) + KMS-axis discriminator (rung 24) | 141 | FLZ reference table; exact Rindler spectral function; ED-locked large-q bridge; chord engine validated 3 independent ways, thermal 2-pt vs ED to 0.6–1.6% |

Each folder is self-contained (README, solver, cross-check, `pytest`). Install the pinned
scientific stack once from the repo root (`pip install -e ".[dev]"`), then run `pytest -q`
inside any module (`-m 'not slow'` skips long ED sweeps).

## Findings so far (short version; details in `road_map.md`)

**Survives scrutiny — often sharper than the paper stated:**
- Singlet self-averaging (Eq. 3.27): rate exactly $2/\binom{N}{p}$. But this is an
  elementary i.i.d.-$\chi^2$ identity for a sum of squared Gaussians — true regardless of
  any diagram, so it is *not* itself a test of the paper's diagrammatic self-averaging
  mechanism.
- Symmetry violation (Eq. 3.28) — the load-bearing claim, and the real test of the
  mechanism: **verified exactly** in the double-scaling regime via closed-form Wick
  combinatorics — $\mathrm{Var}(B_{jk})=4\sigma^4\binom{N-2}{p-1}$, so
  $\ln\mathrm{rms}=-\tfrac14\sqrt N\ln N$, beating $e^{-a\sqrt N}$ for every $a$
  (`exact_wick.py`; previously an ED-unreachable "honest wall").

**Corrected (bookkeeping, not structure):** the four numbered corrections below, plus
$M_2=-\sum\dot\psi^\dagger\dot\psi$ (paper's $+$ sign is wrong), CP must be the *unitary*
particle-hole map, and raw $M_{n\ge2}$ are CP-mixtures — CP lives in resolved correlators.

**New data:** first CP-resolved $M_n$ spectral functions (`mn_spectroscopy.py`); first
't Hooft decay-constant table (all $f_n>0$ in the $\phi(0^+)>0$ convention); CP-resolved
spectral-moment invariants along fixed $p$ (`moments_pipeline.py`);
U(N)-class self-averaging exact + charging curves with GLS statistics
(`qed_campaign.py`, rung 19); baryon-sector first contact — exact 2-point collapse,
$M_B$ extensivity, edge-curvature hierarchy (`baryons.py`, rung 23); **first
exhibition of CK's early/late free-product transition** — $\mathrm{spec}(E{+}F(t))$
lands on the arcsine law with a rigid-spectrum $\sim1/\mathrm{dim}$ floor, onset at
the dissipation time (`freeness.py`, rung 27; the $\ln N$ scrambling-time
separation is beyond ED).

**New calibration (rung 13, `duality/BOOST.md`):** one discrete mass appears in
boost frequency as a *continuum with mass-encoding nodes* — peak-hunting is the
wrong observable, and the ω-symmetric $T_B=\infty$ hologram vs KMS-asymmetric bulk
response ("tomperature" map) is the precisely-stated gating unknown.

**Audited (rung 12):** the systematics behind §3.5 were never published anywhere —
`exact_wick.py` is the only complete prefactor-exact derivation we're aware of and
*confirms* the papers at leading order; the λ factor-2 clash is internal to the Susskind
corpus; the 2511.10907 entropy correction leaves §5 untouched.

**New external constraint (2026-07-23, `duality/ANTISCRAMBLING.md`):**
Cui–Kolchmeyer 2607.13665 (cites the target paper) argue dS observer OTOCs
must *anti-scramble*; Harlow–Zhao 2607.14215 propose bounded-spectrum QM with
Euclidean-folded correlators as the mechanism and pose the SYK test as open.
This repo ran both: DSSYK∞ **scrambles** at every temperature (the growing OTOC
term has $\mathrm{Re}\,a=-1$) — the magnitude *growth* is ED-confirmed, while the
scrambling *sign* is analytic (Streicher-transcribed), inheriting an
empirically-fixed crossed-operator ordering rather than an independent ED sign
measurement. So the flat dictionary cannot be a dS observer at OTOC level,
sharpening the tomperature gap to "must flip the sign" — plus first SYK fold data
(mechanics work; the Lorentzian response tracks the continuation) and the fixed-λ
scan verdict on the Euclidean $N$-trend: **a genuine $N$-obstruction, not a
λ-artifact** — the folded connected part is not $1/N$-suppressed, so HZ's
"continuation commutes with the semiclassical limit" fails at ED sizes
(unfolded control $N$-stable; the chord-side fold has since settled it — next block).
Time-ordered physics — the entire spectrum-match program — is untouched.

**The $N=\infty$ chord verdicts (2026-08-16, `duality/CHORD.md`):** the exact
double-scaled theory is now implemented as finite linear algebra (chord transfer
matrix + matter contours, validated three independent ways; thermal 2-pt vs ED to
0.6–1.6%). Three results: **(rung 18)** at $T_B=\infty$, every spectral channel is a
*structureless continuum* at every λ, converging to the semiclassical
sech$^{2\Delta}$ transform linearly in λ — no tower emerges; **(rung 15)** sharper:
the H-derived singlet bilinear channels $A_n=\sum_i\psi_i(\mathrm{ad}_H)^n\psi_i$
are *exactly conserved* at $N=\infty$ (μ₂ = μ₄ = 0 to machine precision at every λ;
$A_1=-2pH$ is the exact anchor) — **the meson tower is a 1/N effect invisible at
every order of the double-scaled expansion**, and the ED bench shows it turning on
along fixed $p$, λ→0 — exactly the 't Hooft corner. The match lives strictly outside
the strict double-scaling limit; **(rung 26 decider #2)** the HZ Euclidean fold fails
*in the double-scaled theory itself*: the folded connected part stays $O(1)$ as λ→0
(no 1/N suppression), the $N$-scaled folded ratio runs away from the continued
closed form $\propto1/\lambda$ — closing the fold as an anti-scrambling mechanism
for DSSYK∞. Plus **(rung 21)** the bulk vacuum-energy half: LM's exact condensate
transcribed and verified against the solver in the chiral corner (first
$\alpha\to-1$ data: $2\lambda_0\to(2/\pi)\sqrt{a/3}$, measured NLO $\simeq\sqrt a$);
the entire 2d counterpart of §5's fine-tuning drama is one finite number,
$\varepsilon(0)=-0.0407\,N_cg^2$ — nothing to tune. And **(rung 24,
`duality/DISCRIMINATOR.md`)** the KMS/OTOC axis now discriminates the three
competing duals: sine-dilaton's fake-temperature black-hole reading matches the
measured data as-is; MSS flat-space owes a sign-flipping map (the fold is now
excluded); NV engineers dS-KMS and so faces the CK no-go head-on, pending a
computable OTOC sign.

**Still untested:** the headline spectrum match itself — now precisely located: the
tower must be extracted along fixed $p$, λ→0 (a $1/N$ effect from the double-scaled
viewpoint), through the rung-17 $O(1/q)$ residues or direct fixed-$p$ methods; the
$N=\infty$ chord side is closed.

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
know the masses…", promised for a future paper — a live race). The first deliverable —
discrete-vs-continuum at $N=\infty$ — is now **delivered** (`duality/CHORD.md`): the exact
chord answer is *neither* — the singlet channels are conserved at strict double scaling,
so the tower is a $1/N$ effect living along fixed $p$, λ→0. Remaining ladder (rungs in
`road_map.md`): the $O(1/q)$ residues on the large-$q$ skeleton (17 proper) and/or direct
fixed-$p$ methods at growing $p$, through the boost-to-mass map (13, sign-gated by rung
25). The race is now for the *subleading* computation — the one the authors promised.

## Blocked / out of scope

DSSYK∞ = JT-dS itself (assumed by the authors), the $1{+}1\to3{+}1$ uplift, and the
$\bar g^2(N)$ fine-tuning curve — asserted with no functional form to compute against.
