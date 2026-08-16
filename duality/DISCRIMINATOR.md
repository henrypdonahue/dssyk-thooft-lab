# The KMS/OTOC axis: discriminating the three competing duals

*Rung 24 (KMS-axis note, 2026-08-16). Companion to `ANTISCRAMBLING.md`
(rungs 25–26) and `CHORD.md` (rungs 15/18). Sources fetched and read for
this note: Narovlansky–Verlinde 2310.16994 (§2.3, §3), Blommaert–Mertens–
Papalini 2404.03535 (§1), Cui–Kolchmeyer 2607.13665 and Harlow–Zhao
2607.14215 (repo root). Repo-side numbers quoted here are asserted from
committed JSONs in `test_discriminator.py`.*

## The axis

Cui–Kolchmeyer: a dS observer's OTOC must **anti-scramble** (negative
eikonal phase, their $c<0$), and the MSS chaos bound forces $c>0$
(scrambling) for **any unitary QM in a KMS state at the dS temperature** —
so a dS hologram must break one of: unitarity, KMS-at-$\beta_{dS}$, or the
remaining MSS hypotheses. This repo's measured/exact facts on that axis:

- **Rung 25:** under the flat dictionary DSSYK$_\infty$ scrambles at every
  temperature — the growing OTOC term has $\mathrm{Re}\,a(v) = -1$
  exactly, rate $\lambda_L = 2\pi v/\beta$.
- **Rung 26 (both deciders):** the Harlow–Zhao Euclidean fold — the one
  proposed *mechanism* for converting a bounded-spectrum scrambler into an
  anti-scrambler — fails quantitatively in DSSYK: the folded connected
  part is not $1/N$-suppressed at ED sizes (decider #1) **nor at
  $N=\infty$, where it grows toward the semiclassical corner** (decider
  #2, `chord_fold.json`).

## Where each program sits

**1. MSS flat-space ('t Hooft dual, 2607.05678 — this repo's target).**
The hologram is the $T_B=\infty$ maximal-entropy state: *not* KMS at the
dS temperature, so the CK no-go is dodged in letter. But then the
dictionary owes an explanation of how a dS observer's thermal physics
(and its $c<0$) emerges — the "tomperature" map. Rungs 25/26 sharpen the
debt: the map must flip the OTOC sign, and the HZ fold cannot be that map
in this model. **Status: open falsification pressure; the OTOC axis is
where this proposal is weakest.**

**2. Narovlansky–Verlinde doubled DSSYK (2310.16994).** Their physical
operators are dressed and **shifted by $\pm i\beta_{dS}/4$, with
$\beta_{dS}$ "adjusted... to match the KMS property of the correlation
functions in the de Sitter vacuum"** (their §2.3). The construction
*engineers* KMS at the dS temperature — exactly the CK hypothesis. So the
NV program faces the no-go head-on: unitary QM + their KMS structure ⇒
MSS scrambling, while their dS interpretation needs anti-scrambling —
unless the equal-energy constraint / dressing (their crossed-product-like
structure, §2.3) breaks an MSS hypothesis. **The discriminating
computation, not yet published by anyone: the OTOC sign of NV's dressed
$\mathcal O^\pm_\Delta$.** The chord engine here evaluates ordinary
DSSYK contours; the NV objects need the two-sided constrained ($H_L=H_R$)
sector — adjacent tooling, flagged as the sharpest next computation on
this axis.

**3. Sine-dilaton gravity (2404.03535).** Their bulk is a **black hole**
whose "Lorentzian spacetime probed by matter is a smooth black hole with
Hawking temperature equal to the fake temperature
$\beta_{BH} = 2\pi/\sin\theta$", with "all matter correlators... thermal
at the fake temperature" (their §1), and DSSYK as a *sub-maximal
scrambler*. That is the scrambling (AdS-sign) side of CK's dichotomy —
**the only one of the three programs whose stated thermal structure
matches this repo's measured OTOC data with no additional map**: fake-
temperature KMS + scrambling is precisely $\lambda_L = 2\pi v/\beta =
2\pi/\beta_{\rm fake}$ with $\mathrm{Re}\,a = -1$ (rung 25; the identity
$\beta_{\rm fake}\lambda_L = 2\pi$ ties their (1.7) to this repo's $v$
machinery and is asserted in the tests). CK pressure attaches only to the
program's further dS3 reinterpretation (their §1 relation to Verlinde's
proposal), not to the black-hole reading.

## Verdict of the note

The OTOC-sign/KMS axis is now a working discriminator, with data:
DSSYK$_\infty$ as measured here is a fake-temperature scrambler
(sine-dilaton's reading). Making it a *de Sitter* hologram requires either
a sign-flipping state/time map that no current proposal supplies (MSS
flat-space; HZ's fold is now excluded in-model), or a demonstration that
constrained dressing evades the chaos bound (NV; computable, unpublished).
The meson-tower axis (rung 18) and this axis are independent: a proposal
could pass one and fail the other. Current score on this axis:
sine-dilaton consistent, MSS flat-space in debt, NV pending a computable
sign.

*Caveat: NV and sine-dilaton statements above are quoted from their
papers, not re-derived here; the repo-side facts (Re a = −1, the fold
verdicts, β_fake·λ_L = 2π) are the tested content.*
