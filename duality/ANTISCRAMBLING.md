# Anti-scrambling vs DSSYK∞: the sign test and the Euclidean fold

*Rungs 25–26 (2026-07-23). Code: `antiscrambling.py` (closed forms, verified
transcription), `syk-self-averaging/fold_bench.py` (ED laboratory,
`fold_bench.json`), asserted in `test_antiscrambling.py`.*

## Why this exists

Two papers of July 2026 turned the OTOC sector into a sharp external
constraint on the duality this repo stress-tests:

- **Cui–Kolchmeyer** ([2607.13665](https://arxiv.org/abs/2607.13665)): dS
  shockwaves at the cosmological horizon have a **negative** eikonal phase
  (time *advance*). Writing the semiclassical OTOC as
  $\mathrm{OTOC/TOC} = 1 - i\,c\,e^{T} h_Ah_B e^{\cdots}/(\sinh\cdot\sinh)$
  (their Eq. 2.74), AdS has $c>0$ (scrambling), dS needs $c<0$
  (*anti-scrambling*). The MSS chaos bound forces $c>0$ for any unitary QM
  in a KMS state at the dS temperature — dS violates the KMS assumption, the
  Hartle–Hawking state fails to be a trace on the observer's crossed-product
  algebra, and observer-centric static patch holography is challenged.
  2607.05678 (this repo's target) is their ref. [12]. Its $T_B=\infty$
  hologram is **not** a KMS state at the dS temperature, so it dodges the
  letter of the no-go — but then *something* in the dictionary must produce
  $c<0$.
- **Harlow–Zhao** ([2607.14215](https://arxiv.org/abs/2607.14215)): same
  bulk physics (observer 2-pt thermal/KMS, 4-pt anti-scrambles,
  $F(w)\to F(-w)$), plus a proposed QM mechanism: a system with spectrum
  **bounded above and below**, evaluated on a **Euclidean fold** — their
  continuation (6.9), $\tau^{QM} = (-\tau_1,\,2\pi n-\tau_2,\,2\pi-\tau_3,\,
  2\pi(n{+}1)-\tau_4)$, which flips $w\to-w$. They pose the test verbatim:
  *"It would be interesting to study the validity of this assumption in a
  concrete model such as the SYK model."* Finite-$N$ SYK is bounded both
  ways; every folded correlator is exactly computable by ED. That test is
  what `fold_bench.py` runs.

## Result 1 (rung 25, settled): DSSYK∞ scrambles — the flat dictionary cannot anti-scramble

Streicher's exact large-$q$ closed forms (all temperatures, transcription
from the arXiv LaTeX source, conventions bridge $\mathcal J=\sqrt2\,p$
ED-locked) give, at the symmetric crossed configuration
$\theta=(3\pi/2,\pi/2,\pi,0)$ with Lorentzian $\theta$-time on the
$(\theta_1,\theta_2)$ pair:

$$\mathcal F/F_d \sim a(v)\,e^{v\,\theta_L},\qquad
  a(v) = -\frac{e^{-i\pi v/2}}{\cos(\pi v/2)},\qquad
  \boxed{\ \mathrm{Re}\,a(v) = -1\ \text{ exactly, at every temperature}\ }$$

with rate $\lambda_L = 2\pi v/\beta \to 2\mathcal J$ as $T_B\to\infty$: the
scrambling (AdS/time-delay, CK's $c>0$) sign, **unsuppressed at the
$T_B=\infty$ hologram point**.

What ED adds, and what it does not. ED independently confirms the OTOC
*magnitude growth* at all orders in $1/N,1/q$: $\mathrm{Re}\,\mathcal F/F_d$
moves monotonically away from its time-ordered value at every measured
$(N,\beta\mathcal J)$, by well over $100\sigma$ (drop/err $\approx130$–$200$;
`test_bench_scrambling_sign_verdict`). The *sign* that makes this
"scrambling" rather than "anti-scrambling," however, is **analytic**
(Streicher-transcribed): it rides on the crossed-operator ordering, which is
fixed empirically to match Streicher (the opposite order is off by ×25; see
caveats). So ED gives an independent measurement of the growth, not an
independent measurement of the sign — do not read this as ED alone deciding
scrambling vs. anti-scrambling.

**Consequence.** Under the flat identification (boundary time = observer
proper time), DSSYK∞ has the wrong OTOC sign to be CK's dS observer — at
any temperature, including $T_B=\infty$. This is falsification pressure on
assumption #1 of 2607.05678 *unless* the dictionary contains a nontrivial
time/state map (fold, fake-disk complexification, tomperature) that flips
the sign. It sharpens rung 13's "tomperature map" gap from "unstated" to
"must anti-scramble."

## Result 2 (rung 26, first data): the fold works mechanically; whether it commutes with large N is genuinely open

`fold_bench.py` evaluates the HZ fold (6.9, $n=1$) of the same crossed
configuration: the ordered-trace correlator is an **entire function** of the
four $\theta$'s at finite dimension, so the continuation is unambiguous —
freeze the operator sequence at its standard-configuration path order, move
the arguments (this produces the $e^{+\tau H}$ segments; finite because the
SYK spectrum is bounded — exactly HZ's requirement). Findings at
$p=6,\ \beta\mathcal J\in\{2,0.5\},\ N\in\{12,14,16\}$:

1. **Mechanics work.** Every folded correlator (2-pt and 4-pt) is finite and
   well-behaved; the folded 2-pt rises above 1 (Euclidean growth), the
   qualitative signature of the continued closed form.
2. **The Lorentzian fold response tracks the continuation.** The unfolded
   ED OTOC is exactly real (configuration symmetry); the *folded* ED
   correlator develops an $\mathrm{Im}\,\mathcal F/F_d$ that is nonzero,
   negative, and monotone — matching the analytically continued closed form
   to ~10% at the best point ($N{=}12$, $\beta\mathcal J{=}2$) and in
   sign/trend everywhere. The fold breaks the reality symmetry *exactly the
   way the continuation predicts*.
3. **The open question, made quantitative.** Unfolded observables converge
   toward the closed forms as $N$ grows at fixed $p$ (2-pt deviation ratios
   $0.35\to0.44\to0.52$ for $N=12\to14\to16$). The *folded* Euclidean 4-pt
   starts near the continued form at $N{=}12$ (ratio $\approx1.05$) but
   drifts **away** as $N$ grows ($\approx1.21$ at $N{=}14$, $\approx1.42$ at
   $N{=}16$); near the continuation's poles (the folded 2-pt at
   $\theta=-\pi$, the small-$v$ folded 4-pt) the ED values converge slowly
   or not visibly. Two readings, not distinguishable at ED sizes: (a)
   $\lambda=p^2/N$ corrections simply behave differently under folds (still
   $\lambda\ge2.25$ here), or (b) a genuine Stokes-type obstruction to HZ's
   "continuation commutes with the semiclassical limit" assumption is
   developing. **This is the first data on their SYK question either way.**

## What would decide it

- A fixed-$\lambda$ scan ($p=4$, $N=8..24$, sector-resolved ED) separates
  $\lambda$-corrections from genuine $N$-trends.
- Chord/exact methods at $N=\infty$ (the rung-18 toolchain) could evaluate
  the folded 4-pt directly in the double-scaled theory — the fold is a
  contour statement, and DSSYK chord correlators are analytic in the chord
  angles. That would answer HZ exactly, and is adjacent to the headline
  tower extraction.
- The negative-temperature variant (HZ 6.12) is the $H\to-H$ reflection
  (tested as such): for the SYK *ensemble*, $H$ and $-H$ are equidistributed
  at $p\equiv0\bmod4$, so disorder-averaged (6.9) and (6.12) coincide —
  the two prescriptions are ensemble-equivalent here, one less ambiguity.

## Standing caveats (stated once, apply throughout)

Large-$q$ closed forms at ED sizes carry an $O(1/p,\lambda)$ floor measured
in-line (2-pt deviation ratios ~0.3–0.55, TOC ~0.5–0.8 of the closed form):
all 4-pt comparisons are sign/trend-level, each read against the same-size
unfolded calibration. The connected ratio uses
$N(\langle F\rangle - \langle G\rangle\langle G\rangle)$, which includes the
disorder-replica piece (that is what the ladder resums). The crossed-path
sign convention is fixed empirically (the alternative is off by ×25) and the
fold configuration in the JSON is asserted equal to `fold_map` of the
crossed configuration, tying ED and closed forms together.
