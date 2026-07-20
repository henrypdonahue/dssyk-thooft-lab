# The large-q anchor: what the λ→0 corner actually contains (rung 17)

`largeq_anchor.py` encodes the equation-verified large-q SYK closed forms
(Maldacena–Stanford 1604.07818; Streicher 1911.10171; Choi–Mezei–Sárosi
1912.00004; Gross–Rosenhaus 1702.08016) with the exact conventions bridge to this
repo ($\mathcal J = \sqrt2\,p$), validated against the ED bench
(`syk-self-averaging/largeq_bench.json`).

## The structural finding

The large-q singlet 4-point function is known in closed form at **all**
temperatures, and its Sommerfeld–Watson structure is meromorphic with poles only at
$m=0$ (Hamiltonian exchange — the *entire* time-ordered correlator, at every
temperature) and $m=\pm v$ (the scramblon, OTOC only). The would-be bilinear tower
— the quantization sets $\mathcal M^e/\mathcal M^o$ (CMS 3.3), the all-temperature
analog of $k(h)=1$ — **cancels identically at leading order**, and its couplings
are $c_n^2 = O(1/q^2)$ (GR 2.15). At strict $q=\infty$, $k_c(h)=2/[h(h-1)]$ has
only $h=2$.

**Consequence: the λ→0 anchor of the 't Hooft match cannot come from the leading
large-q 4-point function.** The meson tower first appears at $O(1/q)$, organized on
the exact-in-$\beta\mathcal J$ skeleton $\mathcal M^e/\mathcal M^o$ — whose
finite-temperature residues are unpublished. Computing them **is** rung 17 proper
(and no longer looks optional: it is where the entire quantitative content of the
duality's $\lambda\to0$ corner lives).

## Encouraging structure

- $\mathcal M^e/\mathcal M^o$ **interleave** — even/odd alternation, exactly the
  CP-alternation slot of the 't Hooft tower (`duality/dictionary.py`).
- At $T=\infty$ ($v\to0$) the sets degenerate to integers in $2\pi/\beta$ units —
  physical-frequency spacing diverges as $\beta\to0$, consistent with `BOOST.md`:
  masses are not boost-frequency peak positions.
- The ED bench confirms the transcription end-to-end: $\langle H^2\rangle$ matches
  the bridge exactly, and $G(t)=\mathrm{sech}(\sqrt2\,p\,t)^{2/p}$ matches the ED
  curves with no free parameters at the $O(1/p)$ level (measured deviations
  0.07–0.10 at $p=4$ vs $1/p=0.25$; the deviation *plateaus* rather than
  vanishing as $\lambda\to0$ — the intrinsic $1/p$ error survives, and at
  $\lambda\sim1$ the $\lambda$- and $1/p$-corrections partially cancel).

## En-route finding

The "obvious" Majorana singlet $\tfrac iN\sum\psi_a\dot\psi_a$ is exactly conserved:
$\sum_a\psi_a[H,\psi_a]=-2pH$ (the Majorana counterpart of $M_1=-(p/2)iH$; asserted
in the tests). The first dynamical Majorana tower member is $\sum\dot\psi\dot\psi$
— used by the bench.
