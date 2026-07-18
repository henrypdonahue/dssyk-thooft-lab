# Numerical test of DSSYK self-averaging (§3.5)

A brute-force exact-diagonalization test of the **self-averaging claims** in
Section 3.5 of Miyashita–Sekino–Susskind, *"Holograms and Standard Models"*
([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)) — the section that
underwrites the whole "which disorder instance is our universe?" story with a
diagrammatic, partly hand-waved argument. This is the *falsifiable* corner of
the paper: either the suppression is as advertised, or it is weaker, which would
crack Section 3.

## The two claims

For the double-scaled SYK model at infinite Boltzmann temperature, with
$\langle W\rangle=\mathrm{Tr}(W)/\dim$ (Eq. 3.4):

| | claim | eq. |
|---|---|---|
| **Singlet self-averaging** | $\mathrm{Var}\langle W\rangle \sim p!/N^{p}$ for $O(N)$-invariant $W$ | 3.27 |
| **Symmetry violation** | non-singlet (adjoint) expectations $\le e^{-a\sqrt N}$ | 3.28 |

Model (Eq. 3.2, Majorana): $H=i^{p/2}\sum_{i_1<\dots<i_p}J_{i_1\dots i_p}\psi_{i_1}\cdots\psi_{i_p}$,
couplings $\mathrm{Var}(J)=p!/N^{p-1}$ (Eq. 3.8), double-scaling $\lambda=p^2/N$ fixed.

## The subtlety that shapes the whole test

**At fixed $p$, $p!/N^p$ is only a power law $N^{-p}$.** The paper's *striking*
super-exponential suppression ($e^{-aN}$, $e^{-a\sqrt N}$) appears only in the
**double-scaling limit** $p=\sqrt{\lambda N}$, where $p$ grows with $N$. Exact
diagonalization lives at $\dim=2^{N/2}$ and caps out near $N\sim18$–$20$, which
**cannot** reach large $p$ and large $N$ together — the same wall the headline
spectrum-match computation hits.

So the test is split by what each tool can honestly establish:

- **A. The mechanism, exactly & analytically.** The singlet second moment
  $m_2=\mathrm{Tr}(H^2)/\dim=\sum_c J_c^2$ has relative variance **exactly
  $2/\binom{N}{p}$** — no diagonalization needed. This is literally
  *one over the number of couplings*. It equals $2\,p!/N^p$ divided by
  $\prod_{k=1}^{p-1}(1-k/N)\approx e^{-p^2/2N}$ — i.e. times a factor
  $\approx e^{+p^2/2N}$: an $O(1/N)$ correction at *fixed*
  $p$, but an $O(\lambda)$ **order-one factor** in the double-scaling regime
  $\lambda=p^2/N$ that this section is actually about. Followed along $p=\sqrt N$
  to $N=640$ it drops to $\sim10^{-46}$ — super-exponential, the promised
  $p!/N^p\to e^{-aN}$ tail.
  *Caveat on what this proves:* $m_2=\sum_c J_c^2$ is exact, so $2/\binom{N}{p}$
  is an elementary statistical identity (sum of i.i.d. Gaussian squares) that
  holds **regardless of the paper's diagrammatic argument**. Matching it by ED is
  a code-correctness check ($\mathrm{Tr}\,H^2=\sum J^2$), not a test of the
  mechanism. The non-trivial tests are B and C.
- **B. Higher singlet moments ($m_4,m_6$) by ED.** Do they self-average at the
  same combinatorial rate? (Yes — same $\sim1/\binom{N}{p}$ trend, steepening
  slightly; these are the first genuinely non-tautological singlet checks.)
- **C. Symmetry violation by ED.** The off-diagonal of the $H$-dressed adjoint
  bilinear $B_{jk}=\mathrm{Tr}(\psi_j H\psi_k H)/\dim$ has zero ensemble mean for
  $j\ne k$; its typical size measures how much one instance breaks $O(N)$. We
  track its decay with $N$ and $p$.

## Results (summary)

- **Eq. 3.27 — mechanism confirmed and quantified (for the $H^{2k}$ singlets).**
  The singlet self-averaging *rate* is exactly $2/\binom{N}{p}$ (verified: an
  exact identity for $m_2$; ED-confirmed at 2–9%, consistent with the
  $\sqrt{2/(n-1)}\approx10\%$ sampling noise of $n=200$ realizations), which
  along the double-scaling line is super-exponential ($\sim10^{-46}$ at $N=640$;
  the *leading* form is $\sim e^{-\frac12\sqrt N\ln N}$ — subleading
  $O(\sqrt N)$ terms are sizable, e.g. a factor $\sim e^{-24}$ at $N=640$).
  Higher moments $m_4,m_6$ self-average at the same combinatorial rate. What is *not* established: a fully general bounded
  singlet observable (only powers of $H$ are probed), and the $m_2$ identity is
  not itself a test of the paper's diagrammatics (see caveat above).
- **Eq. 3.28 — consistent, but at the edge of ED's reach.** The trustworthy
  signal $\mathrm{rms\_off}$ of the adjoint bilinear decays with $N$ and steepens
  with $p$, as required. Distinguishing the precise $e^{-a\sqrt N}$ form from a
  steep power law needs the double-scaling regime ($p\sim\sqrt N$, large $N$)
  that ED cannot reach; the accessible window ($N\le18$) constrains the trend but
  cannot pin the exponent. (The dimensionless ratio $\mathrm{rms\_off}/\mathrm{rms\_diag}$
  has a spurious bump at $p=N/2$ where its denominator collapses — a metric
  artifact, not physics; read $\mathrm{rms\_off}$.) An honest reach limit,
  **not** a failure of the claim.

Full numbers in [`results.txt`](results.txt); figure in `self_averaging.png`.

## Validation (this is what makes the numbers trustworthy)

The construction is checked against the textbook signature before any physics is
measured — see `validate_syk.py`:

- **Clifford algebra & Hermiticity** exact to machine precision.
- **$[H,P]=0$** with fermion parity $P$, and the **parity-sector subtlety**:
  mixing the two blocks fakes Poisson statistics; resolving parity is mandatory.
- **Random-matrix Bott periodicity**: the mean gap ratio $\langle r\rangle$
  cycles GOE→GUE→GSE→GUE with $N\bmod 8$ — all three RMT classes reproduced
  across four $N\bmod8$ values ($N=10,14,18$→GUE; $12,20$→GSE; $16$→GOE).
  Note: the $N\bmod8$ table used here is the $p\equiv0\pmod4$ classification
  (guarded in `validate_syk.py`).

## Files

| file | purpose |
|------|---------|
| `syk.py` | sparse Majorana + complex SYK Hamiltonians (Jordan–Wigner), parity/charge operators |
| `validate_syk.py` | construction validation via $N\bmod8$ RMT level statistics |
| `self_averaging.py` | the §3.5 measurement (A/B/C above) + verdict |
| `plot_self_averaging.py` | the two-panel figure |
| `test_syk.py` | asserting `pytest` suite |
| `results.txt` | saved run output |
| `requirements.txt` | pinned deps |

## Usage

```bash
pip install -r requirements.txt
python3 validate_syk.py         # RMT Bott-periodicity validation (~30 s)
python3 self_averaging.py       # the full self-averaging test (~7 min)
python3 plot_self_averaging.py  # figure
pytest -q -m 'not slow'         # fast asserts
```

## Honest scope

This tests §3.5, a *supporting* claim, not the duality's headline (that's the
't Hooft spectrum match, which is theory-blocked on both routes). What this
module settles: the self-averaging **mechanism** is real and quantitatively
$1/\binom{N}{p}$; the super-exponential regime is a double-scaling statement that
ED confirms analytically for the leading moment and supports — without fully
resolving — for the subleading and symmetry-violating observables.
