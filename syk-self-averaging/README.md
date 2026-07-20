# DSSYK self-averaging, the $M_n$ tower, and the QED sector

Tests of the *falsifiable* corner of Miyashita–Sekino–Susskind
([arXiv:2607.05678](https://arxiv.org/abs/2607.05678)): §3.5's self-averaging
claims, plus the operator content of the U(N)/Dirac sector its "standard model"
story needs. Model (Eq. 3.2): $H=i^{p/2}\sum J_{i_1\dots i_p}\psi_{i_1}\cdots\psi_{i_p}$,
$\mathrm{Var}(J)=p!/N^{p-1}$ (Eq. 3.8), $\lambda=p^2/N$, $\langle W\rangle=\mathrm{Tr}(W)/\dim$.

| claim | eq. | verdict |
|---|---|---|
| singlet $\mathrm{Var}\langle W\rangle\sim p!/N^p$ | 3.27 | confirmed, rate exactly $2/\binom{N}{p}$ |
| adjoint violations $\le e^{-a\sqrt N}$ | 3.28 | **verified exactly in double scaling** (`exact_wick.py`) |

## Results

**Eq. 3.27.** The singlet rate is exactly $2/\binom{N}{p}$ — one over the number
of couplings; equals the paper's $p!/N^p$ times $2\prod(1-k/N)^{-1}$, an $O(1)$
factor in double scaling. Along $p=\sqrt N$: $\sim10^{-46}$ by $N=640$ (leading
form $e^{-\frac12\sqrt N\ln N}$; subleading terms are sizable). *Caveats:* the
$m_2$ identity holds regardless of the paper's diagrammatics (its ED match is a
code check); the non-tautological content is $m_4,m_6$ self-averaging at the
same rate (confirmed by ED at $N\le18$, within the ~10% noise of 200 draws),
and only $H^{2k}$ singlets are probed.

**Eq. 3.28 — settled beyond ED.** $B_{jk}=\mathrm{Tr}(\psi_jH\psi_kH)/\dim$ is
quadratic in the couplings, so Wick gives closed forms at any $(N,p)$:
$$\mathrm{Var}(B_{jk})=4\sigma^4\binom{N-2}{p-1},\qquad
  \mathbb{E}[B_{jj}]=\sigma^2\binom{N}{p}\Big(1-\tfrac{2p}{N}\Big),$$
verified against literal contraction enumeration (exact) and ED (~1%). Along
$p=\sqrt N$: $\ln\mathrm{rms}=-\tfrac14\sqrt N\ln N\,(1+o(1))$ — beats
$e^{-a\sqrt N}$ for every $a$, verified in the regime the claim is about
(scanned to $N=1600$). At fixed $p$ (the match regime) the asymptote is the
power law $N^{-(p-1)/2}$. $\mathbb{E}[B_{jj}]=0$ at $p=N/2$ *derives* the
ratio artifact flagged below.

**$M_n$ tower and CP, measured** (`dirac.py`, `mn_spectroscopy.py`). Exact and
asserted: $[H,Q]=0$; $M_0=Q$; $M_1=-(p/2)iH$ (paper's $M_1=H$ fails at $O(1)$);
$M_2=-\sum\dot c^\dagger\dot c$ (paper's $+$ sign is wrong);
$M_n^\dagger=(-1)^nM_n$ only for $n\le2$. The paper's CP must be the **unitary**
particle-hole map, $C c_i C^\dagger=(-1)^{N_c-1}c_i^\dagger$ ($M_1\propto iH$
rules out antiunitary). CP$(M_0)=-1$, CP$(M_1)=+1$ exact on C-invariant
instances; raw $M_{n\ge2}$ are CP-*mixtures* (15–68% wrong channel), so
CP$=(-1)^{n+1}$ lives in CP-resolved correlators — `mn_spectroscopy.py`
computes them (first data: conserved $n=0,1$ machine-exact at $\omega=0$;
broad continua at $\lambda=1.6$; not the $\lambda\to0$ match regime).

Numbers: [`results.txt`](results.txt)/`results.json`; figures
`self_averaging.png`, `mn_spectra.png`.

## Why the numbers are trustworthy

- Instance-level ED caps at $N\sim18$–20 ($\dim 2^{N/2}$) and cannot reach the
  double-scaling regime — every double-scaling statement here is exact
  combinatorics, not extrapolation; ED-window fits carry error bars and ΔAIC
  verdicts that say when models are indistinguishable.
- Construction validated before physics: Clifford algebra and $[H,P]=0$ to
  machine precision; $N\bmod8$ RMT Bott periodicity reproduced (GUE/GSE/GOE
  across $N=10..20$; parity-sector resolution is mandatory, $p\equiv0\bmod4$
  table guarded).
- Fast bitwise builders are equivalence-tested against independent sparse
  constructions; read `rms_off`, not the $p\sim N/2$-artifacted ratio.

## Files

| file | purpose |
|------|---------|
| `syk.py` | sparse JW operators — slow reference builders |
| `pauli_strings.py` | bitwise Hamiltonian assembly (10–100×, equivalence-tested) |
| `dirac.py` | complex SYK: $[H,Q]=0$, $M_n$ tower, unitary CP machinery |
| `exact_wick.py` | closed-form Eq. 3.28 statistics, beyond ED |
| `mn_spectroscopy.py` | CP-resolved $M_n$ spectral functions |
| `validate_syk.py` | RMT level-statistics validation |
| `self_averaging.py` | the §3.5 measurement + verdict |
| `plot_self_averaging.py` | figure (reads `results.json`) |
| `test_syk.py` | 24-test asserting suite |

## Usage

```bash
pip install -r requirements.txt
python3 validate_syk.py         # RMT validation (~7 s)
python3 self_averaging.py       # §3.5 measurement (~2 min; writes results.*)
python3 exact_wick.py           # exact Eq. 3.28 + double-scaling scan (~1 s)
python3 dirac.py                # Dirac/QED identities + CP report (~5 s)
python3 mn_spectroscopy.py      # Mn spectral functions (~1 min; --big = Nc 12)
python3 plot_self_averaging.py  # figure
pytest -q -m 'not slow'         # fast asserts (~4 s)
```

## Scope

This settles §3.5 (both equations, one exactly) and the operator-level
prerequisites of the U(N) story. It does not touch the headline spectrum match;
the CP-resolved spectroscopy is the pipeline for it, pending the boost-to-mass
dictionary and the $\lambda\to0$ regime (road_map rungs 13, 15–18).
