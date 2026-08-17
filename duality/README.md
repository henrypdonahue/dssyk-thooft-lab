# The comparison dictionary

Every convention needed to compare a DSSYK∞ number against the 't Hooft
yardstick, encoded once, unit-tested, and built **before** the DSSYK side of the
comparison exists — so that when a number finally appears there is exactly one
place an index offset, factor of 2, or normalization can hide.

| encoded fact | source |
|---|---|
| $n_{\text{paper}} = n_{\text{FLZ}} + 2$ (photon/graviton relabeling) | correction #1 |
| CP $=(-1)^{n+1}$, consistent with FLZ wavefunction parity through the offset | paper §4.2–4.3 |
| $M_1 = -(p/2)\,iH$ (`m1_coefficient`) | correction #3, tested in `syk-self-averaging` |
| paper $\lambda = p^2/N$ **vs** chord $\lambda_{\text{chord}} = 2p^2/N$, $q=e^{-\lambda_{\text{chord}}}$ | the factor-2 landmine |
| $\bar\alpha = \bar g^2 N = p^2$, $\bar g^2 = \lambda$, $\tau = p^2 m_q^2$ | paper §4.1–4.2 |

## The fingerprints (what a DSSYK tower must reproduce)

Built from `thooft-target/reference_spectrum.json`; both are **parameter-free**
(no units, no fitted constants):

1. **Mass-squared ratios** $M_n^2/M_\eta^2$: 1, 2.3794, 3.7285, 5.0892, …
2. **Interleave splittings**: a *single* alternating-CP trajectory with
   near-unit gaps 1.0167, 0.9944, 1.0029, … → 1. Pairwise-degenerate CP
   doublets would falsify the identification. (This also fixes the reading of
   the paper's "two degenerate Regge trajectories… on top of one another":
   equal slopes and interleaving — adjacent levels are *not* degenerate.)

Also in this module:

- [`BOOST.md`](BOOST.md) + `boost_mass.py` — rung-13 boost-to-mass
  calibration: one mass = a boost-frequency *continuum with nodes*; the
  "tomperature" map is the stated open question.
- [`LARGEQ.md`](LARGEQ.md) + `largeq_anchor.py` — rung-17 large-q
  $\lambda\to0$ anchor (conventions bridge $\mathcal J=\sqrt2\,p$, ED-locked);
  the meson tower cancels at leading order, so the anchor lives at $O(1/q)$.
- [`ANTISCRAMBLING.md`](ANTISCRAMBLING.md) + `antiscrambling.py` —
  rung-25/26 response to Cui–Kolchmeyer 2607.13665 and Harlow–Zhao
  2607.14215: DSSYK∞ *scrambles* at every temperature
  ($\mathrm{Re}\,a=-1$ exactly), so the flat dictionary cannot produce dS
  anti-scrambling; plus the first SYK data on the HZ Euclidean-fold
  prescription and the fixed-λ scan verdict — the folded Euclidean drift
  is a genuine $N$-obstruction, not a λ-artifact (ED lab:
  `../syk-self-averaging/fold_bench.py`, `--scan`).
- [`CHORD.md`](CHORD.md) + `chord.py`, `chord_moments.py`, `chord_fold.py`
  — **the exact $N=\infty$ chord engine** (rung-18 toolchain): transfer
  matrix + matter-chord contour evaluator, validated against a brute-force
  chord-diagram enumerator, the published closed forms, a chord-free
  Krawtchouk microscopy, and finite-$N$ ED (thermal 2-pt to 0.6–1.6%).
  Deliverables: the rung-18 λ→0 spectral scan (**structureless continua;
  and the $A_n$ singlet channels are exactly conserved at $N=\infty$** —
  the tower is invisible at every order of the double-scaled expansion;
  `chord_spectra.json`, `chord_moments.json`), and the rung-26 decider #2
  (**the HZ fold fails in the double-scaled theory itself**;
  `chord_fold.json`). Figure: `plot_chord.py`.
- [`DISCRIMINATOR.md`](DISCRIMINATOR.md) — rung-24 KMS/OTOC-axis note:
  where the three competing duals (MSS flat-space, Narovlansky–Verlinde,
  sine-dilaton) sit against the CK no-go and this repo's measured data;
  sine-dilaton's fake-temperature reading is the one consistent with the
  data as-is ($\beta_{\rm fake}\lambda_L=2\pi$ tested).

```bash
python3 dictionary.py     # print the yardstick in both conventions
python3 antiscrambling.py # the OTOC sign table + fold map
python3 chord.py --spectra   # rung-18 lambda->0 scan  (~4 min)
python3 chord_moments.py --run  # rung-15 exact moments (~1 s)
python3 chord_fold.py        # rung-26 decider #2       (~30 s)
pytest -q                 # 151 tests: conventions, boost, large-q,
                          # anti-scrambling, chord engine + verdicts
```
