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

Also in this module: [`BOOST.md`](BOOST.md) + `boost_mass.py` — the rung-13
boost-to-mass calibration (one mass = a boost-frequency *continuum with
nodes*; the "tomperature" map is the stated open question);
[`LARGEQ.md`](LARGEQ.md) + `largeq_anchor.py` — the rung-17 large-q
$\lambda\to0$ anchor (conventions bridge $\mathcal J=\sqrt2\,p$, ED-locked;
the meson tower cancels at leading order, so the anchor lives at $O(1/q)$);
and [`ANTISCRAMBLING.md`](ANTISCRAMBLING.md) + `antiscrambling.py` — the
rung-25/26 response to Cui–Kolchmeyer 2607.13665 and Harlow–Zhao 2607.14215:
DSSYK∞ *scrambles* at every temperature ($\mathrm{Re}\,a=-1$ exactly), so
the flat dictionary cannot produce dS anti-scrambling, and the first SYK
data on the HZ Euclidean-fold prescription (ED lab:
`../syk-self-averaging/fold_bench.py`).

```bash
python3 dictionary.py     # print the yardstick in both conventions
python3 antiscrambling.py # the OTOC sign table + fold map
pytest -q                 # 38 tests: conventions, boost, large-q, anti-scrambling
```
