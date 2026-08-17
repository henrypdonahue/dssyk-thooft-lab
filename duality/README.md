# The comparison dictionary

Every convention needed to compare a DSSYK∞ number against the 't Hooft
yardstick, encoded once and unit-tested — built **before** the DSSYK side
existed, so an index offset, factor of 2, or normalization has exactly one
place to hide.

| encoded fact | source |
|---|---|
| n_paper = n_FLZ + 2 (photon/graviton relabeling) | correction #1 |
| CP = (−1)^(n+1), consistent with FLZ parity through the offset | paper §4.2–4.3 |
| M₁ = −(p/2) iH (`m1_coefficient`) | correction #3, tested in `syk-self-averaging` |
| paper λ = p²/N **vs** chord λ = 2p²/N, q = exp(−λ_chord) | the factor-2 landmine |
| ᾱ = ḡ²N = p², ḡ² = λ, τ = p² mq² | paper §4.1–4.2 |

## The fingerprints (what a DSSYK tower must reproduce)

Built from `thooft-target/reference_spectrum.json`. Both are parameter-free.

1. **Mass-squared ratios** Mn²/Mη²: 1, 2.3794, 3.7285, 5.0892, …
2. **Interleave splittings**: one alternating-CP trajectory with near-unit
   gaps 1.0167, 0.9944, 1.0029, … → 1. Degenerate CP doublets would
   falsify the identification. (This also settles the paper's "two
   degenerate Regge trajectories" wording: equal slopes, interleaved —
   adjacent levels are *not* degenerate.)

## Also in this module

- [BOOST.md](BOOST.md) + `boost_mass.py` — boost-to-mass calibration.
  One mass gives a boost-frequency *continuum with nodes*, not a peak.
  The "tomperature" map is the stated open question.
- [LARGEQ.md](LARGEQ.md) + `largeq_anchor.py` — the large-q, λ → 0 anchor
  (conventions bridge 𝒥 = √2 p, ED-locked). The meson tower cancels at
  leading order; the anchor lives at O(1/q).
- [ANTISCRAMBLING.md](ANTISCRAMBLING.md) + `antiscrambling.py` — the OTOC
  sign test and the Euclidean fold. DSSYK∞ scrambles at every temperature
  (Re a = −1 exactly), so the flat dictionary cannot give dS
  anti-scrambling. The fold does not repair it: a genuine N-obstruction at
  ED sizes (`../syk-self-averaging/fold_bench.py --scan`) and a
  λ-obstruction at N = ∞ (`chord_fold.py`).
- [CHORD.md](CHORD.md) + `chord.py`, `chord_moments.py`, `chord_fold.py` —
  **the exact N = ∞ chord engine**: transfer matrix plus matter-chord
  contours, validated four independent ways (enumerator, closed forms,
  Krawtchouk microscopy, finite-N ED to 0.6–1.6%). Results: structureless
  spectra at every λ; the singlet channels carry zero propagating weight
  at N = ∞; the fold fails in the double-scaled theory itself.
  Figure: `plot_chord.py`.
- [DISCRIMINATOR.md](DISCRIMINATOR.md) — where the three competing duals
  (MSS flat-space, Narovlansky–Verlinde, sine-dilaton) sit against the
  chaos-bound no-go and this repo's data. Sine-dilaton's black-hole
  reading is qualitatively consistent with the data; its quantitative
  θ ↔ v tie is untested.
- `signflip_probe.py` — the first candidate dS map that works: the
  antipodal half-fake-period shift flips the OTOC sign (exact in the
  closed form; exact theory: λ ≤ 1.7) and makes the 2-pt sector exactly
  thermal at the tomperature. Necessary, not sufficient; it is NV's
  dressing, geometrically.

```bash
python3 dictionary.py        # print the yardstick in both conventions
python3 antiscrambling.py    # the OTOC sign table + fold map
python3 chord.py --spectra   # the lambda->0 spectral scan (~4 min)
python3 chord_moments.py --run  # exact N = inf moments    (~1 s)
python3 chord_fold.py        # the fold at N = inf         (~30 s)
python3 signflip_probe.py    # the antipodal map            (~1 min)
python3 signflip_probe.py --deep  # its deep campaign       (~20 s)
pytest -q                    # 165 tests
```
