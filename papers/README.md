# papers/ — local reference copies (not tracked)

The papers this repo engages with, kept here for convenience. The files
are gitignored: the arXiv default license does not permit third-party
redistribution, so a clone gets this README only. Fetch them with:

```bash
cd papers
curl -sLO https://arxiv.org/pdf/2607.05678v1   # Miyashita-Sekino-Susskind, the target paper
curl -sLO https://arxiv.org/pdf/2607.13665v1   # Cui-Kolchmeyer, anti-scrambling algebra
curl -sLO https://arxiv.org/pdf/2607.14215v1   # Harlow-Zhao, anti-scrambling + Euclidean folds
curl -sL  https://arxiv.org/e-print/2607.05678v1 -o src.tar.gz   # target-paper LaTeX source
mkdir -p arXiv-2607.05678v1 && tar -xzf src.tar.gz -C arXiv-2607.05678v1 && rm src.tar.gz
```

Line-number citations in the code (e.g. `qed_campaign.py`, `baryons.py`)
refer to `arXiv-2607.05678v1/HologramsSM.tex`.
