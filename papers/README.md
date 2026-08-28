# papers/ — local reference copies, not tracked

This directory holds the papers the repo engages with. The files are
gitignored. The arXiv default license does not permit third-party
redistribution, so a clone gets this README and nothing else.

Fetch them like this:

```bash
cd papers
curl -sL -o 2607.05678v1.pdf https://arxiv.org/pdf/2607.05678v1   # Miyashita-Sekino-Susskind, the target paper
curl -sL -o 2607.13665v1.pdf https://arxiv.org/pdf/2607.13665v1   # Cui-Kolchmeyer, the anti-scrambling algebra
curl -sL -o 2607.14215v1.pdf https://arxiv.org/pdf/2607.14215v1   # Harlow-Zhao, anti-scrambling and Euclidean folds
curl -sL https://arxiv.org/e-print/2607.05678v1 -o src.tar.gz     # the target paper's LaTeX source
mkdir -p arXiv-2607.05678v1 && tar -xzf src.tar.gz -C arXiv-2607.05678v1 && rm src.tar.gz
```

Line-number citations in the code, for example in `qed_campaign.py` and
`baryons.py`, refer to `arXiv-2607.05678v1/HologramsSM.tex`.
