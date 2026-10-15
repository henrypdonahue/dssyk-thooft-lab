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
curl -sL -o 2607.14215v2.pdf https://arxiv.org/pdf/2607.14215v2   # its v2 (2026-07-31), the version the notes cite
curl -sL https://arxiv.org/e-print/2607.05678v1 -o src.tar.gz     # the target paper's LaTeX source
mkdir -p arXiv-2607.05678v1 && tar -xzf src.tar.gz -C arXiv-2607.05678v1 && rm src.tar.gz
```

Versions as of 2026-10-04 (2607.05678 re-checked 2026-10-14): 2607.05678 is still v1; 2607.14215 has a v2
(minor corrections; the technical note cites v2 §6); the companion paper
2506.18054 is at v3, which corrects its Section 5 statements about the
't Hooft and flat-space limits.

Line-number citations in the code, for example in `qed_campaign.py` and
`baryons.py`, refer to `arXiv-2607.05678v1/HologramsSM.tex`.
