# papers/ — local reference copies, not tracked

This directory holds every paper the repository cites, plus two
DSSYK–de Sitter background papers (2206.01083, 2303.00792). The PDFs are
gitignored, because the arXiv default license does not permit
redistributing third-party papers, so a clone gets this README and
`fetch.sh` only. To rebuild the folder:

```bash
bash papers/fetch.sh      # downloads the PDFs below and the target paper's LaTeX source
```

Line-number citations in the code, for example in `qed_campaign.py` and
`baryons.py`, refer to `arXiv-2607.05678v1/HologramsSM.tex`.

## The proposal under test

| arXiv | Authors | Title | Why it is here |
|---|---|---|---|
| [2607.05678](https://arxiv.org/abs/2607.05678) | Miyashita, Sekino, Susskind | Holograms and Standard Models | The target paper (v1) |
| [2506.18054](https://arxiv.org/abs/2506.18054) | Miyashita, Sekino, Susskind | DSSYK at infinite temperature: the flat-space limit and the 't Hooft model | The companion: the λ = ḡ² dictionary, the flat- and 't Hooft-limit definitions, the Sₙ tower (v3) |
| [2501.09423](https://arxiv.org/abs/2501.09423) | Sekino, Susskind | Double-scaled SYK, QCD, and the flat space limit of de Sitter space | The self-averaging counting the target defers to; λ = 2p²/N there |
| [2012.12326](https://arxiv.org/abs/2012.12326) | Susskind | Electric Forces in the Charged SYK Model | The U(1) Coulomb-force picture behind M₀ = Q |

## DSSYK and de Sitter holography

| arXiv | Authors | Title | Why it is here |
|---|---|---|---|
| [2206.01083](https://arxiv.org/abs/2206.01083) | Lin, Susskind | Infinite Temperature's Not So Hot | Tomperature; DSSYK∞ as the dS static patch |
| [2303.00792](https://arxiv.org/abs/2303.00792) | Susskind | De Sitter Space has no Chords. Almost Everything is Confined. | Confinement to the stretched horizon; quasinormal-mode time scales |
| [2401.08555](https://arxiv.org/abs/2401.08555) | Rahman, Susskind | Infinite Temperature is Not So Infinite: The Many Temperatures of de Sitter Space | Which scales stay finite in the semiclassical limit |
| [2407.12988](https://arxiv.org/abs/2407.12988) | Rahman, Susskind | p-Chords, Wee-Chords, and de Sitter Space | Argues the NV and DSSYK chord algebras conflict |
| [2310.16994](https://arxiv.org/abs/2310.16994) | Narovlansky, Verlinde | Double-scaled SYK and de Sitter holography | The dS₃ dual; m² = 4Δ(1−Δ); the ±iβ/4 dressing; dS quasinormal modes |
| [2404.03535](https://arxiv.org/abs/2404.03535) | Blommaert, Mertens, Papalini | The dilaton gravity hologram of double-scaled SYK | The sine-dilaton dual; fake temperature |
| [2505.08116](https://arxiv.org/abs/2505.08116) | Okuyama | de Sitter JT gravity from double-scaled SYK | A further route to dS JT |
| [2607.13665](https://arxiv.org/abs/2607.13665) | Cui, Kolchmeyer | A de Sitter Anti-Scrambling Algebra | The anti-scrambling requirement (`duality/ANTISCRAMBLING.md`) |
| [2607.14215](https://arxiv.org/abs/2607.14215) | Harlow, Zhao | Anti-scrambling and euclidean folds from observer correlators in de Sitter space | The Euclidean fold; the notes cite v2 §6 (v1 and v2 kept) |
| [1904.01911](https://arxiv.org/abs/1904.01911) | Maldacena, Turiaci, Yang | Two dimensional Nearly de Sitter gravity | Near-Nariai reduction (`duality/jt_ds_vacuum.py`) |
| [2203.00700](https://arxiv.org/abs/2203.00700) | Svesko, Verheijden, Verlinde, Visser | Quasi-local energy and microcanonical entropy in two-dimensional nearly de Sitter gravity | The JT-dS first-law sign |

## Double-scaled SYK: chords and limits

| arXiv | Authors | Title | Why it is here |
|---|---|---|---|
| [1811.02584](https://arxiv.org/abs/1811.02584) | Berkooz, Isachenkov, Narovlansky, Torrents | Towards a full solution of the large N double-scaled SYK model | The chord solution the engines implement |
| [2407.09396](https://arxiv.org/abs/2407.09396) | Berkooz, Mamroud | A cordial introduction to double scaled SYK | The review whose conventions the code follows (Eq. 2.30 vertex) |
| [2006.13983](https://arxiv.org/abs/2006.13983) | Berkooz, Narovlansky, Raj | Complex Sachdev-Ye-Kitaev model in the double scaling limit | Dirac chords (`duality/dirac_chord.py`) |
| [2403.05980](https://arxiv.org/abs/2403.05980) | Berkooz, Brukner, Jia, Mamroud | Path integral for chord diagrams and chaotic-integrable transitions in double scaled SYK | The chord path integral |
| [2403.01950](https://arxiv.org/abs/2403.01950) | Berkooz, Brukner, Jia, Mamroud | From Chaos to Integrability in Double Scaled Sachdev-Ye-Kitaev Model via a Chord Path Integral | The small-λ transition |
| [2601.09801](https://arxiv.org/abs/2601.09801) | Aguilar-Gutierrez, Das, Erdmenger, Xian | Probing the chaos to integrability transition in double-scaled SYK | The same transition, probed |
| [2302.11882](https://arxiv.org/abs/2302.11882) | Narayan, S | SYK Model with global symmetries in the double scaling limit | U(M)-symmetric chords (`duality/dirac_chord.py`) |
| [2512.07715](https://arxiv.org/abs/2512.07715) | Forste, Kruse, Natu | Grand canonical vs. canonical Krylov complexity in double-scaled complex SYK model | Fixed-charge Dirac chord conventions |

## SYK at large interaction order

| arXiv | Authors | Title | Why it is here |
|---|---|---|---|
| [1604.07818](https://arxiv.org/abs/1604.07818) | Maldacena, Stanford | Remarks on the Sachdev-Ye-Kitaev model | The bilinear spectrum and the finite-T ladder kernel |
| [1702.08016](https://arxiv.org/abs/1702.08016) | Gross, Rosenhaus | The Bulk Dual of SYK: Cubic Couplings | Large-q OPE coefficients, O(1/q) against O(1/q²) |
| [1911.10171](https://arxiv.org/abs/1911.10171) | Streicher | SYK Correlators for All Energies | The large-q OTOC closed forms (gauge: his App. A) |
| [1912.00004](https://arxiv.org/abs/1912.00004) | Choi, Mezei, Sárosi | Exact four point function for large q SYK from Regge theory | The 𝓜ᵉ/𝓜ᵒ quantization sets |
| [1610.03816](https://arxiv.org/abs/1610.03816) | García-García, Verbaarschot | Spectral and thermodynamic properties of the Sachdev-Ye-Kitaev model | The N mod 8 random-matrix classes |

## The 't Hooft model

| arXiv | Authors | Title | Why it is here |
|---|---|---|---|
| [0905.2280](https://arxiv.org/abs/0905.2280) | Fateev, Lukyanov, Zamolodchikov | On mass spectrum in 't Hooft's 2D model of mesons | The α = 0 tables, sum rules and d± the target is checked against |
| [2409.11324](https://arxiv.org/abs/2409.11324) | Litvinov, Meshcheriakov | Meson mass spectrum in QCD2 't Hooft's model | The α ≠ 0 anchors (`external_anchors.json`) |
| [2312.15598](https://arxiv.org/abs/2312.15598) | Ambrosino, Komatsu | 2d QCD and integrability. Part I. 't Hooft model | The TQ-Baxter recasting |
| [2406.11078](https://arxiv.org/abs/2406.11078) | Ambrosino, Komatsu | 2d QCD and Integrability, Part II: Generalized QCD | Its continuation |

## Versions

As of 2026-10-14: 2607.05678 is still v1; 2607.14215 has a v2 (minor
corrections; the technical note cites v2 §6); 2506.18054 is at v3, which
corrects its §5 statements about the 't Hooft and flat-space limits and
redefines ḡ² = g²/m_q² (its footnote 9). Titles and authors were checked
against INSPIRE on 2026-10-14.

`EMAIL-draft-local.md` is a local, gitignored outreach draft (stale).
