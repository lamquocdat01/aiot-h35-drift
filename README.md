# aiot-h35-drift

Code, pre-registration and derived numbers for the manuscript

> **Actual versus nominal false-alarm rates of drift monitors for low-cost gas sensors, measured against a reference analyzer** (Dat Lam Quoc, FPT University; under review).

The study measures the actual probability that sequential drift monitors raise a false recalibration alarm on calibration residuals of low-cost gas sensors with a co-located reference analyzer (UCI Air Quality), separates the effect of the daily cycle from that of autocorrelation, and compares corrections (self-normalised blocks, an adaptive block plug-in monitor, a frozen long-run-variance CUSUM) under a pre-registered protocol.

## Contents
| Path | What |
|---|---|
| `prereg/PREREG_H35.md` | **Original** pre-registration (Vietnamese), tagged `prereg-H35-v1` before any confirmatory run |
| `prereg/AMENDMENT_LOG.md` | **Original** amendment log (Vietnamese), 8 amendments, each tagged before the runs it governs (`prereg-H35-v1a` … `v1e`) |
| `prereg/PREREG_H35.sha256` | SHA-256 of the two original files (current versions) |
| `prereg/*_en.md` | English translations, added at release; **not** covered by the hash. The Vietnamese files are the registered documents |
| `code/` | Simulation, analysis, table and figure scripts (Python) |
| `results/*.json` | Raw result files of each experiment (summary statistics with Wilson intervals; no dataset rows) |
| `results/numbers.json` | Every number used in the manuscript, built from the result files |
| `results/fig/` | Figures generated from `numbers.json` |

Commit dates in the history are the original dates of work: the pre-registration and each amendment precede the corresponding results. This public repository is a filtered copy of a private working repository (only `code/`, `prereg/`, `results/` were kept); commit hashes therefore differ from the private copy, while file contents and the SHA-256 of the pre-registration are unchanged.

### Checking the pre-registration hash
`.gitattributes` keeps `prereg/` byte-exact on every platform. Then:
```
cd prereg && sha256sum -c PREREG_H35.sha256
git show prereg-H35-v1:prereg/PREREG_H35.md | sha256sum      # same hash as registered on 21-11-2026
```

## Data (not included)
| Dataset | Source | Licence | Use |
|---|---|---|---|
| UCI Air Quality (#360) | https://archive.ics.uci.edu/dataset/360/air+quality | CC BY 4.0; "exclusively for research purposes" | pre-registered decisions (4 referenced channels) |
| UCI Gas Sensor Array Drift at Different Concentrations (#270) | https://archive.ics.uci.edu/dataset/270 | CC BY 4.0; "exclusively for research purposes" | exploratory, batch level only |
| Intel Berkeley Research Lab | http://db.csail.mit.edu/labdata/labdata.html | no stated licence — not redistributed | exploratory |

```
python code/get_data.py <dataset_dir>
export THS_DATASETS=<dataset_dir>          # Windows: set THS_DATASETS=<dataset_dir>
```
`get_data.py` was written at release time and has not been run end to end by the author (the original runs used local copies with the same layout).

## Environment
Python 3.11; `pip install -r requirements.txt` (numpy, scipy, pandas, matplotlib, python-dateutil; pymupdf only for `check_render.py`).

## Reproducing
### Numbers, tables and figures from the stored results (seconds, no data needed)
```
python -I code/build_numbers.py results results/numbers.json
python code/paper_figs.py results/numbers.json results/fig
python code/make_tex.py results/numbers.json <manuscript_dir>     # LaTeX macros and tables used by the manuscript
```

### Re-running the experiments (needs the data; seeds are fixed inside the scripts)
Runtimes were not recorded systematically; the Monte-Carlo experiments take from minutes to several hours each on a desktop CPU.

| Result file | Command | Manuscript |
|---|---|---|
| `aq_meta.json` | `python code/aq_meta.py results/aq_meta.json` | Table 1 |
| `c1_e1_aq.json`, `c1_e1_syn.json`, `c1_e1_intel.json` | `python code/c1_run.py e1_aq|e1_syn|e1_intel results/c1_<exp>.json` | Table 3, Figure 1 (RQ1); Table 5 (F1, F2); out-of-scope nulls |
| `pilot_B.json`, `s1b_a_decompose.json`, `q0_deseason_ctm.json` | `python code/pilot_B.py <AirQualityUCI.csv> results/pilot_B.json` (same pattern for the others) | Table 4 (RQ2) |
| `c1_e2_delay.json`, `c1_e2_arl0.json` | `python code/c1_run.py e2_delay|e2_arl0 …` | Table 6, Figure 2, H3 |
| `c1_e3.json`, `c1_e4.json`, `c1_e5.json` | `python code/c1_run.py e3|e4|e5 …` | sensitivity grid, assumption (E), H7 |
| `c1_e7.json` | `python code/c1_e78.py e7 results/c1_e7.json` | Tables 5–6, Figure 3 |
| `c1_e8.json` | `python code/c1_e78.py e8 results/c1_e8.json` | Table 7, Figure 4 |
| `s30_oos.json` | `python code/s30_oos.py results/s30_oos.json` | out-of-sample check, Figure 3 (open markers) |
| `e9_ref_noise.json` | `python code/e9_ref_noise.py results/e9_ref_noise.json` | Table 8 |
| `q2_adaptive.json`, `q2_selfcheck.json` | `python code/q2_adaptive.py|q2_selfcheck.py <AirQualityUCI.csv> …` | exploratory checks of Theorem 3 (before pre-registration) |
| `pilot_A.json`, `s1b_thm1_check.json`, `s1b_thm2_sim.json` | see usage line at the top of each script | exploratory (before pre-registration) |
| `s1b_prop0_cex.json` | produced inline during an exploratory step; no separate script | numbers.json key only, not used in the manuscript tables |

Seeds: 20261121–25 (main), 20261200–04 (threshold calibration, E7), 20261300–04 (out-of-sample, Amendment 7), 20261400+ (reference noise, E9). `results/c1_e3_INVALID_run1.json` is kept for traceability (Amendment 4).

`code/` also contains the literature and venue screening tools used while planning the study (`h35api.py`, `novelty_scan.py`, `make_refs_bib.py`, `verify_refs.py`, `venue_g2.py`, `split_gfa.py`, `build_novelty_csv.py`, `dataset_check.py`); they query OpenAlex, Crossref and DataCite, read an optional API key from a local `.env` that is not part of this repository, and are not needed to reproduce the results.

## Use of generative AI
The code was written and run with Claude Opus 5.5 (Anthropic) via Claude Code, under the author's direction; see the manuscript's Methods and Acknowledgements.

## Licence and citation
Code: MIT (`LICENSE`). The datasets keep their own licences. Please cite the manuscript; `CITATION.cff` describes this repository.
