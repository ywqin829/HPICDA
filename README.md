# HPICDA — Human Pulmonary Infection Conditional Dependency Atlas

Analysis code for the HPICDA resource manuscript: a baseline-fitness-aware
re-analysis of 39 CRISPR loss-of-function screens of human respiratory
infection, and an independent holdout adjudication of the candidates it
produces.

The manuscript's central results are:

1. The composite (HG3) score does not validate: at chance as a predictor of
   independent holdout validation on the frozen pre-registered layer
   (gene-level AUROC 0.524) and below chance on the SC2-03-inclusive layer
   (0.447), while a plain count of the pathogens in which a gene recurred
   reaches 0.686 on the frozen layer.
2. Cellular context, not intervention direction, is the pair-level determinant
   of replication (context-matched pair, OR 2.58, 95% CI 1.11–6.02).
3. An estimand-aware reconciliation of the screens into 852,216 gene × screen
   calls, separated into conditional dependency, baseline-fitness confounding
   and signal-state alteration, with a prospectively frozen set of 353
   candidates.

This repository contains **code only**. The data tables and the figures are
distributed separately (see below).

## Repository contents

```
scripts/       analysis pipeline, 01–60 (Python, plus 2 R QC/overlay scripts)
R/             publication figure scripts (7 figures + shared theme)
runner/        checkpoint / pipeline driver (run_pipeline.py, PIPELINE_INDEX.yaml)
environment_lock.txt   pinned environment
```

## Data

The manuscript's supplementary tables (S1–S20b) are archived on Zenodo:

> **Zenodo data archive:** DOI [`10.5281/zenodo.23357807`](https://doi.org/10.5281/zenodo.23357807)

The archive contains `tables/` — 34 TSV files (S1–S20b) plus `TABLES_INDEX.md`
and `TABLES_INDEX.json`.

The pipeline's intermediate data are **not** bundled, here or on Zenodo: the
CRISPR screens (GEO), the DepMap 26Q1 release and the pharmacological resources
(ChEMBL / Open Targets / Orphanet) are external inputs, re-fetched by the
pipeline. See the manuscript Methods and `environment_lock.txt` for sources and
versions.

Every script resolves its working root from the `HPICDA_ROOT` environment
variable (default: the current directory), so point it at wherever you keep the
working tree:

```bash
export HPICDA_ROOT=/path/to/HPICDA        # bash
setx HPICDA_ROOT "D:\path\to\HPICDA"      # cmd (new shells)
```

## Requirements

- Python 3 (see `environment_lock.txt` for the pinned package versions)
- R ≥ 4.4 for `R/` (ggplot2, patchwork, data.table, ggrepel, svglite, ragg) and
  for `scripts/31_*.R`, `38_*.R`, `39_*.R`

## Running

```bash
python runner/run_pipeline.py --list      # list the 48 indexed steps and checkpoint state
python runner/run_pipeline.py --verify    # re-hash every checkpoint output
```

Reproducibility is staged, not single-shot: each stage is committed behind a
SHA-256 checkpoint and the full pipeline was never run end-to-end unattended.
That is stated in the manuscript as a limitation, not a defect.

To regenerate the publication figures:

```bash
cd <repository root>
Rscript R/fig1_study_architecture.R
# … fig2_landscape.R … fig7_holdout.R
```

Figure scripts write PDFs to `25_manuscript_package/figures/`. Run them from the
repository root (the shared theme is sourced as `R/hpicda_theme.R`).

## Citation

If you use this code, please cite the manuscript (details to follow) and the
Zenodo code/data archives. See `CITATION.cff`.

## License

MIT (see `LICENSE`). The Zenodo data archive is released under CC BY 4.0.
