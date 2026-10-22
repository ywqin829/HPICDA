# HPICDA — Human Pulmonary Infection Conditional Dependency Atlas

Analysis code for the HPICDA resource manuscript: a baseline-fitness-aware
re-analysis of 39 CRISPR loss-of-function screens of human respiratory
infection, and an independent holdout adjudication of the candidates it
produces.

The manuscript's central results are:

1. Same-pathogen replication is low (55 of 654 evaluable gene × holdout pairs,
   8.4%), and the frozen HG3 composite is at chance as a global predictor of
   independent validation on the pre-registered layer (gene-level AUROC 0.524)
   and below chance on the SC2-03-inclusive layer (0.447).
2. An unweighted pathogen count initially reached AUROC 0.686, but this
   advantage was not robust after 5,000 equal-opportunity resamples. The number
   or frequency of independent discovery units hit retained a positive paired
   advantage over HG3, whereas pathogen breadth did not. Cellular context is
   exploratory because it cannot be separated from holdout identity.
3. An estimand-aware reconciliation of the screens into 852,216 gene × screen
   calls, separated into conditional dependency, baseline-fitness confounding
   and signal-state alteration, with a prospectively frozen set of 353
   candidates.

This repository contains **code only**. The data tables and the figures are
distributed separately (see below).

## Repository contents

```
scripts/       analysis pipeline, 01–64 plus the revision audit script 66
docs/          post-hoc module-specificity audit and decision record
R/             publication figure scripts (7 figures + shared theme)
runner/        checkpoint / pipeline driver (run_pipeline.py, PIPELINE_INDEX.yaml)
environment_lock.txt   pinned environment
```

## Data

The manuscript's supplementary tables and complete revision analyses are
archived on Zenodo:

> **Zenodo data archive:** DOI [`10.5281/zenodo.23550165`](https://doi.org/10.5281/zenodo.23550165)

The archive contains `tables/` — the complete source tables S1–S21b, including
the full `S5_per_screen_calls.tsv.gz`, plus `TABLES_INDEX.md` and
`TABLES_INDEX.json`. `revision_analyses/` contains the opportunity-adjusted and
equal-opportunity benchmark outputs.

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

The revision-stage analysis and packaging scripts are
`scripts/61_opportunity_balanced_benchmark.py` through
`scripts/64_build_combined_submission_pdf.py`. They use the same
`HPICDA_ROOT` convention as the rest of the pipeline.

`scripts/66_reaudit_module_holdout_specificity.py` records why the proposed
module-by-holdout specificity atlas was not added to the manuscript. The
original proposal counted non-evaluable cross-pathogen rows in its denominator;
the corrected within-holdout analysis does not survive multiplicity
correction.

Figure scripts write PDFs to `25_manuscript_package/figures/`. Run them from the
repository root (the shared theme is sourced as `R/hpicda_theme.R`).

## Citation

If you use this code, please cite the manuscript (details to follow) and the
Zenodo code/data archives. See `CITATION.cff`.

## License

MIT (see `LICENSE`). The Zenodo data archive is released under CC BY 4.0.
