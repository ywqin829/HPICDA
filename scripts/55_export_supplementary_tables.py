#!/usr/bin/env python
"""Export Supplementary Tables S1-S19 for the HPICDA resource manuscript.

Each table is written as a UTF-8 (no BOM) tab-separated file with a stable
column order.  S5 (the full gene x screen call table) is additionally written
as .tsv.gz because of its size.  A machine-readable index and a hash manifest
(S18) are produced alongside.

Table bands: S1-S6 inputs and calls, S7-S7c cross-pathogen evidence,
S8-S11c context / genetics / modules, S12-S17b frozen candidates and
translational annotation, S19-S19f the holdout benchmark release (S19f being the
two-layer version of the ranking benchmark), S20-S20b the context-versus-composite
comparison, and S18 last as the reproducibility manifest over the headline
artifacts.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import shutil
from pathlib import Path

ROOT = Path(
    os.environ.get("HPICDA_ROOT", Path(__file__).resolve().parents[1])
).resolve()
OUT = ROOT / "25_manuscript_package" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

# (id, output_name, caption, [source files])
TABLES = [
    ("S1", "S1_screen_audit",
     "Registry and estimand audit of every functional infection screen considered for the atlas.",
     ["03_functional_screens/screen_estimand_registry.tsv"]),
    ("S2", "S2_endpoint_ontology",
     "Controlled vocabulary for endpoints, directions and pooling rules.",
     ["02_ontology/hpicda_ontology_v1.tsv"]),
    ("S3", "S3_raw_qc",
     "Raw-data and reanalysis quality-control summary per study.",
     ["03_functional_screens/stage3_qc_summary.tsv"]),
    ("S4", "S4_depmap_baseline",
     "DepMap 26Q1 Chronos baseline fitness per gene in the pulmonary cell-line panel.",
     ["04_depmap_baseline/baseline_chronos_by_background.tsv"]),
    ("S5", "S5_per_screen_calls",
     "Complete gene x screen conditional-dependency call table (852k rows; delivered gzipped).",
     ["05_conditional_dependency/per_screen_calls_v2.tsv"]),
    ("S5b", "S5b_per_screen_call_counts",
     "Call-class counts per screen (companion summary to S5).",
     ["05_conditional_dependency/per_screen_call_counts_v2.tsv"]),
    ("S6", "S6_within_pathogen_replication",
     "Within-pathogen replication tiers across independence clusters.",
     ["06_within_pathogen/within_pathogen_replication_v2.tsv"]),
    ("S7", "S7_cross_pathogen_recurrence",
     "Cross-pathogen meta-analysis per gene (Fisher combined evidence, with leave-one-independence-cluster-out stability support).",
     ["07_cross_pathogen/convergence_meta_v2_gene.tsv"]),
    ("S7b", "S7b_loso_stability",
     "Leave-one-independence-cluster-out stability of the cross-pathogen convergence estimates (14 clusters).",
     ["07_cross_pathogen/convergence_v2_loso.tsv"]),
    ("S7c", "S7c_track_ab_sensitivity",
     "Track A / Track B comparability sensitivity analyses per gene: one-vote-per-independent-unit (rank-based) and direction-only scoring, with per-gene convergence status under each.",
     ["24_closing_analyses/track_ab_sensitivity_gene.tsv"]),
    ("S8", "S8_pulmonary_datasets",
     "Human pulmonary response datasets used for transcriptomic context and directional concordance.",
     ["09_transcriptomic_context/stage8_dataset_selection.tsv"]),
    ("S8b", "S8b_per_gene_concordance",
     "Per-candidate transcriptomic concordance summary across the pulmonary datasets.",
     ["09_transcriptomic_context/stage8_per_gene_summary.tsv"]),
    ("S9", "S9_singlecell_expression",
     "Single-cell expression summary of frozen candidates across the human lung reference.",
     ["10_singlecell_context/stage9_candidate_lung_expression_summary.tsv"]),
    ("S10", "S10_human_genetics",
     "Human loss-of-function constraint (gnomAD LOEUF) and genetics evidence tier per candidate.",
     ["08_human_genetics/stage10_genetics_evidence.tsv"]),
    ("S11", "S11_module_coessentiality",
     "Module-level DepMap co-essentiality with permutation background statistics.",
     ["11_modules/stage11_module_coessentiality.tsv"]),
    ("S11b", "S11b_intact_ppi_support",
     "IntAct physical interaction support within curated modules.",
     ["11_modules/stage11_intact_support.tsv"]),
    ("S11c", "S11c_contradictory_pairs",
     "Stratified opposing-direction calls and their baseline-essentiality confounding status.",
     ["11_modules/stage11_contradictory_pairs_stratified.tsv"]),
    ("S12", "S12_frozen_candidates",
     "The 353 pre-registered candidates with all frozen annotation fields.",
     ["09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2_annotated.tsv"]),
    ("S13", "S13_holdout_verdicts",
     "Per gene x holdout adjudication verdicts.",
     ["20_stage13_holdout_adjudication/stage13_gene_holdout_verdicts.tsv"]),
    ("S13b", "S13b_holdout_manifest",
     "Pre-registered holdout screens and their independence rationale.",
     ["09_candidate_freeze/holdout_manifest_v2.tsv"]),
    ("S14", "S14_pharmacology",
     "Open Targets-derived tractability annotation per candidate.",
     ["21_stage14_druggability/stage14_druggability_axis.tsv"]),
    ("S15", "S15_launched_curation",
     "Machine first-pass curation of launched / clinical-stage compounds mapped to candidates; classifications remain provisional pending source-level human review.",
     ["21_stage14_druggability/stage14_launched_curation.tsv"]),
    ("S16", "S16_safety_tolerability",
     "Human tolerability axis (LOEUF, pLI, Genebass, HPA) per candidate.",
     ["16_revision_analyses/16_revision_r12_tolerability/tolerability_axis.tsv"]),
    ("S17", "S17_contradiction_matrix",
     "Per-pathogen direction contradiction matrix for frozen candidates.",
     ["09_candidate_freeze/candidate_contradiction_matrix_v2.tsv"]),
    ("S17b", "S17b_pareto_priority",
     "Stage 15 multi-axis Pareto prioritisation of frozen candidates.",
     ["22_stage15_pareto/stage15_pareto_priority.tsv"]),
    ("S19", "S19_benchmark_rankings",
     "Honest benchmark at the gene level: HG3 composite score versus four naive consensus features and a best-rank percentile, over the 333 evaluable candidates (AUROC, AUPRC, top-K validated rates, odds ratios).",
     ["24_closing_analyses/benchmark_naive_vs_hpicda.tsv"]),
    ("S19b", "S19b_benchmark_per_gene",
     "Per-gene benchmark input: validation status, evaluability and the five ranking features for all 353 frozen candidates.",
     ["24_closing_analyses/benchmark_gene_level_data.tsv"]),
    ("S19c", "S19c_topk_calibration",
     "Top-K holdout-PASS enrichment at the pair and gene levels, computed separately for the frozen adjudication and for the SC2-03-inclusive analysis set.",
     ["24_closing_analyses/holdout_topk_calibration.tsv"]),
    ("S19d", "S19d_context_determinants",
     "Stratified holdout PASS rate by cellular-context match and intervention direction (source of Fig. 7d).",
     ["24_closing_analyses/context_determinants_stratified.tsv"]),
    ("S19e", "S19e_archetype_case_studies",
     "Case-study evidence for the three archetypes (ACE2 anchor, heparan-sulfate module, V-ATPase module) across discovery, holdout and baseline-fitness axes.",
     ["24_closing_analyses/case_studies_ace2_hs_vatpase.tsv"]),
    ("S19f", "S19f_benchmark_two_layer",
     "Two-layer honest benchmark at the gene level: AUROC (with bootstrap 95% CI) of the HG3 composite score and the five naive consensus features over the 333 evaluable candidates, computed separately on the frozen pre-registered outcome (51/333 validated) and on the SC2-03-inclusive outcome (132/333). Source of the two-layer Fig. 7c; the SC2-03-inclusive column is carried alongside the published values of S19 as a cross-check.",
     ["24_closing_analyses/benchmark_two_layer.tsv"]),
    ("S20", "S20_context_vs_composite_gene",
     "Head-to-head comparison of the HG3 composite score, cellular-context exposure, intervention direction and a plain pathogen count as predictors of independent validation at the gene level, over the 333 evaluable candidates, computed separately on the frozen pre-registered outcome and on the SC2-03-inclusive outcome (AUROC with bootstrap CI, sign-free AUROC, and paired differences against the composite).",
     ["24_closing_analyses/headtohead_frozen_gene.tsv"]),
    ("S20b", "S20b_context_vs_composite_pair",
     "Pair-level logistic models of holdout PASS on cellular-context match, intervention direction and Bio-Tier 2 status, marginally and with holdout fixed effects, in the frozen pre-registered layer and in the SC2-03-inclusive layer.",
     ["24_closing_analyses/headtohead_frozen_pair.tsv"]),
    ("S18", "S18_cleanroom_hashes",
     "SHA-256 reproducibility manifest of the frozen inputs and headline outputs.",
     []),  # generated below
]

HEADLINE_HASH_FILES = [
    "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv",
    "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2_annotated.tsv",
    "05_conditional_dependency/per_screen_calls_v2.tsv",
    "07_cross_pathogen/convergence_meta_v2_gene.tsv",
    "07_cross_pathogen/convergence_meta_v2_module.tsv",
    "08_human_genetics/stage10_genetics_evidence.tsv",
    "11_modules/stage11_module_coessentiality.tsv",
    "20_stage13_holdout_adjudication/stage13_gene_holdout_verdicts.tsv",
    "21_stage14_druggability/stage14_druggability_axis.tsv",
    "22_stage15_pareto/stage15_pareto_priority.tsv",
    "24_closing_analyses/benchmark_naive_vs_hpicda.tsv",
    "24_closing_analyses/benchmark_gene_level_data.tsv",
    "24_closing_analyses/holdout_topk_calibration.tsv",
    "24_closing_analyses/benchmark_two_layer.tsv",
    "24_closing_analyses/context_determinants_stratified.tsv",
    "24_closing_analyses/track_ab_sensitivity_gene.tsv",
    "24_closing_analyses/case_studies_ace2_hs_vatpase.tsv",
    "24_closing_analyses/headtohead_frozen_gene.tsv",
    "24_closing_analyses/headtohead_frozen_pair.tsv",
    "24_closing_analyses/context_vs_composite_gene.tsv",
    "07_cross_pathogen/convergence_v2_loso.tsv",
    "09_candidate_freeze/stage12_freeze_v2_runinfo.json",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_tsv(src: Path, dst: Path) -> tuple[int, int]:
    """Copy a TSV verbatim (normalising to LF, UTF-8 without BOM)."""
    rows = 0
    ncol = 0
    with open(src, "r", encoding="utf-8", newline="") as fin, \
         open(dst, "w", encoding="utf-8", newline="\n") as fout:
        rd = csv.reader(fin, delimiter="\t")
        wr = csv.writer(fout, delimiter="\t", lineterminator="\n")
        for i, row in enumerate(rd):
            if i == 0:
                ncol = len(row)
            wr.writerow(row)
            rows += 1
    return rows - 1, ncol


index = []
for tid, name, caption, sources in TABLES:
    if tid == "S18":
        continue
    src = ROOT / sources[0]
    if not src.exists():
        print(f"  !! missing {src}")
        continue
    if tid == "S5":  # large: gzip only
        dst = OUT / f"{name}.tsv.gz"
        with open(src, "rb") as fin, gzip.open(dst, "wb", compresslevel=6) as fout:
            shutil.copyfileobj(fin, fout, 1 << 20)
        with open(src, "r", encoding="utf-8", newline="") as fin:
            rows = sum(1 for _ in fin) - 1
            fin.seek(0)
            ncol = len(next(csv.reader(fin, delimiter="\t")))
        size_mb = dst.stat().st_size / 1e6
        print(f"  {tid:4s} {name:34s} {rows:>8,d} x {ncol:>2d}  ({size_mb:.1f} MB gz)")
        index.append(dict(id=tid, file=dst.name, caption=caption, rows=rows,
                          cols=ncol, gzip=True))
        continue
    dst = OUT / f"{name}.tsv"
    rows, ncol = copy_tsv(src, dst)
    print(f"  {tid:4s} {name:34s} {rows:>8,d} x {ncol:>2d}")
    index.append(dict(id=tid, file=dst.name, caption=caption, rows=rows, cols=ncol))

# --- S18: clean-room reproducibility manifest -------------------------------
man = []
for rel in HEADLINE_HASH_FILES:
    p = ROOT / rel
    if not p.exists():
        continue
    man.append(dict(artifact=rel, sha256=sha256(p), bytes=p.stat().st_size))
s18 = OUT / "S18_cleanroom_hashes.tsv"
with open(s18, "w", encoding="utf-8", newline="\n") as fh:
    wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
    wr.writerow(["artifact", "sha256", "bytes"])
    for r in man:
        wr.writerow([r["artifact"], r["sha256"], r["bytes"]])
captions = {t[0]: t[2] for t in TABLES}
index.append(dict(id="S18", file=s18.name, caption=captions["S18"],
                  rows=len(man), cols=3))
print(f"  S18  {'S18_cleanroom_hashes':34s} {len(man):>8,d} x  3")

with open(OUT / "TABLES_INDEX.json", "w", encoding="utf-8") as fh:
    json.dump(index, fh, indent=2, ensure_ascii=False)

# --- write a readable index -------------------------------------------------
lines = ["# Supplementary tables S1-S20b", "",
         "All tables are UTF-8, tab-separated, LF line endings.",
         "S5 is delivered gzipped because it contains 852k rows.", ""]
for r in index:
    lines += [f"## {r['id']} - `{r['file']}`",
              "", r["caption"], "",
              f"- Rows (excluding header): {r['rows']:,}",
              f"- Columns: {r['cols']}", ""]
(OUT / "TABLES_INDEX.md").write_text("\n".join(lines), encoding="utf-8")

print(f"\n  wrote {len(index)} tables -> {OUT}")
