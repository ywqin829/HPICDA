#!/usr/bin/env python3
"""HPICDA Stage 3: write per-study QC_REPORT.md + summary table + checkpoint."""
import csv, json, sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, "/workspace/HPICDA/runner")
from checkpoint import write_checkpoint

ROOT = Path("/workspace/HPICDA/03_functional_screens")
NOW = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")

qc = json.load(open(ROOT / "stage3_qc_raw.json"))
qc_by_id = {r["study_id"]: r for r in qc}

# Per-study QC verdicts (12 checks assessed from metrics + Stage 2 validation)
verdicts = {
    "GSE210689": {
        "status": "PASS_WITH_WARNINGS",
        "rationale": "Coverage 97.5% but mean depth ~1.2 reads/guide; control-replicate Spearman 0.13 "
                     "(severe bottleneck/drift); >50% zero-count guides in UI/Lat; no FDR<0.25 hit. "
                     "GPR156 directionally recovered (pos p=0.0038, rank 63/18,059), consistent with "
                     "published MMGT1-GPR156 axis -> no conflict with publication. UNDERPOWERED: "
                     "contributes to synthesis only via directional/vote-counting evidence, never as "
                     "a powered negative.",
        "flags": ["LOW_DEPTH", "LOW_REPLICATE_CORRELATION", "UNDERPOWERED"],
    },
    "GSE178699": {
        "status": "PASS_HIGH_CONFIDENCE",
        "rationale": "Coverage 99.6%; control-replicate Spearman 0.75; ACE2 pos_rank#1 in all 3 cell "
                     "lines (FDR 0.0008-0.005), TMPRSS2 #2 in Calu-3; KMT2C/KDM6A nominal. "
                     "Author method qGI; MAGeCK reanalysis concordant on anchor genes.",
        "flags": [],
    },
    "GSE111166": {
        "status": "PASS_HIGH_CONFIDENCE",
        "rationale": "Coverage 99.6%; SLC35A1 pos_rank#1 in all 3 contrasts (p=2.3e-7); progressive "
                     "enrichment across rounds verified at guide level (sequence-keyed merge). "
                     "ORCS-1150 is a bottlenecked single-replicate survival screen: RRA scores used "
                     "as primary evidence per HG2 rule (median LFC sign unreliable).",
        "flags": ["BOTTLENECK_SCREEN_ORCS-1150"],
    },
    "Li2020": {
        "status": "PASS_WITH_WARNINGS",
        "rationale": "Coverage 99.8%; WDR7 #1 (FDR 0.0003), SLC35A1/ATP6V1A/ATP6AP1/CMTR1/CCDC115/"
                     "TMEM199 all recovered at author-concordant ranks. WARNING: source data are "
                     "log2RPM, reverse-transformed to pseudo-counts at assumed 10M depth "
                     "(documented approximation); single-replicate bins (P5_6 vs P4, P6 vs P4).",
        "flags": ["PSEUDO_COUNTS", "SINGLE_REPLICATE_BINS"],
    },
    "GSE175666": {
        "status": "PASS_HIGH_CONFIDENCE",
        "rationale": "Coverage 99.9%; control-replicate Spearman 0.92; AP1G1 #1 + ACE2 #2 "
                     "(FDR 0.0007 each) exactly match author zscore sheet; CRISPRa arm ACE2 "
                     "neg_rank#1 (sign-inverted at Stage 5); Caco-2-ACE2 arm ACE2 #1.",
        "flags": [],
    },
    "PRJNA1281304": {
        "status": "PASS_HIGH_CONFIDENCE",
        "rationale": "DESIGN DISCOVERY: screens A and B used GeCKOv2 half-libraries A and B "
                     "respectively (99%/97% of reads on the cognate half-library) -> independent "
                     "half-library replicates, NOT guide-level replicates. Combined-library MAGeCK "
                     "run superseded by per-half-library runs (6 contrasts). Heparan-sulfate pathway "
                     "replicates across half-libraries: SLC35B2 viability A rank3/B rank4 "
                     "(FDR<0.002 both), GFPNEG A rank3/B rank2; EXT2 viability A rank1; EXT1 "
                     "viability B rank5; B3GALT6 viability B rank8. Coverage 99.8%; ~31% "
                     "empty-vector reads excluded (expected for this amplicon design).",
        "flags": ["HALF_LIBRARY_DESIGN"],
    },
    "GSE154751": {
        "status": "LIST_ONLY",
        "rationale": "Context-only screen (Pseudomonas ExoU toxin sensitivity, not a discovery "
                     "pathogen). Coverage 99.1%; single library-bank control (no true replicates); "
                     "DNAJC5 top hit (FDR 0.035). Contributes context only, never core synthesis.",
        "flags": ["CONTEXT_ONLY", "NO_REPLICATES"],
    },
}

# Write per-study QC_REPORT.md
for sid, v in verdicts.items():
    r = qc_by_id.get(sid, {})
    lines = [
        f"# QC REPORT — {sid}",
        "",
        f"**Status: {v['status']}**",
        "",
        f"date: {NOW}",
        "",
        "## Metrics",
        "",
        f"- guides: {r.get('n_guides', 'NA')} | genes: {r.get('n_genes', 'NA')} | "
        f"median guides/gene: {r.get('median_guides_per_gene', 'NA')}",
        f"- library coverage (>=1 read in >=1 sample): {r.get('coverage', 'NA')}",
        f"- library sizes: {json.dumps(r.get('lib_sizes', {}))}",
        f"- zero-count fraction: {json.dumps(r.get('zero_frac', {}))}",
        f"- control replicate Spearman: {r.get('control_rep_spearman', 'NA')}",
        f"- treated-arm Gini: {json.dumps(r.get('gini_treated', {}))}",
        f"- known controls: {json.dumps(r.get('known_controls', {}))}",
        "",
        "## 12-check assessment",
        "",
        "1. library coverage: see metrics above",
        "2. replicate correlation: see control Spearman (RSV: half-library design, gene-level replication)",
        "3. non-targeting guides: present in libraries; used by MAGeCK null where applicable",
        "4. known controls: see known_controls recovery above",
        "5. guide concordance: multi-guide agreement via MAGeCK RRA (validated on anchor genes)",
        "6. selection bottleneck: see Gini + zero fraction; ORCS-1150 flagged bottlenecked",
        "7. screen saturation: genome-wide libraries, median guides/gene as above",
        "8. baseline abundance: control-arm library sizes as above",
        "9. shared controls: none across studies (HG2-audited); study-internal controls only",
        "10. cell identity: per registry (Cellosaurus-mapped at Stage 0)",
        "11. endpoint direction: per HG1-approved adjudication (+ AM2-10 for ORCS-1207)",
        f"12. author vs reanalysis consistency: {r.get('author_note', 'see Stage 2 report')}",
        "",
        "## Rationale",
        "",
        v["rationale"],
        "",
        f"flags: {', '.join(v['flags']) if v['flags'] else 'none'}",
        "",
    ]
    d = ROOT / sid
    d.mkdir(exist_ok=True)
    (d / "QC_REPORT.md").write_text("\n".join(lines))
    print(f"wrote {sid}/QC_REPORT.md [{v['status']}]")

# Track B group QC (author stats only -> LIST_ONLY or PASS by author evidence)
tb_lines = [
    "# QC REPORT — Track B ORCS author-statistic screens (24)",
    "",
    "**Status: PASS_WITH_WARNINGS (group)**",
    "",
    f"date: {NOW}",
    "",
    "Track B screens contribute author gene-level statistics, not reanalysis. QC limited to:",
    "- table completeness: 23/24 genome-wide (18-21k genes); ORCS-1175 curated subset (1,190)",
    "- known-anchor recovery: TMEM106B/VAC14 top in ORCS-1220/1385; ACE2 RRA rank 8 in ORCS-1431;",
    "  ACE2 pos 1.8e-28 in ORCS-1704; PIK3C3 FDR 0.0008 in ORCS-1821; CTSL pos_rank 2 (MOESM7)",
    "- score heterogeneity (MAGeCK score / Z-score / RRA rank / Log2FC / FDR) handled by per-screen",
    "  rank normalization at Stage 5",
    "",
    "WARNING: author pipelines differ; effect directions harmonized via HG1 adjudication.",
    "GSE177545 (ORCS-1821..1827): no raw data exists; Track B is the ceiling for these screens.",
    "",
]
(ROOT / "TrackB_ORCS").mkdir(exist_ok=True)
(ROOT / "TrackB_ORCS" / "QC_REPORT.md").write_text("\n".join(tb_lines))
print("wrote TrackB_ORCS/QC_REPORT.md [PASS_WITH_WARNINGS]")

# Summary TSV
with open(ROOT / "stage3_qc_summary.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["study_id", "status", "n_guides", "n_genes", "coverage",
                "control_rep_spearman", "flags"])
    for sid, v in verdicts.items():
        r = qc_by_id.get(sid, {})
        w.writerow([sid, v["status"], r.get("n_guides", ""), r.get("n_genes", ""),
                    r.get("coverage", ""), r.get("control_rep_spearman", ""),
                    ";".join(v["flags"])])
    w.writerow(["TrackB_ORCS_24screens", "PASS_WITH_WARNINGS", "", "", "", "", "AUTHOR_STATS_ONLY"])
print("wrote stage3_qc_summary.tsv")

# Stage 3 report
n_pass_hi = sum(1 for v in verdicts.values() if v["status"] == "PASS_HIGH_CONFIDENCE")
n_pass_warn = sum(1 for v in verdicts.values() if v["status"] == "PASS_WITH_WARNINGS")
report = f"""STAGE: 3
STATUS: DONE
TIMESTAMP_UTC: {NOW}
INPUTS:
  - 03_functional_screens/raw/*/counts*.txt (Track A count matrices)
  - 03_functional_screens/raw/*/*.gene_summary.txt (18 Track A MAGeCK runs)
  - 03_functional_screens/trackB_orcs/*.gene_scores.tsv (24 Track B tables)
OUTPUTS:
  - 03_functional_screens/<study>/QC_REPORT.md x 7 Track A studies + TrackB_ORCS group report
  - 03_functional_screens/stage3_qc_summary.tsv
  - 03_functional_screens/stage3_qc_raw.json
  - RSV per-half-library MAGeCK reruns: raw/PRJNA1281304/RSV_lib[AB]_*.gene_summary.txt (6)
KEY_COUNTS:
  trackA_studies_qc: 7
  PASS_HIGH_CONFIDENCE: {n_pass_hi} (GSE178699, GSE111166, GSE175666, PRJNA1281304)
  PASS_WITH_WARNINGS: {n_pass_warn} (GSE210689 low-depth, Li2020 pseudo-counts)
  LIST_ONLY: 1 (GSE154751 context-only)
  FAIL: 0
  ENDPOINT_UNRESOLVED: 0
GATES:
  HG1: APPROVED
  HG2: APPROVED
HUMAN_GATE_STATUS: no pending gate
ESTIMAND_REGISTRY_STATUS: FROZEN (amendments file only)
WARNINGS:
  - RSV half-library design discovered (screens A/B = GeCKOv2 lib A/B); combined-library run
    superseded by 6 per-half-library contrasts; registry amendment AM2-11 required at Stage 5
    accounting (each RSV condition = 2 independent half-library sub-screens)
  - GSE210689 underpowered: directional/vote-counting evidence only
  - ORCS-1150 bottlenecked: RRA-over-LFC rule (HG2)
DOWNGRADES:
  - GSE210689: powered-evidence -> directional-evidence (QC-driven)
NEXT_STAGE: 4 (DepMap 26Q1 baseline freeze)
"""
(ROOT.parent / "19_logs" / "stage3_report.txt").write_text(report)
print("wrote stage3_report.txt")

# Checkpoint
outputs = [ROOT / "stage3_qc_summary.tsv", ROOT / "stage3_qc_raw.json",
           ROOT.parent / "19_logs" / "stage3_report.txt"]
outputs += sorted(ROOT.glob("*/QC_REPORT.md"))
outputs += sorted(ROOT.glob("raw/PRJNA1281304/RSV_lib*.gene_summary.txt"))
outputs += sorted(ROOT.glob("raw/PRJNA1281304/counts_lib*.txt"))
write_checkpoint(
    stage_id="stage3",
    script_path="/workspace/HPICDA/scripts/14_stage3_reports.py",
    input_paths=[str(ROOT / "screen_estimand_registry.tsv"), str(ROOT / "stage3_qc_raw.json")],
    config_path="/workspace/HPICDA/project_config.yaml",
    protocol_path="/workspace/HPICDA/00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md",
    output_paths=[str(p) for p in outputs],
    upstream_stage_ids=["stage2"],
    extra={"verdicts": {k: v["status"] for k, v in verdicts.items()}},
)
print("checkpoint written: stage3")
