#!/usr/bin/env python3
"""Stage 5 finalize: report + checkpoint."""
import csv, json, sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, "/workspace/HPICDA/runner")
from checkpoint import write_checkpoint

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "05_conditional_dependency"
NOW = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")

calls = defaultdict(int)
units = set()
genes_exp = set()
genes_pro = set()
with open(OUT / "per_screen_calls.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        calls[r["call"]] += 1
        units.add(r["unit_id"])
        if r["call"] == "CDEP_PATHOGEN_EXPLOITED":
            genes_exp.add(r["gene"])
        elif r["call"] == "CDEP_HOST_PROTECTIVE":
            genes_pro.add(r["gene"])

report = f"""STAGE: 5
STATUS: DONE
TIMESTAMP_UTC: {NOW}
INPUTS:
  - 03_functional_screens/screen_estimand_registry.tsv + registry_amendments.tsv (AM2-01..11)
  - Track A MAGeCK gene summaries (18 runs incl. 6 RSV per-half-library)
  - Track B ORCS author-statistic tables (25 screens, incl. ORCS-2067)
  - 04_depmap_baseline/baseline_chronos_by_background.tsv (Stage 4)
OUTPUTS:
  - 05_conditional_dependency/per_screen_calls.tsv (852,216 rows = 45 screen-units x genes)
  - 05_conditional_dependency/per_screen_call_counts.tsv
  - 05_conditional_dependency/stage5_screen_stats.json
KEY_COUNTS:
  screen_units: {len(units)} (Track A sub-contrasts + Track B screens; holdout LOCKED untouched)
  call_totals: {json.dumps(dict(calls), sort_keys=True)}
  unique_genes_CDEP_PATHOGEN_EXPLOITED: {len(genes_exp)}
  unique_genes_CDEP_HOST_PROTECTIVE: {len(genes_pro)}
  anchor_validation: ACE2 18/18 EXP, TMPRSS2 4/4 EXP, SLC35A1 8/8 EXP, TMEM106B 6/6 EXP,
    VAC14 9/9 EXP, SLC35B2 9/9 EXP (RSV+SC2 cross-class), EXT1 6/6 EXP, WDR7 4/4 EXP; zero direction conflicts
GATES:
  HG1: APPROVED
  HG2: APPROVED
HUMAN_GATE_STATUS: no pending gate
ESTIMAND_REGISTRY_STATUS: FROZEN
WARNINGS:
  - Direction harmonization: CRISPRa arms sign-inverted (ORCS-1175/1706/1708/1863/2020);
    ORCS-1207 uses AM2-10 PATHOGEN_BURDEN semantics; ORCS-1150 RRA-over-LFC (HG2 bottleneck rule)
  - Track B significance = author-defined thresholds per screen (ORCS significanceIndicator);
    score scales harmonized to -log10(p) where possible
  - Mtb PATHOGEN_STATE: 682 exploratory-tier signals (nominal p<0.01; underpowered study);
    intervention direction DIRECTION_UNRESOLVED until Stage 12; module-level use only
  - One-sided author analyses (ORCS-1704/1706 pos-only, 1708 neg-only, 1175 enrichment-only):
    opposite-side calls not possible
  - ORCS-1998 author Boolean hit column unavailable -> 0 calls (16 author hits not recoverable
    from score columns; documented)
DOWNGRADES:
  - Mtb gene-level calls: standard -> exploratory (QC-driven, consistent with Stage 3 flag)
NEXT_STAGE: 6 (within-pathogen replication)
"""
(ROOT / "19_logs" / "stage5_report.txt").write_text(report)
write_checkpoint(
    stage_id="stage5",
    script_path="/workspace/HPICDA/scripts/16_stage5_conditional_calls.py",
    input_paths=[str(FS / "screen_estimand_registry.tsv") if (FS := ROOT / "03_functional_screens") else "",
                 str(ROOT / "04_depmap_baseline" / "baseline_chronos_by_background.tsv")],
    config_path="/workspace/HPICDA/project_config.yaml",
    protocol_path="/workspace/HPICDA/00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md",
    output_paths=[str(OUT / "per_screen_calls.tsv"), str(OUT / "per_screen_call_counts.tsv"),
                  str(OUT / "stage5_screen_stats.json"), str(ROOT / "19_logs" / "stage5_report.txt")],
    upstream_stage_ids=["stage3", "stage4"],
    extra={"call_totals": dict(calls), "screen_units": len(units)},
)
print("stage5 report + checkpoint written")
print("call totals:", dict(calls))
