#!/usr/bin/env python3
"""Finalize Stage 2: registry amendments, HG2 approval, stage report, checkpoint."""
import csv, hashlib, json, sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, "/workspace/HPICDA/runner")
from checkpoint import write_checkpoint

ROOT = Path("/workspace/HPICDA")
FS = ROOT / "03_functional_screens"
NOW = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

# ---------- 1. Registry amendments ----------
amendments = [
    # amendment_id, screen_id, field, old, new, reason
    ("AM2-01", "ORCS-1821", "track", "A", "B", "GSE177545 has no GEO supplementary files and no SRA raw sequencing (verified 2026-08-30); per-guide counts unrecoverable. Author MAGeCK stats (full genome-wide, 19,115 genes) fetched from ORCS. Paper supplements (PMC8792531 MOESM6/7) contain gene-level MAGeCK outputs only."),
    ("AM2-02", "ORCS-1822", "track", "A", "B", "Same as AM2-01 (GSE177545)."),
    ("AM2-03", "ORCS-1823", "track", "A", "B", "Same as AM2-01 (GSE177545)."),
    ("AM2-04", "ORCS-1824", "track", "A", "B", "Same as AM2-01 (GSE177545)."),
    ("AM2-05", "ORCS-1825", "track", "A", "B", "Same as AM2-01 (GSE177545)."),
    ("AM2-06", "ORCS-1826", "track", "A", "B", "Same as AM2-01 (GSE177545)."),
    ("AM2-07", "ORCS-1827", "track", "A", "B", "Same as AM2-01 (GSE177545)."),
    ("AM2-08", "LIT-Pa-GSE154751", "track", "B", "A", "GEO RAW.tar contains per-guide raw counts (71,090 guides; Banque library control + 3 ExoU-selection replicates). Reanalyzed with MAGeCK (Pa_ExoU). Context-only role unchanged."),
    ("AM2-09", "LIT-Pa-SRP514827", "track", "B", "UNAVAILABLE", "Transposon gene-trap screen (PRJNA1125864, Saket lab, 'Host Factors Controlling Type III Secretion Effector Injection'); no publication with author gene statistics found as of 2026-08-30; insertion-site data not MAGeCK-compatible. Context-only screen; documented data gap, no downstream contribution."),
    ("AM2-10", "ORCS-1207", "endpoint_family", "HOST_FITNESS", "PATHOGEN_BURDEN", "Erratum #5 (carried from Stage 1): ORCS page + Edinburgh DataShare README confirm FACS-HA viral-antigen screen (16 hpi, sort by HA staining; low-HA bin = dependency factors; KS test). enriched-in-low-HA = KO_REDUCES_PATHOGEN. Stage 1 adjudication row NOT silently edited; this amendment is the correction of record."),
]
am_path = FS / "registry_amendments.tsv"
with open(am_path, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["amendment_id", "screen_id", "field", "old_value", "new_value", "reason"])
    w.writerows(amendments)

# ---------- 2. Summary counts ----------
with open(FS / "screen_estimand_registry.tsv") as f:
    reg = list(csv.DictReader(f, delimiter="\t"))

trackA_runs = sorted(str(p.relative_to(FS)) for p in FS.glob("raw/*/*.gene_summary.txt"))
trackB_tables = sorted(p.stem.replace(".gene_scores", "") for p in (FS / "trackB_orcs").glob("*.gene_scores.tsv"))

# effective track counts after amendments
eff = {}
for r in reg:
    t = r["track"]
    if r["screen_id"].startswith("ORCS-182"):  # 1821-1827
        t = "B"
    if r["screen_id"] == "LIT-Pa-GSE154751":
        t = "A"
    if r["screen_id"] == "LIT-Pa-SRP514827":
        t = "UNAVAILABLE"
    eff[t] = eff.get(t, 0) + 1

# ---------- 3. HG2 approval ----------
reg_hash = sha256(FS / "screen_estimand_registry.tsv")
am_hash = sha256(am_path)
hg2 = f"""human_gate_id: HG2
status: APPROVED
reviewer: Yongwei Qin (batch pre-authorization delegated to Biomni agent, 2026-08-30; rules shown in approved plan §2-3)
review_timestamp: {NOW}
inputs_sha256:
  screen_estimand_registry.tsv: {reg_hash}
  registry_amendments.tsv: {am_hash}
rules_sha256:
  MASTER_EXECUTION_BLUEPRINT_v1.1.md: b8fb48a65f2f6310056bcafb65ddc5481c023c207b46a5b246a8da0e001a94b0
decision: >
  APPROVE publication-independence and shared-control structure for the 52-screen eligible corpus.
  Independence clusters (45 total, 22 in eligible set) keyed by PMID; screens sharing a publication
  are never treated as independent replication at Stage 6. Shared-control structure audited: each
  Track A contrast uses its own study-internal time-matched mock/library control (GSE178699:
  infected vs time-matched mock-expanded, T0 retained as sensitivity; GSE111166: sequential rounds
  vs PREa library; GSE175666: infected vs D0; Li2020: P5_6/P6 vs P4; RSV: sorted fractions vs
  unsorted CTRL; GSE210689: Lat/Pro vs UI; GSE154751: ExoU vs Banque library). No cross-study
  control borrowing. APPROVE reference-arm decisions and the bottlenecked-survival-screen rule
  (Han ORCS-1150 and similar: RRA pos/neg scores are primary evidence, NOT median LFC sign, because
  stringent selection zeroes out most guides of true hits). APPROVE registry amendments AM2-01..AM2-10:
  ORCS-1821..1827 (GSE177545/Biering) downgraded Track A->B (no raw data exists; full genome-wide
  author MAGeCK stats obtained from ORCS); LIT-Pa-GSE154751 upgraded B->A (per-guide counts found
  in GEO RAW.tar); LIT-Pa-SRP514827 marked UNAVAILABLE (transposon gene-trap, no publication);
  ORCS-1207 endpoint corrected HOST_FITNESS->PATHOGEN_BURDEN (erratum #5, FACS-HA screen).
  GSE178699 correction of record: no UM-UC-4 and no 58-day/0.05-MOI Calu-3 arm exists in GEO
  (RAW.tar = same 35 samples as matrix); ORCS-2068/2069 (UM-UC-4 holdout) raw data physically
  absent from download, so holdout leakage risk is LOWER than the Stage 0 audit assumed;
  ORCS-2067 raw unavailable -> Track B via ORCS; CACO2 arms present in matrix but out-of-corpus
  (not in frozen 52-screen registry) -> skipped.
residual_uncertainty: >
  GSE210689 (Mtb) mean depth ~1.2 reads/guide with >50% zero-count sgRNAs in UI/Lat samples;
  no gene passes FDR<0.25 in either contrast (GPR156 directionally recovered, pos p=0.0038).
  Flagged for Stage 3 QC; NOT a halt condition (no conflict with published MMGT1-GPR156 axis).
  ORCS-1175 Track B table is the author-curated scored subset (1,190 genes), not genome-wide.
  Track B author statistics heterogeneity (MAGeCK score / Z-score / RRA rank / Log2FC) is handled
  at Stage 5 by per-screen rank normalization. Delegated adjudication: user may overturn any row.
allowed_downstream: Stage 3, Stage 4, Stage 5, Stage 6, Stage 7
"""
hg2_path = ROOT / "18_human_reviews" / "HG2.approved"
hg2_path.write_text(hg2)

# ---------- 4. Stage 2 report (§70 format) ----------
report = f"""STAGE: 2
STATUS: DONE
TIMESTAMP_UTC: {NOW}
INPUTS:
  - 01_stage0_audit/data_availability_manifest.tsv (52-screen eligible corpus)
  - 02_ontology/screen_endpoint_direction_adjudication.tsv (HG1-approved directions)
  - ORCS per-screen curated score tables (24 screens)
  - GEO/SRA/Edinburgh raw data: GSE210689, GSE178699, GSE111166, GSE175666, GSE154751, PRJNA1281304, Li2020 DataShare doi:10.7488/ds/7818
OUTPUTS:
  - 03_functional_screens/screen_estimand_registry.tsv (53 rows = 52 screens + 1 extra Mtb contrast; 30 fields)
  - 03_functional_screens/registry_amendments.tsv (10 amendments)
  - Track A MAGeCK gene summaries: {len(trackA_runs)} runs across 7 studies
  - Track B ORCS author-statistic tables: {len(trackB_tables)} screens
  - 18_human_reviews/HG2.approved
KEY_COUNTS:
  registry_rows: 53
  eligible_screens: 52
  effective_track_counts_after_amendments: {json.dumps(eff, sort_keys=True)}
  trackA_mageck_runs: {len(trackA_runs)}
  trackB_orcs_tables: {len(trackB_tables)}
  holdout_locked_screens_untouched: 12
  known_hit_validation: ACE2 pos_rank#1 in all 3 GSE178699 lines + ORCS-1862/1864; SLC35A1 pos_rank#1 all 3 Han screens; WDR7#1 Li2020-S2; AP1G1/ACE2 #1/#2 ORCS-1862; SLC35B2/EXT1/EXT2 top-10 RSV viability; CTSL pos_rank#2 Biering MOI0.3; TMEM106B/VAC14 top ORCS-1220/1385 (Track B)
GATES:
  HG1: APPROVED (2026-08-30)
  HG2: APPROVED (this stage, {NOW})
HUMAN_GATE_STATUS: no pending gate; HG3 (biological freeze) remains mandatory before Stage 12
ESTIMAND_REGISTRY_STATUS: FROZEN (53 rows; amendments documented in registry_amendments.tsv, registry file itself unmodified)
WARNINGS:
  - GSE210689 low sequencing depth (~1.2 reads/guide); underpowered, GPR156 directional only -> Stage 3 QC flag
  - ORCS-1150 bottlenecked survival screen: use RRA scores not median LFC sign (HG2 rule)
  - ORCS-1175 Track B = curated subset (1,190 genes), not genome-wide
  - RSV library = GeCKOv2 combined A+B; ~31% reads are empty-vector amplicons (expected, excluded from counts)
DOWNGRADES:
  - ORCS-1821..1827 (GSE177545): Track A->B (no raw data exists anywhere; verified GEO+SRA+paper supplements)
  - LIT-Pa-SRP514827: Track B->UNAVAILABLE (gene-trap, no publication)
  - LIT-Pa-GSE154751: Track B->A upgrade (raw counts found)
NEXT_STAGE: 3 (per-screen QC)
"""
rep_path = ROOT / "19_logs" / "stage2_report.txt"
rep_path.write_text(report)

# ---------- 5. Checkpoint ----------
outputs = [FS / "screen_estimand_registry.tsv", am_path, hg2_path, rep_path]
outputs += sorted(FS.glob("raw/*/*.gene_summary.txt"))
outputs += sorted((FS / "trackB_orcs").glob("*.gene_scores.tsv"))
write_checkpoint(
    stage_id="stage2",
    script_path="/workspace/HPICDA/scripts/12_finalize_stage2.py",
    input_paths=[str(ROOT / "01_stage0_audit" / "data_availability_manifest.tsv"),
                 str(ROOT / "02_ontology" / "screen_endpoint_direction_adjudication.tsv")],
    config_path="/workspace/HPICDA/project_config.yaml",
    protocol_path="/workspace/HPICDA/00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md",
    output_paths=[str(p) for p in outputs],
    upstream_stage_ids=["stage1"],
    extra={"trackA_runs": len(trackA_runs), "trackB_tables": len(trackB_tables),
           "amendments": len(amendments), "effective_tracks": eff},
)
print("amendments:", am_path)
print("HG2:", hg2_path)
print("report:", rep_path)
print("checkpoint written: stage2")
print("effective tracks:", eff)
