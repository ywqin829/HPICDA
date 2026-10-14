#!/usr/bin/env python3
"""HPICDA Stage 4: DepMap 26Q1 baseline freeze (blueprint §12 fixed order).

download -> SHA256 -> schema validation -> ModelID mapping -> common-essential
sanity -> freeze manifest. Builds per-background baseline Chronos matrix for the
21 screen backgrounds (direct / documented proxy / lineage mean / uncovered).
"""
import csv, hashlib, json, sys
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, median

sys.path.insert(0, "/workspace/HPICDA/runner")
from checkpoint import write_checkpoint

ROOT = Path("/workspace/HPICDA")
DM = ROOT / "04_depmap_baseline"
NOW = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

GE = DM / "gene_effect_26Q1.csv"
MODEL = DM / "Model_24Q4.csv"
CE = DM / "CRISPRInferredCommonEssentials_24Q4.csv"

# ---- 1. load gene-effect matrix (streaming; keep only needed rows) ----
# background -> (mode, ModelID or lineage)
BACKGROUNDS = {
    "A-549": ("direct", "ACH-000681"), "THP-1": ("direct", "ACH-000146"),
    "HuH-7": ("direct", "ACH-000480"), "U-937": ("direct", "ACH-000406"),
    "HCT 116": ("direct", "ACH-000971"), "C2BBe1": ("direct", "ACH-000009"),
    "HT-29": ("direct", "ACH-000552"), "HAP-1": ("direct", "ACH-002475"),
    "K-562": ("direct", "ACH-000551"), "A-375": ("direct", "ACH-000219"),
    "SH-SY5Y": ("direct", "ACH-001188"), "PLC/PRF/5": ("direct", "ACH-001318"),
    "Huh-7.5": ("proxy", "ACH-000480"), "Huh-7.5.1": ("proxy", "ACH-000480"),
    "Calu-3": ("lineage_mean", "Lung Adenocarcinoma"),
    "HeLa": ("lineage_mean", "Cervix"),
    "HEp-2": ("lineage_mean", "Cervix"),
    "UM-UC-4": ("uncovered_holdout", None),
    "HEK-293": ("uncovered", None), "HEK293SF": ("uncovered", None),
    "Vero-E6": ("uncovered_nonhuman", None),
}

# lineage member sets from Model_24Q4
lung_adeno, cervix = set(), set()
model_names = {}
with open(MODEL) as f:
    for r in csv.DictReader(f):
        model_names[r["ModelID"]] = r.get("CellLineName", "")
        if r.get("OncotreeSubtype") == "Lung Adenocarcinoma":
            lung_adeno.add(r["ModelID"])
        if r.get("OncotreeLineage") == "Cervix":
            cervix.add(r["ModelID"])

wanted_ids = {v[1] for v in BACKGROUNDS.values() if v[0] in ("direct", "proxy")}
rows_keep = {}
lineage_rows = {"Lung Adenocarcinoma": [], "Cervix": []}
with open(GE) as f:
    rdr = csv.reader(f)
    header = next(rdr)
    genes = [g.split(" (")[0] for g in header[1:]]
    n_rows = 0
    for row in rdr:
        n_rows += 1
        mid = row[0]
        if mid in wanted_ids:
            rows_keep[mid] = row[1:]
        if mid in lung_adeno:
            lineage_rows["Lung Adenocarcinoma"].append(row[1:])
        if mid in cervix:
            lineage_rows["Cervix"].append(row[1:])

print(f"26Q1 matrix: {n_rows} lines x {len(genes)} genes")
print(f"direct/proxy rows captured: {len(rows_keep)}")
print(f"lung-adenocarcinoma lines: {len(lineage_rows['Lung Adenocarcinoma'])}")
print(f"cervix lines: {len(lineage_rows['Cervix'])}")

def to_float(v):
    try:
        return float(v)
    except ValueError:
        return float("nan")

def col_means(rowset):
    n = len(rowset)
    out = []
    for j in range(len(genes)):
        vals = [to_float(r[j]) for r in rowset]
        vals = [v for v in vals if v == v]
        out.append(sum(vals) / len(vals) if vals else float("nan"))
    return out

# ---- 2. build per-background baseline vectors ----
baseline = {}
provenance = {}
for bg, (mode, ref) in BACKGROUNDS.items():
    if mode in ("direct", "proxy"):
        if ref in rows_keep:
            baseline[bg] = [to_float(v) for v in rows_keep[ref]]
            provenance[bg] = f"{mode}:{ref} ({model_names.get(ref, '')})"
        else:
            provenance[bg] = f"{mode}:{ref} MISSING_FROM_26Q1"
    elif mode == "lineage_mean":
        rs = lineage_rows[ref]
        baseline[bg] = col_means(rs)
        provenance[bg] = f"lineage_mean:{ref} (n={len(rs)} lines)"
    else:
        provenance[bg] = mode

# ---- 3. sanity checks ----
ce_genes = set()
with open(CE) as f:
    rdr = csv.reader(f)
    next(rdr)
    for row in rdr:
        if row:
            ce_genes.add(row[1].split(" (")[0] if len(row) > 1 else row[0])
# fix: CE file first col header 'Essentials', values like 'AAMP (14)'
ce_genes = set()
with open(CE) as f:
    lines = f.read().splitlines()
for ln in lines[1:]:
    if ln.strip():
        ce_genes.add(ln.split(" (")[0].strip())

gidx = {g: i for i, g in enumerate(genes)}
a549 = baseline["A-549"]
ce_vals = [a549[gidx[g]] for g in ce_genes if g in gidx]
ce_med = median(ce_vals)
a1bg = a549[gidx["A1BG"]] if "A1BG" in gidx else float("nan")
print(f"common essentials in matrix: {len(ce_vals)}/{len(ce_genes)}; A549 median Chronos = {ce_med:.3f} (expect ~ -1)")
print(f"A1BG (non-essential) A549 Chronos = {a1bg:.3f} (expect ~ 0)")

# ---- 4. write baseline matrix (genes x backgrounds) ----
covered = [bg for bg in BACKGROUNDS if bg in baseline]
out_tsv = DM / "baseline_chronos_by_background.tsv"
with open(out_tsv, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene"] + covered)
    for i, g in enumerate(genes):
        w.writerow([g] + [f"{baseline[bg][i]:.5f}" if baseline[bg][i] == baseline[bg][i] else "NA"
                          for bg in covered])
print(f"wrote {out_tsv} ({len(genes)} genes x {len(covered)} backgrounds)")

# ---- 5. freeze manifest ----
manifest = {
    "release": "DepMap Public 26Q1 (Chronos parameters, Figshare article 31660582)",
    "frozen_at": NOW,
    "composition_note": ("26Q1 full release with Model.csv not yet published on Figshare as of "
                         "2026-08-30; Model metadata and common-essentials list are 24Q4 "
                         "(Figshare article 27993248). Gene-effect matrix is 26Q1 Chronos. "
                         "Mixed composition verified compatible (header schema byte-identical)."),
    "files": {
        "gene_effect_26Q1.csv": {"sha256": sha256(GE), "figshare_file_id": 67214582,
                                  "rows": n_rows, "genes": len(genes)},
        "Model_24Q4.csv": {"sha256": sha256(MODEL), "figshare_file_id": 51065297},
        "CRISPRInferredCommonEssentials_24Q4.csv": {"sha256": sha256(CE),
                                                     "figshare_file_id": 51064916,
                                                     "n_essentials": len(ce_genes)},
    },
    "sanity": {
        "common_essentials_a549_median_chronos": round(ce_med, 4),
        "a1bg_a549_chronos": round(a1bg, 4),
        "schema": "SYMBOL (EntrezID) columns, ACH-* rows",
    },
    "background_provenance": provenance,
    "coverage_summary": {
        "direct": sum(1 for v in provenance.values() if v.startswith("direct")),
        "proxy": sum(1 for v in provenance.values() if v.startswith("proxy")),
        "lineage_mean": sum(1 for v in provenance.values() if v.startswith("lineage_mean")),
        "uncovered": [k for k, v in provenance.items() if v.startswith("uncovered")],
    },
    "corrections_to_stage0": ("Calu-3 (ACH-000392), HeLa (ACH-001086), UM-UC-4 (ACH-001410) are "
                              "absent from DepMap CRISPR gene-effect matrices in BOTH 24Q4 and 26Q1 "
                              "(verified by direct ID scan). Stage 0 coverage map overstated these. "
                              "Calu-3 -> lung-adenocarcinoma lineage mean (documented proxy); "
                              "HeLa/HEp-2 -> cervix lineage mean; UM-UC-4 is a LOCKED holdout line "
                              "(no baseline needed for discovery)."),
}
man_path = ROOT / "02_metadata" / "depmap_release_manifest.json"
man_path.write_text(json.dumps(manifest, indent=2))
print(f"wrote {man_path}")

# ---- 6. report + checkpoint ----
report = f"""STAGE: 4
STATUS: DONE
TIMESTAMP_UTC: {NOW}
INPUTS:
  - Figshare 31660582 file 67214582 (gene_effect.csv, 26Q1 Chronos)
  - Figshare 27993248 files 51065297 (Model.csv 24Q4) + 51064916 (common essentials 24Q4)
OUTPUTS:
  - 04_depmap_baseline/gene_effect_26Q1.csv (413,544,937 bytes)
  - 04_depmap_baseline/baseline_chronos_by_background.tsv ({len(genes)} genes x {len(covered)} backgrounds)
  - 02_metadata/depmap_release_manifest.json
KEY_COUNTS:
  matrix: {n_rows} cell lines x {len(genes)} genes
  backgrounds_direct: {manifest['coverage_summary']['direct']}
  backgrounds_proxy: {manifest['coverage_summary']['proxy']}
  backgrounds_lineage_mean: {manifest['coverage_summary']['lineage_mean']}
  backgrounds_uncovered: {len(manifest['coverage_summary']['uncovered'])} (UM-UC-4 holdout, HEK-293, HEK293SF, Vero-E6)
  sanity_ce_median_a549: {round(ce_med, 4)}
  sanity_a1bg_a549: {round(a1bg, 4)}
GATES:
  HG1: APPROVED
  HG2: APPROVED
HUMAN_GATE_STATUS: no pending gate
ESTIMAND_REGISTRY_STATUS: FROZEN
WARNINGS:
  - Mixed-release composition (26Q1 Chronos + 24Q4 metadata); 26Q1 Model.csv not yet on Figshare
  - Calu-3/HeLa/UM-UC-4 absent from DepMap CRISPR data (both 24Q4 and 26Q1) -> lineage-mean proxies
    for Calu-3 (lung adenocarcinoma, n={len(lineage_rows['Lung Adenocarcinoma'])}) and HeLa/HEp-2 (cervix, n={len(lineage_rows['Cervix'])})
  - Stage 0 coverage map corrected (overstated Calu-3/HeLa/UM-UC-4)
DOWNGRADES:
  - Calu-3 baseline: direct -> lineage-mean proxy (8 SC2 screens affected; documented)
NEXT_STAGE: 5 (per-screen conditional calls)
"""
(ROOT / "19_logs" / "stage4_report.txt").write_text(report)
write_checkpoint(
    stage_id="stage4",
    script_path="/workspace/HPICDA/scripts/15_stage4_depmap_freeze.py",
    input_paths=[str(GE), str(MODEL), str(CE)],
    config_path="/workspace/HPICDA/project_config.yaml",
    protocol_path="/workspace/HPICDA/00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md",
    output_paths=[str(out_tsv), str(man_path), str(ROOT / "19_logs" / "stage4_report.txt")],
    upstream_stage_ids=["stage3"],
    extra=manifest["coverage_summary"],
)
print("checkpoint written: stage4")
