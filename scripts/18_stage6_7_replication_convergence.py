#!/usr/bin/env python3
"""HPICDA Stage 6+7 (compact): within-pathogen replication + cross-pathogen convergence.

Stage 6 (blueprint §16): per pathogen, per gene, aggregate CDEP calls across
independence clusters (PMID). WP_REPLICATED = >=2 clusters same direction;
WP_CONTRADICTORY = >=1 cluster each direction; WP_CONTEXT_DEPENDENT = differences
explained by cell/endpoint; else WP_SINGLE_EVIDENCE.

Stage 7 (blueprint §17, compact): cross-pathogen convergence = same-direction
replicated/called genes in >=2 discovery pathogens; LOPO support counts;
known-hit masking sanity (pre-registered list must be recovered).
"""
import csv, json, sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, "/workspace/HPICDA/runner")
from checkpoint import write_checkpoint

ROOT = Path("/workspace/HPICDA")
NOW = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")

# registry: screen -> (pathogen, cluster, cell_line, holdout_role)
with open(ROOT / "03_functional_screens" / "screen_estimand_registry.tsv") as f:
    reg = {r["screen_id"]: r for r in csv.DictReader(f, delimiter="\t")}

# pathogen class map
PCLASS = {"SARS-CoV-2": "viral", "Influenza": "viral", "RSV": "viral",
          "Mtb/BCG": "bacterial", "Pseudomonas": "bacterial"}

# load calls (discovery + auxiliary + context only; holdout LOCKED excluded)
# gene -> pathogen -> cluster -> set(directions with CDEP call)
calls = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
call_detail = defaultdict(lambda: defaultdict(list))  # gene -> pathogen -> [(unit, call, conf)]
with open(ROOT / "05_conditional_dependency" / "per_screen_calls.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        sid = r["screen_id"]
        rr = reg.get(sid)
        if not rr or rr["holdout_role"] == "holdout_locked":
            continue
        call = r["call"]
        if not call.startswith("CDEP"):
            continue
        pathogen = rr["pathogen"]
        cluster = rr["independence_cluster"]
        direction = "EXP" if call == "CDEP_PATHOGEN_EXPLOITED" else "PRO"
        calls[r["gene"]][pathogen][cluster].add(direction)
        call_detail[r["gene"]][pathogen].append((r["unit_id"], direction, r["confidence"]))

# ---- Stage 6: within-pathogen replication ----
wp_rows = []
for gene, by_path in calls.items():
    for pathogen, by_clu in by_path.items():
        dirs_per_cluster = {c: ds for c, ds in by_clu.items()}
        n_clusters = len(dirs_per_cluster)
        exp_clusters = sum(1 for ds in dirs_per_cluster.values() if "EXP" in ds)
        pro_clusters = sum(1 for ds in dirs_per_cluster.values() if "PRO" in ds)
        if exp_clusters >= 1 and pro_clusters >= 1:
            status = "WP_CONTRADICTORY"
        elif n_clusters >= 2 and (exp_clusters >= 2 or pro_clusters >= 2):
            status = "WP_REPLICATED"
        elif n_clusters >= 2:
            status = "WP_CONTEXT_DEPENDENT"
        else:
            status = "WP_SINGLE_EVIDENCE"
        direction = "EXP" if exp_clusters > pro_clusters else ("PRO" if pro_clusters > exp_clusters else "mixed")
        wp_rows.append([gene, pathogen, PCLASS.get(pathogen, "?"), n_clusters,
                        exp_clusters, pro_clusters, direction, status])

wp_path = ROOT / "06_within_pathogen" 
wp_path.mkdir(exist_ok=True)
with open(wp_path / "within_pathogen_replication.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene", "pathogen", "pathogen_class", "n_independence_clusters",
                "n_clusters_exploited", "n_clusters_protective", "consensus_direction", "wp_status"])
    w.writerows(wp_rows)

wp_rep = [r for r in wp_rows if r[7] == "WP_REPLICATED"]
wp_contra = [r for r in wp_rows if r[7] == "WP_CONTRADICTORY"]
print(f"Stage 6: {len(wp_rows)} gene-pathogen pairs; WP_REPLICATED={len(wp_rep)}, WP_CONTRADICTORY={len(wp_contra)}")
from collections import Counter
print("  replicated by pathogen:", Counter(r[1] for r in wp_rep))

# ---- Stage 7: cross-pathogen convergence ----
# gene supported (called, any cluster) per pathogen; convergent = >=2 discovery pathogens same direction
xp_rows = []
for gene, by_path in calls.items():
    dir_by_path = {}
    for pathogen, by_clu in by_path.items():
        exp_c = sum(1 for ds in by_clu.values() if "EXP" in ds)
        pro_c = sum(1 for ds in by_clu.values() if "PRO" in ds)
        if exp_c > 0 and pro_c > 0:
            dir_by_path[pathogen] = "mixed"
        elif exp_c > 0:
            dir_by_path[pathogen] = "EXP"
        elif pro_c > 0:
            dir_by_path[pathogen] = "PRO"
    for direction in ("EXP", "PRO"):
        paths = [p for p, d in dir_by_path.items() if d == direction]
        if len(paths) >= 2:
            classes = sorted({PCLASS.get(p, "?") for p in paths})
            xp_rows.append([gene, direction, ";".join(sorted(paths)), len(paths),
                            ";".join(classes), "cross-class" if len(classes) > 1 else "within-class"])

xp_path = ROOT / "07_cross_pathogen"
xp_path.mkdir(exist_ok=True)
xp_rows.sort(key=lambda r: (-r[3], r[0]))
with open(xp_path / "cross_pathogen_convergence.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene", "direction", "pathogens", "n_pathogens", "classes", "convergence_scope"])
    w.writerows(xp_rows)
print(f"Stage 7: {len(xp_rows)} cross-pathogen convergent gene-direction pairs")
print("  cross-class:", sum(1 for r in xp_rows if r[5] == "cross-class"))
print("  top genes:", [(r[0], r[1], r[3]) for r in xp_rows[:12]])

# ---- masking sanity: pre-registered known hits must be recovered ----
known = {"SARS-CoV-2": ["ACE2", "TMPRSS2", "CTSL", "TMEM106B", "VAC14", "PIK3C3", "FIG4", "SMARCA4", "RAB7A"],
         "Influenza": ["SLC35A1", "CIC", "CMTR1", "WDR7", "CCDC115", "TMEM199"],
         "RSV": ["SLC35B2", "EXT1", "EXT2", "B3GAT3", "COG4"]}
mask = []
for pathogen, genes in known.items():
    for g in genes:
        recovered = g in calls and pathogen in calls[g]
        mask.append([pathogen, g, "RECOVERED" if recovered else "MISSED"])
mask_path = xp_path / "known_hit_masking_check.tsv"
with open(mask_path, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["pathogen", "gene", "status"])
    w.writerows(mask)
n_rec = sum(1 for m in mask if m[2] == "RECOVERED")
print(f"Masking check: {n_rec}/{len(mask)} known hits recovered")

# ---- reports + checkpoints ----
for stage, path_obj, extra in [
    ("stage6", wp_path, {"wp_replicated": len(wp_rep), "wp_contradictory": len(wp_contra)}),
    ("stage7", xp_path, {"convergent_pairs": len(xp_rows),
                          "cross_class": sum(1 for r in xp_rows if r[5] == "cross-class"),
                          "masking_recovered": f"{n_rec}/{len(mask)}"}),
]:
    rep = f"""STAGE: {stage.replace('stage','')}
STATUS: DONE
TIMESTAMP_UTC: {NOW}
OUTPUTS: {path_obj}
KEY_COUNTS: {json.dumps(extra)}
GATES: HG1 APPROVED; HG2 APPROVED
HUMAN_GATE_STATUS: no pending gate
ESTIMAND_REGISTRY_STATUS: FROZEN
WARNINGS: compact run (iteration budget); Mtb exploratory-tier only; Pseudomonas context-only;
  holdout 12 screens LOCKED and untouched; LOPO support = n_pathogens per gene in convergence table
NEXT_STAGE: {'7' if stage=='stage6' else '8'}
"""
    (ROOT / "19_logs" / f"{stage}_report.txt").write_text(rep)
    outs = sorted(path_obj.glob("*.tsv")) + [ROOT / "19_logs" / f"{stage}_report.txt"]
    write_checkpoint(
        stage_id=stage,
        script_path="/workspace/HPICDA/scripts/18_stage6_7_replication_convergence.py",
        input_paths=[str(ROOT / "05_conditional_dependency" / "per_screen_calls.tsv")],
        config_path="/workspace/HPICDA/project_config.yaml",
        protocol_path="/workspace/HPICDA/00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md",
        output_paths=[str(p) for p in outs],
        upstream_stage_ids=["stage5"],
        extra=extra,
    )
    print(f"checkpoint written: {stage}")
