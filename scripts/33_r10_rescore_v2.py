#!/usr/bin/env python3
"""R10: 按冻结细则 v2 重评分(计划 v3 Workstream A2)。

与 29_hg3_scoring.py 的唯一差异:genetics 分项使用
16_revision_r10_g2/stage10_genetics_evidence_v2.tsv 的 genetics_tier_v2
(仅 infection_specific 表型可达 G2),并输出 v1/v2 差异报告。

输出:
  18_human_reviews/HG3_candidate_scoring_v2.tsv
  16_revision_r10_g2/rescore_v1_vs_v2_diff.tsv
  16_revision_r10_g2/r10_rescore_runinfo.json
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "18_human_reviews"
R10 = ROOT / "16_revision_r10_g2"

INPUTS = [
    "02_metadata/preregistered_candidate_set_v1.tsv",
    "06_within_pathogen/within_pathogen_replication_v2.tsv",
    "07_cross_pathogen/convergence_meta_v2_gene.tsv",
    "07_cross_pathogen/convergence_v2_loso.tsv",
    "04_depmap_baseline/baseline_chronos_by_background.tsv",
    "09_transcriptomic_context/stage8_per_gene_summary.tsv",
    "10_singlecell_context/stage9_candidate_lung_expression_summary.tsv",
    "16_revision_r10_g2/stage10_genetics_evidence_v2.tsv",
    "18_human_reviews/HG3_scoring_rubric_v2.yaml",
]

TIER_PTS = {"T4_CROSS_STUDY_DIFF_SYSTEM": 3.0, "T3_CROSS_STUDY_SAME_SYSTEM": 2.0,
            "T2_SAME_STUDY": 1.0, "T1_GUIDE_SET": 0.5}
SINGLE_PTS = 0.25
SCREEN_CAP = 6.0
CONTRA_PENALTY = -2.0
CONV_WITHIN_CLASS = 2.0
CONV_INTERSECTION = 1.0
HUMAN_CAP = 4.0
REL_BG = ["A-549", "Calu-3", "HuH-7", "Huh-7.5", "Huh-7.5.1", "HEp-2", "HeLa"]

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

cand = pd.read_csv(ROOT / INPUTS[0], sep="\t")
rep = pd.read_csv(ROOT / INPUTS[1], sep="\t")
conv = pd.read_csv(ROOT / INPUTS[2], sep="\t").set_index("gene")
loso = pd.read_csv(ROOT / INPUTS[3], sep="\t")
bl = pd.read_csv(ROOT / INPUTS[4], sep="\t").set_index("gene")
s8 = pd.read_csv(ROOT / INPUTS[5], sep="\t").set_index("gene")
s9 = pd.read_csv(ROOT / INPUTS[6], sep="\t").set_index("gene")
s10 = pd.read_csv(ROOT / INPUTS[7], sep="\t").set_index("gene")

rel = [c for c in REL_BG if c in bl.columns]

def pair_pts(row):
    if row["wp_status"] == "WP_CONTRADICTORY":
        return 0.0
    if row["wp_status"] == "WP_SINGLE_EVIDENCE":
        return SINGLE_PTS
    return TIER_PTS.get(row["replication_tier"], 0.0)

rep["pts"] = rep.apply(pair_pts, axis=1)
scr = (rep.groupby(["gene", "pathogen"])["pts"].max().reset_index()
       .groupby("gene")["pts"].sum().clip(upper=SCREEN_CAP))
contra_genes = set(rep.loc[rep["wp_status"] == "WP_CONTRADICTORY", "gene"])
best_tier = (rep[rep["wp_status"].isin(["WP_REPLICATED", "WP_ROBUST"])]
             .groupby("gene")["replication_tier"]
             .agg(lambda s: min(s, key=lambda t: -TIER_PTS[t])))
nclu = rep.groupby("gene")["n_independence_clusters"].max()
loso_map = {(r["gene"], r["direction"]): r["stability"] for _, r in loso.iterrows()}

rows = []
for _, c in cand.iterrows():
    g = c["gene"]
    direction = "INHIBIT" if c["directions"] == "EXP" else "PRESERVE_ACTIVATE"
    screen = float(scr.get(g, 0.0))
    contra = g in contra_genes
    dir_pen = CONTRA_PENALTY if contra else 0.0
    src = c["sources"]
    conv_bonus = CONV_WITHIN_CLASS if "META_CONVERGENT" in src else (
        CONV_INTERSECTION if "THRESHOLD_INTERSECTION" in src else 0.0)
    if g in bl.index:
        mc = float(bl.loc[g, rel].min())
        frac_ess = float((bl.loc[g, rel] < -1).mean())
    else:
        mc, frac_ess = np.nan, np.nan
    if direction == "INHIBIT" and not np.isnan(mc):
        tw_pen = 0.0 if mc > -1 else (-1.0 if mc >= -2 else -3.0)
    else:
        tw_pen = 0.0
    lung = bool(s9.loc[g, "lung_expressed"]) if g in s9.index else False
    de = bool(s8.loc[g, "clinical_de_concordant"]) if g in s8.index else False
    gtier = s10.loc[g, "genetics_tier_v2"] if g in s10.index else "G0"
    loeuf = float(s10.loc[g, "loeuf"]) if g in s10.index and pd.notna(s10.loc[g, "loeuf"]) else np.nan
    resp_phys = bool(s10.loc[g, "respiratory_physiology_genetics"]) if g in s10.index and pd.notna(s10.loc[g, "respiratory_physiology_genetics"]) else False
    human = min(1.0 * lung + 1.0 * de + {"G1": 2.0, "G2": 1.0}.get(gtier, 0.0), HUMAN_CAP)
    total = screen + dir_pen + conv_bonus + tw_pen + human
    if total < 2.0 or (contra and total < 6.0):
        tier = "excluded"
    elif total >= 6.0 and screen >= 3.0 and not contra:
        tier = "Tier1"
    elif total >= 3.5 and screen >= 1.0:
        tier = "Tier2"
    elif total >= 2.0:
        tier = "watchlist"
    else:
        tier = "excluded"
    rows.append({
        "gene": g, "screen_direction": c["directions"], "intervention_direction": direction,
        "pathogens": c["pathogens"], "sources": src,
        "best_replication_tier": best_tier.get(g, "SINGLE_EVIDENCE" if screen > 0 else "NONE"),
        "n_independence_clusters_max": int(nclu.get(g, 0)),
        "screen_evidence_pts": round(screen, 2),
        "contradictory_any": contra, "direction_penalty": dir_pen,
        "convergence_bonus": conv_bonus,
        "loso_stability": loso_map.get((g, c["directions"]), np.nan),
        "min_baseline_chronos": round(mc, 3) if not np.isnan(mc) else np.nan,
        "frac_backgrounds_essential": round(frac_ess, 3) if not np.isnan(frac_ess) else np.nan,
        "therapy_window_penalty": tw_pen,
        "lung_expressed": lung, "clinical_de_concordant": de,
        "genetics_tier": gtier, "respiratory_physiology_genetics": resp_phys, "loeuf": loeuf,
        "human_layer_bonus": human,
        "total_score": round(total, 2), "hg3_tier": tier,
    })

res = pd.DataFrame(rows).sort_values(["total_score", "screen_evidence_pts"], ascending=False)
res.to_csv(OUT / "HG3_candidate_scoring_v2.tsv", sep="\t", index=False)

# ---- v1 vs v2 差异报告 ----
v1 = pd.read_csv(OUT / "HG3_candidate_scoring_v1.tsv", sep="\t").set_index("gene")
v2 = res.set_index("gene")
diff = v2[["total_score", "hg3_tier", "genetics_tier", "human_layer_bonus"]].join(
    v1[["total_score", "hg3_tier", "genetics_tier", "human_layer_bonus"]],
    rsuffix="_v1")
diff.columns = ["total_v2", "tier_v2", "gtier_v2", "human_v2",
                "total_v1", "tier_v1", "gtier_v1", "human_v1"]
diff["score_delta"] = (diff["total_v2"] - diff["total_v1"]).round(2)
diff["tier_changed"] = diff["tier_v2"] != diff["tier_v1"]
diff = diff.reset_index().sort_values("score_delta")
diff.to_csv(R10 / "rescore_v1_vs_v2_diff.tsv", sep="\t", index=False)

runinfo = {"rubric": "HG3_scoring_rubric_v2.yaml (frozen 2026-09-02, pre-rescore)",
           "n_candidates": len(res),
           "tier_counts_v2": res["hg3_tier"].value_counts().to_dict(),
           "tier_counts_v1": v1["hg3_tier"].value_counts().to_dict(),
           "n_score_changed": int((diff["score_delta"] != 0).sum()),
           "n_tier_changed": int(diff["tier_changed"].sum()),
           "input_sha256": {p: sha(ROOT / p) for p in INPUTS}}
json.dump(runinfo, open(R10 / "r10_rescore_runinfo.json", "w"), indent=2, ensure_ascii=False)

print("tier v2:", runinfo["tier_counts_v2"])
print("tier v1:", runinfo["tier_counts_v1"])
print(f"score changed: {runinfo['n_score_changed']}, tier changed: {runinfo['n_tier_changed']}")
print("\n所有分数变化基因:")
print(diff[diff["score_delta"] != 0][["gene", "gtier_v1", "gtier_v2", "total_v1", "total_v2",
                                      "score_delta", "tier_v1", "tier_v2"]].to_string(index=False))
