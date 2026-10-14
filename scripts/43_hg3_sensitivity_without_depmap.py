#!/usr/bin/env python3
"""专家意见 #2:HG3 去 DepMap 敏感性分析(2026-09-04)。

关键事实(先在输出中披露):
  HG3 评分细则 v2 中,DepMap 衍生输入只有 therapy_window_penalty
  (相关背景最小基线 Chronos,仅对 INHIBIT 方向施加,取值 0/-1/-3,即只罚不加)。
  Stage 11 的 DepMap 共必需性仅用于模块级解释,从未给单基因加分。
  因此不存在"DepMap 抬高候选"的方向;本敏感性分析回答的是:
  哪些候选被 DepMap 必需性罚分压住,去除该分项后 Tier 是否改变。

方法:total_nodepmap = total_score - therapy_window_penalty(即罚分归零),
层级规则与冻结细则 v2 逐字相同。冻结文件不修改。

输出:
  18_human_reviews/HG3_sensitivity_without_depmap.tsv
  18_human_reviews/HG3_sensitivity_without_depmap_summary.json
"""
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "18_human_reviews"
sys.path.insert(0, str(ROOT / "runner"))
from checkpoint import write_checkpoint

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

sc = pd.read_csv(OUT / "HG3_candidate_scoring_v2.tsv", sep="\t")
fz = pd.read_csv(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", sep="\t",
                 usecols=["gene", "bio_tier"])
hs_genes = set(pd.read_csv(ROOT / "20_stage13_holdout_adjudication/stage13_gene_holdout_verdicts.tsv",
                           sep="\t", usecols=["gene", "is_hs"])
               .query("is_hs == True")["gene"])

def tier_of(total, screen, contra):
    if total < 2.0 or (contra and total < 6.0):
        return "excluded"
    if total >= 6.0 and screen >= 3.0 and not contra:
        return "Tier1"
    if total >= 3.5 and screen >= 1.0:
        return "Tier2"
    if total >= 2.0:
        return "watchlist"
    return "excluded"

sc["total_nodepmap"] = (sc["total_score"] - sc["therapy_window_penalty"]).round(2)
sc["tier_nodepmap"] = [tier_of(t, s, c) for t, s, c in
                       zip(sc["total_nodepmap"], sc["screen_evidence_pts"], sc["contradictory_any"])]
sc["tier_changed"] = sc["tier_nodepmap"] != sc["hg3_tier"]
sc = sc.merge(fz, on="gene", how="left")
sc["is_hs"] = sc["gene"].isin(hs_genes)

# 排序键与评分脚本一致
sc = sc.sort_values(["total_score", "screen_evidence_pts"], ascending=False).reset_index(drop=True)
sc["rank_v2"] = sc.index + 1
sc["rank_nodepmap"] = (sc["total_nodepmap"].rank(method="first", ascending=False)
                       .astype(int))
top10_v2 = set(sc.nsmallest(10, "rank_v2")["gene"])
top10_nd = set(sc.nsmallest(10, "rank_nodepmap")["gene"])

cols = ["gene", "bio_tier", "hg3_tier", "tier_nodepmap", "tier_changed",
        "total_score", "therapy_window_penalty", "total_nodepmap", "score_delta",
        "rank_v2", "rank_nodepmap", "screen_evidence_pts", "min_baseline_chronos",
        "frac_backgrounds_essential", "intervention_direction", "is_hs"]
sc["score_delta"] = (sc["total_nodepmap"] - sc["total_score"]).round(2)
sc[cols].to_csv(OUT / "HG3_sensitivity_without_depmap.tsv", sep="\t", index=False)

t1 = sc[sc["hg3_tier"] == "Tier1"]
hs = sc[sc["is_hs"]]
changed = sc[sc["tier_changed"]]
summary = {
    "analysis": "HG3 sensitivity without DepMap (expert item #2)",
    "date": "2026-09-04",
    "framing_disclosure": ("HG3 细则 v2 中 DepMap 仅经 therapy_window_penalty 进入评分"
                           "(0/-1/-3,只罚不加,仅 INHIBIT 方向);Stage 11 DepMap 共必需性"
                           "仅用于模块级解释,从未给单基因加分。因此无候选被 DepMap 抬高;"
                           "本分析检验罚分移除后 Tier 稳定性。"),
    "n_candidates": len(sc),
    "n_with_depmap_penalty": int((sc["therapy_window_penalty"] < 0).sum()),
    "n_tier_changed": int(sc["tier_changed"].sum()),
    "tier1_stability": {"n_tier1": len(t1),
                        "n_remain_tier1": int((t1["tier_nodepmap"] == "Tier1").sum()),
                        "rate": float((t1["tier_nodepmap"] == "Tier1").mean())},
    "top10_stability": {"top10_v2": sorted(top10_v2),
                        "n_remain_top10": len(top10_v2 & top10_nd),
                        "rate": len(top10_v2 & top10_nd) / 10},
    "hs_cluster": {"genes": sorted(hs["gene"]),
                   "n_tier_changed": int(hs["tier_changed"].sum()),
                   "all_stable": bool(~hs["tier_changed"].any())},
    "single_gene_check": {g: {"tier_v2": r["hg3_tier"], "tier_nodepmap": r["tier_nodepmap"],
                              "tw_pen": r["therapy_window_penalty"],
                              "total_v2": r["total_score"], "total_nodepmap": r["total_nodepmap"]}
                          for g, r in sc[sc["gene"].isin(["FER", "PTK2B"])].set_index("gene").iterrows()},
    "tier_changes": changed[["gene", "hg3_tier", "tier_nodepmap", "therapy_window_penalty",
                             "total_score", "total_nodepmap"]].to_dict("records"),
    "input_sha256": {"HG3_candidate_scoring_v2": sha(OUT / "HG3_candidate_scoring_v2.tsv"),
                     "freeze_v2": sha(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv")},
}
json.dump(summary, open(OUT / "HG3_sensitivity_without_depmap_summary.json", "w"),
          indent=2, ensure_ascii=False)

write_checkpoint("hg3_sensitivity_nodepmap",
                 str(ROOT / "scripts/43_hg3_sensitivity_without_depmap.py"),
                 [str(OUT / "HG3_candidate_scoring_v2.tsv"),
                  str(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv")],
                 str(ROOT / "18_human_reviews/HG3_scoring_rubric_v2.yaml"),
                 str(ROOT / "00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md"),
                 [str(OUT / "HG3_sensitivity_without_depmap.tsv"),
                  str(OUT / "HG3_sensitivity_without_depmap_summary.json")],
                 ["hg3_amendment_v2_rescore"])

print(f"候选 {len(sc)};带 DepMap 罚分 {summary['n_with_depmap_penalty']};Tier 变化 {summary['n_tier_changed']}")
print(f"Tier1 稳定 {summary['tier1_stability']['n_remain_tier1']}/{summary['tier1_stability']['n_tier1']};"
      f"Top10 稳定 {summary['top10_stability']['n_remain_top10']}/10;"
      f"HS 簇变化 {summary['hs_cluster']['n_tier_changed']}")
print("Tier 变化明细:")
print(changed[["gene", "hg3_tier", "tier_nodepmap", "therapy_window_penalty",
               "total_score", "total_nodepmap"]].to_string(index=False))
