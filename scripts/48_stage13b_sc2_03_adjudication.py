#!/usr/bin/env python3
"""Stage 13b:HO-V-SC2-03(PRJEB72772,Stahl et al. 2024)补充裁决(2026-09-04,收尾项)。

该留存条目在 HG4 前已预注册并完成 R14 情境分层(冻结),仅计算被推迟。
本次按冻结规则补齐,不修改 Stage 13 任何冻结输出(故命名 13b,补充性质)。

数据来源与标准(与 ORCS 留存同标准):
  作者基因级结果(boutroslab/Supp_Stahl_2024):Calu3/A549_results_crispr.txt
  作者 HIT 原样:significant 列(Calu-3 FDR<0.10,A549 FDR<0.05,作者各自阈值)
  分数:log2fc;朝向 +1(正=存活富集=依赖因子),锚点 ACE2 +5.34(Calu3)/+4.95(A549)
  原始计数矩阵 SARS-CRISPR-Deep.count.txt 已存档(可复算),本次不重跑 MAGeCK。

冻结规则(与 scripts/40 逐字一致):
  expected_sign: EXP→+1, PRO→−1
  concordant = author_hit 且 sign 一致;reversed = author_hit 且 sign 相反
  裁决: n_rev>0 且 n_con==0 → FAIL(同病原)/XPATH_REVERSED(跨病原)
        n_con>0 → PASS/PASS_XPATH;跨病原非命中 → NOT_EVALUABLE_XPATH
        同病原非命中:context_matched True→FAIL,False→CONTEXT_BOUNDARY,缺失→NOT_STRATIFIED
  情境分层:16_revision_r14_diagnostics/r14_holdout_context_stratification.tsv(冻结,含 HO-V-SC2-03)

输出(23_holdout_sc2_03_stahl/):
  stage13b_sc2_03_gene_screen_calls.tsv / stage13b_sc2_03_verdicts.tsv / stage13b_sc2_03_summary.json
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace/HPICDA")
DATA = ROOT / "23_holdout_sc2_03_stahl"
sys.path.insert(0, str(ROOT / "runner"))
from checkpoint import write_checkpoint

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
HS = ["B3GAT3", "B4GALT7", "EXT1", "EXT2", "EXTL3", "SLC35B2"]
POS_CTRL = ["ACE2", "TMPRSS2"]

cand = pd.read_csv(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv", sep="\t")
cand = cand.rename(columns={"directions": "direction"})
cand["disc_pathogens"] = cand["pathogens"].fillna("").astype(str)

strat = pd.read_csv(ROOT / "16_revision_r14_diagnostics/r14_holdout_context_stratification.tsv", sep="\t")
strat_map = {(r.holdout_id, r.gene): (bool(r.context_matched), bool(r.lab_independent))
             for r in strat.itertuples()}

frames = []
for arm in ["Calu3", "A549"]:
    d = pd.read_csv(DATA / f"{arm}_results_crispr.txt", sep="\t")
    frames.append(pd.DataFrame({
        "gene": d["gene_name"].astype(str).str.strip(),
        "author_hit": d["significant"].astype(str).str.strip().eq("significant"),
        "score": pd.to_numeric(d["log2fc"], errors="coerce"),
        "screen_id": f"STAHL2024-{arm}"}))
calls = pd.concat(frames, ignore_index=True)
calls = calls.merge(cand[["gene", "direction", "disc_pathogens"]], on="gene", how="inner")
calls["expected_sign"] = np.where(calls["direction"] == "EXP", 1,
                                  np.where(calls["direction"] == "PRO", -1, 0))
calls["oriented_sign"] = np.sign(calls["score"]).fillna(0).astype(int)  # orientation = +1
calls["concordant"] = calls["author_hit"] & (calls["oriented_sign"] == calls["expected_sign"]) & (calls["expected_sign"] != 0)
calls["reversed"] = calls["author_hit"] & (calls["oriented_sign"] == -calls["expected_sign"]) & (calls["expected_sign"] != 0)
calls["hit_sign_unknown"] = calls["author_hit"] & (calls["oriented_sign"] == 0)
calls["holdout_id"] = "HO-V-SC2-03"
calls["same_pathogen"] = calls["disc_pathogens"].str.lower().str.contains("sars-cov-2")
calls.to_csv(DATA / "stage13b_sc2_03_gene_screen_calls.tsv", sep="\t", index=False)

verdicts = []
for g, sub in calls.groupby("gene"):
    n_con = int(sub["concordant"].sum()); n_rev = int(sub["reversed"].sum())
    n_uns = int(sub["hit_sign_unknown"].sum())
    cm, li = strat_map.get(("HO-V-SC2-03", g), (None, None))
    same_p = bool(sub["same_pathogen"].iloc[0])
    if n_rev > 0 and n_con == 0 and n_uns == 0:
        v = "FAIL" if same_p else "XPATH_REVERSED"
    elif n_con > 0 or n_uns > 0:
        v = "PASS" if same_p else "PASS_XPATH"
    elif not same_p:
        v = "NOT_EVALUABLE_XPATH"
    elif cm is True:
        v = "FAIL"
    elif cm is False:
        v = "CONTEXT_BOUNDARY"
    else:
        v = "NOT_STRATIFIED"
    verdicts.append({"holdout_id": "HO-V-SC2-03", "gene": g, "n_screens": len(sub),
                     "n_concordant": n_con, "n_reversed": n_rev, "n_hit_unsigned": n_uns,
                     "same_pathogen": same_p, "context_matched": cm, "lab_independent": li,
                     "verdict": v})
verd = pd.DataFrame(verdicts)
fz = pd.read_csv(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", sep="\t",
                 usecols=["gene", "bio_tier", "intervention_direction"])
verd = verd.merge(fz, on="gene", how="left")
verd["is_hs"] = verd["gene"].isin(HS)
verd.to_csv(DATA / "stage13b_sc2_03_verdicts.tsv", sep="\t", index=False)

ev = verd[verd["verdict"].isin(["PASS", "FAIL", "CONTEXT_BOUNDARY"])]
def rate(sub):
    return {"n": len(sub), "pass": int((sub["verdict"] == "PASS").sum()),
            "rate": round(float((sub["verdict"] == "PASS").mean()), 4) if len(sub) else None}
summary = {
    "stage": "Stage13b_HO-V-SC2-03_supplementary", "date": "2026-09-04",
    "data_source": "Stahl et al. 2024 author gene-level results (boutroslab/Supp_Stahl_2024);"
                   " author hit calls as-is (Calu-3 FDR<0.10, A549 FDR<0.05); orientation +1 anchored on ACE2",
    "note": "补充裁决,不修改 Stage 13 冻结输出;分层用冻结 R14 文件(含 HO-V-SC2-03)",
    "n_candidates_matched": len(verd),
    "verdict_counts": verd["verdict"].value_counts().to_dict(),
    "same_pathogen_evaluable": rate(ev),
    "by_tier": {t: rate(g) for t, g in ev.groupby("bio_tier")},
    "by_context_matched": {str(k): rate(g) for k, g in ev.groupby("context_matched")},
    "hs_cluster": verd[verd["is_hs"]][["gene", "verdict", "n_concordant", "n_reversed",
                                       "context_matched"]].to_dict("records"),
    "positive_controls": verd[verd["gene"].isin(POS_CTRL)][["gene", "verdict", "n_concordant",
                                                            "n_reversed"]].to_dict("records"),
    "biotier2_verdicts": verd[verd["bio_tier"] == "Bio-Tier 2"][["gene", "verdict"]].to_dict("records"),
    "input_sha256": {"calu3_results": sha(DATA / "Calu3_results_crispr.txt"),
                     "a549_results": sha(DATA / "A549_results_crispr.txt"),
                     "candidate_set": sha(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv"),
                     "stratification_frozen": sha(ROOT / "16_revision_r14_diagnostics/r14_holdout_context_stratification.tsv")},
}
json.dump(summary, open(DATA / "stage13b_sc2_03_summary.json", "w"), indent=2, ensure_ascii=False)

write_checkpoint("stage13b_sc2_03",
                 str(ROOT / "scripts/48_stage13b_sc2_03_adjudication.py"),
                 [str(DATA / "Calu3_results_crispr.txt"), str(DATA / "A549_results_crispr.txt"),
                  str(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv"),
                  str(ROOT / "16_revision_r14_diagnostics/r14_holdout_context_stratification.tsv")],
                 None, str(ROOT / "00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md"),
                 [str(DATA / "stage13b_sc2_03_gene_screen_calls.tsv"),
                  str(DATA / "stage13b_sc2_03_verdicts.tsv"),
                  str(DATA / "stage13b_sc2_03_summary.json")],
                 ["stage13_holdout_adjudication"])

print(f"匹配候选 {len(verd)}/353;裁决分布: {summary['verdict_counts']}")
print(f"同病原可裁决: {summary['same_pathogen_evaluable']}")
print(f"按情境: {json.dumps(summary['by_context_matched'])}")
print(f"按层级: {json.dumps(summary['by_tier'])}")
print("HS 簇:", [(r['gene'], r['verdict']) for r in summary['hs_cluster']])
print("阳性对照:", [(r['gene'], r['verdict']) for r in summary['positive_controls']])
print("Bio-Tier 2:", [(r['gene'], r['verdict']) for r in summary['biotier2_verdicts']])
