#!/usr/bin/env python3
"""专家意见 #1 + #3:冻结表派生标注(2026-09-04)。

治理约束:BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv 为 SHA256 冻结件,不修改;
本脚本产出派生标注表(全列 + 3 个新列),并在 runinfo 中记录冻结件哈希以证明未动。

新列定义(规则显式化,可复现):
  benchmark_class:
    benchmark_positive_control = 基因出现在 R14 预定义已知命中真值集
      (07_cross_pathogen/known_hit_masking_check_v2.tsv,20 病原体-基因对)
      或 Stage 13 阳性对照(ACE2, TMPRSS2;scripts/40 POS_CTRL)
    discovery_candidate = 其余
  benchmark_source: 来源明细(R14_known_hit:<pathogen>;stage13_pos_ctrl)
  screen_effect_direction(筛选表型方向,来自预注册候选集 directions 列):
    EXP -> LOF_protects_host      (敲除降低感染 = 宿主依赖/促病毒因子)
    PRO -> LOF_sensitizes_host    (敲除增加感染 = 限制因子)
  therapeutic_hypothesis(治疗假说,与筛选方向解耦):
    PRESERVE_ACTIVATE -> preserve_activate_not_direct_inhibition_target
    INHIBIT 且 min_chronos < -2  -> inhibition_candidate_window_unresolved_baseline_essentiality
    INHIBIT 且 -2<=min_chronos<=-1 -> inhibition_candidate_moderate_essentiality_caution
    INHIBIT 且 min_chronos 缺失   -> inhibition_candidate_window_unknown_no_depmap
    其余 INHIBIT                  -> inhibition_candidate

输出:
  09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2_annotated.tsv
  09_candidate_freeze/freeze_v2_annotation_runinfo.json
"""
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/workspace/HPICDA")
sys.path.insert(0, str(ROOT / "runner"))
from checkpoint import write_checkpoint

FREEZE = ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv"
FREEZE_SHA256_EXPECTED = "1c00f7cc694b3dff15dd671f1288e347726e0b69f0af1657892e71bd02d8d4dc"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

actual = sha(FREEZE)
assert actual == FREEZE_SHA256_EXPECTED, f"冻结 v2 哈希不符:{actual}"
print(f"冻结 v2 哈希校验通过(只读,未修改):{actual[:16]}…")

fz = pd.read_csv(FREEZE, sep="\t")
cand = pd.read_csv(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv", sep="\t",
                   usecols=["gene", "directions"])
truth = pd.read_csv(ROOT / "07_cross_pathogen/known_hit_masking_check_v2.tsv", sep="\t")
POS_CTRL = ["ACE2", "TMPRSS2"]  # scripts/40 冻结定义

truth_map = truth.groupby("gene")["pathogen"].agg(lambda s: ";".join(sorted(set(s)))).to_dict()

def benchmark(row):
    src = []
    if row["gene"] in truth_map:
        src.append(f"R14_known_hit:{truth_map[row['gene']]}")
    if row["gene"] in POS_CTRL:
        src.append("stage13_pos_ctrl")
    return pd.Series({
        "benchmark_class": "benchmark_positive_control" if src else "discovery_candidate",
        "benchmark_source": ";".join(src)})

fz = fz.merge(cand.rename(columns={"directions": "screen_direction_raw"}), on="gene", how="left")
fz = pd.concat([fz, fz.apply(benchmark, axis=1)], axis=1)

DIR_MAP = {"EXP": "LOF_protects_host", "PRO": "LOF_sensitizes_host"}
fz["screen_effect_direction"] = fz["screen_direction_raw"].map(DIR_MAP)

def hypothesis(row):
    if row["intervention_direction"] == "PRESERVE_ACTIVATE":
        return "preserve_activate_not_direct_inhibition_target"
    mc = row["min_baseline_chronos"]
    if pd.isna(mc):
        return "inhibition_candidate_window_unknown_no_depmap"
    if mc < -2:
        return "inhibition_candidate_window_unresolved_baseline_essentiality"
    if mc <= -1:
        return "inhibition_candidate_moderate_essentiality_caution"
    return "inhibition_candidate"

fz["therapeutic_hypothesis"] = fz.apply(hypothesis, axis=1)

out_tsv = ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2_annotated.tsv"
fz.to_csv(out_tsv, sep="\t", index=False)

bt2 = fz[fz["bio_tier"] == "Bio-Tier 2"]
runinfo = {
    "analysis": "freeze v2 derived annotation (expert items #1 + #3)",
    "date": "2026-09-04",
    "governance": ("冻结件未修改;SHA256 校验见 freeze_v2_sha256_verified。"
                   "本表为派生标注,不进入任何已冻结评分。"),
    "freeze_v2_sha256_verified": actual,
    "n_candidates": len(fz),
    "benchmark_class_counts": fz["benchmark_class"].value_counts().to_dict(),
    "biotier2_benchmark_split": bt2["benchmark_class"].value_counts().to_dict(),
    "biotier2_benchmark_genes": sorted(bt2.loc[bt2["benchmark_class"] == "benchmark_positive_control", "gene"]),
    "biotier2_discovery_genes": sorted(bt2.loc[bt2["benchmark_class"] == "discovery_candidate", "gene"]),
    "screen_effect_direction_counts": fz["screen_effect_direction"].value_counts().to_dict(),
    "therapeutic_hypothesis_counts": fz["therapeutic_hypothesis"].value_counts().to_dict(),
    "input_sha256": {"freeze_v2": actual,
                     "candidate_set": sha(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv"),
                     "r14_truth": sha(ROOT / "07_cross_pathogen/known_hit_masking_check_v2.tsv")},
}
json.dump(runinfo, open(ROOT / "09_candidate_freeze/freeze_v2_annotation_runinfo.json", "w"),
          indent=2, ensure_ascii=False)

write_checkpoint("freeze_v2_annotation",
                 str(ROOT / "scripts/45_freeze_annotation_addon.py"),
                 [str(FREEZE), str(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv"),
                  str(ROOT / "07_cross_pathogen/known_hit_masking_check_v2.tsv")],
                 None,
                 str(ROOT / "00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md"),
                 [str(out_tsv), str(ROOT / "09_candidate_freeze/freeze_v2_annotation_runinfo.json")],
                 ["stage12_candidate_freeze_v2"])

print(f"\nbenchmark_class 全体:{runinfo['benchmark_class_counts']}")
print(f"Bio-Tier 2 拆分:{runinfo['biotier2_benchmark_split']}")
print(f"  benchmark:{runinfo['biotier2_benchmark_genes']}")
print(f"  discovery:{runinfo['biotier2_discovery_genes']}")
print(f"\nscreen_effect_direction:{runinfo['screen_effect_direction_counts']}")
print(f"therapeutic_hypothesis:{runinfo['therapeutic_hypothesis_counts']}")
