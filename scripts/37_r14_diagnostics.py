#!/usr/bin/env python3
"""R14: 诊断与披露包(计划 v3 Workstream A6)。

组成:
  1. 召回率零假设:13/20 的超几何期望 + 二项 CI + 阴性对照假阳性率 + 循环性披露
  2. 按读出类型的秩次分布诊断(rank_pct 均匀性,KS 检验)
  3. LOSO 稳定性 × 独立实验室簇数联合报告
  4. 门盲法披露(HG1/HG2 评审时基因名可见)
  5. 留存上下文分层矩阵 + 反向假设预注册

输出: 16_revision_r14_diagnostics/
  r14_recall_null.json
  r14_readout_rank_diagnostics.tsv / .json
  r14_loso_cluster_report.tsv
  r14_gate_blinding_disclosure.md
  r14_holdout_context_stratification.tsv
  r14_reverse_hypothesis_preregistration.md
  r14_runinfo.json
"""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "16_revision_r14_diagnostics"
OUT.mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
t0 = time.time()

# ---------- 1. 召回率零假设 ----------
truth = pd.read_csv(ROOT / "07_cross_pathogen/known_hit_masking_check_v2.tsv", sep="\t")
calls = pd.read_csv(ROOT / "05_conditional_dependency/per_screen_calls_v2.tsv", sep="\t",
                    usecols=["gene", "call"])
cand = pd.read_csv(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv", sep="\t")
N = calls["gene"].nunique()                      # 被测基因宇宙
K = len(truth)                                   # truth set 大小 = 20
n = len(cand)                                    # 候选数 = 353
k = int((truth["status"] == "RECOVERED").sum())  # 13
expected = n * K / N
p_hyper = stats.hypergeom.sf(k - 1, N, K, n)
ci = stats.binomtest(k, K).proportion_ci(confidence_level=0.95, method="wilsoncc")

# 阴性对照:操作性非必需基因(相关背景 Chronos 全 > -0.5)
bl = pd.read_csv(ROOT / "04_depmap_baseline/baseline_chronos_by_background.tsv", sep="\t").set_index("gene")
REL = ["A-549", "Calu-3", "HuH-7", "Huh-7.5", "Huh-7.5.1", "HEp-2", "HeLa"]
rel = [c for c in REL if c in bl.columns]
noness = set(bl.index[(bl[rel] > -0.5).all(axis=1)])
cand_genes = set(cand["gene"])
fpr_noness = len(cand_genes & noness) / len(noness)
tpr = k / K
# 循环性:truth set 中 4 个 RSV 基因(EXT1/EXT2/SLC35B2/B3GAT3)来自同一屏文献
circ = truth[truth["pathogen"] == "RSV"]["gene"].tolist()
recall_null = {
    "recall": f"{k}/{K}", "binomial_ci95": [round(ci.low, 3), round(ci.high, 3)],
    "gene_universe_N": int(N), "candidate_n": n,
    "hypergeometric_expected": round(expected, 2), "hypergeometric_p": float(p_hyper),
    "negative_control": {"definition": "相关背景 Chronos 均 > -0.5 的操作性非必需基因",
                         "n_nonessential": len(noness), "n_in_candidates": len(cand_genes & noness),
                         "fpr": round(fpr_noness, 5), "tpr": round(tpr, 3),
                         "baseline_rate_random": round(n / N, 5)},
    "circularity_disclosure": f"truth set 中 RSV 基因 {circ} 来自与语料相同的屏文献(Sunshine2025 A-549),"
                              f"其'召回'非独立验证;独立召回应仅计 SC2/IAV 部分",
    "sc2_iav_only": {"k": int((truth[truth.pathogen != 'RSV']['status'] == 'RECOVERED').sum()),
                     "K": int((truth.pathogen != 'RSV').sum())},
}
json.dump(recall_null, open(OUT / "r14_recall_null.json", "w"), indent=2, ensure_ascii=False)
print("召回零假设:", json.dumps(recall_null, indent=1)[:600], flush=True)

# ---------- 2. 按读出类型的秩次分布 ----------
calls2 = pd.read_csv(ROOT / "05_conditional_dependency/per_screen_calls_v2.tsv", sep="\t",
                     usecols=["gene", "endpoint_family", "rank_pct", "call", "track"])
diag_rows = []
for (ef, tr), sub in calls2.groupby(["endpoint_family", "track"]):
    rp = sub["rank_pct"].dropna()
    rp = rp[(rp >= 0) & (rp <= 1)]
    if len(rp) < 100:
        continue
    ks = stats.kstest(rp, "uniform")
    sig = sub["call"].isin(["CDEP_PATHOGEN_EXPLOITED", "CDEP_HOST_PROTECTIVE"])
    diag_rows.append({"endpoint_family": ef, "track": tr, "n": len(sub),
                      "rank_pct_median": float(rp.median()),
                      "ks_stat_vs_uniform": float(ks.statistic), "ks_p": float(ks.pvalue),
                      "frac_called": float(sig.mean())})
diag = pd.DataFrame(diag_rows)
diag.to_csv(OUT / "r14_readout_rank_diagnostics.tsv", sep="\t", index=False)
json.dump({"interpretation": "rank_pct 相对均匀分布的 KS 偏差反映该读出类型的信号富集程度;"
                             "frac_called 跨读出类型差异大 → 置换/阈值需按读出分层(已在 Stage 5 v2 处理)",
           "by_readout": diag.to_dict("records")},
          open(OUT / "r14_readout_rank_diagnostics.json", "w"), indent=2, ensure_ascii=False)
print("读出诊断:", diag.to_string(index=False), flush=True)

# ---------- 3. LOSO × 簇数 ----------
loso = pd.read_csv(ROOT / "07_cross_pathogen/convergence_v2_loso.tsv", sep="\t")
rep = pd.read_csv(ROOT / "06_within_pathogen/within_pathogen_replication_v2.tsv", sep="\t")
nclu = rep.groupby("gene")["n_independence_clusters"].max().rename("n_independence_clusters_max")
labs = rep.groupby("gene")["supporting_clusters"].agg(
    lambda s: len(set(";".join(s.dropna()).split(";")))).rename("n_supporting_clusters_total")
lr = loso.merge(nclu, on="gene", how="left").merge(labs, on="gene", how="left")
lr.to_csv(OUT / "r14_loso_cluster_report.tsv", sep="\t", index=False)
print(f"LOSO 报告: {len(lr)} 行;stability=1.0 且簇数>=2: {int(((lr.stability==1.0)&(lr.n_independence_clusters_max>=2)).sum())}", flush=True)

# ---------- 4. 门盲法披露 ----------
(OUT / "r14_gate_blinding_disclosure.md").write_text("""# 门盲法披露(R14,2026-09-02)

## 事实
- HG1(语料冻结)与 HG2(参照臂决策)评审时,评审者可见基因名与屏身份。
  门审并非盲态:评审材料包含屏级 QC 与已知命中检查(如 ACE2/TMPRSS2)。
- HG3 评分细则在 Stage 8–11 结果生成**前**冻结(2026-08-31),评分脚本在结果
  最终化后才运行 —— 该环节为分析盲态(规则先于数据)。
- Stage 12 冻结后,候选身份对已解锁分析公开;holdout(12 屏)保持锁定,
  仅身份/角色可见(accession_only_no_results)。

## 残余风险
- HG1/HG2 阶段的非盲态可能引入"知名基因宽容"偏差(例如对已知命中的 QC 失败
  更倾向保留)。缓解:QC 规则为预注册阈值,逐屏适用,无基因级例外记录。
- 本披露不改变任何已冻结产物,仅记录评审过程的盲态边界。

## 结论
分析盲态(规则先于数据)在 HG3 评分环节成立;门审环节为非盲态,已披露。
""")

# ---------- 5. 留存上下文分层 + 反向假设预注册 ----------
ho = pd.read_csv(ROOT / "09_candidate_freeze/holdout_manifest.tsv", sep="\t")
# 每候选的发现语境(细胞系 + 簇)
ctx = rep.groupby("gene").agg(
    discovery_cell_lines=("cell_lines", lambda s: ";".join(sorted(set(";".join(s.dropna()).split(";"))))),
    discovery_pathogens=("pathogen", lambda s: ";".join(sorted(set(s)))),
    n_clusters=("n_independence_clusters", "max")).reset_index()
# 留存屏的细胞背景(从 manifest rationale 提取,人工标注)
ho_bg = {"HO-V-SC2-01": "Huh-7.5", "HO-V-SC2-02": "UM-UC-4;HCT 116",
         "HO-V-SC2-03": "A-549;Calu-3", "HO-V-IAV-01": "A-549;HeLa",
         "HO-V-IAV-02": "HEK293SF", "HO-V-RSV-01": "HEp-2",
         "HO-B-Mtb-01": "THP-1", "HO-B-Mtb-02": "murine-Hoxb8",
         "HO-CLIN-01": "clinical", "HO-SC-01": "single_cell"}
ho["cell_background"] = ho["holdout_id"].map(ho_bg)
path_map = {"SARS-CoV-2": "SARS-CoV-2", "Influenza A": "Influenza", "RSV": "RSV", "Mtb/BCG": "Mtb/BCG"}
rows = []
for _, h in ho.iterrows():
    hp = path_map.get(h["pathogen"], h["pathogen"])
    for _, c in ctx.iterrows():
        if hp not in c["discovery_pathogens"]:
            continue
        disc_lines = set(c["discovery_cell_lines"].split(";")) if c["discovery_cell_lines"] else set()
        ho_lines = set(h["cell_background"].split(";"))
        matched = bool(disc_lines & ho_lines)
        rows.append({"holdout_id": h["holdout_id"], "holdout_bg": h["cell_background"],
                     "gene": c["gene"], "discovery_cell_lines": c["discovery_cell_lines"],
                     "context_matched": matched,
                     "lab_independent": "Independent lab" in str(h["independence_from_discovery"])})
strat = pd.DataFrame(rows)
strat.to_csv(OUT / "r14_holdout_context_stratification.tsv", sep="\t", index=False)
print(f"留存分层矩阵: {len(strat)} 基因×留存对", flush=True)

hs_genes = ["EXT1", "EXT2", "EXTL3", "B4GALT7", "B3GAT3", "SLC35B2"]
(OUT / "r14_reverse_hypothesis_preregistration.md").write_text(f"""# 反向假设预注册(R14,2026-09-02;HG4 解锁前冻结)

## 背景
评审 §1:HS 基因簇({', '.join(hs_genes)})的跨病原体收敛可能由永生化细胞系
语境驱动,而非真实的广谱宿主依赖。RSV 在分化原代气道上皮中的已知受体为
CX3CR1/IGF1R 而非 HS。

## 预注册预测(解锁前冻结,不得事后修改)
1. **H1(语境驱动)**:HS 基因在**语境匹配**留存屏(发现背景 = 永生化系,
   如 HO-V-SC2-01 Huh-7.5)中 PASS,而在**语境不匹配**留存屏
   (HO-V-SC2-02 UM-UC-4/HCT116 等非肝非肺背景)中 FAIL。
   若观察到该模式 → 支持"细胞语境驱动收敛",HS 簇降级为语境特异发现。
2. **H2(真实广谱)**:HS 基因在匹配与不匹配留存屏中均 PASS
   → 反对语境伪迹假说,支持真实跨语境依赖。
3. **判定规则**:按 §34 的 PASS/CONTEXT_BOUNDARY/FAIL 词汇表;
   语境匹配/不匹配按 r14_holdout_context_stratification.tsv 的
   context_matched 列;统计为 2×2 Fisher 精确检验(匹配/不匹配 × PASS/非PASS),
   显著性阈值 α=0.05,单尾(方向:匹配 PASS 率 > 不匹配 PASS 率)。
4. **对照基因集**:ACE2/TMPRSS2(已知真实依赖,预期匹配与不匹配中
   按已知细胞趋向性分布)与 44 个 ORCS top5% 频繁命中者(预期普遍 PASS,
   作为频繁命中者对照)。

## 禁止
- 不得根据留存结果事后修改 context_matched 定义或基因集构成。
- 任何偏离本预注册的分析须标记为探索性。
""")

json.dump({"revision": "R14", "date": "2026-09-02", "runtime_min": round((time.time() - t0) / 60, 1),
           "input_sha256": {"known_hit_masking_check_v2": sha(ROOT / "07_cross_pathogen/known_hit_masking_check_v2.tsv"),
                            "holdout_manifest": sha(ROOT / "09_candidate_freeze/holdout_manifest.tsv"),
                            "per_screen_calls_v2": sha(ROOT / "05_conditional_dependency/per_screen_calls_v2.tsv")}},
          open(OUT / "r14_runinfo.json", "w"), indent=2, ensure_ascii=False)
print(f"\n总耗时 {(time.time()-t0)/60:.1f} min")
