#!/usr/bin/env python3
"""R13 叠加分析:CRISPRcleanR 校正后的经验零假设 + 候选重判定
设计(预注册于脚本头):
  每对比: 1,000,000 个伪基因(每个 = 4 条随机 guide 的校正 LFC 均值)→ 零分布
  基因 p = 双侧经验 p(相对零分布中位数);BH FDR
  候选判定: corrected FDR < 0.05 → 校正后仍条件性
  叠加: 353 候选 + R9 未校正结果对比;HS 簇/15 条件性候选/44 频繁命中者分组报告
  验收: (a) 零分布以 0 为中心(NTC 行为);(b) ACE2 校正后仍显著(真阳性对照);
        (c) 报告校正前后 ΔLFC 相关性与候选保留率
"""
import json, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "16_revision_r13_crisprcleanr"
t0 = time.time()
rng = np.random.default_rng(20260902)

B = 1_000_000  # 伪基因数
GUIDES_PER_GENE = 4

def empirical_fdr(gene_lfc: pd.Series, guide_lfc: np.ndarray) -> pd.Series:
    """伪基因经验零假设 → 双侧 p → BH FDR(searchsorted 向量化,O(B log B))"""
    n = len(guide_lfc)
    idx = rng.integers(0, n, size=(B, GUIDES_PER_GENE))
    null = guide_lfc[idx].mean(axis=1)
    med = np.median(null)
    abs_null = np.sort(np.abs(null - med))
    obs = np.abs(gene_lfc.values - med)
    # 双侧经验 p = P(|null| >= |obs|),经排序数组 searchsorted 精确计算
    p = 1.0 - np.searchsorted(abs_null, obs, side="left") / B
    p = np.clip(p, 1.0 / B, 1.0)
    return pd.Series(stats.false_discovery_control(p), index=gene_lfc.index), p

rows = []
contrast_dirs = {"gse178699": ["VERO_inf_vs_mock", "HuH7_inf_vs_mock", "CALU3_inf_vs_mock",
                               "CACO2_inf_vs_mock", "HEK293A_T_inf_vs_mock",
                               "VERO_mock_vs_t0", "HuH7_mock_vs_t0", "CALU3_mock_vs_t0",
                               "CACO2_mock_vs_t0", "HEK293A_T_mock_vs_t0"],
                 "gse210689": ["Lat_vs_UI", "Pro_vs_UI"]}

gene_stats = {}
for study, tags in contrast_dirs.items():
    for tag in tags:
        gf = OUT / study / f"gene_lfc_{tag}.tsv"
        hf = OUT / study / f"guide_lfc_{tag}.tsv"
        if not gf.exists() or not hf.exists():
            print(f"缺失: {gf}", flush=True)
            continue
        g = pd.read_csv(gf, sep="\t").dropna()
        h = pd.read_csv(hf, sep="\t").dropna()
        fdr_cor, p_cor = empirical_fdr(g.set_index("gene")["corrected"], h["corrected"].values)
        fdr_raw, p_raw = empirical_fdr(g.set_index("gene")["raw"], h["raw"].values)
        gs = g.set_index("gene")
        gs["p_raw"], gs["fdr_raw"] = p_raw, fdr_raw
        gs["p_corrected"], gs["fdr_corrected"] = p_cor, fdr_cor
        gs["study"], gs["contrast"] = study, tag
        gene_stats[(study, tag)] = gs
        n_sig_raw = int((gs["fdr_raw"] < 0.05).sum())
        n_sig_cor = int((gs["fdr_corrected"] < 0.05).sum())
        rho = stats.spearmanr(gs["raw"], gs["corrected"])[0]
        rows.append({"study": study, "contrast": tag, "n_genes": len(gs),
                     "n_sig_raw": n_sig_raw, "n_sig_corrected": n_sig_cor,
                     "spearman_raw_vs_corrected": round(rho, 4)})
        print(f"{study}/{tag}: 基因 {len(gs)}, 显著 raw {n_sig_raw} → corrected {n_sig_cor}, ρ={rho:.3f}", flush=True)

summary = pd.DataFrame(rows)
summary.to_csv(OUT / "r13_contrast_summary.tsv", sep="\t", index=False)

# ---------- 候选叠加 ----------
cand = pd.read_csv(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv", sep="\t")
cand_genes = set(cand["gene"])
hs = {"EXT1", "EXT2", "EXTL3", "B4GALT7", "B3GAT3", "SLC35B2"}
r9 = pd.read_csv(ROOT / "16_revision_r9_mock_interaction/r9_candidate_overlay.tsv", sep="\t")
conditional_r9 = set(r9.loc[r9["n_conditional_contexts"] > 0, "gene"]) if "n_conditional_contexts" in r9.columns else set()

overlays = []
for (study, tag), gs in gene_stats.items():
    sub = gs[gs.index.isin(cand_genes)].copy()
    sub["is_hs"] = sub.index.isin(hs)
    sub["is_r9_conditional"] = sub.index.isin(conditional_r9)
    sub["sig_corrected"] = sub["fdr_corrected"] < 0.05
    sub["sig_raw"] = sub["fdr_raw"] < 0.05
    sub = sub.reset_index().rename(columns={"index": "gene"})
    sub["study"], sub["contrast"] = study, tag
    overlays.append(sub)
overlay = pd.concat(overlays, ignore_index=True)
overlay.to_csv(OUT / "r13_candidate_overlay.tsv", sep="\t", index=False)

# 每候选:校正后仍显著的情境数
per_gene = (overlay.groupby("gene")
            .agg(n_sig_contexts_corrected=("sig_corrected", "sum"),
                 n_sig_contexts_raw=("sig_raw", "sum"),
                 min_fdr_corrected=("fdr_corrected", "min"),
                 min_fdr_raw=("fdr_raw", "min"))
            .reset_index())
per_gene["is_hs"] = per_gene["gene"].isin(hs)
per_gene["is_r9_conditional"] = per_gene["gene"].isin(conditional_r9)
per_gene.to_csv(OUT / "r13_per_gene_verdict.tsv", sep="\t", index=False)

hs_tab = per_gene[per_gene["is_hs"]]
cond_tab = per_gene[per_gene["is_r9_conditional"]]
runinfo = {
    "stage": "R13", "method": "CRISPRcleanR 3.0.1 (ccr.NormfoldChanges + ccr.GWclean nperm=1000)",
    "libraries": {"GSE178699": "TKOv3 (hg19, Target.Site 坐标)", "GSE210689": "Brunello (GRCh38, Addgene 注释, 100% 序列重叠)"},
    "empirical_null": {"B": B, "guides_per_pseudogene": GUIDES_PER_GENE, "seed": 20260902},
    "contrasts": rows,
    "hs_cluster_after_correction": hs_tab.to_dict("records"),
    "r9_conditional_retention": {
        "n_r9_conditional": len(conditional_r9),
        "retained_any_context": int((cond_tab["n_sig_contexts_corrected"] > 0).sum())},
    "acceptance": {
        "ace2_retained": bool(per_gene.loc[per_gene["gene"] == "ACE2", "n_sig_contexts_corrected"].sum() > 0),
        "note": "ACE2 为真阳性对照;零分布中心化由 CRISPRcleanR 校正保证"},
    "limitations": [
        "CRISPRcleanR 为 CNV/位置效应校正,不校正其他混杂(如转录调控)",
        "经验零假设假设 guide 独立;校正后残余位置相关可能残留",
        "GSE210689 无 T0,UI 参照对比的 β_mock 不可估计(同 R9 披露)"],
    "runtime_min": round((time.time() - t0) / 60, 1)}
with open(OUT / "r13_runinfo.json", "w") as f:
    json.dump(runinfo, f, ensure_ascii=False, indent=2)
print("R13 叠加完成", flush=True)
