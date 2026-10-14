#!/usr/bin/env python3
"""Stage 15 复合分权重敏感性(2026-09-04,收尾项)。

等权复合分已披露为局限(milestone v3 局限 #6)。本分析检验排序对权重选择的稳健性:
  方案 A:留一法 —— 每次去掉一个维度,其余 4 维等权(5 套)
  方案 B:Dirichlet(α=1)随机权重 2,000 套 —— 排名分布
指标:与等权基线排名的 Spearman ρ;top-14(=前沿 1 大小)重叠;逐基因排名区间。
Pareto 前沿本身不依赖权重(非支配排序),本敏感性仅针对可读性复合分。

输出:
  22_stage15_pareto/stage15_weight_sensitivity.tsv
  22_stage15_pareto/stage15_weight_sensitivity.json
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "22_stage15_pareto"
sys.path.insert(0, str(ROOT / "runner"))
from checkpoint import write_checkpoint

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
DIMS = ["x_evidence", "x_holdout", "x_specificity", "x_tolerability", "x_druggability"]

d = pd.read_csv(OUT / "stage15_pareto_priority.tsv", sep="\t")
Z = (d[DIMS] - d[DIMS].min()) / (d[DIMS].max() - d[DIMS].min() + 1e-12)
base_score = Z.mean(axis=1)
base_rank = base_score.rank(ascending=False, method="min")
top14_base = set(d.loc[base_rank <= 14, "gene"])

rng = np.random.default_rng(20260904)
schemes = {f"leave_out_{dim}": np.array([0.0 if x == dim else 0.25 for x in DIMS]) for dim in DIMS}
W_rand = rng.dirichlet(np.ones(5), size=2000)

rank_mat = np.zeros((len(d), len(schemes) + len(W_rand)))
cols = list(schemes) + [f"dirichlet_{i}" for i in range(len(W_rand))]
for j, name in enumerate(list(schemes)):
    s = Z.to_numpy() @ schemes[name]
    rank_mat[:, j] = pd.Series(s).rank(ascending=False, method="min").to_numpy()
for j in range(len(W_rand)):
    s = Z.to_numpy() @ W_rand[j]
    rank_mat[:, len(schemes) + j] = pd.Series(s).rank(ascending=False, method="min").to_numpy()

# 方案级汇总
scheme_rows = []
for j, name in enumerate(cols):
    r = rank_mat[:, j]
    rho = spearmanr(base_rank, r).statistic
    top14 = set(d.loc[r <= 14, "gene"])
    scheme_rows.append({"scheme": name, "spearman_vs_equal": round(float(rho), 4),
                        "top14_overlap": len(top14 & top14_base)})
sch = pd.DataFrame(scheme_rows)

# 逐基因汇总
per_gene = pd.DataFrame({
    "gene": d["gene"], "bio_tier": d["bio_tier"], "pareto_front": d["pareto_front"],
    "composite_equal_weight": base_score.round(4),
    "rank_equal_weight": base_rank.astype(int),
    "rank_min_all_schemes": rank_mat.min(axis=1).astype(int),
    "rank_max_all_schemes": rank_mat.max(axis=1).astype(int),
    "rank_mean_dirichlet": rank_mat[:, len(schemes):].mean(axis=1).round(1),
    "rank_p05_dirichlet": np.percentile(rank_mat[:, len(schemes):], 5, axis=1).round(1),
    "rank_p95_dirichlet": np.percentile(rank_mat[:, len(schemes):], 95, axis=1).round(1),
}).sort_values("rank_equal_weight")
per_gene.to_csv(OUT / "stage15_weight_sensitivity.tsv", sep="\t", index=False)

loo = sch[sch["scheme"].str.startswith("leave_out")]
diri = sch[sch["scheme"].str.startswith("dirichlet")]
key_genes = ["ACE2", "B4GALT7", "FER", "SLC35B2", "EXT1", "ACSL3", "KCNA3"]
summary = {
    "analysis": "Stage 15 composite weight sensitivity", "date": "2026-09-04",
    "note": "Pareto 前沿为非支配排序,不依赖权重;敏感性仅针对等权复合分。",
    "n_genes": len(d), "n_schemes": len(cols),
    "leave_one_out": {"min_spearman": float(loo["spearman_vs_equal"].min()),
                      "min_top14_overlap": int(loo["top14_overlap"].min()),
                      "detail": loo.to_dict("records")},
    "dirichlet_2000": {"min_spearman": float(diri["spearman_vs_equal"].min()),
                       "p05_spearman": float(diri["spearman_vs_equal"].quantile(0.05)),
                       "min_top14_overlap": int(diri["top14_overlap"].min()),
                       "p05_top14_overlap": float(diri["top14_overlap"].quantile(0.05))},
    "key_gene_rank_ranges": per_gene[per_gene["gene"].isin(key_genes)].to_dict("records"),
    "input_sha256": {"stage15_pareto_priority": sha(OUT / "stage15_pareto_priority.tsv")},
}
json.dump(summary, open(OUT / "stage15_weight_sensitivity.json", "w"), indent=2, ensure_ascii=False)

write_checkpoint("stage15_weight_sensitivity",
                 str(ROOT / "scripts/46_stage15_weight_sensitivity.py"),
                 [str(OUT / "stage15_pareto_priority.tsv")], None,
                 str(ROOT / "00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md"),
                 [str(OUT / "stage15_weight_sensitivity.tsv"),
                  str(OUT / "stage15_weight_sensitivity.json")],
                 ["stage15_pareto"])

print(f"留一法:min Spearman={summary['leave_one_out']['min_spearman']}, "
      f"min top14 重叠={summary['leave_one_out']['min_top14_overlap']}/14")
print(f"Dirichlet×2000:min Spearman={summary['dirichlet_2000']['min_spearman']:.4f}, "
      f"P5={summary['dirichlet_2000']['p05_spearman']:.4f}, "
      f"min top14 重叠={summary['dirichlet_2000']['min_top14_overlap']}/14")
print("\n关键基因排名区间:")
print(per_gene[per_gene["gene"].isin(key_genes)].to_string(index=False))
