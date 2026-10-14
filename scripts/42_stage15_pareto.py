#!/usr/bin/env python3
"""Stage 15 — Pareto prioritization with tolerability dimension (Plan v3 C3).
Dimensions (all maximized) over Bio-Tier 2+3 candidates:
  evidence      = freeze v2 total_score
  holdout       = Stage 13 same-pathogen n_pass
  specificity   = 1 - ORCS hit frequency (non-promiscuity)
  tolerability  = 1 - human_intolerant (R12; binary, permissive definition disclosed)
  druggability  = 0 none / 1 tractable / 2 clinical-stage / 3 launched (Stage 14)
Non-dominated sorting (Kung et al. simple O(n^2) fine at n=267).
"""
import json, os, datetime
import pandas as pd
import numpy as np

BASE = "/workspace/HPICDA"
OUT = f"{BASE}/22_stage15_pareto"
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(f"{BASE}/21_stage14_druggability/stage14_druggability_axis.tsv", sep="\t")
d = df[df["bio_tier"].astype(str).str.contains("Bio-Tier 2|Bio-Tier 3", na=False)].copy()

d["x_evidence"] = pd.to_numeric(d["total_score"], errors="coerce").fillna(0)
d["x_holdout"] = pd.to_numeric(d["n_pass"], errors="coerce").fillna(0)
d["x_specificity"] = 1 - pd.to_numeric(d["orcs_hit_freq_all"], errors="coerce").fillna(0)
d["x_tolerability"] = 1 - pd.to_numeric(d.get("human_intolerant", 0), errors="coerce").fillna(0)
d["x_druggability"] = np.select(
    [d["n_launched"].fillna(0) > 0, d["max_phase_rank"].fillna(0) >= 1,
     (d["sm_tractable"].fillna("").str.len() > 0) | (d["ab_tractable"].fillna("").str.len() > 0)],
    [3, 2, 1], default=0)
DIMS = ["x_evidence", "x_holdout", "x_specificity", "x_tolerability", "x_druggability"]

X = d[DIMS].to_numpy(float)
n = len(d)
dominated = np.zeros(n, bool)
fronts = []
remaining = np.arange(n)
Xw = X.copy()
# non-dominated sort
idx_all = np.arange(n)
cur = idx_all.copy()
while len(cur):
    keep = []
    for i in cur:
        xi = X[i]
        dom = False
        for j in cur:
            if i == j:
                continue
            xj = X[j]
            if np.all(xj >= xi) and np.any(xj > xi):
                dom = True
                break
        if not dom:
            keep.append(i)
    fronts.append(keep)
    cur = np.array([i for i in cur if i not in keep])
d["pareto_front"] = 0
for fi, fr in enumerate(fronts, 1):
    d.loc[d.index[fr], "pareto_front"] = fi

# composite (equal weights, min-max normalized) as readable secondary ranking
Z = (d[DIMS] - d[DIMS].min()) / (d[DIMS].max() - d[DIMS].min() + 1e-12)
d["composite_equal_weight"] = Z.mean(axis=1).round(4)
d = d.sort_values(["pareto_front", "composite_equal_weight"], ascending=[True, False])
cols = ["gene", "bio_tier", "intervention_direction", "intervention_question", "pareto_front",
        "composite_equal_weight", "x_evidence", "x_holdout", "x_specificity", "x_tolerability",
        "x_druggability", "druggability_class", "launched_drugs", "loeuf", "mendelian_disease",
        "orcs_promiscuous_top5pct", "n_evaluable", "n_pass", "n_fail", "n_cb", "n_xpath_pass"]
cols = [c for c in cols if c in d.columns]
d[cols].to_csv(f"{OUT}/stage15_pareto_priority.tsv", sep="\t", index=False)

front_sizes = {f"front_{i+1}": len(fr) for i, fr in enumerate(fronts)}
f1 = d[d["pareto_front"] == 1]
summ = {"stage": "Stage15_pareto", "date": "2026-09-03", "n_ranked": n, "dimensions": DIMS,
        "front_sizes": front_sizes,
        "front1_genes": f1[["gene", "bio_tier", "composite_equal_weight"]].to_dict("records"),
        "front1_tier2": f1[f1["bio_tier"].astype(str).str.contains("Bio-Tier 2")]["gene"].tolist(),
        "notes": ["tolerability 维度为二值(人类不可耐受=0),R12 定义为宽松口径,仅作分诊不作否决",
                  "druggability 0-3 分级;launched 注释含间接/宽松映射(如 ACE2-lisinopril),需人工复核",
                  "PRO/PRESERVE 方向候选不进入干预性开发(评审 §7),Pareto 排序仅作证据展示"]}
with open(f"{OUT}/stage15_summary.json", "w") as f:
    json.dump(summ, f, indent=2, default=str)
with open(f"{OUT}/stage15_runinfo.json", "w") as f:
    json.dump({"script": os.path.abspath(__file__), "timestamp": datetime.datetime.now().isoformat(),
               "inputs": ["21_stage14_druggability/stage14_druggability_axis.tsv"]}, f, indent=2)
print(json.dumps({"n": n, "fronts": front_sizes,
                  "front1": [r["gene"] for r in summ["front1_genes"]][:40],
                  "front1_tier2": summ["front1_tier2"]}, indent=2))
