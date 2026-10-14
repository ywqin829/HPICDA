#!/usr/bin/env python3
"""专家意见 #8:Stage 13 top-K 富集(2026-09-04)。

问题:按冻结 v2 排序的 top-K 候选,在独立留存屏中的 PASS 率是否显著高于其余候选?
估计量:同病原体可裁决基因×留存对(verdict ∈ PASS/FAIL/CONTEXT_BOUNDARY),
PASS vs 非 PASS(FAIL+CB),top-K vs 其余,Fisher 精确检验(单侧,富集方向)。
同时报告基因级指标:top-K 中至少 1 个 PASS 的基因比例(限于有可裁决对的基因)。
排序键与评分脚本一致:total_score 降序,screen_evidence_pts 降序破冰。
不改任何冻结文件;作者 HIT 原样,不涉及调参。

输出:
  20_stage13_holdout_adjudication/stage13_topk_enrichment.tsv
  20_stage13_holdout_adjudication/stage13_topk_enrichment.json
"""
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "20_stage13_holdout_adjudication"
sys.path.insert(0, str(ROOT / "runner"))
from checkpoint import write_checkpoint

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

fz = pd.read_csv(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", sep="\t",
                 usecols=["gene", "total_score", "screen_evidence_pts", "bio_tier"])
fz = fz.sort_values(["total_score", "screen_evidence_pts"], ascending=False).reset_index(drop=True)
fz["rank"] = fz.index + 1

verd = pd.read_csv(OUT / "stage13_gene_holdout_verdicts.tsv", sep="\t")
ev = verd[verd["verdict"].isin(["PASS", "FAIL", "CONTEXT_BOUNDARY"])].copy()
ev["is_pass"] = ev["verdict"] == "PASS"
ev = ev.merge(fz[["gene", "rank", "bio_tier"]], on="gene", how="left")

rows = []
for K in [10, 20, 50, 100]:
    top = set(fz.loc[fz["rank"] <= K, "gene"])
    e = ev.copy()
    e["in_topk"] = e["gene"].isin(top)
    a = int(((e["in_topk"]) & (e["is_pass"])).sum())
    b = int(((e["in_topk"]) & (~e["is_pass"])).sum())
    c = int(((~e["in_topk"]) & (e["is_pass"])).sum())
    d = int(((~e["in_topk"]) & (~e["is_pass"])).sum())
    OR, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    # 基因级:有可裁决对的 top-K 基因中,至少 1 个 PASS 的比例
    g = e.groupby("gene")["is_pass"].agg(["sum", "count"])
    gk = g[g.index.isin(top)]
    gr = g[~g.index.isin(top)]
    rows.append({
        "K": K, "pairs_topk": a + b, "pass_topk": a,
        "pass_rate_topk": round(a / (a + b), 4) if a + b else None,
        "pairs_rest": c + d, "pass_rest": c,
        "pass_rate_rest": round(c / (c + d), 4) if c + d else None,
        "odds_ratio": round(OR, 3) if OR == OR else None,
        "fisher_p_one_sided": f"{p:.3e}",
        "genes_topk_evaluable": len(gk),
        "genes_topk_with_pass": int((gk["sum"] > 0).sum()),
        "gene_recall_topk": round(float((gk["sum"] > 0).mean()), 4) if len(gk) else None,
        "genes_rest_evaluable": len(gr),
        "gene_recall_rest": round(float((gr["sum"] > 0).mean()), 4) if len(gr) else None,
    })

res = pd.DataFrame(rows)
res.to_csv(OUT / "stage13_topk_enrichment.tsv", sep="\t", index=False)

summary = {
    "analysis": "Stage 13 top-K enrichment (expert item #8)",
    "date": "2026-09-04",
    "estimand": ("同病原体可裁决对(PASS/FAIL/CONTEXT_BOUNDARY)中 PASS 富集;"
                 "top-K 按冻结 v2 total_score 排序;Fisher 单侧(富集)"),
    "n_evaluable_pairs": len(ev),
    "n_pass_total": int(ev["is_pass"].sum()),
    "n_genes_evaluable": int(ev["gene"].nunique()),
    "ranking": "freeze v2 total_score desc, screen_evidence_pts desc tie-break",
    "results": rows,
    "input_sha256": {"freeze_v2": sha(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv"),
                     "stage13_verdicts": sha(OUT / "stage13_gene_holdout_verdicts.tsv")},
}
json.dump(summary, open(OUT / "stage13_topk_enrichment.json", "w"), indent=2, ensure_ascii=False)

write_checkpoint("stage13_topk_enrichment",
                 str(ROOT / "scripts/44_stage13_topk_enrichment.py"),
                 [str(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv"),
                  str(OUT / "stage13_gene_holdout_verdicts.tsv")],
                 None,
                 str(ROOT / "00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md"),
                 [str(OUT / "stage13_topk_enrichment.tsv"),
                  str(OUT / "stage13_topk_enrichment.json")],
                 ["stage13_holdout_adjudication"])

print(res.to_string(index=False))
