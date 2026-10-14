"""HPICDA script 52 — 留存 top-K 校准主结果(2026-09-05,投稿封口 A3).

问题:HPICDA 冻结排序的顶端是否在真正独立留存中富集可复现宿主依赖?这是主结果而非内部 QC。
设计:
  - 排序 = freeze v2 total_score 降序,screen_evidence_pts 降序打破并列(与 Stage 13 已报告一致;
    已验证 freeze v2 total_score 与 HG3 total_score 完全一致,为留存解锁前冻结,无泄漏)。
  - 结局 = 配对级(gene×holdout 同病原可裁决对,verdict ∈ {PASS,FAIL,CONTEXT_BOUNDARY})PASS 富集;
    同时给基因级 recall(top-K 基因中 ≥1 PASS 的比例)。
  - 两套分析集:Stage-13-only(与 2026-09-04 已报告数字连续)与 Stage-13+13b(完整主结果)。
  - 统计:Fisher 单侧(富集)p;OR 95% CI 用 log-OR(Woolf)法,零格加 0.5 连续性校正。
输出:holdout_topk_calibration.tsv(长表,figure-ready)+ holdout_topk_calibration.json。
"""
import json
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from pathlib import Path

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "24_closing_analyses"
OUT.mkdir(exist_ok=True)
K_LIST = [10, 20, 50, 100]
EVAL = {"PASS", "FAIL", "CONTEXT_BOUNDARY"}

# ---- 冻结排序(全部 353 候选) ----
fr = pd.read_csv(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2_annotated.tsv",
                 sep="\t", usecols=["gene", "total_score", "screen_evidence_pts"])
rank = (fr.sort_values(["total_score", "screen_evidence_pts"], ascending=False)
        .reset_index(drop=True))
rank["rank"] = np.arange(1, len(rank) + 1)

# ---- 留存裁决(Stage 13 + 13b) ----
v13 = pd.read_csv(ROOT / "20_stage13_holdout_adjudication/stage13_gene_holdout_verdicts.tsv", sep="\t")
v13b = pd.read_csv(ROOT / "23_holdout_sc2_03_stahl/stage13b_sc2_03_verdicts.tsv", sep="\t")
v13["analysis_set"] = "stage13_only"
both = pd.concat([v13, v13b], ignore_index=True)


def eval_pairs(df):
    """同病原可裁决对,标注是否 PASS。"""
    e = df[df["same_pathogen"] & df["verdict"].isin(EVAL)].copy()
    e["is_pass"] = (e["verdict"] == "PASS").astype(int)
    return e[["holdout_id", "gene", "is_pass"]]


def logor_ci(a, b, c, d):
    """a=pass_topk b=fail_topk c=pass_rest d=fail_rest;log-OR(Woolf)95% CI,零格加 0.5。"""
    aa, bb, cc, dd = a, b, c, d
    if min(aa, bb, cc, dd) == 0:
        aa, bb, cc, dd = aa + 0.5, bb + 0.5, cc + 0.5, dd + 0.5
    lo = np.log((aa * dd) / (bb * cc))
    se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return np.exp(lo - 1.96 * se), np.exp(lo + 1.96 * se)


def topk_table(pairs, set_name):
    pairs = pairs.merge(rank[["gene", "rank"]], on="gene", how="left")
    rows = []
    for K in K_LIST:
        topk_genes = set(rank.head(K)["gene"])
        tp = pairs[pairs["gene"].isin(topk_genes)]
        rp = pairs[~pairs["gene"].isin(topk_genes)]
        a = int(tp["is_pass"].sum()); b = int(len(tp) - a)
        c = int(rp["is_pass"].sum()); d = int(len(rp) - c)
        OR, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        lo, hi = logor_ci(a, b, c, d)
        # 基因级 recall
        g_top = tp.groupby("gene")["is_pass"].max()
        g_rest = rp.groupby("gene")["is_pass"].max()
        rows.append({
            "analysis_set": set_name, "K": K,
            "pairs_topk": len(tp), "pass_topk": a,
            "pass_rate_topk": round(a / len(tp), 4) if len(tp) else np.nan,
            "pairs_rest": len(rp), "pass_rest": c,
            "pass_rate_rest": round(c / len(rp), 4) if len(rp) else np.nan,
            "odds_ratio": round(OR, 3), "OR_CI_lo": round(lo, 3), "OR_CI_hi": round(hi, 3),
            "fisher_p_one_sided": f"{p:.2e}",
            "genes_topk_evaluable": int(len(g_top)),
            "genes_topk_with_pass": int(g_top.sum()),
            "gene_recall_topk": round(g_top.mean(), 4) if len(g_top) else np.nan,
            "genes_rest_evaluable": int(len(g_rest)),
            "gene_recall_rest": round(g_rest.mean(), 4) if len(g_rest) else np.nan,
        })
    return rows


rows = (topk_table(eval_pairs(v13), "stage13_only")
        + topk_table(eval_pairs(v13b), "sc2_03_only")
        + topk_table(eval_pairs(both), "stage13_plus_13b"))
res = pd.DataFrame(rows)
res.to_csv(OUT / "holdout_topk_calibration.tsv", sep="\t", index=False)

json.dump({"script": "52_holdout_topk_calibration.py", "date": "2026-09-05",
           "estimand": "同病原可裁决对(gene×holdout)中 PASS 富集;排序为 freeze v2 total_score(=HG3,留存前冻结)",
           "ranking": "freeze v2 total_score desc, screen_evidence_pts desc tie-break",
           "n_evaluable_pairs_stage13": int(len(eval_pairs(v13))),
           "n_evaluable_pairs_stage13_plus_13b": int(len(eval_pairs(both))),
           "ci_method": "log-OR (Woolf) 95% CI, 零格 0.5 连续性校正",
           "base_rate_heterogeneity_note": ("各留存基线 PASS 率差异大(SC2-02 2.4% / SC2-03 32.2%);"
                                            "合并配对级分析会被高基线留存主导而稀释富集对比,"
                                            "故同时报告 stage13_only / sc2_03_only / 合并三套以透明呈现;"
                                            "基因级 recall 对该异质性更稳健"),
           "outputs": ["holdout_topk_calibration.tsv"]},
          open(OUT / "holdout_topk_calibration.json", "w"), indent=2, ensure_ascii=False)

pd.set_option("display.width", 200)
print(res.to_string(index=False))
