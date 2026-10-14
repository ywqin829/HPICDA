"""HPICDA script 51 — 情境依赖决定因素(2026-09-05,投稿封口 A2).

问题:同病原可裁决的基因×留存对中,哪些因素决定 PASS vs FAIL/CONTEXT_BOUNDARY?
数据:Stage 13(7 个留存)+ Stage 13b(SC2-03)的同病原可裁决行。
方法:logistic 回归 PASS ~ context_matched + intervention_direction + bio_tier + holdout_id,
报告 OR 与 95% CI;另给 context_matched × direction 分层 PASS 率表。
"""
import json
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from pathlib import Path

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "24_closing_analyses"
OUT.mkdir(exist_ok=True)

v13 = pd.read_csv(ROOT / "20_stage13_holdout_adjudication/stage13_gene_holdout_verdicts.tsv", sep="\t")
v13b = pd.read_csv(ROOT / "23_holdout_sc2_03_stahl/stage13b_sc2_03_verdicts.tsv", sep="\t")
cols = ["holdout_id", "gene", "verdict", "context_matched", "intervention_direction", "bio_tier"]
v = pd.concat([v13[cols + ["is_hs"]], v13b[cols + ["is_hs"]]], ignore_index=True)
ev = v[v["verdict"].isin(["PASS", "FAIL", "CONTEXT_BOUNDARY"])].copy()
ev["y"] = (ev["verdict"] == "PASS").astype(int)
ev["direction"] = ev["intervention_direction"].fillna("UNKNOWN")
ev["tier2"] = (ev["bio_tier"] == "Bio-Tier 2").astype(int)
ev["cm"] = ev["context_matched"].astype(int)
print(f"可裁决行 {len(ev)};PASS {ev['y'].sum()} ({ev['y'].mean():.1%});"
      f"情境匹配 {ev['cm'].sum()} 行")

# 分层表
strat = (ev.groupby(["cm", "direction"])["y"]
         .agg(["sum", "count", "mean"]).reset_index()
         .rename(columns={"sum": "n_pass", "count": "n_eval", "mean": "pass_rate"}))
strat["pass_rate"] = strat["pass_rate"].round(4)
strat.to_csv(OUT / "context_determinants_stratified.tsv", sep="\t", index=False)

# 按留存条目分层
byh = (ev.groupby("holdout_id")["y"].agg(["sum", "count", "mean"]).reset_index()
       .rename(columns={"sum": "n_pass", "count": "n_eval", "mean": "pass_rate"}))
byh["pass_rate"] = byh["pass_rate"].round(4)
byh.to_csv(OUT / "context_determinants_by_holdout.tsv", sep="\t", index=False)

def or_table(m, model_name):
    t = pd.DataFrame({
        "model": model_name,
        "term": m.params.index,
        "logOR": m.params.values.round(4),
        "OR": np.exp(m.params.values).round(3),
        "CI_lo": np.exp(m.conf_int()[0].values).round(3),
        "CI_hi": np.exp(m.conf_int()[1].values).round(3),
        "p": m.pvalues.values,
    })
    t["p"] = t["p"].map(lambda x: f"{x:.2e}")
    return t

# 模型 1(边际关联,主推断之一):不含 holdout 固定效应,干净收敛
m_marg = smf.logit("y ~ cm + direction + tier2", data=ev).fit(disp=0)

# 模型 2(holdout 内保守估计):含 holdout 固定效应。
# HO-V-IAV-02 为 0/27 PASS(完全分离),会使该系数趋向 -inf 并阻碍收敛;
# 该留存对 PASS/FAIL 对比不提供信息,故从 FE 模型剔除(其 0% 见 by_holdout 表)。
SEP = "HO-V-IAV-02"
ev_fe = ev[ev["holdout_id"] != SEP].copy()
m_fe = smf.logit("y ~ cm + direction + tier2 + C(holdout_id)", data=ev_fe).fit(disp=0)

or_tab = pd.concat([or_table(m_marg, "marginal_no_holdoutFE"),
                    or_table(m_fe, f"within_holdoutFE_excl_{SEP}")], ignore_index=True)
or_tab.to_csv(OUT / "context_determinants_logistic.tsv", sep="\t", index=False)

json.dump({"script": "51_context_determinants.py", "date": "2026-09-05",
           "n_evaluable_rows": int(len(ev)), "n_pass": int(ev["y"].sum()),
           "separation_note": (f"{SEP} 为 0/{int((ev['holdout_id']==SEP).sum())} PASS(完全分离),"
                               "已从固定效应模型剔除;其描述性 0% 见 by_holdout 表"),
           "model_note": ("marginal_no_holdoutFE = cm/direction/tier2 的边际关联;"
                          "within_holdoutFE = 控制留存基线后的保守(holdout 内)估计,为更严谨的主结果"),
           "outputs": ["context_determinants_logistic.tsv",
                       "context_determinants_stratified.tsv",
                       "context_determinants_by_holdout.tsv"]},
          open(OUT / "context_determinants_runinfo.json", "w"), indent=2, ensure_ascii=False)
print("=== 边际模型(无 holdout FE)===")
print(or_table(m_marg, "marginal").to_string(index=False))
print(f"\n=== holdout 内模型(剔除 {SEP})===")
print(or_table(m_fe, "within_FE").to_string(index=False))
print("\n=== 分层表 ===")
print(strat.to_string(index=False))
