"""HPICDA script 53 — 三案例结局表(2026-09-05,投稿封口 A4).

问题:用三个代表性案例讲清 HPICDA 的核心哲学——并非每个复发 CRISPR hit 都是治疗靶点。
  - ACE2        :robust_anchor(稳健锚点,跨情境可复现)
  - HS 模块 6 基因:context_dependent(情境依赖;RSV 保留、SC2 不泛化)
  - V-ATPase 9 亚基:confounded_module(基线必需性/混杂模块,机制真实但干预靶点差)
每个基因汇总:排序/层级/方向、各留存裁决(Stage13+13b)、ORCS 混杂、基线必需性、
治疗窗、成药性、Pareto 前沿,并赋解读标签。输出 case_studies_ace2_hs_vatpase.tsv + .md。
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "24_closing_analyses"
OUT.mkdir(exist_ok=True)

CASES = (["ACE2"]
         + ["B3GAT3", "EXT1", "EXT2", "EXTL3", "B4GALT7", "SLC35B2"]
         + ["ATP6AP1", "ATP6AP2", "ATP6V0B", "ATP6V0C", "ATP6V1A",
            "ATP6V1B2", "ATP6V1E1", "ATP6V1G1", "ATP6V1H"])
GROUP = {"ACE2": "robust_anchor"}
for g in ["B3GAT3", "EXT1", "EXT2", "EXTL3", "B4GALT7", "SLC35B2"]:
    GROUP[g] = "context_dependent_HS"
for g in ["ATP6AP1", "ATP6AP2", "ATP6V0B", "ATP6V0C", "ATP6V1A",
          "ATP6V1B2", "ATP6V1E1", "ATP6V1G1", "ATP6V1H"]:
    GROUP[g] = "confounded_module_VATPase"

fr = pd.read_csv(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2_annotated.tsv", sep="\t")
p15 = pd.read_csv(ROOT / "22_stage15_pareto/stage15_pareto_priority.tsv", sep="\t")
v13 = pd.read_csv(ROOT / "20_stage13_holdout_adjudication/stage13_gene_holdout_verdicts.tsv", sep="\t")
v13b = pd.read_csv(ROOT / "23_holdout_sc2_03_stahl/stage13b_sc2_03_verdicts.tsv", sep="\t")
v = pd.concat([v13, v13b], ignore_index=True)
EVAL = {"PASS", "FAIL", "CONTEXT_BOUNDARY"}

# 冻结排序名次
rank = (fr.sort_values(["total_score", "screen_evidence_pts"], ascending=False)
        .reset_index(drop=True))
rank["rank"] = np.arange(1, len(rank) + 1)
rank_map = dict(zip(rank["gene"], rank["rank"]))

rows = []
for g in CASES:
    f = fr[fr["gene"] == g]
    p = p15[p15["gene"] == g]
    vv = v[v["gene"] == g]
    ev = vv[vv["same_pathogen"] & vv["verdict"].isin(EVAL)]
    n_pass = int((ev["verdict"] == "PASS").sum())
    n_fail = int((ev["verdict"] == "FAIL").sum())
    n_cb = int((ev["verdict"] == "CONTEXT_BOUNDARY").sum())
    pass_holds = ",".join(ev.loc[ev["verdict"] == "PASS", "holdout_id"].tolist()) or "none"
    # 各留存一行式裁决串
    vstr = ";".join(f"{r.holdout_id}:{r.verdict}" for r in ev.itertuples()) or "no_evaluable"
    orcs_freq = vv["orcs_hit_freq_all"].dropna().max() if "orcs_hit_freq_all" in vv else np.nan
    row = {
        "gene": g, "case_group": GROUP[g],
        "rank_by_total_score": rank_map.get(g),
        "total_score": round(float(f["total_score"].iloc[0]), 3) if len(f) else np.nan,
        "bio_tier": f["bio_tier"].iloc[0] if len(f) else None,
        "intervention_direction": f["intervention_direction"].iloc[0] if len(f) else None,
        "benchmark_class": f["benchmark_class"].iloc[0] if len(f) and "benchmark_class" in f else None,
        "pathogens": f["pathogens"].iloc[0] if len(f) else None,
        "n_evaluable_pairs": int(len(ev)), "n_pass": n_pass, "n_fail": n_fail, "n_cb": n_cb,
        "pass_holdouts": pass_holds,
        "frac_backgrounds_essential": round(float(f["frac_backgrounds_essential"].iloc[0]), 3) if len(f) else np.nan,
        "min_baseline_chronos": round(float(f["min_baseline_chronos"].iloc[0]), 3) if len(f) else np.nan,
        "therapy_window_penalty": round(float(f["therapy_window_penalty"].iloc[0]), 3) if len(f) else np.nan,
        "orcs_hit_freq_all": round(float(orcs_freq), 4) if pd.notna(orcs_freq) else np.nan,
        "orcs_promiscuous_top5pct": (bool(p["orcs_promiscuous_top5pct"].iloc[0]) if len(p) and pd.notna(p["orcs_promiscuous_top5pct"].iloc[0]) else None),
        "pareto_front": int(p["pareto_front"].iloc[0]) if len(p) and pd.notna(p["pareto_front"].iloc[0]) else None,
        "druggability_class": p["druggability_class"].iloc[0] if len(p) else None,
        "x_tolerability": round(float(p["x_tolerability"].iloc[0]), 3) if len(p) and pd.notna(p["x_tolerability"].iloc[0]) else np.nan,
        "holdout_verdicts": vstr,
    }
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(OUT / "case_studies_ace2_hs_vatpase.tsv", sep="\t", index=False)

# ---- 简短 md 摘要 ----
def line(g):
    r = df[df["gene"] == g].iloc[0]
    return (f"| {g} | {r['case_group']} | rank {r['rank_by_total_score']} | "
            f"{r['bio_tier']} | {r['n_pass']}P/{r['n_fail']}F/{r['n_cb']}CB | "
            f"PASS: {r['pass_holdouts']} | ess={r['frac_backgrounds_essential']} | "
            f"ORCS_top5%={r['orcs_promiscuous_top5pct']} | {r['druggability_class']} |")

md = ["# HPICDA 三案例结局表(2026-09-05)", "",
      "核心哲学:并非每个复发 CRISPR hit 都是治疗靶点。三类结局——稳健锚点 / 情境依赖 / 混杂模块。", "",
      "| 基因 | 案例 | 排名 | 层级 | 留存裁决(P/F/CB) | PASS 留存 | 基线必需性 | ORCS混杂 | 成药性 |",
      "|---|---|---|---|---|---|---|---|---|"]
for g in CASES:
    md.append(line(g))
md += ["", "## 解读",
       "- **ACE2(robust_anchor)**:HG3 total_score 与 SLC35A1 并列第 1(8.00);stage15 复合分数在全部 2005 种权重方案下稳居第 1;在全部 3 个 SC2 留存(SC2-01/02/03)均 PASS(3P/0F/0CB);基线必需性 0;已有上市药物;方向明确(INHIBIT/受体阻断)。",
       "- **HS 模块(context_dependent)**:RSV 保留强信号,但 SC2 不泛化(SC2-03 0/6 PASS,EXT2 真反向)——RSV 偏好 + 细胞情境依赖,非广谱轴。",
       "- **V-ATPase(confounded_module)**:模块机制真实(内体酸化),但高基线必需性 + ORCS 高频混杂 → 治疗窗差,不宜作干预靶点。"]
(OUT / "case_studies_ace2_hs_vatpase.md").write_text("\n".join(md), encoding="utf-8")

json.dump({"script": "53_case_studies.py", "date": "2026-09-05", "n_case_genes": len(CASES),
           "groups": {"robust_anchor": ["ACE2"],
                      "context_dependent_HS": ["B3GAT3", "EXT1", "EXT2", "EXTL3", "B4GALT7", "SLC35B2"],
                      "confounded_module_VATPase": ["ATP6AP1", "ATP6AP2", "ATP6V0B", "ATP6V0C",
                                                    "ATP6V1A", "ATP6V1B2", "ATP6V1E1", "ATP6V1G1", "ATP6V1H"]},
           "outputs": ["case_studies_ace2_hs_vatpase.tsv", "case_studies_ace2_hs_vatpase.md"]},
          open(OUT / "case_studies_runinfo.json", "w"), indent=2, ensure_ascii=False)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
print(df[["gene", "case_group", "rank_by_total_score", "bio_tier", "n_pass", "n_fail", "n_cb",
          "pass_holdouts", "frac_backgrounds_essential", "orcs_promiscuous_top5pct",
          "druggability_class"]].to_string(index=False))
