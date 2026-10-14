#!/usr/bin/env python3
"""R9: mock 臂交互项重算(计划 v3 Workstream A1)— MAGeCK RRA 实现。

方法说明(替代自研置换,理由记录于 runinfo):
  自研 guide 置换零分布在 B=10,000 时最小 p≈1e-4,对 ~1.8 万基因的 BH FDR
  分辨率不足(除 VERO 外全部对比 FDR<0.05 基因数为 0,属方法分辨率假象)。
  改用领域标准 MAGeCK RRA(负二项 guide 模型 + α-RRA 基因级聚合 + 经验零分布),
  与语料既有分析(ORCS-2066/2070/2071、mtb_lat_vs_ui)工具链一致。

对比设计:
  GSE178699(5 细胞系):
    Δ      = mageck test  感染(3) vs 未感染平行培养 mock(3)   → 交互项(感染特异选择)
    β_mock = mageck test  mock(3) vs T0(1)                     → 平行培养生长必需性
  GSE210689(Mtb THP-1,无 T0):
    Δ_Lat  = Lat(3) vs UI(3);Δ_Pro = Pro(3) vs UI(3)           → 对比内建 mock 校正

条件性重定义: Δ 显著(neg 或 pos FDR<0.05)且 β_mock 不显著(两侧 FDR≥0.05)。
  GSE210689 无 T0,β_mock 不可估计 → 仅 Δ 判定,如实披露。
验收: 核糖体基因(RPL/RPS)应 β_mock 显著而被新定义排除。

输出: 16_revision_r9_mock_interaction/
  mageck/                           每对比的 MAGeCK 原始输出
  r9_per_gene_interaction.tsv       每基因 × 对比:LFC/FDR/判定
  r9_candidate_overlay.tsv          353 候选覆盖
  r9_calibration.json               校准 + 核糖体验收
  r9_runinfo.json
"""
import hashlib
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "16_revision_r9_mock_interaction"
MK = OUT / "mageck"
MK.mkdir(parents=True, exist_ok=True)
MAGECK = "/workspace/.conda/envs/crispr/bin/mageck"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
t0 = time.time()

# ---------- 构建 MAGeCK 输入 ----------
raw = pd.read_csv(ROOT / "03_functional_screens/raw/GSE178699/GSE178699_covid_screens_readCounts.txt.gz", sep="\t")
raw = raw.rename(columns={"GENE_CLONE": "sgRNA", "GENE": "Gene"})

CELL_LINES = {
    "VERO":      dict(t0=["VERO.T0"], mock=["VERO.T7A", "VERO.T7B", "VERO.T7C"],
                      inf=["VERO.T16A.SARS.IC90", "VERO.T16B.SARS.IC90", "VERO.T16C.SARS.IC90"]),
    "HuH7":      dict(t0=["HuH7.T0"], mock=["HuH7.T14A", "HuH7.T14B", "HuH7.T14C"],
                      inf=["HuH7.T15A.SARS2", "HuH7.T15B.SARS2", "HuH7.T15C.SARS2"]),
    "CALU3":     dict(t0=["CALU3.T0"], mock=["CALU3.T25A", "CALU3.T25B", "CALU3.T25C"],
                      inf=["CALU3.T43A.lowMOI.SARS", "CALU3.T43B.lowMOI.SARS", "CALU3.T43C.lowMOI.SARS"]),
    "CACO2":     dict(t0=["CACO2.T0"], mock=["CACO2.T39A", "CACO2.T39B", "CACO2.T39C"],
                      inf=["CACO2.T37A.MOI10", "CACO2.T37B.MOI10", "CACO2.T37C.MOI10"]),
    "HEK293A_T": dict(t0=["HEK293A_T.T0"], mock=["HEK293A_T.T12A", "HEK293A_T.T12B", "HEK293A_T.T12C"],
                      inf=["HEK293A_T.T12A.MOI2", "HEK293A_T.T12B.MOI2", "HEK293A_T.T12C.MOI2"]),
}

import os
ENV = dict(os.environ)
ENV["PATH"] = "/workspace/.conda/envs/crispr/bin:" + ENV["PATH"]  # RRA 可执行文件所在

def run_mageck(count_file, treat, ctrl, prefix):
    cmd = [MAGECK, "test", "-k", str(count_file), "-t", ",".join(treat), "-c", ",".join(ctrl),
           "-n", str(prefix), "--norm-method", "median", "--remove-zero", "none",
           "--gene-lfc-method", "alphamean"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=ENV)
    if r.returncode != 0:
        raise RuntimeError(f"mageck failed for {prefix}:\n{r.stderr[-2000:]}")
    return Path(f"{prefix}.gene_summary.txt")

def load_gene_summary(gs):
    df = pd.read_csv(gs, sep="\t")
    df.columns = [c.replace("|", "_") for c in df.columns]
    out = pd.DataFrame({
        "gene": df["id"],
        "neg_fdr": df["neg_fdr"], "pos_fdr": df["pos_fdr"],
        "neg_lfc": df["neg_lfc"], "pos_lfc": df["pos_lfc"],
        "neg_p": df["neg_p-value"], "pos_p": df["pos_p-value"],
    })
    # 方向化:取显著侧;两侧均不显著时取 |lfc| 较大侧记录
    neg_sig = out["neg_fdr"] < 0.05
    pos_sig = out["pos_fdr"] < 0.05
    out["direction"] = np.select([neg_sig & ~pos_sig, pos_sig & ~neg_sig],
                                 ["depleted", "enriched"], default="ns")
    out["lfc"] = np.where(neg_sig, out["neg_lfc"], out["pos_lfc"])
    out["fdr_min"] = out[["neg_fdr", "pos_fdr"]].min(axis=1)
    return out

summaries = {}

# ---------- GSE178699 ----------
for cl, d in CELL_LINES.items():
    cols = d["t0"] + d["mock"] + d["inf"]
    cf = MK / f"counts_{cl}.txt"
    raw[["sgRNA", "Gene"] + cols].to_csv(cf, sep="\t", index=False)
    # Δ: 感染 vs mock
    gs = run_mageck(cf, d["inf"], d["mock"], MK / f"delta_{cl}")
    summaries[f"GSE178699_{cl}_delta"] = load_gene_summary(gs)
    # β_mock: mock vs T0
    gs = run_mageck(cf, d["mock"], d["t0"], MK / f"bmock_{cl}")
    summaries[f"GSE178699_{cl}_beta_mock"] = load_gene_summary(gs)
    print(f"{cl}: delta sig={int((summaries[f'GSE178699_{cl}_delta']['fdr_min']<0.05).sum())}, "
          f"bmock sig={int((summaries[f'GSE178699_{cl}_beta_mock']['fdr_min']<0.05).sum())}", flush=True)

# ---------- GSE210689 ----------
mtb_cf = ROOT / "03_functional_screens/raw/GSE210689/counts_mageck.txt"
for name, treat in [("delta_Lat", ["Lat1", "Lat2", "Lat3"]), ("delta_Pro", ["Pro1", "Pro2", "Pro3"])]:
    gs = run_mageck(mtb_cf, treat, ["UI1", "UI2", "UI3"], MK / f"mtb_{name}")
    summaries[f"GSE210689_THP1_{name}"] = load_gene_summary(gs)
    print(f"{name}: sig={int((summaries[f'GSE210689_THP1_{name}']['fdr_min']<0.05).sum())}", flush=True)

# ---------- 合并宽表 ----------
genes_all = sorted(set().union(*[set(s["gene"]) for s in summaries.values()]))
wide = pd.DataFrame({"gene": genes_all})
for name, s in summaries.items():
    wide = wide.merge(s[["gene", "lfc", "fdr_min", "direction"]].rename(
        columns={"lfc": f"lfc_{name}", "fdr_min": f"fdr_{name}", "direction": f"dir_{name}"}),
        on="gene", how="left")

# 条件性判定
for cl in CELL_LINES:
    d, m = f"fdr_GSE178699_{cl}_delta", f"fdr_GSE178699_{cl}_beta_mock"
    wide[f"conditional_{cl}"] = (wide[d] < 0.05) & ~(wide[m] < 0.05)
for name in ["delta_Lat", "delta_Pro"]:
    wide[f"conditional_THP1_{name}"] = wide[f"fdr_GSE210689_THP1_{name}"] < 0.05
cond_cols = [c for c in wide.columns if c.startswith("conditional_")]
wide["n_conditional_contexts"] = wide[cond_cols].sum(axis=1)
wide.to_csv(OUT / "r9_per_gene_interaction.tsv", sep="\t", index=False)

# ---------- 验收:核糖体基因 ----------
ribo = wide["gene"].str.match(r"^(RPL|RPS)\d")
accept = {}
for cl in CELL_LINES:
    m = f"fdr_GSE178699_{cl}_beta_mock"
    sub = wide[ribo & wide[m].notna()]
    accept[f"ribosomal_n_{cl}"] = int(len(sub))
    accept[f"ribosomal_beta_mock_sig_frac_{cl}"] = round(float((sub[m] < 0.05).mean()), 3) if len(sub) else None
    accept[f"ribosomal_conditional_frac_{cl}"] = round(float(sub[f"conditional_{cl}"].mean()), 3) if len(sub) else None

# ---------- 353 候选覆盖 ----------
cand = pd.read_csv(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv", sep="\t")
overlay = cand[["gene", "directions", "pathogens"]].merge(wide, on="gene", how="left")
overlay.to_csv(OUT / "r9_candidate_overlay.tsv", sep="\t", index=False)

calib = {"acceptance_ribosomal": accept,
         "n_sig_per_contrast": {k: int((v["fdr_min"] < 0.05).sum()) for k, v in summaries.items()}}
json.dump(calib, open(OUT / "r9_calibration.json", "w"), indent=2)
json.dump({"revision": "R9", "date": "2026-09-02", "method": "MAGeCK 0.5.9.5 RRA (alphamean LFC)",
           "method_substitution": "自研 guide 置换(B=1e4)FDR 分辨率不足,改用领域标准 MAGeCK RRA;"
                                  "与原计划'guide 级置换零分布'同为 guide 感知零分布,但具解析/经验 FDR",
           "runtime_min": round((time.time() - t0) / 60, 1),
           "input_sha256": {
               "GSE178699_counts": sha(ROOT / "03_functional_screens/raw/GSE178699/GSE178699_covid_screens_readCounts.txt.gz"),
               "GSE210689_counts": sha(ROOT / "03_functional_screens/raw/GSE210689/counts_mageck.txt")},
           "note": "GSE210689 无 T0;UI 平行臂使对比内建 mock 校正;β_mock 不可估计"},
          open(OUT / "r9_runinfo.json", "w"), indent=2, ensure_ascii=False)

print(f"\n总耗时 {(time.time()-t0)/60:.1f} min")
print("每对比显著基因数:", json.dumps(calib["n_sig_per_contrast"], indent=1))
print("核糖体验收:", json.dumps(accept, indent=1))
