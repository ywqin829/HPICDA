#!/usr/bin/env python3
"""R12: 人类可耐受性轴(计划 v3 Workstream A4)— 优化版(单次 primekg 读取)。

错配定义(预注册,运行前固定):
  depmap_conditional = 相关背景最小 Chronos > -1
  human_intolerant   = LOEUF < 0.35 或 pLI >= 0.9 或 primekg 疾病基因
                       或 genebass LoF 负担 p<1e-6 表型数 >= 1
  mismatch           = depmap_conditional AND human_intolerant
"""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "16_revision_r12_tolerability"
DL = Path("/mnt/datalake")
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
t0 = time.time()

# ---------- 1. gnomAD 约束 ----------
gn = pd.read_csv(OUT / "raw/gnomad.v4.1.constraint_metrics.tsv", sep="\t", low_memory=False)
gn = gn[gn["canonical"] == True]
gn = gn[["gene", "lof.oe", "lof.oe_ci.upper", "lof.pLI", "lof.z_score", "mis.oe", "mis.z_score"]]
gn.columns = ["gene", "lof_oe", "loeuf", "pli", "lof_z", "mis_oe", "mis_z"]
gn = gn.groupby("gene", as_index=False).first()
print(f"gnomAD: {len(gn)} 基因 ({(time.time()-t0)/60:.1f} min)", flush=True)

# ---------- 2. primekg 单次读取(disease_protein + bioprocess_protein)----------
pk = pd.read_csv(DL / "primekg/primekg.csv",
                 usecols=["relation", "x_name", "y_name", "y_source"], dtype=str)
dg = pk[pk["relation"] == "disease_protein"]
bp = pk[pk["relation"] == "bioprocess_protein"][["x_name", "y_name"]]
del pk
print(f"primekg 读取完成 ({(time.time()-t0)/60:.1f} min);disease_protein {len(dg)}, bioprocess {len(bp)}", flush=True)
# 方向判定:protein 是 gene
if dg["y_source"].str.contains("MONDO", na=False).mean() > 0.5:
    dg = dg.rename(columns={"x_name": "gene", "y_name": "disease"})
else:
    dg = dg.rename(columns={"y_name": "gene", "x_name": "disease"})
mendel = dg.groupby("gene")["disease"].agg(lambda s: ";".join(sorted(set(s))[:10])).reset_index()
mendel["n_diseases"] = mendel["disease"].str.count(";") + 1
bp.columns = ["term", "gene"]
del dg

# ---------- 3. genebass pLoF 负担(内存优化:pyarrow 过滤后读入) ----------
import pyarrow.parquet as pq
import pyarrow.compute as pc
_pq = DL / "genebass/genebass_pLoF_all_filtered.parquet"
# 第一遍:只读 gene + Pvalue_Burden(两列,~1.5 GB),向量化聚合
gb_min = pd.read_parquet(_pq, columns=["gene", "Pvalue_Burden"]).dropna(subset=["Pvalue_Burden"])
gb_min["is_sig"] = gb_min["Pvalue_Burden"] < 1e-6
agg = gb_min.groupby("gene", sort=False).agg(
    genebass_min_p=("Pvalue_Burden", "min"),
    genebass_n_sig_pheno=("is_sig", "sum")).reset_index()
agg["genebass_n_sig_pheno"] = agg["genebass_n_sig_pheno"].astype(int)
del gb_min
import gc; gc.collect()
# 第二遍:pyarrow 谓词下推只读显著行(p<1e-6 极少),取表型名
sig_t = pq.read_table(_pq, columns=["gene", "pheno_description"],
                      filters=[("Pvalue_Burden", "<", 1e-6)]).to_pandas()
sig = (sig_t.groupby("gene")["pheno_description"]
       .agg(lambda s: ";".join(sorted(set(s.dropna()))[:5])).reset_index()
       .rename(columns={"pheno_description": "genebass_sig_phenos"}))
agg = agg.merge(sig, on="gene", how="left")
del sig_t; gc.collect()
print(f"genebass: {len(agg)} 基因 ({(time.time()-t0)/60:.1f} min)", flush=True)

# ---------- 4. HPA ----------
hpa = pd.read_csv(DL / "human_protein_atlas/proteinatlas_subset.tsv", sep="\t",
                  usecols=["Gene", "RNA tissue specificity", "RNA tissue distribution",
                           "RNA tissue specificity score", "Disease involvement"])
hpa.columns = ["gene", "hpa_tissue_specificity", "hpa_tissue_distribution",
               "hpa_specificity_score", "hpa_disease_involvement"]

# ---------- 5. 全基因组背景 ----------
gw = gn.merge(mendel[["gene", "n_diseases"]], on="gene", how="left") \
       .merge(agg[["gene", "genebass_min_p", "genebass_n_sig_pheno"]], on="gene", how="left")
gw["n_diseases"] = gw["n_diseases"].fillna(0)
gw["genebass_n_sig_pheno"] = gw["genebass_n_sig_pheno"].fillna(0)
gw["human_intolerant"] = ((gw["loeuf"] < 0.35) | (gw["pli"] >= 0.9) |
                          (gw["n_diseases"] > 0) | (gw["genebass_n_sig_pheno"] >= 1))
gw.to_csv(OUT / "tolerability_genomewide.tsv", sep="\t", index=False)
print(f"全基因组: {len(gw)};intolerant {gw['human_intolerant'].mean():.3f} ({(time.time()-t0)/60:.1f} min)", flush=True)

# ---------- 6. 353 候选轴 ----------
cand = pd.read_csv(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv", sep="\t")
bl = pd.read_csv(ROOT / "04_depmap_baseline/baseline_chronos_by_background.tsv", sep="\t").set_index("gene")
REL = ["A-549", "Calu-3", "HuH-7", "Huh-7.5", "Huh-7.5.1", "HEp-2", "HeLa"]
rel = [c for c in REL if c in bl.columns]
axis = (cand[["gene", "directions", "pathogens"]].merge(gw, on="gene", how="left")
        .merge(hpa, on="gene", how="left")
        .merge(agg[["gene", "genebass_sig_phenos"]], on="gene", how="left")
        .merge(mendel[["gene", "disease"]], on="gene", how="left"))
axis["min_baseline_chronos"] = [float(bl.loc[g, rel].min()) if g in bl.index else np.nan for g in axis["gene"]]
axis["depmap_conditional"] = axis["min_baseline_chronos"] > -1
axis["mismatch"] = axis["depmap_conditional"] & (axis["human_intolerant"] == True)
axis.to_csv(OUT / "tolerability_axis.tsv", sep="\t", index=False)

# ---------- 7. 错配量化 ----------
dep = pd.read_csv(ROOT / "04_depmap_baseline/baseline_chronos_by_background.tsv", sep="\t")
dep["depmap_conditional"] = dep[rel].min(axis=1) > -1
gwm = gw.merge(dep[["gene", "depmap_conditional"]], on="gene", how="inner")
gwm["mismatch"] = gwm["depmap_conditional"] & gwm["human_intolerant"]
mm_rate_gw = float(gwm["mismatch"].mean())
mm_rate_cand = float(axis["mismatch"].mean())

from scipy.stats import fisher_exact
targets = {"glycosylation": bp[bp["term"].str.contains("glycosylation", case=False, na=False)]["gene"].unique(),
           "extracellular_matrix": bp[bp["term"].str.contains("extracellular matrix", case=False, na=False)]["gene"].unique(),
           "heparan_sulfate": bp[bp["term"].str.contains("heparan sulfate", case=False, na=False)]["gene"].unique()}
enrich = {}
mm_genes = set(gwm.loc[gwm["mismatch"], "gene"])
bg_genes = set(gwm["gene"])
for name, genes in targets.items():
    gs = set(genes) & bg_genes
    a = len(mm_genes & gs); b = len(mm_genes - gs)
    c = len((bg_genes - mm_genes) & gs); d = len(bg_genes - mm_genes - gs)
    if a + c > 0:
        orr, p = fisher_exact([[a, b], [c, d]])
        enrich[name] = {"n_term_genes": len(gs), "n_mismatch_in_term": a,
                        "odds_ratio": float(orr), "p": float(p)}
    else:
        enrich[name] = {"n_term_genes": len(gs), "n_mismatch_in_term": 0, "odds_ratio": None, "p": None}

json.dump({"mismatch_rate_genomewide": mm_rate_gw, "mismatch_rate_candidates": mm_rate_cand,
           "n_mismatch_genomewide": int(gwm["mismatch"].sum()),
           "n_mismatch_candidates": int(axis["mismatch"].sum()),
           "pathway_enrichment_in_mismatch": enrich,
           "definitions": {"depmap_conditional": "min Chronos(相关背景) > -1",
                           "human_intolerant": "LOEUF<0.35 | pLI>=0.9 | primekg 疾病基因 | genebass LoF p<1e-6"}},
          open(OUT / "r12_mismatch_analysis.json", "w"), indent=2, ensure_ascii=False)
json.dump({"revision": "R12", "date": "2026-09-02", "runtime_min": round((time.time() - t0) / 60, 1),
           "input_sha256": {"gnomad": sha(OUT / "raw/gnomad.v4.1.constraint_metrics.tsv")},
           "datalake_sources": ["primekg disease_protein/bioprocess_protein", "genebass pLoF", "HPA subset"],
           "limitation": "IMPC 小鼠敲除致死性未纳入(数据湖 mousemine 无活力表型);孟德尔病核实见 D3 引文核实"},
          open(OUT / "r12_runinfo.json", "w"), indent=2, ensure_ascii=False)

print(f"\n错配率: 全基因组 {mm_rate_gw:.3f} (n={int(gwm['mismatch'].sum())}), 候选 {mm_rate_cand:.3f} (n={int(axis['mismatch'].sum())})")
print("通路富集:", json.dumps(enrich, indent=1))
print("\n候选错配基因:")
mm = axis[axis["mismatch"] == True]
print(mm[["gene", "loeuf", "pli", "n_diseases", "genebass_n_sig_pheno", "min_baseline_chronos"]].to_string(index=False))
print(f"\n总耗时 {(time.time()-t0)/60:.1f} min")
