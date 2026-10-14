#!/usr/bin/env python3
"""R11: ORCS 频繁命中者协变量(计划 v3 Workstream A3)。

数据: BIOGRID-ORCS 人类全量屏包(1,952 屏,2.0.18;每屏全基因 HIT=YES/NO)。
计算:
  每基因 hit_freq_all  = 命中屏数 / 被测屏数(全部 1,952 屏)
  每基因 hit_freq_excl = 同上,但排除注册表中 49 个呼吸道感染 ORCS 屏
  promiscuous_flag     = hit_freq_all 处于全基因组 top 5%
输出:
  16_revision_r11_orcs/orcs_gene_hit_frequency.tsv   全基因组频率表
  16_revision_r11_orcs/r11_candidate_orcs.tsv        353 候选并表
  16_revision_r11_orcs/r11_runinfo.json
验收: 已知 promiscuous/必需通路(核糖体 RPL/RPS)频率显著高于背景。
"""
import hashlib
import json
import re
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace/HPICDA")
RAW = ROOT / "16_revision_r11_orcs/raw"
OUT = ROOT / "16_revision_r11_orcs"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
t0 = time.time()

# 注册表中的呼吸道感染 ORCS 屏 ID
reg = pd.read_csv(ROOT / "03_functional_screens/screen_estimand_registry.tsv", sep="\t")
resp_ids = set()
for c in ["screen_id", "study_id"]:
    for v in reg[c].astype(str):
        m = re.match(r"ORCS-(\d+)", v)
        if m:
            resp_ids.add(int(m.group(1)))
print(f"呼吸道感染屏 ID(注册表): {len(resp_ids)} 个")

files = sorted(RAW.glob("BIOGRID-ORCS-SCREEN_*.screen.tab.txt"))
print(f"屏文件: {len(files)}")

hit_all, assayed_all = Counter(), Counter()
hit_excl, assayed_excl = Counter(), Counter()
n_skipped = 0
for i, f in enumerate(files):
    sid = int(re.search(r"SCREEN_(\d+)", f.name).group(1))
    try:
        df = pd.read_csv(f, sep="\t", usecols=["OFFICIAL_SYMBOL", "HIT"], dtype=str)
    except Exception:
        n_skipped += 1
        continue
    df = df.dropna(subset=["OFFICIAL_SYMBOL"])
    df = df[df["HIT"].isin(["YES", "NO"])]
    genes = df["OFFICIAL_SYMBOL"]
    hits = genes[df["HIT"] == "YES"]
    assayed_all.update(genes)
    hit_all.update(hits)
    if sid not in resp_ids:
        assayed_excl.update(genes)
        hit_excl.update(hits)
    if (i + 1) % 400 == 0:
        print(f"  {i+1}/{len(files)} ({(time.time()-t0)/60:.1f} min)", flush=True)

all_genes = sorted(set(assayed_all) | set(assayed_excl))
res = pd.DataFrame({
    "gene": all_genes,
    "orcs_n_screens": [assayed_all.get(g, 0) for g in all_genes],
    "orcs_n_hits": [hit_all.get(g, 0) for g in all_genes],
    "orcs_n_screens_excl_resp": [assayed_excl.get(g, 0) for g in all_genes],
    "orcs_n_hits_excl_resp": [hit_excl.get(g, 0) for g in all_genes],
})
res["orcs_hit_freq_all"] = res["orcs_n_hits"] / res["orcs_n_screens"]
res["orcs_hit_freq_excl_resp"] = res["orcs_n_hits_excl_resp"] / res["orcs_n_screens_excl_resp"]
# top 5% promiscuous(在 ≥50 屏被测的基因中)
elig = res["orcs_n_screens"] >= 50
thr = res.loc[elig, "orcs_hit_freq_all"].quantile(0.95)
res["orcs_promiscuous_top5pct"] = elig & (res["orcs_hit_freq_all"] >= thr)
res.to_csv(OUT / "orcs_gene_hit_frequency.tsv", sep="\t", index=False)

# 验收:核糖体基因
ribo = res[res["gene"].str.match(r"^(RPL|RPS)\d", na=False)]
bg = res[elig]
acc = {"ribosomal_n": int(len(ribo)),
       "ribosomal_median_freq": float(ribo["orcs_hit_freq_all"].median()),
       "background_median_freq": float(bg["orcs_hit_freq_all"].median()),
       "ribosomal_top5pct_frac": float(ribo["orcs_promiscuous_top5pct"].mean()),
       "background_top5pct_frac": float(bg["orcs_promiscuous_top5pct"].mean()),
       "promiscuous_threshold": float(thr)}

# 353 候选并表
cand = pd.read_csv(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv", sep="\t")
cov = cand[["gene", "directions", "pathogens"]].merge(res, on="gene", how="left")
cov.to_csv(OUT / "r11_candidate_orcs.tsv", sep="\t", index=False)
n_cov = cov["orcs_n_screens"].notna().sum()

json.dump({"revision": "R11", "date": "2026-09-02",
           "n_screen_files": len(files), "n_skipped": n_skipped,
           "n_resp_excluded": len(resp_ids),
           "n_genes": len(res), "candidates_covered": f"{n_cov}/353",
           "acceptance": acc, "runtime_min": round((time.time() - t0) / 60, 1),
           "input_sha256": {"orcs_tarball": sha(RAW / "BIOGRID-ORCS-ALL-homo_sapiens.screens.tar.gz")},
           "note": "索引文件无法下载(命名不符);排除版基于注册表 49 个呼吸道感染 ORCS 屏"},
          open(OUT / "r11_runinfo.json", "w"), indent=2, ensure_ascii=False)

print(f"\n总耗时 {(time.time()-t0)/60:.1f} min;基因数 {len(res)};候选覆盖 {n_cov}/353")
print("验收:", json.dumps(acc, indent=1))
print("\n候选中 top5% promiscuous:")
print(cov[cov["orcs_promiscuous_top5pct"] == True][["gene", "orcs_hit_freq_all", "orcs_hit_freq_excl_resp"]]
      .sort_values("orcs_hit_freq_all", ascending=False).to_string(index=False))
