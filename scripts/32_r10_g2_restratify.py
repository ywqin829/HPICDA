#!/usr/bin/env python3
"""R10: G2 表型重分层(计划 v3 Workstream A2)。

输入:
  08_human_genetics/stage10_genetics_evidence.tsv   (v1 分层 + 逐条证据 JSON)
  16_revision_r10_g2/phenotype_classification_v2.yaml (冻结的 disease_id 分类)

输出:
  16_revision_r10_g2/stage10_genetics_evidence_v2.tsv
    新增列:
      g2_class            infection_specific / respiratory_physiology / none
      genetics_tier_v2    仅 infection_specific 的 L2G>=0.5 证据可达 G2
      respiratory_physiology_genetics  披露 JSON(不计分证据)
      g2_evidence_v2      计分证据 JSON(infection_specific 子集)
  16_revision_r10_g2/r10_restratify_runinfo.json
"""
import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "16_revision_r10_g2"
L2G_THRESHOLD = 0.5  # 与 v1 冻结值一致(27_stage10_human_genetics.py)

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

cls = yaml.safe_load(open(OUT / "phenotype_classification_v2.yaml"))
id2class = {}
for cname, cdef in cls["classes"].items():
    for did in (cdef.get("disease_ids") or {}):
        id2class[did] = cname

df = pd.read_csv(ROOT / "08_human_genetics/stage10_genetics_evidence.tsv", sep="\t")

rows = []
n_unclassified = {}
for _, r in df.iterrows():
    g2ev = json.loads(r["g2_evidence"]) if pd.notna(r["g2_evidence"]) else []
    g3ev = json.loads(r["g3_evidence"]) if pd.notna(r["g3_evidence"]) else []
    # 分类每条证据
    def classify(evs):
        out = {"infection_specific": [], "respiratory_physiology": [], "other": []}
        for e in evs:
            c = id2class.get(e["disease_id"], "other")
            if c == "other":
                n_unclassified[e["disease_id"]] = n_unclassified.get(e["disease_id"], 0) + 1
            out[c].append(e)
        return out
    g2c = classify(g2ev)
    g3c = classify(g3ev)
    # v2 G2: 仅 infection_specific 且 score>=0.5(v1 已按 >=0.5 入 g2_evidence,此处再校验)
    g2_v2 = [e for e in g2c["infection_specific"] if e.get("score") is not None and e["score"] >= L2G_THRESHOLD]
    resp_phys = g2c["respiratory_physiology"] + g3c["respiratory_physiology"]
    if g2_v2:
        tier_v2, g2_class = "G2", "infection_specific"
    elif resp_phys:
        # 有呼吸生理证据但不计分:层级按剩余证据重判
        tier_v2 = "G3" if (g3c["infection_specific"] or g3c["other"] or g2c["other"]) else "G0"
        g2_class = "respiratory_physiology"
    else:
        tier_v2 = "G3" if (g3ev or g2c["other"]) else "G0"
        g2_class = "none" if not g2ev else "other"
    rows.append({
        "gene": r["gene"], "ensembl": r["ensembl"], "loeuf": r["loeuf"],
        "tier_v1": r["tier"], "genetics_tier_v2": tier_v2, "g2_class": g2_class,
        "g2_evidence_v2": json.dumps(g2_v2) if g2_v2 else "",
        "respiratory_physiology_genetics": json.dumps(resp_phys) if resp_phys else "",
        "direction_conflict": r["direction_conflict"], "notes": r["notes"] if pd.notna(r["notes"]) else "",
    })

res = pd.DataFrame(rows)
res.to_csv(OUT / "stage10_genetics_evidence_v2.tsv", sep="\t", index=False)

changed = res[res["tier_v1"] != res["genetics_tier_v2"]]
runinfo = {
    "revision": "R10", "date": "2026-09-02",
    "l2g_threshold": L2G_THRESHOLD,
    "classification_file": "16_revision_r10_g2/phenotype_classification_v2.yaml",
    "tier_counts_v1": res["tier_v1"].value_counts().to_dict(),
    "tier_counts_v2": res["genetics_tier_v2"].value_counts().to_dict(),
    "n_tier_changed": int(len(changed)),
    "changed_genes": changed[["gene", "tier_v1", "genetics_tier_v2", "g2_class"]].to_dict("records"),
    "n_unclassified_disease_ids": n_unclassified,
    "input_sha256": {
        "stage10_genetics_evidence.tsv": sha(ROOT / "08_human_genetics/stage10_genetics_evidence.tsv"),
        "phenotype_classification_v2.yaml": sha(OUT / "phenotype_classification_v2.yaml"),
    },
}
json.dump(runinfo, open(OUT / "r10_restratify_runinfo.json", "w"), indent=2, ensure_ascii=False)

print("tier v1:", runinfo["tier_counts_v1"])
print("tier v2:", runinfo["tier_counts_v2"])
print(f"tier changed: {len(changed)}")
print(changed[["gene", "tier_v1", "genetics_tier_v2", "g2_class"]].to_string(index=False))
print("unclassified disease_ids:", n_unclassified)
