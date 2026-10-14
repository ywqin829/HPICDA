#!/usr/bin/env python3
"""Stage 14 上市药注释机器初筛(2026-09-04,收尾项;专家/报告已披露需人工复核)。

对 druggability_class == launched_drug_exists 的 20 个基因逐条做药理学初筛分类:
  direct_target            药物已知主要/确立 MOA 直接作用于该基因产物
  direct_secondary         有直接作用报道但为次要/脱靶/混杂活性
  indirect_pathway         同通路/底物层面,不直接作用于该基因产物
  endogenous_ligand        内源性配体,非可干预药物
  likely_misannotation     疑似注释错误(经典陷阱)
附加标记:direction_mismatch(基因为 PRESERVE_ACTIVATE,抑制型药物方向相反)、
tier_context(ARCHIVE/REJECT 降级)、tolerability(基线必需)。
分类依据为药理学常识性证据的机器初筛,最终转化声明前仍需人工复核签字。

输出:
  21_stage14_druggability/stage14_launched_curation.tsv
  21_stage14_druggability/stage14_launched_curation.json
"""
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "21_stage14_druggability"
sys.path.insert(0, str(ROOT / "runner"))
from checkpoint import write_checkpoint

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

# 机器初筛分类表:(gene, drug) -> (class, rationale)
CUR = {
    ("ACE2", "lisinopril"): ("likely_misannotation", "ACE 抑制剂作用于 ACE(ACE1),不抑制 ACE2;同源酶混淆为经典注释陷阱"),
    ("ACE2", "moexipril"): ("likely_misannotation", "同上,ACE≠ACE2"),
    ("ACSL3", "icosapent"): ("indirect_pathway", "二十碳五烯酸乙酯为酰基-CoA 合成酶底物层面,非 ACSL3 抑制剂"),
    ("DGAT2", "niacin"): ("indirect_pathway", "烟酸为 NAD 前体/维生素,非 DGAT2 直接抑制剂"),
    ("EP300", "curcumin"): ("direct_secondary", "姜黄素对 p300 HAT 有报道抑制,但为混杂天然产物(PAINS 类),不可作靶向药物"),
    ("GNAS", "forskolin"): ("indirect_pathway", "福斯高林直接激活腺苷酸环化酶(AC),GNAS 为 AC 上游 G 蛋白;通路相邻非直接"),
    ("PDE6A", "caffeine"): ("indirect_pathway", "非选择性 PDE 抑制,PDE6A 为视网膜特异亚型;抑制即眼毒性机制,非治疗策略"),
    ("PDE6A", "dipyridamole"): ("indirect_pathway", "同上,非选择性 PDE 抑制"),
    ("PDE6A", "pentoxifylline"): ("indirect_pathway", "同上,非选择性 PDE 抑制"),
    ("PTK2B", "leflunomide"): ("direct_secondary", "主要 MOA 为 DHODH 抑制(特立氟胺);对 PTK2B/FAK2 酪氨酸激酶为报道的次要活性"),
    ("ROCK2", "fasudil"): ("direct_target", "法舒地尔为确立的 ROCK 抑制剂(日本上市,脑血管痉挛)"),
    ("ROCK2", "ripasudil"): ("direct_target", "利帕地尔为确立的 ROCK2 抑制剂(上市,青光眼)"),
    ("CTSL", "boceprevir"): ("direct_secondary", "HCV NS3 蛋白酶抑制剂,有对组织蛋白酶 L 脱靶抑制的报道;可共价结合半胱氨酸蛋白酶"),
    ("CTSL", "felbinac"): ("likely_misannotation", "NSAID(COX),与 CTSL 关联疑似宽松映射,需人工复核"),
    ("GRID1", "L-glutamic-acid"): ("endogenous_ligand", "谷氨酸为谷氨酸受体内源性配体,非可干预药物"),
    ("P4HB", "ribostamycin"): ("likely_misannotation", "氨基糖苷类抗生素与 PDI 关联弱,疑似宽松映射,需人工复核"),
    ("SLC5A1", "canagliflozin"): ("direct_target", "卡格列净对 SGLT1(SLC5A1)有实质抑制(双重 SGLT1/2);选择性警示:该类主要靶向 SLC5A2"),
    ("SLC5A1", "dapagliflozin"): ("direct_secondary", "主要靶向 SGLT2,对 SGLT1 活性弱"),
    ("SLC5A1", "empagliflozin"): ("direct_secondary", "主要靶向 SGLT2,对 SGLT1 活性弱"),
    ("TOP1", "SN-38"): ("direct_target", "伊立替康活性代谢物,确立 TOP1 毒剂"),
    ("TOP1", "aclarubicin"): ("direct_target", "确立 TOP1/II 抑制剂(抗肿瘤药)"),
    ("TOP1", "irinotecan"): ("direct_target", "确立 TOP1 抑制剂(上市抗肿瘤药)"),
    ("TOP1", "topotecan"): ("direct_target", "确立 TOP1 抑制剂(上市抗肿瘤药)"),
    ("CHUK", "acetylcysteine"): ("indirect_pathway", "氧化还原调节剂,对 IKK/NF-κB 为间接影响"),
    ("CHUK", "aminosalicylate"): ("indirect_pathway", "水杨酸类对 IKK 有弱报道,非直接确立 MOA"),
    ("CHUK", "mesalazine"): ("indirect_pathway", "同上,NF-κB 通路层面间接"),
    ("CHUK", "sulfasalazine"): ("indirect_pathway", "有 IKK 抑制报道但弱且非确立主要 MOA"),
    ("KCNA3", "dalfampridine"): ("direct_target", "4-氨基吡啶为广谱 Kv 通道阻滞剂,含 Kv1.3(KCNA3);非选择性"),
    ("KDM1A", "tranylcypromine"): ("direct_target", "机制性 KDM1A(LSD1)抑制剂(MAO 骨架),确立直接作用"),
    ("ROCK1", "fasudil"): ("direct_target", "同 ROCK2,确立 ROCK 抑制剂"),
    ("ROCK1", "ripasudil"): ("direct_target", "同 ROCK2,确立 ROCK 抑制剂"),
    ("CD44", "hyaluronic-acid"): ("endogenous_ligand", "透明质酸为 CD44 内源性配体,非抑制剂"),
    ("TRDMT1", "pentamidine"): ("likely_misannotation", "抗肺孢子虫药与 TRDMT1 关联疑似宽松映射,需人工复核"),
    ("ATP6V1A", "alendronate"): ("direct_secondary", "双膦酸盐主要 MOA 为 FPPS 抑制;对破骨细胞 V-ATPase 有报道的次要作用"),
    ("ATP6V1A", "etidronic-acid"): ("direct_secondary", "同上"),
    ("ATP6V1A", "tiludronate"): ("direct_secondary", "同上"),
}

d = pd.read_csv(OUT / "stage14_druggability_axis.tsv", sep="\t")
launched = d[d["druggability_class"] == "launched_drug_exists"]

rows = []
for _, r in launched.iterrows():
    g = r["gene"]
    for drug in str(r["launched_drugs"]).split(";"):
        drug = drug.strip()
        cls, rat = CUR.get((g, drug), ("needs_review", "未在初筛表中,需人工复核"))
        flags = []
        if r["intervention_direction"] == "PRESERVE_ACTIVATE" and cls.startswith("direct"):
            flags.append("direction_mismatch(基因为限制因子,抑制型药物方向相反)")
        if r["bio_tier"] in ("ARCHIVE", "REJECT"):
            flags.append(f"tier_context:{r['bio_tier']}(不进入转化优先级)")
        if pd.notna(r["min_baseline_chronos"]) and r["min_baseline_chronos"] < -2:
            flags.append("baseline_essential(治疗窗未决)")
        rows.append({"gene": g, "bio_tier": r["bio_tier"],
                     "intervention_direction": r["intervention_direction"],
                     "drug": drug, "curation_class": cls, "rationale": rat,
                     "flags": ";".join(flags)})

cur = pd.DataFrame(rows)
cur.to_csv(OUT / "stage14_launched_curation.tsv", sep="\t", index=False)

summ = {
    "analysis": "Stage 14 launched-drug annotation machine curation (first pass)",
    "date": "2026-09-04",
    "status": "机器初筛,转化声明前需人工复核签字",
    "n_genes": int(cur["gene"].nunique()), "n_pairs": len(cur),
    "class_counts": cur["curation_class"].value_counts().to_dict(),
    "direct_target_pairs": cur[cur["curation_class"] == "direct_target"][["gene", "drug"]].to_dict("records"),
    "actionable_leads_after_curation": [
        {"gene": "ROCK1/ROCK2", "drugs": "fasudil, ripasudil", "note": "确立直接抑制;Bio-Tier 3"},
        {"gene": "KDM1A", "drugs": "tranylcypromine", "note": "直接抑制但 direction_mismatch(PRESERVE_ACTIVATE),抑制方向与治疗假说相反"},
        {"gene": "SLC5A1", "drugs": "canagliflozin", "note": "直接但选择性警示(主靶 SGLT2)"},
        {"gene": "KCNA3", "drugs": "dalfampridine", "note": "直接但非选择性 Kv 阻滞"},
        {"gene": "PTK2B", "drugs": "leflunomide", "note": "次要活性;主 MOA 为 DHODH"},
        {"gene": "CTSL", "drugs": "boceprevir", "note": "脱靶报道"},
    ],
    "confirmed_traps": ["ACE2-lisinopril/moexipril(ACE≠ACE2)"],
    "input_sha256": {"stage14_axis": sha(OUT / "stage14_druggability_axis.tsv")},
}
json.dump(summ, open(OUT / "stage14_launched_curation.json", "w"), indent=2, ensure_ascii=False)

write_checkpoint("stage14_launched_curation",
                 str(ROOT / "scripts/47_stage14_launched_curation.py"),
                 [str(OUT / "stage14_druggability_axis.tsv")], None,
                 str(ROOT / "00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md"),
                 [str(OUT / "stage14_launched_curation.tsv"),
                  str(OUT / "stage14_launched_curation.json")],
                 ["stage14_druggability"])

print(f"{summ['n_genes']} 基因 {summ['n_pairs']} 药物对")
print("分类:", json.dumps(summ["class_counts"], ensure_ascii=False))
print("\n初筛后可行动线索:")
for l in summ["actionable_leads_after_curation"]:
    print(f"  {l['gene']} — {l['drugs']}:{l['note']}")
