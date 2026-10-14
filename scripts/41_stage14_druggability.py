#!/usr/bin/env python3
"""Stage 14 — druggability + tolerability integration (post-HG4, Plan v3 C2). v2
- Open Targets tractability (all modalities, no-arg field; batched aliases)
- Broad Repurposing Hub (datalake) for approved/clinical compounds per target
- Integrated with freeze v2 direction (§35 inhibit/preserve), R12 tolerability,
  R11 promiscuity, Stage 13 holdout rollup.
HS cluster reframed: pathogen-HS interface druggable, biosynthesis enzymes not.
"""
import json, os, time, datetime
import pandas as pd
import numpy as np
import requests

BASE = "/workspace/HPICDA"
OUT = f"{BASE}/21_stage14_druggability"
os.makedirs(OUT, exist_ok=True)
CACHE = f"{OUT}/cache_ot_ensembl.json"
GQL = "https://api.platform.opentargets.org/api/v4/graphql"
REHUB = "/mnt/datalake/broad_drug_repurposing_hub/broad_repurposing_hub_phase_moa_target_info.parquet"

frz = pd.read_csv(f"{BASE}/09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", sep="\t")
tol = pd.read_csv(f"{BASE}/16_revision_r12_tolerability/tolerability_axis.tsv", sep="\t")
r11 = pd.read_csv(f"{BASE}/16_revision_r11_orcs/r11_candidate_orcs.tsv", sep="\t")
s13 = pd.read_csv(f"{BASE}/20_stage13_holdout_adjudication/stage13_per_gene_rollup.tsv", sep="\t")
genes = frz["gene"].tolist()

# ---- 1. Ensembl mapping (reuse valid cache; search only unmapped) ----
ens_map = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
# migrate from v1 cache if present
old = f"{OUT}/cache_ot.json"
if os.path.exists(old) and not ens_map:
    ens_map = {g: r.get("ensembl") for g, r in json.load(open(old)).items()}
SEARCH_Q = """query($q:String!){ search(queryString:$q, entityNames:["target"], page:{index:0,size:8}){ hits{ id name entity } } }"""
def gql(query, variables=None, tries=3):
    for t in range(tries):
        try:
            r = requests.post(GQL, json={"query": query, "variables": variables or {}}, timeout=60)
            if r.status_code == 200:
                return r.json()
        except Exception:
            pass
        time.sleep(1.0 * (t + 1))
    return None
n_search = 0
for g in genes:
    if g not in ens_map:
        s = gql(SEARCH_Q, {"q": g})
        hits = (((s or {}).get("data") or {}).get("search") or {}).get("hits") or []
        exact = [h for h in hits if h.get("name", "").upper() == g.upper() and h.get("entity") == "target"]
        ens_map[g] = exact[0]["id"] if exact else None
        n_search += 1
        time.sleep(0.1)
json.dump(ens_map, open(CACHE, "w"))
mapped = {g: e for g, e in ens_map.items() if e}

# ---- 2. batched tractability ----
BATCH = 20
ens2gene = {e: g for g, e in mapped.items()}
ens_list = list(mapped.values())
tract = {}
for i in range(0, len(ens_list), BATCH):
    chunk = ens_list[i:i + BATCH]
    q = "query{ " + " ".join(
        f'g{j}: target(ensemblId:"{e}"){{ tractability{{ modality label value }} }}'
        for j, e in enumerate(chunk)) + " }"
    d = gql(q)
    data = (d or {}).get("data") or {}
    for j, e in enumerate(chunk):
        t = data.get(f"g{j}") or {}
        tr = t.get("tractability") or []
        tract[e] = {"sm": [x["label"] for x in tr if x.get("modality") == "SM" and x.get("value")],
                    "ab": [x["label"] for x in tr if x.get("modality") == "AB" and x.get("value")]}
    time.sleep(0.15)

# ---- 3. Repurposing Hub known-compound layer ----
hub = pd.read_parquet(REHUB)
hub["target"] = hub["target"].fillna("")
hub_exp = hub.assign(gene=hub["target"].str.split("|")).explode("gene")
hub_exp["gene"] = hub_exp["gene"].str.strip()
hub_c = hub_exp[hub_exp["gene"].isin(genes)]
PHASE_RANK = {"launched": 4, "phase 3": 3, "phase 2": 2, "phase 1": 1, "preclinical": 0}
def phase_rank(p):
    p = str(p).strip().lower()
    for k, v in PHASE_RANK.items():
        if k in p:
            return v
    return 0
hub_c = hub_c.copy()
hub_c["phase_rank"] = hub_c["clinical_phase"].map(phase_rank)
hub_agg = hub_c.groupby("gene").agg(
    n_compounds=("pert_iname", "nunique"),
    n_launched=("phase_rank", lambda s: int((s == 4).sum())),
    max_phase_rank=("phase_rank", "max"),
    moas=("moa", lambda s: ";".join(sorted(set(x for x in s if isinstance(x, str)))[:300]))).reset_index()
launched = (hub_c[hub_c["phase_rank"] == 4].groupby("gene")["pert_iname"]
            .apply(lambda s: ";".join(sorted(set(s))[:300])).reset_index().rename(columns={"pert_iname": "launched_drugs"}))
hub_agg = hub_agg.merge(launched, on="gene", how="left")

# ---- 4. integrate ----
ot = pd.DataFrame([{"gene": ens2gene[e], "ensembl": e,
                    "sm_tractable": ";".join(tract.get(e, {}).get("sm", [])),
                    "ab_tractable": ";".join(tract.get(e, {}).get("ab", []))} for e in ens_list])
tol_keep = ["gene"] + [c for c in ("loeuf", "pli", "human_intolerant", "mismatch_depmap_human",
                                   "mendelian_disease", "n_mendelian_diseases") if c in tol.columns]
df = (frz.merge(ot, on="gene", how="left")
         .merge(hub_agg, on="gene", how="left")
         .merge(tol[tol_keep], on="gene", how="left")
         .merge(r11[["gene", "orcs_hit_freq_all", "orcs_promiscuous_top5pct"]], on="gene", how="left")
         .merge(s13[["gene", "n_evaluable", "n_pass", "n_fail", "n_cb", "n_xpath_pass", "n_xpath_reversed"]],
                on="gene", how="left"))
df["sm_tractable_any"] = df["sm_tractable"].fillna("").str.len() > 0
df["ab_tractable_any"] = df["ab_tractable"].fillna("").str.len() > 0
df["druggability_class"] = np.select(
    [df["n_launched"].fillna(0) > 0, df["max_phase_rank"].fillna(0) >= 1,
     df["sm_tractable_any"] | df["ab_tractable_any"]],
    ["launched_drug_exists", "clinical_stage_compound", "tractable_no_clinical_compound"],
    default="no_tractability_evidence")
df["intervention_question"] = np.where(
    df["intervention_direction"].astype(str).str.upper().str.contains("INHIB"),
    "inhibit", "preserve_or_activate_not_for_intervention")
df.to_csv(f"{OUT}/stage14_druggability_axis.tsv", sep="\t", index=False)

HS = ["B3GAT3", "B4GALT7", "EXT1", "EXT2", "EXTL3", "SLC35B2"]
summ = {
    "stage": "Stage14_druggability_v2", "date": "2026-09-03",
    "n_genes": len(df), "n_ot_mapped": len(mapped),
    "n_with_tractability": int((df["sm_tractable_any"] | df["ab_tractable_any"]).sum()),
    "n_with_hub_compound": int(df["n_compounds"].notna().sum()),
    "n_launched": int((df["n_launched"].fillna(0) > 0).sum()),
    "druggability_by_tier": {str(k): v for k, v in df.groupby("bio_tier")["druggability_class"].apply(lambda s: s.value_counts().to_dict()).items()},
    "tier2_table": df[df["bio_tier"].astype(str).str.contains("Bio-Tier 2", na=False)][
        ["gene", "intervention_direction", "druggability_class", "n_launched", "launched_drugs",
         "sm_tractable", "ab_tractable", "moas"]].to_dict("records"),
    "hs_cluster": df[df["gene"].isin(HS)][["gene", "druggability_class", "sm_tractable", "ab_tractable",
                                           "n_launched", "launched_drugs", "moas"]].to_dict("records"),
    "inhibit_vs_preserve": df["intervention_question"].value_counts().to_dict(),
}
with open(f"{OUT}/stage14_summary.json", "w") as f:
    json.dump(summ, f, indent=2, default=str)
with open(f"{OUT}/stage14_runinfo.json", "w") as f:
    json.dump({"script": os.path.abspath(__file__), "timestamp": datetime.datetime.now().isoformat(),
               "sources": ["Open Targets Platform GraphQL (tractability, all modalities)",
                           "Broad Drug Repurposing Hub (datalake parquet): phase/MOA/target"],
               "inputs": ["BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", "tolerability_axis.tsv",
                          "r11_candidate_orcs.tsv", "stage13_per_gene_rollup.tsv"]}, f, indent=2)
print(json.dumps({"mapped": len(mapped), "n_tract": summ["n_with_tractability"],
                  "n_hub": summ["n_with_hub_compound"], "n_launched": summ["n_launched"],
                  "by_tier": summ["druggability_by_tier"], "hs": summ["hs_cluster"]},
                 indent=2, default=str)[:5000])
