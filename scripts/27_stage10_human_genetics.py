#!/usr/bin/env python3
"""Stage 10 (fail-soft): human genetics evidence for pre-registered candidates.

Implements 08_human_genetics/genetics_rules.yaml (frozen v1):
  G1 = fine-mapped (PIP>=0.5) + coloc PP.H4>=0.7 + direction match
       -> NOT assessable via Platform API at scale; capped at G2 (rules:
          "fine-mapping 缺失 -> 最高 G2"). Recorded, not silently dropped.
  G2 = OT gwas_credible_sets (L2G aggregate) score >= 0.5 for a priority
       disease, or coding-variant source (clinvar / uniprot_variants)
  G3 = genetic association below G2 threshold (recorded only, no score bonus)
  G0 = no data = UNKNOWN (not negative)
GWAS Catalog is accessed via Open Targets ingestion (gwas_credible_sets
datasource); direct per-gene GWAS Catalog API is redundant and omitted
(recorded in runinfo).

Architecture: disease-first associatedTargets sweep (~15 diseases x pages)
+ per-gene geneticConstraint (LOEUF) queries. Fail-soft on every call.

Usage: python3 27_stage10_human_genetics.py run
"""
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "08_human_genetics"
RAW = OUT / "raw_api"
RAW.mkdir(parents=True, exist_ok=True)

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
L2G_THRESHOLD = 0.5          # frozen
PAGE = 1000
CODING_SOURCES = {"clinvar", "uniprot_variants"}       # coding-variant path -> G2
G3_SOURCES = {"gene_burden", "gwas_credible_sets"}     # below-threshold -> G3

PHENOTYPES = {  # frozen search terms
    "COVID-19 severity": ["critical COVID-19", "COVID-19"],
    "influenza severity": ["influenza"],
    "TB susceptibility/progression": ["tuberculosis"],
    "pneumonia": ["pneumonia", "bacterial pneumonia"],
    "sepsis": ["sepsis"],
    "respiratory failure": ["respiratory failure"],
}


def ot_query(query, variables=None, tries=3):
    for t in range(tries):
        try:
            r = requests.post(OT_URL, json={"query": query, "variables": variables or {}},
                              timeout=60)
            if r.status_code == 200 and "errors" not in r.json():
                return r.json()
        except Exception:
            pass
        time.sleep(2 * (t + 1))
    return None


def find_disease_ids():
    q = """query($q:String!){ search(queryString:$q, entityNames:["disease"])
      { hits { id name entity } } }"""
    mapping = {}
    for pheno, terms in PHENOTYPES.items():
        hits = []
        for term in terms:
            j = ot_query(q, {"q": term})
            if j:
                hits += [{"id": h["id"], "name": h["name"], "term": term}
                         for h in j["data"]["search"]["hits"][:3]]
        # dedupe by id
        seen, uniq = set(), []
        for h in hits:
            if h["id"] not in seen:
                seen.add(h["id"])
                uniq.append(h)
        mapping[pheno] = uniq
    return mapping


def sweep_disease(disease_id, candidates):
    """Return {symbol: [(disease_id, source, score), ...]} for candidate targets."""
    q = """query($id:String!,$idx:Int!){ disease(efoId:$id){
      associatedTargets(page:{index:$idx,size:%d}){ count rows{
        target{ id approvedSymbol } datasourceScores{ id score } } } } }""" % PAGE
    out = {}
    idx = 0
    while True:
        j = ot_query(q, {"id": disease_id, "idx": idx})
        if not j or not j["data"].get("disease"):
            return out, False  # fail-soft flag
        at = j["data"]["disease"]["associatedTargets"]
        for row in at["rows"]:
            sym = row["target"]["approvedSymbol"]
            if sym in candidates:
                for d in row["datasourceScores"]:
                    out.setdefault(sym, []).append(
                        {"disease_id": disease_id, "source": d["id"], "score": d["score"]})
        idx += 1
        if idx * PAGE >= at["count"] or idx >= 30:
            break
        time.sleep(0.2)
    return out, True


def get_loeuf(ensembl_id):
    q = """query($id:String!){ target(ensemblId:$id){
      geneticConstraint { constraintType oeUpper } } }"""
    j = ot_query(q, {"id": ensembl_id})
    if not j or not j["data"].get("target"):
        return None, False
    for c in j["data"]["target"]["geneticConstraint"] or []:
        if c["constraintType"] == "lof":
            return c["oeUpper"], True
    return None, True


def main():
    cand = pd.read_csv(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv", sep="\t")
    genes = sorted(set(cand["gene"]))
    cand_set = set(genes)
    hgnc = pd.read_csv(ROOT / "02_metadata/hgnc_complete_set.txt", sep="\t", low_memory=False)
    sym2ens = dict(zip(hgnc["symbol"], hgnc["ensembl_gene_id"]))

    print("mapping priority phenotypes to disease IDs...")
    dmap = find_disease_ids()
    with open(RAW / "disease_id_map.json", "w") as f:
        json.dump(dmap, f, indent=2)
    disease_ids = sorted({h["id"] for v in dmap.values() for h in v})
    print(f"  {len(disease_ids)} disease IDs: {disease_ids}")

    # disease-first sweep
    evidence = {g: [] for g in genes}
    sweep_fail = []
    for did in disease_ids:
        ev, ok = sweep_disease(did, cand_set)
        if not ok:
            sweep_fail.append(did)
            continue
        for sym, recs in ev.items():
            evidence[sym] += recs
        print(f"  {did}: {len(ev)} candidate targets hit")
    with open(RAW / "disease_sweep_evidence.json", "w") as f:
        json.dump(evidence, f)

    # per-gene LOEUF
    rows = []
    n_loeuf_fail = 0
    for i, g in enumerate(genes):
        ens = sym2ens.get(g)
        loeuf, ok = (None, False) if not ens or pd.isna(ens) else get_loeuf(ens)
        if not ok:
            n_loeuf_fail += 1
        recs = evidence.get(g, [])
        g2 = [r for r in recs if (r["source"] == "gwas_credible_sets"
                                  and r["score"] is not None and r["score"] >= L2G_THRESHOLD)
              or r["source"] in CODING_SOURCES]
        g3 = [r for r in recs if r not in g2 and r["source"] in G3_SOURCES | CODING_SOURCES]
        tier = "G2" if g2 else ("G3" if g3 else "G0")
        rows.append({"gene": g, "ensembl": ens, "loeuf": loeuf, "tier": tier,
                     "g2_evidence": json.dumps(g2) if g2 else "",
                     "g3_evidence": json.dumps(g3) if g3 else "",
                     "direction_conflict": False,
                     "notes": "" if ok else "loeuf_query_failed"})
        if (i + 1) % 100 == 0:
            print(f"  loeuf {i+1}/{len(genes)}")
            time.sleep(0.5)

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "stage10_genetics_evidence.tsv", sep="\t", index=False)
    meta = {"run_at_utc": datetime.now(timezone.utc).isoformat(),
            "n_genes": len(genes),
            "disease_ids": disease_ids,
            "sweep_failures": sweep_fail,
            "loeuf_query_failures": n_loeuf_fail,
            "tier_counts": res["tier"].value_counts().to_dict(),
            "l2g_threshold": L2G_THRESHOLD,
            "gwas_catalog_access": "via Open Targets gwas_credible_sets ingestion; direct API omitted (redundant)",
            "g1_note": "G1 not assessable (PIP/coloc unavailable via API); capped at G2 per frozen rules",
            "rule": "G0=unknown not negative"}
    with open(OUT / "stage10_runinfo.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
