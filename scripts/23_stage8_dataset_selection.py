#!/usr/bin/env python3
"""Stage 8 step 1: GEO dataset selection for the human relevance validation layer.

Per human_layer_rules.yaml: human samples, infection vs control (or severity),
bulk expression, n>=10, <=3 datasets per pathogen, selection logged for audit.
Queries GEO via E-utilities; outputs candidate metadata table for selection.
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "08_human_genetics/../09_transcriptomic_context"
OUT = ROOT / "09_transcriptomic_context"
OUT.mkdir(exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

QUERIES = {
    "Influenza": '(influenza[Title]) AND "Homo sapiens"[Organism] AND ("expression profiling by high throughput sequencing"[Filter] OR "expression profiling by array"[Filter]) AND (infected[Title/Abstract] OR infection[Title/Abstract] OR patient[Title/Abstract])',
    "SARS-CoV-2": '(COVID-19[Title] OR SARS-CoV-2[Title]) AND "Homo sapiens"[Organism] AND "expression profiling by high throughput sequencing"[Filter] AND (patient[Title/Abstract] OR blood[Title/Abstract] OR lung[Title/Abstract] OR airway[Title/Abstract])',
    "RSV": '(respiratory syncytial virus[Title] OR RSV[Title]) AND "Homo sapiens"[Organism] AND ("expression profiling by high throughput sequencing"[Filter] OR "expression profiling by array"[Filter]) AND (infection[Title/Abstract] OR patient[Title/Abstract] OR infant[Title/Abstract])',
    "Mtb": '(tuberculosis[Title]) AND "Homo sapiens"[Organism] AND ("expression profiling by high throughput sequencing"[Filter] OR "expression profiling by array"[Filter]) AND (blood[Title/Abstract] OR patient[Title/Abstract])',
}

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "HPICDA/1.0 (research)"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode()

rows = []
for pathogen, q in QUERIES.items():
    url = f"{EUTILS}/esearch.fcgi?db=gds&retmax=40&retmode=json&term={urllib.parse.quote(q)}"
    res = json.loads(get(url))
    ids = res.get("esearchresult", {}).get("idlist", [])
    print(f"{pathogen}: {len(ids)} hits")
    time.sleep(0.4)
    if not ids:
        continue
    url2 = f"{EUTILS}/esummary.fcgi?db=gds&retmode=json&id={','.join(ids)}"
    summ = json.loads(get(url2))["result"]
    time.sleep(0.4)
    for gid in ids:
        r = summ.get(gid, {})
        if r.get("entrytype") != "GSE":
            continue
        rows.append({
            "pathogen": pathogen,
            "gse": r.get("accession", ""),
            "title": r.get("title", ""),
            "summary": (r.get("summary", "") or "")[:300],
            "n_samples": r.get("n_samples", ""),
            "gdstype": r.get("gdstype", ""),
            "pdat": r.get("pdat", ""),
        })

df = pd.DataFrame(rows).drop_duplicates("gse")
df["n_samples"] = pd.to_numeric(df["n_samples"], errors="coerce")
df = df.sort_values(["pathogen", "n_samples"], ascending=[True, False])
df.to_csv(OUT / "stage8_geo_candidates.tsv", sep="\t", index=False)
for p in QUERIES:
    sub = df[df["pathogen"] == p].head(8)
    print(f"\n=== {p} top candidates ===")
    print(sub[["gse", "n_samples", "gdstype", "pdat", "title"]].to_string(index=False, max_colwidth=70))
