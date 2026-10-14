#!/usr/bin/env python3
"""R1: per-screen gene-universe audit + HGNC harmonization map.

Outputs:
  05_conditional_dependency/gene_universe_audit.tsv   (one row per screen-unit)
  05_conditional_dependency/gene_symbol_map.tsv       (raw label -> HGNC approved symbol + class)
"""
import pandas as pd
import re
from pathlib import Path

ROOT = Path("/workspace/HPICDA")
CALLS = ROOT / "05_conditional_dependency/per_screen_calls.tsv"
HGNC = ROOT / "02_metadata/hgnc_complete_set.txt"
REG = ROOT / "03_functional_screens/screen_estimand_registry.tsv"
OUT_AUDIT = ROOT / "05_conditional_dependency/gene_universe_audit.tsv"
OUT_MAP = ROOT / "05_conditional_dependency/gene_symbol_map.tsv"

# ---- HGNC reference ----
hgnc = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
approved = set(hgnc.loc[hgnc["status"] == "Approved", "symbol"])
# alias / previous symbol -> approved symbol (Approved entries only)
alias2sym = {}
for _, r in hgnc[hgnc["status"] == "Approved"].iterrows():
    sym = r["symbol"]
    for col in ("alias_symbol", "prev_symbol"):
        v = r.get(col)
        if isinstance(v, str) and v:
            for a in v.split("|"):
                a = a.strip()
                if a and a not in approved:
                    alias2sym.setdefault(a, sym)
sym2locus = dict(zip(hgnc["symbol"], hgnc["locus_group"]))

MIR = re.compile(r"^(hsa-mir|hsa-let|mmu-mir)", re.I)
NTC = re.compile(r"^(NonTargeting|NTC|NEGCONTROL)", re.I)

def classify(g):
    g = str(g)
    if NTC.match(g):
        return "nontargeting_control", None
    if MIR.match(g):
        return "miRNA", None
    if g in approved:
        return "hgnc_" + sym2locus.get(g, "other").replace(" ", "_"), g
    if g in alias2sym:
        s = alias2sym[g]
        return "hgnc_alias_" + sym2locus.get(s, "other").replace(" ", "_"), s
    return "unmapped", None

# ---- load calls ----
df = pd.read_csv(CALLS, sep="\t", dtype=str)
reg = pd.read_csv(REG, sep="\t", dtype=str)
lib_map = {}
for _, r in reg.iterrows():
    lib_map[r["screen_id"]] = r.get("perturbation_type", "")

# ---- map every unique label once ----
uniq = df["gene"].unique().tolist()
rows = []
for g in uniq:
    cls, sym = classify(g)
    rows.append({"gene_raw": g, "gene_hgnc": sym if sym else "", "gene_class": cls})
gmap = pd.DataFrame(rows)
gmap.to_csv(OUT_MAP, sep="\t", index=False)

lut_cls = dict(zip(gmap["gene_raw"], gmap["gene_class"]))
lut_sym = dict(zip(gmap["gene_raw"], gmap["gene_hgnc"]))
df["gene_class"] = df["gene"].map(lut_cls)
df["gene_hgnc"] = df["gene"].map(lut_sym)

# ---- per-unit audit ----
aud = []
for unit, g in df.groupby("unit_id"):
    n = len(g)
    pc = g["gene_class"].str.contains("protein-coding", na=False).sum()
    mir = (g["gene_class"] == "miRNA").sum()
    ntc = (g["gene_class"] == "nontargeting_control").sum()
    unmapped = (g["gene_class"] == "unmapped").sum()
    other = n - pc - mir - ntc - unmapped
    sid = g["screen_id"].iloc[0]
    aud.append({
        "unit_id": unit,
        "screen_id": sid,
        "pathogen": g["pathogen"].iloc[0],
        "cell_line": g["cell_line"].iloc[0],
        "track": g["track"].iloc[0],
        "n_rows": n,
        "n_protein_coding": int(pc),
        "n_miRNA": int(mir),
        "n_nontargeting": int(ntc),
        "n_other_hgnc": int(other),
        "n_unmapped": int(unmapped),
        "pct_protein_coding": round(100 * pc / n, 2),
    })
aud = pd.DataFrame(aud).sort_values("unit_id")
aud.to_csv(OUT_AUDIT, sep="\t", index=False)

print("=== gene_universe_audit.tsv ===")
print(aud.to_string(index=False))
print()
print("total rows:", len(df), "| sum n_rows:", aud["n_rows"].sum())
print("class totals:", pd.Series(list(lut_cls.values())).value_counts().to_dict())
print("unique raw labels:", len(uniq), "| mapped to HGNC:", (gmap['gene_hgnc'] != '').sum())
