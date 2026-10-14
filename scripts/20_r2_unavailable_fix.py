#!/usr/bin/env python3
"""R2: UNAVAILABLE encoding fix + versioned v2 re-issue of Stage 5 outputs.

Changes vs v1:
  - adds gene_hgnc / gene_class columns (from R1 map)
  - ORCS-1998 (author Boolean hit column unrecoverable): all rows call=UNAVAILABLE,
    confidence='unavailable', intervention_implication='unknown'
  - scans for any other unit with zero informative calls and reports them
  - recomputes call counts + known-hit masking check with UNAVAILABLE excluded
    from denominators
Outputs (v2, v1 files untouched):
  05_conditional_dependency/per_screen_calls_v2.tsv
  05_conditional_dependency/per_screen_call_counts_v2.tsv
  05_conditional_dependency/stage5_screen_stats_v2.json
  07_cross_pathogen/known_hit_masking_check_v2.tsv
Appends AM2-12 to registry_amendments.tsv.
"""
import json
import pandas as pd
from pathlib import Path

ROOT = Path("/workspace/HPICDA")
FS = ROOT / "05_conditional_dependency"

df = pd.read_csv(FS / "per_screen_calls.tsv", sep="\t", dtype=str)
gmap = pd.read_csv(FS / "gene_symbol_map.tsv", sep="\t", dtype=str).fillna("")
lut_sym = dict(zip(gmap["gene_raw"], gmap["gene_hgnc"]))
lut_cls = dict(zip(gmap["gene_raw"], gmap["gene_class"]))
df["gene_hgnc"] = df["gene"].map(lut_sym)
df["gene_class"] = df["gene"].map(lut_cls)

# --- UNAVAILABLE reclassification ---
UNAVAIL_UNITS = {"ORCS-1998": "author Boolean hit column unrecoverable from ORCS table; "
                              "continuous FDR min=0.66 conflicts with 16 author-reported hits; "
                              "excluded from all denominators (AM2-12)"}
for u in UNAVAIL_UNITS:
    m = df["unit_id"] == u
    df.loc[m, "call"] = "UNAVAILABLE"
    df.loc[m, "confidence"] = "unavailable"
    df.loc[m, "intervention_implication"] = "unknown"

# --- scan: units with zero informative calls (informational only) ---
informative = {"CDEP_PATHOGEN_EXPLOITED", "CDEP_HOST_PROTECTIVE",
               "BASELINE_FITNESS_CONFOUNDED", "TOXIN_SPECIFIC", "SIGNAL_STATE_ALTERING"}
zero_units = []
for unit, g in df.groupby("unit_id"):
    if unit in UNAVAIL_UNITS:
        continue
    if not g["call"].isin(informative).any():
        zero_units.append(unit)
print("units with zero informative calls (excluding UNAVAILABLE):", zero_units)

# --- v2 outputs ---
df.to_csv(FS / "per_screen_calls_v2.tsv", sep="\t", index=False)

counts = df.groupby(["unit_id", "call"]).size().unstack(fill_value=0)
counts.to_csv(FS / "per_screen_call_counts_v2.tsv", sep="\t")
stats = {u: g["call"].value_counts().to_dict() for u, g in df.groupby("unit_id")}
with open(FS / "stage5_screen_stats_v2.json", "w") as f:
    json.dump(stats, f, indent=2)

# --- masking check v2 (UNAVAILABLE screens excluded from evidence) ---
mask = pd.read_csv(ROOT / "07_cross_pathogen/known_hit_masking_check.tsv", sep="\t")
called = df[df["call"].isin(informative)]
rec = []
for _, r in mask.iterrows():
    hit = called[(called["pathogen"] == r["pathogen"]) & (called["gene"] == r["gene"])]
    rec.append("RECOVERED" if len(hit) else "MISSED")
mask["status_v1"] = mask["status"]
mask["status"] = rec
mask["note"] = "v2: UNAVAILABLE screens (ORCS-1998) excluded from evidence base"
mask.to_csv(ROOT / "07_cross_pathogen/known_hit_masking_check_v2.tsv", sep="\t", index=False)
print("masking v2:", mask["status"].value_counts().to_dict())
changed = mask[mask["status"] != mask["status_v1"]]
print("changed by v2:", len(changed))

# --- amendment AM2-12 ---
am = pd.read_csv(ROOT / "03_functional_screens/registry_amendments.tsv", sep="\t", dtype=str)
new = pd.DataFrame([{
    "amendment_id": "AM2-12",
    "screen_id": "ORCS-1998",
    "field": "call_encoding",
    "old_value": "all NO_FUNCTIONAL_SIGNAL (0 calls)",
    "new_value": "UNAVAILABLE (whole screen)",
    "rationale": UNAVAIL_UNITS["ORCS-1998"],
}])
# align columns defensively
for c in am.columns:
    if c not in new.columns:
        new[c] = ""
new = new[am.columns]
am2 = pd.concat([am, new], ignore_index=True)
am2.to_csv(ROOT / "03_functional_screens/registry_amendments.tsv", sep="\t", index=False)
print("amendments now:", len(am2))
print("call totals v2:", df["call"].value_counts().to_dict())
