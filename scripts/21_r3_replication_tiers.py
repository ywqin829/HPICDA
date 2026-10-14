#!/usr/bin/env python3
"""R3: four-tier replication redefinition + Stage 6 reclassification.

Tier vocabulary (frozen in replication_tier_definitions.yaml):
  T1_GUIDE_SET        same experiment, same cell line/virus, different guide-set
                      (e.g. GeCKOv2 half-library A vs B) -> technical robustness
  T2_SAME_STUDY       same independence cluster, different screen arms/conditions
                      -> within-study robustness, NOT independent replication
  T3_CROSS_STUDY_SAME_SYSTEM  different independence cluster, same cell line
  T4_CROSS_STUDY_DIFF_SYSTEM  different independence cluster AND different cell line

Inputs : 05_conditional_dependency/per_screen_calls_v2.tsv (R2),
         03_functional_screens/screen_estimand_registry.tsv
Outputs: 06_within_pathogen/within_pathogen_replication_v2.tsv
         06_within_pathogen/replication_tier_definitions.yaml
         07_cross_pathogen/convergent_gene_cluster_support.tsv
"""
import pandas as pd
from pathlib import Path

ROOT = Path("/workspace/HPICDA")
df = pd.read_csv(ROOT / "05_conditional_dependency/per_screen_calls_v2.tsv", sep="\t", dtype=str)
reg = pd.read_csv(ROOT / "03_functional_screens/screen_estimand_registry.tsv", sep="\t", dtype=str)

# screen_id -> metadata (first occurrence; arms share cluster/cell line)
meta = reg.drop_duplicates("screen_id").set_index("screen_id")
unit_meta = {}
for unit, g in df.groupby("unit_id"):
    sid = g["screen_id"].iloc[0]
    cl = meta["independence_cluster"].get(sid, "")
    cell = g["cell_line"].iloc[0]
    # RSV half-library tag: libA vs libB inside unit_id
    half = ""
    if "|libA" in unit or unit.endswith("|libA"):
        half = "libA"
    elif "|libB" in unit or unit.endswith("|libB"):
        half = "libB"
    unit_meta[unit] = {"cluster": cl, "cell_line": cell, "half": half}

# informative directional calls, protein-coding harmonized genes only
INF = {"CDEP_PATHOGEN_EXPLOITED": "EXP", "CDEP_HOST_PROTECTIVE": "PRO"}
d = df[df["call"].isin(INF)].copy()
d["dirn"] = d["call"].map(INF)
d["gene_h"] = d["gene_hgnc"].where(d["gene_hgnc"].fillna("") != "", d["gene"])
d = d[d["gene_class"].str.contains("protein-coding", na=False)]

rows = []
for (gene, path), g in d.groupby(["gene_h", "pathogen"]):
    units = g["unit_id"].unique().tolist()
    clusters = sorted({unit_meta[u]["cluster"] for u in units})
    cells = sorted({unit_meta[u]["cell_line"] for u in units})
    halves = {unit_meta[u]["half"] for u in units} - {""}
    # direction per cluster (a cluster supports a direction if any unit in it shows it)
    cl_dir = {}
    for u in units:
        c = unit_meta[u]["cluster"]
        cl_dir.setdefault(c, set()).update(g.loc[g["unit_id"] == u, "dirn"])
    exp_clusters = sorted([c for c, s in cl_dir.items() if "EXP" in s])
    pro_clusters = sorted([c for c, s in cl_dir.items() if "PRO" in s])
    both = [c for c, s in cl_dir.items() if len(s) > 1]
    n_cl = len(clusters)

    if both or (exp_clusters and pro_clusters):
        status, tier, cons = "WP_CONTRADICTORY", "", "CONFLICT"
    else:
        cons = "EXP" if exp_clusters else "PRO"
        if n_cl >= 2:
            # highest tier: T4 if any supporting pair of clusters differs in cell line
            sup_cells = sorted({unit_meta[u]["cell_line"] for u in units
                                if unit_meta[u]["cluster"] in (exp_clusters or pro_clusters)})
            tier = "T4_CROSS_STUDY_DIFF_SYSTEM" if len(sup_cells) > 1 else "T3_CROSS_STUDY_SAME_SYSTEM"
            status = "WP_REPLICATED"
        elif len(halves) == 2:
            tier, status = "T1_GUIDE_SET", "WP_ROBUST"
        elif len(units) >= 2:
            tier, status = "T2_SAME_STUDY", "WP_ROBUST"
        else:
            tier, status = "", "WP_SINGLE_EVIDENCE"
    rows.append({"gene": gene, "pathogen": path,
                 "pathogen_class": g["pathogen"].iloc[0] and
                 {"Influenza": "viral", "SARS-CoV-2": "viral", "RSV": "viral",
                  "Mtb/BCG": "bacterial", "Pseudomonas": "bacterial"}.get(path, ""),
                 "n_independence_clusters": n_cl,
                 "n_clusters_exploited": len(exp_clusters),
                 "n_clusters_protective": len(pro_clusters),
                 "consensus_direction": cons,
                 "replication_tier": tier,
                 "wp_status": status,
                 "supporting_units": ";".join(sorted(units)),
                 "supporting_clusters": ";".join(clusters),
                 "cell_lines": ";".join(cells)})

out = pd.DataFrame(rows).sort_values(["pathogen", "wp_status", "gene"])
out.to_csv(ROOT / "06_within_pathogen/within_pathogen_replication_v2.tsv", sep="\t", index=False)

print("=== wp_status v2 ===")
print(out["wp_status"].value_counts().to_string())
print("\n=== tier breakdown (replicated/robust only) ===")
print(out.loc[out["replication_tier"] != "", "replication_tier"].value_counts().to_string())
print("\n=== replicated by pathogen ===")
print(out[out["wp_status"] == "WP_REPLICATED"].groupby("pathogen").size().to_string())
print("\n=== RSV T1 guide-set robust pairs ===")
rsv = out[(out["pathogen"] == "RSV") & (out["replication_tier"] == "T1_GUIDE_SET")]
print("n =", len(rsv))
print(rsv[["gene", "consensus_direction"]].head(20).to_string(index=False))

# --- per-gene cluster support for the 40 threshold-intersection genes ---
conv = pd.read_csv(ROOT / "07_cross_pathogen/cross_pathogen_convergence.tsv", sep="\t")
sup_rows = []
for _, r in conv.iterrows():
    gene = r["gene"]
    tot_clusters = set()
    per_path = []
    for p in r["pathogens"].split(";"):
        sub = out[(out["gene"] == gene) & (out["pathogen"] == p)]
        if len(sub):
            cls = set(sub["supporting_clusters"].iloc[0].split(";"))
            tier = sub["replication_tier"].iloc[0]
        else:
            cls, tier = set(), ""
        tot_clusters |= cls
        per_path.append(f"{p}:{len(cls)}cl/{tier or 'single'}")
    sup_rows.append({"gene": gene, "direction": r["direction"],
                     "pathogens": r["pathogens"],
                     "total_independence_clusters": len(tot_clusters),
                     "per_pathogen_support": "; ".join(per_path)})
sup = pd.DataFrame(sup_rows)
sup.to_csv(ROOT / "07_cross_pathogen/convergent_gene_cluster_support.tsv", sep="\t", index=False)
print("\n=== convergent gene cluster support (40) ===")
print(sup[["gene", "direction", "total_independence_clusters", "per_pathogen_support"]].to_string(index=False))

# --- freeze tier vocabulary ---
yaml = """# R3: replication tier vocabulary (frozen 2026-08-31, plan v2)
tiers:
  T1_GUIDE_SET:
    definition: same experiment, same cell line and virus, different guide-set (e.g. GeCKOv2 half-library A vs B)
    interpretation: technical/guide-level robustness; NOT independent replication
  T2_SAME_STUDY:
    definition: same independence cluster (PMID/lab), different screen arms or conditions
    interpretation: within-study robustness; NOT independent replication
  T3_CROSS_STUDY_SAME_SYSTEM:
    definition: different independence clusters, same cell line
    interpretation: independent replication, shared-system confounds remain
  T4_CROSS_STUDY_DIFF_SYSTEM:
    definition: different independence clusters AND different cell lines
    interpretation: strongest replication tier available in this corpus
rules:
  - WP_REPLICATED requires >=2 independence clusters with consistent direction (T3 or T4)
  - T1/T2 evidence is reported as WP_ROBUST, never counted as replication
  - direction conflict across any evidence unit -> WP_CONTRADICTORY
  - UNAVAILABLE screens excluded from all evidence bases
"""
(ROOT / "06_within_pathogen/replication_tier_definitions.yaml").write_text(yaml)
print("\ntier vocabulary frozen.")
