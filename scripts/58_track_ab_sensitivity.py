#!/usr/bin/env python3
"""Track A / Track B comparability sensitivity analyses.

Motivation: the primary convergence statistic combines within-unit one-sided
rank p-values across units into a cluster Stouffer, then a pathogen Stouffer,
then a cross-pathogen Fisher. That chain assumes the per-unit statistics are
commensurable. Track B screens contribute author-reported gene-level results
whose hit-calling thresholds differ from study to study, so a reviewer can
reasonably ask whether the combination is being driven by calibration rather
than by biology.

Two analyses that do not rely on that assumption:

  A. one-vote-per-independent-unit (rank-based)
     Within each unit a gene is a direction-concordant hit if its signed score
     is positive and falls in the top Q of that unit's positive tail. Each
     INDEPENDENCE CLUSTER then casts a single vote (does any unit in the
     cluster call the gene). No p-value from any study enters the statistic,
     and no cluster can outvote another by contributing more screens.

  B. direction-only
     No threshold, no magnitude, no p-value. A gene is scored purely by the
     sign of its effect in each unit. A pathogen supports the gene when the
     majority of its units have the expected sign.

Both are compared against the same convergence gate (>=2 pathogens AND >=2
independence clusters) and against the frozen candidate set and the holdout
verdicts, so that the question "do the principal candidates survive?" is
answered on data rather than on argument.

Inputs
  05_conditional_dependency/per_screen_calls_v2.tsv
  03_functional_screens/screen_estimand_registry.tsv
  09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv
  20_stage13_holdout_adjudication/stage13_gene_holdout_verdicts.tsv

Outputs
  24_closing_analyses/track_ab_sensitivity_gene.tsv
  24_closing_analyses/track_ab_sensitivity_summary.json
"""
import os
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get(
    "HPICDA_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "24_closing_analyses"
Q = 0.01  # per-unit top fraction of the positive tail that counts as a hit

# ------------------------------------------------------------------ load
df = pd.read_csv(ROOT / "05_conditional_dependency/per_screen_calls_v2.tsv",
                 sep="\t", dtype=str)
df["score"] = pd.to_numeric(df["score"], errors="coerce")
d = df[(df["pathogen"].isin(["Influenza", "SARS-CoV-2", "RSV"]))
       & (df["endpoint_family"].isin({"HOST_FITNESS", "PATHOGEN_BURDEN"}))
       & (df["call"] != "UNAVAILABLE")
       & (df["gene_class"].str.contains("protein-coding", na=False))
       & (df["gene_hgnc"].fillna("") != "")].copy()
d = d[d["score"].notna()]
# expected direction: KO protects (enriched) is the pathogen-exploited sense
d["x"] = np.where(d["direction_lof"] == "enriched", d["score"], -d["score"])

reg = pd.read_csv(ROOT / "03_functional_screens/screen_estimand_registry.tsv",
                  sep="\t", dtype=str).drop_duplicates("screen_id").set_index("screen_id")

units = sorted(d["unit_id"].unique())
scr_of_unit = d.groupby("unit_id")["screen_id"].first()
unit_cluster = {u: reg["independence_cluster"].get(scr_of_unit[u], "") for u in units}
unit_pathogen = d.groupby("unit_id")["pathogen"].first().to_dict()
clusters = sorted(set(unit_cluster.values()))
pathogens = sorted(set(unit_pathogen.values()))
print(f"units={len(units)} clusters={len(clusters)} pathogens={pathogens}")

U = sorted(d["gene_hgnc"].unique())
nU = len(U)
print(f"gene universe = {nU}")

# --------------------------------------------------- per-unit signed matrix
uidx = {u: i for i, u in enumerate(units)}
gidx = {g: i for i, g in enumerate(U)}
X = np.full((len(units), nU), np.nan)
for u, g, x in zip(d["unit_id"], d["gene_hgnc"], d["x"]):
    X[uidx[u], gidx[g]] = x

# Analysis A: which genes are hits in each unit (top Q of the positive tail)
HIT = np.zeros((len(units), nU), dtype=bool)
for i in range(len(units)):
    row = X[i]
    ok = np.isfinite(row) & (row > 0)
    n_ok = int(ok.sum())
    if n_ok == 0:
        continue
    k = max(1, int(np.ceil(n_ok * Q)))
    idx = np.where(ok)[0]
    keep = idx[np.argsort(-row[idx])[:k]]
    HIT[i, keep] = True

# cluster-level vote: a cluster votes if any of its units calls the gene
cl_units = {c: [uidx[u] for u in units if unit_cluster[u] == c] for c in clusters}
pa_units = {p: [uidx[u] for u in units if unit_pathogen[u] == p] for p in pathogens}

voteA = np.zeros((len(clusters), nU), dtype=bool)
for ci, c in enumerate(clusters):
    voteA[ci] = HIT[cl_units[c]].any(axis=0)

cl_of = {c: [i for i, c2 in enumerate(clusters) if c2 == c] for c in clusters}
pa_of_cluster = {c: unit_pathogen[[u for u in units if unit_cluster[u] == c][0]]
                 for c in clusters}
pa_clusters = {p: [ci for ci, c in enumerate(clusters) if pa_of_cluster[c] == p]
               for p in pathogens}

nA_cluster = voteA.sum(axis=0)
nA_path = np.array([sum(voteA[ci].any() for ci in pa_clusters[p]) for p in pathogens])
A_path_count = np.zeros(nU, dtype=int)
for p in pathogens:
    A_path_count += np.array([voteA[pa_clusters[p]].any(axis=0)]).reshape(nU)
passA = (A_path_count >= 2) & (nA_cluster >= 2)

# Analysis B: direction only. Per unit, sign of x. Per pathogen, majority sign.
posB = np.zeros((len(pathogens), nU))
negB = np.zeros((len(pathogens), nU))
for pi, p in enumerate(pathogens):
    sub = X[pa_units[p]]
    posB[pi] = np.nansum(sub > 0, axis=0)
    negB[pi] = np.nansum(sub < 0, axis=0)
B_path_majority = posB > negB
B_path_count = B_path_majority.sum(axis=0)
# cluster-level concordance for the second half of the gate
voteB = np.zeros((len(clusters), nU), dtype=bool)
for ci, c in enumerate(clusters):
    sub = X[cl_units[c]]
    voteB[ci] = np.nansum(sub > 0, axis=0) > np.nansum(sub < 0, axis=0)
nB_cluster = voteB.sum(axis=0)
passB = (B_path_count >= 2) & (nB_cluster >= 2)

n_units_pos = (X > 0).sum(axis=0)
n_units_any = np.isfinite(X).sum(axis=0)
frac_pos = np.divide(n_units_pos, n_units_any, out=np.zeros(nU), where=n_units_any > 0)

# ------------------------------------------------------------ frozen sets
frozen = pd.read_csv(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv",
                     sep="\t")
frozen_genes = set(frozen["gene"])
tier2 = set(frozen.loc[frozen["bio_tier"] == "Bio-Tier 2", "gene"])
conv = pd.read_csv(ROOT / "07_cross_pathogen/convergence_v2_loso.tsv", sep="\t")
conv_genes = set(conv["gene"])

gene_df = pd.DataFrame({
    "gene": U,
    "voteA_clusters": nA_cluster,
    "voteA_pathogens": A_path_count,
    "passA_one_vote_per_unit": passA,
    "directionB_pathogens": B_path_count,
    "directionB_clusters": nB_cluster,
    "passB_direction_only": passB,
    "units_evaluated": n_units_any,
    "frac_units_positive": frac_pos,
})
gene_df["in_frozen_353"] = gene_df["gene"].isin(frozen_genes)
gene_df["bio_tier_2"] = gene_df["gene"].isin(tier2)
gene_df["convergent_21"] = gene_df["gene"].isin(conv_genes)
gene_df.to_csv(OUT / "track_ab_sensitivity_gene.tsv", sep="\t", index=False)

# ------------------------------------------------- holdout enrichment check
verdicts = pd.read_csv(ROOT / "20_stage13_holdout_adjudication"
                       / "stage13_gene_holdout_verdicts.tsv", sep="\t")
EV = {"PASS", "FAIL", "CONTEXT_BOUNDARY"}
vv = verdicts[verdicts["verdict"].isin(EV)].copy()
vv["passed"] = vv["verdict"] == "PASS"
pair = vv.groupby("gene").agg(n_pair=("passed", "size"),
                              n_pass=("passed", "sum")).reset_index()


def topk_enrichment(ranking_series, label, Ks=(10, 20, 50)):
    """pair-level PASS enrichment for the top-K genes of a ranking."""
    from scipy.stats import fisher_exact
    m = pd.DataFrame({"gene": U, "rank_score": ranking_series}).merge(pair, on="gene",
                                                                     how="inner")
    m = m.sort_values("rank_score", ascending=False, kind="mergesort")
    tot_p, tot_ok = int(m["n_pair"].sum()), int(m["n_pass"].sum())
    out = []
    for K in Ks:
        top = m.head(K)
        a, b_ = int(top["n_pass"].sum()), int(top["n_pair"].sum()) - int(top["n_pass"].sum())
        c_ = tot_ok - a
        d_ = tot_p - int(top["n_pair"].sum()) - c_
        OR, p = fisher_exact([[a, b_], [c_, d_]])
        out.append({"scheme": label, "K": K, "pairs_topk": int(top["n_pair"].sum()),
                    "pass_topk": a, "pass_rate_topk": round(a / int(top["n_pair"].sum()), 4),
                    "odds_ratio": round(OR, 3), "fisher_p_one_sided": f"{p:.3e}"})
    return out


enrich = []
enrich += topk_enrichment(gene_df["voteA_clusters"].values.astype(float),
                          "one_vote_per_independent_unit")
enrich += topk_enrichment(gene_df["directionB_pathogens"].values.astype(float)
                          + gene_df["directionB_clusters"].values.astype(float) / 1000.0,
                          "direction_only")
hg3 = frozen.set_index("gene")["total_score"]
enrich += topk_enrichment(pd.Series(hg3.reindex(U).values, index=range(nU)),
                          "HG3_total_score")

summary = {
    "analysis": "Track A/B comparability sensitivity (one-vote-per-unit; direction-only)",
    "date": "2026-09-13",
    "per_unit_hit_definition": f"top {Q:.0%} of the positive tail of the signed score",
    "units": len(units), "clusters": len(clusters), "pathogens": pathogens,
    "gene_universe": nU,
    "n_passA_one_vote_per_unit": int(passA.sum()),
    "n_passB_direction_only": int(passB.sum()),
    "frozen_353": len(frozen_genes),
    "convergent_21": len(conv_genes),
    "bio_tier_2": len(tier2),
    "recovery": {
        "A_all_353_frozen": int(gene_df.loc[gene_df["in_frozen_353"], "passA_one_vote_per_unit"].sum()),
        "A_12_biotier2": int(gene_df.loc[gene_df["bio_tier_2"], "passA_one_vote_per_unit"].sum()),
        "A_21_convergent": int(gene_df.loc[gene_df["convergent_21"], "passA_one_vote_per_unit"].sum()),
        "B_all_353_frozen": int(gene_df.loc[gene_df["in_frozen_353"], "passB_direction_only"].sum()),
        "B_12_biotier2": int(gene_df.loc[gene_df["bio_tier_2"], "passB_direction_only"].sum()),
        "B_21_convergent": int(gene_df.loc[gene_df["convergent_21"], "passB_direction_only"].sum()),
    },
    "overlap_between_A_and_B": int((passA & passB).sum()),
    "overlap_frozen_vs_A": int((gene_df["in_frozen_353"] & passA).sum()),
    "overlap_frozen_vs_B": int((gene_df["in_frozen_353"] & passB).sum()),
    "topk_enrichment": enrich,
}
with open(OUT / "track_ab_sensitivity_summary.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2, ensure_ascii=False)

print(json.dumps({k: v for k, v in summary.items() if k != "topk_enrichment"},
                 indent=2, ensure_ascii=False))
print()
print(pd.DataFrame(enrich).to_string(index=False))
