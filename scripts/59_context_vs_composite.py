#!/usr/bin/env python3
"""Head-to-head: does cellular context or composite scoring predict replication?

Motivation
----------
The manuscript reports two numbers side by side that invite a comparison it never
made: holdout PASS is much higher in context-matched strata, and the composite
HG3 score has a gene-level AUROC of 0.447.  "Context outweighs composite
scoring" is a *comparative* claim and needs the two predictors evaluated on the
same genes against the same outcome.  This script does exactly that, and it
reports the ways in which the comparison is fragile.

Design
------
Outcome (gene level, 333 evaluable genes, matching the benchmark definition):
    validated = the gene passed in at least one same-pathogen holdout,
    including the supplementary SC2-03 arm.

Predictors, all evaluated on that same set of genes:
    HG3_total_score      the composite the resource proposes
    n_pathogens_hit      the simplest possible consensus count
    context_exposure     fraction of the gene's evaluable same-pathogen pairs
                         whose holdout matched the discovery cellular context
    direction            gene intervention direction (preserve/activate)
    n_pairs              number of evaluable pairs (a nuisance variable: a gene
                         with more pairs has more chances to be validated)

Three honesty checks are built in, because the pairwise context variable is
weakly identified:
  1. context_matched is constant within four of the six holdouts, so the
     between-holdout estimate of its effect is confounded with holdout identity.
     The within-holdout fixed-effect estimate is reported, and it is identified
     almost entirely by HO-V-SC2-01 and HO-V-SC2-03.
  2. Leave-one-holdout-out refits of the pair-level context effect.
  3. The holdout-level view (n = 6) is reported as descriptive only.

Outputs
  24_closing_analyses/context_vs_composite_gene.tsv
  24_closing_analyses/context_vs_composite_summary.json
"""
from __future__ import annotations

import os

import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = Path(os.environ.get(
    "HPICDA_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "24_closing_analyses"
RNG = np.random.default_rng(20260913)

VERDICTS = ("PASS", "FAIL", "CONTEXT_BOUNDARY")

# --------------------------------------------------------------- pair level
v13 = pd.read_csv(ROOT / "20_stage13_holdout_adjudication/stage13_gene_holdout_verdicts.tsv",
                  sep="\t")
v13b = pd.read_csv(ROOT / "23_holdout_sc2_03_stahl/stage13b_sc2_03_verdicts.tsv", sep="\t")
cols = ["holdout_id", "gene", "verdict", "context_matched", "same_pathogen",
        "intervention_direction", "bio_tier"]
pairs = pd.concat([v13[cols], v13b[cols]], ignore_index=True)
pairs = pairs[pairs["verdict"].isin(VERDICTS) & (pairs["same_pathogen"] == True)].copy()
pairs["y"] = (pairs["verdict"] == "PASS").astype(int)
pairs["cm"] = pairs["context_matched"].astype(int)
pairs["direction"] = pairs["intervention_direction"].fillna("UNKNOWN")
pairs["tier2"] = (pairs["bio_tier"] == "Bio-Tier 2").astype(int)
print(f"pair level: {len(pairs)} same-pathogen evaluable pairs, "
      f"{pairs['y'].sum()} PASS ({pairs['y'].mean():.1%})")

# --------------------------------------------------------------- gene level
bench = pd.read_csv(ROOT / "24_closing_analyses/benchmark_gene_level_data.tsv", sep="\t")
frozen = pd.read_csv(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", sep="\t")

ctx = (pairs.groupby("gene")
       .agg(n_pairs=("y", "size"), n_cm=("cm", "sum"), n_pass=("y", "sum"))
       .reset_index())
ctx["context_exposure"] = ctx["n_cm"] / ctx["n_pairs"]

g = bench.merge(ctx[["gene", "n_pairs", "n_cm", "context_exposure"]], on="gene", how="left")
g = g.merge(frozen[["gene", "intervention_direction"]], on="gene", how="left")
g = g[(g["evaluable"] == True)].copy()
g["context_exposure"] = g["context_exposure"].fillna(0.0)
g["preserve_activate"] = (g["intervention_direction"] == "PRESERVE_ACTIVATE").astype(int)
g["validated_i"] = g["validated"].astype(int)
g["HG3_total_score"] = g["total_score"]
print(f"gene level: {len(g)} evaluable genes, {g['validated_i'].sum()} validated "
      f"({g['validated_i'].mean():.1%})")
print(f"  distinct context_exposure values: {g['context_exposure'].nunique()}")
print(f"  genes with >=1 context-matched pair: {(g['n_cm'] > 0).sum()}")

COVS = ["HG3_total_score", "n_pathogens_hit", "n_units_hit", "hit_frequency",
        "context_exposure", "preserve_activate", "n_pairs"]

# ------------------------------------------------ AUROC with bootstrap CIs
def auroc_ci(y, s, n=4000):
    y = np.asarray(y); s = np.asarray(s, dtype=float)
    if len(np.unique(y)) < 2 or np.allclose(s, s[0]):
        return float("nan"), float("nan"), float("nan")
    pt = roc_auc_score(y, s)
    idx = np.arange(len(y)); boot = []
    for _ in range(n):
        b = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        boot.append(roc_auc_score(y[b], s[b]))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return pt, lo, hi


rows = []
y = g["validated_i"].values
for c in COVS:
    pt, lo, hi = auroc_ci(y, g[c].values)
    rows.append({"model": "single_variable", "predictor": c, "auroc": round(pt, 4),
                 "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                 "includes_chance": bool(lo <= 0.5 <= hi) if lo == lo else None})

# --------------------------------------- multi-predictor models + 5-fold CV
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

MODELS = {
    "n_pairs_only": ["n_pairs"],
    "HG3_only": ["HG3_total_score"],
    "context_only": ["context_exposure"],
    "direction_only": ["preserve_activate"],
    "HG3_plus_context": ["HG3_total_score", "context_exposure"],
    "HG3_plus_context_plus_direction": ["HG3_total_score", "context_exposure", "preserve_activate"],
    "naive_count": ["n_pathogens_hit"],
}
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
X_all = g[["HG3_total_score", "context_exposure", "preserve_activate",
           "n_pathogens_hit", "n_pairs"]].astype(float)
for name, feats in MODELS.items():
    X = X_all[feats].values
    pred = np.zeros(len(g))
    for tr, te in cv.split(X, y):
        pipe = make_pipeline(StandardScaler(),
                             LogisticRegression(max_iter=2000, C=1.0))
        pipe.fit(X[tr], y[tr])
        pred[te] = pipe.predict_proba(X[te])[:, 1]
    pt, lo, hi = auroc_ci(y, pred, n=2000)
    in_sample = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000,
                                                                  C=1.0)).fit(X, y)
    ins = roc_auc_score(y, in_sample.predict_proba(X)[:, 1])
    rows.append({"model": name, "predictor": "+".join(feats),
                 "auroc": round(pt, 4), "ci_lo": round(lo, 4),
                 "ci_hi": round(hi, 4), "in_sample_auroc": round(ins, 4),
                 "includes_chance": bool(lo <= 0.5 <= hi) if lo == lo else None})

# ------------------------- paired AUROC differences against HG3 (bootstrap)
# NOTE ON WHAT THIS NUMBER MEANS.  Each model here is a logistic fit on a single
# feature, and a single-feature logistic is monotone in that feature with the
# sign chosen to maximise the likelihood.  For the HG3 score that means the fit
# FLIPS THE SIGN (as-used AUROC 0.447 -> fitted 0.553).  So these deltas are
# SIGN-FREE comparisons: they ask whether HG3 carries as much rank information
# as the alternative once its direction is free, not whether the score as
# published beats it.  The as-used comparison (0.447 vs 0.574) is in
# scripts/60_headtohead_frozen.py, which also reports the frozen-only layer.
def paired_diff(feats_a, feats_b, n=3000):
    Xa, Xb = X_all[feats_a].values, X_all[feats_b].values
    pa = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(Xa, y).predict_proba(Xa)[:, 1]
    pb = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(Xb, y).predict_proba(Xb)[:, 1]
    d = roc_auc_score(y, pb) - roc_auc_score(y, pa)
    idx = np.arange(len(y)); boot = []
    for _ in range(n):
        b = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        boot.append(roc_auc_score(y[b], pb[b]) - roc_auc_score(y[b], pa[b]))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return d, lo, hi


diffs = []
for label, a, b in [
    ("context_minus_HG3", ["HG3_total_score"], ["context_exposure"]),
    ("direction_minus_HG3", ["HG3_total_score"], ["preserve_activate"]),
    ("context_plus_dir_minus_HG3", ["HG3_total_score"],
     ["context_exposure", "preserve_activate"]),
    ("naive_count_minus_HG3", ["HG3_total_score"], ["n_pathogens_hit"]),
    ("context_minus_naive_count", ["n_pathogens_hit"], ["context_exposure"]),
]:
    d, lo, hi = paired_diff(a, b)
    diffs.append({"comparison": label, "delta_auroc": round(d, 4),
                  "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                  "excludes_zero": bool(lo > 0 or hi < 0)})

# --------------------------------------------- pair-level context effect
ev = pairs.copy()
m_marg = smf.logit("y ~ cm + direction + tier2", data=ev).fit(disp=0)
SEP = "HO-V-IAV-02"
ev_fe = ev[ev["holdout_id"] != SEP].copy()
m_fe = smf.logit("y ~ cm + direction + tier2 + C(holdout_id)", data=ev_fe).fit(disp=0)


def grab(m, term):
    return {"logOR": float(m.params[term]),
            "OR": float(np.exp(m.params[term])),
            "ci_lo": float(np.exp(m.conf_int().loc[term, 0])),
            "ci_hi": float(np.exp(m.conf_int().loc[term, 1])),
            "p": float(m.pvalues[term])}


pair_level = {"marginal_cm": grab(m_marg, "cm"),
              "marginal_direction": grab(m_marg, "direction[T.PRESERVE_ACTIVATE]"),
              "within_holdoutFE_cm": grab(m_fe, "cm"),
              "within_holdoutFE_direction": grab(m_fe, "direction[T.PRESERVE_ACTIVATE]")}

# leave-one-holdout-out refits of the pair-level context OR
loho = []
for h in sorted(ev["holdout_id"].unique()):
    sub = ev[ev["holdout_id"] != h]
    try:
        mm = smf.logit("y ~ cm + direction + tier2", data=sub).fit(disp=0)
        loho.append({"holdout_dropped": h, "OR_cm": round(float(np.exp(mm.params["cm"])), 3),
                     "ci_lo": round(float(np.exp(mm.conf_int().loc["cm", 0])), 3),
                     "ci_hi": round(float(np.exp(mm.conf_int().loc["cm", 1])), 3),
                     "p": float(mm.pvalues["cm"])})
    except Exception as exc:  # separation in a small drop-out refit
        loho.append({"holdout_dropped": h, "OR_cm": None, "note": str(exc)[:80]})

# cm identifiability: within-holdout variation
within = (ev.groupby("holdout_id")["cm"].agg(["mean", "size"]).reset_index()
          .rename(columns={"mean": "cm_fraction", "size": "n_pair"}))
within["cm_varies_within"] = ev.groupby("holdout_id")["cm"].nunique().values > 1

# ------------------------------------------------ holdout-level descriptive
hold = (ev.groupby("holdout_id")
        .agg(n_pair=("y", "size"), pass_rate=("y", "mean"),
             cm_fraction=("cm", "mean"),
             pa_fraction=("direction", lambda s: (s == "PRESERVE_ACTIVATE").mean()))
        .reset_index())

# ------------------------------------------------------------------ write
res = pd.DataFrame(rows)
res.to_csv(OUT / "context_vs_composite_gene.tsv", sep="\t", index=False)

summary = {
    "analysis": "Head-to-head comparison of cellular context and composite scoring as predictors of replication",
    "date": "2026-09-13",
    "outcome": "gene-level validated (>=1 same-pathogen PASS incl. SC2-03), 333 evaluable genes",
    "n_genes_evaluable": int(len(g)),
    "n_genes_validated": int(g["validated_i"].sum()),
    "gene_level_auroc": rows,
    "paired_auroc_differences": diffs,
    "pair_level_logistic": pair_level,
    "pair_level_note": ("context_matched is constant within four of the six holdouts, so the "
                        "marginal OR is confounded with holdout identity; the fixed-effect "
                        "estimate is identified almost entirely by HO-V-SC2-01 and HO-V-SC2-03"),
    "context_identifiability": within.to_dict(orient="records"),
    "leave_one_holdout_out_cm_OR": loho,
    "holdout_level_descriptive": hold.to_dict(orient="records"),
}
with open(OUT / "context_vs_composite_summary.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2, ensure_ascii=False)

print("\n=== gene-level single-variable and model AUROC (5-fold CV) ===")
print(res.to_string(index=False))
print("\n=== paired AUROC differences (bootstrap 95% CI; SIGN-FREE, see note) ===")
print(pd.DataFrame(diffs).to_string(index=False))
print("\n=== pair-level context / direction effects ===")
for k, v in pair_level.items():
    print(f"  {k:28s} OR {v['OR']:6.3f} [{v['ci_lo']:.2f}, {v['ci_hi']:.2f}]  p={v['p']:.2e}")
print("\n=== cm identifiability ===")
print(within.to_string(index=False))
print("\n=== leave-one-holdout-out OR_cm ===")
print(pd.DataFrame(loho).to_string(index=False))
print("\n=== holdout level (descriptive, n=6) ===")
print(hold.round(3).to_string(index=False))
