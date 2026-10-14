#!/usr/bin/env python3
"""Head-to-head, second pass: does the finding survive without SC2-03?

Why this script exists
----------------------
``59_context_vs_composite.py`` compared cellular context and the HG3 composite on
the gene-level benchmark whose outcome is *SC2-03-inclusive*.  Both the predictor
(``context_exposure``, a fraction computed over the gene's evaluable pairs) and
the outcome can be concentrated in the same supplementary arm, so a positive
result there is not by itself evidence that context generalises.  The same worry
applies to the pair-level context effect, and the interaction test made it
concrete: the headline context x direction cell is 77.2% SC2-03 pairs.

This script answers one question only: **which of the reported predictors survive
when the supplementary SC2-03 arm is removed entirely**, i.e. on the frozen,
pre-registered layer alone.  It is deliberately the same machinery as 59 so the
numbers are comparable, and it reports both layers side by side.

Outcomes
  outcome_all     validated = >=1 same-pathogen holdout PASS, frozen + SC2-03
  outcome_frozen  validated = >=1 same-pathogen holdout PASS, frozen only
Predictors (two context definitions, because context_exposure depends on the
pair set it is computed over)
  HG3_total_score, context_exposure_all, context_exposure_frozen,
  preserve_activate, n_pathogens_hit, n_pairs_all, n_pairs_frozen

Outputs
  24_closing_analyses/headtohead_frozen_gene.tsv
  24_closing_analyses/benchmark_two_layer.tsv
  24_closing_analyses/headtohead_frozen_summary.json
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
SEP = "HO-V-IAV-02"          # separation in the within-holdout FE fit


def auc_point(y, s):
    """Point AUROC only.  Must not bootstrap: it is called inside bootstrap loops."""
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    if len(np.unique(y)) < 2 or np.allclose(s, s[0]):
        return float("nan")
    return float(roc_auc_score(y, s))


def auroc(y, s, n=4000):
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    pt = auc_point(y, s)
    if pt != pt:
        return float("nan"), float("nan"), float("nan")
    idx = np.arange(len(y))
    boot = []
    for _ in range(n):
        b = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        boot.append(roc_auc_score(y[b], s[b]))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return pt, float(lo), float(hi)


# --------------------------------------------------------------- pair levels
v13 = pd.read_csv(ROOT / "20_stage13_holdout_adjudication/stage13_gene_holdout_verdicts.tsv", sep="\t")
v13b = pd.read_csv(ROOT / "23_holdout_sc2_03_stahl/stage13b_sc2_03_verdicts.tsv", sep="\t")
COLS = ["holdout_id", "gene", "verdict", "context_matched", "same_pathogen",
        "intervention_direction", "bio_tier"]


def prep(df):
    d = df[COLS].copy()
    d = d[d["verdict"].isin(VERDICTS) & (d["same_pathogen"] == True)]
    d["y"] = (d["verdict"] == "PASS").astype(int)
    d["cm"] = d["context_matched"].astype(int)
    d["direction"] = d["intervention_direction"].fillna("UNKNOWN")
    d["tier2"] = (d["bio_tier"] == "Bio-Tier 2").astype(int)
    return d


p_frozen = prep(v13)
p_all = pd.concat([p_frozen, prep(v13b)], ignore_index=True)
print(f"frozen pairs {len(p_frozen)}  PASS {p_frozen['y'].sum()} ({p_frozen['y'].mean():.1%})")
print(f"all pairs    {len(p_all)}  PASS {p_all['y'].sum()} ({p_all['y'].mean():.1%})")


def fit_pair(d, fe, label):
    sub = d[d["holdout_id"] != SEP] if fe else d
    rhs = "cm + direction + tier2" + (" + C(holdout_id)" if fe else "")
    try:
        m = smf.logit(f"y ~ {rhs}", data=sub).fit(disp=0)
    except Exception as exc:
        return {"layer": label, "fe": fe, "error": str(exc)[:120]}
    out = {"layer": label, "fe": fe, "n": int(len(sub)), "n_pass": int(sub["y"].sum())}
    for term, key in [("cm", "cm"), ("direction[T.PRESERVE_ACTIVATE]", "preserve_activate"),
                      ("tier2", "tier2")]:
        if term in m.params.index:
            out[key] = {"OR": round(float(np.exp(m.params[term])), 3),
                        "ci_lo": round(float(np.exp(m.conf_int().loc[term, 0])), 3),
                        "ci_hi": round(float(np.exp(m.conf_int().loc[term, 1])), 3),
                        "p": float(m.pvalues[term])}
    return out


pair_models = [fit_pair(p_frozen, False, "frozen"), fit_pair(p_frozen, True, "frozen"),
               fit_pair(p_all, False, "all"), fit_pair(p_all, True, "all")]

# number of holdouts in which cm varies (identifiability)
ident = []
for lbl, d in [("frozen", p_frozen), ("all", p_all)]:
    for h, sub in d.groupby("holdout_id"):
        ident.append({"layer": lbl, "holdout": h, "n": int(len(sub)),
                      "cm_frac": round(float(sub["cm"].mean()), 3),
                      "cm_varies": bool(sub["cm"].nunique() > 1),
                      "pass_rate": round(float(sub["y"].mean()), 3)})

# --------------------------------------------------------------- gene level
bench = pd.read_csv(OUT / "benchmark_gene_level_data.tsv", sep="\t")
frozen_tbl = pd.read_csv(ROOT / "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", sep="\t")


def gene_ctx(d, suffix):
    c = (d.groupby("gene").agg(**{f"n_pairs_{suffix}": ("y", "size"),
                                  f"n_cm_{suffix}": ("cm", "sum")}).reset_index())
    c[f"context_exposure_{suffix}"] = c[f"n_cm_{suffix}"] / c[f"n_pairs_{suffix}"]
    return c


ctx_f = gene_ctx(p_frozen, "frozen")
ctx_a = gene_ctx(p_all, "all")
val_f = p_frozen.groupby("gene")["y"].max().rename("validated_frozen").reset_index()
val_a = p_all.groupby("gene")["y"].max().rename("validated_all").reset_index()

g = (bench
     .merge(ctx_f, on="gene", how="left")
     .merge(ctx_a, on="gene", how="left")
     .merge(val_f, on="gene", how="left")
     .merge(val_a, on="gene", how="left")
     .merge(frozen_tbl[["gene", "intervention_direction"]], on="gene", how="left"))
g["validated_frozen"] = g["validated_frozen"].fillna(False)
g["validated_all"] = g["validated_all"].fillna(False)
g["preserve_activate"] = (g["intervention_direction"] == "PRESERVE_ACTIVATE").astype(int)
g["HG3_total_score"] = g["total_score"]

ev_all = g[g["evaluable"] == True].copy()
ev_frozen = g[(g["n_pairs_frozen"] > 0)].copy()
print(f"evaluable (benchmark flag): {len(ev_all)}  validated_all {int(ev_all['validated_all'].sum())}")
print(f"evaluable (>=1 frozen pair): {len(ev_frozen)}  validated_frozen {int(ev_frozen['validated_frozen'].sum())}")
print(f"  of the benchmark-evaluable genes, {(ev_all['n_pairs_frozen'] > 0).sum()} have a frozen pair")
print(f"  agreement validated_all vs bench.validated: "
      f"{int((ev_all['validated_all'] == ev_all['validated']).sum())}/{len(ev_all)}")

rows = []
for setname, sub, outcome in [("all_layer", ev_all, "validated_all"),
                              ("frozen_layer", ev_frozen, "validated_frozen")]:
    y = sub[outcome].values.astype(int)
    preds = ["HG3_total_score", "context_exposure_all", "context_exposure_frozen",
             "preserve_activate", "n_pathogens_hit", "n_pairs_all", "n_pairs_frozen"]
    for p in preds:
        if p not in sub:
            continue
        s = sub[p].astype(float).values
        if np.all(np.isnan(s)):
            continue
        s = np.nan_to_num(s, nan=0.0)
        pt, lo, hi = auroc(y, s)
        # sign-free information content: max(a, 1-a) is not a valid AUROC, so also
        # fit a sign-flipped copy to show the information is or is not there.
        pt_f, lo_f, hi_f = auroc(y, -s)
        rows.append({"analysis_set": setname, "outcome": outcome, "n": len(sub),
                     "n_pos": int(y.sum()), "prevalence": round(float(y.mean()), 4),
                     "predictor": p, "auroc": round(pt, 4),
                     "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                     "auroc_sign_flipped": round(pt_f, 4),
                     "includes_chance": bool(lo <= 0.5 <= hi)})

# paired AUROC difference against HG3 within each analysis set
for setname, sub, outcome in [("all_layer", ev_all, "validated_all"),
                              ("frozen_layer", ev_frozen, "validated_frozen")]:
    y = sub[outcome].values.astype(int)
    idx = np.arange(len(y))

    def series(p):
        return np.nan_to_num(sub[p].astype(float).values, nan=0.0)

    hg3 = series("HG3_total_score")
    for p in ["context_exposure_all", "context_exposure_frozen", "preserve_activate",
              "n_pathogens_hit"]:
        other = series(p)
        d = auc_point(y, other) - auc_point(y, hg3)
        boots = []
        for _ in range(3000):
            b = RNG.choice(idx, len(idx), replace=True)
            if len(np.unique(y[b])) < 2:
                continue
            a1 = auc_point(y[b], other[b])
            a2 = auc_point(y[b], hg3[b])
            if a1 == a1 and a2 == a2:
                boots.append(a1 - a2)
        lo, hi = np.percentile(boots, [2.5, 97.5])
        rows.append({"analysis_set": setname, "outcome": outcome, "n": len(sub),
                     "n_pos": int(y.sum()), "predictor": f"{p} - HG3_total_score",
                     "auroc": round(float(d), 4), "ci_lo": round(float(lo), 4),
                     "ci_hi": round(float(hi), 4),
                     "excludes_zero": bool(lo > 0 or hi < 0)})

res = pd.DataFrame(rows)
res.to_csv(OUT / "headtohead_frozen_gene.tsv", sep="\t", index=False)

# ------------------------------------------------ two-layer benchmark release
# The published benchmark table (benchmark_naive_vs_hpicda.tsv) carries only the
# SC2-03-inclusive layer, and Fig. 7c inherited that single layer.  The manuscript
# reports the composite on both layers, so the figure needs both.  Recompute the
# six ranking features on the frozen outcome as well.
FEATURES = ["HPICDA_HG3_total_score", "hit_frequency", "n_units_hit", "n_pathogens_hit",
            "n_screens_hit", "best_rank_pct"]
COL_OF = {"HPICDA_HG3_total_score": "HG3_total_score"}
for f in FEATURES:
    COL_OF.setdefault(f, f)
# best_rank_pct is a percentile: smaller is better, so its ranking score is its
# negative.  That is the sign convention script 50 uses for the published table,
# and the frozen-layer AUROC must use the same one to be comparable.
DESC = {"best_rank_pct": False}
bench_rank = pd.read_csv(OUT / "benchmark_naive_vs_hpicda.tsv", sep="\t")
ALL_AUROC = bench_rank.set_index("ranking")["auroc"].to_dict()
tl = []
for feat in FEATURES:
    col = COL_OF[feat]
    sign = -1.0 if not DESC.get(feat, True) else 1.0
    row = {"ranking": feat}
    for lbl, sub, outcome in [("all", ev_all, "validated_all"),
                              ("frozen", ev_frozen, "validated_frozen")]:
        y = sub[outcome].values.astype(int)
        s = sign * sub[col].astype(float).values
        pt, lo, hi = auroc(y, s)
        row[f"n_evaluable_{lbl}"] = int(len(sub))
        row[f"n_validated_{lbl}"] = int(y.sum())
        row[f"auroc_{lbl}"] = round(pt, 4)
        row[f"ci_lo_{lbl}"] = round(lo, 4)
        row[f"ci_hi_{lbl}"] = round(hi, 4)
    row["auroc_all_published"] = ALL_AUROC[feat]   # take S19's number as-is
    tl.append(row)
tl = pd.DataFrame(tl)
tl.to_csv(OUT / "benchmark_two_layer.tsv", sep="\t", index=False)
print("\n=== two-layer benchmark (gene level, 333 evaluable) ===")
print(tl.to_string(index=False))
_chk = (tl["auroc_all"] - tl["auroc_all_published"]).abs().max()
assert _chk < 1.5e-3, f"SC2-03-inclusive layer does not reproduce S19 (max diff {_chk:.4f})"
print(f"  SC2-03-inclusive layer reproduces S19 to {_chk:.4f} AUROC; "
      f"the column is kept as a cross-check, not a replacement.")

# pair-level models, flattened for release
prow = []
for m in pair_models:
    r = {"layer": m.get("layer"), "fixed_effects": bool(m.get("fe")), "n": m.get("n"),
         "n_pass": m.get("n_pass")}
    for k in ("cm", "preserve_activate", "tier2"):
        v = m.get(k)
        if isinstance(v, dict):
            r[f"{k}_OR"] = v["OR"]
            r[f"{k}_ci_lo"] = v["ci_lo"]
            r[f"{k}_ci_hi"] = v["ci_hi"]
            r[f"{k}_p"] = v["p"]
    prow.append(r)
pd.DataFrame(prow).to_csv(OUT / "headtohead_frozen_pair.tsv", sep="\t", index=False)

summary = {
    "analysis": "Does the context-vs-composite comparison survive removal of the supplementary SC2-03 arm?",
    "date": "2026-09-13",
    "pair_level_n": {"frozen": int(len(p_frozen)), "all": int(len(p_all))},
    "pair_level_models": pair_models,
    "holdout_identifiability": ident,
    "gene_level": rows,
    "note": ("context_matched is constant within most holdouts, so only the within-holdout "
             "fixed-effect estimate is identified; in the frozen layer that estimate rests on "
             "SC2-01 alone, because SC2-03 is not part of the frozen adjudication."),
}
with open(OUT / "headtohead_frozen_summary.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2, ensure_ascii=False)

print("\n=== pair-level models ===")
print(pd.DataFrame([{k: (v if not isinstance(v, dict) else json.dumps(v)) for k, v in m.items()}
                    for m in pair_models]).to_string(index=False))
print("\n=== holdout identifiability ===")
print(pd.DataFrame(ident).to_string(index=False))
print("\n=== gene-level AUROC / paired differences ===")
print(res.to_string(index=False))
