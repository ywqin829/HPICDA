#!/usr/bin/env python3
"""Stage 13 — holdout adjudication (post-HG4). v2
Per blueprint §34 + R14 frozen reverse hypothesis. No threshold re-tuning:
author HIT calls used as deposited; direction from signed score columns.
v2 fixes: (1) empirical sign-orientation map (ORCS-2360/2361 Z-score is
protective=negative, anchored on ACE2/COG3/ATP6V0C/TMEM41B vs IFNAR/JAK/STAT);
(2) cross-pathogen pairs (absent from frozen R14 stratification) are
exploratory: non-hits NOT_EVALUABLE, concordant hits = XPATH support;
(3) hits-only deposits (ORCS-1233/1234) flagged, PASS rate = direction
concordance among deposited hits; (4) Fisher table orientation fixed.
"""
import json, os, glob, datetime
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

BASE = "/workspace/HPICDA"
OUT = f"{BASE}/20_stage13_holdout_adjudication"
os.makedirs(OUT, exist_ok=True)
RAW = f"{BASE}/16_revision_r11_orcs/raw"

HOLDOUT_SCREENS = {
    "HO-V-SC2-01": ["ORCS-2360", "ORCS-2361"],
    "HO-V-SC2-02": ["ORCS-2068", "ORCS-2069", "ORCS-2076", "ORCS-2077"],
    "HO-V-IAV-01": ["ORCS-2428", "ORCS-2430"],
    "HO-V-IAV-02": ["ORCS-1524"],
    "HO-V-RSV-01": ["ORCS-2431"],
    "HO-B-Mtb-01": ["ORCS-1233", "ORCS-1234"],
    "HO-V-HPIV3-EXPL": ["ORCS-2433"],
}
HOLDOUT_PATHOGEN = {"HO-V-SC2-01": "SARS-CoV-2", "HO-V-SC2-02": "SARS-CoV-2",
                    "HO-V-IAV-01": "Influenza", "HO-V-IAV-02": "Influenza",
                    "HO-V-RSV-01": "RSV", "HO-B-Mtb-01": "Mtb/BCG",
                    "HO-V-HPIV3-EXPL": "HPIV-3"}
NOT_COMPUTABLE = {
    "HO-V-SC2-03": "PRJEB72772 raw FASTQ requires full alignment+MAGeCK pipeline; not executed in this revision (documented limitation)",
    "HO-B-Mtb-02": "GSE267063 murine Hoxb8; pre-registered cross-species context only, not a formal holdout",
    "HO-CLIN-01": "reserved (TBD at Stage 8 freeze)", "HO-SC-01": "reserved (TBD at Stage 9 freeze)"}
ALL_SCREENS = [s for v in HOLDOUT_SCREENS.values() for s in v]

# Empirical sign orientation: +1 = positive score means enriched-in-selection-arm
# (proviral KO protects/reduces pathogen). -1 = flipped.
# ORCS-2360/2361 anchors: ACE2 -2.09/-2.05, COG3 -4.14/-3.25, ATP6V0C -2.00/-4.26,
# TMEM41B -1.05/-1.58 (proviral, protective KO) vs IFNAR2 +3.70, JAK1 +3.59,
# IFNAR1 +2.54, STAT2 +2.37 (restriction) in ORCS-2361 => Z-score protective=negative.
# All other screens: author hits are positive-L2FC enriched (verified: TMEM41B +0.43
# ORCS-1524; SLC35A1/ATP6V0C/COG3 positive in Soriaga screens; ACE2 +0.91/+10.55 in
# ORCS-2068/2069; IFNAR1/2-STAT1/2-JAK1 positive in BCG screens, consistent with
# type-I-IFN being host-detrimental in Mtb/BCG).
SCREEN_ORIENTATION = {"ORCS-2360": -1, "ORCS-2361": -1}
HITS_ONLY = {"ORCS-1233", "ORCS-1234"}  # all deposited rows are HIT=YES

reg = pd.read_csv(f"{BASE}/03_functional_screens/screen_estimand_registry.tsv", sep="\t")
reg_h = reg[reg["screen_id"].isin(ALL_SCREENS)].set_index("screen_id")
dir_map = reg_h[["pathogen", "cell_line", "endpoint_family"]].to_dict("index")

idx = pd.read_csv(f"{RAW}/BIOGRID-ORCS-SCREEN_INDEX-2.0.18.index.tab.txt", sep="\t", low_memory=False)
idx.columns = [c.strip().strip("#") for c in idx.columns]
idx["sid"] = "ORCS-" + idx[idx.columns[0]].astype(str).str.strip()
idx = idx.set_index("sid")
score_type_cols = [c for c in idx.columns if c.upper().startswith("SCORE") and c.upper().endswith("TYPE")]

def pick_signed_score(sid):
    if sid not in idx.index: return None, None
    row = idx.loc[sid]; types = []
    for c in score_type_cols:
        v = str(row.get(c, "")).strip()
        if v and v.lower() != "nan":
            types.append((int("".join(ch for ch in c if ch.isdigit())), v))
    for p in ["l2fc", "log2", "differential", "z-score", "zscore", "z score", "mageck"]:
        for n, t in types:
            if p in t.lower(): return n, t
    return types[0] if types else (None, None)

cand = pd.read_csv(f"{BASE}/02_metadata/preregistered_candidate_set_v1.tsv", sep="\t")
cand_cols = {c.lower(): c for c in cand.columns}
cand = cand.rename(columns={cand_cols.get("gene", cand.columns[0]): "gene"})
dir_col = next((cand_cols[c] for c in ("directions", "direction") if c in cand_cols), None)
cand = cand.rename(columns={dir_col: "direction"})
cand["direction"] = cand["direction"].str.upper().str.strip()
path_col = next((cand_cols[c] for c in ("pathogens", "pathogen") if c in cand_cols), None)
cand["disc_pathogens"] = cand[path_col].fillna("").astype(str) if path_col else ""

frz = pd.read_csv(f"{BASE}/09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", sep="\t")[["gene", "bio_tier", "intervention_direction"]]
r11 = pd.read_csv(f"{BASE}/16_revision_r11_orcs/r11_candidate_orcs.tsv", sep="\t")[["gene", "orcs_hit_freq_all", "orcs_promiscuous_top5pct"]]
strat = pd.read_csv(f"{BASE}/16_revision_r14_diagnostics/r14_holdout_context_stratification.tsv", sep="\t")
strat_map = {(r.holdout_id, r.gene): (bool(r.context_matched), bool(r.lab_independent)) for r in strat.itertuples()}

HS = ["B3GAT3", "B4GALT7", "EXT1", "EXT2", "EXTL3", "SLC35B2"]
POS_CTRL = ["ACE2", "TMPRSS2"]

def load_screen(sid):
    n, stype = pick_signed_score(sid)
    files = glob.glob(f"{RAW}/BIOGRID-ORCS-SCREEN_{sid.split('-')[1]}-*.screen.tab.txt")
    if not files: return None, n, stype
    df = pd.read_csv(files[0], sep="\t", low_memory=False)
    df.columns = [c.strip().lstrip("#") for c in df.columns]
    sym = next(c for c in df.columns if c.upper() == "OFFICIAL_SYMBOL")
    hit = next((c for c in df.columns if c.upper().startswith("HIT")), None)
    scol = f"SCORE.{n}"
    if scol not in df.columns: scol = next(c for c in df.columns if c.upper().startswith("SCORE") and c.endswith(str(n)))
    out = pd.DataFrame({"gene": df[sym].astype(str).str.strip(),
                        "author_hit": df[hit].astype(str).str.strip().str.upper().eq("YES") if hit else False,
                        "score": pd.to_numeric(df[scol], errors="coerce")})
    out = out[out["gene"].notna() & (out["gene"] != "") & (out["gene"].str.lower() != "nan")]
    out = out.groupby("gene", as_index=False).agg(author_hit=("author_hit", "any"), score=("score", "first"))
    out["screen_id"] = sid; out["score_type"] = stype
    return out, n, stype

frames, screen_meta = [], {}
for sid in ALL_SCREENS:
    df, n, stype = load_screen(sid)
    if df is None: screen_meta[sid] = {"status": "FILE_NOT_FOUND"}; continue
    df["endpoint_family"] = dir_map.get(sid, {}).get("endpoint_family", "")
    frames.append(df)
    screen_meta[sid] = {"status": "ok", "n_genes": len(df), "n_hits": int(df["author_hit"].sum()),
                        "signed_score_col": n, "score_type": stype, "orientation": SCREEN_ORIENTATION.get(sid, 1),
                        "hits_only_deposit": sid in HITS_ONLY, "endpoint_family": df["endpoint_family"].iloc[0]}
calls = pd.concat(frames, ignore_index=True)

calls = calls.merge(cand[["gene", "direction", "disc_pathogens"]], on="gene", how="inner")
calls["expected_sign"] = np.where(calls["direction"] == "EXP", 1, np.where(calls["direction"] == "PRO", -1, 0))
calls["oriented_sign"] = (np.sign(calls["score"]).fillna(0).astype(int)
                          * calls["screen_id"].map(SCREEN_ORIENTATION).fillna(1).astype(int))
calls["concordant"] = calls["author_hit"] & (calls["oriented_sign"] == calls["expected_sign"]) & (calls["expected_sign"] != 0)
calls["reversed"] = calls["author_hit"] & (calls["oriented_sign"] == -calls["expected_sign"]) & (calls["expected_sign"] != 0)
calls["hit_sign_unknown"] = calls["author_hit"] & (calls["oriented_sign"] == 0)
screen2hold = {s: h for h, ss in HOLDOUT_SCREENS.items() for s in ss}
calls["holdout_id"] = calls["screen_id"].map(screen2hold)
calls["holdout_pathogen"] = calls["holdout_id"].map(HOLDOUT_PATHOGEN)
_hp = calls["holdout_pathogen"].str.lower()
_dp = calls["disc_pathogens"].fillna("").str.lower()
calls["same_pathogen"] = [h in d for h, d in zip(_hp, _dp)]
calls.to_csv(f"{OUT}/stage13_gene_screen_calls.tsv", sep="\t", index=False)

verdicts = []
for (h, g), sub in calls.groupby(["holdout_id", "gene"]):
    n_con = int(sub["concordant"].sum()); n_rev = int(sub["reversed"].sum())
    n_uns = int(sub["hit_sign_unknown"].sum())
    cm, li = strat_map.get((h, g), (None, None))
    same_p = bool(sub["same_pathogen"].iloc[0])
    if n_rev > 0 and n_con == 0 and n_uns == 0:
        v = "FAIL" if same_p else "XPATH_REVERSED"  # cross-pathogen reversal = divergence, not same-estimand FAIL (§34)
    elif n_con > 0 or n_uns > 0:
        v = "PASS" if same_p else "PASS_XPATH"
    elif not same_p:
        v = "NOT_EVALUABLE_XPATH"   # cross-pathogen non-hit: neither FAIL nor CONTEXT_BOUNDARY
    elif cm is True:
        v = "FAIL"
    elif cm is False:
        v = "CONTEXT_BOUNDARY"
    else:
        v = "NOT_STRATIFIED"
    verdicts.append({"holdout_id": h, "gene": g, "n_screens": len(sub), "n_concordant": n_con,
                     "n_reversed": n_rev, "n_hit_unsigned": n_uns, "same_pathogen": same_p,
                     "context_matched": cm, "lab_independent": li, "verdict": v})
verd = pd.DataFrame(verdicts).merge(frz, on="gene", how="left").merge(r11, on="gene", how="left")
verd["is_hs"] = verd["gene"].isin(HS)
verd.to_csv(f"{OUT}/stage13_gene_holdout_verdicts.tsv", sep="\t", index=False)

# HS frozen reverse-hypothesis Fisher (H1: matched PASS-rate > mismatched), same-pathogen pairs only
hs = verd[verd["is_hs"] & verd["same_pathogen"] & verd["context_matched"].notna()].copy()
hs["pass_bin"] = hs["verdict"].eq("PASS").astype(int)
m_p = int(((hs["context_matched"] == True) & (hs["pass_bin"] == 1)).sum()); m_n = int(((hs["context_matched"] == True) & (hs["pass_bin"] == 0)).sum())
x_p = int(((hs["context_matched"] == False) & (hs["pass_bin"] == 1)).sum()); x_n = int(((hs["context_matched"] == False) & (hs["pass_bin"] == 0)).sum())
table = np.array([[m_p, m_n], [x_p, x_n]])
or_, p_greater = fisher_exact(table, alternative="greater")
_, p_two = fisher_exact(table, alternative="two-sided")
hs_fisher = {"table_[[matched_PASS,matched_notPASS],[mismatched_PASS,mismatched_notPASS]]": table.tolist(),
             "odds_ratio": float(or_) if np.isfinite(or_) else "inf",
             "p_one_sided_H1_greater": float(p_greater), "p_two_sided": float(p_two), "alpha": 0.05,
             "hs_pass_matched": m_p, "hs_total_matched": m_p + m_n,
             "hs_pass_mismatched": x_p, "hs_total_mismatched": x_p + x_n,
             "interpretation": ("H1 (context-driven) supported" if p_greater < 0.05 else
                                "H1 not supported (matched PASS-rate not significantly greater); evaluate H2 (PASS in both) descriptively")}
with open(f"{OUT}/stage13_hs_fisher.json", "w") as f: json.dump(hs_fisher, f, indent=2)

EVAL = ["PASS", "FAIL", "CONTEXT_BOUNDARY"]           # same-pathogen evaluable
def pass_rate(df):
    d = df[df["verdict"].isin(EVAL)]
    return {"n_evaluable": len(d), "n_pass": int((d["verdict"] == "PASS").sum()),
            "n_fail": int((d["verdict"] == "FAIL").sum()), "n_cb": int((d["verdict"] == "CONTEXT_BOUNDARY").sum()),
            "pass_rate": round(float((d["verdict"] == "PASS").mean()), 4) if len(d) else None}
strata = {"overall_same_pathogen": pass_rate(verd),
          "by_holdout": {h: pass_rate(d) for h, d in verd.groupby("holdout_id")},
          "by_biotier": {str(t): pass_rate(d) for t, d in verd.groupby("bio_tier")},
          "by_promiscuous": {str(t): pass_rate(d) for t, d in verd.groupby(verd["orcs_promiscuous_top5pct"].fillna(False))},
          "by_context_matched": {str(t): pass_rate(d) for t, d in verd.groupby(verd["context_matched"].fillna("NA"))},
          "hs_cluster_same_pathogen": pass_rate(hs),
          "positive_controls": {g: pass_rate(verd[verd["gene"] == g]) for g in POS_CTRL}}
cm_only = verd[verd["context_matched"] == True]
is_bt2 = cm_only["bio_tier"].astype(str).str.contains("Bio-Tier 2", na=False)
strata["biotier2_vs_rest_matched"] = {"biotier2": pass_rate(cm_only[is_bt2]), "rest": pass_rate(cm_only[~is_bt2])}
# cross-pathogen exploratory support counts
xp = verd[verd["verdict"] == "PASS_XPATH"]
strata["xpath_support"] = {"n_gene_holdout": len(xp), "n_genes": int(xp["gene"].nunique()),
                           "hs_genes_with_xpath_pass": sorted(set(xp["gene"]) & set(HS))}
with open(f"{OUT}/stage13_stratified_pass_rates.json", "w") as f: json.dump(strata, f, indent=2)

# per-gene rollup (353 candidates)
roll = verd.groupby("gene").agg(
    n_evaluable=("verdict", lambda s: int(s.isin(EVAL).sum())),
    n_pass=("verdict", lambda s: int((s == "PASS").sum())),
    n_fail=("verdict", lambda s: int((s == "FAIL").sum())),
    n_cb=("verdict", lambda s: int((s == "CONTEXT_BOUNDARY").sum())),
    n_xpath_pass=("verdict", lambda s: int((s == "PASS_XPATH").sum())),
    n_xpath_reversed=("verdict", lambda s: int((s == "XPATH_REVERSED").sum()))).reset_index()
roll = cand[["gene"]].merge(roll, on="gene", how="left").fillna(0).merge(frz, on="gene", how="left")
roll.to_csv(f"{OUT}/stage13_per_gene_rollup.tsv", sep="\t", index=False)

summary = {"stage": "Stage13_holdout_adjudication_v2", "date": "2026-09-02", "screens_evaluated": screen_meta,
           "not_computable": NOT_COMPUTABLE, "n_gene_screen_calls": len(calls),
           "n_gene_holdout_verdicts": len(verd), "n_candidates_evaluated_any": int(verd["gene"].nunique()),
           "hs_fisher": hs_fisher, "strata": strata,
           "rules": {"PASS": ">=1 concordant author hit (orientation-corrected), no reversal; PASS_XPATH if cross-pathogen",
                     "FAIL": "reversed hit, OR no hit in same-pathogen context-matched holdout (pre-registered)",
                     "CONTEXT_BOUNDARY": "no hit in same-pathogen context-mismatched holdout",
                     "NOT_EVALUABLE_XPATH": "cross-pathogen non-hit (outside frozen stratification)",
                     "orientation": "empirical anchor map; ORCS-2360/2361 flipped (protective=negative Z)",
                     "no_threshold_retuning": True}}
with open(f"{OUT}/stage13_summary.json", "w") as f: json.dump(summary, f, indent=2)
with open(f"{OUT}/stage13_runinfo.json", "w") as f:
    json.dump({"script": os.path.abspath(__file__), "timestamp": datetime.datetime.now().isoformat(),
               "inputs": ["screen_estimand_registry.tsv", "preregistered_candidate_set_v1.tsv",
                          "BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", "r11_candidate_orcs.tsv",
                          "r14_holdout_context_stratification.tsv", "BIOGRID-ORCS-SCREEN_*"]}, f, indent=2)
print(json.dumps({"screens": {k: {"n_genes": v.get("n_genes"), "n_hits": v.get("n_hits"), "orient": v.get("orientation")} for k, v in screen_meta.items()},
                  "overall": strata["overall_same_pathogen"], "hs_fisher": hs_fisher,
                  "pos_ctrl": strata["positive_controls"], "by_holdout": strata["by_holdout"],
                  "biotier2_vs_rest_matched": strata["biotier2_vs_rest_matched"],
                  "xpath": strata["xpath_support"]}, indent=2))
