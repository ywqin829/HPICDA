#!/usr/bin/env python3
"""Stage 12: biological candidate freeze (blueprint §31), post-HG3-approval.

Freeze basis ONLY: functional conditional dependency, cross-pathogen replication,
baseline fitness, human pulmonary relevance, human genetics, mechanistic plausibility.
ChEMBL/DrugCentral/druggability results are NOT read for gene selection (compliance
recorded in runinfo). Druggability (Stage 13/14) stays locked until
HUMAN_GATE_BIOLOGICAL_FREEZE.

Bio-Tier mapping (documented, derived from approved HG3 scoring + blueprint §32):
  Bio-Tier 1: cross-class functional -> EMPTY (viral-only discovery corpus;
              cross-class convergence not computable; breadth label recorded)
  Bio-Tier 2: HG3 Tier1 v2 (12)  - viral cross-pathogen functional, human context,
              direction resolved, no contradiction
  Bio-Tier 3: HG3 Tier2 v2 - single-pathogen / weaker functional, reasonable
              human context
  ARCHIVE:    HG3 watchlist v2
  REJECT:     HG3 excluded v2 with reject_reason

Outputs (09_candidate_freeze/):
  BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv      (all 353 scored candidates + disposition)
  candidate_direction_adjudication_v2.tsv  (frozen candidates; blueprint §31.1 fields)
  candidate_contradiction_matrix_v2.tsv    (353 x per-pathogen consensus direction)
  holdout_manifest_v2.tsv                  (identity/role only; unlock requires HG4)
  stage12_freeze_v2_runinfo.json           (sha256, timestamp, git commit, protocol hash)
"""
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "09_candidate_freeze"
OUT.mkdir(exist_ok=True)
TS = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

score = pd.read_csv(ROOT / "18_human_reviews/HG3_candidate_scoring_v2.tsv", sep="\t")
rep = pd.read_csv(ROOT / "06_within_pathogen/within_pathogen_replication_v2.tsv", sep="\t")
conv = pd.read_csv(ROOT / "07_cross_pathogen/convergence_meta_v2_gene.tsv", sep="\t").set_index("gene")
hints = pd.read_csv(ROOT / "11_modules/stage11_module_crossclass_hints.tsv", sep="\t")
reg = pd.read_csv(ROOT / "03_functional_screens/screen_estimand_registry.tsv", sep="\t", dtype=str)
hold = pd.read_csv(ROOT / "01_stage0_audit/holdout_inventory.tsv", sep="\t")

# ---- convergent-module membership (rebuild identically to scripts 22/22b/28) ----
hgnc = pd.read_csv(ROOT / "02_metadata/hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False)
entrez2sym = dict(zip(hgnc["entrez_id"], hgnc["symbol"]))
calls = pd.read_csv(ROOT / "05_conditional_dependency/per_screen_calls_v2.tsv", sep="\t", dtype=str)
U = set(calls.loc[calls["gene_class"].str.contains("protein-coding", na=False)
                  & (calls["gene_hgnc"].fillna("") != ""), "gene_hgnc"])
modules = {}
with open(ROOT / "02_metadata/NCBI2Reactome_All_Levels.txt") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6 or p[5] != "Homo sapiens":
            continue
        sym = entrez2sym.get(p[0])
        if sym and sym in U:
            modules.setdefault(f"Reactome:{p[3]}", set()).add(sym)
CURATED = {
    "CURATED:V-ATPase proton pump": {"ATP6V1A","ATP6V1B1","ATP6V1B2","ATP6V1C1","ATP6V1D","ATP6V1E1",
        "ATP6V1F","ATP6V1G1","ATP6V1G2","ATP6V1H","ATP6V0A1","ATP6V0A2","ATP6V0A4","ATP6V0B",
        "ATP6V0C","ATP6V0D1","ATP6V0D2","ATP6V0E1","ATP6V0E2","ATP6AP1","ATP6AP2","VMA21","TCIRG1"},
    "CURATED:COG complex": {"COG1","COG2","COG3","COG4","COG5","COG6","COG7","COG8"},
    "CURATED:heparan-sulfate biosynthesis": {"EXT1","EXT2","EXTL1","EXTL2","EXTL3","B3GALT6","B4GALT7",
        "B3GAT3","XYLT1","XYLT2","SLC35B2","NDST1","NDST2","HS2ST1","HS3ST1","HS3ST6","HS6ST1","HS6ST2","GLCE"},
}
for k, v in CURATED.items():
    modules[k] = v & U
modules = {k: sorted(v) for k, v in modules.items() if 5 <= len(v) <= 500}

mgate = pd.read_csv(ROOT / "07_cross_pathogen/convergence_meta_v2_module_gated.tsv", sep="\t")
conv_mods = set(mgate.loc[mgate["convergent_EXP"] | mgate["convergent_PRO"], "module"])
gene2convmods = {}
for m in conv_mods:
    for g in modules.get(m, []):
        gene2convmods.setdefault(g, []).append(m)

# cross-class HINT carriers (gene carries exploratory bacterial signal in a convergent module)
gene2hint = {}
for _, r in hints[hints["cross_class_hint"]].iterrows():
    for g in str(r["bacterial_member_signals"]).split(";"):
        if g and g != "nan":
            gene2hint.setdefault(g, []).append(r["module"])

# ---- registry lookups ----
reg_ep = dict(zip(reg["screen_id"], reg["endpoint_family"]))
def endpoint_of_units(units):
    eps = set()
    for u in str(units).split(";"):
        sid = u.split("|")[-1]
        ep = reg_ep.get(sid) or reg_ep.get(u)
        if ep:
            eps.add(ep)
    return ";".join(sorted(eps)) if eps else "unknown"

# ---- Bio-Tier assignment ----
BIO = {"Tier1": "Bio-Tier 2", "Tier2": "Bio-Tier 3", "watchlist": "ARCHIVE", "excluded": "REJECT"}
rows = []
for _, r in score.iterrows():
    g = r["gene"]
    bt = BIO[r["hg3_tier"]]
    reason = ""
    if bt == "REJECT":
        if r["contradictory_any"]:
            reason = "intervention direction contradictory (WP_CONTRADICTORY record; total<6)"
        else:
            reason = "insufficient evidence (total<2)"
        if r["therapy_window_penalty"] <= -3:
            reason += "; baseline fitness partially explains (min Chronos<-2)"
    breadth = "WITHIN_CLASS_VIRAL_ONLY (cross-class not computable: viral-only discovery corpus)"
    if g in gene2hint:
        breadth += "; cross-class HINT carrier in " + ";".join(gene2hint[g]) + " (exploratory bacterial signal, HINT level only)"
    rows.append({
        "gene": g, "intervention_direction": r["intervention_direction"],
        "bio_tier": bt, "cross_class_breadth_label": breadth,
        "hg3_tier": r["hg3_tier"], "total_score": r["total_score"],
        "screen_evidence_pts": r["screen_evidence_pts"],
        "best_replication_tier": r["best_replication_tier"],
        "n_independence_clusters_max": r["n_independence_clusters_max"],
        "pathogens": r["pathogens"], "sources": r["sources"],
        "convergence_bonus": r["convergence_bonus"], "loso_stability": r["loso_stability"],
        "min_baseline_chronos": r["min_baseline_chronos"],
        "frac_backgrounds_essential": r["frac_backgrounds_essential"],
        "therapy_window_penalty": r["therapy_window_penalty"],
        "lung_expressed": r["lung_expressed"], "clinical_de_concordant": r["clinical_de_concordant"],
        "genetics_tier": r["genetics_tier"], "loeuf": r["loeuf"],
        "contradictory_any": r["contradictory_any"],
        "convergent_module_memberships": ";".join(sorted(gene2convmods.get(g, []))),
        "reject_reason": reason, "freeze_timestamp": TS,
    })
freeze = pd.DataFrame(rows)
freeze.to_csv(OUT / "BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv", sep="\t", index=False)

# ---- direction adjudication (frozen candidates only) ----
rep_g = rep.groupby("gene")
def direction_confidence(r):
    if r["screen_evidence_pts"] >= 3 and r["convergence_bonus"] >= 2 and \
       (pd.isna(r["loso_stability"]) or r["loso_stability"] >= 0.9):
        return "high"
    if r["screen_evidence_pts"] >= 3 or r["convergence_bonus"] >= 2:
        return "medium"
    return "low"

adj_rows = []
for _, r in freeze[freeze["bio_tier"].isin(["Bio-Tier 2", "Bio-Tier 3"])].iterrows():
    g = r["gene"]
    gexp = r["intervention_direction"] == "INHIBIT"
    gre = rep_g.get_group(g) if g in rep_g.groups else None
    units = ";".join(sorted(set(";".join(gre["supporting_units"].dropna()).split(";")))) if gre is not None else ""
    ko_eff = ("KO/LOF enriched under infection selection (host-protective KO)" if gexp
              else "KO/LOF depleted under infection selection (host-harming KO)")
    pat_eff = ("gene KO reduces pathogen fitness/burden proxy -> pathogen exploits gene"
               if gexp else "gene KO increases pathogen fitness/burden proxy -> gene restricts pathogen")
    mc = r["min_baseline_chronos"]
    host_eff = (f"min baseline Chronos {mc:.2f} across pulmonary-relevant backgrounds; "
                f"frac essential {r['frac_backgrounds_essential']:.2f}") if pd.notna(mc) else \
               "baseline Chronos unavailable (unknown, not zero-risk)"
    ev = [f"screen {r['screen_evidence_pts']:.2f} pts (best {r['best_replication_tier']}, "
          f"{r['n_independence_clusters_max']} clusters)"]
    if r["convergence_bonus"] >= 2:
        q = conv.loc[g, "q_EXP" if gexp else "q_PRO"] if g in conv.index else np.nan
        ev.append(f"meta-convergent q={q:.2e}" + (f", LOSO {r['loso_stability']:.2f}"
                  if pd.notna(r["loso_stability"]) else ""))
    elif r["convergence_bonus"] == 1:
        ev.append("threshold-intersection candidate (pre-FDR vocabulary)")
    if r["clinical_de_concordant"]:
        ev.append("clinical-DE concordant (Stage 8)")
    if r["lung_expressed"]:
        ev.append("lung-expressed (Stage 9)")
    if r["genetics_tier"] in ("G1", "G2"):
        ev.append(f"human genetics {r['genetics_tier']} (Stage 10)")
    adj_rows.append({
        "gene": g, "screen_id": units,
        "KO_or_perturbation_effect": ko_eff, "pathogen_effect": pat_eff,
        "host_effect": host_eff, "endpoint_family": endpoint_of_units(units),
        "baseline_fitness": mc if pd.notna(mc) else "unknown",
        "required_intervention": "INHIBIT" if gexp else "PRESERVE_OR_ACTIVATE",
        "evidence_supporting_direction": "; ".join(ev),
        "contradictory_evidence": "none (no WP_CONTRADICTORY record)",
        "direction_confidence": direction_confidence(r),
        "adjudication_status": "RESOLVED",
        "adjudicator": "Biomni agent per HG3_scoring_rubric_v1 + blueprint §31.1 (HG3 approved 2026-09-01)",
        "freeze_timestamp": TS,
    })
adj = pd.DataFrame(adj_rows)
adj.to_csv(OUT / "candidate_direction_adjudication_v2.tsv", sep="\t", index=False)

# ---- contradiction matrix (all 353 candidates x pathogen) ----
PATHOGENS = ["SARS-CoV-2", "Influenza", "RSV", "Mtb", "BCG", "Pseudomonas"]
cm_rows = []
for g in score["gene"]:
    row = {"gene": g}
    ncontra = 0
    gre = rep_g.get_group(g) if g in rep_g.groups else None
    for p in PATHOGENS:
        cell = "no_evidence"
        if gre is not None:
            sub = gre[gre["pathogen"] == p]
            if len(sub):
                s = sub.iloc[0]
                cell = f"{s['wp_status']}:{s['consensus_direction']}"
                if s["wp_status"] == "WP_CONTRADICTORY":
                    ncontra += 1
        row[p] = cell
    row["n_contradictory_pathogens"] = ncontra
    cm_rows.append(row)
cm = pd.DataFrame(cm_rows)
cm.to_csv(OUT / "candidate_contradiction_matrix_v2.tsv", sep="\t", index=False)

# ---- holdout manifest (identity/role only) ----
hold_out = hold.copy()
hold_out["manifest_note"] = ("identity/role only; no candidate outcomes; "
                             "unlock requires HUMAN_GATE_HOLDOUT_UNLOCK (HG4) per blueprint §33.1")
hold_out.to_csv(OUT / "holdout_manifest_v2.tsv", sep="\t", index=False)

# ---- git commit (freeze provenance) ----
gitignore = ROOT / ".gitignore"
if not gitignore.exists():
    gitignore.write_text("*.csv\n*.npz\n*.npy\n*.gz\n*.log\n__pycache__/\n"
                         "09_transcriptomic_context/raw/\n"
                         "05_conditional_dependency/per_screen_calls.tsv\n"
                         "05_conditional_dependency/per_screen_calls_v2.tsv\n")
def git(*args):
    return subprocess.run(["git", "-C", str(ROOT)] + list(args),
                          capture_output=True, text=True)
if not (ROOT / ".git").exists():
    git("init", "-q")
    git("config", "user.email", "biomni-agent@phylo.local")
    git("config", "user.name", "Biomni agent")
git("add", "-A")
cmts = git("commit", "-q", "-m", f"Stage 12: biological candidate freeze (HG3 approved {TS})")
commit = git("rev-parse", "HEAD").stdout.strip() or "unavailable"

# ---- runinfo ----
runinfo = {
    "stage": "12", "freeze_timestamp_utc": TS,
    "git_commit": commit,
    "protocol_sha256": sha(ROOT / "00_protocol/MASTER_EXECUTION_BLUEPRINT_v1.1.md"),
    "hg3_approval": "18_human_reviews/HG3.approved (2026-09-01)",
    "freeze_basis_only": ["functional conditional dependency", "cross-pathogen replication",
                          "baseline fitness", "human pulmonary relevance", "human genetics",
                          "mechanistic plausibility"],
    "druggability_compliance": "ChEMBL/DrugCentral/druggability results NOT read for gene selection; Stage 13/14 locked until HUMAN_GATE_BIOLOGICAL_FREEZE",
    "bio_tier_rule": {"Bio-Tier 1": "cross-class functional - EMPTY (viral-only discovery corpus; cross-class not computable)",
                      "Bio-Tier 2": "HG3 Tier1", "Bio-Tier 3": "HG3 Tier2",
                      "ARCHIVE": "HG3 watchlist", "REJECT": "HG3 excluded"},
    "bio_tier_counts": freeze["bio_tier"].value_counts().to_dict(),
    "direction_confidence_counts": adj["direction_confidence"].value_counts().to_dict(),
    "output_sha256": {p.name: sha(p) for p in sorted(OUT.glob("*.tsv"))},
    "input_sha256": {p: sha(ROOT / p) for p in [
        "18_human_reviews/HG3_candidate_scoring_v2.tsv",
        "18_human_reviews/HG3_scoring_rubric_v1.yaml",
        "06_within_pathogen/within_pathogen_replication_v2.tsv",
        "07_cross_pathogen/convergence_meta_v2_gene.tsv",
        "07_cross_pathogen/convergence_meta_v2_module_gated.tsv",
        "11_modules/stage11_module_crossclass_hints.tsv",
        "03_functional_screens/screen_estimand_registry.tsv",
        "01_stage0_audit/holdout_inventory.tsv"]},
}
json.dump(runinfo, open(OUT / "stage12_freeze_v2_runinfo.json", "w"), indent=2)

print("freeze timestamp:", TS)
print("git commit:", commit)
print("bio_tier counts:", runinfo["bio_tier_counts"])
print("direction confidence:", runinfo["direction_confidence_counts"])
print("\nBio-Tier 2 (frozen primary candidates):")
print(freeze[freeze["bio_tier"] == "Bio-Tier 2"][["gene", "intervention_direction",
      "pathogens", "total_score"]].to_string(index=False))
print("\nwrote", OUT)
