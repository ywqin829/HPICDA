#!/usr/bin/env python3
"""Stage 11: module-level convergence synthesis + contradictory-pair interpretation.

Components (plan v2):
  A. Gated module convergence table (from R4b) + cross-class attempt:
     all discovery screens are viral -> screen-based cross-class convergence is
     NOT computable; we record cross-class HINTS = convergent viral module with
     >=1 member gene carrying exploratory Mtb/BCG signal (PATHOGEN_STATE/aux)
     or Mtb expression concordance (Stage 8). Labeled as hint, never convergence.
  B. WP_CONTRADICTORY interpretation: stratify by cell-line overlap, cluster
     pattern, direction split, baseline essentiality.
  C. DepMap co-essentiality of convergent modules (independent validation):
     mean pairwise Pearson r of Chronos profiles among members vs size-matched
     random modules (10k draws).
  D. IntAct physical-interaction support (fail-soft): edges among members vs
     size-matched random sets.
Usage: python3 28_stage11_modules_contradictions.py
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "11_modules"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20260831)

# ---------------- A. module convergence + cross-class hints ----------------
mod = pd.read_csv(ROOT / "07_cross_pathogen/convergence_meta_v2_module_gated.tsv", sep="\t")
conv = mod[mod["convergent_EXP"] | mod["convergent_PRO"]].copy()
print(f"gated convergent modules: {len(conv)}")

# module membership (rebuild identically to script 22/22b)
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

# exploratory bacterial signal per gene (Mtb PATHOGEN_STATE / BCG aux / Pa context)
bac = calls[calls["pathogen"].isin(["Mtb", "BCG", "Pseudomonas"])]
bac_sig = set(bac.loc[pd.to_numeric(bac["score"], errors="coerce") >= 2.0, "gene_hgnc"])
hint_rows = []
for _, r in conv.iterrows():
    mem = set(modules.get(r["module"], []))
    hits = sorted(mem & bac_sig)
    hint_rows.append({"module": r["module"], "direction": "EXP" if r["convergent_EXP"] else "PRO",
                      "n_genes": r["n_genes"], "q_EXP": r["q_EXP"], "q_PRO": r["q_PRO"],
                      "frac_members_essential": r["frac_members_essential"],
                      "bacterial_member_signals": ";".join(hits),
                      "cross_class_hint": len(hits) > 0})
hints = pd.DataFrame(hint_rows)
# specificity context (descriptive only; frozen hint definition unchanged):
# base rate of exploratory bacterial signal among all genes / all modules
bac_in_U = bac_sig & U
gene_base = len(bac_in_U) / len(U)
mod_base = sum(1 for m in modules.values() if set(m) & bac_sig) / len(modules)
hints["frac_members_bac_signal"] = [
    len(set(modules.get(r["module"], [])) & bac_sig) / max(len(modules.get(r["module"], [])), 1)
    for _, r in conv.iterrows()]
hints.to_csv(OUT / "stage11_module_crossclass_hints.tsv", sep="\t", index=False)
print(f"cross-class hints: {int(hints['cross_class_hint'].sum())} modules "
      f"(HINT level only - no bacterial discovery screen exists)")
print(f"  specificity context: bac-signal genes {len(bac_in_U)}/{len(U)} ({gene_base:.3f}); "
      f"fraction of ALL modules with >=1 bac member: {mod_base:.3f}")

# ---------------- B. WP_CONTRADICTORY interpretation ----------------
rep = pd.read_csv(ROOT / "06_within_pathogen/within_pathogen_replication_v2.tsv", sep="\t")
contra = rep[rep["wp_status"] == "WP_CONTRADICTORY"].copy()
bl = pd.read_csv(ROOT / "04_depmap_baseline/baseline_chronos_by_background.tsv", sep="\t").set_index("gene")
rel = [c for c in ["A-549", "Calu-3", "HuH-7", "Huh-7.5", "Huh-7.5.1", "HEp-2", "HeLa"] if c in bl.columns]
contra["min_chronos"] = [bl.loc[g, rel].min() if g in bl.index else np.nan for g in contra["gene"]]
contra["same_cell_line"] = contra["cell_lines"].fillna("").apply(
    lambda s: len(set(s.split(";"))) == 1)
contra["pattern"] = np.where(
    (contra["n_clusters_exploited"] > 0) & (contra["n_clusters_protective"] > 0),
    "EXP_in_some_PRO_in_others", "other")
contra["essentiality_confound_suspect"] = contra["min_chronos"] < -1
summ = (contra.groupby(["pathogen", "pattern", "same_cell_line", "essentiality_confound_suspect"])
        .size().reset_index(name="n_pairs"))
summ.to_csv(OUT / "stage11_contradictory_pairs_stratified.tsv", sep="\t", index=False)
contra.to_csv(OUT / "stage11_contradictory_pairs_annotated.tsv", sep="\t", index=False)
print(f"\nWP_CONTRADICTORY pairs: {len(contra)}; stratification:")
print(summ.to_string(index=False))

# ---------------- C. DepMap co-essentiality ----------------
print("\nloading DepMap 26Q1 gene effect matrix...")
ge = pd.read_csv(ROOT / "04_depmap_baseline/gene_effect_26Q1.csv", index_col=0)
ge.columns = [c.split(" ")[0] for c in ge.columns]
ge = ge.T.groupby(level=0).mean().T  # dedupe columns if any
# NOTE: no global dropna — with 18.5k genes every cell line has >=1 NA.
# Per-module listwise deletion is done inside coess_fast instead.
print(f"  matrix: {ge.shape} (per-module listwise NA deletion)")

pool = sorted(g for g in U if g in ge.columns)
pidx = {g: i for i, g in enumerate(pool)}
G = ge[pool].to_numpy(dtype=np.float32)  # cell lines x pool genes
print(f"  pool genes in DepMap: {len(pool)}")

def coess_fast(member_idx):
    """Mean pairwise Pearson r of Chronos profiles (listwise-complete cell lines)."""
    if len(member_idx) < 5:
        return np.nan
    sub = G[:, member_idx]
    sub = sub[~np.isnan(sub).any(axis=1)]
    if sub.shape[0] < 50:
        return np.nan
    r = np.corrcoef(sub.T.astype(np.float64))
    iu = np.triu_indices(len(member_idx), 1)
    return float(np.nanmean(r[iu]))

target_modules = sorted(set(conv["module"]) | set(CURATED))
co_rows = []
for i, name in enumerate(target_modules):
    mem = modules.get(name, [])
    midx = [pidx[g] for g in mem if g in pidx]
    r_obs = coess_fast(midx)
    if np.isnan(r_obs):
        continue
    # size-matched background (2,000 draws)
    k = min(len(mem), len(pool))
    bg = np.array([coess_fast(rng.choice(len(pool), size=k, replace=False))
                   for _ in range(2000)])
    bg = bg[~np.isnan(bg)]
    p_emp = (np.sum(bg >= r_obs) + 1) / (len(bg) + 1)
    co_rows.append({"module": name, "n_members": len(mem), "n_in_depmap": len(midx),
                    "mean_pairwise_r": r_obs, "bg_mean": float(np.mean(bg)),
                    "bg_sd": float(np.std(bg)), "p_empirical": p_emp,
                    "convergent": name in set(conv["module"])})
    if (i + 1) % 50 == 0:
        print(f"  co-essentiality: {i + 1}/{len(target_modules)} modules done", flush=True)
co = pd.DataFrame(co_rows)
co.to_csv(OUT / "stage11_module_coessentiality.tsv", sep="\t", index=False)
print(co.head(30).to_string(index=False))

# ---------------- D. IntAct interaction support (fail-soft, best-effort) ----------------
import requests
def intact_intra_edges(members, tries=2):
    """Count unique intra-module interactions via IntAct findInteractions per member.
    Partner identified by regex on member symbols in interaction record text.
    Best-effort (text match); fail-soft returns None."""
    memset = set(members)
    edges = set()
    for g in members:
        j = None
        for t in range(tries):
            try:
                r = requests.get(f"https://www.ebi.ac.uk/intact/ws/interaction/findInteractions/{g}",
                                 params={"page": 0, "pageSize": 100}, timeout=20)
                if r.status_code == 200:
                    j = r.json()
                    break
            except Exception:
                pass
            time.sleep(2)
        if j is None:
            return None
        for it in j.get("content", []):
            txt = json.dumps(it)
            for other in memset - {g}:
                if f'"{other}"' in txt or f" {other} " in txt or f"({other})" in txt:
                    edges.add(tuple(sorted((g, other))))
        time.sleep(0.15)
    return edges

ia_rows = []
intact_ok = True
for name in [m for m in target_modules if m.startswith("CURATED")]:
    mem = modules.get(name, [])
    edges = intact_intra_edges(mem)
    if edges is None:
        intact_ok = False
        ia_rows.append({"module": name, "status": "intact_unreachable", "n_members": len(mem)})
        continue
    ia_rows.append({"module": name, "status": "ok", "n_members": len(mem),
                    "n_intra_module_edges": len(edges),
                    "edges": ";".join(f"{a}-{b}" for a, b in sorted(edges))})
ia = pd.DataFrame(ia_rows)
ia.to_csv(OUT / "stage11_intact_support.tsv", sep="\t", index=False)
print("\nIntAct:", ia[["module", "status", "n_intra_module_edges"]].to_string(index=False))

meta = {"convergent_modules": len(conv),
        "cross_class_hints": int(hints["cross_class_hint"].sum()),
        "bac_gene_base_rate": round(gene_base, 4),
        "bac_module_base_rate": round(mod_base, 4),
        "contradictory_pairs": len(contra),
        "coessentiality_bg_draws": 2000,
        "coessentiality_modules_tested": len(co),
        "intact_ok": intact_ok,
        "note": "cross-class = HINT level only; no bacterial discovery screen exists"}
json.dump(meta, open(OUT / "stage11_runinfo.json", "w"), indent=2)
print("\nwrote", OUT)
