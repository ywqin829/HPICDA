#!/usr/bin/env python3
"""R4b: module-level convergence gate (post-hoc, uses saved B=10000 null).

Script 22 writes module Fisher p/q but no gate columns. This script:
  - recomputes observed per-unit module z deterministically (same code path)
  - loads the saved permutation null (convergence_v2_null.npz)
  - adds per-pathogen module p-values, per-cluster support counts, and the
    convergence gate: BH-FDR<0.05 (directions jointly) AND >=2 pathogens with
    nominal p<0.01 AND >=2 independence clusters with cluster-level nominal
    p<0.05 (module analog of the gene-level cluster rule; definition frozen here)
  - module LOSO stability over the 14 clusters
Output: 07_cross_pathogen/convergence_meta_v2_module_gated.tsv
        07_cross_pathogen/convergence_v2_module_loso.tsv
Usage: python3 22b_r4_module_gate.py   (run AFTER full 22_r4_convergence_v2.py)
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path("/workspace/HPICDA")

# ---------------- load (identical to script 22) ----------------
df = pd.read_csv(ROOT / "05_conditional_dependency/per_screen_calls_v2.tsv", sep="\t", dtype=str)
df["score"] = pd.to_numeric(df["score"], errors="coerce")
d = df[(df["pathogen"].isin(["Influenza", "SARS-CoV-2", "RSV"]))
       & (df["endpoint_family"].isin({"HOST_FITNESS", "PATHOGEN_BURDEN"}))
       & (df["call"] != "UNAVAILABLE")
       & (df["gene_class"].str.contains("protein-coding", na=False))
       & (df["gene_hgnc"].fillna("") != "")].copy()
d["x"] = np.where(d["direction_lof"] == "enriched", d["score"], -d["score"])

reg = pd.read_csv(ROOT / "03_functional_screens/screen_estimand_registry.tsv", sep="\t", dtype=str)
meta = reg.drop_duplicates("screen_id").set_index("screen_id")
units = sorted(d["unit_id"].unique())
unit_cluster = {u: meta["independence_cluster"].get(d.loc[d["unit_id"] == u, "screen_id"].iloc[0], "")
                for u in units}
unit_pathogen = {u: d.loc[d["unit_id"] == u, "pathogen"].iloc[0] for u in units}
clusters = sorted(set(unit_cluster.values()))
pathogens = sorted(set(unit_pathogen.values()))
cl_units = {c: [i for i, u in enumerate(units) if unit_cluster[u] == c] for c in clusters}
pa_clusters = {p: [c for c in clusters if any(unit_pathogen[u] == p and unit_cluster[u] == c
                                              for u in units)] for p in pathogens}
U = sorted(d["gene_hgnc"].unique())
gidx = {g: i for i, g in enumerate(U)}
nU = len(U)

unit_x = np.full((len(units), nU), np.nan)
unit_present = np.zeros((len(units), nU), dtype=bool)
for i, u in enumerate(units):
    g = d[d["unit_id"] == u]
    ii = g["gene_hgnc"].map(gidx).to_numpy()
    unit_x[i, ii] = g["x"].to_numpy()
    unit_present[i, ii] = True

hgnc = pd.read_csv(ROOT / "02_metadata/hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False)
entrez2sym = dict(zip(hgnc["entrez_id"], hgnc["symbol"]))
modules = {}
with open(ROOT / "02_metadata/NCBI2Reactome_All_Levels.txt") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6 or p[5] != "Homo sapiens":
            continue
        sym = entrez2sym.get(p[0])
        if sym and sym in gidx:
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
    modules[k] = v & set(U)
modules = {k: sorted(v) for k, v in modules.items() if 5 <= len(v) <= 500}
mod_names = sorted(modules)
mod_idx = [np.array([gidx[g] for g in modules[m]]) for m in mod_names]
nM = len(mod_names)
print(f"units={len(units)} clusters={len(clusters)} universe={nU} modules={nM}", flush=True)

# ---------------- observed stats (deterministic) ----------------
def unit_z(xmat):
    z = np.zeros_like(xmat)
    for i in range(xmat.shape[0]):
        pres = unit_present[i]
        vals = xmat[i][pres]
        order = np.argsort(-vals, kind="mergesort")
        ranks = np.empty(len(vals))
        ranks[order] = np.arange(1, len(vals) + 1)
        p = ranks / (len(vals) + 1.0)
        z[i][pres] = norm.ppf(1.0 - np.clip(p, 1e-16, 1 - 1e-16))
    return z

def module_unit_z(zmat):
    """competitive diff-in-means per unit -> (n_units, nM) module z."""
    mz = np.zeros((len(units), nM))
    for i in range(len(units)):
        pres = unit_present[i]
        zp = zmat[i][pres]
        tot, cnt = zp.sum(), pres.sum()
        for mi, idx in enumerate(mod_idx):
            m = idx[pres[idx]]
            if len(m) < 5:
                continue
            mz[i, mi] = (zmat[i][m].sum() - (tot - zmat[i][m].sum()) / max(cnt - len(m), 1) * len(m)) \
                        / np.sqrt(len(m))
    return mz

def combine_full(mz, drop_cluster=None):
    """(n_units, nM) -> fisher, per-pathogen p, per-cluster z."""
    zc = np.zeros((len(clusters), mz.shape[1]))
    for ci, c in enumerate(clusters):
        if c == drop_cluster:
            continue
        iu = cl_units[c]
        zc[ci] = mz[iu].sum(axis=0) / np.sqrt(len(iu))
    pp = np.ones((len(pathogens), mz.shape[1]))
    for pi, p in enumerate(pathogens):
        cs = [clusters.index(c) for c in pa_clusters[p] if c != drop_cluster]
        if not cs:
            continue
        zp = zc[cs].sum(axis=0) / np.sqrt(len(cs))
        pp[pi] = 1.0 - norm.cdf(zp)
    pp = np.clip(pp, 1e-300, 1.0)
    fisher = -2.0 * np.log(pp).sum(axis=0)
    return fisher, pp, zc

def bh_fast(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    ranked = p[o] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[o] = np.clip(ranked, 0, 1)
    return q

z_obs = unit_z(unit_x)
mz_obs = module_unit_z(z_obs)
fE, ppE, zcE = combine_full(mz_obs)
fP, ppP, zcP = combine_full(-mz_obs)

# ---------------- empirical p from saved null ----------------
null = np.load(ROOT / "07_cross_pathogen/convergence_v2_null.npz")
srt_m = np.sort(null["module"].ravel())[::-1]
B = null["module"].shape[0]
print(f"loaded null: B={B}", flush=True)

def emp_p(obs, srt):
    c = np.searchsorted(-srt, -obs)
    return (c + 1) / (srt.size + 1)

pE, pP = emp_p(fE, srt_m), emp_p(fP, srt_m)
qb = bh_fast(np.concatenate([pE, pP]))
qE, qP = qb[:nM], qb[nM:]

# cluster support: cluster-level one-sided nominal p<0.05
cl_pE = 1.0 - norm.cdf(zcE)   # (n_clusters, nM)
cl_pP = 1.0 - norm.cdf(zcP)
nclE = (cl_pE < 0.05).sum(axis=0)
nclP = (cl_pP < 0.05).sum(axis=0)
npatE = (ppE < 0.01).sum(axis=0)
npatP = (ppP < 0.01).sum(axis=0)

convE = (qE < 0.05) & (npatE >= 2) & (nclE >= 2)
convP = (qP < 0.05) & (npatP >= 2) & (nclP >= 2)

base = pd.read_csv(ROOT / "07_cross_pathogen/convergence_meta_v2_module.tsv", sep="\t")
base = base.set_index("module")
if "frac_members_essential" not in base.columns:  # pilot-null compatibility
    base["frac_members_essential"] = np.nan
out = pd.DataFrame({
    "module": mod_names, "n_genes": [len(m) for m in mod_idx],
    "p_EXP": pE, "q_EXP": qE, "p_PRO": pP, "q_PRO": qP,
    "n_pathogens_p01_EXP": npatE, "n_pathogens_p01_PRO": npatP,
    "n_clusters_p05_EXP": nclE, "n_clusters_p05_PRO": nclP,
    "convergent_EXP": convE, "convergent_PRO": convP,
    "p_SC2_EXP": ppE[pathogens.index("SARS-CoV-2")],
    "p_IAV_EXP": ppE[pathogens.index("Influenza")],
    "p_RSV_EXP": ppE[pathogens.index("RSV")],
    "p_SC2_PRO": ppP[pathogens.index("SARS-CoV-2")],
    "p_IAV_PRO": ppP[pathogens.index("Influenza")],
    "p_RSV_PRO": ppP[pathogens.index("RSV")],
    "frac_members_essential": [base["frac_members_essential"].get(m, np.nan) for m in mod_names],
}).sort_values("q_EXP")
out.to_csv(ROOT / "07_cross_pathogen/convergence_meta_v2_module_gated.tsv", sep="\t", index=False)

# ---------------- module LOSO ----------------
rows = []
fdr_pass = np.where(convE | convP)[0]
for c in clusters:
    fEc, _, _ = combine_full(mz_obs, drop_cluster=c)
    fPc, _, _ = combine_full(-mz_obs, drop_cluster=c)
    pEc, pPc = emp_p(fEc, srt_m), emp_p(fPc, srt_m)
    qc = bh_fast(np.concatenate([pEc, pPc]))
    for mi in fdr_pass:
        was = "EXP" if convE[mi] else "PRO"
        still = (qc[:nM] < 0.05)[mi] if was == "EXP" else (qc[nM:] < 0.05)[mi]
        rows.append({"module": mod_names[mi], "direction": was,
                     "dropped_cluster": c, "still_fdr_pass": bool(still)})
loso = pd.DataFrame(rows)
if len(loso):
    loso_sum = loso.groupby(["module", "direction"])["still_fdr_pass"].agg(["sum", "count", "mean"]).reset_index()
    loso_sum.columns = ["module", "direction", "n_pass", "n_runs", "stability"]
else:
    loso_sum = pd.DataFrame(columns=["module", "direction", "n_pass", "n_runs", "stability"])
loso_sum.to_csv(ROOT / "07_cross_pathogen/convergence_v2_module_loso.tsv", sep="\t", index=False)

print(f"\nconvergent modules EXP={convE.sum()} PRO={convP.sum()}")
print(out[convE | convP].head(40).to_string(index=False))
print("\ncurated modules:")
print(out[out["module"].str.startswith("CURATED")].to_string(index=False))
json.dump({"B": int(B), "gate": "BH-FDR<0.05 AND >=2 pathogens p<0.01 AND >=2 clusters p<0.05",
           "n_convergent_EXP": int(convE.sum()), "n_convergent_PRO": int(convP.sum())},
          open(ROOT / "07_cross_pathogen/convergence_v2_module_gate_runinfo.json", "w"), indent=2)
