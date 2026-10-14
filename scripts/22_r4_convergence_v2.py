#!/usr/bin/env python3
"""R4: cross-pathogen convergence statistics v2 (gene + module level).

Design (plan v2 R4, "standard" strictness):
  universe : viral units (SC2/IAV/RSV), HOST_FITNESS/PATHOGEN_BURDEN endpoints,
             UNAVAILABLE excluded, protein-coding HGNC symbols only
  statistic: signed score x = +score if direction_lof==enriched else -score
             (enriched = KO protects = EXP direction; CRISPRa sign-inverted at parse)
  combine  : within-unit one-sided rank p -> cluster Stouffer -> pathogen Stouffer
             -> cross-pathogen Fisher (df=2k)
  null     : B stratified permutations (within-unit, within essentiality decile of
             mapped background; library size preserved by within-unit permutation;
             expression stratum unavailable -> AM2-13 deviation); pooled null across
             genes for p resolution; PRO direction reuses EXP null (sign symmetry)
  calls    : convergence requires BH-FDR<0.05 (both directions corrected together)
             AND >=2 pathogens with pathogen-level nominal p<0.01
             AND >=2 independence clusters total
  LOSO     : leave-one-cluster-out; per FDR-passing gene, fraction of runs staying
             FDR-passing
  modules  : Reactome (5-500 universe genes) + 3 curated complex modules
             (V-ATPase, COG, HS/GAG biosynthesis); competitive diff-in-means of
             unit-level z; same combination + permutation null; BH-FDR
Outputs:
  07_cross_pathogen/convergence_meta_v2_gene.tsv
  07_cross_pathogen/convergence_meta_v2_module.tsv
  07_cross_pathogen/convergence_v2_loso.tsv
  07_cross_pathogen/convergence_v2_runinfo.json
  07_cross_pathogen/convergence_v2_null.npz
Usage: python3 22_r4_convergence_v2.py [B]
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path("/workspace/HPICDA")
B_TOTAL = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
BLOCK = 100
rng = np.random.default_rng(20260831)

# ---------------- load ----------------
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
unit_cell = {u: d.loc[d["unit_id"] == u, "cell_line"].iloc[0] for u in units}
clusters = sorted(set(unit_cluster.values()))
pathogens = sorted(set(unit_pathogen.values()))
cl_units = {c: [i for i, u in enumerate(units) if unit_cluster[u] == c] for c in clusters}
pa_clusters = {p: [c for c in clusters if any(unit_pathogen[u] == p and unit_cluster[u] == c
                                              for u in units)] for p in pathogens}
print(f"units={len(units)} clusters={len(clusters)} pathogens={pathogens}", flush=True)

U = sorted(d["gene_hgnc"].unique())
gidx = {g: i for i, g in enumerate(U)}
nU = len(U)
print("universe:", nU, flush=True)

bl = pd.read_csv(ROOT / "04_depmap_baseline/baseline_chronos_by_background.tsv", sep="\t").set_index("gene")

def strata_for(unit):
    cell = unit_cell[unit]
    if cell in bl.columns:
        ch = pd.Series([bl[cell].get(g, np.nan) for g in U]).fillna(0.0)
        return pd.qcut(ch.rank(method="first"), 10, labels=False).to_numpy()
    return np.zeros(nU, dtype=int)

unit_x = np.full((len(units), nU), np.nan)
unit_present = np.zeros((len(units), nU), dtype=bool)
unit_strata = np.zeros((len(units), nU), dtype=int)
for i, u in enumerate(units):
    g = d[d["unit_id"] == u]
    ii = g["gene_hgnc"].map(gidx).to_numpy()
    unit_x[i, ii] = g["x"].to_numpy()
    unit_present[i, ii] = True
    unit_strata[i] = strata_for(u)

# ---------------- modules ----------------
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
print("modules:", nM, flush=True)

# ---------------- core stats ----------------
def unit_z(xmat):
    """(n_units, nU) signed scores -> (n_units, nU) one-sided EXP z."""
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

def combine(zmat, drop_cluster=None):
    """zmat (n_units, nU) -> (fisher_EXP, per-pathogen p (n_pathogens, nU))."""
    zc = np.zeros((len(clusters), zmat.shape[1]))
    for ci, c in enumerate(clusters):
        if c == drop_cluster:
            continue
        iu = cl_units[c]
        zc[ci] = zmat[iu].sum(axis=0) / np.sqrt(len(iu))
    pp = np.ones((len(pathogens), zmat.shape[1]))
    for pi, p in enumerate(pathogens):
        cs = [clusters.index(c) for c in pa_clusters[p] if c != drop_cluster]
        if not cs:
            continue
        zp = zc[cs].sum(axis=0) / np.sqrt(len(cs))
        pp[pi] = 1.0 - norm.cdf(zp)
    pp = np.clip(pp, 1e-300, 1.0)
    fisher = -2.0 * np.log(pp).sum(axis=0)
    return fisher, pp

def module_fisher(zmat):
    """competitive diff-in-means per unit, then same combine. Returns (nM,) fisher."""
    mz = np.zeros((len(units), nM))
    for i in range(len(units)):
        z = zmat[i]
        pres = unit_present[i]
        zp = z[pres]
        tot, cnt = zp.sum(), pres.sum()
        for mi, idx in enumerate(mod_idx):
            m = idx[pres[idx]]
            if len(m) < 5:
                continue
            mz[i, mi] = (z[m].sum() - (tot - z[m].sum()) / max(cnt - len(m), 1) * len(m)) \
                        / np.sqrt(len(m))
    f, _ = combine(mz)
    return f

def bh_fast(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    ranked = p[o] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[o] = np.clip(ranked, 0, 1)
    return q

t0 = time.time()
z_obs = unit_z(unit_x)
obs_fisher, obs_pp = combine(z_obs)
obs_fisher_pro, obs_pp_pro = combine(-z_obs)
obs_mod = module_fisher(z_obs)
obs_mod_pro = module_fisher(-z_obs)
print(f"observed in {time.time()-t0:.1f}s", flush=True)

# ---------------- permutation ----------------
null_gene, null_mod = [], []
t0 = time.time()
done = 0
while done < B_TOTAL:
    b = min(BLOCK, B_TOTAL - done)
    fg = np.zeros((b, nU))
    fm = np.zeros((b, nM))
    for bi in range(b):
        xp = np.full_like(unit_x, np.nan)
        for i in range(len(units)):
            pres = unit_present[i]
            sp = unit_strata[i][pres]
            keys = rng.random(pres.sum()) + sp * 10.0
            xp[i][pres] = unit_x[i][pres][np.argsort(keys, kind="mergesort")]
        zp_ = unit_z(xp)
        fg[bi], _ = combine(zp_)
        fm[bi] = module_fisher(zp_)
    null_gene.append(fg)
    null_mod.append(fm)
    done += b
    el = time.time() - t0
    print(f"  perms {done}/{B_TOTAL} ({el:.0f}s eta {el/done*(B_TOTAL-done):.0f}s)", flush=True)

null_gene = np.concatenate(null_gene)
null_mod = np.concatenate(null_mod)
srt_g = np.sort(null_gene.ravel())[::-1]
srt_m = np.sort(null_mod.ravel())[::-1]

def emp_p(obs, srt):
    c = np.searchsorted(-srt, -obs)
    return (c + 1) / (srt.size + 1)

p_exp = emp_p(obs_fisher, srt_g)
p_pro = emp_p(obs_fisher_pro, srt_g)   # symmetric null
p_mod = emp_p(obs_mod, srt_m)
p_mod_pro = emp_p(obs_mod_pro, srt_m)

q_both = bh_fast(np.concatenate([p_exp, p_pro]))
q_exp, q_pro = q_both[:nU], q_both[nU:]
q_mod_both = bh_fast(np.concatenate([p_mod, p_mod_pro]))
q_mod, q_mod_pro = q_mod_both[:nM], q_mod_both[nM:]

# convergence gate: >=2 pathogens nominal p<0.01 AND >=2 clusters
n_cl_gene = {}
for g in U:
    pass  # filled below from v2 replication table
rep = pd.read_csv(ROOT / "06_within_pathogen/within_pathogen_replication_v2.tsv", sep="\t")
gene2ncl = rep.groupby("gene")["n_independence_clusters"].max().to_dict()

n_pat_exp = (obs_pp < 0.01).sum(axis=0)
n_pat_pro = (obs_pp_pro < 0.01).sum(axis=0)
ncl = np.array([gene2ncl.get(g, 1) for g in U])
conv_exp = (q_exp < 0.05) & (n_pat_exp >= 2) & (ncl >= 2)
conv_pro = (q_pro < 0.05) & (n_pat_pro >= 2) & (ncl >= 2)

# baseline essentiality annotation (therapy-window / confound flag)
rel_cells = ["A-549", "Calu-3", "HuH-7", "Huh-7.5", "Huh-7.5.1", "HEp-2", "HeLa"]
rel_cells = [c for c in rel_cells if c in bl.columns]
bl_sub = bl[rel_cells]
min_ch = np.array([bl_sub.loc[g].min() if g in bl_sub.index else np.nan for g in U])
frac_ess = np.array([(bl_sub.loc[g] < -1).mean() if g in bl_sub.index else np.nan for g in U])

gene_tab = pd.DataFrame({
    "gene": U,
    "fisher_EXP": obs_fisher, "p_EXP": p_exp, "q_EXP": q_exp,
    "n_pathogens_p01_EXP": n_pat_exp,
    "fisher_PRO": obs_fisher_pro, "p_PRO": p_pro, "q_PRO": q_pro,
    "n_pathogens_p01_PRO": n_pat_pro,
    "n_independence_clusters_max": ncl,
    "convergent_EXP": conv_exp, "convergent_PRO": conv_pro,
    "min_baseline_chronos": min_ch,
    "frac_backgrounds_essential": frac_ess,
    "p_SC2_EXP": obs_pp[pathogens.index("SARS-CoV-2")],
    "p_IAV_EXP": obs_pp[pathogens.index("Influenza")],
    "p_RSV_EXP": obs_pp[pathogens.index("RSV")],
    "p_SC2_PRO": obs_pp_pro[pathogens.index("SARS-CoV-2")],
    "p_IAV_PRO": obs_pp_pro[pathogens.index("Influenza")],
    "p_RSV_PRO": obs_pp_pro[pathogens.index("RSV")],
}).sort_values(["convergent_EXP", "q_EXP"], ascending=[False, True])
gene_tab.to_csv(ROOT / "07_cross_pathogen/convergence_meta_v2_gene.tsv", sep="\t", index=False)

mod_ess = [float(np.nanmean([(bl_sub.loc[g] < -1).mean() if g in bl_sub.index else np.nan
                             for g in modules[m]])) for m in mod_names]
mod_tab = pd.DataFrame({
    "module": mod_names, "n_genes": [len(m) for m in mod_idx],
    "fisher_EXP": obs_mod, "p_EXP": p_mod, "q_EXP": q_mod,
    "fisher_PRO": obs_mod_pro, "p_PRO": p_mod_pro, "q_PRO": q_mod_pro,
    "frac_members_essential": mod_ess,
}).sort_values("q_EXP")
mod_tab.to_csv(ROOT / "07_cross_pathogen/convergence_meta_v2_module.tsv", sep="\t", index=False)

# ---------------- LOSO ----------------
loso_rows = []
fdr_pass = np.where(conv_exp | conv_pro)[0]
for c in clusters:
    fE, _ = combine(z_obs, drop_cluster=c)
    fP, _ = combine(-z_obs, drop_cluster=c)
    pE = emp_p(fE, srt_g)
    pP = emp_p(fP, srt_g)
    qb = bh_fast(np.concatenate([pE, pP]))
    passE = qb[:nU] < 0.05
    passP = qb[nU:] < 0.05
    for gi in fdr_pass:
        gene = U[gi]
        was = "EXP" if conv_exp[gi] else "PRO"
        still = passE[gi] if was == "EXP" else passP[gi]
        loso_rows.append({"gene": gene, "direction": was, "dropped_cluster": c,
                          "still_fdr_pass": bool(still)})
loso = pd.DataFrame(loso_rows)
loso_sum = loso.groupby(["gene", "direction"])["still_fdr_pass"].agg(["sum", "count", "mean"]).reset_index()
loso_sum.columns = ["gene", "direction", "n_pass", "n_runs", "stability"]
loso_sum.to_csv(ROOT / "07_cross_pathogen/convergence_v2_loso.tsv", sep="\t", index=False)

np.savez_compressed(ROOT / "07_cross_pathogen/convergence_v2_null.npz",
                    gene=null_gene, module=null_mod)
json.dump({"B": B_TOTAL, "units": len(units), "clusters": len(clusters),
           "universe_protein_coding": nU, "n_modules": nM,
           "strata": "within-unit essentiality decile (baseline Chronos of mapped background); "
                     "expression stratum unavailable (DepMap expression not frozen) - AM2-13",
           "null": "pooled across genes/modules; PRO reuses EXP null by sign symmetry",
           "gate": "BH-FDR<0.05 (directions jointly corrected) AND >=2 pathogens nominal p<0.01 "
                   "AND >=2 independence clusters",
           "n_convergent_EXP": int(conv_exp.sum()), "n_convergent_PRO": int(conv_pro.sum())},
          open(ROOT / "07_cross_pathogen/convergence_v2_runinfo.json", "w"), indent=2)

print("\n=== RESULTS ===")
print("convergent EXP:", int(conv_exp.sum()), "| PRO:", int(conv_pro.sum()))
print(gene_tab[gene_tab["convergent_EXP"] | gene_tab["convergent_PRO"]].head(60).to_string(index=False))
print("\n=== modules FDR<0.05 ===")
print(mod_tab[(mod_tab["q_EXP"] < 0.05) | (mod_tab["q_PRO"] < 0.05)].to_string(index=False))
print("\n=== curated modules ===")
print(mod_tab[mod_tab["module"].str.startswith("CURATED")].to_string(index=False))
print("\nLOSO stability summary:")
print(loso_sum["stability"].describe().to_string())
