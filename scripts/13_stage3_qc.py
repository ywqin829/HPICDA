#!/usr/bin/env python3
"""HPICDA Stage 3: per-screen QC (blueprint §11, 12 checks).

Reads each Track A study's MAGeCK count matrix + gene summaries, computes QC
metrics, assigns per-screen status, writes 03_functional_screens/<study>/QC_REPORT.md
and a machine-readable stage3_qc_summary.tsv.
"""
import csv, json, math
from pathlib import Path
from statistics import median

ROOT = Path("/workspace/HPICDA/03_functional_screens")
OUT_SUMMARY = ROOT / "stage3_qc_summary.tsv"

def read_counts(path):
    """Return (guides, genes, samples, matrix as dict guide->[counts])."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[2:]
        guides, genes, rows = [], [], {}
        for line in f:
            p = line.rstrip("\n").split("\t")
            guides.append(p[0]); genes.append(p[1])
            rows[p[0]] = [float(x) for x in p[2:]]
    return guides, genes, samples, rows

def spearman(x, y):
    n = len(x)
    rx = rankdata(x); ry = rankdata(y)
    mx = sum(rx)/n; my = sum(ry)/n
    cov = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
    vx = sum((a-mx)**2 for a in rx); vy = sum((b-my)**2 for b in ry)
    return cov/math.sqrt(vx*vy) if vx > 0 and vy > 0 else float("nan")

def rankdata(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0]*len(v)
    i = 0
    while i < len(v):
        j = i
        while j+1 < len(v) and v[order[j+1]] == v[order[i]]:
            j += 1
        avg = (i+j)/2.0 + 1
        for k in range(i, j+1):
            r[order[k]] = avg
        i = j+1
    return r

def gini(vals):
    v = sorted(vals)
    n = len(v); s = sum(v)
    if s == 0: return 0.0
    cum = 0.0
    for i, x in enumerate(v, 1):
        cum += i*x
    return (2*cum)/(n*s) - (n+1)/n

def qc_study(study_id, count_path, gene_summaries, known_pos, known_note,
             control_cols, treat_cols_map, lib_size, author_note=""):
    guides, genes, samples, rows = read_counts(count_path)
    n_guides = len(guides)
    n_genes = len(set(genes))
    # per-sample stats
    lib_sizes = {s: 0 for s in samples}
    zero_frac = {}
    for j, s in enumerate(samples):
        col = [rows[g][j] for g in guides]
        lib_sizes[s] = int(sum(col))
        zero_frac[s] = sum(1 for c in col if c == 0)/n_guides
    # coverage = guides with >=1 read in >=1 sample
    detected = sum(1 for g in guides if any(rows[g][j] > 0 for j in range(len(samples))))
    coverage = detected/n_guides
    # replicate correlation on control columns (log10(x+1))
    rep_corr = None
    if len(control_cols) >= 2:
        idx = [samples.index(c) for c in control_cols[:2]]
        x = [math.log10(rows[g][idx[0]]+1) for g in guides]
        y = [math.log10(rows[g][idx[1]]+1) for g in guides]
        rep_corr = round(spearman(x, y), 4)
    # bottleneck: gini on treated
    gini_t = {}
    for tname, tcols in treat_cols_map.items():
        j = samples.index(tcols[0])
        gini_t[tname] = round(gini([rows[g][j] for g in guides]), 4)
    # known positive controls from first gene summary
    known_rec = {}
    if gene_summaries:
        gs = gene_summaries[0]
        with open(gs) as f:
            gsum = {r["id"]: r for r in csv.DictReader(f, delimiter="\t")}
        for kg in known_pos:
            if kg in gsum:
                r = gsum[kg]
                known_rec[kg] = {"pos_rank": r["pos|rank"], "pos_fdr": r["pos|fdr"]}
    # guides/gene
    from collections import Counter
    gc = Counter(genes)
    med_gpg = median(gc.values())
    return {
        "study_id": study_id, "n_guides": n_guides, "n_genes": n_genes,
        "median_guides_per_gene": med_gpg, "library_size_expected": lib_size,
        "coverage": round(coverage, 4), "lib_sizes": lib_sizes,
        "zero_frac": {k: round(v, 4) for k, v in zero_frac.items()},
        "control_rep_spearman": rep_corr, "gini_treated": gini_t,
        "known_controls": known_rec, "known_note": known_note,
        "author_note": author_note,
    }

# ---- study configurations ----
studies = []

studies.append(dict(
    study_id="GSE210689", count_path=ROOT/"raw/GSE210689/counts_mageck.txt",
    gene_summaries=[ROOT/"raw/GSE210689/mtb_pro_vs_ui.gene_summary.txt"],
    known_pos=["GPR156", "MMGT1"], known_note="MMGT1-GPR156 lipid-droplet axis (PMID 37269834)",
    control_cols=["UI1", "UI2", "UI3"],
    treat_cols_map={"Pro": ["Pro1"], "Lat": ["Lat1"]},
    lib_size=72618,
    author_note="LOW DEPTH ~1.2 reads/guide; >50% zero-count in UI/Lat; no FDR<0.25 hit; GPR156 directional (pos p=0.0038).",
))

studies.append(dict(
    study_id="GSE178699", count_path=ROOT/"raw/GSE178699/ORCS-2066_counts.txt",
    gene_summaries=[ROOT/"raw/GSE178699/ORCS-2066.gene_summary.txt"],
    known_pos=["ACE2", "TMPRSS2", "KMT2C", "KDM6A"], known_note="ACE2/TMPRSS2 entry factors",
    control_cols=["CALU3.T25A", "CALU3.T25B", "CALU3.T25C"], treat_cols_map={"SARS": ["CALU3.T43A.lowMOI.SARS"]},
    lib_size=None, author_note="qGI-based author method; MAGeCK reanalysis. ACE2 #1 all 3 lines.",
))

studies.append(dict(
    study_id="GSE111166", count_path=ROOT/"raw/GSE111166/counts_mageck.txt",
    gene_summaries=[ROOT/"raw/GSE111166/ORCS-1148.gene_summary.txt"],
    known_pos=["SLC35A1", "CIC", "PIGN", "TRIM23"], known_note="SLC35A1 sialic-acid transporter",
    control_cols=["PREa"], treat_cols_map={"Rd2": ["Rd2_rep1"], "Rd5": ["Rd5_rep1"]},
    lib_size=63950, author_note="Sequential screen; ORCS-1150 bottlenecked (RRA over LFC rule).",
))

studies.append(dict(
    study_id="Li2020", count_path=ROOT/"raw/Li2020_Edinburgh/ORCS-1207-S2_counts.txt",
    gene_summaries=[ROOT/"raw/Li2020_Edinburgh/ORCS-1207-S2.gene_summary.txt"],
    known_pos=["WDR7", "SLC35A1", "ATP6V1A", "ATP6AP1", "CMTR1", "CCDC115", "TMEM199"],
    known_note="V-ATPase + SLC35A1 (FACS-HA screen)",
    control_cols=["P4"], treat_cols_map={"P6": ["P6"]},
    lib_size=None, author_note="log2RPM reverse-transformed to pseudo-counts (documented approximation).",
))

studies.append(dict(
    study_id="GSE175666", count_path=ROOT/"raw/GSE175666/ORCS-1862_counts.txt",
    gene_summaries=[ROOT/"raw/GSE175666/ORCS-1862.gene_summary.txt"],
    known_pos=["ACE2", "AP1G1", "TMPRSS2"], known_note="ACE2/AP1G1 bidirectional screen",
    control_cols=["B-SS_D0_r1_run1", "D-SS_D0_R2_run1"], treat_cols_map={"D4": ["E-SS_D4_MOI01_r1_run1"]},
    lib_size=None, author_note="Author zscore sheet matches reanalysis (AP1G1/ACE2 #1/#2).",
))

studies.append(dict(
    study_id="PRJNA1281304", count_path=ROOT/"raw/PRJNA1281304/counts_mageck.txt",
    gene_summaries=[ROOT/"raw/PRJNA1281304/RSV_viability.gene_summary.txt"],
    known_pos=["SLC35B2", "EXT1", "EXT2", "B3GAT3", "COG4"], known_note="heparan-sulfate + COG (RSV entry)",
    control_cols=["CTRL_A", "CTRL_B"], treat_cols_map={"viability": ["viability_A"], "GFP_NEG": ["A_GFP_NEG"]},
    lib_size=123411, author_note="GeCKOv2 A+B combined; ~31% empty-vector reads excluded.",
))

studies.append(dict(
    study_id="GSE154751", count_path=ROOT/"raw/GSE154751/counts_mageck.txt",
    gene_summaries=[ROOT/"raw/GSE154751/Pa_ExoU.gene_summary.txt"],
    known_pos=["DNAJC5"], known_note="ExoU toxin sensitivity (context-only)",
    control_cols=["Banque"], treat_cols_map={"ExoU": ["ExoU_1"]},
    lib_size=71090, author_note="Single control replicate (library bank); context-only screen.",
))

results = []
for st in studies:
    try:
        r = qc_study(**st)
        results.append(r)
        print(f"OK {r['study_id']}: cov={r['coverage']}, rep_rho={r['control_rep_spearman']}, "
              f"guides={r['n_guides']}, genes={r['n_genes']}")
    except Exception as e:
        print(f"FAIL {st['study_id']}: {e}")
        results.append({"study_id": st["study_id"], "error": str(e)})

with open(ROOT / "stage3_qc_raw.json", "w") as f:
    json.dump(results, f, indent=2)
print("saved stage3_qc_raw.json")
