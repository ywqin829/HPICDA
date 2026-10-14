#!/usr/bin/env python3
"""Script 24: freeze pre-registered candidate set v1 (plan v2, blindness rule).

MUST run AFTER R4 full run (convergence_meta_v2_gene.tsv, B=10000) and
BEFORE Stage 8 run-mode / Stage 9 / Stage 10 (human-layer analyses).

Candidate set = union of:
  (a) genes with wp_status == WP_REPLICATED in within_pathogen_replication_v2.tsv
  (b) genes in cross_pathogen_convergence.tsv (threshold-intersection, v1)
  (c) genes with convergent_EXP or convergent_PRO in convergence_meta_v2_gene.tsv

Output: 02_metadata/preregistered_candidate_set_v1.tsv
"""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "02_metadata/preregistered_candidate_set_v1.tsv"

def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    rep = pd.read_csv(ROOT / "06_within_pathogen/within_pathogen_replication_v2.tsv", sep="\t")
    thr = pd.read_csv(ROOT / "07_cross_pathogen/cross_pathogen_convergence.tsv", sep="\t")
    meta = pd.read_csv(ROOT / "07_cross_pathogen/convergence_meta_v2_gene.tsv", sep="\t")

    # (a) WP_REPLICATED genes (T3/T4 only per tier definitions)
    a = rep[rep["wp_status"] == "WP_REPLICATED"]
    genes_a = set(a["gene"])
    dir_a = a.groupby("gene")["consensus_direction"].agg(lambda s: ";".join(sorted(set(s))))
    pat_a = a.groupby("gene")["pathogen"].agg(lambda s: ";".join(sorted(set(s))))

    # (b) threshold-intersection convergent genes (v1)
    genes_b = set(thr["gene"])
    dir_b = thr.groupby("gene")["direction"].agg(lambda s: ";".join(sorted(set(s))))
    pat_b = thr.groupby("gene")["pathogens"].agg(lambda s: ";".join(sorted(set(s))))

    # (c) R4 meta-analysis convergent genes (FDR gate already applied by script 22)
    c = meta[meta["convergent_EXP"] | meta["convergent_PRO"]]
    genes_c = set(c["gene"])
    dir_c = c.groupby("gene").apply(
        lambda d: ";".join(sorted(set(
            (["EXP"] if d["convergent_EXP"].any() else []) +
            (["PRO"] if d["convergent_PRO"].any() else [])))),
        include_groups=False)

    all_genes = sorted(genes_a | genes_b | genes_c)
    rows = []
    for g in all_genes:
        src = []
        if g in genes_a: src.append("WP_REPLICATED")
        if g in genes_b: src.append("THRESHOLD_INTERSECTION")
        if g in genes_c: src.append("META_CONVERGENT")
        dirs, pats = set(), set()
        for d in (dir_a.get(g), dir_b.get(g), dir_c.get(g)):
            if isinstance(d, str):
                dirs.update(d.split(";"))
        for p in (pat_a.get(g), pat_b.get(g)):
            if isinstance(p, str):
                pats.update(p.split(";"))
        rows.append({
            "gene": g,
            "sources": ";".join(src),
            "n_sources": len(src),
            "directions": ";".join(sorted(dirs)),
            "pathogens": ";".join(sorted(pats)),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)

    manifest = {
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "rule": "union of WP_REPLICATED + THRESHOLD_INTERSECTION + META_CONVERGENT genes",
        "inputs": {
            "within_pathogen_replication_v2.tsv": sha256(ROOT / "06_within_pathogen/within_pathogen_replication_v2.tsv"),
            "cross_pathogen_convergence.tsv": sha256(ROOT / "07_cross_pathogen/cross_pathogen_convergence.tsv"),
            "convergence_meta_v2_gene.tsv": sha256(ROOT / "07_cross_pathogen/convergence_meta_v2_gene.tsv"),
        },
        "n_candidates": len(out),
        "n_by_source": {
            "WP_REPLICATED": len(genes_a),
            "THRESHOLD_INTERSECTION": len(genes_b),
            "META_CONVERGENT": len(genes_c),
        },
        "blindness": "frozen before any Stage 8/9/10 human-layer result was computed",
    }
    with open(ROOT / "02_metadata/preregistered_candidate_set_v1.manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps(manifest, indent=2))
    print("\nby n_sources:")
    print(out["n_sources"].value_counts().sort_index().to_string())

if __name__ == "__main__":
    main()
