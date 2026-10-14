#!/usr/bin/env python3
"""Stage 9 (required-lightweight): candidate-gene expression in human lung
single-cell reference (CZ CELLxGENE Census, lung primary data incl. HLCA).

Two-pass memory-safe design (worker-0 has 16 GB):
  pass 1: obs only (soma_joinid, cell_type, dataset_id) for lung primary cells
  pass 2: stratified sample <= MAX_CELLS_PER_TYPE per cell type, then fetch X
          for sampled joinids x candidate genes only
Pre-registered rule: "lung_expressed" if >=1 cell type has pct>=10% AND mean>0.
Missing from reference = unknown (not negative). Upgrade-only.

Usage:
  python3 26_stage9_hlca_expression.py test
  python3 26_stage9_hlca_expression.py run
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace/HPICDA")
OUT = ROOT / "10_singlecell_context"
OUT.mkdir(exist_ok=True)

MAX_CELLS_PER_TYPE = 2000
PCT_THRESHOLD = 10.0
RNG_SEED = 20260831
CENSUS_VERSION = "2025-11-08"
OBS_FILTER = "tissue_general == 'lung' and is_primary_data == True"


def main(mode):
    import cellxgene_census

    if mode == "test":
        genes = ["ACE2", "TMPRSS2", "CTSL"]
    else:
        cand = pd.read_csv(ROOT / "02_metadata/preregistered_candidate_set_v1.tsv", sep="\t")
        genes = sorted(set(cand["gene"]))
    print(f"candidates: {len(genes)}", flush=True)

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        if mode == "test":
            ds = census["census_info"]["datasets"].read().concat().to_pandas()
            lung = ds[ds["dataset_title"].str.contains("lung|HLCA", case=False, na=False)]
            cols = [c for c in ["dataset_id", "dataset_title", "dataset_total_cell_count"] if c in ds.columns]
            print(lung[cols].head(10).to_string())
            return

        # pass 1: obs only
        print("pass 1: fetching obs...", flush=True)
        obs = cellxgene_census.get_obs(
            census, "Homo sapiens", value_filter=OBS_FILTER,
            column_names=["cell_type", "dataset_id"])
        print(f"  {len(obs)} lung primary cells, {obs['cell_type'].nunique()} cell types", flush=True)

        # stratified sample
        rng = np.random.default_rng(RNG_SEED)
        keep_idx = (obs.groupby("cell_type", group_keys=False)
                    .apply(lambda d: d.sample(n=min(len(d), MAX_CELLS_PER_TYPE),
                                              random_state=rng), include_groups=False)
                    .index)
        sampled = obs.loc[sorted(keep_idx)]
        joinids = (sampled["soma_joinid"] if "soma_joinid" in sampled.columns
                   else sampled.index).to_numpy()
        print(f"pass 2: fetching X for {len(joinids)} sampled cells x {len(genes)} genes...", flush=True)

        var_filter = "feature_name in [" + ",".join(f"'{g}'" for g in genes) + "]"
        adata = cellxgene_census.get_anndata(
            census, organism="Homo sapiens",
            obs_coords=joinids, var_value_filter=var_filter,
            column_names={"obs": ["cell_type"]})
    print(f"fetched {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    import scipy.sparse as sp
    X = sp.csr_matrix(adata.X)
    var_names = np.array(adata.var["feature_name"].astype(str))
    ct = adata.obs["cell_type"].to_numpy()
    rows = []
    for c in pd.unique(ct):
        idx = np.where(ct == c)[0]
        if len(idx) < 20:
            continue
        sub = X[idx]
        pct = np.asarray((sub > 0).sum(axis=0)).ravel() / sub.shape[0] * 100.0
        mean = np.asarray(sub.mean(axis=0)).ravel()
        for j, g in enumerate(var_names):
            rows.append({"cell_type": c, "gene": g, "n_cells": int(sub.shape[0]),
                         "pct_expressing": pct[j], "mean_expression": mean[j]})
    mat = pd.DataFrame(rows)
    mat.to_csv(OUT / "stage9_hlca_candidate_expression_by_celltype.tsv", sep="\t", index=False)

    summ = (mat.groupby("gene")
            .agg(max_pct=("pct_expressing", "max"),
                 max_mean=("mean_expression", "max"),
                 n_types_expressed=("pct_expressing", lambda s: int((s >= PCT_THRESHOLD).sum())))
            .reset_index())
    summ["lung_expressed"] = (summ["max_pct"] >= PCT_THRESHOLD) & (summ["max_mean"] > 0)
    missing = sorted(set(genes) - set(summ["gene"]))
    summ = pd.concat([summ, pd.DataFrame({"gene": missing, "lung_expressed": pd.NA})],
                     ignore_index=True)
    summ.to_csv(OUT / "stage9_candidate_lung_expression_summary.tsv", sep="\t", index=False)

    ds_counts = obs["dataset_id"].value_counts()
    meta = {"run_at_utc": datetime.now(timezone.utc).isoformat(),
            "census": CENSUS_VERSION, "max_cells_per_type": MAX_CELLS_PER_TYPE,
            "pct_threshold": PCT_THRESHOLD, "seed": RNG_SEED,
            "n_lung_primary_cells": int(len(obs)),
            "n_cell_types": int(obs["cell_type"].nunique()),
            "n_datasets": int(len(ds_counts)),
            "n_sampled_cells": int(len(joinids)),
            "n_candidates": len(genes), "n_missing_from_reference": len(missing),
            "n_lung_expressed": int(summ["lung_expressed"].eq(True).sum())}
    with open(OUT / "stage9_runinfo.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "test")
