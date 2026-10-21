#!/usr/bin/env python3
"""Opportunity-balanced audit of the HPICDA validation benchmark.

The original gene-level endpoint is "at least one same-pathogen PASS". Genes
with more evaluable holdout pairs have more chances to pass. This script:

1. quantifies that opportunity bias;
2. repeats the gene-level benchmark after randomly selecting exactly one
   same-pathogen holdout pair per gene;
3. reports whether the simple consensus features retain their apparent
   advantage over the frozen HG3 score.

The script is sensitivity analysis only. It does not change the frozen
candidate set, holdout verdicts, or HG3 scores.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from sklearn.metrics import roc_auc_score
from statsmodels.formula.api import logit


ROOT = Path(
    os.environ.get("HPICDA_ROOT", Path(__file__).resolve().parents[1])
).resolve()
OUT = ROOT / "24_closing_analyses"
TABLE_DIR = ROOT / "25_manuscript_package/tables"
RNG = np.random.default_rng(20261021)
N_RESAMPLES = 5000


def fit_logit_table(frame: pd.DataFrame, formulas: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for formula in formulas:
        model = logit(formula, data=frame).fit(disp=0)
        for term in model.params.index:
            if term == "Intercept":
                continue
            ci = np.exp(model.conf_int().loc[term])
            rows.append(
                {
                    "model": formula,
                    "term": term,
                    "n": int(model.nobs),
                    "aic": float(model.aic),
                    "odds_ratio": float(np.exp(model.params[term])),
                    "ci_lo": float(ci.iloc[0]),
                    "ci_hi": float(ci.iloc[1]),
                    "p": float(model.pvalues[term]),
                }
            )
    return pd.DataFrame(rows)


def percentile_summary(values: np.ndarray) -> dict[str, float]:
    lo, median, hi = np.percentile(values, [2.5, 50, 97.5])
    return {
        "median": float(median),
        "ci_lo": float(lo),
        "ci_hi": float(hi),
        "n_values": int(len(values)),
    }


def main() -> None:
    verdicts = pd.read_csv(
        ROOT / "25_manuscript_package/tables/S13_holdout_verdicts.tsv",
        sep="\t",
    )
    benchmark = pd.read_csv(
        OUT / "benchmark_gene_level_data.tsv",
        sep="\t",
    )

    evaluable_verdicts = {"PASS", "FAIL", "CONTEXT_BOUNDARY"}
    pairs = verdicts[
        (verdicts["same_pathogen"].astype(str).str.lower() == "true")
        & verdicts["verdict"].isin(evaluable_verdicts)
    ].copy()
    pairs["pass"] = (pairs["verdict"] == "PASS").astype(int)
    pairs = pairs.sort_values(["gene", "holdout_id"]).reset_index(drop=True)

    features = benchmark.set_index("gene")
    genes = sorted(pairs["gene"].unique())
    if len(genes) != 333:
        raise RuntimeError(f"Expected 333 benchmark-evaluable genes, found {len(genes)}")

    gene_rows: list[dict[str, object]] = []
    grouped_pairs: list[np.ndarray] = []
    for gene in genes:
        block = pairs[pairs["gene"] == gene]
        grouped_pairs.append(block["pass"].to_numpy(dtype=np.int8))
        row = features.loc[gene].to_dict()
        row.update(
            {
                "gene": gene,
                "any_pass": int(block["pass"].max()),
                "n_same_pathogen_pairs": int(len(block)),
                "n_same_pathogen_pass": int(block["pass"].sum()),
            }
        )
        gene_rows.append(row)

    gene_table = pd.DataFrame(gene_rows)
    gene_table["log1p_pairs"] = np.log1p(gene_table["n_same_pathogen_pairs"])

    rho_pathogen, p_pathogen = spearmanr(
        gene_table["n_pathogens_hit"], gene_table["n_same_pathogen_pairs"]
    )
    rho_units, p_units = spearmanr(
        gene_table["n_units_hit"], gene_table["n_same_pathogen_pairs"]
    )

    model_table = fit_logit_table(
        gene_table,
        [
            "any_pass ~ n_pathogens_hit",
            "any_pass ~ n_pathogens_hit + log1p_pairs",
            "any_pass ~ n_units_hit",
            "any_pass ~ n_units_hit + log1p_pairs",
            "any_pass ~ hit_frequency",
            "any_pass ~ hit_frequency + log1p_pairs",
            "any_pass ~ n_screens_hit + log1p_pairs",
            "any_pass ~ n_pathogens_hit + n_same_pathogen_pairs",
            "any_pass ~ n_units_hit + n_same_pathogen_pairs",
        ],
    )

    score_columns = {
        "HG3_total_score": gene_table["total_score"].to_numpy(dtype=float),
        "n_pathogens_hit": gene_table["n_pathogens_hit"].to_numpy(dtype=float),
        "n_units_hit": gene_table["n_units_hit"].to_numpy(dtype=float),
        "n_screens_hit": gene_table["n_screens_hit"].to_numpy(dtype=float),
        "hit_frequency": gene_table["hit_frequency"].to_numpy(dtype=float),
        "best_rank_pct": -gene_table["best_rank_pct"].to_numpy(dtype=float),
    }

    replicate_aucs: dict[str, list[float]] = {name: [] for name in score_columns}
    paired_differences: dict[str, list[float]] = {
        "n_pathogens_hit_minus_HG3": [],
        "n_units_hit_minus_HG3": [],
        "hit_frequency_minus_HG3": [],
    }
    top10_rows: list[dict[str, float]] = []

    for _ in range(N_RESAMPLES):
        sampled_outcome = np.array(
            [RNG.choice(group, size=1)[0] for group in grouped_pairs],
            dtype=np.int8,
        )
        if len(np.unique(sampled_outcome)) < 2:
            continue
        aucs = {
            name: float(roc_auc_score(sampled_outcome, score))
            for name, score in score_columns.items()
        }
        for name, value in aucs.items():
            replicate_aucs[name].append(value)
        for name in paired_differences:
            predictor = name.removesuffix("_minus_HG3")
            paired_differences[name].append(aucs[predictor] - aucs["HG3_total_score"])

        top = np.argsort(-score_columns["HG3_total_score"], kind="mergesort")[:10]
        top_pass = int(sampled_outcome[top].sum())
        top_rate = top_pass / 10
        rest_pass = int(sampled_outcome.sum()) - top_pass
        odds_ratio, p_value = fisher_exact(
            [
                [top_pass, 10 - top_pass],
                [rest_pass, len(sampled_outcome) - 10 - rest_pass],
            ]
        )
        top10_rows.append(
            {
                "top10_pass_rate": top_rate,
                "odds_ratio": float(odds_ratio),
                "p": float(p_value),
            }
        )

    resample_summary = {
        "n_resamples_requested": N_RESAMPLES,
        "n_resamples_used": len(replicate_aucs["HG3_total_score"]),
        "n_genes": len(genes),
        "n_evaluable_pairs": int(len(pairs)),
        "n_pass_pairs": int(pairs["pass"].sum()),
        "opportunity_association": {
            "n_pathogens_hit_vs_n_pairs_spearman": float(rho_pathogen),
            "n_pathogens_hit_vs_n_pairs_p": float(p_pathogen),
            "n_units_hit_vs_n_pairs_spearman": float(rho_units),
            "n_units_hit_vs_n_pairs_p": float(p_units),
        },
        "equal_opportunity_auroc": {
            name: percentile_summary(np.array(values))
            for name, values in replicate_aucs.items()
        },
        "equal_opportunity_paired_differences": {
            name: {
                **percentile_summary(np.array(values)),
                "fraction_at_or_below_zero": float(np.mean(np.array(values) <= 0)),
            }
            for name, values in paired_differences.items()
        },
        "equal_opportunity_top10_HG3": {
            name: percentile_summary(np.array([row[name] for row in top10_rows]))
            for name in ("top10_pass_rate", "odds_ratio", "p")
        },
        "interpretation_rule": (
            "The original any-PASS endpoint gives genes with more same-pathogen pairs "
            "more opportunities to pass. Equal-opportunity resampling selects one pair "
            "per gene per replicate. A predictor is regarded as retaining evidence only "
            "when the equal-opportunity paired-difference interval excludes zero."
        ),
    }

    gene_table.to_csv(
        OUT / "61_opportunity_balanced_gene_table.tsv", sep="\t", index=False
    )
    model_table.to_csv(
        OUT / "61_opportunity_adjusted_models.tsv", sep="\t", index=False
    )
    pd.DataFrame(top10_rows).to_csv(
        OUT / "61_equal_opportunity_top10_replicates.tsv", sep="\t", index=False
    )
    with (OUT / "61_opportunity_balanced_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(resample_summary, handle, indent=2, ensure_ascii=False)

    paired_table = pd.DataFrame(
        [
            {
                "comparison": name,
                **percentile_summary(np.array(values)),
                "fraction_at_or_below_zero": float(np.mean(np.array(values) <= 0)),
            }
            for name, values in paired_differences.items()
        ]
    )
    paired_table.to_csv(
        OUT / "61_equal_opportunity_paired_differences.tsv",
        sep="\t",
        index=False,
    )
    model_table.to_csv(
        TABLE_DIR / "S21a_opportunity_adjusted_models.tsv",
        sep="\t",
        index=False,
    )
    summary_rows: list[dict[str, object]] = []
    for name, values in replicate_aucs.items():
        stats = percentile_summary(np.array(values))
        summary_rows.append(
            {
                "section": "equal_opportunity_auroc",
                "metric": name,
                "median": stats["median"],
                "ci_lo": stats["ci_lo"],
                "ci_hi": stats["ci_hi"],
                "n": stats["n_values"],
                "fraction_at_or_below_zero": "",
            }
        )
    for name, values in paired_differences.items():
        values_array = np.array(values)
        stats = percentile_summary(values_array)
        summary_rows.append(
            {
                "section": "paired_difference_vs_HG3",
                "metric": name,
                "median": stats["median"],
                "ci_lo": stats["ci_lo"],
                "ci_hi": stats["ci_hi"],
                "n": stats["n_values"],
                "fraction_at_or_below_zero": float(np.mean(values_array <= 0)),
            }
        )
    for name in ("top10_pass_rate", "odds_ratio", "p"):
        values_array = np.array([row[name] for row in top10_rows])
        stats = percentile_summary(values_array)
        summary_rows.append(
            {
                "section": "equal_opportunity_top10_HG3",
                "metric": name,
                "median": stats["median"],
                "ci_lo": stats["ci_lo"],
                "ci_hi": stats["ci_hi"],
                "n": stats["n_values"],
                "fraction_at_or_below_zero": "",
            }
        )
    pd.DataFrame(summary_rows).to_csv(
        TABLE_DIR / "S21b_equal_opportunity_summary.tsv",
        sep="\t",
        index=False,
    )

    print(json.dumps(resample_summary, indent=2, ensure_ascii=False))
    print("\nAdjusted models:")
    print(model_table.to_string(index=False))


if __name__ == "__main__":
    main()
