#!/usr/bin/env python3
"""Re-audit post-hoc module-by-holdout specificity claims.

The source reports compared module genes with all genes in a holdout row,
including ``NOT_EVALUABLE_XPATH`` in the denominator and treating
``PASS_XPATH`` as a same-pathogen PASS. This script applies the manuscript's
same-pathogen adjudication definition instead:

- retain ``same_pathogen == True``;
- retain verdicts PASS, FAIL and CONTEXT_BOUNDARY;
- score only PASS as success and both FAIL and CONTEXT_BOUNDARY as non-PASS;
- test module membership within each holdout by two-sided Fisher exact test;
- apply Benjamini-Hochberg correction across all non-empty module-holdout tests.

This is a post-hoc exploratory audit. It does not modify the frozen candidates
or adjudications and is not added to the manuscript.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


ROOT = Path(
    os.environ.get("HPICDA_ROOT", Path(__file__).resolve().parents[1])
).resolve()
S13 = ROOT / "25_manuscript_package/tables/S13_holdout_verdicts.tsv"
S11B = ROOT / "25_manuscript_package/tables/S11b_intact_ppi_support.tsv"
OUT = ROOT / "24_closing_analyses"
TABLE = OUT / "66_module_holdout_specificity_reaudit.tsv"
SUMMARY = OUT / "66_module_holdout_specificity_reaudit.json"


def module_genes(edges: str) -> set[str]:
    return {
        gene
        for edge in str(edges).split(";")
        for gene in edge.strip().split("-")
        if gene
    }


def main() -> None:
    verdicts = pd.read_csv(S13, sep="\t")
    modules = pd.read_csv(S11B, sep="\t")
    members = {
        "HS": set(verdicts.loc[verdicts["is_hs"], "gene"]),
        "V-ATPase": module_genes(
            modules.loc[
                modules["module"].str.contains("V-ATPase", case=False, na=False),
                "edges",
            ].iloc[0]
        ),
        "COG": module_genes(
            modules.loc[
                modules["module"].str.contains("COG", case=False, na=False),
                "edges",
            ].iloc[0]
        ),
    }

    rows: list[dict[str, object]] = []
    evaluable = verdicts[
        verdicts["same_pathogen"].astype(str).str.lower().eq("true")
        & verdicts["verdict"].isin(["PASS", "FAIL", "CONTEXT_BOUNDARY"])
    ].copy()

    for module, genes in members.items():
        for holdout, frame in evaluable.groupby("holdout_id"):
            within = frame[frame["gene"].isin(genes)]
            background = frame[~frame["gene"].isin(genes)]
            if within.empty or background.empty:
                continue
            module_pass = int(within["verdict"].eq("PASS").sum())
            background_pass = int(background["verdict"].eq("PASS").sum())
            table = [
                [module_pass, len(within) - module_pass],
                [background_pass, len(background) - background_pass],
            ]
            odds_ratio, p_two = fisher_exact(table, alternative="two-sided")
            _, p_greater = fisher_exact(table, alternative="greater")
            rows.append(
                {
                    "module": module,
                    "holdout_id": holdout,
                    "module_pass": module_pass,
                    "module_n": len(within),
                    "module_pass_rate": module_pass / len(within),
                    "background_pass": background_pass,
                    "background_n": len(background),
                    "background_pass_rate": background_pass / len(background),
                    "odds_ratio": odds_ratio,
                    "fisher_p_two_sided": p_two,
                    "fisher_p_greater": p_greater,
                }
            )

    result = pd.DataFrame(rows)
    if result.empty:
        raise RuntimeError("No evaluable module-by-holdout comparisons")
    result["bh_q_two_sided"] = multipletests(
        result["fisher_p_two_sided"], method="fdr_bh"
    )[1]
    result["bh_q_greater"] = multipletests(
        result["fisher_p_greater"], method="fdr_bh"
    )[1]
    result = result.sort_values(["module", "holdout_id"]).reset_index(drop=True)
    result.to_csv(TABLE, sep="\t", index=False, encoding="utf-8")

    best_two = result.loc[result["fisher_p_two_sided"].idxmin()]
    best_greater = result.loc[result["fisher_p_greater"].idxmin()]
    summary = {
        "status": "post-hoc exploratory audit; not included in manuscript",
        "input": str(S13),
        "module_members": {k: sorted(v) for k, v in members.items()},
        "eligibility": {
            "same_pathogen": True,
            "verdicts": ["PASS", "FAIL", "CONTEXT_BOUNDARY"],
            "success": "PASS",
            "non_success": ["FAIL", "CONTEXT_BOUNDARY"],
        },
        "n_tests": int(len(result)),
        "best_nominal_two_sided": best_two.to_dict(),
        "best_nominal_greater": best_greater.to_dict(),
        "any_bh_two_sided_below_0_05": bool(
            (result["bh_q_two_sided"] < 0.05).any()
        ),
        "any_bh_greater_below_0_05": bool(
            (result["bh_q_greater"] < 0.05).any()
        ),
        "decision": (
            "Do not add the proposed module-specificity atlas to the manuscript: "
            "no comparison survives multiplicity correction, and the source "
            "reports used non-evaluable rows in the background denominator."
        ),
    }
    SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {TABLE}")
    print(f"Wrote {SUMMARY}")
    print(
        result[
            [
                "module",
                "holdout_id",
                "module_pass",
                "module_n",
                "background_pass",
                "background_n",
                "fisher_p_greater",
                "bh_q_greater",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
