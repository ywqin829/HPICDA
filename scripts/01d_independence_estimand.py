#!/usr/bin/env python3
"""HPICDA Stage 0 — publication_independence.tsv + screen_estimand_precheck.tsv.

Derived from the Stage 0 manifest. Screen-level metadata only.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

MANIFEST = "01_stage0_audit/data_availability_manifest.tsv"
INDEP_OUT = "01_stage0_audit/publication_independence.tsv"
ESTIMAND_OUT = "01_stage0_audit/screen_estimand_precheck.tsv"


def main() -> None:
    rows = list(csv.DictReader(open(MANIFEST), delimiter="\t"))

    # ---- publication independence ----
    clusters = defaultdict(list)
    for r in rows:
        clusters[r["independence_cluster"]].append(r)

    ind_fields = [
        "independence_cluster", "pmid", "n_screens", "pathogens", "cell_lines",
        "host_species", "shares_control_within_cluster", "independent_of_other_clusters",
        "nonhuman_excluded", "notes",
    ]
    with open(INDEP_OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ind_fields, delimiter="\t")
        w.writeheader()
        for cid, rs in sorted(clusters.items()):
            human = [r for r in rs if r["host_species"] == "Homo sapiens"]
            shared = any(r["shared_control_group"].startswith("yes") for r in rs)
            nonhuman = any(r["host_species"] != "Homo sapiens" for r in rs)
            w.writerow({
                "independence_cluster": cid,
                "pmid": ";".join(sorted({r["pmid"] for r in rs if r["pmid"]})),
                "n_screens": len(rs),
                "pathogens": ";".join(sorted({r["pathogen"] for r in rs})),
                "cell_lines": ";".join(sorted({r["cell_line"] for r in rs})),
                "host_species": ";".join(sorted({r["host_species"] for r in rs})),
                "shares_control_within_cluster": "yes" if shared else "no/unknown",
                # True by construction here: clusters are keyed by distinct
                # publication/accession. Cross-lab sample sharing cannot be
                # excluded from metadata alone -> verify at Stage 2 (HG2).
                "independent_of_other_clusters": "yes_metadata_level_verify_HG2",
                "nonhuman_excluded": "yes" if nonhuman and not human else "n/a",
                "notes": f"{len(human)} human screens",
            })

    # ---- screen estimand precheck ----
    est_fields = [
        "study_id", "pathogen", "cell_line", "comparison_arm", "reference_arm",
        "effect_metric", "effect_sign_definition", "SE_provenance",
        "meta_eligible", "estimand_risk_flags",
    ]
    with open(ESTIMAND_OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=est_fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            if r["host_species"] != "Homo sapiens":
                continue
            whole = r["whole_infection_yes_no"] == "yes"
            toxin = r["toxin_only_yes_no"] == "yes"
            strat = r["selection_strategy"]
            if "FACS" in strat:
                comp = "sorted infection/reporter-defined population"
                ref = "unsorted or antigen/reporter-negative population"
            elif "survival" in strat:
                comp = "survivors after selection"
                ref = "initial pool / uninfected control"
            else:
                comp = "selected population (TBD at Stage 2)"
                ref = "reference population (TBD at Stage 2)"
            metric = "MAGeCK/RRA/L2FC (typical; confirm per-screen at Stage 2)"
            sign = ("positive enrichment in comparison arm = host factor whose loss "
                    "is protective/dependency-reducing; sign convention to be frozen "
                    "per-screen in screen_estimand_registry at Stage 2")
            risks = []
            if toxin:
                risks.append("TOXIN_ONLY_not_whole_infection")
            if not whole:
                risks.append("endpoint_not_whole_infection")
            if r["genome_wide_yes_no"] != "yes":
                risks.append("not_genome_wide")
            if r["raw_counts_available"].startswith("verify"):
                risks.append("raw_counts_unverified")
            if r["uninfected_baseline_available"].startswith("verify"):
                risks.append("baseline_unverified")
            meta = ("yes" if (whole and r["genome_wide_yes_no"] == "yes"
                              and not toxin and not risks) else
                    "conditional" if (whole and r["genome_wide_yes_no"] == "yes") else "no")
            w.writerow({
                "study_id": r["study_id"], "pathogen": r["pathogen"],
                "cell_line": r["cell_line"], "comparison_arm": comp,
                "reference_arm": ref, "effect_metric": metric,
                "effect_sign_definition": sign,
                "SE_provenance": ("recomputable_from_raw_counts"
                                  if r["raw_counts_available"].startswith("yes")
                                  else "publication_stats_only_or_unverified"),
                "meta_eligible": meta,
                "estimand_risk_flags": ";".join(risks),
            })

    print(f"Wrote {INDEP_OUT} ({len(clusters)} clusters)")
    print(f"Wrote {ESTIMAND_OUT}")


if __name__ == "__main__":
    main()
