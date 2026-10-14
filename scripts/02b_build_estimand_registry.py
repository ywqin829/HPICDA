#!/usr/bin/env python3
"""HPICDA Stage 2a: build screen_estimand_registry.tsv (blueprint v1.1 §10, 26 fields).

One row per interpretable contrast/estimand. Holdout screens are registered at
metadata level only (holdout_role=holdout_locked) — no gene-level data ingestion
until HG4. meta_eligible is machine-derived from track + SE provenance rules.
"""
import csv, os, hashlib, datetime

ROOT = "/workspace/HPICDA"
OUT = f"{ROOT}/03_functional_screens"
os.makedirs(OUT, exist_ok=True)

man = list(csv.DictReader(open(f"{ROOT}/01_stage0_audit/data_availability_manifest.tsv"), delimiter="\t"))
adj = {r["study_id"]: r for r in csv.DictReader(open(f"{ROOT}/02_ontology/screen_endpoint_direction_adjudication.tsv"), delimiter="\t")}
pre = {r["study_id"]: r for r in csv.DictReader(open(f"{ROOT}/01_stage0_audit/screen_estimand_precheck.tsv"), delimiter="\t")}

HOLDOUT = {"ORCS-2360", "ORCS-2361", "ORCS-2068", "ORCS-2069", "ORCS-2076", "ORCS-2077",
           "ORCS-2428", "ORCS-2430", "ORCS-1524", "ORCS-2431", "ORCS-1233", "ORCS-1234"}
CONTEXT = {"LIT-Pa-GSE154751", "LIT-Pa-SRP514827"}
AUXILIARY = {"ORCS-2020"}  # Influenza B near-neighbor
EXCLUDED = {"ORCS-2433"}   # HPIV-3

# Track A raw-count routes verified against GEO on 2026-08-30 (PMID links from gds esummary):
#   GSE175666 = PMID 35879413 (Rebendenne Nat Genet, bidirectional Gattinara+Calabrese)
#   GSE177545 = PMID 35086559 (Biering, HEK293T/ACE2 MOI series + OC43)
#   GSE178699 = PMID 36597481 (Prelli Bozzo multi-model TKOv3; GEO lacks PMID link,
#               design match: Calu-3/HuH-7/HEK293+ACE2+TMPRSS2 + Vero, 35 samples)
# Stage 0 viral audit had mis-attributed GSE177545 to 36597481 -> erratum #4.
TRACK_A = {  # by study_id (literature screens without PMID-keyed route)
    "LIT-Mtb-GSE210689": "GSE210689 counts matrix (GEO)",
    "LIT-RSV-PRJNA1281304-viability": "SRA PRJNA1281304 (~8 GB, 8 runs)",
    "LIT-RSV-PRJNA1281304-early": "SRA PRJNA1281304 (~8 GB, 8 runs)",
}
TRACK_A_BY_PMID = {
    "37269834": "GSE210689 counts (GEO)",
    "29642015": "GSE111166 counts (GEO)",
    "31919360": "Edinburgh doi:10.7488/ds/7818 counts",
    "35879413": "GSE175666 counts (GEO)",
    "35086559": "GSE177545 counts (GEO)",
    "36597481": "GSE178699 counts (GEO)",
}

def role(sid):
    if sid in EXCLUDED: return "excluded"
    if sid in HOLDOUT: return "holdout_locked"
    if sid in CONTEXT: return "context"
    if sid in AUXILIARY: return "auxiliary"
    return "discovery"

def track(sid, pmid):
    if sid in TRACK_A: return "A", TRACK_A[sid]
    if pmid in TRACK_A_BY_PMID and TRACK_A_BY_PMID[pmid]: return "A", TRACK_A_BY_PMID[pmid]
    return "B", "author gene-level statistics (ORCS/publication)"

rows = []
for r in man:
    if not (r["host_species"] == "Homo sapiens" and r["genome_wide_yes_no"] == "yes"
            and r["whole_infection_yes_no"] == "yes"):
        continue
    sid = r["study_id"]
    a = adj.get(sid, {})
    p = pre.get(sid, {})
    rl = role(sid)
    if rl == "excluded":
        continue
    tr, tr_src = track(sid, r["pmid"])
    if rl == "holdout_locked":
        tr, tr_src = "LOCKED", "holdout: accession/role only until HG4"

    # effect metric + SE provenance
    if tr == "A":
        metric = "MAGeCK RRA score (recomputed from raw counts); L2FC as secondary"
        se_def = "RRA rank-based (no parametric SE); MLE beta-score SE where design supports"
        se_prov = "recomputed_this_project"
        meta = "conditional"
        meta_reason = "trackA_pending_QC (Stage 3)"
    elif tr == "B":
        metric = "author gene-level statistic (per ORCS score_summary; e.g. MAGeCK RRA / CasTLE / L2FC)"
        se_def = "as published; heterogeneous across studies"
        se_prov = "publication_stats_only"
        meta = "conditional"
        meta_reason = "author_stats_heterogeneous; usable for recurrence/permutation, formal meta only after harmonization"
    else:
        metric = "LOCKED"
        se_def = "LOCKED"
        se_prov = "LOCKED"
        meta = "no"
        meta_reason = "holdout_locked_until_HG4"

    # arms
    if a.get("endpoint_family") == "HOST_FITNESS":
        comp, ref = "infected, selection-surviving cell population", "pre-selection library (T0/plasmid) or uninfected control (per study)"
        sel = "survival under infection-induced cytopathic death"
    elif a.get("endpoint_family") == "PATHOGEN_BURDEN":
        comp, ref = "viral-reporter-low/negative FACS fraction", "reporter-high fraction or unsorted infected control"
        sel = "FACS on viral reporter/antigen"
    elif a.get("endpoint_family") == "PATHOGEN_STATE":
        comp, ref = "state-reporter-high FACS fraction (active-replication or persistence, per contrast)", "uninfected-baseline library and/or alternate-state fraction"
        sel = "FACS on bacterial phenotypic-state reporter"
    elif a.get("endpoint_family") == "TOXIN_SENSITIVITY":
        comp, ref = "survivors after 3 rounds of PA14 ExoU intoxication", "unselected starting population"
        sel = "survival under toxin-driven death"
    else:
        comp, ref = "selected population (per publication)", "reference population (per publication)"
        sel = a.get("endpoint_exact", "")

    # GSE210689 has two state contrasts -> emit two estimand rows
    contrasts = [("", "")]
    if sid == "LIT-Mtb-GSE210689":
        contrasts = [("active_replication", "active-replication reporter sort vs uninfected baseline"),
                     ("stressed_nonreplicating", "stressed/non-replicating reporter sort vs uninfected baseline")]

    for ctag, cdesc in contrasts:
        eid = f"EST-{sid.replace('LIT-', '').replace('ORCS-', 'ORCS')}{'_' + ctag if ctag else ''}"
        rows.append({
            "study_id": sid, "screen_id": sid, "pmid": r["pmid"], "pathogen": r["pathogen"],
            "cell_line": r["cell_line"], "cell_type": r["cell_type"] or "see Cellosaurus " + r["cellosaurus_id"],
            "perturbation_type": r["perturbation_type"],
            "comparison_arm": cdesc or comp, "reference_arm": ref if not cdesc else "uninfected baseline (shared within study)",
            "selection_scheme": sel, "endpoint_family": a.get("endpoint_family", ""),
            "endpoint_exact": a.get("endpoint_exact", ""), "timepoint": r["timepoint"] or "per publication (verify Track B)",
            "effect_metric": metric, "effect_sign_definition": a.get("raw_sign_meaning", ""),
            "SE_or_variance_definition": se_def, "SE_provenance": se_prov,
            "unit_of_analysis": "gene (guide-level effects aggregated per gene)",
            "biological_replicate_definition": r["biological_replicates"] or "per publication (verify Track B); GSE210689=3, Sunshine=2",
            "shared_control_cluster": r["independence_cluster"] if r["independence_cluster"] else r["pmid"],
            "independence_cluster": r["independence_cluster"],
            "estimand_id": eid, "meta_eligible": meta, "meta_exclusion_reason": meta_reason,
            "adjudication_status": a.get("adjudication_status", ""),
            "track": tr, "track_source": tr_src, "holdout_role": rl,
            "normalized_direction_enriched": a.get("normalized_direction_enriched", ""),
            "normalized_direction_depleted": a.get("normalized_direction_depleted", ""),
        })

fields = ["study_id", "screen_id", "pmid", "pathogen", "cell_line", "cell_type", "perturbation_type",
          "comparison_arm", "reference_arm", "selection_scheme", "endpoint_family", "endpoint_exact",
          "timepoint", "effect_metric", "effect_sign_definition", "SE_or_variance_definition",
          "SE_provenance", "unit_of_analysis", "biological_replicate_definition",
          "shared_control_cluster", "independence_cluster", "estimand_id", "meta_eligible",
          "meta_exclusion_reason", "adjudication_status", "track", "track_source", "holdout_role",
          "normalized_direction_enriched", "normalized_direction_depleted"]
with open(f"{OUT}/screen_estimand_registry.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

from collections import Counter
print("registry rows:", len(rows))
print("holdout_role:", dict(Counter(r["holdout_role"] for r in rows)))
print("track:", dict(Counter(r["track"] for r in rows)))
print("meta_eligible:", dict(Counter(r["meta_eligible"] for r in rows)))
print("discovery by pathogen:", dict(Counter(r["pathogen"] for r in rows if r["holdout_role"] == "discovery")))
h = hashlib.sha256(open(f"{OUT}/screen_estimand_registry.tsv", "rb").read()).hexdigest()
print("sha256:", h)
