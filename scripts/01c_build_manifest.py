#!/usr/bin/env python3
"""HPICDA Stage 0 — build data_availability_manifest.tsv.

Combines ORCS 2.0.18 census output (02_metadata/core_screens_classified.json)
with literature-derived screens not yet curated in ORCS. Fills the blueprint
section 8 manifest fields. Screen-level metadata only; no gene-level results.

Usage:
    python scripts/01c_build_manifest.py --out 01_stage0_audit/data_availability_manifest.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

# --- classification rules -------------------------------------------------
KNOWN_GW_LIBS = [
    "brunello", "gecko", "tkov3", "avana", "gattinara", "calabrese",
    "dolcetto", "h1 crispri", "edit-r human whole genome", "minlibcas9",
    "liu human crispr knockout", "sabatini", "whitehead", "crisprn (spaan",
    "crisprn (cp1560)", "toronto knockout", "bassik human crispr",
]
TARGETED_HINTS = [
    "targeted", "interactome", "interactor", "surfaceome", "isg", "focused",
    "custom", "414 interferon",
]
# Cell lines whose host species is NOT human despite a human library.
NONHUMAN_CELL_LINES = {"vero-e6", "vero", "vero c1008", "vero-e6|vero c1008"}

CELLOSAURUS = {
    "A-549": "CVCL_0023", "A549-ACE2": "", "A549-AC (hACE2 and Cas9)": "",
    "A549-SunTag ACE2": "", "A549-deltaSTAT1-SunTag": "",
    "A549-ACE2-Cas9-GFP1–10": "",
    "THP-1": "CVCL_0006", "U-937": "CVCL_0007", "Calu-3": "CVCL_0609",
    "Calu1-ACE2": "", "HEp-2": "CVCL_1906", "HeLa": "CVCL_0030",
    "HuH-7": "CVCL_0336", "Huh-7.5": "CVCL_7927", "Huh-7.5.1": "CVCL_E049",
    "Huh-7.5.1-ACE2": "", "HEK293T/ACE2": "", "HEK293T-CRISPRa": "CVCL_0063",
    "HEK293-ACE2-TMPRSS2": "", "HEK-293": "CVCL_0045", "HEK293SF": "CVCL_4U53",
    "Caco-2-ACE2": "", "C2BBe1": "CVCL_1096", "HCT 116": "CVCL_0291",
    "A-375": "CVCL_0132", "UM-UC-4": "CVCL_2756", "HT-29": "CVCL_0320",
    "HAP-1": "CVCL_Y019", "K-562": "CVCL_0004", "SH-SY5Y": "CVCL_0019",
    "PLC/PRF/5": "CVCL_0485", "Vero-E6|VERO C1008": "CVCL_0574",
}

PATHOGEN_CLASS = {
    "SARS-CoV-2": "viral", "Influenza": "viral", "RSV": "viral",
    "Mtb/BCG": "bacterial", "Mtb": "bacterial", "Pseudomonas": "bacterial",
    "S. aureus (PVL toxin)": "bacterial", "S. aureus (aHL toxin)": "bacterial",
    "Legionella (FIREWALLED)": "bacterial_excluded",
}

FIELDS = [
    "study_id", "pmid", "publication_year", "pathogen", "strain",
    "host_species", "cell_line", "cellosaurus_id", "cell_type",
    "perturbation_type", "library", "library_size", "genome_wide_yes_no",
    "endpoint_family", "endpoint_exact", "infection_stage", "timepoint",
    "selection_strategy", "raw_counts_available", "gene_stats_available",
    "supplement_available", "uninfected_baseline_available",
    "biological_replicates", "shared_control_group", "data_access_status",
    "license", "discovery_or_holdout", "notes", "pathogen_class",
    "whole_infection_yes_no", "toxin_only_yes_no", "independence_cluster",
]


def lib_name(s: str) -> str:
    m = re.search(r"Library Name:\s*(.+?)\s*Methodology:", s)
    return m.group(1).strip() if m else s[:60]


def parse_score(s: str):
    hits = re.search(r"Number of Hits:\s*([\d,]+)", s)
    full = re.search(r"Full Dataset Size:\s*([\d,]+)", s)
    f = lambda m: int(m.group(1).replace(",", "")) if m else None  # noqa: E731
    return f(hits), f(full)


def classify_pathogen(cond_name: str):
    c = cond_name.lower()
    if "sars-cov-2" in c:
        return "SARS-CoV-2"
    if "influenza" in c:
        return "Influenza"
    if "respiratory syncytial" in c:
        return "RSV"
    if "tuberculosis" in c or "bovis bcg" in c:
        return "Mtb/BCG"
    if "panton" in c:
        return "S. aureus (PVL toxin)"
    if "legionella" in c:
        return "Legionella (FIREWALLED)"
    return None


def gw_status(lib: str, full):
    """Blueprint section 8: author-declared genome-wide OR >=15,000
    protein-coding genes OR named genome-wide library. A named library
    whose deposited dataset is a small subset does NOT qualify."""
    lib_l = lib.lower()
    named = any(k in lib_l for k in KNOWN_GW_LIBS)
    targeted = any(k in lib_l for k in TARGETED_HINTS)
    if full is not None and full >= 15000:
        return "yes"
    if full is not None and full < 10000:
        return "no"  # small deposited dataset: targeted/subscreen
    if targeted:
        return "no"
    if named:
        return "yes"  # named GW library with 10k-15k or unknown size
    return "no"


def orcs_rows(classified: list[dict]) -> list[dict]:
    out = []
    for o in classified:
        cell = o["cell_line"]
        host = "Homo sapiens"
        if cell.lower() in NONHUMAN_CELL_LINES or "vero" in cell.lower():
            host = "Chlorocebus sabaeus (human library)"
        cond = o["condition"]
        whole = o["whole_infection"]
        toxin = (not whole) and ("toxin" in cond.lower() or "pvl" in cond.lower()
                                 or "spike" in cond.lower() or "pseudovirus" in cond.lower())
        year = ""
        m = re.search(r"\((\d{4})\)", o["dataset_ref"])
        if m:
            year = m.group(1)
        out.append({
            "study_id": f"ORCS-{o['orcs_screen_id']}",
            "pmid": o["pmid"],
            "publication_year": year,
            "pathogen": o["pathogen"].replace(" (PVL toxin)", "").replace(" (FIREWALLED)", ""),
            "strain": re.sub(r"^Virus:\s*", "", cond)[:60],
            "host_species": host,
            "cell_line": cell,
            "cellosaurus_id": CELLOSAURUS.get(cell, ""),
            "cell_type": "",
            "perturbation_type": "CRISPRn" if "CRISPRn" in o["library"] or "Knockout" in o["library"] else (
                "CRISPRa" if "Activation" in o["library"] or "CRISPRa" in o["library"] else (
                    "CRISPRi" if "CRISPRi" in o["library"] or "Inhibition" in o["library"] else "CRISPR")),
            "library": o["library"],
            "library_size": o["full_size"] if o["full_size"] else "",
            "genome_wide_yes_no": gw_status(o["library"], o["full_size"]),
            "endpoint_family": "viral_replication_or_survival" if whole else "entry_or_toxin",
            "endpoint_exact": cond[:80],
            "infection_stage": "",
            "timepoint": "",
            "selection_strategy": "",
            "raw_counts_available": "verify_at_stage2",
            "gene_stats_available": "yes_orcs",
            "supplement_available": "verify_at_stage2",
            "uninfected_baseline_available": "verify_at_stage2",
            "biological_replicates": "",
            "shared_control_group": "",
            "data_access_status": "public_orcs",
            "license": "CC-BY-4.0 (BioGRID ORCS)",
            "discovery_or_holdout": "TBD_HG2",
            "notes": f"ORCS 2.0.18 screen {o['orcs_screen_id']}; hits={o['n_hits']}",
            "pathogen_class": PATHOGEN_CLASS.get(o["pathogen"], ""),
            "whole_infection_yes_no": "yes" if whole else "no",
            "toxin_only_yes_no": "yes" if toxin else "no",
            "independence_cluster": f"PMID{o['pmid']}" if o["pmid"] else f"ORCS{o['orcs_screen_id']}",
        })
    return out


def literature_rows() -> list[dict]:
    """Screens identified by literature search, not (yet) in ORCS 2.0.18."""
    base = {f: "" for f in FIELDS}
    rows = []

    def add(**kw):
        r = dict(base)
        r.update(kw)
        rows.append(r)

    add(study_id="LIT-Mtb-GSE210689", pmid="37269834", publication_year="2023",
        pathogen="Mtb/BCG", strain="M. tuberculosis (reporter strain, Erdman background)",
        host_species="Homo sapiens", cell_line="THP-1", cellosaurus_id="CVCL_0006",
        cell_type="monocyte/macrophage", perturbation_type="CRISPRn",
        library="genome-wide CRISPR (author-declared; name TBD from publication)",
        library_size="", genome_wide_yes_no="yes",
        endpoint_family="bacterial_phenotypic_state",
        endpoint_exact="FACS on Mtb reporter: active replication vs stressed non-replicating",
        infection_stage="active_replication + latency_induction", timepoint="",
        selection_strategy="FACS_sort_reporter",
        raw_counts_available="yes", gene_stats_available="yes_publication",
        supplement_available="yes", uninfected_baseline_available="yes",
        biological_replicates="3", shared_control_group="yes_uninfected_baseline",
        data_access_status="public", license="GEO (public)",
        discovery_or_holdout="TBD_HG2",
        notes="GSE210689; raw count matrix GSE210689_RawCountMatrix.txt.gz + SRA FASTQ; "
              "9 samples (3 uninfected, 3 active, 3 stressed); MMGT1/GPR156 study",
        pathogen_class="bacterial", whole_infection_yes_no="yes", toxin_only_yes_no="no",
        independence_cluster="PMID37269834")

    add(study_id="LIT-RSV-PRJNA1281304-viability", pmid="", publication_year="2025",
        pathogen="RSV", strain="RSV A2-EGFP", host_species="Homo sapiens",
        cell_line="A-549", cellosaurus_id="CVCL_0023", cell_type="lung epithelial",
        perturbation_type="CRISPRn", library="GeCKO v2", library_size="",
        genome_wide_yes_no="yes", endpoint_family="viral_replication_or_survival",
        endpoint_exact="EGFP-negative sorted, 1 wk outgrowth (viability screen)",
        infection_stage="early+outgrowth", timepoint="24h sort + 1wk",
        selection_strategy="FACS_sort_EGFP_neg",
        raw_counts_available="yes_sra", gene_stats_available="yes_publication",
        supplement_available="yes", uninfected_baseline_available="yes_unsorted_control",
        biological_replicates="2", shared_control_group="",
        data_access_status="public", license="SRA (public)",
        discovery_or_holdout="TBD_HG2",
        notes="Sunshine et al. mBio 2025 (10.1128/mbio.01010-25); PRJNA1281304; "
              "scRNA GSE281623; code github.com/Sara-Sunshine/RSV_scRNAseq_Cas9screen_2025",
        pathogen_class="viral", whole_infection_yes_no="yes", toxin_only_yes_no="no",
        independence_cluster="Sunshine2025")

    add(study_id="LIT-RSV-PRJNA1281304-early", pmid="", publication_year="2025",
        pathogen="RSV", strain="RSV A2-EGFP", host_species="Homo sapiens",
        cell_line="A-549", cellosaurus_id="CVCL_0023", cell_type="lung epithelial",
        perturbation_type="CRISPRn", library="GeCKO v2", library_size="",
        genome_wide_yes_no="yes", endpoint_family="viral_replication_or_survival",
        endpoint_exact="EGFP-negative/low sorted 24 hpi (early infection screen)",
        infection_stage="early", timepoint="24h",
        selection_strategy="FACS_sort_EGFP_neg_low",
        raw_counts_available="yes_sra", gene_stats_available="yes_publication",
        supplement_available="yes", uninfected_baseline_available="yes_unsorted_control",
        biological_replicates="2", shared_control_group="yes_shared_with_viability",
        data_access_status="public", license="SRA (public)",
        discovery_or_holdout="TBD_HG2",
        notes="Sunshine et al. mBio 2025; same study as viability screen (shared control)",
        pathogen_class="viral", whole_infection_yes_no="yes", toxin_only_yes_no="no",
        independence_cluster="Sunshine2025")

    add(study_id="LIT-Pa-GSE154751", pmid="", publication_year="2021",
        pathogen="Pseudomonas", strain="P. aeruginosa PA14 (ExoU+)",
        host_species="Homo sapiens", cell_line="A-549", cellosaurus_id="CVCL_0023",
        cell_type="alveolar epithelial", perturbation_type="CRISPRn",
        library="one-component CRISPR-Cas9 library (71,090 gRNAs / 18,053 genes)",
        library_size="18053", genome_wide_yes_no="yes",
        endpoint_family="toxin_mediated_cytotoxicity",
        endpoint_exact="survival after PA14 infection (ExoU-dependent cell death), 3 rounds",
        infection_stage="acute_cytotoxicity", timepoint="90min x3 rounds",
        selection_strategy="survival",
        raw_counts_available="yes_sra", gene_stats_available="verify_at_stage2",
        supplement_available="verify_at_stage2", uninfected_baseline_available="verify_at_stage2",
        biological_replicates="", shared_control_group="",
        data_access_status="public", license="GEO (public)",
        discovery_or_holdout="TBD_HG2",
        notes="GSE154751; whole-bacteria infection but selection is ExoU-toxin-driven death; "
              "borderline whole-infection (flagged for HG2 adjudication)",
        pathogen_class="bacterial", whole_infection_yes_no="yes", toxin_only_yes_no="no",
        independence_cluster="GSE154751")

    add(study_id="LIT-Pa-SRP514827", pmid="", publication_year="2024",
        pathogen="Pseudomonas", strain="P. aeruginosa (T3SS+)",
        host_species="Homo sapiens", cell_line="HAP-1", cellosaurus_id="CVCL_Y019",
        cell_type="haploid fibroblast-like", perturbation_type="gene_trap_mutagenesis",
        library="TnPuro gene-trap (haploid)", library_size="",
        genome_wide_yes_no="yes",
        endpoint_family="effector_injection",
        endpoint_exact="resistance to T3SS effector injection (FACS)",
        infection_stage="intoxication", timepoint="",
        selection_strategy="FACS_sort_T3SS_resistant",
        raw_counts_available="yes_sra", gene_stats_available="verify_at_stage2",
        supplement_available="verify_at_stage2", uninfected_baseline_available="verify_at_stage2",
        biological_replicates="2", shared_control_group="",
        data_access_status="public", license="SRA (public)",
        discovery_or_holdout="TBD_HG2",
        notes="SRP514827; NOT CRISPR (gene-trap); endpoint is T3SS injection not replication",
        pathogen_class="bacterial", whole_infection_yes_no="yes", toxin_only_yes_no="no",
        independence_cluster="SRP514827")

    add(study_id="LIT-Sa-aHL-srep24242", pmid="27052225", publication_year="2016",
        pathogen="S. aureus (aHL toxin)", strain="alpha-hemolysin (purified toxin)",
        host_species="Homo sapiens", cell_line="U-937", cellosaurus_id="CVCL_0007",
        cell_type="myeloid", perturbation_type="CRISPRn",
        library="GeCKO v1 (19,050 genes + 1,864 miRNAs)", library_size="19050",
        genome_wide_yes_no="yes", endpoint_family="toxin_mediated_cytotoxicity",
        endpoint_exact="survival after alpha-hemolysin intoxication",
        infection_stage="intoxication", timepoint="",
        selection_strategy="survival",
        raw_counts_available="verify_at_stage2", gene_stats_available="yes_publication",
        supplement_available="yes", uninfected_baseline_available="n/a_toxin",
        biological_replicates="", shared_control_group="",
        data_access_status="public", license="journal (Sci Rep)",
        discovery_or_holdout="TBD_HG2",
        notes="Sci Rep 2016 srep24242; TOXIN-ONLY (no live bacteria); ADAM10/SYS1/ARFRP1/TSPAN14",
        pathogen_class="bacterial", whole_infection_yes_no="no", toxin_only_yes_no="yes",
        independence_cluster="PMID27052225")

    # Murine Mtb screens — recorded for the bacterial-arm audit, NOT eligible
    # for the human host atlas.
    add(study_id="LIT-Mtb-GSE267063-murine", pmid="", publication_year="2024",
        pathogen="Mtb/BCG", strain="M. tuberculosis", host_species="Mus musculus",
        cell_line="Hoxb8-immortalized macrophages", cellosaurus_id="",
        cell_type="macrophage", perturbation_type="CRISPRn",
        library="genome-wide (murine)", library_size="", genome_wide_yes_no="yes",
        endpoint_family="bacterial_growth_restriction",
        endpoint_exact="FACS on Mtb-induced cytotoxicity / intracellular growth",
        infection_stage="intracellular_growth", timepoint="",
        selection_strategy="FACS_sort",
        raw_counts_available="yes", gene_stats_available="yes_publication",
        supplement_available="yes", uninfected_baseline_available="verify_at_stage2",
        biological_replicates="", shared_control_group="",
        data_access_status="public", license="GEO (public)",
        discovery_or_holdout="EXCLUDED_nonhuman",
        notes="GSE267063/PRJNA1109426; Simwela et al. Nat Commun 2024; GID/CTLH; "
              "MURINE host - excluded from human atlas, listed for audit completeness",
        pathogen_class="bacterial", whole_infection_yes_no="yes", toxin_only_yes_no="no",
        independence_cluster="Simwela2024")

    add(study_id="LIT-Mtb-EMTAB15739-murine", pmid="", publication_year="2024",
        pathogen="Mtb/BCG", strain="M. tuberculosis", host_species="Mus musculus",
        cell_line="conditionally-immortalized macrophages", cellosaurus_id="",
        cell_type="macrophage", perturbation_type="CRISPRn",
        library="Brie (murine genome-wide)", library_size="", genome_wide_yes_no="yes",
        endpoint_family="host_inflammatory_response",
        endpoint_exact="TNF/iNOS regulator expression in infected macrophages",
        infection_stage="innate_response", timepoint="",
        selection_strategy="FACS_sort",
        raw_counts_available="verify_at_stage2", gene_stats_available="verify_at_stage2",
        supplement_available="verify_at_stage2", uninfected_baseline_available="verify_at_stage2",
        biological_replicates="", shared_control_group="",
        data_access_status="public", license="ArrayExpress (public)",
        discovery_or_holdout="EXCLUDED_nonhuman",
        notes="E-MTAB-15739/PRJEB100798; MURINE host - excluded, listed for audit completeness",
        pathogen_class="bacterial", whole_infection_yes_no="yes", toxin_only_yes_no="no",
        independence_cluster="EMTAB15739")

    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--classified", default="02_metadata/core_screens_classified.json")
    ap.add_argument("--out", default="01_stage0_audit/data_availability_manifest.tsv")
    args = ap.parse_args()

    classified = json.loads(Path(args.classified).read_text())
    rows = orcs_rows(classified) + literature_rows()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Wrote {out} with {len(rows)} screens")


if __name__ == "__main__":
    main()
