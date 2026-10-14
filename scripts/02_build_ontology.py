#!/usr/bin/env python3
"""HPICDA Stage 1: freeze unified ontology + per-screen endpoint/direction adjudication.

Blueprint v1.1 §9/§9.1. Reads Stage 0 manifest + ORCS census metadata; emits:
  02_ontology/hpicda_ontology_v1.tsv                      (controlled vocabulary)
  02_ontology/screen_endpoint_direction_adjudication.tsv  (§9.1 per-screen sign-off)
  02_ontology/stage1_manifest_corrections.tsv             (errata found during Stage 1)

Direction semantics are normalized to loss-of-function (KO) equivalents; CRISPRa/CRISPRi
raw signs are inverted and flagged. Screens whose direction cannot be reliably interpreted
are marked DIRECTION_UNRESOLVED (never force-classified).
"""
import csv, json, hashlib, datetime, os, re

ROOT = "/workspace/HPICDA"
MANIFEST = f"{ROOT}/01_stage0_audit/data_availability_manifest.tsv"
CENSUS = f"{ROOT}/02_metadata/orcs_census_screens.tsv"
OUTDIR = f"{ROOT}/02_ontology"
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------- ontology ---
ONTOLOGY = [
    # term_type, term, definition, pooling_rule
    ("endpoint_family", "ENTRY", "Pathogen entry/attachment into host cell (incl. pseudovirus/RBD binding)", "never pool with non-ENTRY"),
    ("endpoint_family", "INTRACELLULAR_REPLICATION", "Pathogen replication inside host cell measured directly", "pool only within same family+metric"),
    ("endpoint_family", "PATHOGEN_BURDEN", "Pathogen amount/activity via antigen, reporter, or molecular readout (incl. PRF reporter)", "pool only within same family+metric"),
    ("endpoint_family", "HOST_FITNESS", "Host survival/proliferation under infection (CPE-death protection)", "pool only within same family+metric"),
    ("endpoint_family", "HOST_DAMAGE", "Host damage readout short of death (e.g. membrane damage)", "pool only within same family+metric"),
    ("endpoint_family", "IMMUNE_OUTPUT", "Immune signaling/cytokine output of infected host cells", "pool only within same family+metric"),
    ("endpoint_family", "PATHOGEN_STATE", "Pathogen phenotypic state (e.g. Mtb active replication vs persistence reporter)", "never pool with fitness/burden families"),
    ("endpoint_family", "TOXIN_SENSITIVITY", "Sensitivity to a defined pathogen toxin/effector (toxin-driven death)", "never pool with whole-infection families"),
    ("endpoint_family", "OTHER_DEFINED", "Any other precisely defined endpoint (e.g. T3SS effector injection)", "case-by-case only"),
    ("endpoint_family", "UNRESOLVED", "Endpoint cannot be reliably interpreted", "excluded from pooling"),
    ("direction", "KO_REDUCES_PATHOGEN", "Loss of gene reduces pathogen burden/replication/function (pro-pathogen dependency factor)", "LOF-normalized"),
    ("direction", "KO_INCREASES_PATHOGEN", "Loss of gene increases pathogen burden/replication (restriction factor)", "LOF-normalized"),
    ("direction", "KO_PROTECTS_HOST", "Loss of gene protects host from infection-induced death/damage", "LOF-normalized"),
    ("direction", "KO_SENSITIZES_HOST", "Loss of gene sensitizes host to infection-induced death/damage", "LOF-normalized"),
    ("direction", "KO_REDUCES_HOST_DAMAGE", "Loss of gene reduces host damage readout", "LOF-normalized"),
    ("direction", "KO_INCREASES_HOST_DAMAGE", "Loss of gene increases host damage readout", "LOF-normalized"),
    ("direction", "KO_REDUCES_IMMUNE_OUTPUT", "Loss of gene reduces immune output", "LOF-normalized"),
    ("direction", "KO_INCREASES_IMMUNE_OUTPUT", "Loss of gene increases immune output", "LOF-normalized"),
    ("direction", "KO_ALTERS_PATHOGEN_STATE", "Loss of gene shifts pathogen phenotypic state (direction of shift recorded per screen)", "LOF-normalized"),
    ("direction", "DIRECTION_UNRESOLVED", "Direction cannot be reliably interpreted", "excluded from directional synthesis"),
    ("hdt_direction", "INHIBIT_TARGET_BENEFICIAL", "Inhibiting the host target is the candidate-beneficial direction", "requires direction adjudication"),
    ("hdt_direction", "PRESERVE_OR_ACTIVATE_TARGET_BENEFICIAL", "Preserving/activating the host target is the candidate-beneficial direction", "requires direction adjudication"),
    ("hdt_direction", "INTERVENTION_DIRECTION_UNRESOLVED", "Intervention direction cannot be assigned from available evidence", "default until adjudicated"),
]

# ------------------------------------------------- per-study adjudication ----
# Keyed by (pmid, pathogen-group). Designs grounded in Stage 0 audits + ORCS
# phenotype/condition metadata + paper titles (PMID 37269834 verified by fetch).
SURV = dict(  # viral survival / CPE-protection screen (CRISPRn)
    endpoint_family="HOST_FITNESS",
    raw_sign_meaning="enriched guides = KO cells surviving virus-induced cytopathic death; depleted = KO cells sensitized",
    normalized_direction_enriched="KO_PROTECTS_HOST",
    normalized_direction_depleted="KO_SENSITIZES_HOST",
    implied_pathogen_effect_enriched="KO_REDUCES_PATHOGEN",
    hdt_interpretability="direct: enriched-tail dependencies are INHIBIT candidates; depleted-tail restriction factors are PRESERVE/ACTIVATE candidates",
    perturbation_sign_inversion="no",
)
CRISPRA_SURV = dict(
    endpoint_family="HOST_FITNESS",
    raw_sign_meaning="enriched guides = CRISPRa overexpression cells surviving virus-induced death (gain-of-function protection)",
    normalized_direction_enriched="KO_SENSITIZES_HOST",   # LOF-equivalent: gene is restriction factor
    normalized_direction_depleted="KO_PROTECTS_HOST",
    implied_pathogen_effect_enriched="KO_INCREASES_PATHOGEN",
    hdt_interpretability="inverted: CRISPRa-enriched genes are restriction factors -> PRESERVE/ACTIVATE candidates",
    perturbation_sign_inversion="yes",
)
FACS_VIRAL_REPORTER = dict(
    endpoint_family="PATHOGEN_BURDEN",
    raw_sign_meaning="enriched guides = KO cells sorted into viral-reporter-low/negative fraction (reduced viral replication)",
    normalized_direction_enriched="KO_REDUCES_PATHOGEN",
    normalized_direction_depleted="KO_INCREASES_PATHOGEN",
    implied_pathogen_effect_enriched="KO_REDUCES_PATHOGEN",
    hdt_interpretability="direct: enriched-tail are INHIBIT candidates",
    perturbation_sign_inversion="no",
)
MTB_STATE = dict(
    endpoint_family="PATHOGEN_STATE",
    raw_sign_meaning=("FACS sort of infected THP-1 by Mtb phenotypic-state reporter (active replication vs "
                      "stressed/non-replicating persistence); enriched guides shift bacterial state distribution "
                      "(PMID 37269834, Cell Host Microbe 2023: host regulators of Mtb phenotypes; MMGT1-GPR156 axis)"),
    normalized_direction_enriched="KO_ALTERS_PATHOGEN_STATE",
    normalized_direction_depleted="KO_ALTERS_PATHOGEN_STATE",
    implied_pathogen_effect_enriched="state-shift (active-vs-persistence direction recorded per contrast at Stage 2)",
    hdt_interpretability="complex: state-shift is not burden reduction; HDT direction requires Stage 12 adjudication",
    perturbation_sign_inversion="no",
)
BCG_SURV = dict(
    endpoint_family="HOST_FITNESS",
    raw_sign_meaning="enriched guides = KO cells protected from BCG-infection-induced host cell death (THP-1)",
    normalized_direction_enriched="KO_PROTECTS_HOST",
    normalized_direction_depleted="KO_SENSITIZES_HOST",
    implied_pathogen_effect_enriched="host-death suppression (BCG near-neighbor context; NOT virulent Mtb)",
    hdt_interpretability="direct but near-neighbor: interpret as Mtb-complex context only",
    perturbation_sign_inversion="no",
)
TOXIN = dict(
    endpoint_family="TOXIN_SENSITIVITY",
    raw_sign_meaning="enriched guides = KO cells surviving ExoU-toxin-driven acute cytotoxicity (live PA14, 3 selection rounds)",
    normalized_direction_enriched="KO_REDUCES_HOST_DAMAGE",
    normalized_direction_depleted="KO_INCREASES_HOST_DAMAGE",
    implied_pathogen_effect_enriched="toxin-context only; not whole-infection fitness",
    hdt_interpretability="toxin-context only; excluded from whole-infection synthesis (Stage 0 borderline ruling)",
    perturbation_sign_inversion="no",
)
T3SS = dict(
    endpoint_family="OTHER_DEFINED",
    raw_sign_meaning="enriched guides = gene-trap-disrupted HAP-1 cells resistant to T3SS effector injection (BlaM reporter FACS)",
    normalized_direction_enriched="KO_REDUCES_PATHOGEN",
    normalized_direction_depleted="KO_INCREASES_PATHOGEN",
    implied_pathogen_effect_enriched="reduced T3SS injection function",
    hdt_interpretability="entry/function-context only; gene-trap (not CRISPR); excluded from whole-infection synthesis",
    perturbation_sign_inversion="no",
)

# (pmid, pathogen_group) -> (rule, endpoint_exact, notes)
STUDY_RULES = {
    ("28813663", "Influenza"):   (CRISPRA_SURV, "IAV PR8 H1N1-induced host death; CRISPRa survival", "Heaton 2017"),
    ("29642015", "Influenza"):   (SURV, "IAV H5N1-induced host death; GeCKOv2 survival", "Han 2018; GSE111166 raw counts"),
    ("31919360", "Influenza"):   (SURV, "IAV PR8 H1N1-induced host death; Avana-4 survival", "Li 2020; Edinburgh raw counts"),
    ("32699298", "Influenza"):   (SURV, "IAV PR8 H1N1-induced host death; Brunello survival in HEK293SF", "Sharon 2020"),
    ("36309510", "Influenza"):   (CRISPRA_SURV, "Influenza B-induced host death; Calabrese CRISPRa survival", "Trimarco 2022; IBV near-neighbor - auxiliary for IAV replication counts"),
    ("40372436", "Influenza"):   (SURV, "IAV PR8 H1N1-induced host death; Brunello survival", "Soriaga 2025 multi-virus study"),
    ("40372436", "RSV"):         (SURV, "RSV A2-induced host death; Brunello survival in HEp-2 (15 MOI)", "Soriaga 2025"),
    ("33147445", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; GeCKOv2 survival in A549", "Wang 2020"),
    ("33333024", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; GeCKOv2 survival in Huh-7.5.1", "Daniloski 2021"),
    ("33382968", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; Brunello survival in Huh-7.5", "Wei 2021"),
    ("33686287", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; Brunello survival in HuH-7 (multi-condition)", "Baggen 2021"),
    ("35086559", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; Brunello survival, MOI series in HEK293T/ACE2", "Biering 2022"),
    ("35469023", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; Brunello survival in Calu-3 incl. Alpha/Beta VOCs", "GATA6 study 2022"),
    ("35879412", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; Brunello KO + Calabrese CRISPRa bidirectional in Calu-3", "Rebendenne 2022a; CRISPRa arms sign-inverted per row perturbation_type"),
    ("35879413", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; Gattinara KO + Calabrese CRISPRa bidirectional", "Rebendenne 2022b Nat Genet"),
    ("36597481", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; TKOv3 survival across 4 human backgrounds", "Prelli Bozzo 2023; GSE177545 raw counts"),
    ("36753415", "SARS-CoV-2"):  (FACS_VIRAL_REPORTER, "SARS-CoV-2 programmed ribosomal frameshifting (PRF) reporter in HCT 116", "PRF reporter = viral molecular activity; PATHOGEN_BURDEN not HOST_FITNESS"),
    ("39316623", "SARS-CoV-2"):  (SURV, "SARS-CoV-2-induced host death; Edit-R crRNA survival in Huh-7.5 (+/- IFN condition arm)", "2024 temporal holdout candidate; deposited 12,119-gene subset of named GW library"),
    ("37269834", "Mtb/BCG"):     (MTB_STATE, "Mtb phenotypic-state reporter (active replication vs stressed persistence) in THP-1; FACS sort", "GSE210689; sole strict bacterial discovery screen"),
    ("32970993", "Mtb/BCG"):     (BCG_SURV, "BCG-infection-induced host death; Brunello CRISPRn + Dolcetto CRISPRi in THP-1", "near-neighbor auxiliary; CRISPRi arm sign-equivalent to KO for direction"),
}

def adjudicate_row(r, census_by_id):
    """Return adjudication dict for one manifest row."""
    sid = r["study_id"]
    pmid, path = r["pmid"], r["pathogen"]
    m = re.search(r"ORCS-(\d+)", sid)
    orcs_id = m.group(1) if m else ""
    cen = census_by_id.get(orcs_id, {})

    # --- exclusion: HPIV-3 misclassified as Influenza ---
    if orcs_id == "2433" or "parainfluenza" in cen.get("condition", "").lower():
        return dict(endpoint_family="EXCLUDED", endpoint_exact="Human parainfluenza 3 virus (HPIV-3) screen",
                    raw_sign_meaning="", normalized_direction_enriched="", normalized_direction_depleted="",
                    implied_pathogen_effect_enriched="", hdt_interpretability="",
                    perturbation_sign_inversion="", adjudication_status="EXCLUDED_NON_CORE_PATHOGEN",
                    adjudication_note="ORCS-2433 condition=HPIV-3 (PMID 40372436 multi-virus study); HPIV-3 is not an HPICDA core pathogen; Stage 0 manifest misclassified as Influenza - corrected here")

    # --- literature (non-ORCS) screens ---
    if sid == "LIT-Mtb-GSE210689":
        return dict(MTB_STATE, endpoint_exact=STUDY_RULES[("37269834", "Mtb/BCG")][1],
                    adjudication_status="ADJUDICATED",
                    adjudication_note="PMID 37269834 verified: host regulators of Mtb phenotypes (state reporter)")
    if sid.startswith("LIT-RSV-PRJNA1281304"):
        if sid.endswith("viability"):
            return dict(SURV, endpoint_exact="RSV A2-EGFP-induced host death; GeCKOv2 survival in A549",
                        adjudication_status="ADJUDICATED", adjudication_note="Sunshine 2025 mBio; viability arm")
        return dict(FACS_VIRAL_REPORTER, endpoint_exact="RSV A2-EGFP early infection; FACS EGFP-negative sort in A549",
                    adjudication_status="ADJUDICATED", adjudication_note="Sunshine 2025 mBio; early-infection reporter arm")
    if sid == "LIT-Pa-GSE154751":
        return dict(TOXIN, endpoint_exact="P. aeruginosa PA14 ExoU-toxin acute cytotoxicity; survival selection in A549",
                    adjudication_status="ADJUDICATED_BORDERLINE", adjudication_note="Stage 0 borderline: toxin endpoint; context role only")
    if sid == "LIT-Pa-SRP514827":
        return dict(T3SS, endpoint_exact="P. aeruginosa T3SS effector injection resistance; HAP-1 gene-trap FACS",
                    adjudication_status="ADJUDICATED_BORDERLINE", adjudication_note="gene-trap not CRISPR; context role only")

    # --- ORCS screens via study rules ---
    key = (pmid, path)
    if key not in STUDY_RULES:
        return dict(endpoint_family="UNRESOLVED", endpoint_exact="", raw_sign_meaning="",
                    normalized_direction_enriched="", normalized_direction_depleted="",
                    implied_pathogen_effect_enriched="", hdt_interpretability="",
                    perturbation_sign_inversion="", adjudication_status="DIRECTION_UNRESOLVED",
                    adjudication_note=f"no study rule for pmid={pmid} pathogen={path}")
    rule, exact, note = STUDY_RULES[key]
    out = dict(rule)
    out["endpoint_exact"] = exact
    out["adjudication_note"] = note
    # per-row CRISPRa arm inversion (Rebendenne studies have mixed KO + CRISPRa arms)
    pert = r.get("perturbation_type", "")
    if "CRISPRa" in pert or "activation" in r.get("library", "").lower():
        if rule is SURV:
            out = dict(CRISPRA_SURV)
            out["endpoint_exact"] = exact + " [CRISPRa arm]"
            out["adjudication_note"] = note + "; row-level CRISPRa inversion applied"
    out["adjudication_status"] = "ADJUDICATED"
    return out

def main():
    man = list(csv.DictReader(open(MANIFEST), delimiter="\t"))
    cen = {r["orcs_screen_id"]: r for r in csv.DictReader(open(CENSUS), delimiter="\t")}
    wi = [r for r in man if r["host_species"] == "Homo sapiens" and r["genome_wide_yes_no"] == "yes"
          and r["whole_infection_yes_no"] == "yes"]

    # 1. ontology vocabulary
    with open(f"{OUTDIR}/hpicda_ontology_v1.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["term_type", "term", "definition", "pooling_rule", "protocol_version", "frozen_date"])
        for t in ONTOLOGY:
            w.writerow([*t, "1.1", "2026-08-30"])

    # 2. per-screen adjudication
    fields = ["study_id", "pmid", "pathogen", "cell_line", "perturbation_type", "library",
              "endpoint_family", "endpoint_exact", "raw_sign_meaning",
              "normalized_direction_enriched", "normalized_direction_depleted",
              "implied_pathogen_effect_enriched", "hdt_interpretability",
              "perturbation_sign_inversion", "adjudication_status", "adjudication_note",
              "adjudicator", "adjudication_timestamp"]
    rows_out = []
    for r in wi:
        adj = adjudicate_row(r, cen)
        rows_out.append({
            "study_id": r["study_id"], "pmid": r["pmid"], "pathogen": r["pathogen"],
            "cell_line": r["cell_line"], "perturbation_type": r["perturbation_type"],
            "library": r["library"], **adj,
            "adjudicator": "Biomni agent draft; ratified by delegation (Yongwei Qin, HG1 batch pre-authorization 2026-08-30)",
            "adjudication_timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        })
    with open(f"{OUTDIR}/screen_endpoint_direction_adjudication.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    # 3. errata
    corr = [
        dict(item="ORCS-2433", severity="EXCLUSION",
             finding="condition=Human parainfluenza 3 virus; manifest pathogen=Influenza is a misclassification",
             action="Excluded from discovery corpus (non-core pathogen). Eligible corpus 53->52 screens; Influenza 10->9; SARS-CoV-2 share 35/52=67.3%"),
        dict(item="HO-V-SC2-02 (holdout_inventory.tsv)", severity="ATTRIBUTION",
             finding="ORCS-2076/2077 (HCT 116) attributed to PMID 36597481; ORCS metadata shows PMID 36753415 and phenotype=PRF reporter",
             action="Holdout text corrected: HCT 116 PRF screens = PMID 36753415; UM-UC-4 screens (ORCS-2068/2069) remain PMID 36597481"),
        dict(item="ORCS-2020 (PMID 36309510)", severity="SPECIES_NOTE",
             finding="condition=Influenza B virus (B/Yamagata/16/1988), not Influenza A",
             action="Retained as Influenza-genus near-neighbor; flagged auxiliary - must not count toward IAV within-pathogen replication at Stage 6"),
    ]
    with open(f"{OUTDIR}/stage1_manifest_corrections.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["item", "severity", "finding", "action"], delimiter="\t")
        w.writeheader()
        w.writerows(corr)

    # summary
    from collections import Counter
    st = Counter(r["adjudication_status"] for r in rows_out)
    fam = Counter(r["endpoint_family"] for r in rows_out if str(r["adjudication_status"]).startswith("ADJUDICATED"))
    print("adjudication_status:", dict(st))
    print("endpoint_family (adjudicated):", dict(fam))
    print("rows:", len(rows_out))
    for fn in ["hpicda_ontology_v1.tsv", "screen_endpoint_direction_adjudication.tsv", "stage1_manifest_corrections.tsv"]:
        h = hashlib.sha256(open(f"{OUTDIR}/{fn}", "rb").read()).hexdigest()
        print(f"{fn}  sha256={h}")

if __name__ == "__main__":
    main()
