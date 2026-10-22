#!/usr/bin/env python3
"""Apply the 2026-10-21 HPICDA revision to the copied v2 manuscript."""

from __future__ import annotations

import os
import re
from pathlib import Path


ROOT = Path(
    os.environ.get("HPICDA_ROOT", Path(__file__).resolve().parents[1])
).resolve()
PACKAGE = ROOT / "25_manuscript_package"
V2_DIR = PACKAGE / "revision_v2_2026-10-21"
MANUSCRIPT = V2_DIR / "HPICDA_manuscript_v2.md"


TITLE = (
    "# A pre-registered benchmark identifies transferable and non-transferable "
    "signals in a pulmonary infection CRISPR atlas"
)
RUNNING_TITLE = "**Running title:** Pre-registered CRISPR atlas benchmark"

ABSTRACT = """CRISPR host–pathogen screens have proliferated, but candidate rankings are rarely benchmarked under pre-fixed rules. We present a pre-registered benchmark using 39 CRISPR infection screens covering SARS-CoV-2, influenza A, RSV, *M. tuberculosis*/BCG and *P. aeruginosa*. We assembled the Human Pulmonary Infection Conditional Dependency Atlas (HPICDA), comprising 852,216 gene × screen records reconciled under an estimand-aware scheme that separates conditional dependency, defined against each gene's own baseline cell fitness, from baseline-fitness confounding and signal-state alteration, and froze 353 candidates before adjudication against unused screens. Same-pathogen replication was 55 of 654 evaluable gene × holdout pairs (8.4%). The HG3 composite was at chance on the pre-registered gene-level outcome (AUROC 0.524, 95% CI 0.434–0.617). A top-10 pair-level signal was present (36.0% versus 7.3%), but equal-opportunity resampling gave a wide uncertainty interval. Pathogen count initially reached AUROC 0.686, but breadth was strongly associated with the number of holdout opportunities. In 5,000 equal-opportunity resamples, pathogen count did not robustly exceed HG3 (difference +0.072, 95% interval −0.055 to +0.194), whereas independent discovery units hit did (+0.108, +0.004 to +0.207). Cellular context remained inseparable from holdout identity and is exploratory. Three archetypes span recurrent-hit meanings: a robust anchor (ACE2), a context-dependent module (heparan-sulfate biosynthesis) and a mechanistically real but largely baseline-essential module (V-ATPase). We release the call table, frozen set and adjudication as a benchmark and Pareto-front resource for candidate prioritization."""

IMPORTANCE = """Genome-wide CRISPR screens are now routine in respiratory infection research, but weighted rankings of host dependencies are rarely benchmarked against independent studies. We froze 353 candidates from 39 pulmonary infection screens and evaluated their transferability under pre-registered rules. The same-pathogen replication rate was 8.4%, and the frozen composite did not discriminate global validation. A simple pathogen count appeared superior, but this advantage largely reflected how many validation opportunities each gene received. After giving every gene one holdout opportunity, a count of independent discovery units, rather than pathogen breadth, retained predictive value. The benchmark supports reporting candidate evidence as a membership set and identifies the features that transfer."""

INTRODUCTION_RESULTS = """Three results define the benchmark's output. First, across 654 pre-registered gene × holdout pairs the same-pathogen PASS rate was 8.4%, quantifying replication at the individual-call level. Second, the HG3 composite was at chance as a global predictor of independent validation (AUROC 0.524), while its top ranks were enriched. Pathogen breadth did not account for transfer after equalizing validation opportunities; the count or frequency of independent discovery units hit did. Third, cellular context may contribute to heterogeneity, but it cannot be separated from holdout identity in this design and is therefore reported as exploratory."""

OPPORTUNITY_SECTION = """### Opportunity-balanced benchmarking identifies discovery-unit breadth

The original gene-level endpoint counted a candidate as validated if it passed in at least one same-pathogen holdout. That endpoint is sensitive to the number of evaluable pairs: pathogen breadth was strongly associated with holdout opportunity (Spearman rho = 0.647, *P* = 7.7 × 10⁻⁴¹). We therefore added a revision-stage sensitivity analysis that gives every gene the same number of opportunities. In each of 5,000 replicates, one same-pathogen pair was selected at random for each of the 333 evaluable genes, and the AUROCs of the frozen HG3 score and simple discovery features were recomputed (**Supplementary Tables S21a and S21b**).

Under equal opportunity, the HG3 score remained at chance (median AUROC 0.538, 95% resampling interval 0.441–0.640). Pathogen breadth remained higher than HG3 in the median but did not retain a robust paired advantage (median difference +0.072, 95% interval −0.055 to +0.194). In contrast, the number of independent discovery units hit exceeded HG3 (median AUROC 0.646; difference +0.108, 95% interval +0.004 to +0.207), as did hit frequency (difference +0.108, 95% interval +0.004 to +0.206). The same distinction held in logistic models adjusted for the number of evaluable pairs: pathogen count was not retained (OR 2.38, 95% CI 0.77–7.41, *P* = 0.133), whereas independent-unit count was retained (OR 1.46, 95% CI 1.17–1.81, *P* = 6.8 × 10⁻⁴). The robust consensus feature was therefore the count or frequency of supporting independent discovery units, rather than pathogen breadth.

The top of the HG3 ranking also remained enriched under equal-opportunity resampling, but with wide uncertainty: the top 10 had a median PASS rate of 0.30 (95% interval 0.10–0.60) and a median odds ratio of 7.26 (95% interval 1.68–22.07). This supports a top-heavy prioritization signal without supporting a precise global ordering.

### Cellular context is exploratory and not identifiable from holdout identity

Context-matched holdouts had a higher PASS rate than mismatched holdouts (25.9% versus 5.8%), and a pair-level logistic model with holdout fixed effects estimated an odds ratio of 2.58 (95% CI 1.11–6.02). This estimate cannot be interpreted as a determinant of replication in the present design. The context flag is constant within four of the five frozen holdouts, so the marginal comparison is also a comparison between holdouts, and only HO-V-SC2-01 identifies the effect within a holdout. The largest matched, preserve/activate stratum is driven mainly by the supplementary SC2-03 arm, and the context-by-direction interaction falls from OR 8.54 in the pooled data to 1.12 (95% CI 0.26–4.84) without that arm (**Supplementary Tables S19d and S20b**). We therefore report context as an exploratory source of heterogeneity and not as a validated determinant of replication."""

NEW_BENCHMARK_PARAGRAPH = """On the SC2-03-inclusive outcome, global discrimination was reversed. The HG3 total score achieved an AUROC of 0.447 with an AUPRC of 0.419 against a prevalence of 0.396. The number of pathogens hit achieved 0.574 (AUPRC 0.531) and the number of units hit 0.562. On the frozen outcome the composite was at chance: AUROC 0.524 (95% CI 0.428–0.617, which includes 0.5). The below-chance figure is specific to the supplementary-inclusive outcome and is produced by the 81 candidates that validate only in the permissive SC2-03 arm. A pathogen count reached AUROC 0.686 on the frozen any-PASS outcome, but this endpoint gives genes with more evaluable pairs more chances to pass; opportunity-balanced analysis showed that the pathogen-breadth advantage was not robust. The number of evaluable pairs itself reached AUROC 0.748 and is treated as a nuisance variable."""

NEW_CONTEXT_STRATIFICATION_PARAGRAPH = """First, the context flag is constant within four of the five frozen holdouts—every evaluable pair is context-matched in HO-V-IAV-01 and none is in HO-V-IAV-02, HO-V-RSV-01 or HO-V-SC2-02—so the marginal contrast is in large part a contrast between holdouts. Only HO-V-SC2-01 identifies the effect within a holdout, and there context-matched pairs passed at 20.0% (12/60) against 8.0% (16/201). Second, the stratum that combines context match with a preserve/activate direction and gives the largest contrast, 62.3%, is a pooled figure: 77.2% of its 114 pairs come from the supplementary SC2-03 arm. In the frozen adjudication alone the same stratum is 5 of 26 pairs, or 19.2%, and the context-by-direction interaction the pooled data estimate at odds ratio 8.54 (95% CI 3.41–21.40) falls to 1.12 (95% CI 0.26–4.84, not significant) once SC2-03 is removed (**Supplementary Table S19d**, **S20**). This stratification is exploratory and is not evidence that context and direction act jointly in the pre-registered layer."""

NEW_ADJUDICATION_HEADING = "### Independent adjudication identifies a concentrated top-rank signal"

NEW_TOP_RANK_PARAGRAPH = """Ranking quality was assessed by top-*K* enrichment at the pair level (**Fig. 7b**). The top 10 candidates by frozen total score contributed 25 evaluable pairs, of which 9 passed (36.0%), against 46 of 629 (7.3%) in the remainder: an odds ratio of 7.13 (*P* = 7.8 × 10⁻⁵), with gene recall of 70% among evaluable top-10 candidates. Enrichment decayed rapidly with *K*: odds ratios of 3.91 at *K* = 20, 1.82 at *K* = 50 and 1.53 at *K* = 100, the last two not significant. Under the equal-opportunity analysis below, the top-10 signal had a median odds ratio of 7.26 but a 95% resampling interval of 1.68–22.07, with a wide interval for the P value as well. The result is therefore a limited top-rank signal, not evidence for a stable set of individually robust candidates."""

NEW_QUALIFICATION_PARAGRAPH = """Two qualifications define the limits of global discrimination. First, allowing the composite sign to reverse raises the SC2-03-inclusive AUROC from 0.447 to 0.553, close to the simple unit-count comparator, while on the frozen layer it changes 0.524 to 0.476. Its global discrimination is therefore direction-dependent rather than absent. Second, the composite retains a high hit rate within its top ranks (8 of the top 10 under the reported tie-break; 0.80–1.00 across admissible tie-breaks). This combination supports a top-head prioritization tool rather than a transferable scalar ranking. Weight sensitivity over 2,005 schemes left only ACE2 at rank 1, so candidate identities are reported as a Pareto-front membership set (**Supplementary Table S17b**)."""

NEW_DISCUSSION_SECOND = """The second result concerns the ranking machinery. The HG3 composite is at chance as a global discriminator of independent validation on the pre-registered layer (AUROC 0.524) and below chance on the supplementary-inclusive layer (0.447). The head of the ranking is enriched, but this is a local property: under equal-opportunity resampling, the top-10 PASS rate had a wide interval (0.10–0.60), and the global ordering did not transfer. We interpret HG3 as a top-head prioritization tool rather than a transferable ranking. The construction itself still has value: the estimand registry makes explicit what would otherwise be silently averaged, the baseline-fitness gate removes a class of artefact otherwise indistinguishable from signal, the prospective freeze makes the adjudication meaningful, and the biological-before-drug firewall keeps pharmacological priors out of discovery."""

NEW_DISCUSSION_THIRD = """The third result is a boundary on simple alternatives and a limitation on context. An unweighted pathogen count appeared superior on the original any-PASS endpoint (AUROC 0.686), but pathogen breadth was strongly associated with the number of holdout opportunities. After equal-opportunity resampling, that advantage was not robust (difference +0.072, 95% interval −0.055 to +0.194). The simpler feature that retained value was the count or frequency of independent discovery units hit (difference +0.108, 95% interval +0.004 to +0.207). The evidence therefore supports a simple consensus count, but not pathogen breadth specifically. Cellular context may matter, but this benchmark cannot separate it from holdout identity: the context flag is constant within four of five frozen holdouts, one holdout identifies the within-holdout contrast, and the largest interaction cell is driven by the supplementary arm. Context is reported as exploratory heterogeneity, not as the determinant of replication."""

NEW_PRACTICE_SECOND = """Second, **treat context as a design variable.** Context-matched strata differed descriptively, but the current benchmark cannot isolate context from holdout identity. Future validation should vary cell context within each holdout rather than infer context effects from comparisons between holdouts."""

NEW_PRACTICE_THIRD = """Third, **design the benchmark so that predictors and outcomes are separable.** The original any-PASS endpoint gave genes with more pairs more chances to succeed. Equal-opportunity resampling is a minimum sensitivity check; a purpose-built benchmark should also hold the number and type of opportunities approximately fixed and commit the outcome to one layer before scoring."""

NEW_PRACTICE_FOURTH = """Fourth, **benchmark the score against simple alternatives, and report every layer separately.** Equal-opportunity analysis identifies independent discovery-unit breadth as the feature that transfers, whereas pathogen breadth does not. This distinction is the practical value of a benchmark that compares a score with plausible baselines."""

NEW_LIMITATION_13 = """13. **The benchmark cannot separate cellular context from holdout identity.** The context flag is constant within four of the five frozen holdouts, so the marginal context estimate is confounded with holdout identity; only HO-V-SC2-01 identifies the within-holdout effect. The context-by-direction interaction is estimated almost entirely by the supplementary arm and does not survive its removal. Context effects are therefore exploratory. The original gene-level any-PASS endpoint was also sensitive to the number of evaluable pairs; this was addressed by adjustment and 5,000 equal-opportunity resamples, but the outcome remains discrete and does not itself model all dependence among repeated pairs."""

NEW_LIMITATION_6 = """5. **The composite score is weight-sensitive.** Only ACE2 is stable across the 2,005 weightings tested, and the top-10 pair-level signal is imprecise under opportunity-balanced resampling. Downstream communication therefore uses Pareto-front membership, with benchmark performance reported separately for the frozen and SC2-03-inclusive layers."""

NEW_METHOD_PARAGRAPH = """Because the original gene-level endpoint (validation in at least one same-pathogen holdout) is sensitive to the number of evaluable pairs, opportunity dependence was quantified and addressed in two ways. First, logistic models of the any-PASS outcome included log(1 + number of same-pathogen evaluable pairs); the pathogen-count and unit-count models were fitted separately because those features are alternative ranking rules. Second, an equal-opportunity resampling analysis selected exactly one same-pathogen pair at random for every gene in each of 5,000 replicates (seed 20261021). AUROCs of the frozen HG3 score and the simple discovery features were recomputed in each replicate; paired differences were taken within the same sampled outcome. The median and 2.5th–97.5th percentiles of the resampling distribution are reported. This analysis changes neither the frozen candidates nor the adjudicated holdout verdicts (**Supplementary Tables S21a and S21b**)."""

NEW_CONTEXT_METHOD = """The comparison between cellular context, intervention direction and the composite score is reported in **Supplementary Table S20**. Context exposure is defined per gene as the fraction of its evaluable same-pathogen pairs whose holdout matched the discovery cellular context, computed once over all pairs and once over the frozen pairs only. Pair-level effects were estimated by logistic regression of the PASS indicator on the context flag, intervention direction and Bio-Tier 2 status, marginally and with holdout fixed effects; the fixed-effect fit excludes HO-V-IAV-02, where separation occurs. Because context is constant within four of five frozen holdouts, the marginal and fixed-effect context estimates are exploratory rather than design-identified. Direction was coded preserve/activate versus inhibit. Both layers are reported separately."""

NEW_FIG7 = """**Figure 7 | Holdout adjudication and opportunity-balanced benchmarking.** (**a**) Pre-registered holdout PASS rate by stratum; the dashed line marks the overall rate. (**b**) Top-*K* enrichment of holdout PASS at the pair level (top 10: 9 of 25 pairs, 36%). The top-rank signal is limited and imprecise after opportunity-balanced resampling. (**c**) AUROC of the HG3 composite score against naive consensus features, at the gene level over 333 evaluable candidates, computed separately on the frozen pre-registered outcome and on the SC2-03-inclusive outcome; bars are the frozen pre-registered layer and the hollow marker is the SC2-03-inclusive layer (**Supplementary Table S19f**). (**d**) Equal-opportunity sensitivity for the frozen outcome: one same-pathogen pair was sampled per gene in each of 5,000 replicates. Points and intervals show the median paired AUROC difference from HG3. The apparent pathogen-breadth advantage was not robust (difference +0.072, 95% interval −0.055 to +0.194), whereas independent-unit breadth and hit frequency retained a positive advantage (**Supplementary Tables S21a and S21b**)."""

NEW_DRUG_RESULT = """Candidate–compound links are released in Supplementary Table S15 as database-derived provenance and were not used for candidate selection or quantitative ranking. On the tolerability axis (**Supplementary Table S16**), 237 of 332 evaluable candidates carry a human-intolerance flag under a deliberately permissive definition; this axis is used for triage only and never to exclude a candidate."""

NEW_DEPMAP_METHOD = """Baseline fitness was taken from the DepMap Public 26Q1 Chronos gene-effect matrix (Figshare DOI 10.6084/m9.figshare.31660582.v2) computed across the pulmonary cell-line panel (A-549, Calu-3, HEp-2, Huh-7.5, Huh-7.5.1, and the wider reference panel), Supplementary Table S4. Because a complete 26Q1 Figshare bundle was unavailable when the baseline data were frozen, model annotations and the common-essential list came from DepMap 24Q4 Public (DOI 10.25452/figshare.plus.27993248.v1); the gene-effect matrix is 26Q1 Chronos, and schema compatibility and file hashes are recorded in Supplementary Table S4. Chronos is a cell-population-dynamics model of CRISPR fitness effects [PMID 34930405]; DepMap dependency resources are described in [PMID 28753430, PMID 29083409]."""

NEW_TRACTABILITY_METHOD = """Author-reported launched and clinical-stage compounds were unlocked only after the freeze. Tractability was obtained directly from the Open Targets Platform 26.06 GraphQL API, the release current on 2026-09-03 [PMID 39657122]. Launched and clinical-stage compounds, phase and mechanism annotations were obtained from the Broad Drug Repurposing Hub snapshot used on 2026-09-03 [PMID 28388612]. ChEMBL records [PMID 37933841] and Orphanet annotations [PMID 31527858] reached this analysis only through Open Targets and were not separately queried. Genetics annotation used gnomAD LOEUF [PMID 32461654], Genebass [PMID 36778668] and the Human Protein Atlas. Database-derived candidate–compound links are reported in Supplementary Table S15 for provenance and were not used for candidate selection or quantitative ranking."""

NEW_FIG6 = """**Figure 6 | Baseline-fitness liability and tractability.** (**a**) Baseline cellular-fitness liability landscape for the 343 frozen candidates with a DepMap minimum baseline Chronos score: minimum baseline Chronos against HG3 total score, coloured by biological tier, with selected genes labelled. The dashed line marks the essentiality threshold (−0.5); candidates to its right are non-essential in the cell-line panel and therefore lower-liability, which is not the same as a demonstrated therapeutic window. (**b**) Open Targets tractability class by intervention direction."""

NEW_GOVERNANCE_METHOD = """The governance sequence comprised algorithmically documented release gates HG0–HG4. Each gate produced an owner-authorized approval artefact with input and rule hashes; approvals were rollback-capable, and all superseded hashes were retained. The candidate freeze carries SHA-256 `1c00f7cc694b3dff15dd671f1288e347726e0b69f0af1657892e71bd02d8d4dc`; a full reproducibility manifest is in Supplementary Table S18."""

REFERENCE_ADDITIONS = """24. DepMap Broad Institute. Chronos parameters (Public 26Q1), version 2. *figshare* 2026. DOI: 10.6084/m9.figshare.31660582.v2
25. DepMap Broad Institute. DepMap 24Q4 Public. *figshare+* 2024. DOI: 10.25452/figshare.plus.27993248.v1
26. Buniello A, Suveges D, Cruz-Castillo C, et al. Open Targets Platform: facilitating therapeutic hypotheses building in drug discovery. *Nucleic Acids Res* 2025. PMID: 39657122. DOI: 10.1093/nar/gkae1128
27. Zdrazil B, Felix E, Hunter F, et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. *Nucleic Acids Res* 2024. PMID: 37933841. DOI: 10.1093/nar/gkad1004
28. Corsello SM, Bittker JA, Liu Z, et al. The Drug Repurposing Hub: a next-generation drug library and information resource. *Nat Med* 2017. PMID: 28388612. DOI: 10.1038/nm.4306
29. Nguengang Wakap S, Lambert DM, Olry A, et al. Estimating cumulative point prevalence of rare diseases: analysis of the Orphanet database. *Eur J Hum Genet* 2020. PMID: 31527858. DOI: 10.1038/s41431-019-0508-0"""


def replace_paragraph(text: str, prefix: str, replacement: str) -> str:
    parts = text.split("\n\n")
    hits = [i for i, part in enumerate(parts) if part.strip().startswith(prefix)]
    if len(hits) != 1:
        raise RuntimeError(f"Expected one paragraph beginning {prefix!r}; found {len(hits)}")
    parts[hits[0]] = replacement
    return "\n\n".join(parts)


def replace_section(text: str, heading: str, replacement: str) -> str:
    parts = text.split("\n\n")
    hits = [i for i, part in enumerate(parts) if part.strip() == heading]
    if len(hits) != 1:
        raise RuntimeError(f"Expected one section heading {heading!r}; found {len(hits)}")
    start = hits[0]
    end = len(parts)
    for i in range(start + 1, len(parts)):
        if parts[i].startswith("### ") or parts[i].startswith("## "):
            end = i
            break
    parts[start:end] = [replacement]
    return "\n\n".join(parts)


def replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Expected one occurrence of {old[:80]!r}; found {count}")
    return text.replace(old, new, 1)


def replace_line(text: str, prefix: str, replacement: str) -> str:
    lines = text.split("\n")
    hits = [i for i, line in enumerate(lines) if line.startswith(prefix)]
    if len(hits) != 1:
        raise RuntimeError(f"Expected one line beginning {prefix!r}; found {len(hits)}")
    lines[hits[0]] = replacement
    return "\n".join(lines)


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")

    text = replace_once(
        text,
        "# Independent validation exposes context-dependent reproducibility and the limits of composite scoring across pulmonary infection CRISPR screens",
        TITLE,
    )
    text = replace_once(
        text,
        "**Running title:** Context-dependent reproducibility and composite-score limits",
        RUNNING_TITLE,
    )
    text = replace_paragraph(
        text,
        "Genome-wide perturbation screens of host–pathogen interactions have proliferated",
        ABSTRACT,
    )
    text = replace_paragraph(
        text,
        "Genome-wide CRISPR screens are now run routinely",
        IMPORTANCE,
    )
    text = replace_paragraph(
        text,
        "Three results answer those questions, and two of them are negative.",
        INTRODUCTION_RESULTS,
    )

    text = replace_paragraph(
        text,
        "On the SC2-03-inclusive outcome the result was negative.",
        NEW_BENCHMARK_PARAGRAPH,
    )
    text = replace_paragraph(
        text,
        "First, the context flag is constant within four of the five frozen holdouts",
        NEW_CONTEXT_STRATIFICATION_PARAGRAPH,
    )
    text = replace_once(
        text,
        "### Independent adjudication: weak overall, informative at the head",
        NEW_ADJUDICATION_HEADING,
    )
    text = replace_paragraph(
        text,
        "Ranking quality was assessed by top-*K* enrichment at the pair level",
        NEW_TOP_RANK_PARAGRAPH,
    )
    text = replace_paragraph(
        text,
        "Two further readings qualify the negative result",
        NEW_QUALIFICATION_PARAGRAPH,
    )
    text = replace_section(
        text,
        "### What predicts replication: context, not the composite",
        OPPORTUNITY_SECTION,
    )

    text = replace_paragraph(
        text,
        "The second is that the ranking machinery does not do what it was built to do.",
        NEW_DISCUSSION_SECOND,
    )
    text = replace_paragraph(
        text,
        "The third finding concerns what does predict replication.",
        NEW_DISCUSSION_THIRD,
    )
    text = replace_paragraph(
        text,
        "Second, **treat context as a variable, not as noise.**",
        NEW_PRACTICE_SECOND,
    )
    text = replace_paragraph(
        text,
        "Third, **design the benchmark so that its predictors and its outcome are separable.**",
        NEW_PRACTICE_THIRD,
    )
    text = replace_paragraph(
        text,
        "Fourth, **benchmark the score against trivial alternatives, and report every layer separately.**",
        NEW_PRACTICE_FOURTH,
    )
    text = replace_line(
        text,
        "13. **The benchmark cannot separate cellular context from holdout identity.**",
        NEW_LIMITATION_13,
    )
    text = replace_line(
        text,
        "6. **The composite score is weight-sensitive and does not validate.**",
        NEW_LIMITATION_6,
    )

    text = replace_paragraph(
        text,
        "The HG3 total score was compared with four naive consensus features",
        NEW_METHOD_PARAGRAPH,
    )
    text = replace_paragraph(
        text,
        "The comparison between cellular context, intervention direction and the composite score",
        NEW_CONTEXT_METHOD,
    )
    text = replace_paragraph(
        text,
        "We then audited the drug annotations themselves",
        NEW_DRUG_RESULT,
    )
    text = replace_paragraph(
        text,
        "Baseline fitness was taken from DepMap 26Q1",
        NEW_DEPMAP_METHOD,
    )
    text = replace_paragraph(
        text,
        "Author-reported launched and clinical-stage compounds",
        NEW_TRACTABILITY_METHOD,
    )
    text = replace_paragraph(
        text,
        "**Figure 6 | Baseline-fitness liability and tractability.**",
        NEW_FIG6,
    )
    text = replace_paragraph(
        text,
        "The governance sequence comprised gates HG0–HG4",
        NEW_GOVERNANCE_METHOD,
    )
    reference_23 = (
        "23. Gillespie M, Jassal B, Stephan R, et al. The Reactome pathway "
        "knowledgebase 2022. *Nucleic Acids Res* 2022. PMID: 34788843"
    )
    text = replace_once(text, reference_23, reference_23 + "\n" + REFERENCE_ADDITIONS)
    text = replace_once(
        text,
        "all supplementary tables S1–S20b",
        "all supplementary tables S1–S21b",
    )
    text = replace_paragraph(
        text,
        "**Figure 7 | Holdout adjudication and honest benchmarking.**",
        NEW_FIG7,
    )

    MANUSCRIPT.write_text(text, encoding="utf-8", newline="\n")
    abstract_words = len(re.sub(r"[*_`]", "", ABSTRACT).split())
    importance_words = len(re.sub(r"[*_`]", "", IMPORTANCE).split())
    print(f"Wrote {MANUSCRIPT}")
    print(f"Abstract words: {abstract_words}")
    print(f"Importance words: {importance_words}")
    if abstract_words > 250 or importance_words > 150:
        raise RuntimeError("Word limit exceeded")


if __name__ == "__main__":
    main()
