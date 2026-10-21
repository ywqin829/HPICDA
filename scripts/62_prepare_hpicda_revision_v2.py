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
    "# Frozen validation of a pulmonary infection CRISPR atlas shows limited "
    "cross-screen reproducibility and non-transferable composite weighting"
)
RUNNING_TITLE = "**Running title:** Limited CRISPR atlas reproducibility"

ABSTRACT = """Genome-wide perturbation screens of host–pathogen interactions have proliferated, but their reproducibility is rarely tested under pre-fixed rules. We assembled 39 analysable CRISPR infection screens covering SARS-CoV-2, influenza A, RSV, *M. tuberculosis*/BCG and *P. aeruginosa* into the Human Pulmonary Infection Conditional Dependency Atlas (HPICDA): 852,216 gene × screen records reconciled under an estimand-aware scheme that separates conditional dependency—defined against each gene's own baseline cell fitness—from baseline-fitness confounding and signal-state alteration, yielding a frozen 353-candidate set adjudicated against screens not used for discovery. Same-pathogen validation was weak: 55 of 654 evaluable gene × holdout pairs passed (8.4%). The frozen HG3 composite was at chance on the pre-registered gene-level outcome (AUROC 0.524, 95% CI 0.434–0.617). A top-10 pair-level signal was present (36.0% versus 7.3%), but equal-opportunity resampling gave a wide uncertainty interval. An unweighted pathogen count initially reached AUROC 0.686, but pathogen breadth was strongly associated with the number of holdout opportunities. In 5,000 equal-opportunity resamples selecting one pair per gene, pathogen count did not robustly exceed HG3 (difference +0.072, 95% interval −0.055 to +0.194), whereas the number of independent discovery units hit did (+0.108, +0.004 to +0.207). Cellular context could not be separated from holdout identity and is therefore exploratory. Three archetypes span what a recurrent hit can represent: a robust anchor (ACE2), a context-dependent module (heparan-sulfate biosynthesis) and a mechanistically real but largely baseline-essential module (V-ATPase). We release the call table, frozen set and adjudication as a Pareto-front resource rather than a precise target ranking."""

IMPORTANCE = """Genome-wide CRISPR screens are now routine in respiratory infection research, but weighted rankings of host dependencies are rarely tested against independent studies. We froze 353 candidates from 39 pulmonary infection screens before validation. The same-pathogen replication rate was 8.4%, and the frozen composite score was globally uninformative. A simple pathogen count appeared superior, but this advantage largely reflected how many validation opportunities each gene received. After giving every gene one holdout opportunity, a count of independent discovery units, rather than pathogen breadth, retained predictive value. The results support reporting candidate evidence as a membership set rather than as a precisely transferable ranking."""

INTRODUCTION_RESULTS = """Three results answer those questions, and two of them are negative. First, reproducibility is low even under favourable conditions: across 654 pre-registered gene × holdout pairs the same-pathogen PASS rate is 8.4%, so most individual calls do not reproduce. Second, the HG3 composite that ranks the candidate set is at chance as a global predictor of independent validation (AUROC 0.524), although its top ranks are enriched. The apparent superiority of a simple pathogen count is not retained after equalizing the number of validation opportunities; the simple feature that survives is the count or frequency of independent discovery units hit. Third, cellular context may contribute to heterogeneity, but it cannot be separated from holdout identity in this design and is therefore reported as exploratory rather than as a determinant of replication."""

OPPORTUNITY_SECTION = """### Opportunity-balanced benchmarking removes the apparent pathogen-breadth advantage

The original gene-level endpoint counted a candidate as validated if it passed in at least one same-pathogen holdout. That endpoint is sensitive to the number of evaluable pairs: pathogen breadth was strongly associated with holdout opportunity (Spearman rho = 0.647, *P* = 7.7 × 10⁻⁴¹). We therefore added a revision-stage sensitivity analysis that gives every gene the same number of opportunities. In each of 5,000 replicates, one same-pathogen pair was selected at random for each of the 333 evaluable genes, and the AUROCs of the frozen HG3 score and simple discovery features were recomputed (**Supplementary Tables S21a and S21b**).

Under equal opportunity, the HG3 score remained at chance (median AUROC 0.538, 95% resampling interval 0.441–0.640). Pathogen breadth remained higher than HG3 in the median but did not retain a robust paired advantage (median difference +0.072, 95% interval −0.055 to +0.194). In contrast, the number of independent discovery units hit exceeded HG3 (median AUROC 0.646; difference +0.108, 95% interval +0.004 to +0.207), as did hit frequency (difference +0.108, 95% interval +0.004 to +0.206). The same distinction held in logistic models of the original any-PASS outcome adjusted for the number of evaluable pairs: pathogen count was not retained (OR 2.38, 95% CI 0.77–7.41, *P* = 0.133), whereas the number of independent units hit was retained (OR 1.46, 95% CI 1.17–1.81, *P* = 6.8 × 10⁻⁴). We therefore withdraw the pathogen-breadth claim. The evidence supports a simpler count or frequency of supporting independent discovery units, not a claim that pathogen breadth predicts transfer.

The top of the HG3 ranking also remained enriched under equal-opportunity resampling, but with wide uncertainty: the top 10 had a median PASS rate of 0.30 (95% interval 0.10–0.60) and a median odds ratio of 7.26 (95% interval 1.68–22.07). This supports a top-heavy prioritization signal without supporting a precise global ordering.

### Cellular context is exploratory and not identifiable from holdout identity

Context-matched holdouts had a higher PASS rate than mismatched holdouts (25.9% versus 5.8%), and a pair-level logistic model with holdout fixed effects estimated an odds ratio of 2.58 (95% CI 1.11–6.02). This estimate cannot be interpreted as a determinant of replication in the present design. The context flag is constant within four of the five frozen holdouts, so the marginal comparison is also a comparison between holdouts, and only HO-V-SC2-01 identifies the effect within a holdout. The largest matched, preserve/activate stratum is driven mainly by the supplementary SC2-03 arm, and the context-by-direction interaction falls from OR 8.54 in the pooled data to 1.12 (95% CI 0.26–4.84) without that arm (**Supplementary Tables S19d and S20b**). We therefore report context as an exploratory source of heterogeneity and not as a validated determinant of replication."""

NEW_BENCHMARK_PARAGRAPH = """On the SC2-03-inclusive outcome the result was negative. The HG3 total score achieved an AUROC of 0.447—*below* chance—with an AUPRC of 0.419 against a prevalence of 0.396. The number of pathogens hit achieved 0.574 (AUPRC 0.531) and the number of units hit 0.562. On the frozen outcome the composite is at chance rather than reversed: AUROC 0.524 (95% CI 0.428–0.617, which includes 0.5). The below-chance figure is specific to the supplementary-inclusive outcome and is produced by the 81 candidates that validate only in the permissive SC2-03 arm. A pathogen count reached AUROC 0.686 on the frozen any-PASS outcome, but this endpoint gives genes with more evaluable pairs more chances to pass; the opportunity-balanced analysis below shows that the pathogen-breadth advantage is not robust. The number of evaluable pairs itself reached AUROC 0.748 and is treated as a nuisance variable."""

NEW_CONTEXT_STRATIFICATION_PARAGRAPH = """First, the context flag is constant within four of the five frozen holdouts—every evaluable pair is context-matched in HO-V-IAV-01 and none is in HO-V-IAV-02, HO-V-RSV-01 or HO-V-SC2-02—so the marginal contrast is in large part a contrast between holdouts. Only HO-V-SC2-01 identifies the effect within a holdout, and there context-matched pairs passed at 20.0% (12/60) against 8.0% (16/201). Second, the stratum that combines context match with a preserve/activate direction and gives the largest contrast, 62.3%, is a pooled figure: 77.2% of its 114 pairs come from the supplementary SC2-03 arm. In the frozen adjudication alone the same stratum is 5 of 26 pairs, or 19.2%, and the context-by-direction interaction the pooled data estimate at odds ratio 8.54 (95% CI 3.41–21.40) falls to 1.12 (95% CI 0.26–4.84, not significant) once SC2-03 is removed (**Supplementary Table S19d**, **S20**). This stratification is exploratory and is not evidence that context and direction act jointly in the pre-registered layer."""

NEW_ADJUDICATION_HEADING = "### Independent adjudication: low overall and a limited top-rank signal"

NEW_TOP_RANK_PARAGRAPH = """Ranking quality was assessed by top-*K* enrichment at the pair level (**Fig. 7b**). The top 10 candidates by frozen total score contributed 25 evaluable pairs, of which 9 passed (36.0%), against 46 of 629 (7.3%) in the remainder: an odds ratio of 7.13 (*P* = 7.8 × 10⁻⁵), with gene recall of 70% among evaluable top-10 candidates. Enrichment decayed rapidly with *K*: odds ratios of 3.91 at *K* = 20, 1.82 at *K* = 50 and 1.53 at *K* = 100, the last two not significant. Under the equal-opportunity analysis below, the top-10 signal had a median odds ratio of 7.26 but a 95% resampling interval of 1.68–22.07, with a wide interval for the P value as well. The result is therefore a limited top-rank signal, not evidence for a stable set of individually robust candidates."""

NEW_QUALIFICATION_PARAGRAPH = """Two qualifications constrain the negative result without rescuing the weighting scheme. First, the composite is not noise: allowing its sign to reverse raises the SC2-03-inclusive AUROC from 0.447 to 0.553, close to the simple unit-count comparator, while on the frozen layer it changes 0.524 to 0.476. The failure is therefore in the stated direction and weights rather than in the complete absence of rank information. Second, the composite retains a high hit rate within its top ranks (8 of the top 10 under the reported tie-break; 0.80–1.00 across admissible tie-breaks). The dissociation—sharp at the extreme and uninformative globally—supports a top-head prioritization tool, not a transferable ranking. Weight sensitivity over 2,005 schemes left only ACE2 at rank 1, so candidate identities are reported as a Pareto-front membership set rather than a scalar ordering (**Supplementary Table S17b**)."""

NEW_DISCUSSION_SECOND = """The second result concerns the ranking machinery. The HG3 composite is at chance as a global discriminator of independent validation on the pre-registered layer (AUROC 0.524) and below chance on the supplementary-inclusive layer (0.447). The head of the ranking is enriched, but this is a local property: under equal-opportunity resampling, the top-10 PASS rate had a wide interval (0.10–0.60), and the global ordering did not transfer. We interpret HG3 as a top-head prioritization tool rather than a transferable ranking. The construction itself still has value: the estimand registry makes explicit what would otherwise be silently averaged, the baseline-fitness gate removes a class of artefact otherwise indistinguishable from signal, the prospective freeze makes the adjudication meaningful, and the biological-before-drug firewall keeps pharmacological priors out of discovery."""

NEW_DISCUSSION_THIRD = """The third result is a boundary on simple alternatives and a limitation on context. An unweighted pathogen count appeared superior on the original any-PASS endpoint (AUROC 0.686), but pathogen breadth was strongly associated with the number of holdout opportunities. After equal-opportunity resampling, that advantage was not robust (difference +0.072, 95% interval −0.055 to +0.194). The simpler feature that retained value was the count or frequency of independent discovery units hit (difference +0.108, 95% interval +0.004 to +0.207). The evidence therefore supports a simple consensus count, but not pathogen breadth specifically. Cellular context may matter, but this benchmark cannot separate it from holdout identity: the context flag is constant within four of five frozen holdouts, one holdout identifies the within-holdout contrast, and the largest interaction cell is driven by the supplementary arm. Context is reported as exploratory heterogeneity, not as the determinant of replication."""

NEW_PRACTICE_SECOND = """Second, **treat context as a design variable.** Context-matched strata differed descriptively, but the current benchmark cannot isolate context from holdout identity. Future validation should vary cell context within each holdout rather than infer context effects from comparisons between holdouts."""

NEW_PRACTICE_THIRD = """Third, **design the benchmark so that predictors and outcomes are separable.** The original any-PASS endpoint gave genes with more pairs more chances to succeed. Equal-opportunity resampling is a minimum sensitivity check; a purpose-built benchmark should also hold the number and type of opportunities approximately fixed and commit the outcome to one layer before scoring."""

NEW_PRACTICE_FOURTH = """Fourth, **benchmark the score against trivial alternatives, and report every layer separately.** The equal-opportunity analysis withdraws the apparent pathogen-breadth advantage but supports a simple count or frequency of independent discovery units. That distinction is the difference between a score that has been compared with a plausible baseline and one that has merely been used."""

NEW_LIMITATION_13 = """13. **The benchmark cannot separate cellular context from holdout identity.** The context flag is constant within four of the five frozen holdouts, so the marginal context estimate is confounded with holdout identity; only HO-V-SC2-01 identifies the within-holdout effect. The context-by-direction interaction is estimated almost entirely by the supplementary arm and does not survive its removal. Context effects are therefore exploratory. The original gene-level any-PASS endpoint was also sensitive to the number of evaluable pairs; this was addressed by adjustment and 5,000 equal-opportunity resamples, but the outcome remains discrete and does not itself model all dependence among repeated pairs."""

NEW_LIMITATION_6 = """6. **The composite score is weight-sensitive and does not validate.** Only ACE2 is stable across the 2,005 weightings tested. The top-10 pair-level signal is also imprecise under opportunity-balanced resampling. Downstream communication therefore uses Pareto-front membership rather than the composite score, and its benchmark AUROC is reported on both layers (0.524 frozen, 0.447 SC2-03-inclusive) rather than as a single figure."""

NEW_METHOD_PARAGRAPH = """Because the original gene-level endpoint (validation in at least one same-pathogen holdout) is sensitive to the number of evaluable pairs, opportunity dependence was quantified and addressed in two ways. First, logistic models of the any-PASS outcome included log(1 + number of same-pathogen evaluable pairs); the pathogen-count and unit-count models were fitted separately because those features are alternative ranking rules. Second, an equal-opportunity resampling analysis selected exactly one same-pathogen pair at random for every gene in each of 5,000 replicates (seed 20261021). AUROCs of the frozen HG3 score and the simple discovery features were recomputed in each replicate; paired differences were taken within the same sampled outcome. The median and 2.5th–97.5th percentiles of the resampling distribution are reported. This analysis changes neither the frozen candidates nor the adjudicated holdout verdicts (**Supplementary Tables S21a and S21b**)."""

NEW_CONTEXT_METHOD = """The comparison between cellular context, intervention direction and the composite score is reported in **Supplementary Table S20**. Context exposure is defined per gene as the fraction of its evaluable same-pathogen pairs whose holdout matched the discovery cellular context, computed once over all pairs and once over the frozen pairs only. Pair-level effects were estimated by logistic regression of the PASS indicator on the context flag, intervention direction and Bio-Tier 2 status, marginally and with holdout fixed effects; the fixed-effect fit excludes HO-V-IAV-02, where separation occurs. Because context is constant within four of five frozen holdouts, the marginal and fixed-effect context estimates are exploratory rather than design-identified. Direction was coded preserve/activate versus inhibit. Both layers are reported separately."""

NEW_FIG7 = """**Figure 7 | Holdout adjudication and opportunity-balanced benchmarking.** (**a**) Pre-registered holdout PASS rate by stratum; the dashed line marks the overall rate. (**b**) Top-*K* enrichment of holdout PASS at the pair level (top 10: 9 of 25 pairs, 36%). The top-rank signal is limited and imprecise after opportunity-balanced resampling. (**c**) AUROC of the HG3 composite score against naive consensus features, at the gene level over 333 evaluable candidates, computed separately on the frozen pre-registered outcome and on the SC2-03-inclusive outcome; bars are the frozen pre-registered layer and the hollow marker is the SC2-03-inclusive layer (**Supplementary Table S19f**). (**d**) Equal-opportunity sensitivity for the frozen outcome: one same-pathogen pair was sampled per gene in each of 5,000 replicates. Points and intervals show the median paired AUROC difference from HG3. The apparent pathogen-breadth advantage was not robust (difference +0.072, 95% interval −0.055 to +0.194), whereas independent-unit breadth and hit frequency retained a positive advantage (**Supplementary Tables S21a and S21b**)."""


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
