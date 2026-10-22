# Module-by-holdout specificity re-audit

Date: 2026-10-22

## Decision

Do not add the proposed module-specificity atlas to the manuscript.

The two source reports contain a denominator error and a verdict-definition
error:

1. `NOT_EVALUABLE_XPATH` rows were counted in the background denominator. These
   are genes with no same-pathogen evaluable observation in that holdout and
   should not be treated as failures.
2. `PASS_XPATH` was treated as a successful same-pathogen PASS. The manuscript
   defines PASS only among same-pathogen, evaluable observations.
3. The resulting Fisher tests treated genes as independent and did not correct
   for the multiple module-by-holdout comparisons.

## Correct eligibility rule

- retain `same_pathogen == True`;
- retain verdicts `PASS`, `FAIL` and `CONTEXT_BOUNDARY`;
- count only `PASS` as success;
- count `FAIL` and `CONTEXT_BOUNDARY` as non-success;
- test module membership versus all other evaluable genes within the same
  holdout;
- apply Benjamini-Hochberg correction across all non-empty tests.

## Corrected results

| Module | Holdout | Module PASS | Background PASS | Nominal one-sided P | BH q |
|---|---|---:|---:|---:|---:|
| HS | HO-V-RSV-01 | 5/5 | 4/8 | 0.098 | 0.490 |
| HS | HO-V-SC2-01 | 1/5 | 27/256 | 0.436 | 1.000 |
| HS | HO-V-SC2-02 | 0/6 | 8/322 | 1.000 | 1.000 |
| V-ATPase | HO-V-IAV-01 | 7/10 | 3/15 | 0.018 | 0.182 |
| V-ATPase | HO-V-IAV-02 | 0/10 | 0/17 | 1.000 | 1.000 |
| V-ATPase | HO-V-SC2-01 | 1/8 | 27/253 | 0.602 | 1.000 |
| V-ATPase | HO-V-SC2-02 | 0/10 | 8/318 | 1.000 | 1.000 |
| COG | HO-V-IAV-01 | 1/1 | 9/24 | 0.400 | 1.000 |
| COG | HO-V-IAV-02 | 0/1 | 0/26 | 1.000 | 1.000 |
| COG | HO-V-RSV-01 | 1/1 | 8/12 | 0.692 | 1.000 |

No comparison survives multiplicity correction.

## Why this was not in the original analysis

The original module analysis was designed around co-essentiality and physical
PPI support, not post-hoc module enrichment within each holdout. Adding
holdout-specific enrichment after seeing the results would constitute a new
exploratory analysis. The current frozen-candidate data are also sparse at the
holdout level: most candidate genes have no evaluable same-pathogen row in a
given holdout. A formal specificity claim would require a pre-specified module
enrichment design with an appropriate gene-level background and multiplicity
control.

## Reproduction

Script:

`24_closing_analyses/66_reaudit_module_holdout_specificity.py`

Outputs:

- `24_closing_analyses/66_module_holdout_specificity_reaudit.tsv`
- `24_closing_analyses/66_module_holdout_specificity_reaudit.json`
