# Figure 7 — Independent holdout adjudication and opportunity-balanced
# benchmarking. Shared paths, palette and export helpers come from the theme.
source("R/hpicda_theme.R")

topk <- fread(file.path(ROOT, "20_stage13_holdout_adjudication/stage13_topk_enrichment.tsv"))
bench <- fread(file.path(ROOT, "24_closing_analyses/benchmark_two_layer.tsv"))
opp <- fread(file.path(ROOT, "25_manuscript_package/tables/S21b_equal_opportunity_summary.tsv"))
NUMCOLS <- c("median", "ci_lo", "ci_hi", "n", "fraction_at_or_below_zero")
opp[, (NUMCOLS) := lapply(.SD, as.numeric), .SDcols = NUMCOLS]

## ------------------------------------------- (a) PASS rate across strata
STR <- rbindlist(list(
  data.table(grp = "Overall (same-pathogen)", n = 654, k = 55),
  data.table(grp = "Bio-Tier 2",             n = 29,  k = 9),
  data.table(grp = "Bio-Tier 3",             n = 449, k = 30),
  data.table(grp = "ARCHIVE",                n = 147, k = 11),
  data.table(grp = "REJECT",                 n = 29,  k = 5),
  data.table(grp = "Context-matched",        n = 85,  k = 22),
  data.table(grp = "Context-mismatched",     n = 569, k = 33),
  data.table(grp = "Bio-Tier 2 (matched)",   n = 8,   k = 3),
  data.table(grp = "Rest (matched)",         n = 77,  k = 19),
  data.table(grp = "HS cluster (matched)",   n = 16,  k = 6)))
STR[, p := k / n]
BASE <- STR[grp == "Overall (same-pathogen)", p]
STR[, family := c("All", "Tier", "Tier", "Tier", "Tier", "Context", "Context",
                  "Matched", "Matched", "Matched")]
STR[, grp := factor(grp, levels = rev(grp))]
STR[, fam := factor(family, levels = c("All", "Tier", "Context", "Matched"))]
NCOL <- c(All = HP$ink, Tier = HP$state, Context = HP$protective, Matched = HP$exploited)
pa <- ggplot(STR, aes(p, grp)) +
  geom_vline(xintercept = BASE, linetype = "dashed", linewidth = 0.35, colour = HP$muted) +
  geom_segment(aes(x = 0, xend = p, yend = grp, colour = fam), linewidth = 0.4) +
  geom_point(aes(colour = fam, size = n)) +
  geom_text(aes(label = sprintf("%d/%d = %.0f%%", k, n, p * 100)), hjust = -0.2,
            size = PT(6.6), colour = HP$ink) +
  scale_colour_manual(values = NCOL, guide = "none") +
  scale_size_continuous(range = c(1.5, 3.2), guide = "none") +
  scale_x_continuous(limits = c(0, 0.95), expand = c(0, 0)) +
  labs(x = "Holdout PASS rate", y = NULL,
       title = "Pre-registered holdout: PASS rate by stratum") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 6.9),
        plot.margin = margin(5.5, 3, 5.5, 2))

## ------------------------------------------------- (b) top-K enrichment
topk[, K := factor(K, levels = K)]
pb <- ggplot(topk, aes(K, odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35, colour = HP$muted) +
  geom_col(aes(fill = gene_recall_topk), width = 0.5) +
  geom_text(aes(label = sprintf("OR %.2f\np = %s", odds_ratio,
                                formatC(fisher_p_one_sided, format = "g", digits = 2))),
            vjust = -0.35, size = PT(6.6), colour = HP$ink, lineheight = 1.1) +
  scale_fill_gradientn(colours = HP$seq_warm, limits = c(0, 0.8), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.46))) +
  labs(x = "Top-K of the frozen ranking", y = "PASS enrichment (odds ratio)",
       title = "Top-rank enrichment is limited and imprecise") +
  theme_hp() +
  theme(plot.margin = margin(5.5, 3, 5.5, 2))

## ---------------------------------------- (c) two-layer ranking benchmark
NICE <- c(HPICDA_HG3_total_score = "HPICDA HG3 score",
          hit_frequency = "Hit frequency",
          n_units_hit = "n units hit",
          n_pathogens_hit = "n pathogens hit",
          n_screens_hit = "n screens hit",
          best_rank_pct = "Best rank pct")
bench[, lab := factor(NICE[ranking], levels = rev(unname(NICE)))]
bench[, above := auroc_frozen > 0.5]
bc <- ggplot(bench, aes(auroc_frozen, lab)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.35, colour = HP$muted) +
  geom_col(aes(fill = above), width = 0.52) +
  geom_point(aes(colour = "SC2-03-inclusive"), shape = 21, size = 1.5,
             fill = "white", stroke = 0.35) +
  geom_text(aes(label = sprintf("%.3f", auroc_frozen)), hjust = -0.22,
            size = PT(ANNOT_PT), colour = HP$ink) +
  scale_fill_manual(values = c(`FALSE` = HP$confounded, `TRUE` = HP$protective),
                    guide = "none") +
  scale_colour_manual(values = c(`SC2-03-inclusive` = HP$ink), name = NULL) +
  scale_x_continuous(limits = c(0, 0.86), expand = c(0, 0)) +
  labs(x = "AUROC for independent validation", y = NULL,
       title = "The composite score is at chance, not above it") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.0),
        legend.position = "bottom", legend.key.size = unit(3.0, "mm"),
        legend.text = element_text(size = 6.8),
        legend.margin = margin(t = -1.5, b = -1.5),
        plot.margin = margin(5.5, 3, 5.5, 2))

## -------------------------- (d) equal-opportunity paired AUROC differences
opp_d <- opp[section == "paired_difference_vs_HG3"]
setnames(opp_d, "metric", "comparison")
COMP_LABELS <- c(
  n_pathogens_hit_minus_HG3 = "Pathogen breadth",
  n_units_hit_minus_HG3 = "Independent units hit",
  hit_frequency_minus_HG3 = "Hit frequency")
opp_d[, comparison_label := factor(
  unname(COMP_LABELS[as.character(comparison)]),
  levels = c("Pathogen breadth", "Independent units hit", "Hit frequency"))]
opp_d[, robust := fraction_at_or_below_zero < 0.025]
opp_d[, ypos := match(comparison_label, c("Pathogen breadth", "Independent units hit", "Hit frequency"))]
pd <- ggplot(opp_d, aes(median, ypos, colour = robust)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = HP$muted) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = ypos, yend = ypos), linewidth = 0.45) +
  geom_point(size = 1.8) +
  geom_text(aes(label = sprintf("%+.3f [%.3f, %.3f]", median, ci_lo, ci_hi)),
            vjust = -0.9, size = PT(6.6), colour = HP$ink) +
  scale_colour_manual(values = c(`FALSE` = HP$confounded, `TRUE` = HP$exploited),
                      guide = "none") +
  scale_x_continuous(limits = c(-0.07, 0.23),
                     breaks = seq(-0.05, 0.20, by = 0.05)) +
  scale_y_continuous(breaks = 1:3,
                     labels = c("Pathogen breadth", "Independent units hit", "Hit frequency"),
                     limits = c(0.75, 3.35)) +
  labs(x = "Paired AUROC difference versus HG3", y = NULL,
       title = "Opportunity-balanced comparison") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.0),
        plot.margin = margin(5.5, 3, 5.5, 2))

fig <- (pa | pb) / (bc | pd) +
  plot_layout(heights = c(1.0, 0.80)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold"))

save_pub(fig, "fig7_holdout", 183, 140)
