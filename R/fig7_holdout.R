# Figure 7 — Independent holdout adjudication and honest benchmarking
## Shared figure theme (defines ROOT, the warm palette, PT() and save_pub()).
## Run from the repository root, or set HPICDA_ROOT first (see README).
source("R/hpicda_theme.R")

topk <- fread(file.path(ROOT, "20_stage13_holdout_adjudication/stage13_topk_enrichment.tsv"))
bench <- fread(file.path(ROOT, "24_closing_analyses/benchmark_two_layer.tsv"))
ctxd <- fread(file.path(ROOT, "24_closing_analyses/context_determinants_stratified.tsv"))

## ------------------------------------------- (a) PASS rate across strata (forest)
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

## ------------------------------------------------- (b) top-K enrichment of PASS
topk[, K := factor(K, levels = K)]
PEAK <- topk[which.max(odds_ratio)]
pb <- ggplot(topk, aes(K, odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35, colour = HP$muted) +
  geom_col(aes(fill = gene_recall_topk), width = 0.5) +
  geom_text(aes(label = sprintf("OR %.2f\np = %s", odds_ratio,
                                formatC(fisher_p_one_sided, format = "g", digits = 2))),
            vjust = -0.35, size = PT(6.6), colour = HP$ink, lineheight = 1.1) +
  scale_fill_gradientn(colours = HP$seq_warm, limits = c(0, 0.8), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.46))) +
  labs(x = "Top-K of the frozen ranking", y = "PASS enrichment (odds ratio)",
       title = "Enrichment is confined to the top decile") +
  theme_hp() +
  theme(plot.margin = margin(5.5, 3, 5.5, 2))

## ---------------------------------------- (c) benchmark against naive rankings
## Bars are the frozen, pre-registered layer (the primary analysis); the hollow
## marker is the SC2-03-inclusive layer.  A single-layer panel would carry only
## the supplementary-inclusive AUROCs, which is where the composite looks
## anti-predictive -- on the frozen layer it is merely at chance.
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

## ---------------------------------------------- (d) context determinants of PASS
ctxd[, ctx := factor(ifelse(cm == 1, "Context-matched", "Context-mismatched"),
                     levels = c("Context-mismatched", "Context-matched"))]
ctxd[, dir := factor(ifelse(direction == "INHIBIT", "Inhibit", "Preserve / activate"),
                     levels = c("Inhibit", "Preserve / activate"))]
pd <- ggplot(ctxd, aes(ctx, pass_rate, fill = dir)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62,
           colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d/%d)", pass_rate * 100, n_pass, n_eval)),
            position = position_dodge(width = 0.7), vjust = -0.28, size = PT(6.6),
            colour = HP$ink, lineheight = 1.1) +
  scale_fill_manual(values = c(Inhibit = HP$state, `Preserve / activate` = HP$exploited),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 0.9), expand = c(0, 0)) +
  labs(x = NULL, y = "Holdout PASS rate",
       title = "Matched context + activation direction") +
  theme_hp() +
  theme(axis.text.x = element_text(size = 7.0),
        legend.position = "bottom", legend.key.size = unit(3.0, "mm"),
        legend.text = element_text(size = 6.8),
        plot.margin = margin(5.5, 3, 5.5, 2))

## ------------------------------------------------------------------- compose
fig <- (pa | pb) / (bc | pd) +
  plot_layout(heights = c(1.0, 0.80)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold"))

save_pub(fig, "fig7_holdout", 183, 140)
