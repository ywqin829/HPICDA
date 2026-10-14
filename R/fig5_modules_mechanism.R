# Figure 5 — Shared host modules: co-essentiality, physical support, contradictions
## Shared figure theme (defines ROOT, the warm palette, PT() and save_pub()).
## Run from the repository root, or set HPICDA_ROOT first (see README).
source("R/hpicda_theme.R")

mc <- fread(file.path(ROOT, "11_modules/stage11_module_coessentiality.tsv"))
cp <- fread(file.path(ROOT, "11_modules/stage11_contradictory_pairs_stratified.tsv"))
ia <- fread(file.path(ROOT, "11_modules/stage11_intact_support.tsv"))

## ------------------------------------------------ (a) module co-essentiality
mc[, curated := grepl("^CURATED:", module)]
mc[, lab := sub("^Reactome:", "", module)]
mc[, lab := sub("^CURATED:", "", lab)]
mc[, lab := gsub('"', "", lab)]
## Two-line axis labels, capped at two lines: a third line needs more row pitch
## than 13 rows leave, and unwrapped names push the right-hand column to a
## sliver. Names too long for two lines are elided; full names stay in the
## supplementary module table.
wrap2 <- function(s, width = 28) {
  w <- strwrap(s, width = width)
  if (length(w) > 2) {
    w <- c(w[1], paste0(substr(paste(w[-1], collapse = " "), 1, width - 3), "..."))
  }
  paste(w, collapse = "\n")
}
TOP <- mc[order(-mean_pairwise_r)][1:13]
TOP[, lab := vapply(lab, wrap2, character(1))]
TOP[, lab := factor(lab, levels = rev(lab))]
TOP[, grp := fifelse(curated, "Curated convergent module", "Reactome pathway")]
GCOL <- c("Curated convergent module" = HP$exploited, "Reactome pathway" = HP$state)
bgm <- mc$bg_mean[1]; bgs <- mc$bg_sd[1]
pa <- ggplot(TOP, aes(mean_pairwise_r, lab, fill = grp)) +
  geom_col(width = 0.66) +
  geom_errorbar(aes(xmin = bgm - 3 * bgs, xmax = bgm + 3 * bgs), width = 0,
                orientation = "y", colour = HP$muted, linewidth = 0.35) +
  geom_vline(xintercept = bgm, linewidth = 0.35, colour = HP$muted, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.2f", mean_pairwise_r)), hjust = -0.22,
            size = PT(ANNOT_PT), colour = HP$ink) +
  scale_fill_manual(values = GCOL, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = "Mean pairwise co-essentiality (r)", y = NULL,
       title = "Modules with coherent DepMap co-essentiality") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 6.9),
        legend.position = "bottom", legend.key.size = unit(3.0, "mm"),
        legend.text = element_text(size = 6.8),
        plot.margin = margin(5.5, 3, 5.5, 2))

## ------------------------------------------------------ (b) physical support
ia[, lab := sub("^CURATED:", "", module)]
ia[, short := c("COG complex", "V-ATPase", "Heparan-sulfate\nbiosynthesis")[
        match(lab, c("COG complex", "V-ATPase proton pump", "heparan-sulfate biosynthesis"))]]
ia[, density := n_intra_module_edges / n_members]
ia[, short := factor(short, levels = short[order(density)])]
pb <- ggplot(ia, aes(density, short)) +
  geom_segment(aes(x = 0, xend = density, yend = short), colour = HP$rule, linewidth = 0.35) +
  geom_point(aes(size = n_members), colour = HP$protective) +
  ## two-line label: a single line ("112 edges / 23 members") ran past the panel
  ## edge and was clipped; the label is now narrow enough to sit inside the
  ## expansion that the plot actually reserves.
  geom_text(aes(label = sprintf("%d edges\n%d members", n_intra_module_edges, n_members)),
            hjust = -0.18, vjust = 0.5, lineheight = 0.95,
            size = PT(6.6), colour = HP$ink) +
  scale_size_continuous(range = c(2.0, 3.8), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.85))) +
  labs(x = "PPI edges per member", y = NULL,
       title = "Physical interaction support") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.2, lineheight = 1.15),
        plot.margin = margin(5.5, 3, 5.5, 2))

## --------------------------------------------------- (c) contradictory pairs
## The source table has only four rows and RSV contributes none. A faceted
## stacked bar on a shared linear scale hid the two small strata (n = 1 and
## n = 7) beside n = 331, so the "same cell line" facet read as empty. The
## panel is now a log-scale lollipop: every stratum is visible at its true
## magnitude and nothing is stacked.
cp[, path := fifelse(pathogen == "Influenza", "Influenza A", pathogen)]
cp[, cls := fifelse(essentiality_confound_suspect, "Baseline-essentiality suspect",
                    "Not confounded")]
PA <- cp[, .(n_pairs = sum(n_pairs)), by = .(path, same_cell_line, cls)]
PA[, xlab := fifelse(same_cell_line, "same cell line", "different cell lines")]
PA[, row := factor(sprintf("%s\n(%s)", path, xlab),
                   levels = unique(sprintf("%s\n(%s)", path, xlab)[order(-n_pairs)]) )]
CC <- c("Baseline-essentiality suspect" = HP$confounded, "Not confounded" = HP$state)
N_TOT <- PA[, sum(n_pairs)]
FLOOR <- 0.75
pc <- ggplot(PA, aes(n_pairs, row)) +
  geom_segment(aes(x = FLOOR, xend = n_pairs, yend = row, colour = cls),
               linewidth = 0.4) +
  geom_point(aes(colour = cls, size = n_pairs)) +
  geom_text(aes(label = n_pairs), hjust = -0.45, size = PT(ANNOT_PT),
            colour = HP$ink) +
  scale_colour_manual(values = CC, name = NULL) +
  scale_size_continuous(range = c(1.8, 3.6), guide = "none") +
  scale_x_log10(limits = c(FLOOR, 1400), breaks = c(1, 10, 100, 1000)) +
  labs(x = "Pathogen-gene pairs (log scale)",
       y = NULL,
       title = "Discordant calls by cell-line context",
       subtitle = sprintf("n = %d pairs; RSV contributed none", N_TOT)) +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.0, lineheight = 1.05),
        legend.position = "bottom", legend.direction = "vertical",
        legend.key.size = unit(3.0, "mm"),
        legend.text = element_text(size = 6.8),
        plot.subtitle = element_text(size = 6.6, colour = HP$muted),
        plot.margin = margin(5.5, 3, 5.5, 2))

## ------------------------------------------------------------------- compose
## the 20-row ranking panel sets the figure height; the two right-hand panels are
## given height in proportion to their content, not stretched to match.
right <- pb / pc + plot_layout(heights = c(0.58, 1.42))
fig <- (pa | right) +
  plot_layout(widths = c(1.0, 1.08)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold"))

save_pub(fig, "fig5_modules_mechanism", 183, 132)
