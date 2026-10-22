# Figure 6 — Host-directed tractability and baseline-fitness liability
## Shared theme defines ROOT, the warm palette, PT() and save_pub().
source("R/hpicda_theme.R")

ax <- fread(file.path(ROOT, "21_stage14_druggability/stage14_druggability_axis.tsv"))

## ------------------------------------------------- (a) baseline-fitness plane
## The x axis is DepMap baseline fitness, not a measured therapeutic window: the
## panel reads the same liability direction as Fig. 2, so the annotations name
## the baseline-essentiality side rather than a window.
w <- ax[!is.na(min_baseline_chronos)]
w[, tier := factor(bio_tier, levels = HP$tier_levels)]
lab_genes <- c("ACE2", "EXT1", "SLC35A1", "BAP1", "FER", "CHMP4B")
pa <- ggplot(w, aes(min_baseline_chronos, total_score)) +
  geom_vline(xintercept = ESS, linetype = "dashed", linewidth = 0.35, colour = HP$muted) +
  geom_point(aes(colour = tier, size = tier, alpha = tier)) +
  ## the two side annotations occupy the empty band above y = 8.5 (the highest
  ## data point is at y = 8.0); ylim keeps repelled labels out of that band so
  ## they cannot land on top of the annotations
  geom_text_repel(data = w[gene %in% lab_genes], aes(label = gene), size = PT(6.6),
                  fontface = "bold", colour = HP$ink, min.segment.length = 0,
                  segment.size = 0.25, box.padding = 0.3, max.overlaps = 30, seed = 5,
                  ylim = c(-0.9, 8.5)) +
  ## short two-word labels: at this panel width anything longer than ~10
  ## characters per line runs past the panel edge
  annotate("text", x = ESS - 0.06, y = 9.3, label = "baseline\nessential",
           size = PT(6.6), colour = HP$muted, fontface = "italic", hjust = 1) +
  annotate("text", x = ESS + 0.06, y = 9.3, label = "baseline\nnon-essential",
           size = PT(6.6), colour = HP$muted, fontface = "italic", hjust = 0) +
  scale_colour_manual(values = HP$tier_cols, name = NULL) +
  scale_size_manual(values = c(1.6, 1.0, 0.75, 0.65), guide = "none") +
  scale_alpha_manual(values = c(0.95, 0.65, 0.45, 0.35), guide = "none") +
  coord_cartesian(xlim = c(-3.6, 1.0), ylim = c(-0.8, 9.8), clip = "off") +
  labs(x = "Min. baseline Chronos (more negative = more essential)",
       y = "HG3 total score",
       title = "Baseline cellular-fitness liability\nacross the ranked candidates") +
  theme_hp() +
  ## the legend sits outside the panel: at this panel height an inside legend
  ## would cover ~70 of the 343 points in the lower-right quadrant
  theme(legend.position = "right", legend.key.size = unit(3.0, "mm"),
        legend.spacing.y = unit(0.3, "mm"), legend.text = element_text(size = 6.8),
        plot.margin = margin(5.5, 2, 5.5, 2))

## ---------------------------------------------------- (b) druggability class
DCLS <- c(launched_drug_exists = "Launched drug",
          clinical_stage_compound = "Clinical-stage compound",
          tractable_no_clinical_compound = "Tractable, no clinical compound",
          no_tractability_evidence = "No tractability evidence")
ax[, dcls := factor(DCLS[druggability_class], levels = unname(DCLS))]
DIRL <- c(INHIBIT = "Inhibit", PRESERVE_ACTIVATE = "Preserve / activate")
pb_dt <- ax[, .N, by = .(dcls, intervention_direction)]
pb_dt[, dir := factor(DIRL[intervention_direction], levels = rev(unname(DIRL)))]
DCOL <- c("Launched drug" = HP$exploited, "Clinical-stage compound" = HP$protective,
          "Tractable, no clinical compound" = HP$state,
          "No tractability evidence" = HP$neutral)
## dark labels on the two pale fills
pb_dt[, lab_col := fifelse(dcls %in% c("No tractability evidence",
                                       "Tractable, no clinical compound"),
                           HP$ink, "white")]
pb <- ggplot(pb_dt, aes(N, dir, fill = dcls)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = N, colour = lab_col), position = position_stack(vjust = 0.5),
            size = PT(ANNOT_PT), fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = DCOL, name = NULL) +
  scale_colour_identity() +
  ## a hair of right-hand expansion so the last stack segment is not flush with
  ## the panel edge (which reads as a clipped bar)
  scale_x_continuous(expand = expansion(mult = c(0, 0.025))) +
  labs(x = "Frozen candidates", y = NULL,
       title = "Open Targets tractability") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.2),
        legend.position = "bottom", legend.key.size = unit(3.0, "mm"),
        legend.text = element_text(size = 6.6), legend.spacing.y = unit(0.2, "mm"),
        plot.margin = margin(5.5, 3, 5.5, 2))

## ------------------------------------------------------------------- compose
fig <- pa / pb +
  plot_layout(heights = c(0.62, 0.38)) +
  plot_annotation(tag_levels = c("A", "B")) &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold"))

save_pub(fig, "fig6_tractability", 183, 118)
