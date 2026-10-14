# Figure 6 — Host-directed tractability and baseline-fitness liability
## Shared figure theme (defines ROOT, the warm palette, PT() and save_pub()).
## Run from the repository root, or set HPICDA_ROOT first (see README).
source("R/hpicda_theme.R")

ax <- fread(file.path(ROOT, "21_stage14_druggability/stage14_druggability_axis.tsv"))
cur <- fread(file.path(ROOT, "21_stage14_druggability/stage14_launched_curation.tsv"))

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
       title = "Tractability annotation (Orphanet / ChEMBL / OpenTargets)") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.2),
        legend.position = "bottom", legend.key.size = unit(3.0, "mm"),
        legend.text = element_text(size = 6.6), legend.spacing.y = unit(0.2, "mm"),
        plot.margin = margin(5.5, 3, 5.5, 2))

## ----------------------------------------- (c) launched-compound curation audit
CLS <- c(direct_target = "Direct target",
         direct_secondary = "Direct, secondary MOA",
         indirect_pathway = "Indirect / pathway-adjacent",
         endogenous_ligand = "Endogenous ligand",
         likely_misannotation = "Likely misannotation")
cur[, cls := factor(CLS[curation_class], levels = unname(CLS))]
cc <- cur[, .N, by = cls]
cc[, frac := N / sum(N)]
CCOL <- c("Direct target" = HP$tier2, "Direct, secondary MOA" = HP$tier3,
          "Indirect / pathway-adjacent" = HP$archive,
          "Endogenous ligand" = HP$state, "Likely misannotation" = HP$confounded)
## wrap the two long class names: unwrapped they consumed most of the panel and
## pushed the value labels past the right-hand edge, which clipped them.
cc[, clab := gsub(" / ", " /\n", as.character(cls))]
cc[, clab := gsub(", ", ",\n", clab)]
cc[, clab := factor(clab, levels = clab[order(N)])]
pc <- ggplot(cc, aes(N, clab, fill = cls)) +
  geom_col(width = 0.6) +
  ## labels sit on the real scale; white-on-fill is unreadable at these widths,
  ## so they stay dark and the axis reserves room via the expansion below
  geom_text(aes(label = sprintf("%d (%.0f%%)", N, frac * 100)), hjust = -0.10,
            size = PT(6.6), colour = HP$ink) +
  scale_fill_manual(values = CCOL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.42))) +
  ## title wrapped and axis label shortened: this panel is only ~130 pt wide and
  ## the single-line versions ran off the right-hand page edge
  labs(x = "Gene-drug pairs (n = 36)", y = NULL,
       title = "Human audit of\n\"launched drug\" annotations") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.0, lineheight = 1.05),
        plot.margin = margin(5.5, 3, 5.5, 2))

## ------------------------------------------------------------------- compose
## a used to run the full height of the figure for a scatter whose own data
## ranges are only about 1:2, which left a tall, half-empty left column. Panels
## a and c now share the top row (a scatter left, the five audit bars right) and
## the two-category stacked bar in b runs along the bottom, where its natural
## shape is wide and short. Tags are given explicitly because the reading order
## is a, c, b.
## widths belong to the inner one-row patchwork: on the outer (pa|pc)/pb the
## layout has a single column and an outer widths vector would be ignored.
top <- (pa | pc) + plot_layout(widths = c(1.0, 1.05))
fig <- top / pb +
  plot_layout(heights = c(0.63, 0.37)) +
  plot_annotation(tag_levels = list(c("A", "C", "B"))) &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold"))

save_pub(fig, "fig6_tractability", 183, 118)
