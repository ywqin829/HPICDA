# Figure 4 — Human pulmonary and genetic context of the frozen candidates
## Shared figure theme (defines ROOT, the warm palette, PT() and save_pub()).
## Run from the repository root, or set HPICDA_ROOT first (see README).
source("R/hpicda_theme.R")

d <- load_hpicda()
frozen <- d$frozen
s8   <- fread(file.path(ROOT, "09_transcriptomic_context/stage8_per_gene_summary.tsv"))
s9ct <- fread(file.path(ROOT, "10_singlecell_context/stage9_hlca_candidate_expression_by_celltype.tsv"))
s10  <- fread(file.path(ROOT, "08_human_genetics/stage10_genetics_evidence.tsv"))

## ------------------------------------------------- (a) transcriptomic evidence
## A compact two-row summary, not a full-width banner: the four-category
## composition is one summary row per direction and does not need the width.
e <- merge(frozen[, .(gene, intervention_direction, bio_tier)],
           s8, by = "gene", all.x = TRUE)
e[, same_pathogen_deg := (n_sig_up_same + n_sig_down_same) > 0]
e[, cat := fifelse(clinical_de_concordant,
                   "Clinical DE concordant",
            fifelse(same_pathogen_deg, "Same-pathogen DEG",
            fifelse(any_sig_any_dataset, "Any-dataset DEG", "No transcriptomic signal")))]
CATS <- c("Clinical DE concordant", "Same-pathogen DEG",
          "Any-dataset DEG", "No transcriptomic signal")
e[, cat := factor(cat, levels = CATS)]
pa_dt <- e[, .N, by = .(intervention_direction, cat)]
DIRLAB <- c(INHIBIT = "Inhibit (host dependency)",
            PRESERVE_ACTIVATE = "Preserve / activate")
pa_dt[, dir := factor(DIRLAB[intervention_direction],
                      levels = rev(unname(DIRLAB)))]
COLS_A <- c("Clinical DE concordant" = HP$exploited,
            "Same-pathogen DEG"      = HP$protective,
            "Any-dataset DEG"        = HP$state,
            "No transcriptomic signal" = HP$neutral)
## white labels are unreadable on the pale "no signal" / olive fills
pa_dt[, lab_col := fifelse(cat %in% c("No transcriptomic signal", "Any-dataset DEG"),
                           HP$ink, "white")]
DIRN <- pa_dt[, .(t = sum(N)), by = intervention_direction][, max(t)]
pa <- ggplot(pa_dt, aes(N, dir, fill = cat)) +
  geom_col(width = 0.46, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(N >= 12, N, ""), colour = lab_col),
            position = position_stack(vjust = 0.5),
            size = PT(ANNOT_PT), fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = COLS_A, name = NULL) +
  scale_colour_identity() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, DIRN * 1.02)) +
  scale_y_discrete(labels = function(x) gsub(" \\(.*", "", x)) +
  labs(x = "Frozen candidates (n = 353)", y = NULL,
       ## "corroboration" overstates what a composition count can show: these
       ## datasets supply context and a direction check, not confirmation.
       title = "Transcriptomic context and directional concordance\nin human lung datasets") +
  theme_hp() +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.key.height = unit(3.4, "mm"), legend.key.width = unit(3.6, "mm"),
        legend.text = element_text(size = 6.6),
        axis.text.y = element_text(size = 7.0),
        plot.margin = margin(5.5, 3, 5.5, 2))

## ------------------------------------------------ (b) lung cell-type expression
## Spans the full figure width so the 19 cell-type columns are wide enough to
## carry their three-digit values.
LUNG_CT <- c(
  "pulmonary alveolar type 1 cell" = "AT1",
  "pulmonary alveolar type 2 cell" = "AT2",
  "epithelial cell of lung"        = "Lung epithelium",
  "ciliated cell"                  = "Ciliated",
  "club cell"                      = "Club",
  "lung goblet cell"               = "Goblet",
  "respiratory basal cell"         = "Basal",
  "macrophage"                     = "Macrophage",
  "inflammatory macrophage"        = "Inflam. macrophage",
  "CD14-positive monocyte"         = "Monocyte",
  "myeloid dendritic cell"         = "Myeloid DC",
  "lung pericyte"                  = "Pericyte",
  "bronchus fibroblast of lung"    = "Fibroblast",
  "pulmonary artery endothelial cell" = "Pulm. artery EC",
  "endothelial cell of vascular tree" = "Vascular EC",
  "lung resident memory CD8-positive, alpha-beta T cell" = "TRM CD8",
  "lung resident memory CD4-positive, alpha-beta T cell" = "TRM CD4",
  "IgA plasmablast"                = "IgA plasmablast",
  "mast cell"                      = "Mast")
t2 <- frozen[bio_tier == "Bio-Tier 2", gene]
hb <- s9ct[gene %in% t2 & cell_type %in% names(LUNG_CT)]
hb <- merge(hb, frozen[, .(gene, intervention_direction, total_score)], by = "gene")
hb[, gene_lab := sprintf("%s (%s)", gene,
                         fifelse(intervention_direction == "INHIBIT", "inh", "pre"))]
setorder(hb, -total_score)
GL <- unique(hb$gene_lab)
hb[, gene_lab := factor(gene_lab, levels = rev(GL))]
hb[, ct := factor(LUNG_CT[cell_type], levels = unname(LUNG_CT))]
hb[, pct := pmin(pct_expressing, 60)]
pb <- ggplot(hb, aes(ct, gene_lab, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct >= 2, sprintf("%.0f", pct_expressing), "")),
            size = PT(6.5), colour = ifelse(hb$pct > 14, "white", HP$ink)) +
  scale_fill_gradientn(colours = HP$seq_warm, limits = c(0, 60), na.value = "white",
                       breaks = c(0, 20, 40, 60),
                       name = "% cells\nexpressing", guide = guide_colourbar(
                         barwidth = unit(1.8, "mm"), barheight = unit(13, "mm"))) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Bio-Tier 2 expression across lung cell types") +
  theme_hp() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = DENSE_PT, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7.0),
        legend.position = "right",
        legend.title = element_text(size = 6.8),
        legend.text = element_text(size = 6.6),
        plot.margin = margin(5.5, 2, 5.5, 2))

## ------------------------------------------------------- (c) human genetics
## Boxplot only. An overlaid jitter layer read as two plots stacked on top of
## each other, so the distribution is carried by the box and the group sizes by
## the tick labels.
g <- s10[!is.na(loeuf)]
g <- merge(g, frozen[, .(gene, bio_tier, intervention_direction)], by = "gene")
g[, tier := factor(tier, levels = c("G3", "G2", "G0"))]
TCOL <- c(G3 = HP$exploited, G2 = HP$protective, G0 = HP$neutral)
g[, is_t2 := bio_tier == "Bio-Tier 2"]
## plain character keyed by tier: an earlier version made the label vector a
## factor, and ggplot2 then fell back to the factor codes, printing "1 2 3" as
## the tick labels instead of "G3 (n=..)".
NLAB <- g[, .(lab = sprintf("%s\n(n=%d)", as.character(tier), .N)), by = tier]
set.seed(4)
pc <- ggplot(g, aes(tier, loeuf)) +
  geom_hline(yintercept = 0.35, linetype = "dashed", linewidth = 0.3, colour = HP$muted) +
  geom_boxplot(aes(colour = tier), width = 0.46, outlier.shape = NA, fill = NA,
               linewidth = 0.35) +
  geom_point(data = g[is_t2 == TRUE], aes(tier, loeuf), shape = 21,
             size = 1.7, stroke = 0.35, colour = HP$ink, fill = HP$tier2) +
  geom_text_repel(data = g[is_t2 == TRUE], aes(label = gene), size = PT(6.6),
                  fontface = "bold", colour = HP$ink, min.segment.length = 0,
                  segment.size = 0.25, box.padding = 0.3, max.overlaps = 30, seed = 3) +
  scale_colour_manual(values = TCOL, guide = "none") +
  scale_x_discrete(labels = setNames(NLAB$lab, as.character(NLAB$tier))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.16))) +
  labs(x = "Constraint tier (G3 = most constrained)",
       y = "LOEUF (lower = more constrained)",
       title = "Human loss-of-function constraint") +
  theme_hp() +
  theme(axis.text.x = element_text(size = 6.8, lineheight = 1.1),
        plot.margin = margin(5.5, 3, 5.5, 2))

## ------------------------------------------------------------------- compose
## a and c share the top row at half width each; b spans the full figure width
## so its 19 cell-type columns are wide enough for three-digit values.
## Tags are given explicitly because the panel order is not a, b, c.
top <- (pa | pc) + plot_layout(widths = c(1.0, 1.0))
fig <- top / pb +
  plot_layout(heights = c(0.42, 0.58)) +
  plot_annotation(tag_levels = list(c("A", "C", "B"))) &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold"))

save_pub(fig, "fig4_human_context", 183, 138)
