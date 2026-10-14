# Figure 3 — Cross-pathogen functional dependency map
## Shared figure theme (defines ROOT, the warm palette, PT() and save_pub()).
## Run from the repository root, or set HPICDA_ROOT first (see README).
source("R/hpicda_theme.R")

d <- load_hpicda()
frozen <- d$frozen; conv <- d$conv; modules <- d$modules

## gene-level convergence, restricted to the 353 frozen candidates
m <- merge(conv, frozen[, .(gene, bio_tier, total_score, convergent_module_memberships)],
           by = "gene")
m[, tier_r := match(bio_tier, HP$tier_levels)]
setorder(m, tier_r, -total_score)
m[, gene := factor(gene, levels = rev(gene))]

COLS <- c("p_SC2_EXP", "p_IAV_EXP", "p_RSV_EXP",
          "p_SC2_PRO", "p_IAV_PRO", "p_RSV_PRO")
CMP_LAB <- c(p_SC2_EXP = "SARS-CoV-2\nexploited", p_IAV_EXP = "Influenza A\nexploited",
             p_RSV_EXP = "RSV\nexploited", p_SC2_PRO = "SARS-CoV-2\nprotective",
             p_IAV_PRO = "Influenza A\nprotective", p_RSV_PRO = "RSV\nprotective")
ml <- melt(m, id.vars = c("gene", "bio_tier", "total_score",
                          "convergent_module_memberships"),
           measure.vars = COLS, variable.name = "cmp", value.name = "p")
ml[, value := pmin(pmax(-log10(pmax(p, 1e-12)), 0), 8)]
ml[, cmp := factor(cmp, levels = COLS, labels = CMP_LAB[COLS])]

## ------------------------------------------------- (a) candidate-level heatmap
grid_a <- as.data.table(expand.grid(gene = levels(m$gene), cmp = levels(ml$cmp)))
grid_a <- merge(grid_a, ml[, .(gene, cmp, value, bio_tier)],
                by = c("gene", "cmp"), all.x = TRUE)
grid_a[, gene := factor(gene, levels = levels(m$gene))]
grid_a[, cmp := factor(cmp, levels = levels(ml$cmp))]

## module-membership strip, aligned to the same gene order and drawn inside the
## same panel (continuous x) so it does not consume a patchwork tag
MODS <- c(HS = "CURATED:heparan-sulfate biosynthesis",
          `V-ATPase` = "CURATED:V-ATPase proton pump",
          COG = "CURATED:COG complex")
strip <- rbindlist(lapply(names(MODS), function(nm)
  data.table(gene = m$gene, mod = nm,
             on = as.integer(grepl(MODS[[nm]], m$convergent_module_memberships,
                                   fixed = TRUE)))))
strip[, gene := factor(gene, levels = levels(m$gene))]
strip[, mod := factor(mod, levels = names(MODS))]

GAP <- 0.75; SW <- 1.0
grid_a[, x := as.numeric(cmp)]
strip[, x := length(COLS) + GAP + as.numeric(mod) * SW]

X_BRK <- c(seq_along(COLS), length(COLS) + GAP + seq_along(MODS) * SW)
X_LAB <- c(CMP_LAB[COLS], names(MODS))

rng <- range(grid_a$value, na.rm = TRUE)
pa <- ggplot() +
  geom_tile(data = grid_a, aes(x, gene, fill = value)) +
  geom_tile(data = strip, aes(x, gene, alpha = factor(on)), fill = HP$strip_on,
            width = SW * 0.92) +
  geom_vline(xintercept = length(COLS) + GAP / 2, linewidth = 0.4, colour = HP$rule) +
  geom_hline(yintercept = sum(m$bio_tier == "Bio-Tier 2") + 0.5,
             linewidth = 0.4, colour = HP$rule) +
  scale_fill_gradientn(colours = HP$seq_warm, limits = rng, na.value = "white",
                       name = "-log10 p") +
  scale_alpha_manual(values = c(`0` = 0.10, `1` = 1), guide = "none") +
  scale_x_continuous(expand = c(0, 0), breaks = X_BRK, labels = X_LAB) +
  scale_y_discrete(expand = c(0, 0), breaks = NULL) +
  labs(x = NULL, y = NULL, title = "Candidate convergence map (n = 353)") +
  theme_hp() +
  # heatmap panel: no axis text/ticks (rows are unnamed), but the boxed panel
  # frame is retained so the panel does not read as an unfinished plot
  theme(axis.line = element_line(linewidth = 0.28, colour = HP$ink),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = DENSE_PT, angle = 90, hjust = 0, vjust = 0.5),
        legend.position = "right",
        legend.key.width = unit(1.8, "mm"), legend.key.height = unit(7, "mm"),
        legend.title = element_text(size = 7.0, angle = 90, vjust = 0.5),
        legend.text = element_text(size = 6.6),
        plot.margin = margin(6, 2, 6, 3))

## ---------------------------------------------------- (b) Bio-Tier 2 detail
t2 <- m[bio_tier == "Bio-Tier 2"]
t2[, gene := factor(gene, levels = rev(gene[order(total_score)]))]
g2 <- as.data.table(merge(expand.grid(gene = levels(t2$gene), cmp = levels(ml$cmp)),
                          ml[, .(gene, cmp, value)],
                          by = c("gene", "cmp"), all.x = TRUE))
g2[, gene := factor(gene, levels = levels(t2$gene))]
g2[, cmp := factor(cmp, levels = levels(ml$cmp))]
pb <- ggplot(g2, aes(cmp, gene, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.1f", value))),
            size = PT(DENSE_PT),
            colour = ifelse(!is.na(g2$value) & g2$value > 4.2, "white", HP$ink)) +
  scale_fill_gradientn(colours = HP$seq_warm, limits = rng, na.value = "white",
                       guide = "none") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Bio-Tier 2 primary candidates") +
  theme_hp() +
  theme(axis.line = element_line(linewidth = 0.28, colour = HP$ink),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = DENSE_PT, angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 7.0),
        plot.margin = margin(6, 3, 6, 3))

## -------------------------------------------------- (c) curated module level
cur <- modules[grepl("^CURATED:", module)]
cur[, lab := sub("^CURATED:", "", module)]
setorder(cur, fisher_EXP)
cur[, lab := factor(lab, levels = lab)]
pc <- ggplot(cur, aes(fisher_EXP, lab)) +
  geom_col(fill = HP$state, width = 0.58) +
  geom_text(aes(label = sprintf("%s  (q = %.0e, n = %d)", format(round(fisher_EXP),
            big.mark = ","), q_EXP, n_genes)),
            hjust = -0.04, size = PT(ANNOT_PT), colour = HP$ink) +
  scale_x_continuous(limits = c(0, max(cur$fisher_EXP) * 2.6), expand = c(0, 0)) +
  labs(x = "Fisher combined evidence (EXP)", y = NULL,
       title = "Curated module convergence") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.4))

## ------------------------------------------------------------------- compose
right <- pb / pc + plot_layout(heights = c(1.0, 0.75))
fig <- (pa | right) +
  plot_layout(widths = c(0.72, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold"))

save_pub(fig, "fig3_cross_pathogen_map", 183, 150)
