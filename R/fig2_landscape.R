# Figure 2 — Baseline fitness versus infection-dependency landscape
## Shared figure theme (defines ROOT, the warm palette, PT() and save_pub()).
## Run from the repository root, or set HPICDA_ROOT first (see README).
source("R/hpicda_theme.R")

d <- load_hpicda()
calls <- d$calls

CLASSES <- c("CDEP_PATHOGEN_EXPLOITED", "CDEP_HOST_PROTECTIVE",
             "BASELINE_FITNESS_CONFOUNDED", "SIGNAL_STATE_ALTERING")

obs <- calls[!is.na(baseline_chronos) & !is.na(lfc)]
set.seed(1)
bg <- obs[call == "NO_FUNCTIONAL_SIGNAL"][sample(.N, min(.N, 60000))]
sig <- obs[call %in% CLASSES]
sig[, cls := factor(call, levels = CLASSES)]
cat(sprintf("panel a: background %s sampled, %s signal records\n",
            format(nrow(bg), big.mark = ","), format(nrow(sig), big.mark = ",")))

## ------------------------------------------------------------- (a) landscape
pa <- ggplot() +
  geom_point(data = bg, aes(baseline_chronos, lfc), shape = 16,
             colour = HP$neutral, size = 0.12, alpha = 0.12) +
  geom_point(data = sig, aes(baseline_chronos, lfc, colour = cls), shape = 16,
             size = 0.45, alpha = 0.62) +
  geom_vline(xintercept = ESS, linetype = "dashed", linewidth = 0.4, colour = HP$muted) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = HP$muted) +
  scale_colour_manual(values = HP$call_cols[CLASSES], labels = HP$call_labs[CLASSES],
                      name = NULL, drop = FALSE) +
  annotate("text", x = -4.35, y = 8.4, label = "baseline-essential",
           size = PT(DENSE_PT), colour = HP$muted, fontface = "italic", hjust = 0) +
  annotate("text", x = 1.2, y = 8.4, label = "baseline non-essential",
           size = PT(DENSE_PT), colour = HP$muted, fontface = "italic", hjust = 1) +
  coord_cartesian(xlim = c(-4.5, 1.3), ylim = c(-11.5, 9.2)) +
  labs(x = "Baseline fitness (DepMap Chronos; more negative = more essential)",
       y = "Infection effect (log2 fold-change)",
       title = "Conditional calls occupy non-essential space") +
  theme_hp() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.03),
        legend.justification = c(0, 0), legend.key.size = unit(3.2, "mm"),
        legend.spacing.y = unit(0.4, "mm"),
        legend.text = element_text(size = 6.6))

## --------------------------------------------- (b) baseline-essential fraction
## restrict to the five defined call classes: UNAVAILABLE / TOXIN_SPECIFIC have no
## plotted level and must not appear as an "NA" bar.
ALLC <- c(CLASSES, "NO_FUNCTIONAL_SIGNAL")
frac <- obs[call %in% ALLC, .(pct = mean(baseline_chronos < ESS) * 100), by = call]
frac <- frac[order(match(call, ALLC))]
frac[, cls := factor(call, levels = rev(ALLC))]
YLAB <- c("Pathogen-exploited\nconditional", "Host-protective\nconditional",
          "Baseline-fitness\nconfounded", "Signal-state\naltering",
          "No functional\nsignal")
names(YLAB) <- ALLC
pb <- ggplot(frac, aes(pct, cls, fill = cls)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), hjust = -0.16, size = PT(ANNOT_PT),
            fontface = "bold", colour = HP$ink) +
  geom_vline(xintercept = 50, linetype = "dotted", linewidth = 0.4, colour = HP$rule) +
  scale_fill_manual(values = HP$call_cols, guide = "none") +
  scale_y_discrete(labels = YLAB) +
  scale_x_continuous(limits = c(0, 122), expand = c(0, 0)) +
  labs(x = "Calls that are baseline-essential (%)", y = NULL,
       title = "Conditional calls are rarely essential") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 6.8, lineheight = 1.15),
        plot.margin = margin(6, 4, 6, 3))

## ------------------------------------------------------------------- compose
## Two panels only. A third panel (minimum baseline Chronos against HG3 total
## score for the 353 frozen candidates) used to sit here, but it plotted exactly
## the same 343 gene/value pairs as Fig. 6a. The therapeutic-window reading of
## that scatter belongs with the tractability section, so it is carried by
## Fig. 6a alone and this figure stops at the landscape (a) and its class
## breakdown (b).
fig <- (pa | pb) +
  plot_layout(widths = c(1.0, 0.78)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold"))

save_pub(fig, "fig2_baseline_vs_infection_landscape", 183, 92)
