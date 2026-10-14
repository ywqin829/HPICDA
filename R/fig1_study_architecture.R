# Figure 1 — Study architecture, candidate funnel and screen census
#
# Panel a is RESERVED for a hand-drawn study schematic; see
# 25_manuscript_package/figure_prompts/fig1a_study_schematic_BioRender_prompt.md
# The two data panels below are tagged b and c so the schematic can be dropped in
# as panel a without renumbering anything.
## Shared figure theme (defines ROOT, the warm palette, PT() and save_pub()).
## Run from the repository root, or set HPICDA_ROOT first (see README).
source("R/hpicda_theme.R")

d <- load_hpicda()
calls <- d$calls; frozen <- d$frozen

n_records <- nrow(calls)
n_signal  <- calls[call %in% c("CDEP_PATHOGEN_EXPLOITED", "CDEP_HOST_PROTECTIVE",
                               "SIGNAL_STATE_ALTERING"), .N]
n_cand <- nrow(frozen)
n_t12  <- frozen[bio_tier %in% c("Bio-Tier 2", "Bio-Tier 3"), .N]
n_t2   <- frozen[bio_tier == "Bio-Tier 2", .N]
cat(sprintf("funnel: %s -> %s -> %d -> %d -> %d\n",
            format(n_records, big.mark = ","), format(n_signal, big.mark = ","),
            n_cand, n_t12, n_t2))

## ---------------------------------------------------------------- (b) funnel
STAGES <- c("Gene x screen\ncall records", "Functional\nsignal records",
            "Pre-registered\ncandidates", "HG3-scored\n(Tier 1-2)",
            "Bio-Tier 2\nprimary candidates")
fun <- data.frame(
  stage = factor(STAGES, levels = rev(STAGES)),
  n = c(n_records, n_signal, n_cand, n_t12, n_t2))
# log-scale bars must not be anchored at 0 (dropped by the scale); anchor at the
# axis floor instead, so the bar reads as a magnitude band rather than 0 -> n.
# y stays the stage FACTOR (not as.numeric) so the discrete scale supplies the
# stage names; the rect offsets are computed from the factor code inside aes().
FLOOR <- 5
pb <- ggplot(fun, aes(y = stage)) +
  geom_rect(aes(xmin = FLOOR, xmax = n, ymin = as.numeric(stage) - 0.32,
                ymax = as.numeric(stage) + 0.32, fill = stage)) +
  geom_text(aes(x = n, label = format(n, big.mark = ",")), hjust = -0.14,
            size = PT(ANNOT_PT), fontface = "bold", colour = HP$ink) +
  scale_fill_manual(values = rev(HP$seq_warm[c(3, 4, 5, 6, 8)]), guide = "none") +
  # plain decade labels: plotmath 10^x renders the exponent at 0.7x the tick
  # size, which drops below the 6 pt floor
  scale_x_log10(limits = c(5, 9e6), breaks = 10^(0:6),
                labels = c("1", "10", "100", "1k", "10k", "100k", "1M")) +
  scale_y_discrete(expand = expansion(add = 0.12)) +
  labs(x = "Records / genes (log scale)", y = NULL,
       title = "Candidate funnel") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.0, lineheight = 1.15),
        plot.margin = margin(6, 4, 6, 3))

## ----------------------------------------------------------- (c) screen census
census <- calls[, .(n = uniqueN(screen_id)), by = .(pathogen, track)]
tot <- census[, .(tot = sum(n)), by = pathogen][order(tot)]
census[, pathogen := factor(pathogen, levels = tot$pathogen)]
census[, track := factor(track, levels = c("A", "B"))]
pc <- ggplot(census, aes(n, pathogen, fill = track)) +
  geom_col(width = 0.62) +
  geom_text(data = tot, aes(x = tot + 0.7, y = pathogen, label = tot),
            inherit.aes = FALSE, hjust = 0, size = PT(ANNOT_PT),
            fontface = "bold", colour = HP$ink) +
  scale_fill_manual(values = c(A = HP$exploited, B = HP$protective),
                    labels = c("Track A (raw counts)", "Track B (gene-level)"),
                    name = NULL) +
  scale_x_continuous(limits = c(0, max(tot$tot) * 1.24)) +
  labs(x = "Screens contributing to conditional calls", y = NULL,
       title = "Screen census") +
  theme_hp() +
  theme(axis.text.y = element_text(size = 7.4),
        legend.position = "inside", legend.position.inside = c(0.98, 0.16),
        legend.justification = c(1, 0),
        legend.key.size = unit(3.2, "mm"), legend.spacing.y = unit(0.5, "mm"),
        plot.margin = margin(6, 4, 6, 3))

## ------------------------------------------------------------------- compose
fig <- (pb | pc) +
  plot_layout(widths = c(1.0, 0.95)) +
  plot_annotation(tag_levels = list(c("B", "C"))) &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold"))

save_pub(fig, "fig1_study_architecture", 183, 72)
