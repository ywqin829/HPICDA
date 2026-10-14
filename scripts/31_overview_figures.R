# HPICDA overview figures (funnel + Bio-Tier 2 evidence matrix)
# User-facing summary figures; data from verified Stage 5/6/12 artifacts.
suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

OUT <- "/workspace/HPICDA/15_figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

font <- "Liberation Sans"
theme_set(theme_bw(base_size = 12, base_family = font) +
            theme(panel.grid.minor = element_blank()))

## ---------- Figure 1: candidate funnel ----------
funnel <- data.frame(
  stage = factor(c("逐屏条件依赖记录", "功能信号记录\n(EXP+PRO+STATE)",
                   "预注册候选基因\n(HG3 前冻结)", "HG3 评分达标\n(Tier1+Tier2)",
                   "Bio-Tier 2\n主要候选"),
                 levels = c("逐屏条件依赖记录", "功能信号记录\n(EXP+PRO+STATE)",
                            "预注册候选基因\n(HG3 前冻结)", "HG3 评分达标\n(Tier1+Tier2)",
                            "Bio-Tier 2\n主要候选")),
  n = c(852216, 7765, 353, 270, 18)
)

funnel$stage <- factor(funnel$stage, levels = rev(levels(funnel$stage)))
p1 <- ggplot(funnel, aes(x = stage, y = n)) +
  geom_col(fill = c("#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"),
           width = 0.68) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.12,
            family = font, size = 3.6) +
  scale_y_log10(expand = expansion(mult = c(0, 0.22))) +
  coord_flip() +
  labs(x = NULL, y = "数量(log10 刻度)") +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(OUT, "fig_funnel.svg"), p1, width = 7.2, height = 3.6)
ggsave(file.path(OUT, "fig_funnel.png"), p1, width = 7.2, height = 3.6, dpi = 200)

## ---------- Figure 2: Bio-Tier 2 evidence matrix ----------
fr <- read.delim("/workspace/HPICDA/09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN.tsv",
                 stringsAsFactors = FALSE)
bt2 <- fr[fr$bio_tier == "Bio-Tier 2", ]
bt2 <- bt2[order(-bt2$total_score, bt2$gene), ]
bt2$gene_dir <- ifelse(bt2$intervention_direction == "INHIBIT",
                       paste0(bt2$gene, " (抑制)"), paste0(bt2$gene, " (保护/激活)"))

# evidence columns -> long format, per-column 0-1 scaled fill + raw text label
mk <- function(value, disp, lo, hi) {
  data.frame(value = value, disp = disp,
             scaled = pmax(0, pmin(1, (value - lo) / (hi - lo))),
             stringsAsFactors = FALSE)
}
rows <- list()
add <- function(rows, gene_dir, col, value, disp, lo, hi) {
  d <- mk(value, disp, lo, hi)
  d$gene_dir <- gene_dir; d$evidence <- col
  rows[[length(rows) + 1]] <- d
  rows
}
for (i in seq_len(nrow(bt2))) {
  g <- bt2$gene_dir[i]
  rows <- add(rows, g, "屏幕证据 (0–6)", bt2$screen_evidence_pts[i],
              sprintf("%.1f", bt2$screen_evidence_pts[i]), 0, 6)
  rows <- add(rows, g, "跨病原体收敛 (0–2)", bt2$convergence_bonus[i],
              sprintf("%.0f", bt2$convergence_bonus[i]), 0, 2)
  loso <- bt2$loso_stability[i]
  rows <- add(rows, g, "LOSO 稳定性",
              ifelse(is.na(loso), NA_real_, loso),
              ifelse(is.na(loso), "—", sprintf("%.2f", loso)), 0, 1)
  chr <- bt2$min_baseline_chronos[i]
  rows <- add(rows, g, "基线 Chronos\n(越负越必需)",
              ifelse(is.na(chr), NA_real_, chr),
              ifelse(is.na(chr), "—", sprintf("%.2f", chr)), -2.5, 0.5)
  rows <- add(rows, g, "肺表达", as.numeric(bt2$lung_expressed[i] == "True"),
              ifelse(bt2$lung_expressed[i] == "True", "是", "否"), 0, 1)
  rows <- add(rows, g, "临床 DE 一致", as.numeric(bt2$clinical_de_concordant[i] == "True"),
              ifelse(bt2$clinical_de_concordant[i] == "True", "是", "否"), 0, 1)
  gt <- bt2$genetics_tier[i]
  gv <- c(G0 = 0, G3 = 1 / 3, G2 = 2 / 3, G1 = 1)[gt]
  rows <- add(rows, g, "人类遗传层", as.numeric(gv), gt, 0, 1)
  np <- length(strsplit(bt2$pathogens[i], ";")[[1]])
  short_path <- gsub("Influenza", "IAV", gsub("SARS-CoV-2", "SC2", bt2$pathogens[i]))
  rows <- add(rows, g, "病原体支持", np, short_path, 1, 2)
}
long <- do.call(rbind, rows)
long$gene_dir <- factor(long$gene_dir, levels = rev(bt2$gene_dir))
long$evidence <- factor(long$evidence, levels = c(
  "屏幕证据 (0–6)", "跨病原体收敛 (0–2)", "LOSO 稳定性", "基线 Chronos\n(越负越必需)",
  "肺表达", "临床 DE 一致", "人类遗传层", "病原体支持"))
long$txt_col <- ifelse(!is.na(long$scaled) & long$scaled > 0.55, "white", "black")

p2 <- ggplot(long, aes(x = evidence, y = gene_dir)) +
  geom_tile(aes(fill = scaled), color = "white", linewidth = 0.6) +
  geom_text(aes(label = disp, color = txt_col), family = font, size = 2.9) +
  scale_color_identity() +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", na.value = "#E8E8E8",
                      guide = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 9.5),
        axis.text.y = element_text(size = 9.5),
        panel.grid = element_blank(),
        plot.margin = margin(6, 10, 6, 6))

ggsave(file.path(OUT, "fig_biotier2_evidence.svg"), p2, width = 8.6, height = 6.2)
ggsave(file.path(OUT, "fig_biotier2_evidence.png"), p2, width = 8.6, height = 6.2, dpi = 200)

cat("figures written to", OUT, "\n")
print(bt2[, c("gene", "intervention_direction", "total_score")])
