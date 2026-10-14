# HPICDA publication figure theme (warm palette, Nature-style specifications)
# Source this file at the top of every figure script.
#
# CONTRACT (do not violate in panel-level theme() calls):
#   * Both axis lines are always drawn. Never use axis.line.* = element_blank().
#     Exception: heatmap/tile panels that have no meaningful axis (set
#     axis.ticks = element_blank() only, and keep axis.line).
#   * All text is >= 6.5 pt at final size; body text is 7-8 pt.
#     Never pass a bare number to geom_text(size=); use PT(<points>).

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggrepel)
  library(svglite)
  library(ragg)
  library(grid)
})

## Repository root: the directory that holds this checkout together with the
## data downloaded from Zenodo (see README). Override with the HPICDA_ROOT
## environment variable, or simply run from that directory.
ROOT <- Sys.getenv("HPICDA_ROOT", unset = normalizePath("."))
OUT  <- file.path(ROOT, "25_manuscript_package", "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ggplot2 geom_text(size=) is in mm; font size in points = size * 2.8453.
## PT() converts a point size into the geom_text size argument so that panel
## code reads in real typographic units.
.PT_PER_MM <- 72.27 / 25.4
PT <- function(points) points / .PT_PER_MM

## ------------------------------------------------------------------ typography
BASE_PT   <- 8.0   # axis text / titles / legends  (skill spec: 7-9 pt)
TITLE_PT  <- 9.0   # panel titles
TAG_PT    <- 9.5   # panel tags A/B/C
ANNOT_PT  <- 7.0   # in-panel annotations and value labels
DENSE_PT  <- 6.6   # dense tick labels (heatmap columns, long categories)
TICK_PT   <- 7.2   # axis tick labels

## ---------------------------------------------------------------- warm palette
## Restrained warm family: one signal family (reds/oranges), one neutral family
## (tans/beiges), one accent family (olive/gold). No cold blue/grey/black.
HP <- list(
  # conditional-dependency call classes
  exploited   = "#C0392B",  # pathogen-exploited conditional   (signal)
  protective  = "#E98B2A",  # host-protective conditional      (signal)
  confounded  = "#6E4B21",  # baseline-fitness confounded      (neutral-dark)
  state       = "#A8B545",  # signal-state altering            (accent)
  neutral     = "#E6DCCB",  # no functional signal             (neutral-light)
  # evidence tiers
  tier2 = "#A93226", tier3 = "#E8912A", archive = "#B08850", reject = "#DCCFB8",
  # ink / rules
  ink = "#2B2118", muted = "#8A7B69", rule = "#CBBFAE", grid = "#EFE7DA",
  # warm sequential ramp (heatmaps)
  seq_warm = c("#FFFFE5", "#FFF7BC", "#FEE391", "#FEC44F",
               "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506"),
  # warm diverging-free two-tone for module membership strips
  strip_on = "#B5651D", strip_off = "#F3EBDD"
)
HP$call_cols <- c(CDEP_PATHOGEN_EXPLOITED = HP$exploited,
                  CDEP_HOST_PROTECTIVE    = HP$protective,
                  BASELINE_FITNESS_CONFOUNDED = HP$confounded,
                  SIGNAL_STATE_ALTERING   = HP$state,
                  NO_FUNCTIONAL_SIGNAL    = HP$neutral)
HP$call_labs <- c(CDEP_PATHOGEN_EXPLOITED = "Pathogen-exploited conditional",
                  CDEP_HOST_PROTECTIVE    = "Host-protective conditional",
                  BASELINE_FITNESS_CONFOUNDED = "Baseline-fitness confounded",
                  SIGNAL_STATE_ALTERING   = "Signal-state altering",
                  NO_FUNCTIONAL_SIGNAL    = "No functional signal")
HP$tier_cols <- c(`Bio-Tier 2` = HP$tier2, `Bio-Tier 3` = HP$tier3,
                  ARCHIVE = HP$archive, REJECT = HP$reject)
HP$tier_levels <- c("Bio-Tier 2", "Bio-Tier 3", "ARCHIVE", "REJECT")
ESS <- -0.5  # Chronos threshold for baseline essentiality

## --------------------------------------------------------------------- theme
theme_hp <- function(base_size = BASE_PT, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # both axes always visible -- this is a hard requirement of the figure set.
      # Kept light (0.28 pt): at 183 mm print width a heavier rule reads as a
      # frame rather than an axis.
      axis.line   = element_line(linewidth = 0.28, colour = HP$ink),
      axis.ticks  = element_line(linewidth = 0.28, colour = HP$ink),
      axis.text   = element_text(colour = HP$ink, size = TICK_PT),
      axis.title  = element_text(colour = HP$ink, size = base_size),
      legend.title = element_text(size = base_size),
      legend.text  = element_text(size = base_size - 0.8),
      legend.key.size = unit(3.4, "mm"),
      legend.background = element_blank(),
      strip.text  = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      plot.title  = element_text(size = TITLE_PT, face = "bold",
                                 colour = HP$ink, hjust = 0,
                                 margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 0.5, colour = HP$muted,
                                   hjust = 0, margin = margin(b = 2.5)),
      panel.grid  = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.tag = element_text(size = TAG_PT, face = "bold", family = base_family,
                              colour = HP$ink)
    )
}
theme_set(theme_hp())

## -------------------------------------------------------------- export helper
## PDF only, by request. Vector, Arial, editable text.
save_pub <- function(plot, stem, width_mm = 183, height_mm = 120, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(file.path(OUT, paste0(stem, ".pdf")), plot, width = w, height = h,
         device = cairo_pdf, family = "Arial")
  for (ext in c("svg", "png")) {
    f <- file.path(OUT, paste0(stem, ".", ext))
    if (file.exists(f)) unlink(f)
  }
  cat(sprintf("  wrote %s.pdf  (%.0f x %.0f mm)\n", stem, width_mm, height_mm))
}

## ------------------------------------------------------------------ load data
load_hpicda <- function() {
  list(
    calls  = fread(file.path(ROOT, "05_conditional_dependency/per_screen_calls_v2.tsv"),
                   showProgress = FALSE),
    frozen = fread(file.path(ROOT, "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv"),
                   showProgress = FALSE),
    annot  = fread(file.path(ROOT, "09_candidate_freeze/BIOLOGICAL_CANDIDATES_FROZEN_v2_annotated.tsv"),
                   showProgress = FALSE),
    conv   = fread(file.path(ROOT, "07_cross_pathogen/convergence_meta_v2_gene.tsv"),
                   showProgress = FALSE),
    modules = fread(file.path(ROOT, "07_cross_pathogen/convergence_meta_v2_module.tsv"),
                    showProgress = FALSE)
  )
}
