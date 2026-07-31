#!/usr/bin/env Rscript
################################################################################
# Focal Ranking Test - Cross-Database-Size Comparison
#
# Purpose: Compare the focal ranking test ACROSS unrelated-database sizes
#          (e.g. 10k vs 20k vs 100k). The per-size figures live in
#          plot_ranking_summary.R; this script overlays the sizes so the
#          effect of a larger haystack is visible in one figure.
#
#   Figure A (flagship): Recovery curves - the empirical CDF of the true
#            relative's rank. x = rank threshold N (log10), y = proportion of
#            replicates with rank <= N. Database size = colored lines, tested
#            hypothesis = line type, faceted true_relationship (rows) x
#            loci_set (cols). Dotted verticals mark top 10/50/100/200.
#            "Proportion in top N" is exactly this curve evaluated at N, so the
#            four bars in the per-size figures are four points on these lines.
#
#   Figure B: Threshold degradation - recovery at the four operational cutoffs
#            (top 10/50/100/200) as a function of database size, with Wilson
#            95% binomial CIs. Faceted true_relationship (rows) x
#            loci_set + tested_relationship (cols).
#
#   Tables:  recovery_by_threshold_across_databases.csv (long, with CIs)
#            rank_summary_across_databases.csv (median/mean/min/max rank)
#
# Input:  two or more combined ranking_outcomes CSVs, one per database size
#         (the combine_ranking_outcomes.sh output for each focal_test_* run).
#         Database size is read from the n_database column, not the filename.
#
# Usage:
#   Rscript code/compare_ranking_across_databases.R --out <dir> \
#       output/focal_test_10k/combined_all_N10000_outcomes.csv \
#       output/focal_test_20k/combined_all_N20k_outcomes.csv \
#       output/focal_test_100k/combined_all_N100000_outcomes.csv
#
# Date: 2026-07-31
################################################################################

# ==============================================================================
# SETUP AND DEPENDENCIES
# ==============================================================================

suppressMessages({
  library(tidyverse)
  library(data.table)
  library(scales)
})

log_message <- function(message) {
  cat(sprintf("[%s] %s\n", Sys.time(), message))
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

# --out <dir>  (or --out=<dir>); everything else is treated as an input file.
output_dir <- NULL
input_files <- character(0)
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a == "--out") {
    output_dir <- args[i + 1]; i <- i + 2
  } else if (grepl("^--out=", a)) {
    output_dir <- sub("^--out=", "", a); i <- i + 1
  } else {
    input_files <- c(input_files, a); i <- i + 1
  }
}

if (length(input_files) < 2) {
  stop("Provide at least two outcomes CSVs (one per database size).\n",
       "Usage: Rscript compare_ranking_across_databases.R --out <dir> file1.csv file2.csv [...]")
}
if (is.null(output_dir)) {
  output_dir <- file.path(dirname(input_files[1]), "cross_database_figures")
}

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Input file(s) not found:\n  ", paste(missing_files, collapse = "\n  "))
}

log_message(sprintf("Input files (%d):", length(input_files)))
for (f in input_files) log_message(sprintf("  %s", f))
log_message(sprintf("Output directory: %s", output_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# CONSTANTS AND STYLING
# ==============================================================================

RELATIONSHIP_ORDER <- c("parent_child", "full_siblings", "half_siblings",
                        "cousins", "second_cousins", "unrelated")
RELATIONSHIP_LABELS <- c("Parent-Child", "Full Siblings", "Half Siblings",
                         "Cousins", "Second Cousins", "Unrelated")

LOCI_ORDER <- c("core_13", "identifiler_15", "expanded_20",
                "supplementary", "autosomal_29")
LOCI_LABELS <- c("Core 13", "Identifiler 15", "Expanded 20",
                 "Supplementary", "Autosomal 29")

# Tested hypothesis -> line type (only PC / FS are tested in this design)
tested_linetypes <- c("Parent-Child" = "solid", "Full Siblings" = "22")

# Rank thresholds of interest
THRESHOLDS <- c(10, 50, 100, 200)

# Threshold colours (sequential amber; light = tight cutoff)
threshold_colors <- c("Top 10"  = "#EF9F27", "Top 50"  = "#BA7517",
                      "Top 100" = "#854F0B", "Top 200" = "#412402")

# Human-readable database-size label (10000 -> "10k", 100000 -> "100k")
db_size_label <- function(n) {
  ifelse(n >= 1e6, paste0(round(n / 1e6, 1), "M"),
         ifelse(n >= 1e3, paste0(round(n / 1e3), "k"), as.character(n)))
}

# Wilson score interval for a binomial proportion
wilson_ci <- function(k, n, z = 1.96) {
  p <- k / n
  denom  <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  tibble(lo = pmax(0, centre - half), hi = pmin(1, centre + half))
}

theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title    = element_text(size = rel(1.2), hjust = 0.5),
      plot.subtitle = element_text(size = rel(0.85), hjust = 0.5, color = "gray30"),
      axis.title    = element_text(size = rel(1)),
      axis.text     = element_text(size = rel(0.9)),
      strip.background = element_rect(fill = "gray90", color = "gray60"),
      strip.text       = element_text(size = rel(0.8)),
      legend.title    = element_text(size = rel(1)),
      legend.text     = element_text(size = rel(0.9)),
      legend.position = "bottom",
      legend.box      = "horizontal",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# DATA LOADING
# ==============================================================================

required_cols <- c("original_batch_id", "original_pair_id", "true_relationship",
                   "loci_set", "tested_relationship", "rank", "n_database",
                   "top_200", "top_100", "top_50", "top_10")

read_one <- function(f) {
  dt <- fread(f)
  miss <- setdiff(required_cols, names(dt))
  if (length(miss) > 0) {
    stop(sprintf("%s is missing columns: %s", basename(f), paste(miss, collapse = ", ")))
  }
  nd <- unique(dt$n_database)
  if (length(nd) > 1) {
    log_message(sprintf("WARNING: %s has multiple n_database values (%s) - each row is tagged by its own value.",
                        basename(f), paste(nd, collapse = ", ")))
  }
  dt$source_file <- basename(f)
  dt
}

log_message("Loading outcome files...")
all_raw <- rbindlist(lapply(input_files, read_one), use.names = TRUE, fill = TRUE)

all <- as_tibble(all_raw) %>%
  mutate(
    rank        = as.numeric(rank),
    n_database  = as.numeric(n_database),
    db_size_num = n_database - 1L,              # unrelated pool = candidates - true relative
    true_relationship = factor(true_relationship,
                               levels = RELATIONSHIP_ORDER, labels = RELATIONSHIP_LABELS),
    tested_relationship = factor(tested_relationship,
                                 levels = RELATIONSHIP_ORDER, labels = RELATIONSHIP_LABELS),
    loci_set = factor(loci_set, levels = LOCI_ORDER, labels = LOCI_LABELS)
  ) %>%
  filter(!is.na(rank))

# Ordered database-size factor + a sequential colour ramp (light -> dark)
db_num_levels <- sort(unique(all$db_size_num))
db_lab_levels <- db_size_label(db_num_levels)
all <- all %>% mutate(db_size = factor(db_size_label(db_size_num), levels = db_lab_levels))
db_colors <- setNames(
  colorRampPalette(c("#9FC3EA", "#04233F"))(length(db_lab_levels)),
  db_lab_levels
)

log_message(sprintf("Loaded %s rows across %d database size(s): %s",
                    format(nrow(all), big.mark = ","),
                    length(db_lab_levels), paste(db_lab_levels, collapse = ", ")))

# Sanity: unique focal replicates per (db, relationship), collision-proof key
rep_counts <- all %>%
  distinct(db_size, true_relationship, original_batch_id, original_pair_id) %>%
  count(db_size, true_relationship, name = "n_replicates")
cat("\nUnique focal replicates per database x true relationship:\n")
print(as.data.frame(rep_counts))

# ==============================================================================
# FIGURE A: RECOVERY CURVES (empirical CDF of rank)
# ==============================================================================

log_message("Generating Figure A: recovery curves across database sizes...")

max_n <- max(all$n_database, na.rm = TRUE)

# Log-spaced evaluation grid, with the operational thresholds forced in
n_grid <- sort(unique(round(c(
  1:9,
  10^seq(1, log10(max_n), length.out = 240),
  THRESHOLDS
))))
n_grid <- n_grid[n_grid >= 1 & n_grid <= max_n]

curve_data <- all %>%
  group_by(db_size, true_relationship, tested_relationship, loci_set) %>%
  group_modify(~{
    f <- ecdf(.x$rank)                       # P(rank <= N), handles averaged ties
    tibble(n_threshold = n_grid, prop = f(n_grid))
  }) %>%
  ungroup()

fig_curve <- ggplot(
  curve_data,
  aes(x = n_threshold, y = prop, color = db_size, linetype = tested_relationship)
) +
  geom_vline(xintercept = THRESHOLDS, color = "gray70",
             linetype = "dotted", linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  facet_grid(true_relationship ~ loci_set) +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000, 1e5, 1e6),
    labels = comma_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(values = db_colors, name = "Unrelated database size") +
  scale_linetype_manual(values = tested_linetypes, name = "Tested hypothesis") +
  labs(
    title    = "Focal Ranking Test: Recovery vs. Rank Threshold Across Database Sizes",
    subtitle = "Curve = proportion of true relatives with rank <= N | dotted verticals = top 10/50/100/200",
    x = "Rank threshold N (log scale)",
    y = "Proportion of true relatives with rank \u2264 N"
  ) +
  theme_publication() +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))

ggsave(file.path(output_dir, "recovery_curves_across_databases.pdf"),
       fig_curve, width = 11, height = 8.5, units = "in", bg = "white")
ggsave(file.path(output_dir, "recovery_curves_across_databases.png"),
       fig_curve, width = 11, height = 8.5, units = "in", dpi = 300, bg = "white")
log_message("Figure A saved: recovery_curves_across_databases.png/.pdf")

# ==============================================================================
# FIGURE B: RECOVERY AT EACH THRESHOLD VS DATABASE SIZE (with Wilson CIs)
# ==============================================================================

log_message("Generating Figure B: threshold degradation with database size...")

thresh_summary <- all %>%
  pivot_longer(c(top_10, top_50, top_100, top_200),
               names_to = "threshold", values_to = "in_top_n") %>%
  mutate(threshold = factor(threshold,
                            levels = c("top_10", "top_50", "top_100", "top_200"),
                            labels = c("Top 10", "Top 50", "Top 100", "Top 200"))) %>%
  group_by(db_size, db_size_num, true_relationship, tested_relationship,
           loci_set, threshold) %>%
  summarize(n = n(), k = sum(in_top_n), prop = mean(in_top_n), .groups = "drop") %>%
  bind_cols(wilson_ci(.$k, .$n))

fig_thresh <- ggplot(
  thresh_summary,
  aes(x = db_size, y = prop, color = threshold, group = threshold)
) +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.4) +
  geom_point(size = 1.8) +
  facet_grid(true_relationship ~ loci_set + tested_relationship) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(values = threshold_colors, name = "Rank threshold") +
  labs(
    title    = "Focal Ranking Test: Recovery at Each Threshold vs. Database Size",
    subtitle = "Points = proportion in top N | error bars = Wilson 95% CI | columns = loci set x tested hypothesis",
    x = "Unrelated database size",
    y = "Proportion of true relatives in top N"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "recovery_by_threshold_across_databases.pdf"),
       fig_thresh, width = 12, height = 8.5, units = "in", bg = "white")
ggsave(file.path(output_dir, "recovery_by_threshold_across_databases.png"),
       fig_thresh, width = 12, height = 8.5, units = "in", dpi = 300, bg = "white")
log_message("Figure B saved: recovery_by_threshold_across_databases.png/.pdf")

# ==============================================================================
# SUMMARY TABLES
# ==============================================================================

log_message("Writing summary tables...")

write_csv(
  thresh_summary %>%
    select(db_size, db_size_num, true_relationship, tested_relationship,
           loci_set, threshold, n, k, prop, ci_lo = lo, ci_hi = hi) %>%
    arrange(true_relationship, tested_relationship, loci_set, threshold, db_size_num),
  file.path(output_dir, "recovery_by_threshold_across_databases.csv")
)

rank_summary <- all %>%
  group_by(db_size, db_size_num, true_relationship, tested_relationship, loci_set) %>%
  summarize(
    n_replicates = n(),
    median_rank  = median(rank),
    mean_rank    = mean(rank),
    min_rank     = min(rank),
    max_rank     = max(rank),
    .groups = "drop"
  ) %>%
  arrange(true_relationship, tested_relationship, loci_set, db_size_num)

write_csv(rank_summary, file.path(output_dir, "rank_summary_across_databases.csv"))

cat("\nTop-200 recovery by database size (matched / most-informative reading):\n")
print(as.data.frame(
  thresh_summary %>%
    filter(threshold == "Top 200") %>%
    select(true_relationship, tested_relationship, loci_set, db_size, prop) %>%
    pivot_wider(names_from = db_size, values_from = prop) %>%
    arrange(true_relationship, tested_relationship, loci_set)
))

cat("\nOUTPUT FILES:\n")
cat("   - recovery_curves_across_databases.png/.pdf (Figure A)\n")
cat("   - recovery_by_threshold_across_databases.png/.pdf (Figure B)\n")
cat("   - recovery_by_threshold_across_databases.csv\n")
cat("   - rank_summary_across_databases.csv\n")
cat("\nAll files saved to:", output_dir, "\n")
log_message("Comparison complete!")
