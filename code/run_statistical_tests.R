#!/usr/bin/env Rscript

# =============================================================================
# Statistical Tests: LR Simulation Analysis
# =============================================================================
# Runs all inferential tests across the three analysis datasets and writes
# results to CSV tables for use in the manuscript results section.
#
# SECTION 1 — Matched population (pair-level)
#   Input: combined_LR_match.csv.gz
#   Q1: Does log10(LR) differ by relationship type? (KW + epsilon-squared, per loci set)
#   Q2: Does log10(LR) increase with more loci?     (Spearman, per relationship)
#   Q3: Does population affect log10(LR)?           (KW + epsilon-squared, at autosomal 29)
#
# SECTION 2 — Mismatched population (pair-level)
#   Input: combined_LR_all.rds (filtered to true positives, rel_match == TRUE)
#   Q4: Does using wrong population frequencies inflate log10(LR)?
#       (paired Wilcoxon signed-rank on log10(LR_wrong / LR_correct), per
#        population x tested_population x relationship x loci_set)
#   Q5: Which factors drive LR inflation due to population mismatch?
#       (linear model: log10_ratio ~ true_population + tested_population + loci_set,
#        per relationship type)
#
# SECTION 3 — FPR cutoff analysis (aggregated)
#   Input: proportions_with_classification.csv
#   Q6: Does FPR increase monotonically with loci count?
#       (Cochran-Armitage trend test, per tested_relationship x population x classification)
#   Q7: Do relationship type, population, and loci set affect FPR?
#       (logistic regression: cbind(successes, failures) ~ known_relationship +
#        population + loci_set, per tested_relationship x threshold)
#
# NOTE: With large simulation sample sizes p-values will be near-zero for most
# tests. Effect sizes are the informative statistic throughout.
#
# Usage:
#   Rscript code/run_statistical_tests.R <input_dir> [output_dir]
#
#   input_dir   Full path to analysis directory (e.g., output/lr_analysis_20260410)
#   output_dir  Where to write results CSVs (default: <input_dir>/stats)
# =============================================================================

suppressMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(rstatix)     # kruskal_effsize (epsilon-squared)
}))

log_message <- function(msg) cat(paste0("[", Sys.time(), "] ", msg, "\n"))

# --- Argument parsing ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript code/run_statistical_tests.R <input_dir> [output_dir]")
}

input_dir  <- args[1]
output_dir <- if (length(args) >= 2) args[2] else file.path(input_dir, "stats")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_message(paste("Input directory: ", input_dir))
log_message(paste("Output directory:", output_dir))


# =============================================================================
# SHARED CONSTANTS
# =============================================================================

relationship_order  <- c("parent_child", "full_siblings", "half_siblings",
                         "cousins", "second_cousins", "unrelated")
relationship_labels <- c(
  "parent_child"   = "Parent-Child",
  "full_siblings"  = "Full Siblings",
  "half_siblings"  = "Half Siblings",
  "cousins"        = "Cousins",
  "second_cousins" = "Second Cousins",
  "unrelated"      = "Unrelated"
)

loci_set_order  <- c("core_13", "identifiler_15", "expanded_20",
                     "supplementary", "autosomal_29")
loci_set_labels <- c(
  "core_13"        = "Core 13",
  "identifiler_15" = "Identifiler 15",
  "expanded_20"    = "Expanded 20",
  "supplementary"  = "Supplementary",
  "autosomal_29"   = "Autosomal 29"
)

loci_counts <- c(
  "core_13" = 13, "identifiler_15" = 15, "expanded_20" = 20,
  "supplementary" = 23, "autosomal_29" = 29
)

population_order <- c("AfAm", "Asian", "Cauc", "Hispanic", "all")


# =============================================================================
# SECTION 1: MATCHED POPULATION TESTS
# Input: combined_LR_match.csv.gz (pair-level, strict match)
# =============================================================================

log_message("=== SECTION 1: Matched population tests ===")

match_file <- file.path(input_dir, "combined_LR_match.csv.gz")
if (!file.exists(match_file)) {
  stop(sprintf("File not found: %s", match_file))
}

log_message("Loading combined_LR_match.csv.gz...")
matched <- read_csv(match_file, show_col_types = FALSE) %>%
  filter(is_correct_pop == TRUE, known_relationship == tested_relationship) %>%
  mutate(
    combined_LR  = as.numeric(combined_LR),
    log_LR       = log10(combined_LR),
    relationship = factor(known_relationship,
                          levels = relationship_order,
                          labels = relationship_labels),
    loci_set     = factor(loci_set,
                          levels = loci_set_order,
                          labels = loci_set_labels),
    population   = factor(population, levels = population_order)
  )

log_message(sprintf("Loaded %s matched observations.", format(nrow(matched), big.mark = ",")))

# Data objects matching the objects used in plots_matched_publication.R
main_data    <- matched %>% filter(population == "all")
pop_data_29  <- matched %>% filter(loci_set == "Autosomal 29", population != "all")

# ------------------------------------------------------------------------------
# TEST 1: Relationship type effect per loci set
# Kruskal-Wallis + epsilon-squared ("all" population only)
# ------------------------------------------------------------------------------
log_message("Test 1: KW relationship type effect per loci set...")

kw_rel_test <- main_data %>%
  group_by(loci_set) %>%
  summarise(
    kw_statistic = kruskal.test(log_LR ~ relationship)$statistic,
    kw_df        = kruskal.test(log_LR ~ relationship)$parameter,
    kw_p_value   = kruskal.test(log_LR ~ relationship)$p.value,
    n_pairs      = n(),
    .groups = "drop"
  )

kw_rel_effsize <- main_data %>%
  group_by(loci_set) %>%
  group_modify(~ {
    es <- kruskal_effsize(.x, formula = log_LR ~ relationship, ci = FALSE)
    tibble(epsilon_sq = es$effsize, magnitude = as.character(es$magnitude))
  }) %>%
  ungroup()

test1 <- left_join(kw_rel_test, kw_rel_effsize, by = "loci_set") %>%
  mutate(
    section  = "Matched",
    test     = "Kruskal-Wallis",
    question = "Does log10(LR) differ by relationship type?",
    grouping = as.character(loci_set),
    note     = "All populations combined (population = 'all')"
  )

# ------------------------------------------------------------------------------
# TEST 2: Loci count effect — Spearman correlation per relationship
# ------------------------------------------------------------------------------
log_message("Test 2: Spearman loci count effect...")

median_by_loci <- main_data %>%
  group_by(relationship, loci_set) %>%
  summarise(median_logLR = median(log_LR), .groups = "drop") %>%
  mutate(loci_count = loci_counts[as.character(loci_set)])

test2 <- median_by_loci %>%
  group_by(relationship) %>%
  summarise(
    spearman_rho = tryCatch(
      cor(loci_count, median_logLR, method = "spearman"),
      warning = function(w) NA_real_,
      error   = function(e) NA_real_
    ),
    spearman_p = tryCatch(
      cor.test(loci_count, median_logLR,
               method = "spearman", exact = FALSE)$p.value,
      warning = function(w) NA_real_,
      error   = function(e) NA_real_
    ),
    n_loci_sets  = n(),
    .groups = "drop"
  ) %>%
  mutate(
    section  = "Matched",
    test     = "Spearman correlation",
    question = "Does log10(LR) increase with more loci?",
    grouping = as.character(relationship),
    note     = "rho = Spearman rank correlation; NA indicates zero variance (e.g. Unrelated)"
  )

# ------------------------------------------------------------------------------
# TEST 3: Population effect at Autosomal 29
# Kruskal-Wallis + epsilon-squared per relationship
# ------------------------------------------------------------------------------
log_message("Test 3: KW population effect at Autosomal 29...")

kw_pop_test <- pop_data_29 %>%
  group_by(relationship) %>%
  summarise(
    kw_statistic = kruskal.test(log_LR ~ population)$statistic,
    kw_df        = kruskal.test(log_LR ~ population)$parameter,
    kw_p_value   = kruskal.test(log_LR ~ population)$p.value,
    n_pairs      = n(),
    .groups = "drop"
  )

kw_pop_effsize <- pop_data_29 %>%
  group_by(relationship) %>%
  group_modify(~ {
    es <- kruskal_effsize(.x, formula = log_LR ~ population, ci = FALSE)
    tibble(epsilon_sq = es$effsize, magnitude = as.character(es$magnitude))
  }) %>%
  ungroup()

test3 <- left_join(kw_pop_test, kw_pop_effsize, by = "relationship") %>%
  mutate(
    section  = "Matched",
    test     = "Kruskal-Wallis",
    question = "Does population affect log10(LR)?",
    grouping = as.character(relationship),
    note     = "Autosomal 29 only; named populations (excludes 'all')"
  )

# Save section 1
matched_stats <- bind_rows(
  test1 %>% transmute(section, test, question, grouping,
                      statistic = round(coalesce(kw_statistic), 3),
                      df = kw_df, p_value = kw_p_value,
                      effect_size = round(epsilon_sq, 3),
                      effect_size_type = "epsilon-squared",
                      magnitude, n = n_pairs, note),
  test2 %>% transmute(section, test, question, grouping,
                      statistic = round(spearman_rho, 3),
                      df = NA_real_, p_value = spearman_p,
                      effect_size = round(spearman_rho, 3),
                      effect_size_type = "Spearman rho",
                      magnitude = NA_character_, n = n_loci_sets, note),
  test3 %>% transmute(section, test, question, grouping,
                      statistic = round(kw_statistic, 3),
                      df = kw_df, p_value = kw_p_value,
                      effect_size = round(epsilon_sq, 3),
                      effect_size_type = "epsilon-squared",
                      magnitude, n = n_pairs, note)
)

fwrite(matched_stats, file.path(output_dir, "stats_matched.csv"))
log_message(sprintf("Wrote stats_matched.csv (%d rows)", nrow(matched_stats)))

rm(matched, main_data, pop_data_29)
gc()


# =============================================================================
# SECTION 2: MISMATCHED POPULATION TESTS
# Input: combined_LR_all.rds (pair-level, filtered to true positives)
# =============================================================================

log_message("=== SECTION 2: Mismatched population tests ===")

rds_file <- file.path(input_dir, "combined_LR_all.rds")
if (!file.exists(rds_file)) {
  stop(sprintf("File not found: %s\nRun analyze_lr_outputs.R first.", rds_file))
}

log_message("Loading combined_LR_all.rds...")
all_combined <- readRDS(rds_file)
log_message(sprintf("Loaded %s rows.", format(nrow(all_combined), big.mark = ",")))

# Preprocess — keep only true positive pairs (rel_match == TRUE)
# so LR_wrong and LR_correct are directly comparable per pair
all_combined <- all_combined %>%
  mutate(
    combined_LR  = as.numeric(combined_LR),
    log10_LR     = suppressWarnings(log10(combined_LR)),
    rel_match    = known_relationship == tested_relationship,
    pop_match    = tested_population == population,
    known_relationship  = factor(known_relationship,  levels = relationship_order),
    tested_relationship = factor(tested_relationship, levels = relationship_order),
    loci_set            = factor(loci_set,            levels = loci_set_order),
    population          = factor(population,          levels = population_order)
  ) %>%
  filter(rel_match == TRUE,
         known_relationship %in% c("parent_child", "full_siblings"),
         !is.na(combined_LR))

log_message(sprintf("After filtering to true positives: %s rows.",
                    format(nrow(all_combined), big.mark = ",")))

# Build pair-level log10 ratio: join correct and wrong LRs by pair
correct_lrs <- all_combined %>%
  filter(pop_match == TRUE) %>%
  select(batch_id, pair_id, population, known_relationship,
         loci_set, tested_relationship,
         correct_log10_LR = log10_LR)

wrong_lrs <- all_combined %>%
  filter(pop_match == FALSE) %>%
  select(batch_id, pair_id, population, known_relationship,
         loci_set, tested_relationship, tested_population,
         wrong_log10_LR = log10_LR)

pairs_ratio <- wrong_lrs %>%
  left_join(correct_lrs,
            by = c("batch_id", "pair_id", "population",
                   "known_relationship", "loci_set", "tested_relationship")) %>%
  mutate(log10_ratio = wrong_log10_LR - correct_log10_LR) %>%
  filter(!is.na(log10_ratio), is.finite(log10_ratio))

log_message(sprintf("Pair-level ratio table: %s rows.",
                    format(nrow(pairs_ratio), big.mark = ",")))

# ------------------------------------------------------------------------------
# TEST 4: Does using wrong population frequencies inflate log10(LR)?
# Paired Wilcoxon signed-rank test (H0: median log10_ratio = 0) per
# true_population x tested_population x relationship x loci_set cell.
# Effect size: median log10_ratio (directly interpretable: +1 = 10x inflation)
# ------------------------------------------------------------------------------
log_message("Test 4: Paired Wilcoxon on log10(LR_wrong / LR_correct)...")

test4 <- pairs_ratio %>%
  group_by(population, tested_population, known_relationship, loci_set) %>%
  summarise(
    median_log10_ratio = median(log10_ratio, na.rm = TRUE),
    wilcox_p           = tryCatch(
      wilcox.test(log10_ratio, mu = 0, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    n_pairs            = n(),
    .groups = "drop"
  ) %>%
  mutate(
    section          = "Mismatched population",
    test             = "Wilcoxon signed-rank",
    question         = "Does using wrong population frequencies inflate log10(LR)?",
    effect_size      = round(median_log10_ratio, 3),
    effect_size_type = "median log10(LR_wrong / LR_correct)",
    note             = "+1 = 10x inflation; -1 = 10x deflation; 0 = no effect"
  )

# ------------------------------------------------------------------------------
# TEST 5: Which factors drive LR inflation?
# Linear model: log10_ratio ~ true_population + tested_population + loci_set
# per relationship type. Reports eta-squared (SS_effect / SS_total) as effect size.
# Note: fitted on group medians (one row per cell) to avoid pseudoreplication
# from the large pair-level n inflating F statistics.
# ------------------------------------------------------------------------------
log_message("Test 5: Linear model of log10_ratio by population and loci factors...")

# Summarise to one median per cell to avoid pseudoreplication
cell_medians <- pairs_ratio %>%
  group_by(population, tested_population, known_relationship, loci_set) %>%
  summarise(median_log10_ratio = median(log10_ratio, na.rm = TRUE),
            .groups = "drop")

fit_lm_section <- function(data, relationship) {
  df <- data %>% filter(known_relationship == relationship)
  if (nrow(df) < 10) return(NULL)

  fit  <- lm(median_log10_ratio ~ population + tested_population + loci_set, data = df)
  aov  <- anova(fit)
  ss_total <- sum(aov$`Sum Sq`)

  tibble(
    term             = rownames(aov),
    sum_sq           = aov$`Sum Sq`,
    df_term          = aov$Df,
    f_value          = aov$`F value`,
    p_value          = aov$`Pr(>F)`,
    eta_sq           = aov$`Sum Sq` / ss_total
  ) %>%
    filter(term != "Residuals") %>%
    mutate(known_relationship = relationship)
}

test5 <- bind_rows(lapply(
  c("parent_child", "full_siblings"),
  function(rel) fit_lm_section(cell_medians, rel)
)) %>%
  mutate(
    section          = "Mismatched population",
    test             = "Linear model (ANOVA)",
    question         = "Which factors drive LR inflation from population mismatch?",
    effect_size_type = "eta-squared",
    note             = "Fitted on per-cell medians to avoid pseudoreplication"
  )

# Save section 2
mismatch_stats <- bind_rows(
  test4 %>% transmute(section, test, question,
                      grouping = paste(population, tested_population,
                                       known_relationship, loci_set, sep = " | "),
                      statistic = NA_real_, df = NA_real_,
                      p_value = wilcox_p,
                      effect_size, effect_size_type,
                      magnitude = NA_character_, n = n_pairs, note),
  test5 %>% transmute(section, test, question,
                      grouping = paste(known_relationship, term, sep = " | "),
                      statistic = round(f_value, 3),
                      df = df_term, p_value,
                      effect_size = round(eta_sq, 3),
                      effect_size_type,
                      magnitude = NA_character_, n = NA_real_, note)
)

fwrite(mismatch_stats, file.path(output_dir, "stats_mismatched_population.csv"))
log_message(sprintf("Wrote stats_mismatched_population.csv (%d rows)", nrow(mismatch_stats)))

rm(all_combined, correct_lrs, wrong_lrs, pairs_ratio, cell_medians)
gc()


# =============================================================================
# SECTION 3: FPR CUTOFF TESTS
# Input: proportions_with_classification.csv (aggregated)
# =============================================================================

log_message("=== SECTION 3: FPR cutoff tests ===")

props_file <- file.path(input_dir, "proportions_with_classification.csv")
if (!file.exists(props_file)) {
  stop(sprintf(
    "File not found: %s\nRun prepare_combined_lr_intermediates.R first.", props_file
  ))
}

log_message("Loading proportions_with_classification.csv...")
proportions_all <- fread(props_file) %>%
  filter(tested_population == population) %>%   # matched frequencies only
  mutate(
    known_relationship  = factor(known_relationship,  levels = relationship_order),
    tested_relationship = factor(tested_relationship, levels = relationship_order),
    loci_set            = factor(loci_set,            levels = loci_set_order),
    population          = factor(population,          levels = population_order)
  )

log_message(sprintf("Loaded %d rows.", nrow(proportions_all)))

# ------------------------------------------------------------------------------
# TEST 6: Does FPR increase monotonically with loci count?
# Cochran-Armitage trend test per tested_relationship x population x
# classification x threshold combination.
# Tests the ordered hypothesis: more loci -> higher (or lower) FPR.
# ------------------------------------------------------------------------------
log_message("Test 6: Cochran-Armitage trend tests for loci count effect on FPR...")

threshold_cols <- c(
  "prop_LR_gt_1"    = "Fixed (LR > 1)",
  "prop_LR_gt_10"   = "1% FPR",
  "prop_LR_gt_100"  = "0.1% FPR",
  "prop_LR_gt_1000" = "0.01% FPR"
)

run_ca_trend <- function(data, threshold_col, threshold_label) {
  data %>%
    filter(!is.na(.data[[threshold_col]]),
           tested_relationship %in% c("parent_child", "full_siblings"),
           classification != "True Positive",
           population != "all") %>%
    group_by(tested_relationship, population, classification) %>%
    arrange(n_loci) %>%
    group_modify(~ {
      n_exceed <- round(.x[[threshold_col]] * .x$n_related)
      test_result <- tryCatch(
        prop.trend.test(n_exceed, .x$n_related, score = .x$n_loci),
        error = function(e) NULL
      )
      tibble(
        ca_statistic = if (is.null(test_result)) NA_real_ else unname(test_result$statistic),
        ca_p_value   = if (is.null(test_result)) NA_real_ else test_result$p.value,
        threshold    = threshold_label
      )
    }) %>%
    ungroup()
}

test6 <- bind_rows(mapply(
  run_ca_trend,
  threshold_col   = names(threshold_cols),
  threshold_label = unname(threshold_cols),
  MoreArgs = list(data = proportions_all),
  SIMPLIFY = FALSE
)) %>%
  mutate(
    section          = "FPR cutoff",
    test             = "Cochran-Armitage trend",
    question         = "Does FPR change monotonically with loci count?",
    effect_size_type = NA_character_,
    note             = "False positive rows only; matched population frequencies"
  )

# ------------------------------------------------------------------------------
# TEST 7: Do relationship type, population, and loci set affect FPR?
# Logistic regression: successes/failures ~ known_relationship + population +
# loci_set, per tested_relationship x threshold, using n_related as binomial
# denominator. Reports deviance explained (1 - deviance/null deviance) per term
# via drop1().
# ------------------------------------------------------------------------------
log_message("Test 7: Logistic regression of FPR by relationship, population, loci set...")

run_logistic <- function(data, tested_rel, threshold_col, threshold_label) {
  df <- data %>%
    filter(tested_relationship == tested_rel,
           tested_population == population,
           population != "all",
           classification != "True Positive",
           !is.na(.data[[threshold_col]])) %>%
    mutate(
      successes = round(.data[[threshold_col]] * n_related),
      failures  = n_related - successes
    ) %>%
    filter(successes >= 0, failures >= 0)

  if (nrow(df) < 5) return(NULL)

  fit <- tryCatch(
    glm(cbind(successes, failures) ~ known_relationship + population + loci_set,
        data = df, family = binomial),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  drop_table <- drop1(fit, test = "Chisq")

  tibble(
    term      = rownames(drop_table)[-1],   # drop <none> row
    deviance  = drop_table$Deviance[-1],
    df_term   = drop_table$Df[-1],
    p_value   = drop_table$`Pr(>Chi)`[-1],
    dev_explained = drop_table$Deviance[-1] / fit$null.deviance
  ) %>%
    mutate(
      tested_relationship = tested_rel,
      threshold           = threshold_label
    )
}

test7 <- bind_rows(
  lapply(c("parent_child", "full_siblings"), function(rel) {
    bind_rows(mapply(
      run_logistic,
      threshold_col   = names(threshold_cols),
      threshold_label = unname(threshold_cols),
      MoreArgs = list(data = proportions_all, tested_rel = rel),
      SIMPLIFY = FALSE
    ))
  })
) %>%
  mutate(
    section          = "FPR cutoff",
    test             = "Logistic regression (drop1 LRT)",
    question         = "Do relationship type, population, and loci set affect FPR?",
    effect_size_type = "deviance explained (term deviance / null deviance)",
    note             = "False positive rows only; matched frequencies; excludes 'all' population"
  )

# Save section 3
fpr_stats <- bind_rows(
  test6 %>% transmute(section, test, question,
                      grouping = paste(tested_relationship, population,
                                       classification, threshold, sep = " | "),
                      statistic = round(ca_statistic, 3),
                      df = NA_real_, p_value = ca_p_value,
                      effect_size = NA_real_,
                      effect_size_type, magnitude = NA_character_,
                      n = NA_real_, note),
  test7 %>% transmute(section, test, question,
                      grouping = paste(tested_relationship, threshold, term, sep = " | "),
                      statistic = round(deviance, 3),
                      df = df_term, p_value,
                      effect_size = round(dev_explained, 3),
                      effect_size_type, magnitude = NA_character_,
                      n = NA_real_, note)
)

fwrite(fpr_stats, file.path(output_dir, "stats_fpr_cutoffs.csv"))
log_message(sprintf("Wrote stats_fpr_cutoffs.csv (%d rows)", nrow(fpr_stats)))


# =============================================================================
# SUMMARY
# =============================================================================

cat("\nOUTPUT FILES:\n")
for (f in c("stats_matched.csv", "stats_mismatched_population.csv", "stats_fpr_cutoffs.csv")) {
  fpath <- file.path(output_dir, f)
  if (file.exists(fpath)) {
    cat(sprintf("  %s (%d rows)\n", f, nrow(fread(fpath))))
  }
}
cat(sprintf("\nAll files saved to: %s\n", output_dir))
log_message("All statistical tests complete.")
