#!/usr/bin/env Rscript

#### Plots of proportions exceeding cutoffs for tested parent-child and full-sibling relationships
#### Population Status: MATCHED ✓ (tested_population == population means you're using the correct allele frequencies)
#### Relationship Status: MISMATCHED ✗ (known_relationship != tested_relationship means you're testing the wrong relationship hypothesis)


# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)  # Includes ggplot2, dplyr, tidyr, etc.
  library(data.table) # Required for efficient fread()
  library(scales)     # For number formatting
  library(shades)     # For nice colors
}))

# Source required modules
source("code/module9_combinedLR_stats_functions.R")

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# --- Argument Parsing and Setup ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript code/plots_mismatched.R <input_dir> [output_dir]")
}

input_dir <- args[1]
log_message(paste("Input directory:", input_dir))

if (length(args) >= 2) {
  output_subdir <- args[2]
  output_dir <- file.path("output", output_subdir)
} else {
  timestamp <- format(Sys.time(), "%Y%m%d")
  output_dir <- file.path(input_dir, paste0("mismatched_pop_plots_", timestamp)) 
}
log_message(paste("Output directory:", output_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)



today <- format(Sys.time(), format = "%Y-%m-%d")


# Define factors, colors, and labels
populations <- c("AfAm", "Cauc", "Hispanic", "Asian", "all")

relationship_order <- c("parent_child", "full_siblings", "half_siblings", 
                        "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings", 
                         "Cousins", "Second Cousins", "Unrelated")
#names(relationship_labels) <- relationship_order

loci_set_order <- c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29")
loci_set_labels <- c("Core 13", "Identifiler 15", "Expanded 20", "Supplementary", "Autosomal 29")
#names(loci_set_labels) <- loci_set_order

# Define color palette for populations (frequency sources)
true_pop_colors <- c(
  "AfAm" = "#E41A1C",     # Red
  "Asian" = "#377EB8",    # Blue
  "Cauc" = "#4DAF4A",     # Green
  "Hispanic" = "#984EA3",  # Purple
  "all" = "#FF7F00"       # Orange
)

unnamed_colors <- unname(true_pop_colors)
light_colors_unnamed <- lightness(unnamed_colors, scalefac(1.2))
light_pop_colors <- light_colors_unnamed
names(light_pop_colors) <- names(true_pop_colors)


######## Read in combined_LR files ########
# combined_LR_all.rds made in analyze_lr_outputs.R
all_combined_file <- file.path("output", input_dir, "combined_LR_all.rds")
all_combined <- readRDS(all_combined_file)

######## Make sure LR is numeric ######## 
all_combined <- all_combined %>% mutate(across(c(combined_LR), as.numeric))
all_combined <- all_combined %>% mutate(
  across(c(known_relationship, tested_relationship), 
         ~factor(., levels = relationship_order)),
  loci_set = factor(loci_set, levels = loci_set_order)
)

cutoffs <- calculate_cutoffs(all_combined, c(1, 0.1, 0.01))
cutoffs_file <- file.path(output_dir, "combined_LR_cutoffs.csv")
fwrite(cutoffs, cutoffs_file)

proportions_exceeding_cutoffs <- calculate_proportions_exceeding_cutoffs(all_combined, cutoffs)
proportions_file <- file.path(output_dir, "combined_LR_exceeding_cutoffs.csv")
fwrite(proportions_exceeding_cutoffs, proportions_file)

# ============================================================================
# SECTION 1: CLASSIFICATION AND STATISTICAL ANALYSIS
# ============================================================================

# Add classification column (from simulation_analysis.Rmd)
proportions_all <- proportions_exceeding_cutoffs %>%
  mutate(
    classification = case_when(
      as.character(known_relationship) == as.character(tested_relationship) ~ "True Positive",
      known_relationship == "unrelated" ~ "Unrelated FP",
      TRUE ~ "Related FP"
    ),
    n_loci = case_when(
      loci_set == "core_13" ~ 13,
      loci_set == "identifiler_15" ~ 15,
      loci_set == "expanded_20" ~ 20,
      loci_set == "supplementary" ~ 23,
      loci_set == "autosomal_29" ~ 29
    ),
    prop_LR_gt_1 = proportion_exceeding_fixed,
    prop_LR_gt_10 = proportion_exceeding_1,
    prop_LR_gt_100 = proportion_exceeding_0_1,
    prop_LR_gt_1000 = proportion_exceeding_0_01
  )

# convert to data.table
setDT(proportions_all)

# Save enhanced proportions
fwrite(proportions_all, file.path(output_dir, "proportions_with_classification.csv"))

# --- Statistical Test 1: Chi-square for loci effect ---
log_message("Running chi-square tests for loci effect...")

test_loci_effect <- function(data, known_rel, tested_rel, pop, threshold_col = "prop_LR_gt_1") {
  subset <- data[known_relationship == known_rel & 
                   tested_relationship == tested_rel & 
                   population == pop]
  
  if (nrow(subset) < 2) return(NULL)
  
  n_exceed <- round(subset[[threshold_col]] * subset$n_related)
  n_not_exceed <- subset$n_related - n_exceed
  
  if (any(n_exceed < 0) || any(n_not_exceed < 0)) return(NULL)
  if (sum(n_exceed) == 0 || sum(n_not_exceed) == 0) return(NULL)
  
  contingency <- matrix(c(n_exceed, n_not_exceed), nrow = 2, byrow = TRUE,
                        dimnames = list(c("Exceed", "Not_Exceed"), as.character(subset$loci_set)))
  
  test_result <- tryCatch(
    chisq.test(contingency, simulate.p.value = TRUE, B = 2000),
    error = function(e) NULL
  )
  
  if (is.null(test_result)) return(NULL)
  
  data.frame(
    population = pop,
    known_relationship = known_rel,
    tested_relationship = tested_rel,
    statistic = test_result$statistic,
    p_value = test_result$p.value,
    significant = test_result$p.value < 0.05
  )
}

loci_tests <- rbindlist(lapply(populations, function(pop) {
  rbindlist(lapply(relationship_order, function(known_rel) {
    rbindlist(lapply(c("parent_child", "full_siblings"), function(tested_rel) {
      test_loci_effect(proportions_all, known_rel, tested_rel, pop)
    }))
  }))
}))

fwrite(loci_tests, file.path(output_dir, "stat_test_loci_effect.csv"))

# --- Statistical Test 2: Chi-square for population effect ---
log_message("Running chi-square tests for population effect...")

test_pop_effect <- function(data, known_rel, tested_rel, threshold_col = "prop_LR_gt_1") {
  subset <- data[known_relationship == known_rel & 
                   tested_relationship == tested_rel & 
                   loci_set == "autosomal_29"]
  
  if (nrow(subset) < 2) return(NULL)
  
  n_exceed <- round(subset[[threshold_col]] * subset$n_related)
  n_not_exceed <- subset$n_related - n_exceed
  
  if (any(n_exceed < 0) || any(n_not_exceed < 0)) return(NULL)
  if (sum(n_exceed) == 0 || sum(n_not_exceed) == 0) return(NULL)
  
  contingency <- matrix(c(n_exceed, n_not_exceed), nrow = 2, byrow = TRUE,
                        dimnames = list(c("Exceed", "Not_Exceed"), as.character(subset$population)))
  
  test_result <- tryCatch(
    chisq.test(contingency, simulate.p.value = TRUE, B = 2000),
    error = function(e) NULL
  )
  
  if (is.null(test_result)) return(NULL)
  
  data.frame(
    known_relationship = known_rel,
    tested_relationship = tested_rel,
    statistic = test_result$statistic,
    p_value = test_result$p.value,
    significant = test_result$p.value < 0.05
  )
}

pop_tests <- rbindlist(lapply(relationship_order, function(known_rel) {
  rbindlist(lapply(c("parent_child", "full_siblings"), function(tested_rel) {
    test_pop_effect(proportions_all, known_rel, tested_rel)
  }))
}))

fwrite(pop_tests, file.path(output_dir, "stat_test_population_effect.csv"))

# --- Statistical Test 3: Linear trend analysis ---
log_message("Running linear trend tests...")

trend_tests <- proportions_all[, {
  if (.N < 3) return(NULL)
  
  model <- tryCatch(lm(prop_LR_gt_1 ~ n_loci), error = function(e) NULL)
  if (is.null(model)) return(NULL)
  
  summary_model <- summary(model)
  coef_table <- coef(summary_model)
  
  if (nrow(coef_table) < 2) return(NULL)
  
  list(
    slope = coef_table[2, 1],
    slope_se = coef_table[2, 2],
    p_value = coef_table[2, 4],
    r_squared = summary_model$r.squared,
    trend = ifelse(coef_table[2, 1] > 0, "INCREASING", "DECREASING"),
    significant = coef_table[2, 4] < 0.05
  )
}, by = .(population, known_relationship, tested_relationship)]

fwrite(trend_tests, file.path(output_dir, "stat_test_linear_trend.csv"))

log_message("Statistical tests completed and saved.")

# ============================================================================
# SECTION 2: CLASSIFICATION SUMMARY PLOT (from simulation_analysis.Rmd))
# ============================================================================
log_message("Creating classification summary plot...")

summary_29 <- proportions_all[loci_set == "autosomal_29"]

classification_plot <- ggplot(summary_29, aes(x = known_relationship, y = prop_LR_gt_1, fill = classification)) +
  geom_col(position = "dodge") +
  facet_grid(tested_relationship ~ population,
             labeller = labeller(population = as_labeller(c(
               "AfAm" = "African American",
               "Asian" = "Asian",
               "Cauc" = "Caucasian",
               "Hispanic" = "Hispanic",
               "all" = "All Populations"
             )),
             tested_relationship = as_labeller(c(
               "parent_child" = "Parent-Child",
               "full_siblings" = "Full Siblings"
             )))) +
  scale_fill_manual(values = c("True Positive" = "#2ca02c", 
                               "Related FP" = "#ff7f0e", 
                               "Unrelated FP" = "#d62728")) +
  scale_x_discrete(labels = function(x) str_wrap(relationship_labels[match(x, relationship_order)], width = 8)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Classification Performance at 29 Autosomal Loci",
    subtitle = "Rows = Tested Hypothesis, Columns = Population",
    x = "True Relationship",
    y = "Proportion with LR > 1",
    fill = "Classification"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "classification_summary_29loci.png"), 
       plot = classification_plot, width = 14, height = 6, dpi = 300, bg = "white")

log_message("Classification summary plot saved.")


######## Bar chart plots ######## 
# Open PDF device to save all plots to one PDF
pdf(file.path(output_dir, "population_mismatch_proportions_analysis.pdf"), width = 14, height = 10)

# For relationship mismatch analysis (mimicking the single image)
relationships_to_test <- c("parent_child", "full_siblings")

for (rel in relationships_to_test) {
  # Filter to include ALL relationships (both correct and incorrect)
  rel_data <- proportions_exceeding_cutoffs %>% 
    filter(tested_relationship == rel,           # Testing this relationship
           tested_population == population)      # Using correct population frequencies
  
  # Get counts
  counts_summary <- rel_data %>%
    group_by(population, known_relationship, loci_set) %>%
    summarize(n_pairs = first(n_related), .groups = 'drop')
  
  cat("For", rel, " hypothesis each population x known_relationship x loci set has:", unique(counts_summary$n_pairs), "pairs\n")
  
  # Create long format data for plotting
  proportions_long <- rel_data %>%
    pivot_longer(cols = starts_with("proportion_exceeding"),
                 names_to = "Cutoff_Type", values_to = "Proportion",
                 names_prefix = "proportion_exceeding_")
  
  proportions_long$Cutoff_Type <- factor(proportions_long$Cutoff_Type,
                                         levels = c("fixed", "1", "0_1","0_01"),
                                         labels = c("Fixed (1.00)", "1% FPR", "0.1% FPR", "0.01% FPR"))
  
  proportions_long <- proportions_long %>%
    mutate(
      known_relationship = recode_factor(known_relationship,
                                         "parent_child" = "Parent-Child",
                                         "full_siblings" = "Full Siblings", 
                                         "half_siblings" = "Half Siblings",
                                         "cousins" = "Cousins",
                                         "second_cousins" = "Second Cousins",
                                         "unrelated" = "Unrelated"),
      loci_set = recode_factor(loci_set,
                               "core_13" = "Core 13",
                               "identifiler_15" = "Identifiler 15",
                               "expanded_20" = "Expanded 20", 
                               "supplementary" = "Supplementary",
                               "autosomal_29" = "Autosomal 29"),
      population = factor(population, levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all", "unrelated"))
    )
  
  
  # Plot with known_relationship on x-axis (this creates the variation!)
  p <- ggplot(proportions_long, aes(x = known_relationship, y = Proportion, fill = population)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_grid(Cutoff_Type ~ loci_set, scales = "fixed", labeller = labeller(loci_set = loci_set_labels)) +
    scale_fill_manual(values = true_pop_colors, name = "True Population") +
    labs(
      title = paste("Proportions exceeding likelihood cutoffs:", rel, "Hypothesis"),
      subtitle = paste("Using population-matched allele frequencies"),
      x = "True Relationship",
      y = "Proportion Exceeding Cutoff"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "darkblue"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  print(p)
  
  # Save individual PNG
  filename <- paste0("relationship_analysis_", rel, "_hypot_", today, ".png")
  ggsave(file.path(output_dir, filename), plot = p, width = 14, height = 10, bg = "white", dpi = 300)
}

# ============================================================================
# SECTION 3: DETAILED RATES TABLE
# ============================================================================

log_message("Creating detailed rates table...")

# Create detailed table page for the PDF
rates_table <- proportions_all[
  loci_set == "autosomal_29",
  .(population, known_relationship, tested_relationship, classification, 
    prop_LR_gt_1, prop_LR_gt_10, prop_LR_gt_100, prop_LR_gt_1000)
]

setorder(rates_table, tested_relationship, population, classification, known_relationship)

# Create a text plot for the PDF
par(mar = c(0, 0, 2, 0))
plot.new()
title("Detailed Rates at 29 Autosomal Loci", cex.main = 1.5, font.main = 2)

# Add table as text (basic visualization for PDF)
text_lines <- capture.output(print(rates_table, row.names = FALSE))
n_lines <- length(text_lines)
y_positions <- seq(0.9, 0.1, length.out = min(n_lines, 50))

for (i in seq_along(text_lines[1:min(50, n_lines)])) {
  text(0.5, y_positions[i], text_lines[i], family = "mono", cex = 0.5, adj = 0.5)
}

# Save as separate CSV
fwrite(rates_table, file.path(output_dir, "detailed_rates_29loci.csv"))

log_message("Detailed rates table saved.")

dev.off()

# ============================================================================
# SECTION 4: KEY FINDINGS SUMMARY
# ============================================================================

log_message("Generating key findings summary...")

key_findings_file <- file.path(output_dir, "key_findings_summary.txt")
sink(key_findings_file)

cat("=== KEY FINDINGS SUMMARY ===\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. TRUE POSITIVE RATES (at 29 loci, LR > 1):\n")
tp_rates <- proportions_all[loci_set == "autosomal_29" & 
                              as.character(known_relationship) == as.character(tested_relationship),
                            .(population, tested_relationship, prop_LR_gt_1)]
print(dcast(tp_rates, tested_relationship ~ population, value.var = "prop_LR_gt_1"))

cat("\n2. UNRELATED FALSE POSITIVE RATES (at 29 loci, LR > 1):\n")
unrel_fp <- proportions_all[loci_set == "autosomal_29" & 
                              known_relationship == "unrelated",
                            .(population, tested_relationship, prop_LR_gt_1)]
print(dcast(unrel_fp, tested_relationship ~ population, value.var = "prop_LR_gt_1"))

cat("\n3. HALF-SIBLING AS FULL-SIBLING FALSE ID (at 29 loci):\n")
hs_fp <- proportions_all[loci_set == "autosomal_29" & 
                           known_relationship == "half_siblings" &
                           tested_relationship == "full_siblings",
                         .(population, prop_LR_gt_1, prop_LR_gt_10, prop_LR_gt_100)]
print(hs_fp)

cat("\n4. TREND DIRECTIONS (Full-Siblings hypothesis, half-siblings):\n")
if (nrow(trend_tests) > 0) {
  hs_trend <- trend_tests[tested_relationship == "full_siblings" & known_relationship == "half_siblings"]
  print(hs_trend[, .(population, trend, slope, p_value, significant)])
} else {
  cat("  No trend data available\n")
}

cat("\n5. STATISTICAL TEST SUMMARY:\n")
cat("\nLoci Effect Tests (significant at p < 0.05):\n")
if (nrow(loci_tests) > 0) {
  print(loci_tests[significant == TRUE, .(population, known_relationship, tested_relationship, p_value)])
} else {
  cat("  No significant loci effects detected\n")
}

cat("\nPopulation Effect Tests (significant at p < 0.05):\n")
if (nrow(pop_tests) > 0) {
  print(pop_tests[significant == TRUE, .(known_relationship, tested_relationship, p_value)])
} else {
  cat("  No significant population effects detected\n")
}

cat("\n=== END OF SUMMARY ===\n")

sink()

log_message(paste("Key findings saved to:", key_findings_file))
log_message("All analyses complete!")