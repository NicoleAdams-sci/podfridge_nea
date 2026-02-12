#!/usr/bin/env Rscript

#### Plots of combined LR distribution for tested parent-child and full-sibling relationships

# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)  # Includes ggplot2, dplyr, tidyr, etc.
  library(data.table) # Required for efficient fread()
  library(scales)     # For number formatting
}))

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


# Define factors, colors, and labels
populations <- c("AfAm", "Cauc", "Hispanic", "Asian", "all")

relationship_order <- c("parent_child", "full_siblings", "half_siblings", 
                        "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings", 
                         "Cousins", "Second Cousins", "Unrelated")
names(relationship_labels) <- relationship_order

loci_set_order <- c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29")
loci_set_labels <- c("Core 13", "Identifiler 15", "Expanded 20", "Supplementary", "Autosomal 29")
names(loci_set_labels) <- loci_set_order

# Define color palette for populations (frequency sources)
freq_pop_colors <- c(
  "AfAm" = "#E41A1C",     # Red
  "Asian" = "#377EB8",    # Blue
  "Cauc" = "#4DAF4A",     # Green
  "Hispanic" = "#984EA3",  # Purple
  "all" = "#FF7F00"       # Orange
)


######## Read in combined_LR files ########
# combined_LR_all.rds made in analyze_lr_outputs.R
all_combined_file <- file.path(input_dir, "combined_LR_all.rds")
all_combined <- readRDS(all_combined_file)

# Make sure LR is numeric
all_combined <- all_combined %>% mutate(across(c(combined_LR), as.numeric))
all_combined <- all_combined %>% mutate(
  known_relationship = recode_factor(known_relationship,
                                     "parent_child" = "Parent-Child",
                                     "full_siblings" = "Full Siblings", 
                                     "half_siblings" = "Half Siblings",
                                     "cousins" = "Cousins",
                                     "second_cousins" = "Second Cousins",
                                     "unrelated" = "Unrelated"),
  tested_relationship = recode_factor(tested_relationship,  # Add this line!
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
  population = factor(population, levels = c("AfAm", "Asian", "Cauc", "Hispanic", "all"))
)

######## Relationship mismatch boxplots ########
plot3 <- ggplot(all_combined, aes(x = tested_relationship, y = combined_LR, fill = population)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  facet_grid(known_relationship ~ loci_set) +
  scale_fill_manual(values = freq_pop_colors) +
  labs(
    title = "LR Distributions for Relationship Comparisons",
    x = "Tested Relationship Type",
    y = "combined LR",
    fill = "Population"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_log10()
plot3.name <- paste0("relationship_mismatch_LRboxplot.png")
ggsave(file.path(output_dir, plot3.name), plot = plot3, width = 14, height = 11, bg = "white")


# Add an indication of the correct relationship - rectangles on the diagonal 
plot3.alt <- ggplot(all_combined, aes(x = tested_relationship, y = combined_LR, fill = population)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  facet_grid(known_relationship ~ loci_set) +
  # Add gold/yellow rectangles to highlight correct matches
  geom_rect(data = all_combined %>% filter(tested_relationship == as.character(known_relationship)),
            aes(xmin = as.integer(known_relationship)-0.5, xmax = as.integer(known_relationship)+0.5, ymin = 0.1, ymax = 1e36),
            fill = NA, alpha = 0.2, color = "orange", 
            linetype = "dashed", linewidth = 0.7,
            inherit.aes = FALSE) +
  scale_fill_manual(values = freq_pop_colors) +
  scale_y_log10() +
  labs(
    title = "LR Distributions for Relationship Comparisons",
    subtitle = "Yellow highlights show when tested relationship matches true relationship",
    x = "Tested Relationship Type",
    y = "Combined LR (log scale)",
    fill = "Population"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )  



######## Population mismatch boxplots ########
# Define facet grid
facet_formula <- as.formula("known_relationship ~ loci_set")

# For relationship mismatch analysis (mimicking the single image)
relationships_to_test <- c("Parent-Child", "Full Siblings")

# Open PDF device to save all plots to one PDF
pdf(file.path(output_dir, "mismatched_pops_all_relationships_LRboxplots.pdf"), width = 11, height = 8)

for (rel in relationships_to_test) {
  rel.df <-  all_combined %>% filter(tested_relationship == rel)
  
  for (pop in populations) {
    rel_pop <- rel.df %>% filter(tested_population == pop)
    
    # dist of combined_LR
    box.p <- ggplot(rel_pop, aes(x=population, y=combined_LR, fill=population)) +
      geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7, outlier.size = 0.5) +
      facet_grid(facet_formula, scales = "fixed", labeller = labeller(
        loci_set = loci_set_labels,
        known_relationship = relationship_labels,
        tested_relationship = relationship_labels)) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_fill_manual(values = freq_pop_colors) +
      ggtitle(paste(rel, " tested relationship for", pop, "allele frequencies")) +
      labs(
        x = "Known population",
        y = "Log likelihood ratio (LR)",
        #fill = "Known Population"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        #plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 14),
        #legend.title = element_text(size = 12, margin = margin(t = 10)),
        # legend.text = element_text(size = 10),
        legend.position = "none",
        # legend.key.width = unit(2, "cm"),
        # legend.key.height = unit(0.5, "cm"),
        # legend.box.margin = margin(t = 10),
        strip.text = element_text(size = 12)
      )
    
    # Save to PDF (add to the open PDF file)
    print(box.p)
    
    box.name <- paste0(gsub(" ", "_", rel), "_LRboxplot_", pop, ".png")
    ggsave(file.path(output_dir, box.name), plot = box.p, width = 11, height = 8, bg = "white")
    
  }
}
# Close PDF device
dev.off()