#!/usr/bin/env Rscript

# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)  
}))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript code/plots_individual_locus_NEA.R <input_dir> [output_dir]")
}

input_dir <- args[1]
log_message(paste("Input directory:", input_dir))

# If output directory is provided, use it, otherwise create one based on the input directory
if (length(args) >= 2) {
  output_subdir <- args[2]
  output_dir <- file.path("output", output_subdir)
} else {
  timestamp <- format(Sys.time(), "%Y%m%d")
  output_dir <- file.path("output", paste0("individual_locus_plots_", timestamp))
}
log_message(paste("Output directory:", output_dir))

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define populations
populations <- c("AfAm", "Cauc", "Hispanic", "Asian")

# Define relationship types and ensure correct ordering for plotting
relationship_order <- c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings", "Cousins", "Second Cousins", "Unrelated")

# Function to safely read a CSV file
safe_read_csv <- function(file_path) {
  if (!file.exists(file_path)) {
    log_message(paste("File not found:", file_path))
    return(NULL)
  }
  
  tryCatch({
    log_message(paste("Attempting to read:", file_path))
    data <- read.csv(file_path)
    log_message(paste("Successfully read file with", nrow(data), "rows"))
    return(data)
  }, error = function(e) {
    log_message(paste("Error reading file:", e$message))
    return(NULL)
  })
}

# Load raw genotype data with individual locus information
load_raw_genotype_data <- function() {
  log_message("Loading raw genotype data for individual locus analysis...")
  
  # Initialize empty data frame
  all_raw_data <- data.frame()
  
  # Track how many files were successfully loaded
  files_loaded <- 0
  
  for (pop in populations) {
    log_message(paste("Processing population:", pop))
    
    # Look for raw genotype files in the input directory structure
    # First try the population-specific directory pattern
    pop_dir <- file.path(input_dir, paste0(pop, "_summary"))
    
    if (!dir.exists(pop_dir)) {
      # Try alternative directory patterns
      pop_dir_pattern <- paste0("summary_", pop, "_")
      potential_dirs <- list.dirs(input_dir, recursive = FALSE)
      matching_dirs <- potential_dirs[grepl(pop_dir_pattern, basename(potential_dirs))]
      
      if (length(matching_dirs) > 0) {
        pop_dir <- matching_dirs[1]
        log_message(paste("Found directory using pattern match:", pop_dir))
      } else {
        log_message(paste("No directory found for population:", pop))
        next
      }
    }
    
    # Look for raw genotype files
    # These might be in subdirectories of the population summary directory
    # or we might need to look in the original simulation directories
    raw_files <- list.files(pop_dir, pattern = "sim_processed_genotypes.*\\.csv$", 
                           full.names = TRUE, recursive = TRUE)
    
    if (length(raw_files) == 0) {
      # Try looking in the main input directory for combined files
      raw_files <- list.files(input_dir, 
                             pattern = paste0("sim_processed_genotypes_", pop, "_.*_combined\\.csv$"), 
                             full.names = TRUE, recursive = FALSE)
    }
    
    if (length(raw_files) == 0) {
      log_message(paste("No raw genotype files found for population:", pop))
      next
    }
    
    # Load and combine all files for this population
    for (file_path in raw_files) {
      log_message(paste("Loading file:", basename(file_path)))
      data <- safe_read_csv(file_path)
      
      if (!is.null(data) && nrow(data) > 0) {
        # Check if this file has the required columns
        required_cols <- c("population", "relationship_type", "sim_id", "locus", 
                          "LR_AfAm", "LR_Cauc", "LR_Hispanic", "LR_Asian")
        
        if (all(required_cols %in% names(data))) {
          # Filter for the current population
          pop_data <- data %>% filter(population == pop)
          all_raw_data <- bind_rows(all_raw_data, pop_data)
          files_loaded <- files_loaded + 1
          log_message(paste("Added", nrow(pop_data), "rows for population", pop))
        } else {
          log_message(paste("File missing required columns:", file_path))
          log_message(paste("Available columns:", paste(names(data), collapse = ", ")))
        }
      }
    }
  }
  
  # Check if we have data
  if (nrow(all_raw_data) == 0) {
    log_message("ERROR: No raw genotype data was loaded. Please check the input directory and file patterns.")
    stop("No data loaded. Check log for details.")
  }
  
  log_message(paste("Successfully loaded data from", files_loaded, "files"))
  log_message(paste("Total rows in combined data:", nrow(all_raw_data)))
  
  return(all_raw_data)
}

# Function to transform raw data to long format for plotting
prepare_locus_data <- function(raw_data) {
  log_message("Preparing locus-level data for plotting...")
  
  # Melt the LR columns to create a long format
  lr_cols <- c("LR_AfAm", "LR_Cauc", "LR_Hispanic", "LR_Asian")
  
  long_data <- raw_data %>%
    select(population, relationship_type, sim_id, locus, all_of(lr_cols)) %>%
    pivot_longer(cols = all_of(lr_cols), 
                 names_to = "freq_population", 
                 values_to = "LR", 
                 names_prefix = "LR_") %>%
    mutate(
      is_correct_pop = (population == freq_population),
      relationship_type = factor(relationship_type, 
                                levels = relationship_order,
                                labels = relationship_labels),
      population = factor(population, levels = populations),
      freq_population = factor(freq_population, levels = populations)
    ) %>%
    filter(!is.na(LR), is.finite(LR), LR > 0)  # Remove invalid LR values
  
  log_message(paste("Prepared long-format data with", nrow(long_data), "rows"))
  
  # Load Core Loci Data for ordering
  log_message("Loading core loci data for ordering...")
  core_loci_time <- system.time({
    core_loci <- read.csv("data/core_CODIS_loci.csv")
  
    # Create ordered list of loci based on CODIS sets
    core_13_loci <- core_loci %>% filter(core_13 == 1) %>% pull(locus)
    identifiler_15_additional <- core_loci %>% filter(identifiler_15 == 1, core_13 == 0) %>% pull(locus)
    expanded_20_additional <- core_loci %>% filter(expanded_20 == 1, identifiler_15 == 0) %>% pull(locus)
    supplementary_additional <- core_loci %>% filter(supplementary == 1, expanded_20 == 0) %>% pull(locus)
    
    # Get all loci from the data and find any that aren't in the above sets
    all_data_loci <- sort(unique(long_data$locus))
    autosomal_29_additional <- setdiff(all_data_loci, c(core_13_loci, identifiler_15_additional, 
                                                       expanded_20_additional, supplementary_additional))
    
      # Create the final ordered loci list
      ordered_loci <- c(core_13_loci, identifiler_15_additional, expanded_20_additional, 
                       supplementary_additional, autosomal_29_additional)
      
      # Filter to only include loci that are actually in our data
      ordered_loci <- ordered_loci[ordered_loci %in% all_data_loci]
  })
  log_message(paste("Created ordered loci list with", length(ordered_loci), "loci in", core_loci_time["elapsed"], "seconds."))

  return(list(long_data = long_data, ordered_loci = ordered_loci))
}

# Function to create individual locus plots for each population
create_individual_locus_plots <- function(long_data, ordered_loci) {
  log_message("Creating individual locus plots for each population...")
  
  plot_objects <- list()
  
  # Get unique loci, ordered by CODIS sets
  unique_loci <- ordered_loci[ordered_loci %in% unique(long_data$locus)]
  log_message(paste("Found", length(unique_loci), "unique loci in CODIS order"))
  
  # Define color palette for frequency populations
  freq_pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(freq_pop_colors) <- levels(long_data$freq_population)
  
  # Create plots for each population
  for (pop in populations) {
    log_message(paste("Creating plot for population:", pop))
    
    # Filter data for this population and exclude correct population matches
    pop_data <- long_data %>% 
      filter(population == pop, is_correct_pop == FALSE)
    
    if (nrow(pop_data) == 0) {
      log_message(paste("No mismatched data available for population:", pop))
      next
    }
    
    # Calculate sample sizes for subtitle
    sample_sizes <- pop_data %>%
      group_by(relationship_type) %>%
      summarise(n_samples = n_distinct(sim_id), .groups = 'drop')
    
    total_samples <- sum(sample_sizes$n_samples)
    
    # Create the plot
    tryCatch({
      # Determine number of columns for facet_wrap based on number of loci
      n_loci <- length(unique(pop_data$locus))
      n_cols <- if (n_loci <= 15) 5 else if (n_loci <= 25) 6 else 7
      
      p <- ggplot(pop_data, aes(x = relationship_type, y = LR, fill = freq_population)) +
        geom_boxplot(position = position_dodge(width = 0.85), 
                     alpha = 0.7, 
                     outlier.size = 0.3,
                     size = 0.3) +
        facet_wrap(~ locus, scales = "fixed", ncol = n_cols) +
        scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_fill_manual(values = freq_pop_colors, 
                         name = "Assumed Ancestry\n(Frequency Source)",
                         labels = c("African American", "Caucasian", "Hispanic", "Asian")) +
        labs(
          title = paste("LR Distributions by Individual Locus:", pop, "Population"),
          subtitle = paste("Using Mismatched Population Allele Frequencies (n =", 
                          format(total_samples, big.mark = ","), "individuals)"),
          x = "True Relationship Type",
          y = "Likelihood Ratio (log10 scale)"
        ) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9),
          legend.position = "bottom",
          legend.key.width = unit(1.5, "cm"),
          legend.key.height = unit(0.4, "cm"),
          strip.text = element_text(size = 8, face = "bold"),
          panel.spacing = unit(0.5, "lines"),
          legend.box.margin = margin(t = 10)
        ) +
        guides(fill = guide_legend(nrow = 1))
      
      # Save the plot
      filename <- paste0("individual_locus_lr_", pop, "_population.png")
      
      # Adjust plot dimensions based on number of loci
      plot_width <- if (n_loci <= 15) 16 else if (n_loci <= 25) 18 else 20
      plot_height <- ceiling(n_loci / n_cols) * 2.5 + 4  # Dynamic height based on rows
      
      plot_objects[[paste0("lr_", pop)]] <- p
      ggsave(file.path(output_dir, filename), plot = p, width = plot_width, height = plot_height, bg = "white", dpi = 300)
      
      log_message(paste("Saved plot:", filename))
      
    }, error = function(e) {
      log_message(paste("Error creating plot for population", pop, ":", e$message))
    })
  }
  
  return(plot_objects)
}
  
# Function to create ratio boxplots for each population by locus
create_locus_ratio_boxplots <- function(long_data, ordered_loci) {
  log_message("Creating locus-level ratio boxplots for each population...")
  
  ratio_plot_objects <- list() 
  
  # Prepare ratio data by calculating correct vs wrong population LRs
  correct_lrs <- long_data %>%
    filter(is_correct_pop == TRUE) %>%
    select(population, relationship_type, sim_id, locus, correct_LR = LR)
  
  wrong_lrs <- long_data %>%
    filter(is_correct_pop == FALSE) %>%
    select(population, relationship_type, sim_id, locus, freq_population, wrong_LR = LR)
  
  # Merge and calculate ratios
  ratio_data <- wrong_lrs %>%
    left_join(correct_lrs, by = c("population", "relationship_type", "sim_id", "locus")) %>%
    mutate(ratio = wrong_LR / correct_LR) %>%
    filter(is.finite(ratio))
  
  if (nrow(ratio_data) == 0) {
    log_message("No valid ratio data available for boxplots")
    return(NULL)
  }
  
  # Define color palette for frequency populations
  freq_pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(freq_pop_colors) <- levels(ratio_data$freq_population)
  
  # Create boxplots for each population
  for (pop in populations) {
    log_message(paste("Creating ratio boxplot for population:", pop))
    
    pop_ratio_data <- ratio_data %>% 
      filter(population == pop)
    
    if (nrow(pop_ratio_data) == 0) {
      log_message(paste("No ratio data available for population:", pop))
      next
    }
    
    # Calculate sample sizes for subtitle
    sample_sizes <- pop_ratio_data %>%
      group_by(relationship_type) %>%
      summarise(n_samples = n_distinct(sim_id), .groups = 'drop')
    
    total_samples <- sum(sample_sizes$n_samples)
    
    tryCatch({
      # Order loci by CODIS sets for consistency
      pop_ratio_data$locus <- factor(pop_ratio_data$locus, 
                                    levels = ordered_loci)
      
      # Determine number of columns for facet_wrap based on number of loci
      n_loci <- length(unique(pop_ratio_data$locus))
      n_cols <- if (n_loci <= 15) 5 else if (n_loci <= 25) 6 else 7
      
      p <- ggplot(pop_ratio_data, aes(x = relationship_type, y = ratio, fill = freq_population)) +
        geom_boxplot(position = position_dodge(width = 0.85), 
                     alpha = 0.7, 
                     outlier.size = 0.3,
                     size = 0.3) +
        facet_wrap(~ locus, scales = "fixed", ncol = n_cols) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_fill_manual(values = freq_pop_colors, 
                         name = "Assumed Ancestry\n(Frequency Source)",
                         labels = c("African American", "Caucasian", "Hispanic", "Asian")) +
        labs(
          title = paste("LR Ratios by Individual Locus:", pop, "Population"),
          subtitle = paste("Ratio of Mismatched Population LR to Correct Population LR (n =", 
                          format(total_samples, big.mark = ","), "individuals)"),
          x = "True Relationship Type",
          y = "Wrong LR / Correct LR (log10 scale)"
        ) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9),
          legend.position = "bottom",
          legend.key.width = unit(1.5, "cm"),
          legend.key.height = unit(0.4, "cm"),
          strip.text = element_text(size = 8, face = "bold"),
          panel.spacing = unit(0.5, "lines"),
          legend.box.margin = margin(t = 10)
        ) +
        guides(fill = guide_legend(nrow = 1))
      
      # Save the plot
      filename <- paste0("individual_locus_ratio_boxplot_", pop, "_population.png")
      
      # Adjust plot dimensions based on number of loci
      plot_width <- if (n_loci <= 15) 16 else if (n_loci <= 25) 18 else 20
      plot_height <- ceiling(n_loci / n_cols) * 2.5 + 4  # Dynamic height based on rows
      
      ratio_plot_objects[[paste0("ratio_", pop)]] <- p
      ggsave(file.path(output_dir, filename), plot = p, width = plot_width, height = plot_height, bg = "white", dpi = 300)
      
      log_message(paste("Saved ratio boxplot:", filename))
      
    }, error = function(e) {
      log_message(paste("Error creating ratio boxplot for population", pop, ":", e$message))
    })
  }
  
  return(ratio_plot_objects)
}

# Function to create heatmap plots for each population showing LR ratios by locus
create_locus_heatmap_plots <- function(long_data, ordered_loci) {
  log_message("Creating locus-level heatmap plots for each population...")
  
  heatmap_plot_objects <- list()
  
  # Prepare ratio data by calculating correct vs wrong population LRs
  correct_lrs <- long_data %>%
    filter(is_correct_pop == TRUE) %>%
    select(population, relationship_type, sim_id, locus, correct_LR = LR)
  
  wrong_lrs <- long_data %>%
    filter(is_correct_pop == FALSE) %>%
    select(population, relationship_type, sim_id, locus, freq_population, wrong_LR = LR)
  
  # Merge and calculate ratios
  ratio_data <- wrong_lrs %>%
    left_join(correct_lrs, by = c("population", "relationship_type", "sim_id", "locus")) %>%
    mutate(ratio = wrong_LR / correct_LR) %>%
    filter(is.finite(ratio))
  
  if (nrow(ratio_data) == 0) {
    log_message("No valid ratio data available for heatmaps")
    return(NULL)
  }
  
  # Calculate mean ratios for heatmap
  heatmap_data <- ratio_data %>%
    group_by(population, locus, relationship_type, freq_population) %>%
    summarise(mean_ratio = mean(ratio, na.rm = TRUE), .groups = 'drop') %>%
    filter(is.finite(mean_ratio))
  
  # Calculate global scale limits across ALL populations for consistent scaling
  global_min_ratio <- min(heatmap_data$mean_ratio, na.rm = TRUE)
  global_max_ratio <- max(heatmap_data$mean_ratio, na.rm = TRUE)
  
  log_message(paste("Global ratio range:", round(global_min_ratio, 4), "to", round(global_max_ratio, 4)))
  
  # Define consistent color scale based on global data range
  log_min_val <- log10(global_min_ratio)
  log_max_val <- log10(global_max_ratio)
  
  # Use full color palette and calculate white position
  my_colors_full_palette <- c("blue", "green", "white", "yellow", "red")
  val_for_white_rc <- scales::rescale(0, from = c(log_min_val, log_max_val), to = c(0,1))
  
  if (val_for_white_rc <= 0) { 
    # All ratios > 1, use white-yellow-red
    final_colors <- c("white", "yellow", "red")
    final_values <- c(0, 0.5, 1)
  } else if (val_for_white_rc >= 1) {
    # All ratios < 1, use blue-green-white
    final_colors <- c("blue", "green", "white")
    final_values <- c(0, 0.5, 1)
  } else {
    # Ratios span across 1, use full palette
    final_colors <- my_colors_full_palette
    final_values <- c(0,                            
                      val_for_white_rc * 0.5,       
                      val_for_white_rc,             
                      val_for_white_rc + (1 - val_for_white_rc) * 0.5, 
                      1)                            
    final_values <- sort(unique(final_values))
  }
  
  # Define comprehensive breaks and labels for the legend
  all_breaks <- c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000)
  all_labels <- c("0.001", "0.003", "0.01", "0.03", "0.1", "0.3", "1", "3", "10", "30", "100", "300", "1000")
  
  # Filter breaks to only include those within the actual data range (with some padding)
  valid_breaks <- all_breaks[all_breaks >= global_min_ratio * 0.8 & all_breaks <= global_max_ratio * 1.2]
  valid_labels <- all_labels[all_breaks >= global_min_ratio * 0.8 & all_breaks <= global_max_ratio * 1.2]
  
  # Ensure we always include 1 (white/no change) if it's anywhere near our range
  if (global_min_ratio <= 1 & global_max_ratio >= 1) {
    if (!1 %in% valid_breaks) {
      valid_breaks <- sort(c(valid_breaks, 1))
      valid_labels <- sort(c(valid_labels, "1"))
    }
  }
  
  log_message(paste("Using breaks:", paste(valid_breaks, collapse=", ")))
  
  # Create heatmap for each population
  for (pop in populations) {
    log_message(paste("Creating heatmap for population:", pop))
    
    pop_heatmap_data <- heatmap_data %>% 
      filter(population == pop)
    
    if (nrow(pop_heatmap_data) == 0) {
      log_message(paste("No heatmap data available for population:", pop))
      next
    }
    
    # Ensure complete data for all combinations
    pop_heatmap_data <- pop_heatmap_data %>%
      complete(locus, relationship_type, freq_population)
    
    tryCatch({
      # Order loci by CODIS sets for consistency
      pop_heatmap_data$locus <- factor(pop_heatmap_data$locus, 
                                      levels = ordered_loci)
      
      p <- ggplot(pop_heatmap_data, aes(x = locus, y = relationship_type, fill = mean_ratio)) +
        geom_tile(color = "grey80", size = 0.1) +
        facet_wrap(~ freq_population, ncol = 2) +
        scale_fill_gradientn(
          colors = final_colors,
          values = final_values,
          trans = "log10",
          breaks = valid_breaks,
          labels = valid_labels,
          limits = c(global_min_ratio, global_max_ratio),  # Use global limits
          na.value = "grey75",
          guide = guide_colorbar(
            barwidth = 25,
            barheight = 1.5,
            title.position = "top",
            title.hjust = 0.5,
            label.theme = element_text(size = 8),
            ticks.colour = "black",
            ticks.linewidth = 0.5,
            frame.colour = "black",
            frame.linewidth = 0.5
          )
        ) +
        labs(
          title = paste("Mean LR Ratios by Individual Locus:", pop, "Population"),
          subtitle = "Ratio of Mismatched Population LR to Correct Population LR",
          x = "Locus",
          y = "True Relationship Type",
          fill = "Ratio: LR(Tested Pop) / LR(True Pop) \n(log10 scale, 1 = white = no change)"
        ) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 10, face = "bold", lineheight = 1.2),
          legend.text = element_text(size = 8),
          legend.position = "bottom",
          strip.text = element_text(size = 11, face = "bold", margin = margin(t=3,b=3)),
          panel.spacing = unit(1, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      # Save the heatmap
      filename <- paste0("individual_locus_heatmap_", pop, "_population.png")
      
      # Determine dimensions based on number of loci
      n_loci <- length(unique(pop_heatmap_data$locus))
      plot_width <- max(16, n_loci * 0.4)  # Minimum 16, scale with loci
      plot_height <- 12
      
      heatmap_plot_objects[[paste0("ratio_", pop)]] <- p
      ggsave(file.path(output_dir, filename), 
             plot = p, 
             width = plot_width, 
             height = plot_height, 
             bg = "white", 
             dpi = 300)
      
      log_message(paste("Saved heatmap:", filename))
      
    }, error = function(e) {
      log_message(paste("Error creating heatmap for population", pop, ":", e$message))
    })
  }
  
  return(heatmap_plot_objects)
}

# Function to create a comparison plot showing all populations together
create_comparison_plot <- function(long_data, ordered_loci) {
  log_message("Creating comparison plot with all populations...")
  
  # Filter for mismatched populations only
  comparison_data <- long_data %>% 
    filter(is_correct_pop == FALSE)
  
  if (nrow(comparison_data) == 0) {
    log_message("No mismatched data available for comparison plot")
    return(NULL)
  }
  
  # Calculate median LR by locus and population for ordering
  locus_order <- comparison_data %>%
    group_by(locus) %>%
    summarise(median_lr = median(LR, na.rm = TRUE), .groups = 'drop') %>%
    arrange(median_lr) %>%
    pull(locus)
  
  # Use CODIS ordering instead of median LR ordering
  comparison_data$locus <- factor(comparison_data$locus, levels = ordered_loci)
  
  # Define color palette for populations
  pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(pop_colors) <- levels(comparison_data$population)
  
  tryCatch({
    p <- ggplot(comparison_data, 
                aes(x = locus, y = LR, fill = population)) +
      geom_boxplot(position = position_dodge(width = 0.85), 
                   alpha = 0.7, 
                   outlier.size = 0.2,
                   size = 0.3) +
      facet_wrap(~ relationship_type, scales = "fixed", ncol = 2) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_fill_manual(values = pop_colors,
                       name = "True Population") +
      labs(
        title = "LR Distributions by Locus: All Populations Comparison",
        subtitle = "Using Mismatched Population Allele Frequencies\n(Loci ordered by CODIS sets: Core 13, Identifiler 15, Expanded 20, Supplementary, Autosomal 29)",
        x = "Locus (ordered by CODIS sets)",
        y = "Likelihood Ratio (log10 scale)"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"),
        panel.spacing = unit(0.8, "lines")
      )
    
    ggsave(file.path(output_dir, "individual_locus_lr_all_populations_comparison.png"), 
           plot = p, 
           width = 20, 
           height = 14, 
           bg = "white", 
           dpi = 300)
    
    log_message("Saved comparison plot: individual_locus_lr_all_populations_comparison.png")
    
  }, error = function(e) {
    log_message(paste("Error creating comparison plot:", e$message))
  })
}

# Function to save plot objects to PDF
save_plots_to_pdf <- function(lr_plots, ratio_plots, heatmap_plots) {
  log_message("Creating combined PDF from existing plot objects...")
  
  pdf_filename <- file.path(output_dir, "individual_locus_plots_combined.pdf")
  pdf(pdf_filename, width = 16, height = 12)
  
  tryCatch({
    # Print LR plots
    for (plot_name in names(lr_plots)) {
      print(lr_plots[[plot_name]])
    }
    
    # Print ratio plots
    for (plot_name in names(ratio_plots)) {
      print(ratio_plots[[plot_name]])
    }
    
    # Print heatmap plots
    for (plot_name in names(heatmap_plots)) {
      print(heatmap_plots[[plot_name]])
    }
    
    log_message(paste("Combined PDF saved as:", pdf_filename))
    
  }, finally = {
    dev.off()
  })
}

# Main execution
log_message("Starting individual locus plotting process...")

tryCatch({
  # Load raw genotype data
  raw_data <- load_raw_genotype_data()
  
  # Prepare data for plotting
  prepared_data <- prepare_locus_data(raw_data)
  long_data <- prepared_data$long_data
  ordered_loci <- prepared_data$ordered_loci
  
  # Create individual plots and store objects
  lr_plots <- create_individual_locus_plots(long_data, ordered_loci)
  ratio_plots <- create_locus_ratio_boxplots(long_data, ordered_loci)
  
  # Create heatmap plots
  heatmap_plots <- create_locus_heatmap_plots(long_data, ordered_loci)
  
  # Create comparison plot
  create_comparison_plot(long_data, ordered_loci)
  
  # Save all plots to combined PDF
  save_plots_to_pdf(lr_plots, ratio_plots, heatmap_plots)
  
  log_message("All individual locus plots created successfully")
  
}, error = function(e) {
  log_message(paste("Critical error in main execution:", e$message))
  log_message("Script execution failed. Please check the logs and data files.")
  stop(paste("Critical error:", e$message))
})

log_message("Individual locus plotting complete!")