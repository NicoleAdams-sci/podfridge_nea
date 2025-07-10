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
  stop("Usage: Rscript code/plots_known_vs_tested_NEA.R <input_dir> [output_dir]")
}

input_dir <- args[1]
log_message(paste("Input directory:", input_dir))

# If output directory is provided, use it, otherwise create one based on the input directory
if (length(args) >= 2) {
  output_subdir <- args[2]
  output_dir <- file.path("output", output_subdir)
} else {
  timestamp <- format(Sys.time(), "%Y%m%d")
  output_dir <- file.path("output", paste0("mismatched_pop_plots_", timestamp))
}
log_message(paste("Output directory:", output_dir))

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define populations
populations <- c("AfAm", "Cauc", "Hispanic", "Asian")

# Define relationship types and ensure correct ordering for plotting
relationship_order <- c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated")
relationship_labels <- c("Parent-Child", "Full Siblings", "Half Siblings", "Cousins", "Second Cousins", "Unrelated")

# Define loci sets and ensure correct ordering for plotting
loci_set_order <- c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29")
loci_set_labels <- c("Core 13", "Identifiler 15", "Expanded 20", "Supplementary", "Autosomal 29")

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

# Load data from each population - BUT FILTER ONLY FOR MISMATCHED POPULATIONS
load_combined_data <- function() {
  log_message("Loading combined LR data for all populations...")
  
  # Initialize empty data frames
  all_combined_lrs <- data.frame()
  all_ratio_stats <- data.frame()
  
  # Track how many files were successfully loaded
  files_loaded <- 0
  
  for (pop in populations) {
    log_message(paste("Processing population:", pop))
    
    # Define the population-specific directory
    pop_dir <- file.path(input_dir, paste0(pop, "_summary"))
    
    if (!dir.exists(pop_dir)) {
      log_message(paste("Directory not found:", pop_dir))
      # Try alternative directory pattern without the "_summary" suffix
      pop_dir <- file.path(input_dir, paste0("summary_", pop, "_"))
      if (!dir.exists(pop_dir)) {
        # Try with just the date pattern
        pop_dir_pattern <- paste0("summary_", pop, "_")
        potential_dirs <- list.dirs(input_dir, recursive = FALSE)
        matching_dirs <- potential_dirs[grepl(pop_dir_pattern, basename(potential_dirs))]
        
        if (length(matching_dirs) > 0) {
          pop_dir <- matching_dirs[1]  # Take the first matching directory
          log_message(paste("Found directory using pattern match:", pop_dir))
        } else {
          log_message(paste("No directory found for population:", pop))
          next
        }
      } else {
        log_message(paste("Found directory using alternative pattern:", pop_dir))
      }
    } else {
      log_message(paste("Found directory:", pop_dir))
    }
    
    # Load combined LRs
    combined_lrs_file <- file.path(pop_dir, paste0("sim_summary_genotypes_", pop, ".csv"))
    data <- safe_read_csv(combined_lrs_file)
    
    if (!is.null(data)) {
      # Print a sample of the data to debug
      log_message(paste("Sample columns in combined LRs data:", paste(names(data)[1:min(5, length(names(data)))], collapse=", ")))
      log_message(paste("Unique base_loci_set values:", paste(unique(data$base_loci_set), collapse=", ")))
      
      # We'll now keep both matched and mismatched population data to enable ratio calculations
      # Store a copy of just the mismatched data for plotting functions that need only mismatched data
      mismatched_data <- data %>% filter(is_correct_pop == FALSE)
      
      # Check if we have any mismatched data to report
      if(nrow(mismatched_data) > 0) {
        log_message(paste("Found", nrow(mismatched_data), "rows with mismatched populations out of", nrow(data), "total rows"))
      } else {
        log_message("WARNING: No mismatched population data found in this file")
      }
      
      # Check if we have data for all loci sets
      if (length(unique(data$base_loci_set)) < length(loci_set_order)) {
        log_message(paste("WARNING: Found only", length(unique(data$base_loci_set)), 
                          "loci sets out of", length(loci_set_order), "expected sets"))
        log_message(paste("Missing loci sets:", 
                          paste(setdiff(loci_set_order, unique(data$base_loci_set)), collapse=", ")))
      }
      
      if (nrow(data) > 0) {
        files_loaded <- files_loaded + 1
        all_combined_lrs <- bind_rows(all_combined_lrs, data)
      }
    }
    
    # Load ratio stats if available
    ratio_stats_file <- file.path(pop_dir, paste0("sim_lr_ratio_summary_", pop, ".csv"))
    ratio_data <- safe_read_csv(ratio_stats_file)
    
    if (!is.null(ratio_data)) {
      all_ratio_stats <- bind_rows(all_ratio_stats, ratio_data)
    }
  }
  
  # Check if we have data
  if (nrow(all_combined_lrs) == 0) {
    log_message("ERROR: No mismatched population LR data was loaded. Please check the input directory and file patterns.")
    stop("No data loaded. Check log for details.")
  }
  
  log_message(paste("Successfully loaded data from", files_loaded, "files"))
  log_message(paste("Total rows in combined data:", nrow(all_combined_lrs)))
  log_message(paste("Unique loci sets found:", paste(unique(all_combined_lrs$base_loci_set), collapse=", ")))
  
  # Set factor levels for proper ordering in plots
  log_message("Setting factor levels...")
  all_combined_lrs$relationship_type <- factor(
    all_combined_lrs$relationship_type, 
    levels = relationship_order,
    labels = relationship_labels
  )
  
  # Handle whatever column name is used for loci sets
  if ("base_loci_set" %in% names(all_combined_lrs)) {
    loci_set_column <- "base_loci_set"
  } else if ("loci_set" %in% names(all_combined_lrs)) {
    loci_set_column <- "loci_set"
  } else {
    log_message("WARNING: Could not find loci set column. Expected 'base_loci_set' or 'loci_set'")
    loci_set_column <- NULL
  }
  
  if (!is.null(loci_set_column)) {
    all_combined_lrs[[loci_set_column]] <- factor(
      all_combined_lrs[[loci_set_column]], 
      levels = loci_set_order,
      labels = loci_set_labels
    )
  }
  
  all_combined_lrs$population <- factor(
    all_combined_lrs$population, 
    levels = populations
  )
  
  all_combined_lrs$freq_population <- factor(
    all_combined_lrs$freq_population, 
    levels = populations
  )
  
  if (nrow(all_ratio_stats) > 0) {
    all_ratio_stats$relationship_type <- factor(
      all_ratio_stats$relationship_type, 
      levels = relationship_order,
      labels = relationship_labels
    )
    
    # Same handling for ratio stats data
    if ("base_loci_set" %in% names(all_ratio_stats)) {
      all_ratio_stats$base_loci_set <- factor(
        all_ratio_stats$base_loci_set, 
        levels = loci_set_order,
        labels = loci_set_labels
      )
    } else if ("loci_set" %in% names(all_ratio_stats)) {
      all_ratio_stats$loci_set <- factor(
        all_ratio_stats$loci_set, 
        levels = loci_set_order,
        labels = loci_set_labels
      )
    }
    
    all_ratio_stats$population <- factor(
      all_ratio_stats$population, 
      levels = populations
    )
    
    all_ratio_stats$freq_population <- factor(
      all_ratio_stats$freq_population, 
      levels = populations
    )
  }
  
  # Get tallies of populations and relationship_types
  detailed_tallies <- list(
    # Unique samples (what probably want for reporting sample sizes)
    unique_samples = all_combined_lrs %>% 
      filter(is_correct_pop == TRUE) %>%
      group_by(population, relationship_type) %>% 
      summarise(n_unique_samples = n_distinct(sim_id), .groups = 'drop'),
    
    # Total LR calculations (includes all population frequency combinations)
    total_lr_calculations = all_combined_lrs %>% 
      group_by(population, relationship_type, is_correct_pop) %>% 
      tally() %>%
      pivot_wider(names_from = is_correct_pop, values_from = n, 
                  names_prefix = "lr_calculations_correct_pop_"),
    
    # Summary by loci set (if you want to see how many samples per loci set)
    by_loci_set = all_combined_lrs %>% 
      filter(is_correct_pop == TRUE) %>%
      group_by(population, relationship_type, base_loci_set) %>% 
      summarise(n_unique_samples = n_distinct(sim_id), .groups = 'drop')
  )
  
  return(list(
    combined_lrs = all_combined_lrs,
    ratio_stats = all_ratio_stats,
    population_relationship_tallies = detailed_tallies$unique_samples,
    detailed_tallies = detailed_tallies
  ))
}

# Function to create boxplot of LR distributions comparing frequency sources
plot_lr_distributions_by_freq_source <- function(combined_lrs, population_relationship_tallies) {
  log_message("Creating boxplot of LR distributions by frequency source...")
  
  # Determine which column contains the loci set information
  if ("base_loci_set" %in% names(combined_lrs)) {
    loci_set_column <- "base_loci_set"
  } else if ("loci_set" %in% names(combined_lrs)) {
    loci_set_column <- "loci_set"
  } else {
    log_message("ERROR: Could not find loci set column in data. Expected 'base_loci_set' or 'loci_set'")
    return(NULL)
  }
  
  log_message(paste("Using", loci_set_column, "as loci set column"))
  log_message(paste("Unique values in loci set column:", 
                    paste(levels(combined_lrs[[loci_set_column]]), collapse=", ")))
  
  # Calculate total samples per population from the tallies
  population_totals <- population_relationship_tallies %>%
    group_by(population) %>%
    summarize(total_n = sum(n_unique_samples), .groups = 'drop')
  
  # Create legend labels with population names and total sample counts
  legend_labels <- setNames(
    paste0(population_totals$population, " (n=", format(population_totals$total_n, big.mark = ","), ")"),
    population_totals$population
  )
  
  # Define a color palette for frequency populations
  freq_pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(freq_pop_colors) <- levels(combined_lrs$freq_population)
  
  tryCatch({
    # Create the formula for facet_grid based on the loci set column
    facet_formula <- as.formula(paste("population ~", loci_set_column))
    
    p <- ggplot(combined_lrs, aes(x = relationship_type, y = LR, fill = freq_population)) +
      geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7, outlier.size = 0.5) +
      facet_grid(facet_formula, scales = "fixed") +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_fill_manual(values = freq_pop_colors, labels = legend_labels) +
      labs(
        title = "LR Distributions Using Mismatched Population Allele Frequencies",
        subtitle = "True Population (rows) vs. Frequency Source Population (colors)",
        x = "Relationship Type",
        y = "LR (log scale)",
        fill = "Frequency Population"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12, margin = margin(t = 10)),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.box.margin = margin(t = 10),
        strip.text = element_text(size = 12)
      )
    
    ggsave(file.path(output_dir, "mismatched_pop_lr_boxplot.png"), plot = p, width = 16, height = 12, bg = "white", dpi = 300)
    log_message("LR distributions plot saved to mismatched_pop_lr_boxplot.png")
    
    return(p)
  }, error = function(e) {
    log_message(paste("Error creating boxplot:", e$message))
    return(NULL)
  })
}

# Function to create violin plots showing distribution of LRs by frequency source
plot_violin_lr_by_freq_source <- function(combined_lrs, population_relationship_tallies) {
  log_message("Creating violin plot of LR distributions by frequency source...")
  
  # Determine which column contains the loci set information
  if ("base_loci_set" %in% names(combined_lrs)) {
    loci_set_column <- "base_loci_set"
  } else if ("loci_set" %in% names(combined_lrs)) {
    loci_set_column <- "loci_set"
  } else {
    log_message("ERROR: Could not find loci set column in data. Expected 'base_loci_set' or 'loci_set'")
    return(NULL)
  }
  
  log_message(paste("Using", loci_set_column, "as loci set column for violin plot"))
  
  # Calculate total samples per population from the tallies
  population_totals <- population_relationship_tallies %>%
    group_by(population) %>%
    summarize(total_n = sum(n_unique_samples), .groups = 'drop')
  
  # Create legend labels with population names and total sample counts
  legend_labels <- setNames(
    paste0(population_totals$population, " (n=", format(population_totals$total_n, big.mark = ","), ")"),
    population_totals$population
  )
  
  # Define a color palette for frequency populations
  freq_pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(freq_pop_colors) <- levels(combined_lrs$freq_population)
  
  tryCatch({
    # Create the formula for facet_grid based on the loci set column
    facet_formula <- as.formula(paste("population ~", loci_set_column))
    
    p <- ggplot(combined_lrs, aes(x = relationship_type, y = LR, fill = freq_population)) +
      geom_violin(position = position_dodge(width = 0.85), alpha = 0.7, scale = "width", trim = TRUE) +
      facet_grid(facet_formula, scales = "fixed") +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_fill_manual(values = freq_pop_colors, labels = legend_labels) +
      labs(
        title = "LR Distribution Densities Using Mismatched Population Allele Frequencies",
        subtitle = "True Population (rows) vs. Frequency Source Population (colors)",
        x = "Relationship Type",
        y = "LR (log scale)",
        fill = "Frequency Population"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.text = element_text(size = 12)
      )
    
    ggsave(file.path(output_dir, "mismatched_pop_violin_lr.png"), plot = p, width = 16, height = 12, bg = "white", dpi = 300)
    log_message("Violin plot saved to mismatched_pop_violin_lr.png")
    
    return(p)
  }, error = function(e) {
    log_message(paste("Error creating violin plot:", e$message))
    return(NULL)
  })
}

# Function to plot ratio of wrong vs. correct population LRs
plot_ratio_heatmap <- function(ratio_stats) {
  log_message("Creating heatmap of LR ratios (wrong/correct)...")
  
  if (is.null(ratio_stats) || nrow(ratio_stats) == 0) {
    log_message("No ratio stats data available for heatmap")
    return(NULL)
  }
  
  # Determine which column contains the loci set information
  if ("base_loci_set" %in% names(ratio_stats)) {
    loci_set_column <- "base_loci_set"
  } else if ("loci_set" %in% names(ratio_stats)) {
    loci_set_column <- "loci_set"
  } else {
    log_message("ERROR: Could not find loci set column in ratio stats data. Expected 'base_loci_set' or 'loci_set'")
    return(NULL)
  }
  
  log_message(paste("Using", loci_set_column, "as loci set column for ratio heatmap"))
  
  tryCatch({
    # Filter out NaN and Inf values
    ratio_statsA <- ratio_stats %>% 
      filter(is.finite(mean_ratio))
    
    ratio_stats <- ratio_statsA %>%
      complete(population, freq_population, relationship_type, !!sym(loci_set_column))
    
    # --- Define colors and values for the gradient ---
    my_colors_full_palette <- c("blue", "green", "white", "yellow", "red")
    
    # Calculate range on log10 scale from the filtered data
    log_transformed_ratios <- log10(ratio_stats$mean_ratio)
    log_min_val <- min(log_transformed_ratios, na.rm = TRUE)
    log_max_val <- max(log_transformed_ratios, na.rm = TRUE)
    
    # Determine final colors and values for scale_fill_gradientn
    if (log_min_val == log_max_val) {
      # Only one unique value in the data after filtering and log-transform
      # Pick a single color based on where this value lies relative to log10(1)=0
      if (log_min_val == 0) {
        final_colors <- "white"
      } else if (log_min_val < 0) {
        # Could be blue or green; green is a bit softer if it's not extremely low
        final_colors <- "green" 
      } else { # log_min_val > 0
        final_colors <- "yellow"
      }
      final_values <- NULL # Not needed for a single color / default ggplot behavior
    } else {
      # Calculate the normalized position of "white" (log10(1) = 0)
      # within the log-transformed data range [log_min_val, log_max_val]
      val_for_white_rc <- scales::rescale(0, from = c(log_min_val, log_max_val), to = c(0,1))
      
      if (val_for_white_rc <= 0) { 
        # Center (0) is at or below the minimum of the data's log range
        # (e.g., all data ratios are > 1, so log values are all positive)
        # Scale should go from White (at data min) -> Yellow -> Red (at data max)
        final_colors <- c("white", "yellow", "red")
        final_values <- c(0, 0.5, 1) 
      } else if (val_for_white_rc >= 1) {
        # Center (0) is at or above the maximum of the data's log range
        # (e.g., all data ratios are < 1, so log values are all negative)
        # Scale should go from Blue (at data min) -> Green -> White (at data max)
        final_colors <- c("blue", "green", "white")
        final_values <- c(0, 0.5, 1)
      } else {
        # Center (0) is within the data's log range. Use the full B-G-W-Y-R.
        # White is at val_for_white_rc.
        # Green is halfway between data_min_rescaled (0) and val_for_white_rc.
        # Yellow is halfway between val_for_white_rc and data_max_rescaled (1).
        final_colors <- my_colors_full_palette
        final_values <- c(0,                            # Blue at data_min_rescaled
                          val_for_white_rc * 0.5,       # Green
                          val_for_white_rc,             # White
                          val_for_white_rc + (1 - val_for_white_rc) * 0.5, # Yellow
                          1)                            # Red at data_max_rescaled
        # Ensure final_values are unique and sorted (they should be by this logic if 0 < val_for_white_rc < 1)
        # To be perfectly safe if val_for_white_rc is extremely close to 0 or 1:
        final_values <- sort(unique(final_values))
        if(length(final_values) < length(final_colors)){
          # This case implies some points collapsed, e.g. val_for_white_rc was ~0 or ~1
          # For simplicity, we'll let scale_fill_gradientn handle it if values are not strictly increasing
          # by only using the first color if multiple colors map to the same value.
          # The above if/else if/else structure *should* prevent this issue mostly.
          # A more robust way if points collapse is to also filter final_colors, but it gets complex.
          # The current logic for the three cases (val_for_white_rc <=0, >=1, or in between) simplifies this.
        }
      }
    }
    
    # Create the formula for facet_grid based on the loci set column
    facet_formula <- as.formula(paste("population ~", loci_set_column))
    
    p <- ggplot(ratio_stats, aes(x = freq_population, y = relationship_type, fill = mean_ratio)) +
      geom_tile(color = "grey80", size = 0.1) +
      facet_grid(facet_formula) +
      scale_fill_gradientn(
        colors = final_colors,
        values = final_values, # Normalized positions for the colors 
        trans = "log10",
        breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
        labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
        limits = c(min(ratio_stats$mean_ratio), max(ratio_stats$mean_ratio)), # Use actual data limits
        na.value = "grey75",
        guide = guide_colorbar(
          barwidth = 20,
          barheight = 1.2,
          title.position = "top",
          title.hjust = 0.5,
          label.theme = element_text(size = 9)
        )
      ) +
      labs(
        #title = "Mean Ratio of Mismatched Population LR to Correct Population LR",
        #subtitle = "True Population (rows) vs. Wrong Population Frequency Source (columns)",
        #x = "Wrong Population Frequency Source",
        #y = "Relationship Type",
        #fill = "Wrong LR / Correct LR\n(log scale)",
        x = "Assumed Ancestry (Frequency Source)", # Relabelled for clarity
        y = "True Relationship Type",               # Relabelled
        fill = "Ratio: LR(Tested Pop) / LR(True Pop) \n(log10 scale, 1 = white = no change)" # More descriptive
      ) +
      theme_minimal(base_size = 11) +
      theme(
        #plot.title = element_text(hjust = 0.5, size = 16),
        #plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 12, face = "bold", lineheight = 1.2),
        legend.text = element_text(size = 9),
        legend.position = "bottom",
        #strip.text = element_text(size = 9)
        strip.text = element_text(size = 10, face = "bold", margin = margin(t=3,b=3)), # Facet labels
        panel.spacing = unit(1, "lines"), # Add a bit more space between facets
        panel.grid.major = element_blank(), # Remove major grid lines for cleaner tiles
        panel.grid.minor = element_blank()  # Remove minor grid lines
      )
    
    ggsave(file.path(output_dir, "mismatched_pop_ratio_heatmap.png"), plot = p, width = 16, height = 14, bg="white", dpi = 300)
    log_message("Ratio heatmap saved to mismatched_pop_ratio_heatmap.png")
    
    return(p)
  }, error = function(e) {
    log_message(paste("Error creating ratio heatmap:", e$message))
    return(NULL)
  })
}

# Function to plot ratio distributions as boxplots
plot_ratio_boxplots <- function(combined_lrs, population_relationship_tallies) {
  log_message("Creating ratio distribution boxplots...")
  
  # We need to prepare the data by:
  # 1. Getting the correct population LR for each case
  # 2. Calculating the ratio for each wrong population LR
  
  # Determine which column contains the loci set information
  if ("base_loci_set" %in% names(combined_lrs)) {
    loci_set_column <- "base_loci_set"
  } else if ("loci_set" %in% names(combined_lrs)) {
    loci_set_column <- "loci_set"
  } else {
    log_message("ERROR: Could not find loci set column in data. Expected 'base_loci_set' or 'loci_set'")
    return(NULL)
  }
  
  log_message(paste("Using", loci_set_column, "as loci set column for ratio boxplots"))
  
  tryCatch({
    # Separate the combined LRs into correct and wrong
    log_message("Separating correct and wrong population data for ratio calculation...")
    
    # Check if we have both correct and wrong population data
    if (sum(combined_lrs$is_correct_pop == TRUE) == 0) {
      log_message("ERROR: No correct population data available. Cannot calculate ratios.")
      return(NULL)
    }
    
    if (sum(combined_lrs$is_correct_pop == FALSE) == 0) {
      log_message("ERROR: No wrong population data available. Cannot calculate ratios.")
      return(NULL)
    }
    
    correct_lrs <- combined_lrs %>%
      filter(is_correct_pop == TRUE) %>%
      select(population, relationship_type, sim_id, all_of(loci_set_column), correct_LR = LR)
    
    wrong_lrs <- combined_lrs %>%
      filter(is_correct_pop == FALSE) %>%
      select(population, relationship_type, sim_id, all_of(loci_set_column), freq_population, wrong_LR = LR)
    
    # Merge correct and wrong to calculate ratios
    ratio_data <- wrong_lrs %>%
      left_join(
        correct_lrs,
        by = c("population", "relationship_type", "sim_id", loci_set_column)
      ) %>%
      mutate(ratio = wrong_LR / correct_LR,
      # Calculate percentage change: ((new - old) / old) * 100
      percent_change = ((wrong_LR - correct_LR) / correct_LR) * 100) %>%
      filter(is.finite(ratio), is.finite(percent_change))
    
    if (nrow(ratio_data) == 0) {
      log_message("No valid ratio data available for boxplot")
      return(NULL)
    }
    
    # Calculate total samples per population from the tallies
    population_totals <- population_relationship_tallies %>%
      group_by(population) %>%
      summarize(total_n = sum(n_unique_samples), .groups = 'drop')
    
    # Create legend labels with population names and total sample counts
    legend_labels <- setNames(
      paste0(population_totals$population, " (n=", format(population_totals$total_n, big.mark = ","), ")"),
      population_totals$population
    )
    
    # Define a color palette for frequency populations
    freq_pop_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
    names(freq_pop_colors) <- levels(ratio_data$freq_population)
    
    # Create the formula for facet_grid based on the loci set column
    facet_formula <- as.formula(paste("population ~", loci_set_column))
    
    p1 <- ggplot(ratio_data, aes(x = relationship_type, y = ratio, fill = freq_population)) +
      geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7, outlier.size = 0.5) +
      facet_grid(facet_formula, scales = "fixed") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "goldenrod") +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_fill_manual(values = freq_pop_colors, labels = legend_labels) +
      labs(
        title = "Ratio of Mismatched Population LR to Correct Population LR",
        subtitle = "True Population (rows) vs. Wrong Population Frequency Source (colors)",
        x = "Relationship Type",
        y = "Wrong LR / Correct LR (log scale)",
        fill = "Wrong Population"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.text = element_text(size = 12)
      )
    
    ggsave(file.path(output_dir, "mismatched_pop_ratio_boxplot.png"), plot = p1, width = 16, height = 12,bg = "white", dpi = 300)
    log_message("Ratio boxplot saved to mismatched_pop_ratio_boxplot.png")
    
    # Percent change plot with fixed y scale - squishes all data don't use for final
    p2 <- ggplot(ratio_data, aes(x = relationship_type, y = percent_change, fill = freq_population)) +
      geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7, outlier.size = 0.5) +
      facet_grid(facet_formula, scales = "fixed") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "goldenrod", size = 0.8) +
      scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        breaks = scales::pretty_breaks(n = 6)
      ) +
      scale_fill_manual(values = freq_pop_colors, labels = legend_labels) +
      labs(
        title = "Percentage Change in LR When Using Wrong Population Frequencies",
        subtitle = paste0("True Population (rows) vs. Wrong Population Frequency Source (colors)\n",
                          "Dashed line = no change (0%). "),
        x = "Relationship Type",
        y = "Percentage Change in LR",
        fill = "Wrong Population"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.text = element_text(size = 12)
      )
    ggsave(file.path(output_dir, "mismatched_pop_percent_change_boxplot.png"), 
           plot = p2, width = 16, height = 12, bg = "white", dpi = 300)
    log_message("Percentage change boxplot saved to mismatched_pop_percent_change_boxplot.png")
    
    # Percent change plot with FREE axes -- this isn't working due to facet formula I think
    # p2a <- ggplot(ratio_data, aes(x = relationship_type, y = percent_change, fill = freq_population)) +
    #   geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7, outlier.size = 0.5) +
    #   facet_grid(facet_formula, scales = "free") +
    #   #geom_hline(yintercept = 0, linetype = "dashed", color = "goldenrod", size = 0.8) +
    #   scale_y_continuous(
    #     labels = function(x) paste0(x, "%"),
    #     breaks = scales::pretty_breaks(n = 6)
    #   ) +
    #   scale_fill_manual(values = freq_pop_colors, labels = legend_labels) +
    #   labs(
    #     title = "Percentage Change in LR When Using Wrong Population Frequencies",
    #     subtitle = paste0("True Population (rows) vs. Wrong Population Frequency Source (colors)\n",
    #                       "Dashed line = no change (0%). "),
    #     x = "Relationship Type",
    #     y = "Percentage Change in LR",
    #     fill = "Wrong Population"
    #   ) +
    #   theme_minimal() +
    #   theme(
    #     plot.title = element_text(hjust = 0.5, size = 16),
    #     plot.subtitle = element_text(hjust = 0.5, size = 14),
    #     axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    #     axis.title = element_text(size = 14),
    #     legend.title = element_text(size = 12),
    #     legend.text = element_text(size = 10),
    #     legend.position = "bottom",
    #     strip.text = element_text(size = 12)
    #   )
    # ggsave(file.path(output_dir, "mismatched_pop_percent_change_boxplot_free.png"), 
    #        plot = p2a, width = 16, height = 12, bg = "white", dpi = 300)
    # log_message("Percentage change boxplot saved to mismatched_pop_percent_change_boxplot_free.png")
    # 
    # Separate plot focusing on the main distribution of percent change (zoomed in)
    # Calculate reasonable zoom range based on IQR
    q25 <- quantile(ratio_data$percent_change, 0.25, na.rm = TRUE)
    q75 <- quantile(ratio_data$percent_change, 0.75, na.rm = TRUE)
    iqr <- q75 - q25
    
    # Use 1.5*IQR rule for "typical" outlier detection, but extend a bit more
    zoom_min <- q25 - 3 * iqr
    zoom_max <- q75 + 3 * iqr
    
    # Count points outside zoom range
    zoom_stats <- ratio_data %>%
      summarise(
        pct_outside_zoom = mean(percent_change < zoom_min | percent_change > zoom_max, na.rm = TRUE) * 100,
        n_outside = sum(percent_change < zoom_min | percent_change > zoom_max, na.rm = TRUE),
        total_n = n()
      )
    
    p2b <- ggplot(ratio_data, aes(x = relationship_type, y = percent_change, fill = freq_population)) +
      geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.7, outlier.size = 0.5) +
      facet_grid(facet_formula, scales = "fixed") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
      coord_cartesian(ylim = c(zoom_min, zoom_max)) +
      scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        breaks = scales::pretty_breaks(n = 8)
      ) +
      scale_fill_manual(values = freq_pop_colors, labels = legend_labels) +
      labs(
        title = "Percentage Change in LR - Main Distribution (Zoomed View)",
        subtitle = sprintf("Focused view excluding extreme outliers (%d points, %.1f%% of data, outside range)\nSee companion plot for full data range", 
                           zoom_stats$n_outside, zoom_stats$pct_outside_zoom),
        x = "Relationship Type",
        y = "Percentage Change in LR",
        fill = "Wrong Population"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.text = element_text(size = 12)
      )
    ggsave(file.path(output_dir, "mismatched_pop_percent_change_boxplot_zoom.png"), 
           plot = p2b, width = 16, height = 12, bg = "white", dpi = 300)
    log_message("Percentage change boxplot saved to mismatched_pop_percent_change_boxplot_zoom.png")
    
    
    return(list(
      ratio_plot = p1, 
      percent_change_plot = p2,
      #percent_change_freey = p2a,
      percent_change_zoomed = p2b,
    ))
  }, error = function(e) {
    log_message(paste("Error creating ratio boxplots:", e$message))
    return(NULL)
  })
}

# Main execution
log_message("Starting mismatched population plotting process...")

# Load all data
tryCatch({
  all_data <- load_combined_data()
  
  # Create and save plots
  if (!is.null(all_data$combined_lrs) && nrow(all_data$combined_lrs) > 0) {
    # Create a filtered dataset that includes only mismatched populations for the distribution plots
    mismatched_lrs <- all_data$combined_lrs %>% filter(is_correct_pop == FALSE)
    
    if(nrow(mismatched_lrs) > 0) {
      # Plot LR distributions by frequency source using only mismatched data
      log_message("Creating distribution plots using only mismatched population data...")
      lr_boxplot <- plot_lr_distributions_by_freq_source(mismatched_lrs, all_data$population_relationship_tallies)
      violin_plot <- plot_violin_lr_by_freq_source(mismatched_lrs, all_data$population_relationship_tallies)
    } else {
      log_message("No mismatched population data available for distribution plots")
    }
    
    # Plot ratio distributions using the complete dataset
    ratio_boxplots <- plot_ratio_boxplots(all_data$combined_lrs, all_data$population_relationship_tallies)
  } else {
    log_message("No combined_lrs data available for plotting")
  }
  
  if (!is.null(all_data$ratio_stats) && nrow(all_data$ratio_stats) > 0) {
    # Plot ratio heatmap
    ratio_heatmap <- plot_ratio_heatmap(all_data$ratio_stats)
  } else {
    log_message("No ratio_stats data available for plotting")
  }
  
  log_message("All available plots created successfully")
}, error = function(e) {
  log_message(paste("Critical error in main execution:", e$message))
  log_message("Script execution failed. Please check the logs and data files.")
  stop(paste("Critical error:", e$message))
})

log_message("Mismatched population plotting complete!")