# AMR RESISTOME LATENT PROFILE ANALYSIS (LPA)
# Multi-Country Analysis with 7 Antibiotics
# Description: Latent Profile Analysis of AMR resistome patterns across 
#              multiple countries using selected antibiotic MIC data

# SETUP & LIBRARIES

# Install required packages if missing
required_packages <- c("mclust", "boot", "readxl", "dplyr", "tidyr", "stringr",
                       "ggplot2", "scales", "gridExtra", "knitr", "forcats",
                       "viridis", "ggpubr", "patchwork", "psych", "corrplot",
                       "reshape2", "grid", "cowplot", "RColorBrewer", "openxlsx",
                       "kableExtra", "pheatmap", "ComplexHeatmap", "circlize",
                       "ggridges")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("ComplexHeatmap", "circlize")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

suppressPackageStartupMessages({
  library(mclust)
  library(boot)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(gridExtra)
  library(knitr)
  library(forcats)
  library(viridis)
  library(ggpubr)
  library(patchwork)
  library(psych)
  library(corrplot)
  library(reshape2)
  library(grid)
  library(cowplot)
  library(RColorBrewer)
  library(openxlsx)
  library(kableExtra)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(ggridges)
})

cat("\nAMR RESISTOME LATENT PROFILE ANALYSIS - MULTI-COUNTRY")
cat("\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# CONFIGURATION

# Analysis parameters
MAX_CLASSES <- 5
MIN_ISOLATES_PER_COUNTRY <- 30
MIN_YEAR_COVERAGE <- 0.60
MIN_DATA_COMPLETE <- 0.70
MIC_GE_MULTIPLIER <- 2
MCLUST_MODEL_NAMES <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
RANDOM_SEED <- 123

# Selected antibiotics (7 antibiotics from the data)
SELECTED_ANTIBIOTICS <- c(
  "Amikacin",
  "Amoxycillin clavulanate",
  "Ampicillin",
  "Cefepime",
  "Levofloxacin",
  "Piperacillin tazobactam",
  "Tigecycline"
)

# Country panels for visualization
COUNTRY_PANELS <- list(
  panel_a = c("Costa Rica", "Latvia", "Qatar"),
  panel_b = c("Argentina", "Brazil", "Kuwait", "Russia"),
  panel_c = c("Romania", "Ukraine", "Spain", "Switzerland", "South Africa")
)

# Color palettes
RDBU_PALETTE <- rev(brewer.pal(11, "RdBu"))
CLASS_COLORS <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00")
names(CLASS_COLORS) <- paste0("Class", 1:5)
ICU_COLORS <- c("ICU" = "#D55E00", "Non-ICU" = "#0072B2")
AGE_COLORS <- c("Adult" = "#2E86AB", "Elderly" = "#A23B72", 
                "Pediatric" = "#F18F01", "Unknown" = "#6C757D")

# CORE FUNCTIONS

# Create directory structure for outputs
create_dirs <- function(base_dir) {
  subdirs <- c("tables", "figures/heatmaps", "figures/trajectories",
               "figures/icu", "figures/age", "figures/posterior",
               "models", "reports", "diagnostics")
  for (d in subdirs) {
    dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
  }
}

# Save table to CSV or Excel
save_table <- function(data, filename, output_dir, format = "csv") {
  path <- file.path(output_dir, "tables", filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (format == "csv") write.csv(data, paste0(path, ".csv"), row.names = FALSE)
  else if (format == "xlsx") openxlsx::write.xlsx(data, paste0(path, ".xlsx"))
  message(sprintf("  Table saved: %s", path))
}

# Save figure as PNG and PDF
save_figure <- function(plot, filename, output_dir, subdir = "", 
                        width = 12, height = 10, dpi = 300) {
  path <- file.path(output_dir, "figures", subdir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(path, ".png"), plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(path, ".pdf"), plot, width = width, height = height, dpi = dpi)
  message(sprintf("  Figure saved: %s", path))
}

# Load Excel data
load_data <- function(path, sheet = 1) {
  as.data.frame(read_excel(path, sheet = sheet))
}

# Flag ICU cases
flag_icu <- function(speciality_vec) {
  grepl("ICU", speciality_vec, ignore.case = TRUE)
}

# Parse MIC values with >= handling
parse_mic <- function(x, ge_multiplier = 2) {
  x <- trimws(as.character(x))
  is_ge <- grepl("^>=?", x)
  num <- suppressWarnings(as.numeric(gsub("[<>=]", "", x)))
  num[is_ge] <- num[is_ge] * ge_multiplier
  num
}

# Build MIC matrix from data
build_mic_matrix <- function(df, ge_multiplier = 2) {
  i_cols <- grep("_I$", names(df), value = TRUE)
  mic_cols <- str_remove(i_cols, "_I$")
  mic_cols <- mic_cols[mic_cols %in% names(df)]
  if (length(mic_cols) == 0) stop("No MIC columns found.")
  mat <- df[, mic_cols, drop = FALSE]
  mat[] <- lapply(mat, parse_mic, ge_multiplier = ge_multiplier)
  mat[] <- lapply(mat, function(x) { x[!is.na(x) & x <= 0] <- NA; log2(x) })
  names(mat) <- make.names(names(mat))
  as.data.frame(mat)
}

# Select core antibiotic panel
select_core_panel <- function(df, min_coverage = 0.60, selected_ab = NULL) {
  years <- unique(df$Year)
  
  year_mats <- lapply(years, function(y) {
    df_yr <- df %>% filter(Year == y)
    build_mic_matrix(df_yr)
  })
  names(year_mats) <- as.character(years)
  
  all_ab <- unique(unlist(lapply(year_mats, names)))
  
  if (!is.null(selected_ab)) {
    selected_ab_clean <- make.names(selected_ab)
    all_ab <- all_ab[all_ab %in% selected_ab_clean]
  }
  
  coverage <- sapply(all_ab, function(ab) {
    mean(sapply(year_mats, function(mat) ab %in% names(mat)))
  })
  
  core_ab <- names(coverage)[coverage >= min_coverage]
  
  valid_ab <- sapply(core_ab, function(ab) {
    all(sapply(year_mats, function(mat) {
      if (!ab %in% names(mat)) return(FALSE)
      sum(!is.na(mat[[ab]])) / nrow(mat) >= MIN_DATA_COMPLETE
    }))
  })
  
  core_panel <- core_ab[valid_ab]
  
  cat(sprintf("\nCore panel: %d antibiotics\n", length(core_panel)))
  print(core_panel)
  
  list(core_panel = core_panel, coverage = coverage[core_panel],
       year_matrices = year_mats, years = years)
}

# Impute missing values with column means
impute_missing <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    col_mean <- mean(mat[, j], na.rm = TRUE)
    if (!is.na(col_mean)) mat[is.na(mat[, j]), j] <- col_mean
  }
  mat
}

# Fit LPA model
fit_lpa <- function(mat, k_range = 2:5, seed = 123, verbose = TRUE) {
  set.seed(seed)
  
  results <- list()
  best_bic <- Inf
  best_fit <- NULL
  best_k <- NA
  
  for (k in k_range) {
    if (verbose) cat(sprintf("  G=%d...", k))
    fit <- tryCatch(
      Mclust(data = as.matrix(mat), G = k, modelNames = MCLUST_MODEL_NAMES, verbose = FALSE),
      error = function(e) { if (verbose) cat(sprintf(" failed: %s\n", e$message)); NULL }
    )
    
    if (!is.null(fit) && length(fit$modelName) > 0) {
      if (verbose) cat(sprintf(" done (BIC=%.1f)\n", -fit$bic))
      results[[as.character(k)]] <- fit
      if (-fit$bic < best_bic) {
        best_bic <- -fit$bic
        best_fit <- fit
        best_k <- k
      }
    } else {
      if (verbose) cat(" failed\n")
    }
  }
  
  if (is.null(best_fit)) return(NULL)
  
  list(all_fits = results, best_fit = best_fit, best_k = best_k, best_bic = best_bic)
}

# Extract class profiles
extract_profile <- function(fit, var_names) {
  means <- t(fit$parameters$mean)
  means <- as.data.frame(means)
  names(means) <- var_names
  rownames(means) <- paste0("Class", seq_len(nrow(means)))
  means
}

# Predict classes for new data
predict_classes <- function(fit, newdata) {
  pred <- predict(fit, newdata = as.matrix(newdata))
  list(classification = pred$classification, uncertainty = pred$uncertainty,
       posterior = pred$z)
}

# ANALYZE COUNTRY

# Analyze a single country
analyze_country <- function(df_all, country_name, output_base, selected_ab) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("ANALYZING:", country_name, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  country_dir <- file.path(output_base, country_name)
  create_dirs(country_dir)
  
  df_country <- df_all %>% filter(Country == country_name)
  df_kp <- df_country %>% filter(Species == "Klebsiella pneumoniae", Source == "Blood")
  
  cat(sprintf("  Isolates: %d\n", nrow(df_kp)))
  
  if (nrow(df_kp) < MIN_ISOLATES_PER_COUNTRY) {
    cat(sprintf("  SKIP: %d isolates (< %d)\n", nrow(df_kp), MIN_ISOLATES_PER_COUNTRY))
    return(NULL)
  }
  
  core_info <- select_core_panel(df_kp, min_coverage = MIN_YEAR_COVERAGE,
                                 selected_ab = selected_ab)
  core_panel <- sort(core_info$core_panel)
  
  cat(sprintf("  Core antibiotics: %d\n", length(core_panel)))
  print(core_panel)
  
  if (length(core_panel) < 3) {
    cat(sprintf("  SKIP: %d antibiotics (< 3)\n", length(core_panel)))
    return(NULL)
  }
  
  mat_all <- build_mic_matrix(df_kp)
  mat_core <- mat_all[, intersect(names(mat_all), core_panel), drop = FALSE]
  mat_imp <- impute_missing(mat_core)
  
  k_max <- min(MAX_CLASSES, floor(nrow(mat_imp) / 30))
  k_range <- 2:max(2, min(k_max, MAX_CLASSES))
  
  ref_model <- fit_lpa(mat_imp, k_range = k_range, seed = RANDOM_SEED)
  
  if (is.null(ref_model)) {
    cat("  SKIP: No model converged\n")
    return(NULL)
  }
  
  cat(sprintf("  Classes: %d (%s)\n", ref_model$best_k, ref_model$best_fit$modelName))
  
  log2prof <- extract_profile(ref_model$best_fit, names(mat_imp))
  log2prof <- log2prof[, sort(colnames(log2prof))]
  
  pred <- predict_classes(ref_model$best_fit, mat_imp)
  
  df_kp$class_assignment <- pred$classification
  df_kp$class_uniform <- factor(pred$classification, levels = 1:MAX_CLASSES)
  df_kp$icu_flag <- flag_icu(df_kp$Speciality)
  
  if ("AgeGroup" %in% names(df_kp)) {
    df_kp$age_group <- df_kp$AgeGroup
  } else {
    df_kp$age_group <- NA
  }
  
  class_by_year <- list()
  for (yr in sort(unique(df_kp$Year))) {
    idx <- which(df_kp$Year == yr)
    if (length(idx) > 0) {
      class_by_year[[as.character(yr)]] <- list(
        year = yr,
        n = length(idx),
        classification = pred$classification[idx],
        posterior = pred$posterior[idx, , drop = FALSE]
      )
    }
  }
  
  return(list(
    country = country_name,
    core_panel = core_panel,
    ref_model = ref_model,
    log2prof = log2prof,
    df = df_kp,
    posterior = pred$posterior,
    class_by_year = class_by_year,
    n_classes = ref_model$best_k,
    max_classes = MAX_CLASSES
  ))
}

# VISUALIZATION FUNCTIONS

# Create heatmap for country
create_heatmap <- function(log2prof, country, output_dir, max_classes = 5, vmax = 1.9334) {
  scaled <- t(scale(t(log2prof)))
  ab_order <- sort(colnames(scaled))
  scaled <- scaled[, ab_order, drop = FALSE]
  colnames(scaled) <- gsub("\\.", " ", colnames(scaled))
  
  n_actual <- nrow(scaled)
  uni_mat <- matrix(NA, nrow = max_classes, ncol = ncol(scaled))
  rownames(uni_mat) <- paste0("Class", 1:max_classes)
  colnames(uni_mat) <- colnames(scaled)
  for (i in 1:min(n_actual, max_classes)) {
    uni_mat[i, ] <- as.numeric(scaled[i, ])
  }
  
  df_heat <- melt(uni_mat, varnames = c("Class", "Antibiotic"), value.name = "value")
  df_heat$Class <- factor(df_heat$Class, levels = paste0("Class", 1:max_classes))
  
  p <- ggplot(df_heat, aes(x = Antibiotic, y = Class, fill = value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradientn(
      colors = RDBU_PALETTE,
      limits = c(-vmax, vmax),
      na.value = "grey95",
      name = "Scaled\nlog2(MIC)",
      guide = guide_colorbar(barwidth = 0.8, barheight = 4)
    ) +
    labs(title = country, x = "", y = "") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
    )
  
  save_figure(p, paste0("heatmap_", country), output_dir, "heatmaps", width = 8, height = 6)
  return(p)
}

# Create trajectory plot for country
create_trajectory <- function(class_by_year, country, output_dir, max_classes = 5) {
  class_dist <- do.call(rbind, lapply(names(class_by_year), function(yr) {
    yr_data <- class_by_year[[yr]]
    counts <- table(factor(yr_data$classification, levels = 1:max_classes))
    data.frame(
      year = as.integer(yr),
      class = paste0("Class", 1:max_classes),
      n = as.integer(counts),
      pct = as.numeric(counts) / yr_data$n * 100
    )
  }))
  
  class_dist$class <- factor(class_dist$class, levels = paste0("Class", 1:max_classes))
  
  p <- ggplot(class_dist, aes(x = year, y = pct, color = class, group = class)) +
    geom_line(size = 1.2) +
    geom_point(size = 3, shape = 21, fill = "white", stroke = 0.5) +
    scale_color_manual(values = CLASS_COLORS, name = "Class", drop = FALSE) +
    scale_x_continuous(breaks = unique(class_dist$year)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, max(class_dist$pct, na.rm = TRUE) * 1.1)) +
    labs(title = country, x = "Year", y = "Proportion (%)") +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
    ) +
    guides(color = guide_legend(nrow = 1))
  
  save_figure(p, paste0("trajectory_", country), output_dir, "trajectories", width = 7, height = 5)
  return(p)
}

# Create panel figure
create_panel_figure <- function(country_list, panel_label, all_results, output_dir,
                                type = "heatmap", vmax = 1.9334) {
  cat(sprintf("\nCreating Panel %s\n", panel_label))
  
  avail <- country_list[country_list %in% names(all_results)]
  if (length(avail) == 0) {
    cat("  No countries available\n")
    return(NULL)
  }
  
  all_data <- data.frame()
  
  if (type == "heatmap") {
    for (country in avail) {
      res <- all_results[[country]]
      if (is.null(res$log2prof)) next
      
      scaled <- t(scale(t(res$log2prof)))
      ab_order <- sort(colnames(scaled))
      scaled <- scaled[, ab_order, drop = FALSE]
      colnames(scaled) <- gsub("\\.", " ", colnames(scaled))
      
      uni_mat <- matrix(NA, nrow = 5, ncol = ncol(scaled))
      rownames(uni_mat) <- paste0("Class", 1:5)
      colnames(uni_mat) <- colnames(scaled)
      for (i in 1:min(nrow(scaled), 5)) {
        uni_mat[i, ] <- as.numeric(scaled[i, ])
      }
      
      df <- melt(uni_mat, varnames = c("Class", "Antibiotic"), value.name = "value")
      df$Country <- country
      all_data <- rbind(all_data, df)
    }
    
    if (nrow(all_data) == 0) return(NULL)
    
    all_data$Class <- factor(all_data$Class, levels = paste0("Class", 1:5))
    all_data$Country <- factor(all_data$Country, levels = avail)
    
    n_countries <- length(avail)
    n_cols <- if (n_countries <= 3) n_countries else if (n_countries <= 4) 2 else 3
    n_rows <- ceiling(n_countries / n_cols)
    
    p <- ggplot(all_data, aes(x = Antibiotic, y = Class, fill = value)) +
      geom_tile(color = "white", linewidth = 0.3) +
      facet_wrap(~Country, scales = "free_y", ncol = n_cols) +
      scale_fill_gradientn(
        colors = RDBU_PALETTE,
        limits = c(-vmax, vmax),
        na.value = "grey95",
        name = "Scaled\nlog2(MIC)",
        guide = guide_colorbar(barwidth = 0.6, barheight = 4)
      ) +
      labs(title = paste("Panel", panel_label), x = "", y = "") +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "bold"),
        axis.text.y = element_text(size = 9, face = "bold"),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.spacing = unit(0.3, "cm")
      )
    
    fig_w <- max(8, n_cols * 3.5)
    fig_h <- max(6, n_rows * 3 + 1.5)
    
    save_figure(p, paste0("panel_", panel_label, "_heatmap"), 
                output_dir, "heatmaps", width = fig_w, height = fig_h)
    
  } else if (type == "trajectory") {
    for (country in avail) {
      res <- all_results[[country]]
      if (is.null(res$class_by_year)) next
      
      class_dist <- do.call(rbind, lapply(names(res$class_by_year), function(yr) {
        yr_data <- res$class_by_year[[yr]]
        counts <- table(factor(yr_data$classification, levels = 1:5))
        data.frame(
          year = as.integer(yr),
          class = paste0("Class", 1:5),
          pct = as.numeric(counts) / yr_data$n * 100
        )
      }))
      
      if (nrow(class_dist) > 0) {
        class_dist$Country <- country
        all_data <- rbind(all_data, class_dist)
      }
    }
    
    if (nrow(all_data) == 0) return(NULL)
    
    all_data$class <- factor(all_data$class, levels = paste0("Class", 1:5))
    all_data$Country <- factor(all_data$Country, levels = avail)
    
    n_countries <- length(avail)
    n_cols <- if (n_countries <= 3) n_countries else if (n_countries <= 4) 2 else 3
    n_rows <- ceiling(n_countries / n_cols)
    
    p <- ggplot(all_data, aes(x = year, y = pct, color = class, group = class)) +
      geom_line(size = 1) +
      geom_point(size = 2, shape = 21, fill = "white", stroke = 0.3) +
      facet_wrap(~Country, ncol = n_cols, scales = "free") +
      scale_color_manual(values = CLASS_COLORS, name = "Class", drop = FALSE) +
      scale_x_continuous(breaks = function(x) {
        if (length(unique(x)) > 6) return(unique(x)[round(seq(1, length(unique(x)), length.out = 4))])
        return(unique(x))
      }) +
      scale_y_continuous(labels = function(x) paste0(x, "%"),
                         limits = function(x) c(0, max(x, na.rm = TRUE) * 1.1)) +
      labs(title = paste("Panel", panel_label), x = "Year", y = "Proportion (%)") +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.spacing = unit(0.3, "cm")
      ) +
      guides(color = guide_legend(nrow = 1))
    
    fig_w <- max(8, n_cols * 3.5)
    fig_h <- max(6, n_rows * 3 + 1.5)
    
    save_figure(p, paste0("panel_", panel_label, "_trajectory"), 
                output_dir, "trajectories", width = fig_w, height = fig_h)
  }
  
  return(p)
}

# Create posterior probability plots
create_posterior_plots <- function(res, output_dir) {
  posterior <- res$posterior
  df <- res$df
  n_classes <- res$n_classes
  
  if (is.null(posterior)) return(NULL)
  
  post_long <- as.data.frame(posterior)
  colnames(post_long) <- paste0("Class", 1:ncol(posterior))
  post_long$assigned <- df$class_assignment
  
  post_long <- post_long %>%
    pivot_longer(cols = starts_with("Class"), 
                 names_to = "target", 
                 values_to = "prob")
  
  p <- ggplot(post_long, aes(x = target, y = prob, fill = target)) +
    geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, outlier.size = 0.3) +
    facet_wrap(~assigned, nrow = 1) +
    scale_fill_manual(values = CLASS_COLORS) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
    labs(title = paste("Posterior Probabilities -", res$country),
         x = "Target Class", y = "Probability") +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(face = "bold", size = 8),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3)
    )
  
  save_figure(p, "posterior_combined", output_dir, "posterior", width = 12, height = 5)
  
  summary_table <- post_long %>%
    group_by(assigned, target) %>%
    summarise(
      n = n(),
      mean = mean(prob, na.rm = TRUE),
      median = median(prob, na.rm = TRUE),
      sd = sd(prob, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ round(., 4)))
  
  save_table(summary_table, "posterior_summary", output_dir)
  
  for (assigned in 1:n_classes) {
    sub <- post_long %>% filter(assigned == assigned)
    if (nrow(sub) == 0) next
    
    p2 <- ggplot(sub, aes(x = target, y = prob, fill = target)) +
      geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, outlier.size = 0.3) +
      scale_fill_manual(values = CLASS_COLORS) +
      scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
      labs(title = paste("Class", assigned, "-", res$country),
           x = "Target Class", y = "Probability") +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(face = "bold"),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
      )
    
    save_figure(p2, paste0("posterior_class", assigned), output_dir,
                "posterior", width = 7, height = 5)
  }
  
  return(post_long)
}

# Analyze proportion by group (ICU/age)
analyze_proportion <- function(df, group_var, group_name, output_dir) {
  props <- df %>%
    group_by(class_uniform, !!sym(group_var)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(class_uniform) %>%
    mutate(total = sum(n), pct = n / total * 100) %>%
    filter(!is.na(!!sym(group_var)))
  
  if (nrow(props) == 0) return(NULL)
  
  colors <- if (group_name == "ICU") ICU_COLORS else AGE_COLORS
  
  p1 <- ggplot(props, aes(x = class_uniform, y = pct, fill = as.factor(!!sym(group_var)))) +
    geom_col(color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = ifelse(pct > 5, paste0(round(pct, 1), "%"), "")),
              position = position_stack(vjust = 0.5), size = 3, color = "white") +
    scale_fill_manual(values = colors, name = group_name) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
    labs(title = paste(group_name, "Proportion by Class"), x = "Class", y = "Proportion (%)") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
    )
  
  save_figure(p1, paste0(group_name, "_stacked"), output_dir, 
              tolower(group_name), width = 7, height = 5)
  
  yearly <- df %>%
    group_by(Year, class_uniform, !!sym(group_var)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Year, class_uniform) %>%
    mutate(total = sum(n), pct = n / total * 100) %>%
    filter(!is.na(!!sym(group_var)))
  
  if (nrow(yearly) > 0) {
    p2 <- ggplot(yearly, aes(x = Year, y = pct, color = as.factor(!!sym(group_var)), group = !!sym(group_var))) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      facet_wrap(~class_uniform) +
      scale_color_manual(values = colors, name = group_name) +
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      labs(title = paste(group_name, "Proportion Over Time"), x = "Year", y = "Proportion (%)") +
      theme_minimal(base_size = 10) +
      theme(
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3)
      )
    
    save_figure(p2, paste0(group_name, "_trajectory"), output_dir,
                tolower(group_name), width = 10, height = 6)
    
    save_table(yearly, paste0(group_name, "_yearly"), output_dir)
  }
  
  save_table(props, paste0(group_name, "_proportions"), output_dir)
  
  return(list(props = props, yearly = yearly))
}

# MAIN EXECUTION

# Load data
df_all <- load_data("data/amr_data.xlsx", sheet = 1)
cat(sprintf("\nData loaded: %d rows, %d columns\n", nrow(df_all), ncol(df_all)))

# Verify selected antibiotics
i_cols <- grep("_I$", names(df_all), value = TRUE)
all_mic <- str_remove(i_cols, "_I$")

found_ab <- c()
for (ab in SELECTED_ANTIBIOTICS) {
  ab_clean <- make.names(ab)
  if (ab_clean %in% all_mic) {
    found_ab <- c(found_ab, ab)
    cat(sprintf("  ✅ %s -> FOUND\n", ab))
  } else {
    matches <- all_mic[grepl(gsub(" ", "|", ab), all_mic, ignore.case = TRUE)]
    if (length(matches) > 0) {
      found_ab <- c(found_ab, ab)
      cat(sprintf("  ✅ %s -> FOUND (partial)\n", ab))
    } else {
      cat(sprintf("  ❌ %s -> NOT FOUND\n", ab))
    }
  }
}

SELECTED_AB_CLEAN <- make.names(found_ab)
if (length(SELECTED_AB_CLEAN) < 3) {
  stop("Insufficient antibiotics found!")
}

# Get countries with sufficient isolates
countries <- sort(unique(df_all$Country))
country_n <- df_all %>%
  filter(Species == "Klebsiella pneumoniae", Source == "Blood") %>%
  group_by(Country) %>%
  summarise(n = n(), .groups = "drop")

countries_to_analyze <- country_n %>%
  filter(n >= MIN_ISOLATES_PER_COUNTRY) %>%
  pull(Country)

cat(sprintf("\nAnalyzing %d countries\n", length(countries_to_analyze)))

# Create output directory
OUTPUT_DIR <- "results"
create_dirs(OUTPUT_DIR)

# Analyze each country
all_results <- list()
for (country in countries_to_analyze) {
  tryCatch({
    res <- analyze_country(df_all, country, OUTPUT_DIR, SELECTED_AB_CLEAN)
    if (!is.null(res)) all_results[[country]] <- res
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}

cat(sprintf("\n✅ Analyzed: %d countries\n", length(all_results)))

# Create visualizations
cat("\nCREATING VISUALIZATIONS\n")

# Heatmaps
for (country in names(all_results)) {
  res <- all_results[[country]]
  if (!is.null(res$log2prof)) {
    create_heatmap(res$log2prof, country, file.path(OUTPUT_DIR, country))
  }
}

# Trajectories
for (country in names(all_results)) {
  res <- all_results[[country]]
  if (!is.null(res$class_by_year)) {
    create_trajectory(res$class_by_year, country, file.path(OUTPUT_DIR, country))
  }
}

# Panel figures
panel_labels <- c("A", "B", "C")
for (i in seq_along(COUNTRY_PANELS)) {
  panel_name <- names(COUNTRY_PANELS)[i]
  countries <- COUNTRY_PANELS[[panel_name]]
  label <- panel_labels[i]
  
  create_panel_figure(countries, label, all_results, OUTPUT_DIR, type = "heatmap")
  create_panel_figure(countries, label, all_results, OUTPUT_DIR, type = "trajectory")
}

# ICU & Age analyses
for (country in names(all_results)) {
  res <- all_results[[country]]
  df <- res$df
  country_dir <- file.path(OUTPUT_DIR, country)
  
  if (any(!is.na(df$icu_flag))) {
    analyze_proportion(df, "icu_flag", "ICU", country_dir)
  }
  
  if ("age_group" %in% names(df) && any(!is.na(df$age_group))) {
    analyze_proportion(df, "age_group", "Age_Group", country_dir)
  }
}

# Posterior probability plots
for (country in names(all_results)) {
  res <- all_results[[country]]
  create_posterior_plots(res, file.path(OUTPUT_DIR, country))
}

# SUMMARY

cat("\nANALYSIS COMPLETE!\n")
cat(sprintf("Countries analyzed: %d\n", length(all_results)))
cat(sprintf("Uniform classes: %d (1-5)\n", MAX_CLASSES))

cat("\nAntibiotics Used:\n")
for (ab in SELECTED_AB_CLEAN) {
  orig_ab <- found_ab[which(make.names(found_ab) == ab)]
  if (length(orig_ab) > 0) {
    cat(sprintf("  - %s\n", orig_ab))
  }
}

cat("\nOutput Directory: results/\n")
cat("  - figures/         : All visualizations\n")
cat("  - tables/          : Summary tables\n")
cat("  - models/          : Model outputs\n")
cat("  - reports/         : Analysis reports\n")

cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")