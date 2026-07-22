# ============================================================================
# METHODS
#   Step 1: Fit reference model on pooled data (all years) 
#   Step 2: Classify all isolates using reference model parameters
#   Step 3: BCH correction for ICU association (if applicable)
# ============================================================================

DATA_PATH <- "data/amr_data.xlsx"
SHEET     <- "Sheet1"

# Output directory
OUTPUT_DIR <- "output"

# Parameter analisis
K_RANGE_MAX          <- 5
MAX_MISSING_PROP     <- 0.30
MIN_N_COMPLETE       <- 10
TOP_N_ANTIBIOTIK     <- 10
RESISTANCE_THRESHOLD <- 0.5
MIC_GE_MULTIPLIER    <- 2
MIN_YEAR_COVERAGE    <- 0.70
MCLUST_MODEL_NAMES   <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
RUN_BOOTSTRAP        <- TRUE
R_BOOT               <- 200
RUN_BCH_RQ2          <- TRUE
USE_POSTHOC_SIR_LABELS <- TRUE
RANDOM_SEED          <- 123

# Parameter Three-Step
USE_IMPUTATION      <- TRUE
IMPUTATION_METHOD   <- "mean"
MIN_CLASS_SIZE      <- 10
CONFIDENCE_LEVEL    <- 0.95

pkgs_needed <- c("mclust", "boot", "readxl", "dplyr", "tidyr", "stringr",
                 "ggplot2", "scales", "survey", "gridExtra", "knitr", 
                 "forcats", "viridis", "ggpubr", "patchwork", "psych",
                 "corrplot", "reshape2", "grid", "cowplot", "RColorBrewer",
                 "ggdendro", "dendextend", "openxlsx")  
pkgs_missing <- pkgs_needed[!sapply(pkgs_needed, requireNamespace, quietly = TRUE)]
if (length(pkgs_missing) > 0) {
  message("Installing missing packages: ", paste(pkgs_missing, collapse = ", "))
  install.packages(pkgs_missing)
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
  library(survey)
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
})

create_output_dirs <- function(base_dir) {
  dirs <- c(
    "tables",
    "figures"
  )
  for (d in dirs) {
    dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
  }
  message(sprintf("Output directories created in: %s", base_dir))
}

create_output_dirs(OUTPUT_DIR)

# LOAD DATA 
load_amr_data <- function(path, sheet = 1) {
  as.data.frame(read_excel(path, sheet = sheet))
}

# FLAG ICU vs NON-ICU
flag_icu <- function(speciality_vec) {
  grepl("ICU", speciality_vec, ignore.case = TRUE)
}

# PARSE RAW MIC STRINGS 
parse_mic_numeric <- function(x, ge_multiplier = 2) {
  x <- trimws(as.character(x))
  is_ge <- grepl("^>=?", x)
  num <- suppressWarnings(as.numeric(gsub("[<>=]", "", x)))
  num[is_ge] <- num[is_ge] * ge_multiplier
  num
}

# BUILD LOG2(MIC) MATRIX
build_mic_matrix_log2 <- function(df, ge_multiplier = 2) {
  i_cols <- grep("_I$", names(df), value = TRUE)
  mic_cols <- str_remove(i_cols, "_I$")
  mic_cols <- mic_cols[mic_cols %in% names(df)]
  if (length(mic_cols) == 0) stop("No raw MIC columns found matching the '*_I' antibiotic names.")
  mat <- df[, mic_cols, drop = FALSE]
  mat[] <- lapply(mat, parse_mic_numeric, ge_multiplier = ge_multiplier)
  mat[] <- lapply(mat, function(x) { x[!is.na(x) & x <= 0] <- NA; log2(x) })
  names(mat) <- make.names(names(mat))
  as.data.frame(mat)
}

# SELECT CORE PANEL ACROSS ALL YEARS 
select_core_panel <- function(df, year_col = "Year", min_coverage = 0.70) {
  years <- unique(df[[year_col]])
  
  year_matrices <- lapply(years, function(y) {
    df_yr <- df %>% filter(!!sym(year_col) == y)
    mat <- build_mic_matrix_log2(df_yr, ge_multiplier = MIC_GE_MULTIPLIER)
    return(mat)
  })
  names(year_matrices) <- as.character(years)
  
  all_antibiotics <- unique(unlist(lapply(year_matrices, names)))
  
  coverage <- sapply(all_antibiotics, function(ab) {
    mean(sapply(year_matrices, function(mat) ab %in% names(mat)))
  })
  
  core_antibiotics <- names(coverage)[coverage >= min_coverage]
  
  valid_antibiotics <- sapply(core_antibiotics, function(ab) {
    all(sapply(year_matrices, function(mat) {
      if (!ab %in% names(mat)) return(FALSE)
      prop_complete <- sum(!is.na(mat[[ab]])) / nrow(mat)
      prop_complete >= 0.80
    }))
  })
  
  core_panel <- core_antibiotics[valid_antibiotics]
  
  list(
    core_panel = core_panel,
    coverage = coverage[core_panel],
    all_antibiotics = all_antibiotics,
    year_matrices = year_matrices,
    years = years
  )
}

# IMPUTE MISSING VALUES 
impute_missing <- function(mat, method = "mean") {
  for (j in seq_len(ncol(mat))) {
    if (method == "mean") {
      col_means <- mean(mat[, j], na.rm = TRUE)
      if (is.na(col_means)) next
      mat[is.na(mat[, j]), j] <- col_means
    } else if (method == "median") {
      col_med <- median(mat[, j], na.rm = TRUE)
      if (is.na(col_med)) next
      mat[is.na(mat[, j]), j] <- col_med
    }
  }
  mat
}

# FIT LPA WITH AUTOMATIC BIC SELECTION 
fit_lpa_auto <- function(mat, k_range = 2:5, model_names = MCLUST_MODEL_NAMES, 
                         seed = 123, verbose = TRUE) {
  set.seed(seed)
  
  results <- list()
  best_bic <- Inf
  best_fit <- NULL
  best_k <- NA
  
  for (k in k_range) {
    if (verbose) cat(sprintf("  Fitting G=%d...", k))
    fit <- tryCatch(
      mclust::Mclust(data = as.matrix(mat), G = k, modelNames = model_names, 
                     verbose = FALSE),
      error = function(e) {
        if (verbose) cat(sprintf(" failed: %s\n", conditionMessage(e)))
        NULL
      }
    )
    
    if (!is.null(fit) && length(fit$modelName) > 0) {
      if (verbose) cat(sprintf(" done (BIC=%.1f)\n", -fit$bic))
      results[[as.character(k)]] <- fit
      bic_val <- -fit$bic
      if (bic_val < best_bic) {
        best_bic <- bic_val
        best_fit <- fit
        best_k <- k
      }
    } else {
      if (verbose) cat(" failed\n")
    }
  }
  
  if (is.null(best_fit)) {
    warning("No model converged!")
    return(NULL)
  }
  
  list(
    all_fits = results,
    best_fit = best_fit,
    best_k = best_k,
    best_bic = best_bic
  )
}

# COMPUTE ENTROPY 
compute_entropy_matrix <- function(post) {
  n <- nrow(post); k <- ncol(post)
  if (k <= 1) return(NA_real_)
  ent_i <- -rowSums(post * log(post + 1e-12))
  1 - sum(ent_i) / (n * log(k))
}

# MODEL COMPARISON TABLE
model_selection_table_lpa <- function(fit_list) {
  rows <- lapply(names(fit_list), function(g) {
    fit <- fit_list[[g]]
    if (is.null(fit)) {
      return(data.frame(G = as.integer(g), model = NA, BIC = NA, entropy = NA,
                        min_class_pct = NA, n_used = NA))
    }
    post <- fit$z
    tab  <- table(factor(fit$classification, levels = 1:fit$G))
    data.frame(
      G             = as.integer(g),
      model         = fit$modelName,
      BIC           = -fit$bic,
      entropy       = compute_entropy_matrix(post),
      min_class_pct = min(tab) / sum(tab) * 100,
      n_used        = fit$n
    )
  })
  tbl <- do.call(rbind, rows)
  tbl <- tbl[order(tbl$G), ]
  tbl$BIC_rank <- rank(tbl$BIC, na.last = "keep")
  tbl$selected <- ifelse(tbl$BIC == min(tbl$BIC, na.rm = TRUE), "***", "")
  tbl
}

# POSTERIOR DIAGNOSTICS 
posterior_diagnostics_lpa <- function(fit) {
  post <- fit$z
  assigned <- fit$classification
  ks <- sort(unique(assigned))
  avg_pp <- sapply(ks, function(k) mean(post[assigned == k, k]))
  data.frame(
    class = paste0("Class", ks),
    n = as.integer(table(factor(assigned, levels = ks))),
    pct = as.numeric(table(factor(assigned, levels = ks))) / length(assigned) * 100,
    avg_posterior_prob = avg_pp
  )
}

# EXTRACT MIC PROFILE
extract_mic_profile_log2 <- function(fit, var_names) {
  means <- t(fit$parameters$mean)
  means <- as.data.frame(means)
  names(means) <- var_names
  rownames(means) <- paste0("Class", seq_len(nrow(means)))
  means
}

geometric_mean_mic_profile <- function(log2_profile) {
  out <- as.data.frame(round(2 ^ as.matrix(log2_profile), 3))
  rownames(out) <- rownames(log2_profile)
  out
}

# PREDICT CLASSES USING REFERENCE MODE
predict_lpa_classes <- function(fit, newdata) {
  pred <- predict(fit, newdata = as.matrix(newdata))
  list(
    classification = pred$classification,
    uncertainty = pred$uncertainty,
    posterior = pred$z
  )
}

# POST-HOC S/I/R PROFILE 
compute_posthoc_resistance_profile <- function(df_rows, modal_class, var_names) {
  i_cols <- paste0(var_names, "_I")
  i_cols <- i_cols[i_cols %in% names(df_rows)]
  if (length(i_cols) == 0) return(NULL)
  
  mat_i <- df_rows[, i_cols, drop = FALSE]
  names(mat_i) <- str_remove(names(mat_i), "_I$")
  ks <- sort(unique(modal_class))
  
  prof <- sapply(names(mat_i), function(v) {
    x <- as.character(mat_i[[v]])
    sapply(ks, function(k) {
      xk <- x[modal_class == k]
      xk <- xk[!is.na(xk)]
      if (length(xk) == 0) return(NA_real_)
      mean(xk == "Resistant")
    })
  })
  prof <- as.data.frame(prof)
  rownames(prof) <- paste0("Class", ks)
  prof
}

# BOOTSTRAP FOR CLASS PROPORTIONS 
bootstrap_lpa_se <- function(mat, fit, R = 200, seed = 123) {
  set.seed(seed)
  k <- fit$G
  model_names <- fit$modelName
  
  fit_once <- function(d) {
    fit_tmp <- tryCatch(
      mclust::Mclust(data = as.matrix(d), G = k, modelNames = model_names, 
                     verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit_tmp) || is.null(fit_tmp$parameters)) return(rep(NA_real_, k))
    sort(fit_tmp$parameters$pro)
  }
  
  boot_stat <- function(data, indices) fit_once(data[indices, , drop = FALSE])
  boot_obj <- boot::boot(data = mat, statistic = boot_stat, R = R)
  
  ests <- boot_obj$t
  se <- apply(ests, 2, sd, na.rm = TRUE)
  ci <- t(apply(ests, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  
  data.frame(
    class = paste0("Class", 1:k),
    proportion = round(boot_obj$t0, 3),
    se = round(se, 3),
    ci_lower = round(ci[, 1], 3),
    ci_upper = round(ci[, 2], 3),
    n_valid_reps = colSums(!is.na(ests))
  )
}

# SAVE FUNCTIONS 
save_table <- function(data, filename, output_dir, format = "csv") {
  path <- file.path(output_dir, "tables", filename)
  if (format == "csv") {
    write.csv(data, paste0(path, ".csv"), row.names = FALSE)
  } else if (format == "txt") {
    write.table(data, paste0(path, ".txt"), row.names = FALSE, sep = "\t")
  } else if (format == "rds") {
    saveRDS(data, paste0(path, ".rds"))
  }
  message(sprintf("  Saved: %s", paste0(path, ".", format)))
}

save_figure <- function(plot, filename, output_dir, width = 12, height = 8, dpi = 300) {
  path <- file.path(output_dir, "figures", filename)
  ggsave(paste0(path, ".png"), plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(path, ".pdf"), plot, width = width, height = height, dpi = dpi)
  message(sprintf("  Saved: %s", path))
}

save_model <- function(model, filename, output_dir) {
  path <- file.path(output_dir, "models", filename)
  saveRDS(model, paste0(path, ".rds"))
  message(sprintf("  Saved model: %s", path))
}

# ============================================================================
# LOAD DATA
# ============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("AMR RESISTOME LATENT PROFILE ANALYSIS - THREE-STEP APPROACH\n")
cat("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

df_all <- load_amr_data(DATA_PATH, sheet = SHEET)
cat(sprintf("\nData loaded: %d rows, %d columns.\n", nrow(df_all), ncol(df_all)))

# SELECT CORE PANEL ACROSS ALL YEARS
cat("\n", paste(rep("-", 80), collapse = ""), "\n")
cat("STEP 1: Selecting Core Antibiotic Panel\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

df_kp_blood <- df_all %>% filter(Species == "Klebsiella pneumoniae", Source == "Blood")
cat(sprintf("\nTotal KP + Blood isolates: %d\n", nrow(df_kp_blood)))

years <- sort(unique(df_kp_blood$Year))
group_n <- sapply(years, function(y) sum(df_kp_blood$Year == y))
cat("\nIsolates per year:\n")
print(group_n)

# Select core panel
core_info <- select_core_panel(df_kp_blood, year_col = "Year", min_coverage = MIN_YEAR_COVERAGE)
core_panel <- core_info$core_panel

# Sort antibiotics alphabetically
core_panel <- sort(core_panel)

cat(sprintf("\nCore panel selected: %d antibiotics (coverage >= %.0f%%)\n", 
            length(core_panel), MIN_YEAR_COVERAGE * 100))
cat("Antibiotics (alphabetical):\n")
print(core_panel)

# Save core panel info
save_table(data.frame(
  antibiotic = core_panel,
  coverage = round(core_info$coverage[core_panel], 3)
), "core_panel_info", OUTPUT_DIR)


# SECTION 2: BUILD POOLED DATA AND SELECT VARIABLES
cat("\n", paste(rep("-", 80), collapse = ""), "\n")
cat("STEP 2: Building Pooled Dataset\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Build MIC matrix for all isolates
mat_all <- build_mic_matrix_log2(df_kp_blood, ge_multiplier = MIC_GE_MULTIPLIER)

# Keep only core panel
mat_core <- mat_all[, intersect(names(mat_all), core_panel), drop = FALSE]

cat(sprintf("Pooled dataset: %d isolates, %d antibiotics\n", 
            nrow(mat_core), ncol(mat_core)))

# Check missingness
missing_props <- sapply(mat_core, function(x) mean(is.na(x)))
cat("\nMissingness per antibiotic:\n")
print(round(missing_props, 3))

# Impute missing values
if (USE_IMPUTATION) {
  cat(sprintf("\nImputing missing values using %s imputation...\n", IMPUTATION_METHOD))
  mat_imp <- impute_missing(mat_core, method = IMPUTATION_METHOD)
  cat(sprintf("Total missing values imputed: %d\n", sum(is.na(mat_core))))
} else {
  complete_idx <- complete.cases(mat_core)
  mat_imp <- mat_core[complete_idx, ]
  df_kp_blood <- df_kp_blood[complete_idx, ]
  cat(sprintf("Complete cases: %d out of %d (%.1f%%)\n", 
              nrow(mat_imp), nrow(mat_core), 
              100 * nrow(mat_imp) / nrow(mat_core)))
}

# Save imputed data
save_table(mat_imp, "pooled_imputed_data", OUTPUT_DIR, format = "csv")


# STEP 1 - FIT REFERENCE MODEL
cat("\n", paste(rep("-", 80), collapse = ""), "\n")
cat("STEP 3: Fitting Reference Model on Pooled Data\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Determine number of classes to test
n_total <- nrow(mat_imp)
k_max <- min(K_RANGE_MAX, floor(n_total / 30))
k_range <- 2:max(2, min(k_max, 5))

cat(sprintf("\nTesting models with G = %s\n", paste(k_range, collapse = ", ")))

# Fit reference model
ref_model <- fit_lpa_auto(mat_imp, k_range = k_range, model_names = MCLUST_MODEL_NAMES, 
                          seed = RANDOM_SEED, verbose = TRUE)

if (is.null(ref_model)) {
  stop("No reference model converged. Cannot proceed.")
}

cat(sprintf("\nReference model selected: G=%d (%s), BIC=%.1f\n", 
            ref_model$best_k, ref_model$best_fit$modelName, ref_model$best_bic))

# Model comparison table
model_comp <- model_selection_table_lpa(ref_model$all_fits)
cat("\nModel comparison:\n")
print(model_comp)

save_table(model_comp, "reference_model_comparison", OUTPUT_DIR)

# Extract profiles
vars <- names(mat_imp)
log2prof_ref <- extract_mic_profile_log2(ref_model$best_fit, vars)

# Sort antibiotics alphabetically
log2prof_ref <- log2prof_ref[, sort(colnames(log2prof_ref))]

geoprof_ref <- geometric_mean_mic_profile(log2prof_ref)

# Class diagnostics
class_diag_ref <- posterior_diagnostics_lpa(ref_model$best_fit)
cat("\nClass distribution in reference model:\n")
print(class_diag_ref, row.names = FALSE)

save_table(class_diag_ref, "reference_class_diagnostics", OUTPUT_DIR)

# Save reference model
save_model(ref_model, "reference_model", OUTPUT_DIR)

# STEP 2 - CLASSIFY ISOLATES BY YEAR USING REFERENCE MODEL
cat("\n", paste(rep("-", 80), collapse = ""), "\n")
cat("STEP 4: Classifying Isolates by Year\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Storage for classification results
classification_results <- list()
posthoc_results <- list()
phenotype_results <- list()
bootstrap_results <- list()

for (yr in years) {
  cat(sprintf("\n=== Year %s (n=%d) ===\n", yr, group_n[which(years == yr)]))
  
  # Get year-specific data
  df_yr <- df_kp_blood %>% filter(Year == yr)
  mat_yr <- build_mic_matrix_log2(df_yr, ge_multiplier = MIC_GE_MULTIPLIER)
  mat_yr_core <- mat_yr[, intersect(names(mat_yr), core_panel), drop = FALSE]
  
  # Handle missing values
  if (USE_IMPUTATION) {
    mat_yr_imp <- impute_missing(mat_yr_core, method = IMPUTATION_METHOD)
  } else {
    complete_idx <- complete.cases(mat_yr_core)
    mat_yr_imp <- mat_yr_core[complete_idx, ]
    df_yr <- df_yr[complete_idx, ]
  }
  
  # Check if we have enough data
  if (nrow(mat_yr_imp) < 10) {
    cat(sprintf("  WARNING: Only %d isolates after processing. Skipping.\n", nrow(mat_yr_imp)))
    next
  }
  
  # Predict classes using reference model
  pred <- predict_lpa_classes(ref_model$best_fit, mat_yr_imp)
  
  # Store results
  classification_results[[as.character(yr)]] <- list(
    year = yr,
    n = nrow(mat_yr_imp),
    classification = pred$classification,
    uncertainty = pred$uncertainty,
    posterior = pred$posterior,
    class_proportions = table(pred$classification) / nrow(mat_yr_imp),
    mat_imp = mat_yr_imp
  )
  
  # Class distribution
  class_levels <- 1:ref_model$best_k
  class_counts <- table(factor(pred$classification, levels = class_levels))
  
  uncertainty_by_class <- rep(NA, length(class_levels))
  names(uncertainty_by_class) <- paste0("Class", class_levels)
  
  for (i in seq_along(class_levels)) {
    class_name <- class_levels[i]
    idx <- which(pred$classification == class_name)
    if (length(idx) > 0) {
      uncertainty_by_class[i] <- mean(pred$uncertainty[idx], na.rm = TRUE)
    }
  }
  
  class_dist <- data.frame(
    class = paste0("Class", class_levels),
    n = as.integer(class_counts),
    pct = as.numeric(class_counts) / nrow(mat_yr_imp) * 100,
    avg_uncertainty = round(uncertainty_by_class, 4)
  )
  
  cat("\n  Class distribution:\n")
  print(class_dist, row.names = FALSE)
  
  save_table(class_dist, sprintf("year_%s_class_distribution", yr), OUTPUT_DIR)
  
  # Post-hoc S/I/R profiles
  if (USE_POSTHOC_SIR_LABELS) {
    sir_prof <- compute_posthoc_resistance_profile(df_yr, pred$classification, vars)
    if (!is.null(sir_prof)) {
      posthoc_results[[as.character(yr)]] <- list(
        year = yr,
        sir_profile = sir_prof,
        class_avg_resistance = rowMeans(sir_prof, na.rm = TRUE)
      )
    }
  }
  
  # Bootstrap
  if (RUN_BOOTSTRAP && nrow(mat_yr_imp) >= 30) {
    cat(sprintf("  Running bootstrap (R=%d)...\n", R_BOOT))
    boot_result <- bootstrap_lpa_se(mat_yr_imp, ref_model$best_fit, R = R_BOOT)
    boot_result$year <- yr
    bootstrap_results[[as.character(yr)]] <- boot_result
    cat("\n  Bootstrap 95% CIs for class proportions:\n")
    print(boot_result, row.names = FALSE)
    save_table(boot_result, sprintf("year_%s_bootstrap", yr), OUTPUT_DIR)
  }
}

# COMPILE YEARLY RESULTS
cat("\n", paste(rep("-", 80), collapse = ""), "\n")
cat("STEP 5: Compiling Yearly Results\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Compile class distributions across years
class_dist_all <- do.call(rbind, lapply(names(classification_results), function(yr) {
  res <- classification_results[[yr]]
  class_levels <- 1:ref_model$best_k
  
  class_counts <- table(factor(res$classification, levels = class_levels))
  
  uncertainty_by_class <- rep(NA, length(class_levels))
  names(uncertainty_by_class) <- paste0("Class", class_levels)
  
  for (i in seq_along(class_levels)) {
    class_name <- class_levels[i]
    idx <- which(res$classification == class_name)
    if (length(idx) > 0) {
      uncertainty_by_class[i] <- mean(res$uncertainty[idx], na.rm = TRUE)
    }
  }
  
  data.frame(
    year = yr,
    class = paste0("Class", class_levels),
    n = as.integer(class_counts),
    pct = as.numeric(class_counts) / res$n * 100,
    avg_uncertainty = uncertainty_by_class
  )
}))

class_dist_all$avg_uncertainty[is.na(class_dist_all$avg_uncertainty)] <- 0
save_table(class_dist_all, "class_distribution_all_years", OUTPUT_DIR)

# Summary table
summary_table <- do.call(rbind, lapply(names(classification_results), function(yr) {
  res <- classification_results[[yr]]
  class_factors <- factor(res$classification, levels = 1:ref_model$best_k)
  
  data.frame(
    year = yr,
    n = res$n,
    n_classes_observed = length(unique(res$classification)),
    entropy = compute_entropy_matrix(res$posterior),
    min_class_size = min(table(class_factors)),
    max_class_size = max(table(class_factors))
  )
}))
save_table(summary_table, "yearly_summary", OUTPUT_DIR)

# VISUALIZATION - CLASS DISTRIBUTION
cat("\n", paste(rep("-", 80), collapse = ""), "\n")
cat("STEP 6: Generating Visualizations\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Get class order (Class1, Class2, ...)
class_order <- paste0("Class", 1:ref_model$best_k)

# Class distribution over years (stacked bar) 
if (exists("class_dist_all") && nrow(class_dist_all) > 0) {
  plot_data <- class_dist_all
  plot_data$year <- factor(plot_data$year, levels = sort(unique(plot_data$year)))
  plot_data$class <- factor(plot_data$class, levels = class_order)
  
  p1 <- ggplot(plot_data, aes(x = year, y = pct, fill = class)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = ifelse(pct >= 5, sprintf("%.0f%%", pct), "")),
              position = position_stack(vjust = 0.5),
              size = 3.5, fontface = "bold", color = "white", na.rm = TRUE) +
    scale_fill_manual(values = pub_colors[class_order], name = "Class") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                       labels = percent_format(scale = 1)) +
    labs(
      title = "Latent Profile Class Distribution by Year",
      subtitle = "K. pneumoniae blood isolates",
      x = "Year",
      y = "Class Proportion (%)"
    ) +
    theme_publication(base_size = 12)
  
  save_figure(p1, "class_distribution_stacked", OUTPUT_DIR, width = 10, height = 7)
  
  # Class trajectories
  p2 <- ggplot(plot_data, aes(x = year, y = pct, color = class, group = class)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = pub_colors[class_order], name = "Class") +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    labs(
      title = "Class Proportion Trends Over Time",
      x = "Year",
      y = "Class Proportion (%)"
    ) +
    theme_publication(base_size = 12)
  
  save_figure(p2, "class_trajectories", OUTPUT_DIR, width = 10, height = 6)
}

# Heatmap of MIC profiles (log2 MIC - no scaling) 
if (exists("log2prof_ref") && !is.null(log2prof_ref) && nrow(log2prof_ref) > 0) {
  
  # Ensure antibiotics are sorted alphabetically
  antibiotic_order_heat <- sort(colnames(log2prof_ref))
  
  # Ensure classes are in order (Class1, Class2, ...)
  class_order_heat <- class_order
  
  # Reorder data
  log2prof_ref_ordered <- log2prof_ref[class_order_heat, antibiotic_order_heat, drop = FALSE]
  
  log2prof_ref_matrix <- as.matrix(log2prof_ref_ordered)
  rownames(log2prof_ref_matrix) <- rownames(log2prof_ref_ordered)
  colnames(log2prof_ref_matrix) <- gsub("\\.", " ", colnames(log2prof_ref_ordered))
  
  if (requireNamespace("reshape2", quietly = TRUE)) {
    p3 <- ggplot(reshape2::melt(log2prof_ref_matrix), aes(x = Var2, y = Var1, fill = value)) +
      geom_tile() +
      scale_fill_viridis_c(name = "log2(MIC)") +
      labs(
        title = "MIC Profile Heatmap by Class",
        subtitle = "Antibiotics ordered alphabetically",
        x = "Antibiotic",
        y = "Class"
      ) +
      theme_publication(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    save_figure(p3, "mic_profile_heatmap_log2", OUTPUT_DIR, width = 14, height = 6)
  }
}

# Bootstrap confidence intervals
if (length(bootstrap_results) > 0) {
  boot_plot_data <- do.call(rbind, lapply(names(bootstrap_results), function(yr) {
    res <- bootstrap_results[[yr]]
    res$year <- yr
    res$class <- factor(res$class, levels = class_order)
    return(res)
  }))
  
  if (!is.null(boot_plot_data) && nrow(boot_plot_data) > 0) {
    p6 <- ggplot(boot_plot_data, aes(x = year, y = proportion, color = class)) +
      geom_point(size = 3, position = position_dodge(width = 0.3)) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                    width = 0.2, position = position_dodge(width = 0.3)) +
      scale_color_manual(values = pub_colors[class_order], name = "Class") +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      labs(
        title = "Class Proportions with 95% Bootstrap CIs",
        x = "Year",
        y = "Class Proportion"
      ) +
      theme_publication(base_size = 11)
    
    save_figure(p6, "bootstrap_ci_proportions", OUTPUT_DIR, width = 10, height = 6)
  }
}

#ICU PROPORTION VISUALIZATION
cat("\n", paste(rep("-", 80), collapse = ""), "\n")
cat("STEP 7: ICU Proportion Visualization\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Add ICU flag to data
df_kp_blood$icu_flag <- flag_icu(df_kp_blood$Speciality)

# Get classifications for all isolates
if (!exists("pred_all")) {
  pred_all <- predict_lpa_classes(ref_model$best_fit, mat_imp)
}

# Create dataset with class assignments and ICU status
icu_class_data <- data.frame(
  class = paste0("Class", pred_all$classification),
  icu_flag = df_kp_blood$icu_flag,
  stringsAsFactors = FALSE
)

# Calculate overall ICU proportions
icu_overall <- icu_class_data %>%
  group_by(class, icu_flag) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(class) %>%
  mutate(
    total = sum(count),
    pct = round(count / total * 100, 2)
  ) %>%
  ungroup() %>%
  mutate(
    icu_status = ifelse(icu_flag, "ICU", "Non-ICU")
  ) %>%
  select(class, icu_status, count, total, pct) %>%
  arrange(class, desc(icu_status))

# Ensure class order
icu_overall$class <- factor(icu_overall$class, levels = class_order)

# Stacked Bar Chart
p_icu_stacked <- ggplot(icu_overall, aes(x = class, y = pct, fill = icu_status)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = ifelse(pct >= 5, paste0(round(pct, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold", color = "white") +
  scale_fill_manual(values = pub_icu_colors, name = "ICU Status") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     labels = percent_format(scale = 1),
                     limits = c(0, 100)) +
  labs(
    title = "ICU vs Non-ICU Proportion by Latent Class",
    subtitle = "K. pneumoniae blood isolates, all years combined",
    x = "Latent Class",
    y = "Proportion of Isolates (%)"
  ) +
  theme_publication(base_size = 12)

save_figure(p_icu_stacked, "ICU_Proportion_Stacked", OUTPUT_DIR, width = 8, height = 6)

#EXPORT DATA TO EXCEL
cat("\n", paste(rep("-", 80), collapse = ""), "\n")
cat("STEP 8: Exporting Data to Excel\n")
cat(paste(rep("-", 80), collapse = ""), "\n")