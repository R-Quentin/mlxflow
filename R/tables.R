    ### ---- MODULE 1b: Nice outputs ----
#' Merges pop_params + shrinkage (from indiv_params) into a single display-ready
#'
#' tibble for use with kable / gt.
#'
#' Shrinkage values are joined on the base parameter name (stripping the _pop
#' suffix). Parameters with no random effect (e.g. D3, ka1, V3) will have NA.
#'
#' @param summary_obj output of read_summary(..., full = TRUE)
#' @param shrinkage_method "CD" (Conditional Distribution, default) or "CM" (Conditional Mode)
#' @param add_criteria if TRUE, append OFV / AIC / BICc / BIC rows at the bottom under an "Information Criteria" group
#'
#' @export
build_param_table <- function(summary_obj,
                              shrinkage_method = "CD",
                              add_criteria     = TRUE) {
  
  pop   <- summary_obj$pop_params
  indiv <- summary_obj$indiv_params
  
  if (is.null(pop)) stop("pop_params is NULL - re-run read_summary() with full = TRUE")
  
  method_label <- if (shrinkage_method == "CM") "Conditional Mode" else "Conditional Distribution"
  
  shrink <- if (!is.null(indiv) && nrow(indiv) > 0) {
    indiv |>
      filter(method == method_label) |>
      select(parameter, shrinkage)   # parameter names here have no _pop suffix
  } else {
    NULL
  }
  
  tbl <- pop |>
    mutate(param_base = strip_pop(parameter))
  
  if (!is.null(shrink)) {
    tbl <- tbl |>
      left_join(shrink |> rename(param_base = parameter),
                by = "param_base")
  } else {
    tbl$shrinkage <- NA_real_
  }
  
  
  if (add_criteria) {
    crit <- summary_obj$criteria |>
      rename(parameter = criterion, estimate = value) |>
      mutate(
        type      = "Information Criteria",
        se_sa     = NA_real_, rse_pct   = NA_real_,
        ci_low    = NA_real_, ci_high   = NA_real_,
        shrinkage = NA_real_
      )
    tbl <- bind_rows(tbl, crit)
  }
  
  
  tbl <- tbl |>
    select(any_of(c("parameter", "estimate", "se_sa", "rse_pct", "ci_low", "ci_high", "shrinkage", "type"))) |>
    select(where(function(x) any(!is.na(x)))) # Remove absent (empty) columns
  
  
  tbl
}
#' 2) clean_param_labels
#'
#' Derives clean covariate labels from raw "beta_" parameter names
#'
#' Uses the "*_pop" rows already present in `tbl` as the ground-truth list of real parameter names.
#'
#' @param tbl output of build_param_table()
#' @param label_col name of the new column to create (default "label")
#' @param prefix string prepended to derived covariate labels only (default "~"); use "" to disable
#' @param overwrite if TRUE, replaces `parameter` in place instead of adding a new column (if FALSE, raw Monolix names are preserved, label column is added)
#'
#' @return tbl with an additional (or replacing) clean label column. Rows that are not "beta_" covariate effects are left untouched. /!\\ The logic of logt/t as a prefix is based on Monolix's automatic naming of transformed covariates. If a covariate whose real name begins with a "t" followed by a capital letter is ever introduced, it will have to be explicitly excluded.
#'
#' @export
clean_param_labels <- function(tbl,
                               label_col = "label",
                               prefix    = "~",
                               overwrite = TRUE) {
  
  stopifnot("parameter" %in% names(tbl))
  
  # Ground-truth PK parameter names, from the "*_pop" rows already in tbl,
  # longest first so e.g. "ka_IM" is tried before "ka"
  param_names <- tbl$parameter |>
    keep(~ str_detect(.x, "_pop$")) |>
    str_remove("_pop$") |>
    unique()
  param_names <- param_names[order(-nchar(param_names))]
  
  strip_beta_prefix <- function(x) {
    if (!str_detect(x, "^beta_") || length(param_names) == 0) return(NA_character_)
    body <- str_remove(x, "^beta_")
    hit  <- param_names[str_starts(body, paste0(param_names, "_"))][1]
    if (is.na(hit)) return(NA_character_)
    str_remove(body, paste0("^", hit, "_"))
  }
  
  clean_one <- function(x) {
    cov <- strip_beta_prefix(x)
    if (is.na(cov)) return(x)   # not a beta_ row, or no known parameter matched
    # cov <- clean_cov_labels(cov)
    cov <- str_remove_all(cov,  "_G_[A-Za-z0-9_ ]+") 
    cov <- str_replace_all(cov, "_+", "_")
    cov <- str_remove(cov, "_$")
    cov <- str_remove(cov, "^log")
    cov <- str_remove(cov, "^t(?=[A-Z])")
    # cov <- str_remove(cov, regex("_cov$", ignore_case = TRUE))
    cov <- str_remove_all(cov, "^_+|_+$")
    paste0(prefix, cov)
  }
  
  tbl[[label_col]] <- map_chr(tbl$parameter, clean_one)
  
  if (overwrite) {
    tbl$parameter <- tbl[[label_col]]
    tbl[[label_col]] <- NULL
  }
  
  tbl
}
#' Custom kable styling function for param_table, with RSE and Shrinkage highlightinf depending on threshold
#'
#' @param data dataframe (kbl_obj from build_param_table)
#' @param caption table caption
#' @param cols_to_remove Columns of `data` that are not displayed (default: standard error and 95\% CI bounds).
#' @param labels_mapping mapping of colnames and desired, clean labels
#' @param rse_thresholds list of thresholds for RSE, median and low can be different for fixed and random effect
#' @param shrink_thresholds list of thresholds for shrinkage value
#' @param thresholds_colors list of colours for each flag level Default built-in mapping -> requires minimal arguments but easily overriden if needed
#'
#' @export
kable_param_table <- function(
    data, 
    caption = NULL,
    cols_to_remove = c("se_sa", "ci_low", "ci_high"),
    labels_mapping = c(
      parameter = "Parameter", label = "Parameter", estimate = "Estimate", se_sa = "SE", 
      rse_pct = "RSE (%)", ci_low = "CI 2.5%", ci_high = "CI 97.5%", 
      shrinkage = "Shrinkage (%)*"
    ),
    rse_thresholds = list(high = 100, fe_med = 50, fe_low = 30, re_med = 70, re_low = 50),
    shrink_thresholds = list(high = 100, med = 60, low = 40),
    thresholds_colors = list(high = "firebrick", med = "chocolate3", low = "goldenrod3")
) {
  
  default_color <- "black"
  
  # --- 1. Compute Styles FIRST (While data is still numeric) ---
  rse_colors <- NULL
  rse_bold   <- NULL
  if ("rse_pct" %in% names(data)) {
    # Ensure it's numeric just in case it was converted outside
    rse_num <- as.numeric(data$rse_pct) 
    
    rse_colors <- case_when(
      rse_num > rse_thresholds$high ~ thresholds_colors$high,
      data$type == "Fixed Effects" & rse_num > rse_thresholds$fe_med  ~ thresholds_colors$med,
      data$type == "Fixed Effects" & rse_num > rse_thresholds$fe_low  ~ thresholds_colors$low,
      data$type == "Standard Deviation of the Random Effects" & rse_num > rse_thresholds$re_med ~ thresholds_colors$med,
      data$type == "Standard Deviation of the Random Effects" & rse_num > rse_thresholds$re_low ~ thresholds_colors$low,
      TRUE ~ default_color
    )
    
    rse_bold <- case_when(
      data$type == "Fixed Effects" & rse_num > rse_thresholds$fe_low ~ TRUE,
      data$type == "Standard Deviation of the Random Effects" & rse_num > rse_thresholds$re_low ~ TRUE,
      TRUE ~ FALSE
    )
  }
  
  shrink_colors <- NULL
  shrink_bold   <- NULL
  if ("shrinkage" %in% names(data)) {
    shrink_num <- as.numeric(data$shrinkage)
    
    shrink_colors <- case_when(
      shrink_num > shrink_thresholds$high ~ thresholds_colors$high,
      shrink_num > shrink_thresholds$med  ~ thresholds_colors$med,
      shrink_num > shrink_thresholds$low  ~ thresholds_colors$low,
      TRUE ~ default_color
    )
    
    shrink_bold <- !is.na(shrink_num) & shrink_num > shrink_thresholds$low
  }
  
  # --- 2. Clean up NAs and format strings for display SECOND ---
  display_data <- data |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    mutate(across(everything(), \(x) ifelse(is.na(x) | x == "NA", "\u2014", as.character(x))))
  
  # --- 3. Dynamic Column Filtering & Ordering ---
  all_mapped_cols <- intersect(names(labels_mapping), names(display_data))
  final_raw_cols  <- setdiff(all_mapped_cols, cols_to_remove)
  final_headers   <- labels_mapping[final_raw_cols]
  
  rse_idx    <- which(final_headers == "RSE (%)")
  shrink_idx <- which(final_headers == "Shrinkage (%)*")
  
  # --- 4. Build the Base Kable Object ---
  kbl_obj <- display_data |>
    select(all_of(final_raw_cols)) |>
    kable(
      booktabs  = TRUE,
      col.names = final_headers,
      caption   = caption
    ) |>
    auto_pack_rows(display_data, group_vars = "type") |> 
    column_spec(1, width_min = "4cm") |>
    column_spec(2:length(final_raw_cols), width_min = "1.5cm")
  
  
  # --- 5. Apply Pre-computed Styles Safely ---
  if (!is.null(rse_colors) && length(rse_idx) > 0) {
    kbl_obj <- kbl_obj |> column_spec(column = rse_idx, color = rse_colors, bold = rse_bold)
  }
  
  if (!is.null(shrink_colors) && length(shrink_idx) > 0) {
    kbl_obj <- kbl_obj |> column_spec(column = shrink_idx, color = shrink_colors, bold = shrink_bold)
  }
  
  return(kbl_obj)
}
#' /!\\ Deprecated for master function kable_param_table()
#'
#' apply_pmetrics_styles() : Custom styling function for RSE and Shrinkage depending on threshold
#'
#' @noRd
apply_pmetrics_styles <- function(
    kbl_obj, 
    raw_data, 
    col_names_vector,
    rse_thresholds = list(high = 100, fe_med = 50, fe_low = 30, re_med = 70, re_low = 50),
    shrink_thresholds = list(high = 100, med = 60, low = 40),
    thresholds_colors = list(high = "firebrick", med = "chocolate3", low = "goldenrod3")
) {
  
  # Default parameters
  default_color <- "black"
  
  # --- 1. Dynamic Column Mapping ---
  # Find column positions based on what is currently displayed in kable
  rse_idx    <- which(col_names_vector == "RSE (%)")
  shrink_idx <- which(col_names_vector == "Shrinkage (%)*")
  
  # --- 2. Compute RSE Styling ---
  if ("rse_pct" %in% names(raw_data) && length(rse_idx) > 0) {
    rse_colors <- case_when(
      raw_data$rse_pct > rse_thresholds$high                                                                 ~ thresholds_colors$high,
      raw_data$type == "Fixed Effects" & raw_data$rse_pct > rse_thresholds$fe_med                            ~ thresholds_colors$med,
      raw_data$type == "Fixed Effects" & raw_data$rse_pct > rse_thresholds$fe_low                            ~ thresholds_colors$low,
      raw_data$type == "Standard Deviation of the Random Effects" & raw_data$rse_pct > rse_thresholds$re_med ~ thresholds_colors$med,
      raw_data$type == "Standard Deviation of the Random Effects" & raw_data$rse_pct > rse_thresholds$re_low ~ thresholds_colors$low,
      TRUE ~ default_color
    )
    
    rse_bold <- case_when(
      # data$type == "Information Criteria" | is.na(rse_num) ~ FALSE,
      data$type == "Fixed Effects" & rse_num > rse_thresholds$fe_low ~ TRUE,
      data$type == "Standard Deviation of the Random Effects" & rse_num > rse_thresholds$re_low ~ TRUE,
      TRUE ~ FALSE
    )
    
    
    kbl_obj <- kbl_obj |> 
      column_spec(column = rse_idx, color = rse_colors, bold = rse_bold)
  }
  
  # --- 2. Compute Shrinkage Styles ---
  if ("shrinkage" %in% names(raw_data) && length(shrink_idx) > 0) {
    shrink_colors <- case_when(
      raw_data$shrinkage > shrink_thresholds$high ~ thresholds_colors$high,
      raw_data$shrinkage > shrink_thresholds$med  ~ thresholds_colors$med,
      raw_data$shrinkage > shrink_thresholds$low  ~ thresholds_colors$low,
      TRUE ~ default_color
    )
    
    shrink_bold <- !is.na(raw_data$shrinkage) & raw_data$shrinkage > shrink_thresholds$low
    
    kbl_obj <- kbl_obj |> 
      column_spec(column = shrink_idx, color = shrink_colors, bold = shrink_bold)
  }
  
  
  return(kbl_obj)
}
#' Heatmap of the Fisher correlation matrix
#'
#' Cells with |correlation| >= threshold are marked with a cross.
#' The subtitle turns red when the eigen value ratio exceeds ratio_warn.
#'
#' @param summary_obj output of read_summary(..., full = TRUE)
#' @param diag If TRUE (default), the diagonal is shown and only the lower triangle is drawn.
#' @param threshold absolute correlation flagged as high (default 0.7)
#' @param ratio_warn eigen max/min ratio above which a warning is shown (default 20)
#'
#' @export
plot_fisher_matrix <- function(summary_obj,
                               diag = TRUE, # lower triangle only
                               threshold  = 0.5,
                               ratio_warn = 20) {
  
  mat   <- summary_obj$fisher_matrix
  eigen <- summary_obj$eigen_values
  
  if (is.null(mat)) stop("fisher_matrix is NULL - re-run read_summary() with full = TRUE")
  
  params <- rownames(mat)
  n      <- length(params)
  
  
  # --- Lower triangle only (including diagonal) ---
  long <- mat |>
    rownames_to_column("param_row") |>
    pivot_longer(-param_row, names_to = "param_col", values_to = "corr") |>
    mutate(
      row_idx = match(param_row, params),
      col_idx = match(param_col, params)
    )
  
  if (diag) { # filter lower triangle
    long <- long |>
      filter(row_idx >= col_idx) } 
  
  long <- long |>
    mutate(
      param_row = factor(param_row, levels = rev(params)),
      param_col = factor(param_col, levels = params),
      flag      = abs(corr) >= threshold & row_idx != col_idx
    )
  
  
  # --- Subtitle with eigen value info ---
  subtitle <- if (!is.null(eigen)) {
    warn   <- eigen$ratio > ratio_warn
    status <- if (warn) " \u26a0 high ratio" else " \u2713"
    sprintf("Eigen values \u2014 min: %s  max: %s  ratio: %s%s",
            eigen$eigen_min, eigen$eigen_max, eigen$ratio, status)
  } else {
    ""
  }
  
  subtitle_color <- if (!is.null(eigen) && eigen$ratio > ratio_warn) {
    "firebrick"
  } else {
    "grey40"
  }
  
  thres_legend <- sprintf("|r| \u2265 %s", threshold) # "|r| >= 0.5"
  
  # --- Plot ---
  ggplot(long, aes(x = param_col, y = param_row, fill = corr)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#D6604D",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Correlation"
    ) +
    labs(
      title    = sprintf("Correlation matrix \u2014 %s", summary_obj$model_name),
      subtitle = subtitle,
      x        = NULL,
      y        = NULL
    ) +
    geom_point(
      data        = filter(long, flag),
      aes(x = param_col, y = param_row, color = thres_legend),
      shape       = 4, size = 2, stroke = 1.2,
      inherit.aes = FALSE
    ) +
    # Cross color = black and remove legend title
    scale_color_manual(
      values = c("black"), 
      name   = NULL
    ) +
    # Ensures the legend displays a visible cross (shape = 4).
    guides(
      color = guide_legend(override.aes = list(shape = 4, size = 3, stroke = 1.5))
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y   = element_text(size = 7),
      panel.grid    = element_blank(),
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = subtitle_color)
    )
}
#' Group the rows of a kable table by one or two columns
#'
#' Adds kableExtra row-group headers (`pack_rows()`) to an existing kable. The
#' grouping columns are read from `data`, so they can be left out of the
#' displayed table. With two grouping columns, the second gives italic
#' sub-headers inside each group of the first.
#'
#' @param kbl_obj A kable object (from `knitr::kable()`).
#' @param data The data frame the kable was built from, still holding the
#'   grouping columns; its rows must be in the same order as the table.
#' @param group_vars One or two column names of `data`, outermost group first.
#'
#' @return The kable object with the row-group headers.
#' @export
auto_pack_rows <- function(kbl_obj, data, group_vars) {
  stopifnot(is.data.frame(data))
  stopifnot(length(group_vars) %in% c(1, 2))
  stopifnot(all(group_vars %in% names(data)))
  
  # Add row ID
  df <- data |>
    mutate(.row_id = row_number())
  
  # Level 2 (if any)
  if (length(group_vars) == 2) {
    sub_groups <- df |>
      group_by(across(all_of(group_vars))) |>
      summarise(start = min(.row_id), end = max(.row_id), .groups = "drop")
    
    kbl_obj <- reduce(seq_len(nrow(sub_groups)), .init = kbl_obj, .f = function(kb, i) {
      g <- sub_groups[i, ]
      label <- as.character(g[[group_vars[2]]])
      pack_rows(kb, label, g$start, g$end, bold = FALSE, italic = TRUE)
    })
  }
  
  # Level 1 (mandatory)
  main_groups <- df |>
    group_by(across(all_of(group_vars[1]))) |>
    summarise(start = min(.row_id), end = max(.row_id), .groups = "drop")
  
  kbl_obj <- reduce(seq_len(nrow(main_groups)), .init = kbl_obj, .f = function(kb, i) {
    g <- main_groups[i, ]
    label <- as.character(g[[group_vars[1]]])
    pack_rows(kb, label, g$start, g$end, bold = TRUE)
  })
  
  return(kbl_obj)
}
#' auto_style_cells
#'
#' @noRd
auto_style_cells <- function(kbl_obj, data, target_col, bold_col = NULL, color_col = NULL) {
  # Defensive checks
  stopifnot(is.data.frame(data))
  
  # Find target column index
  if (is.character(target_col)) {
    stopifnot(target_col %in% names(data))
    col_idx <- which(names(data) == target_col)
  } else {
    col_idx <- target_col
  }
  
  # Prepare a list of arguments to pass to column_spec
  spec_args <- list(kable_input = kbl_obj, column = col_idx)
  
  # Process bold formatting
  if (!is.null(bold_col)) {
    stopifnot(bold_col %in% names(data))
    # Replace NAs with 0, then create a TRUE/FALSE vector
    bold_vector <- coalesce(data[[bold_col]], 0)
    spec_args$bold <- (bold_vector == 1)
  }
  
  # Process color formatting
  if (!is.null(color_col)) {
    stopifnot(color_col %in% names(data))
    # Replace NAs with "black" and pass the full vector of colors
    spec_args$color <- coalesce(data[[color_col]], "black")
  }
  
  # Execute column_spec in one single shot with the vectors
  if (!is.null(bold_col) || !is.null(color_col)) {
    kbl_obj <- do.call(column_spec, spec_args)
  }
  
  return(kbl_obj)
}
