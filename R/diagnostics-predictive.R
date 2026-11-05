# Functions: Predictive performances ===========================================
  ## 1. Observations vs Predictions: read, score, plot          ----------------
#' This module reads a model's ObsVsPred ChartsData, computes population and
#'
#' individual prediction-error metrics, and builds a diagnostic scatter plot.
#' It is designed to be called once per model inside a map()/lapply() loop so
#' results from several models can be row-bound into a single comparison
#' table (see obspred_diagnostics() to get the data, metrics table and both
#' plots from a single read_obspred() call).
#' A) Detect the DV name and number of error models from .mlxtran ---------
#' Monolix uses two different conventions for naming ObsVsPred files,
#' depending on whether the model has one or several error models attached
#' to the same observed variable:
#'
#'   - Single error model: [LONGITUDINAL] DEFINITION: directly uses the DV
#'     name, e.g. "CAB_conc = \{distribution=normal, prediction=Cc,
#'     errorModel=combined1(a, b)\}" -> file is "CAB_conc_obsVsPred.txt"
#'
#'   - Multiple error models on the same DV: <DATAFILE> declares
#'     "RPV_conc = \{use=observation, yname=\{'1','2'\}, ...\}" and
#'     [LONGITUDINAL] DEFINITION: uses generic names y1, y2, ... (one per
#'     yname level) -> files are "y1_obsVsPred.txt", "y2_obsVsPred.txt", ...
#'     and must be row-bound under the single DV name "RPV_conc".
#'
#' @return tibble: [dv_name, n_levels, y_labels (list-column)] n_levels > 1 means the generic y1/y2/... file naming applies. Note: Monolix .mlxtran files nest two levels of sections: top-level blocks in angle brackets (<DATAFILE>, <MODEL>, <FIT>, ...) and sub-blocks in square brackets ([FILEINFO], [CONTENT], ...) inside them. The "use=observation" declaration lives in <DATAFILE>'s [CONTENT] sub-block, so the search window below spans the whole <DATAFILE> block.
#'
#' @noRd
detect_dv_from_mlxtran <- function(mlxtran_path) {

  lines <- readLines(mlxtran_path, warn = FALSE)

  datafile_start <- which(str_detect(lines, "^<DATAFILE>"))
  if (length(datafile_start) == 0) stop("No <DATAFILE> section found in .mlxtran.")
  next_sec <- which(str_detect(lines, "^<") & seq_along(lines) > datafile_start[1])[1]
  datafile_lines <- lines[datafile_start[1]:next_sec]

  # Require "observation" to be followed by a comma or closing brace, so
  # "use=observationtype" (a different declaration entirely) is never
  # mistaken for "use=observation".
  obs_lines <- datafile_lines[str_detect(datafile_lines, "use\\s*=\\s*observation\\s*[,}]")]

  if (length(obs_lines) == 0) stop("No 'use=observation' declaration found in <DATAFILE>.")

  map_dfr(obs_lines, function(ln) {

    # ChartsData export filenames (CAB_conc_obsVsPred.txt,
    # concentration_percentiles.txt, ...) are based on Monolix's own
    # "name=" attribute inside the braces, NOT the raw dataset column name
    # on the left of "=". The two happen to be identical in some projects
    # (e.g. "CAB_conc = {use=observation, name=CAB_conc, ...}"), which is
    # why taking the left-hand side alone used to work -- but Monolix's own
    # example projects often name the column generically (e.g. "dv") while
    # name= carries the real observation name used for every export (e.g.
    # "dv = {use=observation, name=concentration, ...}" -> ChartsData files
    # are "concentration_*.txt", not "dv_*.txt"). name= is preferred here
    # whenever present, falling back to the left-hand side only if it
    # genuinely isn't there.
    lhs_name  <- str_match(str_trim(ln), "^(\\w+)\\s*=")[2]
    # (?<![A-Za-z0-9_]) so that "yname=" is never mistaken for "name="
    name_attr <- str_match(ln, "(?<![A-Za-z0-9_])name\\s*=\\s*['\"]?([A-Za-z0-9_.]+)['\"]?\\s*[,}]")[, 2]
    dv_name   <- if (!is.na(name_attr)) name_attr else lhs_name

    yname_match <- str_match(ln, "yname\\s*=\\s*\\{([^}]*)\\}")
    y_labels <- if (!is.na(yname_match[1])) {
      str_extract_all(yname_match[2], "'([^']*)'")[[1]] |> str_remove_all("'")
    } else {
      character(0)
    }

    tibble(dv_name  = dv_name,
          n_levels = max(length(y_labels), 1L),
          y_labels = list(if (length(y_labels) > 0) y_labels else dv_name))
  })
}
#' B) Read the obsVsPred file(s) for one model, with prediction-error   ---
#'
#'    metrics already computed                                        ---
#' Locates ChartsData/ObservationsVsPredictions/, resolves the DV name (from
#' the .mlxtran unless dv overrides it), and reads either:
#'   - "<dv>_obsVsPred.txt"                          (single error model), or
#'   - "y1_obsVsPred.txt", "y2_obsVsPred.txt", ...    (one per declared level)
#' stacking the latter into one tibble with a Y_TYPE column recording which
#' file each row came from, and renaming "y<k>_*" columns to "<dv>_*" so
#' downstream code only ever has to deal with one DV-named column set. Also
#' relabels Y_TYPE / split (via apply_labels(), Shared helpers), resolves
#' indPred (mean or mode, with a per-Y_TYPE fallback when the requested one
#' is unavailable), and computes the prediction-error columns (pPE/ApPE,
#' iPE/AiPE) that obspred_metrics() / plot_obspred() / plot_aipe() below all
#' read from -- so those three (and obspred_diagnostics(), which calls all
#' of them at once) never repeat this logic, only read_obspred() does the
#' real work, following the same read_*() / plot_*() split as the rest of
#' this project. model_name, split (whether split-faceting was requested)
#' and Monolix's own spline overlay (visual_guides, from *_visualGuides.txt)
#' are carried as attributes for the plotting functions to reuse without
#' re-deriving them.
#'
#' @param model_dir Monolix model output folder (same one used elsewhere in this module, containing summary.txt / *.mlxtran / ChartsData)
#' @param dv observed-variable name to use. If NULL (default), it is auto-detected from .mlxtran's <DATAFILE> "use=observation" declaration. Pass it explicitly to skip detection or to disambiguate when a model declares more than one DV.
#' @param ipred_type "mode" (default) or "mean" -- which individual prediction column to prefer as `indPred`.
#' @param split FALSE (default). TRUE relabels/keeps the "split" column (via split_labels) as an extra facet dimension for the downstream plots/metrics; errors if no "split" column exists in the data.
#' @param y_type_labels named lookup vectors passed to apply_labels() (Shared helpers) for renaming/ordering Y_TYPE / split.
#' @param split_labels named lookup vectors passed to apply_labels() (Shared helpers) for renaming/ordering Y_TYPE / split.
#'
#' @return tibble with at least: ID, time, <dv>, popPred, indivPredMean and/or indivPredMode (whichever Monolix exported), Y_TYPE, DV, observed, indPred, ipred_type (per row, reflecting any per-Y_TYPE fallback), pPE, ApPE, pPE_abs, iPE, AiPE, iPE_abs. Every column actually present on disk is preserved as-is -- including model-specific extras such as OCC_Admin / OCC_Route, which are not always present and are not normalized to a common name.
#'
#' @export
read_obspred <- function(model_dir, dv = NULL, ipred_type = c("mode", "mean"),
                         split = FALSE, y_type_labels = NULL, split_labels = NULL) {
  
  ipred_type <- match.arg(ipred_type)
  charts_dir <- file.path(model_dir, "ChartsData", "ObservationsVsPredictions")
  if (!dir.exists(charts_dir)) {
    stop("ChartsData/ObservationsVsPredictions not found in: ", model_dir)
  }
  
  # Resolve DV name / expected file layout (see resolve_dv_layout())
  layout   <- resolve_dv_layout(model_dir, charts_dir, "obsVsPred", dv = dv)
  dv_name  <- layout$dv_name
  n_levels <- layout$n_levels

  # ---- Single error model: read "<dv>_obsVsPred.txt" directly ----------------
  if (n_levels == 1) {
    
    direct_file <- file.path(charts_dir, paste0(dv_name, "_obsVsPred.txt"))
    if (!file.exists(direct_file)) {
      # Fall back to "y1_obsVsPred.txt" if that's what's on disk
      direct_file <- file.path(charts_dir, "y1_obsVsPred.txt")
    }
    if (!file.exists(direct_file)) {
      stop("Could not find an obsVsPred file for DV '", dv_name, "' in ", charts_dir)
    }
    
    out <- read.table(direct_file, sep = ",", comment.char = "", header = TRUE,
                      stringsAsFactors = FALSE) |>
      as_tibble()
    
    # If the file was the generic "y1_..." fallback, rename its y1_* columns
    if (!str_detect(basename(direct_file), stringr::fixed(dv_name))) {
      out <- out |> rename_with(~ str_replace(.x, "^y1", dv_name))
    }
    
    out <- out |>
      mutate(Y_TYPE = dv_name, DV = dv_name)
    
  } else {
    
    # ---- Multiple error models: read y1_*, y2_*, ... and stack -------------
    out <- map_dfr(seq_len(n_levels), function(k) {
      
      f <- file.path(charts_dir, sprintf("y%d_obsVsPred.txt", k))
      if (!file.exists(f)) {
        warning("Expected file not found, skipped: ", f, call. = FALSE)
        return(NULL)
      }
      
      read.table(f, sep = ",", comment.char = "", header = TRUE,
                 stringsAsFactors = FALSE) |>
        as_tibble() |>
        mutate(Y_TYPE = sprintf("y%d", k)) |>
        rename_with(~ str_replace(.x, sprintf("^y%d", k), dv_name))
    }) |>
      mutate(DV = dv_name)
  }
  
  id_col   <- intersect(c("ID", "id"), colnames(out))[1]
  time_col <- intersect(c("time", "TIME", "Time"), colnames(out))[1]
  
  if (!dv_name %in% colnames(out)) {
    # The observed values are normally in a column named after the DV; if
    # not, accept a single generic candidate (y / dv / obs / observation).
    cand <- grep("^(y[0-9]*|dv|obs|observation|observed)$", colnames(out), ignore.case = TRUE, value = TRUE)
    if (length(cand) != 1) {
      stop("Column '", dv_name, "' (observed values) not found in the obsVsPred data. Columns present: ",
          paste(colnames(out), collapse = ", "), call. = FALSE)
    }
    message("Observed values read from column '", cand, "' (no column named '", dv_name, "').")
    out[[dv_name]] <- out[[cand]]
  }

  out <- out |>
    arrange(across(any_of(c(id_col, time_col)))) |>
    mutate(observed = .data[[dv_name]])
  
  # ---- Visual guides (Monolix spline + prediction interval) ------------------
  # Same file-naming convention as obsVsPred itself: "<dv>_visualGuides.txt"
  # for a single error model, or "y1_visualGuides.txt", "y2_visualGuides.txt"
  # ... stacked under Y_TYPE for multiple error models. Read here (using the
  # dv_name / n_levels already resolved above) so plot_obspred() can plot
  # Monolix's own spline instead of refitting a loess smoother.
  read_one_guide <- function(file, y_type) {
    if (!file.exists(file)) return(NULL)
    read.table(file, sep = ",", comment.char = "", header = TRUE,
               stringsAsFactors = FALSE) |>
      as_tibble() |>
      mutate(Y_TYPE = y_type)
  }
  
  visual_guides <- if (n_levels == 1) {
    f <- file.path(charts_dir, paste0(dv_name, "_visualGuides.txt"))
    if (!file.exists(f)) f <- file.path(charts_dir, "y1_visualGuides.txt")
    read_one_guide(f, dv_name)
  } else {
    map_dfr(seq_len(n_levels), function(k) {
      read_one_guide(file.path(charts_dir, sprintf("y%d_visualGuides.txt", k)),
                     sprintf("y%d", k))
    })
  }
  
  model_name <- basename(model_dir)
  
  # Relabel Y_TYPE / split values via the apply_labels() helper defined
  # above (the order in which `labels` is written also controls facet /
  # table row order).
  out$Y_TYPE <- apply_labels(out$Y_TYPE, y_type_labels)
  if (!is.null(visual_guides) && "Y_TYPE" %in% colnames(visual_guides)) {
    visual_guides$Y_TYPE <- apply_labels(visual_guides$Y_TYPE, y_type_labels)
  }
  
  if (split) {
    if (!"split" %in% colnames(out)) {
      warning("split = TRUE was requested but no 'split' column was found in the ",
             "obsVsPred data for ", model_dir, " -- continuing without stratification.",
             call. = FALSE)
      split <- FALSE
    } else if (n_distinct(out$split) <= 1) {
      split <- FALSE   # a single stratum ("All"): nothing to separate
    } else {
      out$split <- apply_labels(out$split, split_labels)
      if (!is.null(visual_guides) && "split" %in% colnames(visual_guides)) {
        visual_guides$split <- apply_labels(visual_guides$split, split_labels)
      }
    }
  }

  # ---- Resolve indivPred column per Y_TYPE, with fallback -------------------
  mean_col <- "indivPredMean"
  mode_col <- "indivPredMode"
  has_mean <- mean_col %in% colnames(out)
  has_mode <- mode_col %in% colnames(out)
  
  if (!has_mean && !has_mode) {
    stop("Neither '", mean_col, "' nor '", mode_col,
         "' found in the obsVsPred data for ", model_dir)
  }
  
  requested_col   <- if (ipred_type == "mode") mode_col else mean_col
  fallback_col    <- if (ipred_type == "mode") mean_col else mode_col
  requested_avail <- if (ipred_type == "mode") has_mode else has_mean
  fallback_avail  <- if (ipred_type == "mode") has_mean else has_mode
  
  # A column can exist but be entirely NA for a given Y_TYPE (this is
  # exactly the "double error model, one level less complete" case you
  # described) -- so the fallback is resolved per Y_TYPE, not globally.
  out <- out |>
    group_by(Y_TYPE) |>
    mutate(
      .requested_ok = requested_avail && !all(is.na(.data[[if (requested_avail) requested_col else mode_col]])),
      indPred    = if (requested_avail) .data[[requested_col]] else NA_real_,
      indPred    = if_else(.requested_ok, indPred,
                           if (fallback_avail) .data[[fallback_col]] else NA_real_),
      ipred_type = if_else(.requested_ok, ipred_type,
                           if_else(fallback_avail, setdiff(c("mean", "mode"), ipred_type), NA_character_))
    ) |>
    ungroup() |>
    select(-.requested_ok)
  
  fallback_msg <- out |>
    distinct(Y_TYPE, ipred_type) |>
    filter(ipred_type != !!ipred_type)
  
  if (nrow(fallback_msg) > 0) {
    message("Note: '", requested_col, "' unavailable for Y_TYPE = ",
            paste(fallback_msg$Y_TYPE, collapse = ", "),
            " -- falling back to '", fallback_col, "' for those rows.")
  }
  
  if (any(is.na(out$indPred))) {
    n_missing <- sum(is.na(out$indPred))
    warning(n_missing, " row(s) have no individual prediction available ",
            "(neither mean nor mode) and will be dropped from individual metrics.",
            call. = FALSE)
  }
  
  # ---- Prediction errors: absolute (PE) and relative (%PE) ------------------
  out <- out |>
    mutate(
      pPE     = (popPred - observed) / observed,
      ApPE    = abs(pPE),
      pPE_abs = popPred - observed,
      iPE     = (indPred - observed) / observed,
      AiPE    = abs(iPE),
      iPE_abs = indPred - observed
    )
  
  attr(out, "visual_guides") <- visual_guides
  attr(out, "model_name")    <- model_name
  attr(out, "split")         <- split
  out
}
#' C) Prediction-error metrics table, obs-vs-pred plot, AiPE plot,      ---
#'
#'    and a master function returning all three plus the data          ---
#' Internal (not meant to be called directly -- use obspred_metrics() /
#' plot_obspred() / plot_aipe() / obspred_diagnostics() below): computes the
#' metrics table from a read_obspred() tibble. Grouping columns: DV x Y_TYPE
#' always, plus split when stratification was requested. RMSE computed by
#' hand (sqrt(mean((pred-obs)^2))) to avoid depending on the Metrics package
#' for a one-line formula.
#'
#' @noRd
.obspred_metrics_table <- function(obs_data, split) {

  # effective split flag as resolved by read_obspred() (it drops split when the
  # data has no split column or a single stratum)
  split <- isTRUE(attr(obs_data, "split"))
  
  group_cols <- if (split) c("DV", "Y_TYPE", "split") else c("DV", "Y_TYPE")
  
  compute_metrics_block <- function(df, pe_col, pe_abs_col, ape_col, pred_col, prediction_label) {
    df |>
      group_by(across(all_of(group_cols))) |>
      summarise(
        N          = sum(!is.na(.data[[pe_col]])),
        ipred_type = if (prediction_label == "Individual predictions") first(ipred_type) else NA_character_,
        MPE        = round(mean(.data[[pe_col]],  na.rm = TRUE) * 100, 1),
        MDPE       = round(median(.data[[pe_col]], na.rm = TRUE) * 100, 1),
        RMSPE      = round(sqrt(mean(.data[[pe_col]]^2, na.rm = TRUE)) * 100, 1),
        RMSE       = round(sqrt(mean(.data[[pe_abs_col]]^2, na.rm = TRUE)), 2),
        MAE        = round(mean(.data[[ape_col]],  na.rm = TRUE) * 100, 1),
        MDAPE      = round(median(.data[[ape_col]], na.rm = TRUE) * 100, 1),
        F10        = round(mean(.data[[ape_col]] < 0.10, na.rm = TRUE) * 100, 1),
        F20        = round(mean(.data[[ape_col]] < 0.20, na.rm = TRUE) * 100, 1),
        F30        = round(mean(.data[[ape_col]] < 0.30, na.rm = TRUE) * 100, 1),
        .groups = "drop"
      ) |>
      mutate(Prediction = prediction_label)
  }
  
  metrics_pop <- compute_metrics_block(obs_data, "pPE", "pPE_abs", "ApPE", "popPred",
                                       "Population predictions")
  metrics_ind <- compute_metrics_block(obs_data |> filter(!is.na(indPred)),
                                       "iPE", "iPE_abs", "AiPE", "indPred",
                                       "Individual predictions")
  
  bind_rows(metrics_pop, metrics_ind) |>
    relocate(DV, Y_TYPE, any_of("split"), Prediction, N, ipred_type)
}
#' Internal: the Observed vs Predicted plot from a read_obspred() tibble
#'
#' Axes follow the classic PopPK convention: predicted on x, observed on y.
#' Facet rows are ordered Population first, Individual second (matches
#' Monolix's own panel order). Facet columns are Y_TYPE and/or split, but
#' only when they actually vary -- a facet column would otherwise just
#' repeat the same single label. The smoother is Monolix's own spline
#' (popPred_spline / indivPred_spline from *_visualGuides.txt, carried as an
#' attribute by read_obspred()), not a locally refit loess -- so the curve
#' shown is the one Monolix actually computed for this model, falling back
#' to loess only if no visualGuides file was found.
#'
#' @noRd
.plot_obspred_from_data <- function(obs_data, split, base_size = 9) {

  split <- isTRUE(attr(obs_data, "split"))   # effective flag, see read_obspred()
  
  model_name    <- attr(obs_data, "model_name")
  visual_guides <- attr(obs_data, "visual_guides")
  pred_levels   <- c("Population prediction", "Individual prediction")
  
  facet_cols <- character(0)
  if (n_distinct(obs_data$Y_TYPE) > 1) facet_cols <- c(facet_cols, "Y_TYPE")
  if (split) facet_cols <- c(facet_cols, "split")
  
  plot_cols <- c(facet_cols, "observed", "popPred", "indPred")
  
  plot_data <- obs_data |>
    select(all_of(plot_cols)) |>
    pivot_longer(cols = c(popPred, indPred),
                 names_to = "pred_type", values_to = "predicted") |>
    mutate(pred_type = recode(pred_type,
                              popPred = "Population prediction",
                              indPred = "Individual prediction"),
           pred_type = factor(pred_type, levels = pred_levels))
  
  axis_max <- max(c(plot_data$observed, plot_data$predicted), na.rm = TRUE)
  
  spline_data <- NULL
  if (!is.null(visual_guides) && nrow(visual_guides) > 0) {
    
    guide_cols <- intersect(facet_cols, colnames(visual_guides))
    
    # If visual_guides varies by a column not in facet_cols -- typically
    # "split", when split = FALSE but this model's ChartsData export was
    # itself stratified (e.g. by route) in the Monolix interface -- using
    # Monolix's own spline as-is would draw several distinct per-stratum
    # curves on top of each other in the same panel. Fall back to a
    # freshly fit smoother on the (pooled) points instead of overlaying
    # mismatched curves, and say so once.
    unmatched_variation <- setdiff(intersect(c("split", "Y_TYPE"), colnames(visual_guides)), facet_cols)
    guide_mismatch <- length(unmatched_variation) > 0 &&
      any(vapply(unmatched_variation, function(col) n_distinct(visual_guides[[col]]) > 1, logical(1)))
    
    if (guide_mismatch) {
      message("Note: Monolix's own spline varies by ", paste(unmatched_variation, collapse = ", "),
              ", which isn't part of the current facet (split = FALSE?) -- falling back to a ",
              "freshly fit smoother instead of overlaying mismatched curves. Set split = TRUE to ",
              "see Monolix's own spline per stratum instead.")
    } else {
      spline_data <- bind_rows(
        visual_guides |>
          transmute(across(all_of(guide_cols)),
                    predicted = popPred_spline_abscissa, observed = popPred_spline,
                    pred_type = factor("Population prediction", levels = pred_levels)),
        visual_guides |>
          transmute(across(all_of(guide_cols)),
                    predicted = indivPred_spline_abscissa, observed = indivPred_spline,
                    pred_type = factor("Individual prediction", levels = pred_levels))
      ) |>
        filter(!is.na(predicted), !is.na(observed))
    }
  }
  
  p <- ggplot(plot_data, aes(x = predicted, y = observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.6) +
    geom_point(alpha = 0.45, size = 1.6, color = "#2166AC")
  
  if (!is.null(spline_data) && nrow(spline_data) > 0) {
    p <- p + geom_line(data = spline_data, aes(x = predicted, y = observed),
                       color = "#B2182B", linewidth = 0.8)
  } else {
    # method = "gam" (mgcv, ships with base R) rather than "loess": loess
    # is prone to "zero-width neighborhood" / near-singularity warnings on
    # the clustered, sparse-sampling x-ranges typical of PK data; a
    # penalized spline handles that far more gracefully with no extra
    # dependency.
    p <- p + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE,
                         color = "#B2182B", linewidth = 0.8)
  }
  
  facet_formula <- if (length(facet_cols) > 0) {
    as.formula(paste("pred_type ~", paste(facet_cols, collapse = " + ")))
  } else {
    pred_type ~ .
  }
  
  p +
    coord_equal(xlim = c(0, axis_max), ylim = c(0, axis_max)) +
    facet_grid(facet_formula) +
    labs(title = sprintf("Observed vs Predicted \u2014 %s", model_name),
         x = "Predicted", y = "Observed") +
    theme_bw(base_size = base_size) +
    theme(strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 12))
}
#' Internal: |AiPE| (individual absolute prediction error) vs observed
#'
#' concentration, from a read_obspred() tibble. Same facet columns as
#' .plot_obspred_from_data() (Y_TYPE / split, only when they vary), but not
#' split by prediction type since AiPE only concerns individual predictions.
#' Censored (BLQ) observations are highlighted the same way as the NPDE
#' scatter plots below (layer_censored_points(), Shared helpers).
#'
#' @noRd
.plot_aipe_from_data <- function(obs_data, split, highlight_censored) {

  split <- isTRUE(attr(obs_data, "split"))   # effective flag, see read_obspred()
  
  model_name <- attr(obs_data, "model_name")
  
  facet_cols <- character(0)
  if (n_distinct(obs_data$Y_TYPE) > 1) facet_cols <- c(facet_cols, "Y_TYPE")
  if (split) facet_cols <- c(facet_cols, "split")
  
  aipe_data <- obs_data |> filter(!is.na(AiPE))
  if ("censored" %in% colnames(aipe_data)) {
    aipe_data$censored_label <- factor(aipe_data$censored, levels = c(0, 1),
                                       labels = c("Above", "Below"))
  }
  
  p <- ggplot(aipe_data, aes(x = observed, y = AiPE)) +
    layer_censored_points(aipe_data, highlight_censored) +
    # method = "gam" rather than "loess": far more robust on sparse/
    # clustered PK x-ranges, no extra dependency (mgcv ships with base R).
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE,
                color = "black", linewidth = 0.5) +
    scale_y_continuous(transform = "log1p", labels = scales::label_percent()) +
    labs(title = sprintf("Absolute individual prediction error vs concentration \u2014 %s", model_name),
         x = "Observed concentration", y = "|Individual PE|",
         fill = "Quantification limit") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12), legend.position = "bottom")
  
  if (length(facet_cols) > 0) {
    p <- p + facet_wrap(as.formula(paste("~", paste(facet_cols, collapse = " + "))), labeller = label_both)
  }
  p
}
#' Prediction-error metrics table for one model
#'
#' Reads a model's obsVsPred data (read_obspred()) and summarises
#' population and individual prediction errors (MPE, MDPE, RMSPE, RMSE,
#' MAE, MDAPE, F10/F20/F30) by DV x Y_TYPE (and split, if requested).
#'
#' @param model_dir,dv,ipred_type,split,y_type_labels,split_labels
#'   same conventions as read_obspred()
#' @return a tibble
#' @export
obspred_metrics <- function(model_dir, dv = NULL, ipred_type = c("mode", "mean"),
                            split = FALSE, y_type_labels = NULL, split_labels = NULL) {
  obs_data <- read_obspred(model_dir, dv = dv, ipred_type = ipred_type, split = split,
                           y_type_labels = y_type_labels, split_labels = split_labels)
  .obspred_metrics_table(obs_data, split = split)
}
#' Bold the acceptable values of a metrics table, colour the others red
#'
#' Applies column_spec() styling to a kableExtra table already built from
#' obspred_metrics()'s (or obspred_diagnostics()$metrics's) output. Columns
#' are matched by NAME against `thresholds`, so any subset of the usual
#' MPE/MDPE/RMSPE/RMSE/MAE/MDAPE/F10/F20/F30 columns can be present -- one
#' that's missing from `data` is silently skipped rather than erroring.
#'
#' A value is "acceptable" when it meets the threshold: |MPE| and |MDPE|
#' at or below their cutoff, RMSPE/MAE/MDAPE at or below theirs, and
#' F10/F20/F30 (fraction of predictions within X%) at or above theirs.
#' Acceptable values are shown in bold; values that miss the threshold are
#' shown in `flagged_color` (not bold); missing values are left plain.
#' Defaults are commonly used PopPK rule-of-thumb cutoffs, not a universal
#' standard -- override via `thresholds` for your own project's criteria.
#'
#' @param kbl_obj        a kableExtra table already built from `data`
#' @param data           the tibble the table was built from -- column
#'                       positions are looked up here, so it must have the
#'                       same columns, in the same order, as what was
#'                       passed to kable() (if you select(-Prediction)
#'                       before kable(), pass that same reduced tibble)
#' @param thresholds     named list overriding any default cutoff (all in
#'                       %): MPE = 15, MDPE = 15, RMSPE = 30, MAE = 30,
#'                       MDAPE = 30, F10 = 50, F20 = 80, F30 = 90
#' @param flagged_color  colour of values that miss their threshold
#'                       ("firebrick")
#'
#' @return kbl_obj with column_spec() styling applied
#' @export
style_metrics_table <- function(kbl_obj, data, thresholds = list(), flagged_color = "firebrick") {

  default_thresholds <- list(
    MPE = 15, MDPE = 15, RMSPE = 30, MAE = 30, MDAPE = 30,
    F10 = 50, F20 = 80, F30 = 90
  )
  thresholds <- modifyList(default_thresholds, thresholds)

  lower_is_bad <- c("F10", "F20", "F30")   # fractions: acceptable when AT OR ABOVE threshold
  abs_flagged  <- c("MPE", "MDPE")         # bias: acceptable on |value|

  for (col in names(thresholds)) {
    if (!col %in% colnames(data)) next   # silently skip columns not present

    idx    <- which(colnames(data) == col)
    values <- data[[col]]
    cutoff <- thresholds[[col]]

    acceptable <- if (col %in% lower_is_bad) {
      values >= cutoff
    } else if (col %in% abs_flagged) {
      abs(values) <= cutoff
    } else {
      values <= cutoff
    }
    acceptable[is.na(acceptable)] <- FALSE
    missed <- !acceptable & !is.na(values)

    kbl_obj <- kbl_obj |>
      column_spec(idx, bold = acceptable, color = if_else(missed, flagged_color, "grey30"))
  }

  kbl_obj
}
#' Observed vs Predicted plot for one model
#'
#' @param model_dir,dv,ipred_type,split,y_type_labels,split_labels
#'   same conventions as read_obspred()
#' @param base_size Base font size of the plots (ggplot2 theme).
#' @return a ggplot object
#' @export
plot_obspred <- function(model_dir, dv = NULL, ipred_type = c("mode", "mean"),
                         split = FALSE, y_type_labels = NULL, split_labels = NULL,
                         base_size = 9) {
  obs_data <- read_obspred(model_dir, dv = dv, ipred_type = ipred_type, split = split,
                           y_type_labels = y_type_labels, split_labels = split_labels)
  .plot_obspred_from_data(obs_data, split = split, base_size = base_size)
}
#' |AiPE| (individual prediction error) vs concentration plot for one model
#'
#' @param model_dir,dv,ipred_type,split,y_type_labels,split_labels
#'   same conventions as read_obspred()
#' @param highlight_censored  highlight BLQ observations (see
#'   layer_censored_points(), Shared helpers). TRUE by default.
#' @return a ggplot object
#' @export
plot_aipe <- function(model_dir, dv = NULL, ipred_type = c("mode", "mean"),
                      split = FALSE, y_type_labels = NULL, split_labels = NULL,
                      highlight_censored = TRUE) {
  obs_data <- read_obspred(model_dir, dv = dv, ipred_type = ipred_type, split = split,
                           y_type_labels = y_type_labels, split_labels = split_labels)
  .plot_aipe_from_data(obs_data, split = split, highlight_censored = highlight_censored)
}
#' All of the above at once, from a single read_obspred() call
#'
#' Convenience master function for when you want the data, the metrics
#' table, and both plots together without read_obspred() running three
#' times over (once per obspred_metrics() / plot_obspred() / plot_aipe()
#' call) -- each of those remains usable on its own when you only need one
#' piece.
#'
#' @param model_dir,dv,ipred_type,split,y_type_labels,split_labels,highlight_censored  same conventions as read_obspred() / plot_aipe()
#' @param base_size Base font size of the plots (ggplot2 theme).
#' @return list(model_name, data, metrics, obspred, aipe)
#' @export
obspred_diagnostics <- function(model_dir, dv = NULL, ipred_type = c("mode", "mean"),
                                split = FALSE, y_type_labels = NULL, split_labels = NULL,
                                highlight_censored = TRUE, base_size = 9) {
  obs_data <- read_obspred(model_dir, dv = dv, ipred_type = ipred_type, split = split,
                           y_type_labels = y_type_labels, split_labels = split_labels)
  list(model_name = attr(obs_data, "model_name"),
       data       = obs_data,
       metrics    = .obspred_metrics_table(obs_data, split = split),
       obspred    = .plot_obspred_from_data(obs_data, split = split, base_size = base_size),
       aipe       = .plot_aipe_from_data(obs_data, split = split, highlight_censored = highlight_censored))
}





  ## 2. Visual Predictive Check & NPDE                          ----------------
#' This module reads a model's VisualPredictiveCheck, ScatterPlotOfTheResiduals
#'
#' and DistributionOfTheResiduals ChartsData, and builds the corresponding
#' diagnostic plots. It reuses the same DV / Y_TYPE auto-detection as section 1
#' (resolve_mlxtran_path() / detect_dv_from_mlxtran()), and adds one extra
#' layer of stratification on top: Monolix's own "split" column (e.g.
#' "#OCC_Route:1", or "#OCC_Route:1#BMI_J0_cat:25-30" when several splits are
#' combined), which parse_split() (Shared helpers) unpacks into one plain
#' column per key so it can be filtered / faceted on like any other variable.
#'
#' plot_vpc() and plot_npde() share the same four conventions, implemented
#' once in Shared helpers (relabel_columns(), resolve_facets(), add_facets(),
#' title_with_filter()) so both functions -- and any future one -- apply them
#' identically:
#'   - filter=     an unquoted dplyr-style expression (e.g. OCC_Route == "IM"),
#'                 evaluated against the Y_TYPE / split-derived columns. There
#'                 is no separate mini-language to learn -- reference whatever
#'                 column names parse_split() produced for this model (i.e.
#'                 whatever you named the split in Monolix). NULL (default) =
#'                 no filter. Matching zero rows is an error, not a silent
#'                 empty plot (see resolve_facets()).
#'   - facet_by / facet_max_levels =  "auto" (default) facets on every Y_TYPE /
#'                 split-derived column that still varies after filter= AND
#'                 has at most facet_max_levels (6 by default) distinct
#'                 values -- a higher-cardinality column (e.g. an
#'                 administration occasion index) is left out of the
#'                 auto-selection instead of exploding into dozens of panels.
#'                 NULL disables faceting; a character vector forces/restricts
#'                 which columns are used regardless of cardinality. See
#'                 resolve_facets() (Shared helpers) for the full logic.
#'   - value_labels = a named list of named character vectors used to
#'                 rename/reorder Y_TYPE and/or split-derived columns via
#'                 apply_labels() (Shared helpers), e.g.
#'                 list(OCC_Route = c("1" = "Oral", "2" = "IM")). Relabeling
#'                 is applied BEFORE filter=, so filter= should reference the
#'                 relabeled values when value_labels renames a column (e.g.
#'                 filter = OCC_Route == "IM", not == "2").
#'   - x_trans / y_trans, xlab / ylab / title -- axis scale transforms (passed
#'                 to scale_*_continuous(transform = )) and labels/title, with
#'                 sensible defaults built from the model name (and the
#'                 filter= expression, if any).
#' A) Read and stack VPC / residuals / density ChartsData for one model ----
#' Generic reader shared by read_vpc(), read_residuals() and
#' read_residuals_pdf() below: resolves the DV name / y1, y2, ... layout
#' exactly like read_obspred() above, then reads either
#' "<dv>_<suffix>.txt" or "y1_<suffix>.txt", "y2_<suffix>.txt", ... stacking
#' the latter under a Y_TYPE column. Unlike obsVsPred, column names in these
#' three ChartsData exports do NOT carry a "y<k>_" / "<dv>_" prefix (e.g.
#' "bins_middles", "npde_simBlq" are the same whatever the DV/Y_TYPE), so no
#' renaming step is needed. When a "split" column is present, it is left
#' untouched and additionally unpacked via parse_split() (Shared helpers)
#' into its component columns.
#'
#' @noRd
read_chartsdata_stacked <- function(model_dir, subdir, suffix, dv = NULL) {
  
  charts_dir <- file.path(model_dir, "ChartsData", subdir)
  if (!dir.exists(charts_dir)) {
    stop("ChartsData/", subdir, " not found in: ", model_dir)
  }
  
  layout   <- resolve_dv_layout(model_dir, charts_dir, suffix, dv = dv)
  dv_name  <- layout$dv_name
  n_levels <- layout$n_levels

  read_one <- function(file, y_type) {
    if (!file.exists(file)) {
      warning("Expected file not found, skipped: ", file, call. = FALSE)
      return(NULL)
    }
    read.table(file, sep = ",", comment.char = "", header = TRUE,
               stringsAsFactors = FALSE) |>
      as_tibble() |>
      mutate(Y_TYPE = y_type)
  }
  
  if (n_levels == 1) {
    f <- file.path(charts_dir, sprintf("%s_%s.txt", dv_name, suffix))
    if (!file.exists(f)) f <- file.path(charts_dir, sprintf("y1_%s.txt", suffix))
    out <- read_one(f, dv_name)
    if (is.null(out)) {
      stop("Could not find a '", suffix, "' file for DV '", dv_name, "' in ", charts_dir)
    }
  } else {
    out <- map_dfr(seq_len(n_levels), function(k) {
      read_one(file.path(charts_dir, sprintf("y%d_%s.txt", k, suffix)), sprintf("y%d", k))
    })
  }
  
  out <- out |> mutate(DV = dv_name)
  
  if ("split" %in% colnames(out)) {
    parsed <- parse_split(out$split)
    # Some exports (e.g. VPC observations) already carry a column named after
    # the split key (a numeric OCC_Route next to split = "#OCC_Route:1"). The
    # parsed character version replaces it, so labels / filters behave the
    # same in every reader and bind_cols() never has to rename duplicates.
    out <- bind_cols(out |> select(-any_of(colnames(parsed))), parsed)
  }
  
  attr(out, "model_name") <- basename(model_dir)
  out
}
#' read_vpc(): ChartsData/VisualPredictiveCheck/*_percentiles.txt -- bin
#'
#' mid-time, empirical 10th/50th/90th percentiles, and the corresponding
#' simulated (theoretical) prediction intervals for those same percentiles.
#'
#' @param model_dir Monolix results folder of the model (ChartsData must have been exported).
#' @param dv Optional name of the observed variable; by default all of them are stacked under a Y_TYPE column.
#'
#' @export
read_vpc <- function(model_dir, dv = NULL) {
  read_chartsdata_stacked(model_dir, "VisualPredictiveCheck", "percentiles", dv = dv)
}
#' read_vpc_observations(): ChartsData/VisualPredictiveCheck/*_observations.txt
#'
#' -- the raw observations behind the VPC, stacked / split-parsed like every
#' other ChartsData reader. Per observation type "<p>" (the DV name, or y1,
#' y2, ... when error models are split) the file holds <p> (raw value),
#' <p>_pc (prediction-corrected), <p>_simBlq and <p>_simBlq_pc (censored
#' values replaced by simulated ones), plus time, timeAfterLastDose,
#' censored, split, ... The value column is chosen per row from its own
#' Y_TYPE, so stacked y1/y2 files work too.
#'   value     : which variant to return as obs_value -- "pc" (default,
#'               matches a prediction-corrected VPC), "raw", "simBlq" or
#'               "simBlq_pc"
#'   value_col : explicit column name overriding `value` for every row
#'   time_col  : "time" (default) or e.g. "timeAfterLastDose"
#' Adds obs_time / obs_value columns.
#'
#' @noRd
read_vpc_observations <- function(model_dir, dv = NULL, value = c("pc", "raw", "simBlq", "simBlq_pc"),
                                  value_col = NULL, time_col = "time") {

  value <- match.arg(value)
  obs   <- read_chartsdata_stacked(model_dir, "VisualPredictiveCheck", "observations", dv = dv)
  cols  <- colnames(obs)

  if (!time_col %in% cols) {
    stop("Time column '", time_col, "' not found in the VPC observations file. Columns present: ",
        paste(cols, collapse = ", "), call. = FALSE)
  }

  if (!is.null(value_col)) {
    if (!value_col %in% cols) {
      stop("value_col '", value_col, "' not found in the VPC observations file. Columns present: ",
          paste(cols, collapse = ", "), call. = FALSE)
    }
    obs$obs_value <- obs[[value_col]]
  } else {
    suffix <- c(raw = "", pc = "_pc", simBlq = "_simBlq", simBlq_pc = "_simBlq_pc")[[value]]
    obs$obs_value <- NA_real_
    for (yt in unique(obs$Y_TYPE)) {
      col <- paste0(yt, suffix)
      if (col %in% cols) {
        rows <- obs$Y_TYPE == yt
        obs$obs_value[rows] <- obs[[col]][rows]
      }
    }
    if (all(is.na(obs$obs_value))) {
      stop("No '<observation>", suffix, "' column found in the VPC observations file for ",
          paste(unique(obs$Y_TYPE), collapse = ", "), ". Columns present: ",
          paste(cols, collapse = ", "), ". Pass value_col= to choose one.", call. = FALSE)
    }
  }

  obs$obs_time <- obs[[time_col]]
  obs
}
#' read_residuals(): ChartsData/ScatterPlotOfTheResiduals/*_residuals.txt --
#'
#' one row per observation, with one (prediction_<type>, <type>,
#' <type>_simBlq) trio per residual family found on disk -- typically pwRes,
#' iwRes_mean, iwRes_mode and npde -- plus "censored" (BLQ flag) and Monolix's
#' own "color" / "filter" export columns. The "_simBlq" variant is the one to
#' plot: for censored (BLQ) observations Monolix computes it by simulation
#' rather than leaving it undefined, unlike the plain (non-simBlq) column.
#'
#' @param model_dir Monolix results folder of the model (ChartsData must have been exported).
#' @param dv Optional name of the observed variable; by default all of them are stacked under a Y_TYPE column.
#'
#' @export
read_residuals <- function(model_dir, dv = NULL) {
  read_chartsdata_stacked(model_dir, "ScatterPlotOfTheResiduals", "residuals", dv = dv)
}
#' read_residuals_pdf(): ChartsData/DistributionOfTheResiduals/*_pdf.txt
#'
#' (empirical density of pwRes / iwRes / npde, per Y_TYPE and split) plus
#' theoreticalGuides.txt (the reference N(0,1) density/cdf curve -- a single
#' file shared by the whole model, not stacked, since it does not depend on
#' Y_TYPE or split).
#'
#' @noRd
read_residuals_pdf <- function(model_dir, dv = NULL) {
  density <- read_chartsdata_stacked(model_dir, "DistributionOfTheResiduals", "pdf", dv = dv)
  theo_file <- file.path(model_dir, "ChartsData", "DistributionOfTheResiduals", "theoreticalGuides.txt")
  theoretical <- if (file.exists(theo_file)) {
    read.table(theo_file, sep = ",", comment.char = "", header = TRUE,
               stringsAsFactors = FALSE) |> as_tibble()
  } else NULL
  list(density = density, theoretical = theoretical, model_name = basename(model_dir))
}
#' Prediction-corrected visual predictive check for one model
#'
#' Reads a model's VisualPredictiveCheck ChartsData (read_vpc()) and plots
#' the empirical 10th/50th/90th percentiles (points + lines) against
#' Monolix's simulated 90% prediction intervals for those same three
#' percentiles (shaded ribbons), against bin mid-time.
#'
#' @param model_dir,dv        same convention as read_obspred()
#' @param filter               unquoted filter expression evaluated on the
#'                             (split-unpacked) data, e.g.
#'                             filter = OCC_Route == "IM". NULL (default) =
#'                             no filter. Errors clearly if it matches zero
#'                             rows (see resolve_facets(), Shared helpers).
#' @param facet_by,facet_max_levels  which columns to facet by --
#'                             "auto" (default), NULL, or a character
#'                             vector -- and the cardinality cap for "auto".
#'                             See resolve_facets() in Shared helpers.
#' @param value_labels         named list of named character vectors used to
#'                             rename/reorder Y_TYPE and/or split-derived
#'                             columns, e.g.
#'                             list(OCC_Route = c("1" = "Oral", "2" = "IM")).
#'                             Applied BEFORE filter=, so filter= should
#'                             reference the relabeled values.
#' @param show_observations    FALSE (default). TRUE adds the raw observations
#'                             (grey, semi-transparent) from the VPC
#'                             observations export, relabeled / filtered
#'                             like the rest. If that file can't be read a
#'                             warning is raised and the VPC is drawn
#'                             without them.
#' @param obs_value            which variant of the observations to plot:
#'                             "pc" (default, prediction-corrected),
#'                             "raw", "simBlq" or "simBlq_pc" (censored
#'                             values replaced by simulated ones)
#' @param obs_value_col        explicit observation column, overriding
#'                             obs_value
#' @param obs_time_col         "time" (default) or e.g. "timeAfterLastDose"
#' @param x_trans,y_trans     passed to scale_x/y_continuous(transform = );
#'                             e.g. y_trans = "log10". NULL = linear.
#' @param xlab,ylab,title    axis / plot titles. NULL (ylab, title) uses
#'                             a sensible default built from the model name
#'                             (and the filter= expression, if any).
#'
#' @param facet_scales Scales of the facets, passed on to facet_wrap(): "fixed", "free", "free_x" or "free_y".
#' @return a ggplot object
#' @export
plot_vpc <- function(model_dir, dv = NULL, filter = NULL,
                     facet_by = "auto", facet_max_levels = 6, facet_scales = "free_x",
                     value_labels = NULL,
                     show_observations = FALSE, obs_value = "pc", obs_value_col = NULL,
                     obs_time_col = "time",
                     x_trans = NULL, y_trans = NULL,
                     xlab = "Time", ylab = NULL, title = NULL) {
  
  filter_quo <- rlang::enquo(filter)
  data       <- read_vpc(model_dir, dv = dv)
  model_name <- attr(data, "model_name")
  
  data <- relabel_columns(data, value_labels)
  
  fixed_cols <- c("bins_middles", "empirical_median", "empirical_lower", "empirical_upper",
                  "theoretical_median_median_pc", "theoretical_median_piLower_pc", "theoretical_median_piUpper_pc",
                  "theoretical_lower_median_pc", "theoretical_lower_piLower_pc", "theoretical_lower_piUpper_pc",
                  "theoretical_upper_median_pc", "theoretical_upper_piLower_pc", "theoretical_upper_piUpper_pc",
                  "split", "DV")
  candidate_facets <- setdiff(colnames(data), fixed_cols)
  
  r <- resolve_facets(data, filter_quo, candidate_facets,
                      facet_by = facet_by, facet_max_levels = facet_max_levels)
  data       <- r$data
  facet_cols <- r$facet_cols
  
  # Optional raw observations (drawn first, so ribbons and lines sit on top)
  p <- ggplot(data, aes(x = bins_middles))
  obs_shown <- FALSE
  if (show_observations) {
    obs <- tryCatch({
      o <- read_vpc_observations(model_dir, dv = dv, value = obs_value,
                                 value_col = obs_value_col, time_col = obs_time_col)
      o <- relabel_columns(o, value_labels)
      if (!rlang::quo_is_null(filter_quo)) o <- dplyr::filter(o, !!filter_quo)
      o
    }, error = function(e) {
      warning("Observations could not be added to the VPC: ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (!is.null(obs) && nrow(obs) > 0) {
      unstratified <- setdiff(facet_cols, colnames(obs))
      if (length(unstratified) > 0) {
        warning("VPC observations carry no '", paste(unstratified, collapse = "', '"),
               "' column: they are shown in every panel.", call. = FALSE)
      }
      p <- p + geom_point(data = obs, aes(x = obs_time, y = obs_value), inherit.aes = FALSE,
                          colour = "grey50", alpha = 0.35, size = 0.9)
      obs_shown <- TRUE
    }
  }
  
  p <- p +
    geom_ribbon(aes(ymin = theoretical_upper_piLower_pc, ymax = theoretical_upper_piUpper_pc),
                fill = "#6494c5", alpha = 0.4) +
    geom_ribbon(aes(ymin = theoretical_lower_piLower_pc, ymax = theoretical_lower_piUpper_pc),
                fill = "#6494c5", alpha = 0.4) +
    geom_ribbon(aes(ymin = theoretical_median_piLower_pc, ymax = theoretical_median_piUpper_pc),
                fill = "#2166ac", alpha = 0.55) +
    geom_line(aes(y = empirical_upper), linewidth = 0.5, linetype = "longdash") +
    geom_point(aes(y = empirical_upper), size = 1.4) +
    geom_line(aes(y = empirical_median), linewidth = 0.5) +
    geom_point(aes(y = empirical_median), size = 1.4) +
    geom_line(aes(y = empirical_lower), linewidth = 0.5, linetype = "longdash") +
    geom_point(aes(y = empirical_lower), size = 1.4)
  
  p <- add_facets(p, facet_cols, scales = facet_scales)
  
  if (!is.null(x_trans)) p <- p + scale_x_continuous(transform = x_trans)
  p <- p + scale_y_continuous(transform = if (is.null(y_trans)) "identity" else y_trans,
                              breaks = scales::pretty_breaks())
  
  default_title <- title_with_filter(
    sprintf("Prediction-corrected visual predictive check \u2014 %s", model_name), filter_quo)
  
  p +
    labs(x = xlab, y = if (is.null(ylab)) "Concentration" else ylab,
         title = if (is.null(title)) default_title else title,
         caption = paste0("Shaded bands: simulated 90% CI of the theoretical 10th/50th/90th percentiles. Points/lines: empirical 10th/50th/90th percentiles (solid/darker = median, dashed/lighter = outer).",
                          if (obs_shown) " Grey points: observations." else "")) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 12),
          plot.caption = element_text(size = 8, color = "grey40"),
          legend.position = "none")
}
#' NPDE / residual diagnostic plots for one model
#'
#' One entry point for four classic residual diagnostics built from a
#' model's ChartsData: residuals vs time, residuals vs population
#' prediction, a Q-Q plot against the theoretical N(0,1), and the residuals'
#' empirical density against that same theoretical curve. Despite the name,
#' it is not limited to npde -- residual_type = lets it read any of the
#' residual families Monolix exports in the same files (pwRes, iwRes_mean,
#' iwRes_mode), since read_residuals() / read_residuals_pdf() above already
#' read all of them at once.
#'
#' @param model_dir,dv        same convention as read_obspred()
#' @param type                 "time", "predictions", "density" or "qq"
#' @param residual_type        "npde" (default), "pwRes", "iwRes_mean" or
#'                             "iwRes_mode"
#' @param highlight_censored   for type = "time" / "predictions" / "qq":
#'                             highlight BLQ observations (white fill) among
#'                             the plotted points, using the "censored"
#'                             column carried by the residuals file. TRUE by
#'                             default; ignored for type = "density".
#' @param filter               unquoted filter expression, e.g.
#'                             filter = OCC_Route == "IM". NULL (default) =
#'                             no filter. Errors clearly if it matches zero
#'                             rows (see resolve_facets(), Shared helpers).
#' @param facet_by,facet_max_levels  which columns to facet by --
#'                             "auto" (default), NULL, or a character
#'                             vector -- and the cardinality cap for "auto"
#'                             (e.g. keeps a 2-level OCC_Route auto-faceted
#'                             while leaving a many-level OCC_Admin
#'                             occasion index out unless asked for). See
#'                             resolve_facets() in Shared helpers.
#' @param value_labels         named list of named character vectors, as in
#'                             plot_vpc(). Applied BEFORE filter=.
#' @param x_trans,y_trans     as in plot_vpc()
#' @param xlab,ylab,title    axis / plot titles. NULL uses a sensible
#'                             default built from the model name / type.
#'
#' @return a ggplot object
#' @export
plot_npde <- function(model_dir, dv = NULL,
                      type = c("time", "predictions", "density", "qq"),
                      residual_type = "npde", highlight_censored = TRUE,
                      filter = NULL, facet_by = "auto", facet_max_levels = 6,
                      value_labels = NULL,
                      x_trans = NULL, y_trans = NULL,
                      xlab = NULL, ylab = NULL, title = NULL) {
  
  type       <- match.arg(type)
  filter_quo <- rlang::enquo(filter)
  model_name <- basename(model_dir)
  
  # ---- Density: distribution vs theoretical N(0,1) --------------------------
  if (type == "density") {
    
    rd          <- read_residuals_pdf(model_dir, dv = dv)
    data        <- relabel_columns(rd$density, value_labels)
    theoretical <- rd$theoretical
    
    abscissa_col <- paste0(residual_type, "_abscissa")
    pdf_col      <- paste0(residual_type, "_pdf")
    if (!all(c(abscissa_col, pdf_col) %in% colnames(data))) {
      stop("residual_type = '", residual_type, "' has no matching column ('",
           abscissa_col, "' / '", pdf_col, "') in the density file.")
    }
    
    density_cols     <- colnames(data)[str_detect(colnames(data), "_abscissa$|_pdf$")]
    candidate_facets <- setdiff(colnames(data), c(density_cols, "split", "DV"))
    
    r <- resolve_facets(data, filter_quo, candidate_facets,
                        facet_by = facet_by, facet_max_levels = facet_max_levels)
    data       <- r$data
    facet_cols <- r$facet_cols
    
    # Bar width = 80% of the (per-group) step between successive abscissa
    # points, so the histogram-style bars don't overlap or leave gaps when
    # the binning isn't perfectly regular. The last bin per group has no
    # "next" point to measure a gap from -- carry the previous gap forward.
    data <- data |>
      group_by(across(any_of(facet_cols))) |>
      arrange(.data[[abscissa_col]]) |>
      mutate(.width = 0.8 * (lead(.data[[abscissa_col]]) - .data[[abscissa_col]])) |>
      fill(.width, .direction = "down") |>
      mutate(.xmin = .data[[abscissa_col]] - .width / 2,
             .xmax = .data[[abscissa_col]] + .width / 2) |>
      ungroup()
    
    # geom_rect() with explicit bounds rather than geom_col(aes(width = )):
    # `width` as an aesthetic is not reliably picked up by every ggplot2
    # version, geom_rect() always renders variable-width bars correctly.
    p <- ggplot(data, aes(xmin = .xmin, xmax = .xmax, ymin = 0, ymax = .data[[pdf_col]])) +
      geom_rect(fill = "grey40", alpha = 0.5)
    
    if (!is.null(theoretical) && nrow(theoretical) > 0) {
      p <- p + geom_line(data = theoretical, aes(x = abscissa, y = pdf),
                         colour = "black", linewidth = 0.5, inherit.aes = FALSE)
    }
    
    p <- add_facets(p, facet_cols)
    if (!is.null(x_trans)) p <- p + scale_x_continuous(transform = x_trans)
    if (!is.null(y_trans)) p <- p + scale_y_continuous(transform = y_trans)
    
    default_title <- title_with_filter(
      sprintf("%s distribution \u2014 %s", toupper(residual_type), model_name), filter_quo)
    return(
      p + labs(x = if (is.null(xlab)) toupper(residual_type) else xlab,
               y = if (is.null(ylab)) "Density" else ylab,
               title = if (is.null(title)) default_title else title,
               caption = "Bars: empirical density. Line: theoretical N(0,1) density.") +
        theme_bw(base_size = 11) +
        theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 12),
              plot.caption = element_text(size = 8, color = "grey40"))
    )
  }
  
  # ---- Time / Predictions / QQ: all read from *_residuals.txt --------------
  data <- relabel_columns(read_residuals(model_dir, dv = dv), value_labels)
  
  pred_col   <- paste0("prediction_", residual_type)
  simblq_col <- paste0(residual_type, "_simBlq")
  if (!simblq_col %in% colnames(data)) {
    stop("residual_type = '", residual_type, "' has no '", simblq_col, "' column in the residuals file.")
  }
  
  residual_cols <- c("ID", "Admin", "time", "censored", "split", "DV", "color", "filter",
                     colnames(data)[str_starts(colnames(data), "prediction_")],
                     colnames(data)[str_ends(colnames(data), "_simBlq")],
                     intersect(colnames(data), c("pwRes", "iwRes_mean", "iwRes_mode", "npde")))
  candidate_facets <- setdiff(colnames(data), residual_cols)
  
  r <- resolve_facets(data, filter_quo, candidate_facets,
                      facet_by = facet_by, facet_max_levels = facet_max_levels)
  data       <- r$data
  facet_cols <- r$facet_cols
  
  data <- data |> dplyr::filter(!is.na(.data[[simblq_col]]))
  if (nrow(data) == 0) {
    stop("No non-missing '", simblq_col, "' values left for residual_type = '", residual_type, "'.")
  }
  
  if ("censored" %in% colnames(data)) {
    data$censored_label <- factor(data$censored, levels = c(0, 1), labels = c("Above", "Below"))
  }
  
  if (type == "qq") {
    
    data <- data |>
      group_by(across(any_of(facet_cols))) |>
      arrange(.data[[simblq_col]]) |>
      mutate(.theoretical_q = qnorm(ppoints(n()))) |>
      ungroup()
    
    p <- ggplot(data, aes(x = .theoretical_q, y = .data[[simblq_col]])) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      layer_censored_points(data, highlight_censored)
    
    p <- add_facets(p, facet_cols)
    
    default_title <- title_with_filter(
      sprintf("%s Q-Q plot \u2014 %s", toupper(residual_type), model_name), filter_quo)
    return(
      p + labs(x = if (is.null(xlab)) "Theoretical quantiles (N(0,1))" else xlab,
               y = if (is.null(ylab)) toupper(residual_type) else ylab,
               title = if (is.null(title)) default_title else title,
               fill = "Quantification limit") +
        theme_bw(base_size = 11) +
        theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 12),
              legend.position = "bottom")
    )
  }
  
  # type %in% c("time", "predictions") -----------------------------------
  x_col <- if (type == "time") "time" else pred_col
  
  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[simblq_col]])) +
    layer_censored_points(data, highlight_censored) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
    # method = "gam" rather than "loess": far more robust on sparse/
    # clustered PK x-ranges, no extra dependency (mgcv ships with base R).
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, colour = "black", linewidth = 0.5)
  
  p <- add_facets(p, facet_cols)
  if (!is.null(x_trans)) p <- p + scale_x_continuous(transform = x_trans)
  if (!is.null(y_trans)) p <- p + scale_y_continuous(transform = y_trans)
  
  default_xlab <- if (type == "time") "Time" else {
    if (residual_type %in% c("npde", "pwRes")) "Population prediction" else "Individual prediction"
  }
  default_title <- title_with_filter(
    sprintf("%s vs %s \u2014 %s", toupper(residual_type),
            if (type == "time") "time" else "predictions", model_name), filter_quo)
  
  p + labs(x = if (is.null(xlab)) default_xlab else xlab,
           y = if (is.null(ylab)) toupper(residual_type) else ylab,
           title = if (is.null(title)) default_title else title,
           fill = "Quantification limit") +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom")
}
