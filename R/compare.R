    ### ---- MODULE 1c: compare_summaries() ------------------------------------
#' Internal helper: order-preserving merge of two parameter NAME VECTORS
#'
#' Core primitive shared by .merge_param_order() (2 models, used by
#' compare_two_models/2()) and .merge_param_order_multi() (N models, used
#' by compare_models/2()).
#'
#' The previous approach tried to RE-DERIVE a parameter/covariate hierarchy
#' from the parameter strings themselves (isolating a "root" via
#' str_split_fixed("_", n = 2)), then fell back to alphabetical order
#' (arrange(..., parameter)) as a tie-break within each root. That fallback
#' is exactly what broke the order: alphabetically "beta_ka2_FBMI" < "ka2_pop"
#' (b < k), so a covariate row got sorted BEFORE the base parameter it
#' belongs to -- which, combined with clean = TRUE turning "beta_ka2_FBMI"
#' into "~FBMI", looks like the covariate is now attached to the wrong
#' parameter even though the join itself was correct.
#'
#' This version never re-derives structure from the names at all (so it
#' doesn't care whether a parameter name contains underscores, e.g.
#' "ka_IM"). It simply keeps the NATIVE row order Monolix already writes in
#' summary.txt -- which already interleaves each base parameter with its own
#' covariate-effect rows correctly. `params1` is treated as the backbone and
#' kept in its own order. Any name present only in `params2` is spliced in
#' right after the nearest preceding name the two vectors have in common (or
#' at the very start if none has been seen yet), preserving `params2`'s own
#' relative order for any run of consecutive new names.
#'
#' @noRd
.merge_order_vec <- function(params1, params2) {
  backbone   <- params1
  anchor     <- NA_character_
  insertions <- list()
  
  for (p in params2) {
    if (p %in% backbone) {
      anchor <- p
    } else {
      key <- if (is.na(anchor)) "__START__" else anchor
      insertions[[key]] <- c(insertions[[key]], p)
    }
  }
  
  ordered <- insertions[["__START__"]]
  for (p in backbone) {
    ordered <- c(ordered, p)
    if (!is.null(insertions[[p]])) ordered <- c(ordered, insertions[[p]])
  }
  ordered
}
#' Internal helper: order-preserving merge of TWO parameter tables
#'
#' Used by compare_two_models() and compare_two_models().
#'
#' @return tibble [type, parameter, monolix_order]
#'
#' @noRd
.merge_param_order <- function(tbl_first, tbl_second, block_levels) {
  map_dfr(block_levels, function(blk) {
    p1 <- tbl_first  |> filter(type == blk) |> pull(parameter)
    p2 <- tbl_second |> filter(type == blk) |> pull(parameter)
    if (length(p1) == 0 && length(p2) == 0) return(NULL)
    tibble(type = blk, parameter = .merge_order_vec(p1, p2))
  }) |>
    mutate(monolix_order = row_number())
}
#' Internal helper: order-preserving merge of N parameter tables
#'
#' Used by compare_models() and compare_models(). Folds
#' .merge_order_vec() across all models in sequence: the first model sets
#' the backbone, and each subsequent model's rows are spliced in in turn --
#' so a covariate present only in the 3rd of 5 models still lands right
#' after the base parameter it belongs to, instead of being pushed wherever
#' alphabetical order happens to put it. The accumulator returned by each
#' fold step already contains every name seen so far (original backbone
#' plus anything already spliced in), so later models are always merged
#' against the full history, not just the first model.
#'
#' @return tibble [type, parameter, monolix_order]
#'
#' @noRd
.merge_param_order_multi <- function(tbls, block_levels) {
  map_dfr(block_levels, function(blk) {
    param_lists <- map(tbls, ~ filter(.x, type == blk) |> pull(parameter))
    param_lists <- keep(param_lists, ~ length(.x) > 0)
    if (length(param_lists) == 0) return(NULL)
    ordered <- reduce(param_lists, .merge_order_vec)
    tibble(type = blk, parameter = ordered)
  }) |>
    mutate(monolix_order = row_number())
}
#' build_param_table() drops the columns that are empty for the model (a partial
#'
#' estimation without standard errors has no se_sa / rse_pct / ci_* columns).
#' Comparisons need the same columns everywhere, so the missing ones are added as NA.
#'
#' @noRd
.ensure_columns <- function(tbl, cols) {
  for (v in setdiff(cols, names(tbl))) tbl[[v]] <- NA_real_
  tbl
}





# Accepts a read_summary(..., full = TRUE) object, a summary.txt path or a model folder
#' @noRd
.as_summary <- function(x) {
  if (is.list(x)) return(x)
  if (dir.exists(x)) x <- file.path(x, "summary.txt")
  read_summary(x, full = TRUE)
}

.block_levels <- c("Fixed Effects", "Standard Deviation of the Random Effects",
                   "Error Model Parameters", "Information Criteria")


#' Compare the parameters of two models side by side
#'
#' Joins the parameter tables of two models, keeping the native Monolix row
#' order, adds the difference (`second - first`) and the relative difference
#' for each requested column, and builds a ready-to-print table.
#'
#' @param first,second The two models: the output of `read_summary(file, full = TRUE)`,
#'   a summary.txt path, or a model folder.
#' @param values Columns to compare, e.g. `"estimate"` or `c("estimate", "rse_pct")`.
#'   Models without standard errors (partial estimation) have `NA` for the
#'   standard-error columns.
#' @param clean If TRUE (default), covariate labels are simplified
#'   (`beta_Cl_WT` becomes `~WT`).
#' @param bold_above Optional named numeric vector: values of a column above the
#'   given threshold are printed in bold in the table, e.g. `c(rse_pct = 30)`
#'   highlights the relative standard errors above 30 percent for both models.
#'   Only the table is affected, not `data`.
#'
#' @return A list with `data` (tibble: `type`, `parameter`, then for each value the
#'   columns `<value>_1`, `<value>_2`, `delta_<value>` and `delta_rel_<value>`) and
#'   `kable` (the formatted table, a kableExtra object you can keep styling).
#' @export
compare_two_models <- function(first, second, values = "estimate", clean = TRUE, bold_above = NULL) {
  first  <- .as_summary(first)
  second <- .as_summary(second)

  # One column per requested value and model, suffixed _1 / _2 up front so a single
  # pair of joins is enough
  tbl_first <- build_param_table(first) |>
    .ensure_columns(values) |>
    dplyr::select(type, parameter, dplyr::all_of(values)) |>
    dplyr::rename_with(~ paste0(.x, "_1"), dplyr::all_of(values))
  tbl_second <- build_param_table(second) |>
    .ensure_columns(values) |>
    dplyr::select(type, parameter, dplyr::all_of(values)) |>
    dplyr::rename_with(~ paste0(.x, "_2"), dplyr::all_of(values))

  # Order-preserving merge of the two parameter lists, then data onto that layout
  master_order <- .merge_param_order(tbl_first, tbl_second, .block_levels)
  compared_tbl <- master_order |>
    dplyr::left_join(tbl_first,  by = c("type", "parameter")) |>
    dplyr::left_join(tbl_second, by = c("type", "parameter")) |>
    dplyr::arrange(monolix_order) |>
    dplyr::select(-monolix_order)

  # Differences, per requested value
  for (v in values) {
    val1 <- as.numeric(compared_tbl[[paste0(v, "_1")]])
    val2 <- as.numeric(compared_tbl[[paste0(v, "_2")]])
    compared_tbl[[paste0("delta_", v)]]     <- round(val2 - val1, 3)
    compared_tbl[[paste0("delta_rel_", v)]] <- scales::label_percent(accuracy = 0.1)((val2 - val1) / val1)
  }

  # Columns grouped by value: v_1, v_2, delta, delta_rel
  ordered_cols <- c("type", "parameter",
                    unlist(lapply(values, function(v)
                      c(paste0(v, c("_1", "_2")), paste0(c("delta_", "delta_rel_"), v)))))
  compared_tbl <- dplyr::select(compared_tbl, dplyr::all_of(ordered_cols))

  # Labels are cleaned LAST: matching is done on the raw parameter names
  if (clean) compared_tbl <- clean_param_labels(compared_tbl, overwrite = TRUE)

  # Table: one block of four columns per value, models marked 1 and 2 (names in the footnote)
  value_label <- function(v) switch(v, estimate = "Estimate", rse_pct = "RSE (%)", v)
  col_names <- c("Parameter", unlist(lapply(values, function(v)
    c(paste0(value_label(v), "\u00b9"), paste0(value_label(v), "\u00b2"), "\u0394", "Relative \u0394"))))
  shown <- compared_tbl |>
    dplyr::select(-type) |>
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                \(x) ifelse(is.na(x) | x == "NA", "\u2014", as.character(x))))
  kable_obj <- shown |>
    knitr::kable(booktabs = TRUE, col.names = col_names, caption = "Parameter changes") |>
    auto_pack_rows(compared_tbl, "type") |>
    kableExtra::footnote(number = c(first$model_name, second$model_name)) |>
    kableExtra::kable_styling(full_width = FALSE, font_size = 11, bootstrap_options = c("striped", "bordered", "hover"))

  # Optional highlighting, e.g. bold_above = c(rse_pct = 30)
  for (v in intersect(names(bold_above), values)) {
    for (suffix in c("_1", "_2")) {
      col <- paste0(v, suffix)
      flag <- suppressWarnings(as.numeric(compared_tbl[[col]])) > bold_above[[v]]
      kable_obj <- kableExtra::column_spec(kable_obj, which(names(compared_tbl) == col) - 1,
                                           bold = !is.na(flag) & flag)
    }
  }

  list(data = compared_tbl, kable = kable_obj)
}


#' Compare the parameters of all the models of a folder
#'
#' Joins the parameter tables of every model found in `dir` (flat or nested
#' layout, see [read_models()]), keeping the native Monolix row order, and
#' builds a ready-to-print kable table. Purely visual: no difference is
#' computed (use [compare_two_models()] for that).
#'
#' @param dir Folder holding one sub-folder per model.
#' @param values Columns to compare: one value (`"estimate"`, one column per
#'   model) or several (`c("estimate", "rse_pct")`, columns grouped by model).
#'   Models without standard errors (e.g. COSSAC runs) show `---` for the
#'   standard-error columns.
#' @param clean If TRUE (default), covariate labels are simplified
#'   (`beta_Cl_WT` becomes `~WT`).
#' @param bold_above Optional named numeric vector: values above the threshold are
#'   printed in bold in the table, e.g. `c(rse_pct = 30)`. With one value the name must
#'   be that value (`values = "rse_pct"`); with several, the columns of that value are
#'   highlighted. Only the table is affected, not `data`.
#'
#' @return A list with `data` (the joined tibble) and `kable` (the formatted table).
#' @export
compare_models <- function(dir, values = "estimate", clean = TRUE, bold_above = NULL) {
  summary_files <- find_summary_files(dir)
  single <- length(values) == 1

  # One table per model; columns are named after the model folder
  all_models_data <- purrr::imap(summary_files, function(summary_file, model_id) {
    df_raw <- build_param_table(read_summary(summary_file, full = TRUE)) |>
      .ensure_columns(values) |>
      dplyr::select(type, parameter, dplyr::all_of(values))
    for (v in values) {
      new_name <- if (single) model_id else paste0(model_id, "_", v)
      df_raw <- dplyr::rename(df_raw, !!new_name := !!rlang::sym(v))
    }
    df_raw
  }) |> purrr::keep(~ !is.null(.))

  # Order-preserving merge across all models, then data onto that layout
  master_order <- .merge_param_order_multi(all_models_data, .block_levels)
  merged_tbl   <- purrr::reduce(all_models_data, function(df1, df2)
    dplyr::full_join(df1, df2, by = c("type", "parameter")))
  final_comp_table <- master_order |>
    dplyr::left_join(merged_tbl, by = c("type", "parameter")) |>
    dplyr::arrange(monolix_order) |>
    dplyr::select(-monolix_order)
  if (clean) final_comp_table <- clean_param_labels(final_comp_table, overwrite = TRUE)

  # Formatted table
  model_names <- setdiff(names(final_comp_table), c("type", "parameter"))
  kable_obj_data <- dplyr::select(final_comp_table, -type)    # displayed columns (after Parameter), to locate the highlighted ones
  kable_obj <- final_comp_table |>
    dplyr::select(-type) |>
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                \(x) ifelse(is.na(x) | x == "NA", "\u2014", as.character(x))))
  if (single) {
    kable_obj <- kable_obj |>
      knitr::kable(booktabs = TRUE, escape = FALSE, col.names = c("Parameter", model_names),
                   caption = paste("Global Comparison of", values, "across models")) |>
      auto_pack_rows(final_comp_table, "type") |>
      kableExtra::kable_styling(full_width = FALSE, font_size = 10, bootstrap_options = c("striped", "hover"))
  } else {
    top_header <- c(" " = 1)
    for (m in names(all_models_data)) top_header[basename(m)] <- length(values)
    kable_obj <- kable_obj |>
      knitr::kable(booktabs = TRUE, escape = FALSE,
                   col.names = c("Parameter", rep(values, length(all_models_data))),
                   caption = "Multi-Metric Global Model Comparison") |>
      kableExtra::add_header_above(top_header, bold = TRUE) |>
      auto_pack_rows(final_comp_table, "type") |>
      kableExtra::kable_styling(full_width = FALSE, font_size = 9, bootstrap_options = c("striped", "hover"))
  }
  kable_obj <- kableExtra::column_spec(kable_obj, 1, bold = TRUE)

  # Optional highlighting, e.g. bold_above = c(rse_pct = 30)
  for (v in intersect(names(bold_above), values)) {
    cols <- if (single) model_names else model_names[endsWith(model_names, paste0("_", v))]
    for (col in cols) {
      flag <- suppressWarnings(as.numeric(final_comp_table[[col]])) > bold_above[[v]]
      kable_obj <- kableExtra::column_spec(kable_obj, which(names(kable_obj_data) == col),
                                           bold = !is.na(flag) & flag)
    }
  }

  list(data = final_comp_table, kable = kable_obj)
}


#' Compare the information criteria of two models
#'
#' OFV, AIC, BICc and BIC of two models side by side, with the difference
#' `second - first`. These criteria are only comparable if the log-likelihood
#' of both models was computed the same way (importance sampling or
#' linearization): a warning is raised otherwise.
#'
#' @param first,second The two models: the output of `read_summary()`, a
#'   summary.txt path, or a model folder.
#' @param labels Names of the two models used as column names (default: the
#'   model names found in the summaries).
#'
#' @return A tibble: `criterion`, one column per model and `delta`
#'   (`second - first`). The likelihood methods are kept in the attribute
#'   `ll_method`.
#' @export
compare_criteria <- function(first, second, labels = NULL) {
  first  <- .as_summary(first)
  second <- .as_summary(second)
  if (is.null(labels)) labels <- c(first$model_name, second$model_name)
  if (anyDuplicated(labels)) labels <- paste0(labels, c(" (1)", " (2)"))

  out <- dplyr::full_join(first$criteria, second$criteria, by = "criterion", suffix = c("_1", "_2"))
  out <- dplyr::rename(out, !!labels[1] := value_1, !!labels[2] := value_2)
  out$delta <- round(out[[labels[2]]] - out[[labels[1]]], 2)

  methods <- c(first$ll_method, second$ll_method)
  if (length(methods) == 2 && !anyNA(methods) && methods[1] != methods[2]) {
    warning("The log-likelihood was computed by different methods (", methods[1], " vs ", methods[2],
            "): the criteria are not comparable.", call. = FALSE)
  }
  attr(out, "ll_method") <- stats::setNames(c(first$ll_method, second$ll_method), labels)
  out
}


#' Compare the prediction-error metrics of two models
#'
#' Puts the metrics of [obspred_diagnostics()] (MPE, RMSE, MAE, F20, ...) of
#' two models in one table: one row per metric, prediction type (population,
#' individual) and observed variable, one column per model and the difference
#' `second - first`.
#'
#' @param first,second The two models, given as: a model folder (its ChartsData is
#'   read), the output of `read_summary()` (the folder of the summary is used), the
#'   output of [obspred_diagnostics()], or its `metrics` table.
#' @param labels Names of the two models used as column names.
#'
#' @return A tibble: the identifying columns (`Y_TYPE` when present,
#'   `Prediction`), `metric`, one column per model and `delta`.
#' @export
compare_obspred_metrics <- function(first, second, labels = c("Model 1", "Model 2")) {
  as_long <- function(x, label) {
    m <- .as_metrics(x)
    keys <- intersect(c("DV", "Y_TYPE", "Prediction"), names(m))
    metrics <- setdiff(names(m), c(keys, "ipred_type"))
    m |>
      dplyr::select(dplyr::all_of(c(keys, metrics))) |>
      dplyr::mutate(dplyr::across(dplyr::all_of(metrics), as.numeric)) |>
      tidyr::pivot_longer(dplyr::all_of(metrics), names_to = "metric", values_to = label)
  }
  a <- as_long(first, labels[1]); b <- as_long(second, labels[2])
  keys <- intersect(setdiff(names(a), c("metric", labels[1])), names(b))
  out <- dplyr::full_join(a, b, by = c(keys, "metric"))
  # keep the order in which the first model lists its rows and metrics
  for (k in c(keys, "metric")) {
    out[[k]] <- factor(out[[k]], levels = unique(c(a[[k]], b[[k]])))
  }
  out <- dplyr::arrange(out, dplyr::across(dplyr::all_of(c(keys, "metric"))))
  for (k in c(keys, "metric")) out[[k]] <- as.character(out[[k]])
  out$delta  <- round(out[[labels[2]]] - out[[labels[1]]], 4)
  out
}


# Metrics table of one model, whatever the way it is given: a metrics table, the output of
# obspred_diagnostics(), a read_summary() object, or a model folder (ChartsData must exist there)
#' @noRd
.as_metrics <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && is.data.frame(x$metrics)) return(x$metrics)
  dir <- if (is.character(x) && length(x) == 1) x
         else if (is.list(x) && !is.null(x$summary_path)) dirname(x$summary_path)
  if (is.null(dir)) {
    stop("Cannot get the prediction-error metrics from this object: give a model folder, ",
         "a read_summary() object, the output of obspred_diagnostics() or its `metrics` table.",
         call. = FALSE)
  }
  if (!dir.exists(dir)) stop("Model folder not found: ", dir, call. = FALSE)
  check_chartsdata(dir)
  obspred_diagnostics(dir, dv = resolve_dv(dir), split = TRUE)$metrics
}
