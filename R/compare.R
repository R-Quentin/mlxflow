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


# What the `delta` argument asks for: TRUE / FALSE, or "abs" and / or "rel"
#' @noRd
.resolve_delta <- function(delta) {
  none <- c(abs = FALSE, rel = FALSE)
  if (is.null(delta) || (is.logical(delta) && length(delta) == 1 && !is.na(delta) && !delta)) return(none)
  if (is.logical(delta) && length(delta) == 1 && !is.na(delta)) return(c(abs = TRUE, rel = TRUE))
  if (is.character(delta) && length(delta) > 0 && all(delta %in% c("abs", "rel"))) {
    return(c(abs = "abs" %in% delta, rel = "rel" %in% delta))
  }
  stop("`delta` must be TRUE, FALSE, \"abs\", \"rel\" or c(\"abs\", \"rel\").", call. = FALSE)
}

# Relative difference as a percentage string; NA when it is not defined (0 / 0, x / 0)
#' @noRd
.fmt_rel_change <- function(x) {
  out <- scales::label_percent(accuracy = 0.1)(x)
  out[!is.finite(x)] <- NA_character_
  out
}

# Column header of a compared value
#' @noRd
.value_label <- function(v) switch(v, estimate = "Estimate", rse_pct = "RSE (%)", shrinkage = "Shrinkage (%)", v)


#' Compare the parameters of two models side by side
#'
#' Joins the parameter tables of two models, keeping the native Monolix row
#' order, adds the difference (`second - first`) and the relative difference
#' for each requested column, and builds a ready-to-print table.
#'
#' @param first,second The two models: the output of `read_summary(file, full = TRUE)`,
#'   a summary.txt path, or a model folder.
#' @param values Columns to compare, e.g. `"estimate"` or `c("estimate", "rse_pct")`
#'   (also `"shrinkage"`, `"se_sa"`, ...). Models without standard errors (partial
#'   estimation) have `NA` for the standard-error columns.
#' @param clean If TRUE (default), covariate labels are simplified
#'   (`beta_Cl_WT` becomes `~WT`).
#' @param bold_above Optional named numeric vector: values of a column above the
#'   given threshold are printed in bold in the table, e.g. `c(rse_pct = 30)`
#'   highlights the relative standard errors above 30 percent for both models.
#'   Only the table is affected, not `data`.
#' @param delta Which differences the **table** shows: `c("abs", "rel")` (default, or
#'   `TRUE`) the absolute difference and the relative one, `"abs"` or `"rel"` only one
#'   of them, `FALSE` none (the values of the two models only). Only the table is
#'   affected: `data` always holds both.
#'
#' @return A list with `data` (tibble: `type`, `parameter`, then for each value the
#'   columns `<value>_1`, `<value>_2`, `delta_<value>` and `delta_rel_<value>`) and
#'   `kable` (the formatted table, a kableExtra object you can keep styling).
#' @export
compare_two_models <- function(first, second, values = "estimate", clean = TRUE, bold_above = NULL,
                               delta = c("abs", "rel")) {
  show_delta <- .resolve_delta(delta)
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
    compared_tbl[[paste0("delta_rel_", v)]] <- .fmt_rel_change((val2 - val1) / val1)
  }

  # Columns grouped by value: v_1, v_2, delta, delta_rel
  ordered_cols <- c("type", "parameter",
                    unlist(lapply(values, function(v)
                      c(paste0(v, c("_1", "_2")), paste0(c("delta_", "delta_rel_"), v)))))
  compared_tbl <- dplyr::select(compared_tbl, dplyr::all_of(ordered_cols))

  # Labels are cleaned LAST: matching is done on the raw parameter names
  if (clean) compared_tbl <- clean_param_labels(compared_tbl, overwrite = TRUE)

  # Table: one block per value, models marked 1 and 2 (names in the footnote), then the
  # differences that were asked for
  shown_cols <- "parameter"
  col_names  <- "Parameter"
  for (v in values) {
    shown_cols <- c(shown_cols, paste0(v, c("_1", "_2")))
    col_names  <- c(col_names, paste0(.value_label(v), c("\u00b9", "\u00b2")))
    if (show_delta[["abs"]]) { shown_cols <- c(shown_cols, paste0("delta_", v));     col_names <- c(col_names, "\u0394") }
    if (show_delta[["rel"]]) { shown_cols <- c(shown_cols, paste0("delta_rel_", v)); col_names <- c(col_names, "Relative \u0394") }
  }
  shown <- compared_tbl |>
    dplyr::select(dplyr::all_of(shown_cols)) |>
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
      kable_obj <- kableExtra::column_spec(kable_obj, which(shown_cols == col),
                                           bold = !is.na(flag) & flag)
    }
  }

  list(data = compared_tbl, kable = kable_obj)
}


#' Compare the parameters of all the models of a folder
#'
#' Joins the parameter tables of every model found in `dir` (flat or nested
#' layout, see [read_models()]), keeping the native Monolix row order, and
#' builds a ready-to-print kable table. By default purely visual; with `delta`
#' the difference of every model to a reference model (the first one) is added.
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
#' @param delta Differences to the reference model, next to each other model: `FALSE`
#'   (default) none, `TRUE` or `c("abs", "rel")` the absolute and the relative
#'   difference, `"abs"` or `"rel"` only one of them. The columns are added to both
#'   `data` (`delta_<column>`, `delta_rel_<column>`) and the table; the reference
#'   model has none.
#' @param reference The reference model of the differences: its position in the
#'   table (default 1, the first folder in alphabetical order) or its folder name.
#'
#' @return A list with `data` (the joined tibble) and `kable` (the formatted table).
#' @export
compare_models <- function(dir, values = "estimate", clean = TRUE, bold_above = NULL,
                           delta = FALSE, reference = 1) {
  show_delta    <- .resolve_delta(delta)
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
  model_ids <- names(all_models_data)

  ref <- if (is.character(reference)) match(reference, model_ids) else as.integer(reference)
  if (length(ref) != 1 || is.na(ref) || ref < 1 || ref > length(model_ids)) {
    stop("`reference` must be a position between 1 and ", length(model_ids),
         " or one of the model folders: ", paste(model_ids, collapse = ", "), call. = FALSE)
  }
  value_col <- function(m, v) if (single) m else paste0(m, "_", v)

  # Order-preserving merge across all models, then data onto that layout
  master_order <- .merge_param_order_multi(all_models_data, .block_levels)
  merged_tbl   <- purrr::reduce(all_models_data, function(df1, df2)
    dplyr::full_join(df1, df2, by = c("type", "parameter")))
  final_comp_table <- master_order |>
    dplyr::left_join(merged_tbl, by = c("type", "parameter")) |>
    dplyr::arrange(monolix_order) |>
    dplyr::select(-monolix_order)
  if (clean) final_comp_table <- clean_param_labels(final_comp_table, overwrite = TRUE)

  # Differences to the reference model, and the layout of the table: per model, per value,
  # the value then its differences (the reference model has none)
  spec <- list()
  for (m in model_ids) {
    for (v in values) {
      col <- value_col(m, v)
      spec[[length(spec) + 1]] <- list(col = col, group = m, value = v, is_value = TRUE,
                                       label = if (single && !any(show_delta)) m else .value_label(v))
      if (m == model_ids[ref] || !any(show_delta)) next
      ref_val <- suppressWarnings(as.numeric(final_comp_table[[value_col(model_ids[ref], v)]]))
      val     <- suppressWarnings(as.numeric(final_comp_table[[col]]))
      if (show_delta[["abs"]]) {
        final_comp_table[[paste0("delta_", col)]] <- round(val - ref_val, 3)
        spec[[length(spec) + 1]] <- list(col = paste0("delta_", col), group = m, value = v,
                                         is_value = FALSE, label = "\u0394")
      }
      if (show_delta[["rel"]]) {
        final_comp_table[[paste0("delta_rel_", col)]] <- .fmt_rel_change((val - ref_val) / ref_val)
        spec[[length(spec) + 1]] <- list(col = paste0("delta_rel_", col), group = m, value = v,
                                         is_value = FALSE, label = "Relative \u0394")
      }
    }
  }
  spec <- dplyr::bind_rows(lapply(spec, tibble::as_tibble))
  final_comp_table <- dplyr::select(final_comp_table, type, parameter, dplyr::all_of(spec$col))

  # Formatted table
  shown_cols <- c("parameter", spec$col)
  kable_obj <- final_comp_table |>
    dplyr::select(dplyr::all_of(shown_cols)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                \(x) ifelse(is.na(x) | x == "NA", "\u2014", as.character(x))))
  if (single && !any(show_delta)) {
    kable_obj <- kable_obj |>
      knitr::kable(booktabs = TRUE, escape = FALSE, col.names = c("Parameter", spec$label),
                   caption = paste("Global Comparison of", values, "across models")) |>
      auto_pack_rows(final_comp_table, "type") |>
      kableExtra::kable_styling(full_width = FALSE, font_size = 10, bootstrap_options = c("striped", "hover"))
  } else {
    top_header <- c(" " = 1)
    for (m in model_ids) top_header[m] <- sum(spec$group == m)
    kable_obj <- kable_obj |>
      knitr::kable(booktabs = TRUE, escape = FALSE, col.names = c("Parameter", spec$label),
                   caption = if (any(show_delta)) {
                     sprintf("Global comparison of the models (changes: model - %s)", model_ids[ref])
                   } else "Multi-Metric Global Model Comparison") |>
      kableExtra::add_header_above(top_header, bold = TRUE) |>
      auto_pack_rows(final_comp_table, "type") |>
      kableExtra::kable_styling(full_width = FALSE, font_size = 9, bootstrap_options = c("striped", "hover"))
  }
  kable_obj <- kableExtra::column_spec(kable_obj, 1, bold = TRUE)

  # Optional highlighting, e.g. bold_above = c(rse_pct = 30)
  for (v in intersect(names(bold_above), values)) {
    for (col in spec$col[spec$is_value & spec$value == v]) {
      flag <- suppressWarnings(as.numeric(final_comp_table[[col]])) > bold_above[[v]]
      kable_obj <- kableExtra::column_spec(kable_obj, which(shown_cols == col),
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
#' two models in one table, with the difference `second - first`.
#'
#' By default (`format = "wide"`) there is one block per observed variable, stratum
#' and prediction type, with three rows (the two models and `"delta"`) and one column
#' per metric, which reads like a published performance table. `format = "long"`
#' gives one row per metric instead, and one column per model.
#'
#' Observed variables and strata are matched by name. When the two models do not
#' name their observed variable the same way but each has only one, they are matched
#' anyway (a message tells so). When only one of the models has strata (the `split`
#' of the Monolix plots), both are compared overall, which needs the model folders (or
#' the output of [obspred_diagnostics()]) to recompute the metrics; a metrics table
#' cannot be pooled and an error explains what to give instead.
#'
#' @param first,second The two models, given as: a model folder (its ChartsData is
#'   read), the output of `read_summary()` (the folder of the summary is used), the
#'   output of [obspred_diagnostics()], or its `metrics` table.
#' @param labels Names of the two models used in the table (default `"Model 1"` and
#'   `"Model 2"`).
#' @param format `"wide"` (default): rows are models, columns are metrics. `"long"`:
#'   one row per metric, one column per model.
#' @param split If TRUE (default), keep the strata of the Monolix `split` when both
#'   models have them: the metrics are then given per stratum. If FALSE, compare the
#'   metrics of all the observations together.
#' @param y_type_labels,split_labels Optional named vectors to rename (and order) the
#'   observed variables and the strata of the models read from a folder, so that the
#'   names of the two models agree, e.g. `y_type_labels = c(y1 = "Assay A")`. See
#'   [read_obspred()].
#'
#' @return A tibble. Wide: the identifying columns (`DV`, `Y_TYPE`, `split` and
#'   `Prediction`, when present), `Model` (the two labels, then `"delta"`), then one column
#'   per metric (`N`, `MPE`, ...). Long: the identifying columns, `metric`, one column
#'   per model and `delta`. The order of the first model is kept.
#' @export
compare_obspred_metrics <- function(first, second, labels = c("Model 1", "Model 2"),
                                    format = c("wide", "long"), split = TRUE,
                                    y_type_labels = NULL, split_labels = NULL) {
  format <- match.arg(format)
  if (length(labels) != 2) stop("`labels` must hold two names, one per model.", call. = FALSE)
  if (identical(labels[1], labels[2])) labels <- paste0(labels, c(" (1)", " (2)"))

  src1 <- .metrics_source(first,  split, y_type_labels, split_labels)
  src2 <- .metrics_source(second, split, y_type_labels, split_labels)
  m1 <- src1$metrics; m2 <- src2$metrics

  # Only one of the two has strata: pooled metrics for both, if they can be recomputed
  if (split && ("split" %in% names(m1)) != ("split" %in% names(m2))) {
    if (!src1$recompute || !src2$recompute) {
      stop("Only one of the two models has strata (a `split` column), and a metrics table ",
           "cannot be pooled. Give model folders or obspred_diagnostics() outputs, or the ",
           "metrics computed the same way for both (split = TRUE or FALSE in ",
           "obspred_diagnostics()).", call. = FALSE)
    }
    message("Only one of the models has strata (split): the metrics of all the observations ",
            "are compared.")
    m1 <- .metrics_source(first,  split = FALSE, y_type_labels, split_labels)$metrics
    m2 <- .metrics_source(second, split = FALSE, y_type_labels, split_labels)$metrics
  }

  # One observed variable in each model, named differently: they are the same thing here
  for (k in intersect(c("DV", "Y_TYPE"), intersect(names(m1), names(m2)))) {
    if (dplyr::n_distinct(m1[[k]]) == 1 && dplyr::n_distinct(m2[[k]]) == 1 && m1[[k]][1] != m2[[k]][1]) {
      message("`", k, "` is '", m1[[k]][1], "' in the first model and '", m2[[k]][1],
              "' in the second: matched as the same observed variable.")
      m2[[k]] <- m1[[k]][1]
    }
  }

  id_cols <- c("DV", "Y_TYPE", "split", "Prediction")
  as_long <- function(m, label) {
    keys    <- intersect(id_cols, names(m))
    metrics <- setdiff(names(m), c(keys, "ipred_type"))
    m |>
      dplyr::select(dplyr::all_of(c(keys, metrics))) |>
      dplyr::mutate(dplyr::across(dplyr::all_of(metrics), as.numeric)) |>
      tidyr::pivot_longer(dplyr::all_of(metrics), names_to = "metric", values_to = label)
  }
  a <- as_long(m1, labels[1]); b <- as_long(m2, labels[2])
  keys <- intersect(setdiff(names(a), c("metric", labels[1])), setdiff(names(b), c("metric", labels[2])))
  out <- dplyr::full_join(a, b, by = c(keys, "metric"))

  # What exists in one model only (a stratum, an observed variable) shows NA: say it
  lonely <- out[out$metric == "N" & (is.na(out[[labels[1]]]) | is.na(out[[labels[2]]])), keys, drop = FALSE]
  if (nrow(lonely) > 0) {
    warning("Rows present in one model only (NA for the other): ",
            paste(apply(lonely, 1, paste, collapse = " / "), collapse = "; "),
            ". Rename the observed variables or strata with `y_type_labels` / `split_labels`.",
            call. = FALSE)
  }

  # keep the order in which the first model lists its rows and metrics
  for (k in c(keys, "metric")) {
    out[[k]] <- factor(out[[k]], levels = unique(c(a[[k]], b[[k]])))
  }
  out <- dplyr::arrange(out, dplyr::across(dplyr::all_of(c(keys, "metric"))))
  for (k in c(keys, "metric")) out[[k]] <- as.character(out[[k]])
  out$delta <- round(out[[labels[2]]] - out[[labels[1]]], 4)
  if (format == "long") return(out)

  out |>
    tidyr::pivot_longer(dplyr::all_of(c(labels, "delta")), names_to = "Model", values_to = "value") |>
    tidyr::pivot_wider(names_from = "metric", values_from = "value") |>
    dplyr::mutate(Model = factor(.data$Model, levels = c(labels, "delta"))) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(keys), \(x) factor(x, levels = unique(x))), .data$Model) |>
    dplyr::mutate(Model = as.character(.data$Model))
}


# The metrics table of one model, whatever the way it is given, and whether it can be computed
# again with another `split`: a model folder (ChartsData must exist there), a read_summary()
# object, or the output of obspred_diagnostics() (its data are kept) can; a metrics table cannot
#' @noRd
.metrics_source <- function(x, split = TRUE, y_type_labels = NULL, split_labels = NULL) {
  if (is.data.frame(x)) return(list(metrics = x, recompute = FALSE))
  if (is.list(x) && is.data.frame(x$data) && "Y_TYPE" %in% names(x$data)) {
    d <- x$data
    attr(d, "split") <- isTRUE(split) && "split" %in% names(d) && dplyr::n_distinct(d$split) > 1
    return(list(metrics = .obspred_metrics_table(d, split = attr(d, "split")), recompute = TRUE))
  }
  if (is.list(x) && is.data.frame(x$metrics)) return(list(metrics = x$metrics, recompute = FALSE))
  dir <- if (is.character(x) && length(x) == 1) x
         else if (is.list(x) && !is.null(x$summary_path)) dirname(x$summary_path)
  if (is.null(dir)) {
    stop("Cannot get the prediction-error metrics from this object: give a model folder, ",
         "a read_summary() object, the output of obspred_diagnostics() or its `metrics` table.",
         call. = FALSE)
  }
  if (!dir.exists(dir)) stop("Model folder not found: ", dir, call. = FALSE)
  check_chartsdata(dir)
  res <- obspred_diagnostics(dir, dv = resolve_dv(dir), split = split,
                             y_type_labels = y_type_labels, split_labels = split_labels)
  list(metrics = res$metrics, recompute = TRUE)
}
