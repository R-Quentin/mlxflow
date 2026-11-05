# Locates the summary.txt of every model sub-folder of `dir`, whatever the layout:
#   flat   : <dir>/<model>/summary.txt            (default Monolix output, backward runs)
#   nested : <dir>/<model>/<model>/summary.txt    (Monolix ModelBuilding / SCM / COSSAC)
# Returns a named character vector: names = model folder names, values = summary paths.
# Stops with an explicit message if nothing is found (instead of failing later in reduce()).
#' @noRd
find_summary_files <- function(dir, layout = c("auto", "flat", "nested")) {
  layout <- match.arg(layout)
  if (!dir.exists(dir)) stop("Directory not found: ", dir, call. = FALSE)
  subs   <- as.character(fs::dir_ls(dir, type = "directory"))
  flat   <- file.path(subs, "summary.txt")
  nested <- file.path(subs, basename(subs), "summary.txt")
  found  <- switch(layout,
    auto   = ifelse(file.exists(flat), flat, ifelse(file.exists(nested), nested, NA_character_)),
    flat   = ifelse(file.exists(flat), flat, NA_character_),
    nested = ifelse(file.exists(nested), nested, NA_character_))
  names(found) <- basename(subs)
  found <- found[!is.na(found)]
  if (length(found) == 0) {
    stop("No summary.txt found in the sub-folders of ", dir,
         " (looked for <model>/summary.txt and <model>/<model>/summary.txt)", call. = FALSE)
  }
  found
}


#' Read every model of a folder into one table
#'
#' Reads the summary.txt of each model sub-folder and returns one row per model
#' (information criteria, covariate relations, ...). The folder layout is
#' detected for each sub-folder: either `<path>/<model>/summary.txt` (default
#' Monolix output, backward runs) or `<path>/<model>/<model>/summary.txt`
#' (Monolix model building runs: SCM, COSSAC).
#'
#' `model_id` is the number carried by the folder name (`ModelBuilding_016`
#' gives 16, `run003_I1_V_WT` gives 3). If the names carry no unique number, ids
#' follow the alphabetical order of the folder names, with a warning; this is
#' harmless for a backward flowchart but not for SCM/COSSAC flowcharts, whose
#' builders read `model_id` as the order in which the models were tested.
#' `model_name` is the folder name.
#'
#' @param path Folder holding one sub-folder per model.
#' @param full If TRUE, also parse the population/individual parameters and the
#'   Fisher matrix of each summary (slower).
#' @param layout `"auto"` (default) detects the layout of each sub-folder;
#'   `"flat"` or `"nested"` force one of the two.
#'
#' @return A tibble with one row per model, ordered by `model_id`. `ll_method` gives how
#'   the log-likelihood was computed (`"is"` importance sampling, `"lin"` linearization);
#'   a warning is raised if the models do not share the same method.
#' @export
read_models <- function(path, full = FALSE, layout = c("auto", "flat", "nested")) {
  layout <- match.arg(layout)
  files  <- find_summary_files(path, layout = layout)
  ids    <- stats::setNames(.assign_model_ids(names(files)), names(files))
  message(sprintf("  %d models found in %s", length(files), path))

  models <- purrr::imap(files, function(summary_path, run_name) {
    res <- .parse_run(summary_path, run_name, full = full)
    if (is.null(res)) return(NULL)
    dplyr::mutate(res, model_id = ids[[run_name]], model_name = run_name)
  }) |>
    purrr::list_rbind() |>
    dplyr::arrange(model_id) |>
    .add_relation_counts()

  # OFV, AIC and BIC from different likelihood methods must not be compared
  methods <- unique(stats::na.omit(models$ll_method))
  if (length(methods) > 1) {
    warning("The log-likelihood of these models was not computed by the same method (",
            paste(methods, collapse = ", "), "): their OFV, AIC and BIC are not comparable.", call. = FALSE)
  }
  models
}



  ## Other imported function ---------------------------------------------------
#' Strip categorical suffixes from covariate labels
#'
#' Turns the labels found in summary.txt into one label per covariate relation, e.g. "Cl_tTabac_G_Current_smoker" becomes "Cl_tTabac". Labels without suffix are returned unchanged.
#'
#' @param x Character vector of covariate labels ("<parameter>_<covariate>" with an optional "_G_<category>" suffix).
#'
#' @return character vector of the same length as `x`.
#'
#' @export
clean_cov_labels <- function(x) {
  x |>
    # Remove "_G_<category>" suffix (Monolix categorical encoding)
    str_remove_all("_G_[A-Za-z0-9_ ]+") |>
    # Collapse any double underscores or trailing underscores left behind
    str_replace_all("_+", "_") |>
    str_remove("_$")
}
#' Reads a summary.txt and returns a one-row tibble
#'
#' @noRd
.parse_run <- function(summary_path, run_name, full = FALSE) {
  
  parsed <- tryCatch(
    read_summary(summary_path, full = full),
    error = function(e) {
      warning(glue("Error reading {run_name}: {e$message}"))
      NULL
    }
  )
  
  if (is.null(parsed)) return(NULL)
  
  model_nm <- if (!is.null(parsed$model_name) && !is.na(parsed$model_name)) {
    parsed$model_name
  } else {
    run_name
  }
  
  pull_scalar <- function(df, crit) {
    val <- df |> filter(criterion == crit) |> pull(value)
    if (length(val) == 0) NA_real_ else val[[1]]
  }
  
  covs <- if (is.null(parsed$covariates)) character(0) else parsed$covariates
  
  tibble(
    model_id   = as.integer(str_extract(run_name, "\\d+")),
    model_name = model_nm,
    OFV        = pull_scalar(parsed$criteria, "OFV"),
    AIC        = pull_scalar(parsed$criteria, "AIC"),
    BICc       = pull_scalar(parsed$criteria, "BICc"),
    BIC        = pull_scalar(parsed$criteria, "BIC"),
    covariates = list(covs),
    n_cov      = length(covs),
    ll_method  = if (is.null(parsed$ll_method)) NA_character_ else parsed$ll_method
  )
}






  ## MODULE 2: read_all_models_()   ============================================
#' NB: .parse_run() (reads one summary.txt into a one-row tibble) lives in
#'
#' model_diagnostic_functions.R, next to read_summary().
#' Assigns a unique integer model_id to each run folder.
#' Folders such as "run003_I1_V_WT", "Backward_003" or "ModelBuilding_016" keep
#' the number they carry (first number in the name). When that number is
#' missing or not unique, ids fall back to the natural alphabetical order of
#' the folder names, from 0, with a warning: this is harmless for the backward
#' flowchart (iterations come from the covariate sets) but NOT for SCM/COSSAC,
#' whose builders read model_id as the chronological order of the runs.
#' model_id is used as the graph vertex name, so it MUST be unique.
#'
#' @noRd
.assign_model_ids <- function(run_names) {
  ids <- as.integer(str_extract(run_names, "\\d+"))
  if (anyNA(ids) || anyDuplicated(ids) > 0) {
    ord <- str_order(run_names, numeric = TRUE)
    ids <- integer(length(run_names))
    ids[ord] <- seq_along(run_names) - 1L
    warning("Folder names do not carry unique numbers: model_id assigned from the ",
            "alphabetical order of the folder names. Do not use this for SCM/COSSAC ",
            "results (their builders assume model_id follows the testing order).",
            call. = FALSE)
  }
  ids
}
#' NB: model_name is taken from the FOLDER name, not from the project name written inside
#'
#' summary.txt: a run is estimated while the parent project is still loaded and is only
#' saved under its own name afterwards, so summary.txt can carry the parent's name.
#' Relation count and iteration index. n_cov counts beta coefficients (one per
#' non-reference level of a categorical covariate), n_rel counts relations, so
#' that a 3-level categorical covariate is ONE relation of 2 coefficients.
#'
#' @noRd
.add_relation_counts <- function(df) {
  df |>
    mutate(n_rel = map_int(covariates, ~ length(unique(clean_cov_labels(.x)))),
           iter  = n_rel)
}

    ### ---- MODULE 2a: read_models() ---



    ### ---- MODULE 2b: read_models() ---



    #### ---- HELPER: clean_cov_labels() ---
# NB: clean_cov_labels() lives in model_diagnostic_functions.R.
