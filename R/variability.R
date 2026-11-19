  ## Testing the addition of variability (IOV or IIV) on the parameters of a model ---------

# Files and folders created by test_variability() in output_dir (same run<NNN>_ naming as
# backward_elimination(), plus its own log). Only these names are ever touched.
#' @noRd
list_variability_outputs <- function(output_dir) {
  list.files(output_dir, full.names = TRUE,
             pattern = "^(run[0-9]{3}_.*|variability_log\\.csv|variability_state\\.rds)$")
}
#' One row of variability_log.csv (run 0 = base model, then one row per tested model)
#' @noRd
variability_log_row <- function(run, iter, parameter, level, fit, elapsed_sec, project,
                                alpha = NA_real_, delta_ofv = NA_real_, df = NA_real_,
                                critical = NA_real_, p_value = NA_real_, significant = NA,
                                rse_new = NA_real_, eligible = NA, status = "ok") {
  data.frame(
    run = run, iter = iter, parameter = parameter, level = level,
    delta_ofv = delta_ofv, df = df, critical_ofv = critical, alpha = alpha,
    p_value = p_value, significant = significant, rse_new = rse_new, eligible = eligible,
    ofv = if (is.null(fit)) NA_real_ else crit_value(fit, "OFV"),
    bicc = if (is.null(fit)) NA_real_ else crit_value(fit, "BICc"),
    elapsed_sec = elapsed_sec, status = status, project = basename(project),
    stringsAsFactors = FALSE
  )
}
#' @noRd
log_variability <- function(output_dir, row) {
  log_file <- file.path(output_dir, "variability_log.csv")
  write.table(round_cols(row), log_file, sep = ";", row.names = FALSE,
              col.names = !file.exists(log_file), append = file.exists(log_file))
}
#' Switch the variability of one parameter on, at one level, in the loaded project
#' (only edits the in-memory model, nothing is run)
#' @noRd
set_variability <- function(level, parameter) {
  args <- stats::setNames(list(stats::setNames(list(TRUE), parameter)), level)
  do.call(.lx, c(list("setIndividualParameterVariability"), args))
  invisible(TRUE)
}
#' Save the project, then compute the charts data. In this order the charts are written to the
#' saved project (checked on Monolix: saveProject() alone writes no charts, saveProject() then
#' computeChartsData() does; a second save is not needed).
#' @noRd
save_with_charts <- function(path, export_charts = TRUE) {
  .lx("saveProject", path)
  if (export_charts) .lx("computeChartsData")
  invisible(path)
}
#' RSE (percent) of the parameter(s) a candidate added, from the standard errors of the loaded
#' project; NA when they were not estimated
#' @noRd
new_parameter_rse <- function(parameters, linearization = FALSE) {
  if (length(parameters) == 0) return(NA_real_)
  tryCatch({
    se  <- .lx("getEstimatedStandardErrors")
    tab <- if (linearization) se$linearization else se$stochasticApproximation
    rse <- suppressWarnings(as.numeric(tab$rse[tab$parameter %in% parameters]))
    if (length(rse) == 0 || all(is.na(rse))) NA_real_ else max(rse, na.rm = TRUE)
  }, error = function(e) NA_real_)
}
#' Estimate one model of the search in its own folder, save it, and compare it with the reference
#' by a likelihood ratio test. The degrees of freedom are the number of population parameters
#' that the change added (1 for the variance of one random effect).
#' @noRd
test_variability_one <- function(parent_project, level, parameter, run, iter, output_dir,
                                 reference_fit, alpha, linearization, scenario, warm_start,
                                 export_charts, max_rse = NULL) {
  t0 <- Sys.time()
  name <- run_name(run, iter, parameter, level)
  path <- file.path(output_dir, paste0(name, ".mlxtran"))

  load_project_quiet(parent_project)
  # Results go straight to this run's folder: the parent's results are never overwritten
  .lx("setProjectSettings", directory = file.path(output_dir, name))
  # Start from what the parent estimated (before the change, while results and model agree)
  if (warm_start) {
    tryCatch(.lx("setInitialEstimatesToLastEstimates"),
             error = function(e) warning("Could not start from the last estimates (", conditionMessage(e),
                                         "): the initial values of the project are used.", call. = FALSE))
  }
  set_variability(level, parameter)

  fit <- run_full_estimation(linearization = linearization, scenario = scenario)
  save_with_charts(path, export_charts)

  df        <- fit$n_params - reference_fit$n_params
  ofv       <- crit_value(fit, "OFV")
  delta_ofv <- ofv - crit_value(reference_fit, "OFV")

  # A likelihood that is not a finite number (importance sampling that breaks down, a diverged
  # estimation) cannot be tested: the model is kept in the table, flagged, and never adopted
  valid <- is.finite(delta_ofv)
  # Adding parameters can only lower the OFV: a rise means the estimation did not do as well
  p_value <- if (valid) 1 - stats::pchisq(max(-delta_ofv, 0), df = max(df, 1)) else NA_real_
  significant <- valid && p_value < alpha

  # Quality of the new random effect: its relative standard error
  new_par <- setdiff(names(fit$pop_params), names(reference_fit$pop_params))
  rse_new <- new_parameter_rse(new_par, linearization)
  precise <- is.null(max_rse) || (is.finite(rse_new) && rse_new <= max_rse)

  list(parameter = parameter, level = level, run = run, df = df, delta_ofv = if (valid) delta_ofv else NA_real_,
       p_value = p_value, significant = significant, eligible = significant && precise,
       rse_new = rse_new, new_parameter = paste(new_par, collapse = ", "), ofv_raw = ofv,
       status = if (valid) "ok" else "invalid: OFV is not a finite number",
       fit = fit, path = path, elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}


#' One row per estimated model, in run order: the base model, then every tested model
#' (failed ones included), with the verdict of the search
#' @noRd
variability_runs_table <- function(base_row, cands, adopted) {
  rows <- lapply(cands, function(x) {
    failed <- isTRUE(x$failed)
    data.frame(run = x$run, iter = x$iter, parameter = x$parameter, level = x$level,
               ofv = if (failed) NA_real_ else crit_value(x$fit, "OFV"),
               bicc = if (failed) NA_real_ else crit_value(x$fit, "BICc"),
               delta_ofv = if (failed) NA_real_ else x$delta_ofv,
               df = if (failed) NA_real_ else x$df,
               p_value = if (failed) NA_real_ else x$p_value,
               significant = if (failed) NA else x$significant,
               rse_new = if (failed) NA_real_ else x$rse_new,
               eligible = if (failed) NA else x$eligible,
               adopted = x$run %in% adopted,
               status = if (failed) paste("error:", x$status) else x$status,
               elapsed_sec = if (failed) NA_real_ else x$elapsed,
               project = basename(x$path), stringsAsFactors = FALSE)
  })
  tibble::as_tibble(do.call(rbind, c(list(base_row), rows)))
}

#' Test the addition of variability (IOV or IIV) on the parameters of a model
#'
#' Adds a random effect to one parameter at a time, at one variability level, and compares
#' every model with the reference by a likelihood ratio test (LRT). With the level of an
#' occasion (e.g. `"OCC"`) it tests **inter-occasion variability (IOV)**; with the level of the
#' individuals (`"id"`) it tests **inter-individual variability (IIV)**, i.e. whether a
#' parameter, fixed so far, should have an `omega`.
#'
#' Two ways of searching:
#' \itemize{
#'   \item `forward = FALSE` (default): one pass. Run 0 is the base model, then one run per
#'     tested parameter, each differing from the base model by that single random effect.
#'   \item `forward = TRUE`: stepwise forward selection, the counterpart of
#'     [backward_elimination()]. At each iteration every remaining parameter is tested on the
#'     current model, the significant one with the lowest p-value is adopted, and the search
#'     goes on from there until nothing is significant (or `max_iter`).
#' }
#'
#' The base model is estimated again first (run 0), so that every OFV comes from the same
#' estimation chain. Like in [backward_elimination()], everything is written in `output_dir`
#' (no sub-folder per model to browse): each model is
#' `run<NNN>_I<k>_<parameter>_<level>.mlxtran` next to its result folder of the same name
#' (`NNN` is never reused, `k` is the iteration), the starting project and its results are left
#' untouched, and `variability_log.csv` holds every test. The folders can be read with
#' [read_models()] and compared with [compare_models()] and [compare_obspred_metrics()].
#'
#' An OFV that is not a finite number (`-Inf`, `NaN`), or that changes by thousands for one added
#' variance, usually means that the importance-sampling likelihood broke down, not that the model got
#' better. Such a model is flagged `invalid` and never adopted; look at the standard error of the
#' likelihood in its `summary.txt`, compare with `linearization = TRUE`, and use `max_rse` to
#' require a well estimated new variance.
#'
#' The LRT of a variance is approximate: the null value (0) is on the boundary of the parameter
#' space, so the test is conservative, and an importance-sampling likelihood is itself
#' noisy. Read the p-values together with the shrinkage and the diagnostics of the models.
#'
#' @param project_file Path of the model (.mlxtran) to start from.
#' @param output_dir Folder of the runs. It must not hold the outputs of a previous run
#'   (see `resume` and `overwrite`).
#' @param parameters Parameters to test (names as in the individual model, e.g. `c("ka", "V")`).
#'   Default: all the parameters of the model. Those that already have a random effect at
#'   `level` are skipped, with a message.
#' @param level Variability level of the new random effect, as named in the model (see
#'   `getIndividualParameterModel()$variability`). Default: the first level that is not the
#'   individual level (`"id"`), which is the usual way to test IOV. A model with that level
#'   only gives an error: to test IIV, set `level = "id"`.
#' @param alpha Level of the LRT to keep a random effect (default 0.05). Only decides what
#'   is adopted when `forward = TRUE`; the p-values are given in any case.
#' @param forward If TRUE, select the random effects stepwise (see above); if FALSE
#'   (default) test each parameter once against the base model.
#' @param max_iter Maximum number of iterations of the forward selection (default 10).
#' @param prune_alpha Only with `forward = TRUE`. After the first iteration (every parameter
#'   against the base model), the parameters whose p-value is not below `prune_alpha` are
#'   dropped from the next iterations instead of being tested again. This is the pre-selection:
#'   with 9 parameters and 3 candidates left it takes 9 + 3 + 2 + 1 estimations instead of 9 + 8 + 7 + ...
#'   Use a value above `alpha` (e.g. 0.1) to be cautious. The risk is to miss an effect that only
#'   appears once another one is in the model. Default `NULL`: every remaining parameter is tested
#'   at every iteration.
#' @param max_rse Maximum relative standard error (percent) of the new random effect (its
#'   `gamma_` or `omega_` parameter) for a model to be adopted. A model whose LRT is significant
#'   but whose new variance is poorly estimated stays in the table (`eligible = FALSE`) and is
#'   not adopted. It needs the standard errors: use a `scenario` that enables
#'   `standardErrorEstimation`. Default `NULL`: no condition on the precision.
#' @param scenario Optional Monolix scenario (a list as returned by `getScenario()`) used to
#'   estimate every model, for instance a fast one without standard errors or a complete one
#'   with all the diagnostics. Default `NULL`: SAEM, conditional distribution sampling and
#'   log-likelihood, which is what the LRT needs. See [backward_elimination()] for what the
#'   scenario must enable.
#' @param linearization If TRUE, the log-likelihood is computed by linearization instead of
#'   importance sampling. Ignored when `scenario` is given (its `linearization` is used).
#' @param warm_start If TRUE (default), each tested model starts its estimation from the last
#'   estimates of the model it derives from (the base model, or the last adopted one) instead of
#'   the initial values of the project. This is not run 0: the base model is always estimated
#'   once at the start, whatever this argument says.
#' @param export_charts If TRUE (default), the charts data (VPC, observations vs predictions,
#'   ...) of every model are computed and saved, which the diagnostic functions need.
#'   FALSE saves time.
#' @param resume If TRUE and `output_dir` holds an interrupted run, continue it from the last
#'   estimated model (what an interrupted estimation left behind is discarded). The settings
#'   of the run (level, parameters, alpha, forward) are those of the interrupted run. If that
#'   run already finished, its stored result is returned without recomputing.
#' @param overwrite If TRUE, delete the outputs of a previous run in `output_dir` (only the
#'   files and folders created by this function) and start over.
#'
#' @return A list, invisibly: `runs` (one row per estimated model: iteration, parameter, level,
#'   OFV, BICc, `delta_ofv` and `p_value` against the reference, `significant`, `rse_new` (RSE of
#'   the new variance), `eligible` (significant and precise enough to be adopted), `adopted`,
#'   `status`; a model whose OFV is not a finite number has `status` "invalid" and no test), `history` (the adopted models, base model first), `final_project` (the model
#'   the search ends on: the base model when nothing is adopted), `stop_reason` and `settings`
#'   (level, alpha, ... of the run). Draw the search with [build_variability_data()] and [plot_flowchart()].
#' @examples
#' \dontrun{
#' init_monolix("iov_project.mlxtran")
#'
#' # Is there inter-occasion variability, and on which parameter?
#' res <- test_variability("iov_project.mlxtran", output_dir = "iov_project/IOV_test", level = "OCC")
#' res$runs
#'
#' # Stepwise selection, with a quick scenario; safe to launch again after an interruption
#' res <- test_variability("iov_project.mlxtran", output_dir = "iov_project/IOV_forward",
#'                         level = "OCC", forward = TRUE, alpha = 0.05, resume = TRUE,
#'                         scenario = list(tasks = c(populationParameterEstimation = TRUE,
#'                                                   conditionalDistributionSampling = TRUE,
#'                                                   logLikelihoodEstimation = TRUE),
#'                                         linearization = FALSE))
#'
#' # Does a parameter need a random effect at all (IIV)?
#' test_variability("model.mlxtran", output_dir = "model/IIV_test", level = "id")
#' }
#' @export
test_variability <- function(project_file,
                             output_dir,
                             parameters    = NULL,
                             level         = NULL,
                             alpha         = 0.05,
                             forward       = FALSE,
                             max_iter      = 10,
                             prune_alpha   = NULL,
                             max_rse       = NULL,
                             scenario      = NULL,
                             linearization = FALSE,
                             warm_start    = TRUE,
                             export_charts = TRUE,
                             resume        = FALSE,
                             overwrite     = FALSE) {

  total_t0 <- Sys.time()
  mlx_console$seen <- character()
  if (!is.null(scenario)) {
    check_scenario(scenario)
    linearization <- isTRUE(scenario$linearization)
  }

  if (!is.null(max_rse) && !(!is.null(scenario) && scenario_task(scenario, "standardErrorEstimation"))) {
    stop("`max_rse` needs the standard errors of every model: give a `scenario` that enables ",
         "standardErrorEstimation.", call. = FALSE)
  }
  if (!is.null(prune_alpha) && !forward) stop("`prune_alpha` is for the forward selection (forward = TRUE).", call. = FALSE)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_dir   <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  project_file <- normalizePath(project_file, winslash = "/", mustWork = TRUE)
  state_file   <- file.path(output_dir, "variability_state.rds")

  # Never silently mix a new run with the outputs of a previous one
  previous <- list_variability_outputs(output_dir)
  if (overwrite && length(previous) > 0) {
    unlink(previous, recursive = TRUE)
    message(sprintf("Removed %d file(s)/folder(s) from a previous run.", length(previous)))
    previous <- character()
  }
  if (length(previous) > 0) {
    if (!resume) {
      stop("output_dir already contains the outputs of a previous run. Use resume = TRUE to ",
           "continue it, or overwrite = TRUE to delete it and start over.", call. = FALSE)
    }
    if (!file.exists(state_file)) {
      stop("output_dir contains outputs but no checkpoint (variability_state.rds), so the run ",
           "cannot be resumed safely. Use overwrite = TRUE to start over.", call. = FALSE)
    }
  }

  if (length(previous) > 0) {
    ## Resume ------------------------------------------------------------------
    state <- readRDS(state_file)
    if (isTRUE(state$finished)) {
      message(sprintf("This run is already complete (%s). Returning the stored result; ",
                      state$result$stop_reason), "use overwrite = TRUE to run it again.")
      return(invisible(state$result))
    }
    settings <- state$settings
    if ((!is.null(level) && !identical(level, settings$level)) ||
        (!is.null(parameters) && !setequal(parameters, settings$parameters)) ||
        (!missing(alpha) && !identical(alpha, settings$alpha)) ||
        (!missing(forward) && !identical(forward, settings$forward)) ||
        (!is.null(prune_alpha) && !identical(prune_alpha, settings$prune_alpha)) ||
        (!is.null(max_rse) && !identical(max_rse, settings$max_rse))) {
      stop("The settings differ from those of the interrupted run (level '", settings$level,
           "', parameters ", paste(settings$parameters, collapse = ", "), ", alpha ", settings$alpha,
           ", forward = ", settings$forward, "). Use the same ones, or overwrite = TRUE to start over.",
           call. = FALSE)
    }
    level <- settings$level; parameters <- settings$parameters
    alpha <- settings$alpha; forward <- settings$forward
    prune_alpha <- settings$prune_alpha; max_rse <- settings$max_rse
    current_project <- state$current_project; reference_fit <- state$reference_fit
    cands <- state$cands; adopted <- state$adopted; history <- state$history
    remaining <- state$remaining; iter_done <- state$iter_done; run <- state$last_run
    clean_incomplete_runs(output_dir, run, log_name = "variability_log.csv")
    load_project_quiet(current_project)
    message(sprintf("Resuming after run %d (iteration %d completed, model: %s).",
                    run, iter_done, basename(current_project)))
    base_row <- state$base_row
  } else {
    ## Start -------------------------------------------------------------------
    # Defaults and checks, from the model itself
    load_project_quiet(project_file)
    model  <- .lx("getIndividualParameterModel")
    levels <- names(model$variability)

    if (is.null(level)) {
      other <- levels[tolower(levels) != "id"]
      if (length(other) == 0) {
        stop("Only one variability level found: '", levels[1], "'. If you are trying to test IIV ",
             "(omega), use level = '", levels[1], "'.", call. = FALSE)
      }
      level <- other[1]
      message("Variability level: '", level, "' (the first level that is not the individual one).")
    }
    if (!level %in% levels) {
      stop("Unknown variability level: ", level, ". Available: ", paste(levels, collapse = ", "), call. = FALSE)
    }
    if (is.null(parameters)) parameters <- model$name
    unknown <- setdiff(parameters, model$name)
    if (length(unknown) > 0) {
      stop("Unknown parameter(s): ", paste(unknown, collapse = ", "),
           ". Available: ", paste(model$name, collapse = ", "), call. = FALSE)
    }
    already <- names(which(unlist(model$variability[[level]])))
    skipped <- intersect(parameters, already)
    if (length(skipped) > 0) {
      message("Skipping (already a random effect on '", level, "'): ", paste(skipped, collapse = ", "))
    }
    parameters <- setdiff(parameters, skipped)
    if (length(parameters) == 0) stop("Nothing to test: every parameter already has a random effect on '", level, "'.", call. = FALSE)
    settings <- list(level = level, parameters = parameters, alpha = alpha, forward = forward,
                     prune_alpha = prune_alpha, max_rse = max_rse)

    # Run 0: the base model, estimated like every candidate
    run       <- 0L
    base_name <- run_name(run, 0)
    current_project <- file.path(output_dir, paste0(base_name, ".mlxtran"))
    message("Run 0 (base model)")
    .lx("setProjectSettings", directory = file.path(output_dir, base_name))
    t0 <- Sys.time()
    reference_fit <- run_full_estimation(linearization = linearization, scenario = scenario)
    base_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    save_with_charts(current_project, export_charts)
    if (!is.finite(crit_value(reference_fit, "OFV"))) {
      stop("The OFV of the base model is not a finite number (", crit_value(reference_fit, "OFV"),
           "): the likelihood cannot be used. Check the estimation of the base model, or try ",
           "linearization = TRUE.", call. = FALSE)
    }
    message(sprintf("  %s fitted in %s: OFV = %.2f, BICc = %.2f", base_name, fmt_duration(base_secs),
                    crit_value(reference_fit, "OFV"), crit_value(reference_fit, "BICc")))
    log_variability(output_dir, variability_log_row(run, 0, NA_character_, level, reference_fit,
                                                     base_secs, current_project, alpha = alpha))
    base_row <- data.frame(run = run, iter = 0L, parameter = NA_character_, level = level,
                           ofv = crit_value(reference_fit, "OFV"), bicc = crit_value(reference_fit, "BICc"),
                           delta_ofv = NA_real_, df = NA_real_, p_value = NA_real_, significant = NA,
                           rse_new = NA_real_, eligible = NA,
                           adopted = NA, status = "ok", elapsed_sec = base_secs,
                           project = basename(current_project), stringsAsFactors = FALSE)
    cands     <- list()      # every tested model, failed ones included
    adopted   <- integer()   # runs adopted by the forward selection
    history   <- list(base_row[, c("iter", "run", "parameter", "level", "ofv", "bicc", "project")])
    remaining <- parameters
    iter_done <- 0L
  }
  # What a resumed run needs is saved after every estimated model
  save_state <- function(finished = FALSE, result = NULL) {
    saveRDS(list(settings = settings, current_project = current_project, reference_fit = reference_fit,
                 cands = cands, adopted = adopted, history = history, remaining = remaining,
                 iter_done = iter_done, last_run = run, base_row = base_row,
                 finished = finished, result = result), state_file)
  }
  save_state(finished = FALSE, result = NULL)

  stop_reason <- if (forward) sprintf("max_iter (%d) reached", max_iter) else "single pass (forward = FALSE)"
  iters <- if (forward) { if (iter_done + 1 <= max_iter) (iter_done + 1):max_iter else integer(0) } else if (iter_done < 1) 1L else integer(0)
  if (!forward && iter_done >= 1) stop_reason <- "single pass (forward = FALSE)"

  for (iter in iters) {
    if (length(remaining) == 0) { stop_reason <- "no parameter left to test"; break }
    tested  <- vapply(Filter(function(x) x$iter == iter, cands), function(x) x$parameter, character(1))
    to_test <- setdiff(remaining, tested)
    message(sprintf("Iteration %d: %d parameter(s) to test", iter, length(to_test)))

    for (p in to_test) {
      run <- run + 1L
      res <- tryCatch(
        test_variability_one(current_project, level, p, run, iter, output_dir, reference_fit,
                             alpha, linearization, scenario, warm_start, export_charts, max_rse),
        error = function(e) {
          warning(run_name(run, iter, p, level), " failed: ", conditionMessage(e), call. = FALSE)
          list(parameter = p, level = level, run = run, failed = TRUE, status = conditionMessage(e),
               path = file.path(output_dir, paste0(run_name(run, iter, p, level), ".mlxtran")))
        })
      res$iter <- iter

      if (isTRUE(res$failed)) {
        log_variability(output_dir, variability_log_row(run, iter, p, level, NULL, NA_real_, res$path,
                                                         alpha = alpha, status = paste("error:", res$status)))
      } else {
        thresh <- critical_ofv(alpha, max(res$df, 1))
        verdict <- if (res$status != "ok") sprintf("INVALID (OFV = %s)", format(res$ofv_raw))
                   else if (res$eligible) "SIGNIFICANT"
                   else if (res$significant) sprintf("significant, not adopted (RSE of %s: %s)", res$new_parameter,
                                                     if (is.finite(res$rse_new)) sprintf("%.0f%% > %g%%", res$rse_new, max_rse) else "not estimated")
                   else "not significant"
        message(sprintf("  run%03d  %s on %s: delta OFV = %s (critical %.2f, df = %d), p = %s -> %s%s [%s]",
                        run, p, level, if (is.finite(res$delta_ofv)) sprintf("%.2f", res$delta_ofv) else "NA",
                        -thresh, as.integer(res$df), fmt_pvalue(res$p_value), verdict,
                        if (isTRUE(res$delta_ofv > 0)) " (OFV increased: estimation not as good as the reference)" else "",
                        fmt_duration(res$elapsed)))
        log_variability(output_dir, variability_log_row(
          run, iter, p, level, res$fit, res$elapsed, res$path, alpha = alpha, delta_ofv = res$delta_ofv,
          df = res$df, critical = thresh, p_value = res$p_value, significant = res$significant,
          rse_new = res$rse_new, eligible = res$eligible, status = res$status))
      }
      cands[[length(cands) + 1]] <- res
      save_state(finished = FALSE, result = NULL)   # a resumed run restarts after this model
    }

    # Adopted: the significant model with the lowest p-value (ties: the lowest OFV)
    this_iter <- Filter(function(x) x$iter == iter && !isTRUE(x$failed), cands)
    sig <- Filter(function(x) isTRUE(x$eligible), this_iter)
    chosen <- if (length(sig) > 0) {
      sig[[order(vapply(sig, function(x) x$p_value, numeric(1)),
                 vapply(sig, function(x) crit_value(x$fit, "OFV"), numeric(1)))[1]]]
    } else NULL

    # Pre-selection: after the first iteration, only the parameters that showed something stay
    if (forward && iter == 1 && !is.null(prune_alpha)) {
      kept    <- vapply(Filter(function(x) isTRUE(x$p_value < prune_alpha), this_iter), function(x) x$parameter, character(1))
      dropped <- setdiff(remaining, kept)
      if (length(dropped) > 0) {
        message(sprintf("Pre-selection (p < %s against the base model): %s kept; not tested again: %s.", prune_alpha,
                        if (length(kept)) paste(kept, collapse = ", ") else "none", paste(dropped, collapse = ", ")))
        remaining <- setdiff(remaining, dropped)
      }
    }
    if (!forward) { iter_done <- iter; break }
    if (is.null(chosen)) {
      stop_reason <- sprintf("iteration %d: no random effect is significant%s (alpha = %s)", iter,
                             if (is.null(max_rse)) "" else " and well estimated", alpha)
      message("Stopping at ", stop_reason, ".")
      iter_done <- iter
      break
    }
    current_project <- chosen$path
    reference_fit   <- chosen$fit
    remaining       <- setdiff(remaining, chosen$parameter)
    adopted         <- c(adopted, chosen$run)
    iter_done       <- iter
    history[[length(history) + 1]] <- data.frame(
      iter = iter, run = chosen$run, parameter = chosen$parameter, level = level,
      ofv = crit_value(chosen$fit, "OFV"), bicc = crit_value(chosen$fit, "BICc"),
      project = basename(chosen$path), stringsAsFactors = FALSE)
    save_state(finished = FALSE, result = NULL)
    message(sprintf("Iteration %d done: added a random effect on %s (run%03d, p = %s), OFV %.2f",
                    iter, chosen$parameter, chosen$run, fmt_pvalue(chosen$p_value), crit_value(chosen$fit, "OFV")))
  }

  runs <- variability_runs_table(base_row, cands, adopted)
  result <- list(runs = runs,
                 history = tibble::as_tibble(do.call(rbind, history)),
                 final_project = current_project,
                 stop_reason = stop_reason,
                 settings = settings[c("level", "alpha", "forward", "prune_alpha", "max_rse")])
  save_state(finished = TRUE, result = result)

  load_project_quiet(current_project)
  message(sprintf("Variability test finished in %s (%s). Final model: %s",
                  fmt_duration(as.numeric(difftime(Sys.time(), total_t0, units = "secs"))),
                  stop_reason, basename(current_project)))
  invisible(result)
}
