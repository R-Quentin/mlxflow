  ## Testing the addition of variability (IOV or IIV) on the parameters of a model ---------

# Files and folders created by test_variability() in output_dir (same run<NNN>_ naming as
# backward_elimination(), plus its own log). Only these names are ever touched.
#' @noRd
list_variability_outputs <- function(output_dir) {
  list.files(output_dir, full.names = TRUE,
             pattern = "^(run[0-9]{3}_.*|variability_log\\.csv)$")
}
#' One row of variability_log.csv (run 0 = base model, then one row per tested model)
#' @noRd
variability_log_row <- function(run, iter, parameter, level, fit, elapsed_sec, project,
                                alpha = NA_real_, delta_ofv = NA_real_, df = NA_real_,
                                critical = NA_real_, p_value = NA_real_, significant = NA,
                                status = "ok") {
  data.frame(
    run = run, iter = iter, parameter = parameter, level = level,
    delta_ofv = delta_ofv, df = df, critical_ofv = critical, alpha = alpha,
    p_value = p_value, significant = significant,
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
#' Estimate one model of the search in its own folder, save it, and compare it with the reference
#' by a likelihood ratio test. The degrees of freedom are the number of population parameters
#' that the change added (1 for the variance of one random effect).
#' @noRd
test_variability_one <- function(parent_project, level, parameter, run, iter, output_dir,
                                 reference_fit, alpha, linearization, scenario, warm_start,
                                 export_charts) {
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
  if (export_charts) .lx("computeChartsData")
  .lx("saveProject", path)

  df        <- fit$n_params - reference_fit$n_params
  delta_ofv <- crit_value(fit, "OFV") - crit_value(reference_fit, "OFV")
  # Adding parameters can only lower the OFV: a rise means the estimation did not do as well
  p_value   <- 1 - stats::pchisq(max(-delta_ofv, 0), df = max(df, 1))
  list(parameter = parameter, level = level, run = run, df = df, delta_ofv = delta_ofv,
       p_value = p_value, significant = p_value < alpha, fit = fit, path = path,
       elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
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
#' estimation chain. Every model is saved as `run<NNN>_I<k>_<parameter>_<level>.mlxtran`
#' with its own result folder (`NNN` is never reused, `k` is the iteration), so the original
#' project and its results are left untouched, and each estimation starts from the estimates of the
#' model it derives from. `variability_log.csv` holds every test. The folders can be read
#' with [read_models()] and compared with [compare_models()] and [compare_obspred_metrics()].
#'
#' The LRT of a variance is approximate: the null value (0) is on the boundary of the parameter
#' space, so the test is conservative, and an importance-sampling likelihood is itself
#' noisy. Read the p-values together with the shrinkage and the diagnostics of the models.
#'
#' @param project_file Path of the model (.mlxtran) to start from.
#' @param output_dir Folder of the runs. It must not hold the outputs of a previous run
#'   (see `overwrite`).
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
#' @param scenario Optional Monolix scenario (a list as returned by `getScenario()`) used to
#'   estimate every model, for instance a fast one without standard errors or a complete one
#'   with all the diagnostics. Default `NULL`: SAEM, conditional distribution sampling and
#'   log-likelihood, which is what the LRT needs. See [backward_elimination()] for what the
#'   scenario must enable.
#' @param linearization If TRUE, the log-likelihood is computed by linearization instead of
#'   importance sampling. Ignored when `scenario` is given (its `linearization` is used).
#' @param warm_start If TRUE (default), each estimation starts from the last estimates of the
#'   model it derives from instead of the initial values of the project.
#' @param export_charts If TRUE (default), the charts data (VPC, observations vs predictions,
#'   ...) of every model are computed and saved, which the diagnostic functions need.
#'   FALSE saves time.
#' @param overwrite If TRUE, delete the outputs of a previous run in `output_dir` (only the
#'   files and folders created by this function) and start over.
#'
#' @return A list, invisibly: `runs` (one row per estimated model: iteration, parameter, level,
#'   OFV, BICc, `delta_ofv` and `p_value` against the reference, `significant`, `adopted`,
#'   `status`), `history` (the adopted models, base model first), `final_project` (the model
#'   the search ends on: the base model when nothing is adopted) and `stop_reason`.
#' @examples
#' \dontrun{
#' init_monolix("iov_project.mlxtran")
#'
#' # Is there inter-occasion variability, and on which parameter?
#' res <- test_variability("iov_project.mlxtran", output_dir = "iov_project/IOV_test", level = "OCC")
#' res$runs
#'
#' # Stepwise selection, with a quick scenario
#' res <- test_variability("iov_project.mlxtran", output_dir = "iov_project/IOV_forward",
#'                         level = "OCC", forward = TRUE, alpha = 0.05,
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
                             scenario      = NULL,
                             linearization = FALSE,
                             warm_start    = TRUE,
                             export_charts = TRUE,
                             overwrite     = FALSE) {

  total_t0 <- Sys.time()
  mlx_console$seen <- character()
  if (!is.null(scenario)) {
    check_scenario(scenario)
    linearization <- isTRUE(scenario$linearization)
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_dir   <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  project_file <- normalizePath(project_file, winslash = "/", mustWork = TRUE)

  # Never silently mix a new run with the outputs of a previous one
  previous <- list_variability_outputs(output_dir)
  if (length(previous) > 0) {
    if (!overwrite) {
      stop("output_dir already contains the outputs of a previous run. Use overwrite = TRUE to ",
           "delete them and start over, or another output_dir.", call. = FALSE)
    }
    unlink(previous, recursive = TRUE)
    message(sprintf("Removed %d file(s)/folder(s) from a previous run.", length(previous)))
  }

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

  # Run 0: the base model, estimated like every candidate
  run       <- 0L
  base_name <- run_name(run, 0)
  current_project <- file.path(output_dir, paste0(base_name, ".mlxtran"))
  message("Run 0 (base model)")
  .lx("setProjectSettings", directory = file.path(output_dir, base_name))
  t0 <- Sys.time()
  reference_fit <- run_full_estimation(linearization = linearization, scenario = scenario)
  base_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (export_charts) .lx("computeChartsData")
  .lx("saveProject", current_project)
  message(sprintf("  %s fitted in %s: OFV = %.2f, BICc = %.2f", base_name, fmt_duration(base_secs),
                  crit_value(reference_fit, "OFV"), crit_value(reference_fit, "BICc")))
  log_variability(output_dir, variability_log_row(run, 0, NA_character_, level, reference_fit,
                                                   base_secs, current_project, alpha = alpha))

  runs <- list(data.frame(run = run, iter = 0L, parameter = NA_character_, level = level,
                          ofv = crit_value(reference_fit, "OFV"), bicc = crit_value(reference_fit, "BICc"),
                          delta_ofv = NA_real_, df = NA_real_, p_value = NA_real_, significant = NA,
                          adopted = NA, status = "ok", elapsed_sec = base_secs,
                          project = basename(current_project), stringsAsFactors = FALSE))
  history <- list(runs[[1]][, c("iter", "run", "parameter", "level", "ofv", "bicc", "project")])
  stop_reason <- "single pass (forward = FALSE)"
  remaining   <- parameters

  for (iter in seq_len(if (forward) max_iter else 1L)) {
    if (length(remaining) == 0) { stop_reason <- "no parameter left to test"; break }
    message(sprintf("Iteration %d: %d parameter(s) to test", iter, length(remaining)))
    candidates <- list()

    for (p in remaining) {
      run <- run + 1L
      res <- tryCatch(
        test_variability_one(current_project, level, p, run, iter, output_dir, reference_fit,
                             alpha, linearization, scenario, warm_start, export_charts),
        error = function(e) {
          warning(run_name(run, iter, p, level), " failed: ", conditionMessage(e), call. = FALSE)
          list(parameter = p, level = level, run = run, failed = TRUE, status = conditionMessage(e),
               path = file.path(output_dir, paste0(run_name(run, iter, p, level), ".mlxtran")))
        })

      if (isTRUE(res$failed)) {
        log_variability(output_dir, variability_log_row(run, iter, p, level, NULL, NA_real_, res$path,
                                                         alpha = alpha, status = paste("error:", res$status)))
        runs[[length(runs) + 1]] <- data.frame(
          run = run, iter = iter, parameter = p, level = level, ofv = NA_real_, bicc = NA_real_,
          delta_ofv = NA_real_, df = NA_real_, p_value = NA_real_, significant = NA, adopted = FALSE,
          status = paste("error:", res$status), elapsed_sec = NA_real_, project = basename(res$path),
          stringsAsFactors = FALSE)
        next
      }

      candidates[[length(candidates) + 1]] <- res
      thresh <- critical_ofv(alpha, max(res$df, 1))
      message(sprintf("  run%03d  %s on %s: delta OFV = %.2f (critical -%.2f, df = %d), p = %s -> %s [%s]",
                      run, p, level, res$delta_ofv, thresh, as.integer(res$df), fmt_pvalue(res$p_value),
                      if (res$significant) "SIGNIFICANT" else "not significant", fmt_duration(res$elapsed)))
      log_variability(output_dir, variability_log_row(
        run, iter, p, level, res$fit, res$elapsed, res$path, alpha = alpha, delta_ofv = res$delta_ofv,
        df = res$df, critical = thresh, p_value = res$p_value, significant = res$significant))
    }

    # Adopted: the significant model with the lowest p-value (ties: the lowest OFV)
    sig <- Filter(function(x) x$significant, candidates)
    chosen <- if (length(sig) > 0) {
      sig[[order(vapply(sig, function(x) x$p_value, numeric(1)),
                 vapply(sig, function(x) crit_value(x$fit, "OFV"), numeric(1)))[1]]]
    } else NULL

    for (cand in candidates) {
      runs[[length(runs) + 1]] <- data.frame(
        run = cand$run, iter = iter, parameter = cand$parameter, level = level,
        ofv = crit_value(cand$fit, "OFV"), bicc = crit_value(cand$fit, "BICc"),
        delta_ofv = cand$delta_ofv, df = cand$df, p_value = cand$p_value, significant = cand$significant,
        adopted = forward && !is.null(chosen) && cand$run == chosen$run, status = "ok",
        elapsed_sec = cand$elapsed, project = basename(cand$path), stringsAsFactors = FALSE)
    }

    if (!forward) break
    if (is.null(chosen)) {
      stop_reason <- sprintf("iteration %d: no random effect is significant (alpha = %s)", iter, alpha)
      message("Stopping at ", stop_reason, ".")
      break
    }
    current_project <- chosen$path
    reference_fit   <- chosen$fit
    remaining       <- setdiff(remaining, chosen$parameter)
    history[[length(history) + 1]] <- data.frame(
      iter = iter, run = chosen$run, parameter = chosen$parameter, level = level,
      ofv = crit_value(chosen$fit, "OFV"), bicc = crit_value(chosen$fit, "BICc"),
      project = basename(chosen$path), stringsAsFactors = FALSE)
    message(sprintf("Iteration %d done: added a random effect on %s (run%03d, p = %s), OFV %.2f",
                    iter, chosen$parameter, chosen$run, fmt_pvalue(chosen$p_value), crit_value(chosen$fit, "OFV")))
    if (iter == max_iter) stop_reason <- sprintf("max_iter (%d) reached", max_iter)
  }

  load_project_quiet(current_project)
  message(sprintf("Variability test finished in %s (%s). Final model: %s",
                  fmt_duration(as.numeric(difftime(Sys.time(), total_t0, units = "secs"))),
                  stop_reason, basename(current_project)))
  invisible(list(runs = tibble::as_tibble(do.call(rbind, runs)),
                 history = tibble::as_tibble(do.call(rbind, history)),
                 final_project = current_project,
                 stop_reason = stop_reason))
}
