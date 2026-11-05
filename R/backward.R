  ## A) Session setup and console control   --------------------------------------

# Lines already shown during the current backward run (see capture_mlx)
mlx_console <- new.env()
#' Evaluate an expression while capturing everything Monolix prints (stdout and
#' message/warning stream). show = "none" discards it, show = "new" prints only
#' the lines that were not already displayed earlier in the run, so a repeated
#' engine notice appears once instead of at every refit. If the expression
#' fails, the captured text is printed before the error is re-raised.
#' @noRd
capture_mlx <- function(expr, show = c("none", "new")) {
  show <- match.arg(show)
  con <- textConnection("mlx_buffer", "w", local = TRUE)
  sink(con); sink(con, type = "message")

  finish <- function() {
    sink(type = "message"); sink()
    lines <- trimws(textConnectionValue(con))
    close(con)
    lines[nzchar(lines)]
  }

  value <- tryCatch(expr, error = function(e) {
    captured <- finish()
    if (length(captured) > 0) message(paste(captured, collapse = "\n"))
    stop(e)
  })

  captured <- finish()
  if (show == "new") {
    fresh <- setdiff(unique(captured), mlx_console$seen)
    mlx_console$seen <- c(mlx_console$seen, fresh)
    if (length(fresh) > 0) message(paste(fresh, collapse = "\n"))
  }
  invisible(value)
}
#' .lx("loadProject") without the data set notices Monolix prints at every (re)read
#' (e.g. "subjects-occasion without observation are removed"). Call
#' .lx("loadProject") directly if you ever need to see them.
#' @noRd
load_project_quiet <- function(path) {
  capture_mlx(.lx("loadProject", path), show = "none")
}
#' Initialize the Lixoft connectors and load a Monolix project.
#' software_path is only passed on if given - .lx("initializeLixoftConnectors")
#' expects a string for "path" and errors out on an explicit NULL.
#' @param project_path Path of the Monolix project (.mlxtran) to load.
#' @param software_path Path of the MonolixSuite installation. NULL (default): lixoftConnectors locates it by itself.
#' @export
init_monolix <- function(project_path, software_path = NULL) {
  if (is.null(software_path)) {
    .lx("initializeLixoftConnectors", software = "monolix", force = TRUE)
  } else {
    .lx("initializeLixoftConnectors", software = "monolix", path = software_path, force = TRUE)
  }
  load_project_quiet(project_path)
  invisible(TRUE)
}



  ## B) Small utilities             --------------------------------------------
#' One criterion (OFV, BICc, ...) of a fit returned by run_full_estimation()
#' @noRd
crit_value <- function(fit, name) {
  v <- suppressWarnings(as.numeric(unlist(fit$criteria[name])))
  if (length(v) == 0) NA_real_ else v[1]
}
#' Display formatting - the wording of a relation is defined here only
#'
#' @noRd
fmt_relation <- function(parameter, covariate) sprintf("%s with %s", parameter, covariate)
fmt_pvalue   <- function(p) format.pval(p, digits = 2, eps = 1e-4)
#' fmt_duration
#'
#' @noRd
fmt_duration <- function(secs) {
  if (secs < 90) sprintf("%.0f s", secs)
  else if (secs < 5400) sprintf("%.1f min", secs / 60)
  else sprintf("%.1f h", secs / 3600)
}
#' Name of a run (project file and result folder): a strictly incremental run
#' number, unique over the whole backward like Monolix's ModelBuilding_NNN, the
#' iteration it belongs to (I<k>) and the relation tested, e.g.
#' "run000_I0_base", "run001_I1_V_WT", "run002_I1_Cl_WT", "run003_I2_Cl_WT".
#' The LRT verdict is never part of the name (it is only known after the run):
#' it is logged in backward_log.csv. Labels are reduced to [A-Za-z0-9_] and
#' truncated - the run number alone guarantees uniqueness, and short names keep
#' the result paths under the Windows 260-character limit.
#' @noRd
run_name <- function(run, iter, parameter = NULL, covariate = NULL) {
  label <- if (is.null(parameter)) "base" else paste(parameter, covariate, sep = "_")
  label <- substr(gsub("[^A-Za-z0-9_]", "_", label), 1, 30)
  sprintf("run%03d_I%d_%s", run, iter, label)
}
#' Rounding applied to what is written in the CSV files (console and returned
#' objects keep full precision): OFV-like values 2 decimals, p-values 3
#' significant digits (so very small p-values do not collapse to 0), times
#' 1 decimal
#' @noRd
round_cols <- function(x) {
  x |>
    mutate(across(any_of(c("delta_ofv", "critical_ofv", "ofv", "bicc")), \(v) round(v, 2)),
           across(any_of("p_value"), \(v) signif(v, 3)),
           across(any_of("elapsed_sec"), \(v) round(v, 1)))
}
#' One row of backward_log.csv (run 0 = base model, then one row per candidate)
#' @noRd
log_row <- function(run, iter, fit, elapsed_sec, project, alpha = NA_real_,
                    parameter = NA_character_, covariate = NA_character_,
                    delta_ofv = NA_real_, df = NA_real_, critical = NA_real_,
                    p_value = NA_real_, significant = NA) {
  data.frame(
    run = run, iter = iter, parameter = parameter, covariate = covariate,
    delta_ofv = delta_ofv, df = df, critical_ofv = critical, alpha = alpha,
    p_value = p_value, significant = significant,
    ofv = crit_value(fit, "OFV"), bicc = crit_value(fit, "BICc"),
    elapsed_sec = elapsed_sec, project = basename(project)
  )
}
#' Append one row to the persistent CSV trace of every tested model - this
#' makes the search auditable after the fact, or while it is still running
#' as a background job
#' @noRd
log_candidate <- function(output_dir, row) {
  log_file <- file.path(output_dir, "backward_log.csv")
  write.table(round_cols(row), log_file, sep = ";", row.names = FALSE,
              col.names = !file.exists(log_file), append = file.exists(log_file))
}
#' One row of the table of ADOPTED models (iteration 0 = base model)
#' @noRd
history_row <- function(run, iter, fit, project,
                        parameter = NA_character_, covariate = NA_character_,
                        delta_ofv = NA_real_, df = NA_real_, p_value = NA_real_,
                        n_candidates_tested = NA_integer_) {
  data.frame(
    iter = iter, run = run, parameter = parameter, covariate = covariate,
    delta_ofv = delta_ofv, df = df, p_value = p_value,
    ofv = crit_value(fit, "OFV"), bicc = crit_value(fit, "BICc"),
    n_candidates_tested = n_candidates_tested, project = basename(project)
  )
}
#' Call the user's diagnostic callback without ever letting it kill a run:
#' hours of estimation must not be lost to a typo in a plotting function
#' @noRd
.call_diagnostic_fun <- function(diagnostic_fun, iter, project_file) {
  if (is.null(diagnostic_fun)) return(invisible(NULL))
  tryCatch(
    diagnostic_fun(iter, project_file),
    error = function(e) warning(sprintf(
      "diagnostic_fun failed at iteration %d (the run continues): %s",
      iter, conditionMessage(e)), call. = FALSE, immediate. = TRUE)
  )
  invisible(NULL)
}
#' Files and folders produced by a previous backward run in output_dir (current
#' naming, plus the "stepNN_*" / backward_steps.csv of earlier versions). Only
#' names created by this script are matched - nothing else in the folder is
#' ever touched.
#' @noRd
list_backward_outputs <- function(output_dir) {
  list.files(output_dir, full.names = TRUE, pattern = paste0(
    "^(run[0-9]{3}_.*|step[0-9]{2}_.*|backward_state\\.rds|backward_log\\.csv|",
    "backward_iterations\\.csv|backward_steps\\.csv)$"))
}
#' When resuming: discard what an interrupted iteration left behind (project
#' files, result folders and log rows of runs after the last completed one)
#' @noRd
clean_incomplete_runs <- function(output_dir, last_run) {
  items <- list.files(output_dir, pattern = "^run[0-9]{3}_", full.names = TRUE)
  stale <- items[as.integer(substr(basename(items), 4, 6)) > last_run]
  unlink(stale, recursive = TRUE)

  log_file <- file.path(output_dir, "backward_log.csv")
  if (file.exists(log_file)) {
    old_log <- read.table(log_file, sep = ";", header = TRUE, stringsAsFactors = FALSE)
    write.table(old_log[old_log$run <= last_run, ], log_file, sep = ";", row.names = FALSE)
  }
  invisible(length(stale))
}



  ## C) Covariate model helpers     --------------------------------------------
#' Tidy view of the current parameter-covariate relations
#' One row per (parameter, covariate) pair with a logical "included" flag
#' @export
get_covariate_relations <- function() {
  cov_model <- .lx("getIndividualParameterModel")$covariateModel
  do.call(rbind, lapply(names(cov_model), function(param) {
    covs <- cov_model[[param]]
    data.frame(
      parameter = param,
      covariate = names(covs),
      included  = as.logical(covs),
      stringsAsFactors = FALSE
    )
  }))
}
#' Print the active parameter-covariate relations of the currently loaded
#' project - useful as a sanity check before launching a long run
#' @export
print_covariate_relations <- function() {
  active <- get_covariate_relations()
  active <- active[active$included, ]
  if (nrow(active) == 0) {
    message("No active covariate relation in the current model.")
  } else {
    message(sprintf("Active covariate relations (%d):", nrow(active)))
    message(paste0("  - ", fmt_relation(active$parameter, active$covariate), collapse = "\n"))
  }
  invisible(active)
}
#' Switch a single parameter-covariate relation on/off in the currently
#' loaded project (only edits the in-memory model, does not run anything)
#' @noRd
set_covariate_relation <- function(parameter, covariate, included) {
  cov_model <- .lx("getIndividualParameterModel")$covariateModel
  cov_model[[parameter]][covariate] <- included
  do.call(.lx, c(list("setCovariateModel"), cov_model))
  invisible(cov_model)
}



  ## D) Full estimation chain       --------------------------------------------
#' Run SAEM, conditional distribution sampling and log-likelihood estimation
#' (importance sampling by default), then collect what is needed to compare
#' two nested models. Engine notices are shown once per backward run.
#' @noRd
run_full_estimation <- function(linearization = FALSE) {
  capture_mlx({
    .lx("runPopulationParameterEstimation")
    .lx("runConditionalDistributionSampling")   # required before importance-sampling LL
    .lx("runLogLikelihoodEstimation", linearization = linearization)
  }, show = "new")

  ll_all <- .lx("getEstimatedLogLikelihood")
  method <- if (linearization) "linearization" else "importanceSampling"
  pop_params <- .lx("getEstimatedPopulationParameters")

  list(
    pop_params = pop_params,
    criteria   = ll_all[[method]],   # named vector: OFV, AIC, BIC, BICc
    n_params   = length(pop_params)  # used to infer LRT degrees of freedom
  )
}



  ## E) Wald pre-screening          --------------------------------------------
#' Rank the currently active parameter-covariate relations by Wald p-value,
#' using the standard errors of the CURRENTLY LOADED model (no refit needed).
#' Only used to decide in which order to LRT-test candidates at each backward
#' step - the accept/reject decision is always made by the LRT, never by the
#' Wald test alone.
#' Requires .lx("runStandardErrorEstimation") to have been run beforehand.
#' @noRd
wald_prescreen <- function() {
  se_tab <- .lx("getEstimatedStandardErrors")$stochasticApproximation  # columns: parameter, se, rse
  pop    <- .lx("getEstimatedPopulationParameters")

  relations <- get_covariate_relations()
  active <- relations[relations$included, ]

  active$wald_p <- vapply(seq_len(nrow(active)), function(i) {
    # beta coefficient(s) for this relation - a categorical covariate with
    # several levels has several "beta_<param>_<covariate>_<level>" rows
    prefix <- sprintf("^beta_%s_%s", active$parameter[i], active$covariate[i])
    rows <- se_tab[grepl(prefix, se_tab$parameter), ]
    if (nrow(rows) == 0) return(NA_real_)
    z <- pop[rows$parameter] / rows$se
    # conservative screening score: the smallest (most significant) p-value
    # among the dummy coefficients of a categorical covariate
    min(2 * (1 - pnorm(abs(z))), na.rm = TRUE)
  }, numeric(1))

  active[order(-active$wald_p), ]   # least significant (safest to drop) first
}



  ## F) Likelihood Ratio Test       --------------------------------------------
#' Remove one parameter-covariate relation from the currently loaded project,
#' refit, and compare with a reference (nested) model via a chi-square LRT.
#' Degrees of freedom are inferred from the change in the number of estimated
#' population parameters, so categorical covariates with several dummy
#' coefficients are handled automatically.
#' @noRd
test_covariate_removal <- function(parameter, covariate, reference_fit, linearization = FALSE) {
  t0 <- Sys.time()
  set_covariate_relation(parameter, covariate, FALSE)
  fit <- run_full_estimation(linearization = linearization)

  df        <- reference_fit$n_params - fit$n_params
  delta_ofv <- fit$criteria[["OFV"]] - reference_fit$criteria[["OFV"]]
  p_value   <- 1 - pchisq(delta_ofv, df = max(df, 1))
  elapsed   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  list(parameter = parameter, covariate = covariate,
       df = df, delta_ofv = delta_ofv, p_value = p_value,
       elapsed = elapsed, fit = fit)
}



  ## G) Backward elimination driver --------------------------------------------
#' Backward covariate elimination for a Monolix project, with per-candidate
#' logging, resume support and optional Wald pre-screening to cut down the
#' number of refits per iteration.
#'
#' The base model is first re-estimated and saved as run000_I0_base, so that
#' the reference OFV is computed exactly like every candidate's and the
#' original project (and its result folder) is left untouched. Every model
#' estimated afterwards gets the next run number, whatever its iteration:
#' "run<NNN>_I<k>_<parameter>_<covariate>" is the candidate where that relation
#' is removed, at iteration k (run numbers never repeat, like Monolix's
#' ModelBuilding_NNN). The chosen candidate of an iteration is never refit a
#' second time - it is promoted directly from the run already performed while
#' testing it, and each refit is warm-started from the saved estimates of the
#' previous iteration.
#'
#' Files written in output_dir: run<NNN>_*.mlxtran (+ result folders),
#' backward_log.csv (every tested model, run 0 = base model),
#' backward_iterations.csv (adopted models only) and backward_state.rds
#' (checkpoint).
#'
#' @param project_file Path of the Monolix project (.mlxtran) holding the covariate model to reduce, e.g. the final model of a COSSAC or SCM search.
#' @param linearization If TRUE, the log-likelihood is computed by linearization instead of importance sampling. Every model of the run is evaluated the same way, so the OFV differences stay comparable.
#' @param alpha               LRT level to KEEP a relation (Monolix's own SCM
#'                             backward step defaults to 0.01). The critical OFV
#'                             increase is critical_ofv(alpha, df), df being the
#'                             number of coefficients removed with the relation
#' @param covariates_to_test  optional "parameter:covariate" vector restricting
#'                             which relations are allowed to be dropped
#' @param diagnostic_fun      optional callback(iteration, project_file), called
#'                             for the base model (iteration 0) and after each
#'                             accepted removal. The model is still loaded in the
#'                             Monolix engine at that point, so connectors getters
#'                             (e.g. getEtaShrinkage(), MonolixSuite 2024+) can be
#'                             used directly. An error inside the callback is
#'                             reported as a warning and does not stop the run.
#' @param output_dir          folder for everything listed above
#' @param fast_screening      if TRUE (default), candidates are tested in
#'                             Wald-ranked order (least significant first) and
#'                             the iteration stops at the first LRT-confirmed
#'                             non-significant one instead of testing every
#'                             active relation. Set to FALSE for an exhaustive
#'                             search (closer to classic SCM), e.g. for a final
#'                             confirmatory run. Both modes adopt the candidate
#'                             with the highest LRT p-value among those tested.
#' @param save_all            if TRUE (default), every candidate is kept on disk.
#'                             If FALSE, every candidate but the one adopted for
#'                             the iteration is deleted once the iteration
#'                             completes: saves disk space, but a flowchart built
#'                             from the folders then only shows the adopted path
#'                             (backward_log.csv still holds every candidate)
#' @param export_charts       if TRUE, calls .lx("computeChartsData") (VPC, scatter
#'                             plots, ...) before saving - .lx("saveProject") alone
#'                             only writes the numeric tables (summary.txt,
#'                             populationParameters.txt, ...), not the charts.
#'                             Adds computation time, especially combined with
#'                             save_all and fast_screening = FALSE
#' @param resume              if TRUE and output_dir holds an interrupted run,
#'                             continue it from its last completed iteration
#'                             (whatever an interrupted iteration left behind is
#'                             discarded). If that run already finished, its
#'                             stored result is returned without recomputing
#' @param overwrite           if TRUE, delete the outputs of any previous run in
#'                             output_dir (only the files and folders created by
#'                             this script, nothing else) and start from scratch
#' @param max_iter            maximum number of iterations
#' @export
backward_elimination <- function(project_file,
                                  alpha = 0.01,
                                  covariates_to_test = NULL,
                                  diagnostic_fun = NULL,
                                  output_dir = tempdir(),
                                  linearization = FALSE,
                                  fast_screening = TRUE,
                                  save_all = TRUE,
                                  export_charts = FALSE,
                                  resume = FALSE,
                                  overwrite = FALSE,
                                  max_iter = 30) {

  total_t0 <- Sys.time()
  mlx_console$seen <- character()

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_dir   <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  project_file <- normalizePath(project_file, winslash = "/", mustWork = TRUE)
  state_file   <- file.path(output_dir, "backward_state.rds")
  steps_file   <- file.path(output_dir, "backward_iterations.csv")

  # Sanity check: never silently mix a new run with the outputs of a previous one
  previous <- list_backward_outputs(output_dir)
  if (overwrite && length(previous) > 0) {
    unlink(previous, recursive = TRUE)
    message(sprintf("Removed %d file(s)/folder(s) from a previous backward run.", length(previous)))
    previous <- character()
  }
  if (length(previous) > 0) {
    if (!resume) {
      stop("output_dir already contains the outputs of a previous backward run. ",
           "Use resume = TRUE to continue it, or overwrite = TRUE to delete it ",
           "and start over.", call. = FALSE)
    }
    if (!file.exists(state_file)) {
      stop("output_dir contains backward outputs but no checkpoint (backward_state.rds), ",
           "so the run cannot be resumed safely. Use overwrite = TRUE to start over.",
           call. = FALSE)
    }
  }

  if (length(previous) > 0) {
    state <- readRDS(state_file)
    if (isTRUE(state$finished)) {
      message(sprintf("This backward run is already complete (%s). Returning the stored result; ",
                      state$result$stop_reason),
              "use overwrite = TRUE to run it again.")
      return(state$result)
    }
    clean_incomplete_runs(output_dir, state$last_run)
    current_project <- state$current_project
    reference_fit   <- state$reference_fit
    history         <- state$history
    last_iter       <- state$iter
    run             <- state$last_run
    load_project_quiet(current_project)
    message(sprintf("Resuming after iteration %d (model: %s).", last_iter, basename(current_project)))
    print_covariate_relations()
  } else {
    # Iteration 0: base model, re-estimated and saved in output_dir
    load_project_quiet(project_file)
    message("Iteration 0 (base model)")
    print_covariate_relations()

    run       <- 0L
    last_iter <- 0L
    base_name <- run_name(run, 0)
    current_project <- file.path(output_dir, paste0(base_name, ".mlxtran"))
    .lx("setProjectSettings", directory = file.path(output_dir, base_name))

    t0 <- Sys.time()
    reference_fit <- run_full_estimation(linearization = linearization)
    base_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (export_charts) .lx("computeChartsData")
    .lx("saveProject", current_project)

    message(sprintf("  %s fitted in %s: OFV = %.2f, BICc = %.2f",
                    base_name, fmt_duration(base_secs), crit_value(reference_fit, "OFV"),
                    crit_value(reference_fit, "BICc")))
    log_candidate(output_dir, log_row(run, 0, reference_fit, base_secs, current_project, alpha = alpha))

    history <- list(history_row(run, 0, reference_fit, current_project))
    write.table(round_cols(do.call(rbind, history)), steps_file, sep = ";", row.names = FALSE)
    saveRDS(list(current_project = current_project, reference_fit = reference_fit,
                 history = history, iter = last_iter, last_run = run, finished = FALSE),
            state_file)
    .call_diagnostic_fun(diagnostic_fun, 0, current_project)
  }

  stop_reason <- sprintf("max_iter (%d) reached", max_iter)
  iters <- if (last_iter + 1 <= max_iter) (last_iter + 1):max_iter else integer(0)

  for (iter in iters) {

    relations <- get_covariate_relations()
    active <- relations[relations$included, ]
    if (!is.null(covariates_to_test)) {
      key <- paste(active$parameter, active$covariate, sep = ":")
      active <- active[key %in% covariates_to_test, ]
    }
    if (nrow(active) == 0) {
      stop_reason <- "no covariate relation left to test"
      break
    }

    # order candidates: Wald pre-screening (fast) or as-is (exhaustive)
    if (fast_screening) {
      capture_mlx(.lx("runStandardErrorEstimation", linearization = FALSE), show = "new")
      active <- wald_prescreen()
      if (!is.null(covariates_to_test)) {
        key <- paste(active$parameter, active$covariate, sep = ":")
        active <- active[key %in% covariates_to_test, ]
      }
    }

    message(sprintf("Iteration %d: %d candidate(s) to test", iter, nrow(active)))
    iter_t0 <- Sys.time()
    candidates <- list()
    chosen <- NULL

    for (i in seq_len(nrow(active))) {
      load_project_quiet(current_project)

      # Each candidate gets the next run number (never reused) and is named
      # after its iteration and the relation tested - the LRT verdict is only
      # known afterward and never goes into the name, look it up in
      # backward_log.csv instead. The results directory is always redirected
      # to this candidate's own folder before estimation, so every run writes
      # straight to its final location: nothing is ever overwritten, and the
      # chosen candidate can be promoted to the iteration's official model
      # below with no extra refit.
      run       <- run + 1L
      cand_name <- run_name(run, iter, active$parameter[i], active$covariate[i])
      cand_path <- file.path(output_dir, paste0(cand_name, ".mlxtran"))
      .lx("setProjectSettings", directory = file.path(output_dir, cand_name))

      res <- test_covariate_removal(active$parameter[i], active$covariate[i],
                                     reference_fit, linearization)
      res$path <- cand_path
      res$run  <- run
      candidates[[i]] <- res

      thresh <- critical_ofv(alpha, max(res$df, 1))
      message(sprintf("  run%03d  %s: delta OFV = %.2f (critical %.2f, df = %d), p = %s -> %s [%s]",
                      run, fmt_relation(res$parameter, res$covariate), res$delta_ofv, thresh,
                      as.integer(res$df), fmt_pvalue(res$p_value),
                      if (res$p_value < alpha) "SIGNIFICANT" else "removable",
                      fmt_duration(res$elapsed)))

      log_candidate(output_dir, log_row(
        run, iter, res$fit, res$elapsed, cand_path, alpha = alpha,
        parameter = res$parameter, covariate = res$covariate,
        delta_ofv = res$delta_ofv, df = res$df, critical = thresh,
        p_value = res$p_value, significant = res$p_value < alpha))

      if (export_charts) .lx("computeChartsData")
      .lx("saveProject", cand_path)   # always saved - this is what lets the chosen
                                # candidate become the iteration's official model
                                # with no separate confirmation refit

      if (fast_screening && res$p_value >= alpha) {
        chosen <- res   # first non-significant candidate found: stop testing
        break
      }
    }

    # Adopted candidate = highest LRT p-value among those tested (in fast mode
    # that is the last one, the only non-significant one). The flowchart
    # builder applies the same rule.
    if (is.null(chosen)) {
      p_values <- vapply(candidates, function(x) x$p_value, numeric(1))
      chosen <- candidates[[which.max(p_values)]]
    }

    stopping <- chosen$p_value < alpha

    # reload whichever model actually carries on (current_project if
    # stopping, chosen's otherwise) FIRST, so Monolix releases any file
    # handle on the other candidates' folders before they get discarded
    if (!stopping) current_project <- chosen$path
    load_project_quiet(current_project)

    if (!save_all) {
      keep_path <- if (stopping) NA_character_ else chosen$path
      for (cand in candidates) {
        if (!is.na(keep_path) && identical(cand$path, keep_path)) next
        unlink(cand$path)
        unlink(file.path(dirname(cand$path), tools::file_path_sans_ext(basename(cand$path))),
               recursive = TRUE)
      }
    }

    if (stopping) {
      stop_reason <- sprintf("iteration %d: the least significant relation (%s) is still significant",
                             iter, fmt_relation(chosen$parameter, chosen$covariate))
      message(sprintf("Stopping at iteration %d: %s is still significant (delta OFV = %.2f > critical %.2f, p = %s).",
                      iter, fmt_relation(chosen$parameter, chosen$covariate), chosen$delta_ofv,
                      critical_ofv(alpha, max(chosen$df, 1)), fmt_pvalue(chosen$p_value)))
      break
    }

    old_ofv <- crit_value(reference_fit, "OFV")
    reference_fit <- chosen$fit   # reuse the fit already computed while testing - no refit needed
    last_iter <- iter
    message(sprintf("Iteration %d done in %s: dropped %s (run%03d, p = %s), OFV %.2f -> %.2f [%d/%d candidates tested]",
                    iter, fmt_duration(as.numeric(difftime(Sys.time(), iter_t0, units = "secs"))),
                    fmt_relation(chosen$parameter, chosen$covariate), chosen$run,
                    fmt_pvalue(chosen$p_value), old_ofv, crit_value(reference_fit, "OFV"),
                    length(candidates), nrow(active)))

    history[[length(history) + 1]] <- history_row(
      chosen$run, iter, reference_fit, current_project,
      parameter = chosen$parameter, covariate = chosen$covariate,
      delta_ofv = chosen$delta_ofv, df = chosen$df, p_value = chosen$p_value,
      n_candidates_tested = length(candidates))
    write.table(round_cols(do.call(rbind, history)), steps_file, sep = ";", row.names = FALSE)

    saveRDS(list(current_project = current_project, reference_fit = reference_fit,
                 history = history, iter = last_iter, last_run = run, finished = FALSE),
            state_file)

    .call_diagnostic_fun(diagnostic_fun, iter, current_project)
  }

  result <- list(final_project = current_project,
                 final_fit     = reference_fit,
                 history       = do.call(rbind, history),
                 stop_reason   = stop_reason)
  saveRDS(list(current_project = current_project, reference_fit = reference_fit,
               history = history, iter = last_iter, last_run = run, finished = TRUE,
               result = result),
          state_file)

  message(sprintf("Backward elimination finished in %s (%s). Final model: %s",
                  fmt_duration(as.numeric(difftime(Sys.time(), total_t0, units = "secs"))),
                  stop_reason, basename(current_project)))
  if (!save_all) {
    message("Note: save_all = FALSE, rejected candidates were deleted. A flowchart built ",
            "from the folders will only show the adopted path; backward_log.csv holds every candidate.")
  }
  result
}



  ## H) Optional COSSAC search before the manual backward check ----------------
#' Run a COSSAC covariate search on the currently loaded project and apply
#' the winning covariate model explicitly (safety net in case the loaded
#' project was not already updated).
#' @param covariates Optional character vector restricting the search to these covariates (default: the model building settings of the project).
#' @param use_linearization If TRUE, the search relies on linearization (faster, less accurate) instead of importance sampling.
#' @export
run_cossac <- function(covariates = NULL, use_linearization = FALSE) {

  mb_settings <- .lx("getModelBuildingSettings")
  mb_settings$strategy <- "cossac"
  mb_settings$useLin   <- use_linearization
  if (!is.null(covariates)) mb_settings$covariates <- covariates
  # inspect str(mb_settings) to see the exact stopping-criterion / threshold
  # fields available in your MonolixSuite version before customizing further

  .lx("runModelBuilding", settings = mb_settings)

  results <- .lx("getModelBuildingResults")
  best <- results[[which(vapply(results, function(x) isTRUE(x$bestModel), logical(1)))[1]]]

  im <- best$individualModels   # parameter x covariate logical matrix
  cov_list <- setNames(lapply(rownames(im), function(p) im[p, ]), rownames(im))
  do.call(.lx, c(list("setCovariateModel"), cov_list))

  list(results = results, best_model = im)
}
