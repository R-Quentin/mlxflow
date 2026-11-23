# Results of a Monolix convergence assessment: SAEM run several times from different initial
# values. The results come from lixoftConnectors::getAssessmentResults(): a list with one
# element per run, each holding populationParameters$convergence (one row per SAEM iteration,
# one column per parameter plus convergenceIndicator) and standardErrors$values.

#' Tidy the results of a convergence assessment
#'
#' Reshapes the output of `lixoftConnectors::getAssessmentResults()` into two
#' long tables: the SAEM trajectory of every parameter in every run, and the
#' final estimate of every parameter in every run (with its standard error and
#' 95 percent interval when available).
#'
#' @param results The object returned by `lixoftConnectors::getAssessmentResults()`
#'   (one element per run), the path of an .rds file holding it, or the `Assessment` folder
#'   written by Monolix (see [read_assessment()]).
#'
#' @return A list with `iterations` (tibble: `run`, `iteration`, `parameter`,
#'   `value`) and `final` (tibble: `run`, `parameter`, `estimate`, `se`, `ymin`, `ymax`).
#' @export
tidy_convergence_assessment <- function(results) {
  results <- .read_assessment(results)
  conv <- lapply(results, function(r) r$populationParameters$convergence)
  if (any(vapply(conv, is.null, logical(1)))) {
    stop("Each run must hold `populationParameters$convergence`: is this the output of getAssessmentResults()?",
         call. = FALSE)
  }

  # Parameter order of the first run, convergence indicator first
  param_levels <- c("convergenceIndicator", setdiff(names(conv[[1]]), "convergenceIndicator"))
  run_levels   <- as.character(seq_along(conv))

  iterations <- purrr::imap(conv, function(df, idx) {
    df |>
      dplyr::mutate(iteration = dplyr::row_number(), run = factor(idx, levels = run_levels)) |>
      tidyr::pivot_longer(cols = -c(iteration, run), names_to = "parameter", values_to = "value")
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(parameter = factor(parameter, levels = param_levels)) |>
    dplyr::select(run, iteration, parameter, value)

  final <- purrr::imap(results, function(res, idx) {
    last <- res$populationParameters$convergence
    last <- last[nrow(last), setdiff(names(last), "convergenceIndicator"), drop = FALSE]
    out  <- tibble::tibble(parameter = names(last), estimate = as.numeric(unlist(last)))
    se   <- res$standardErrors$values
    # getAssessmentResults() gives the standard errors as a plain vector in alphabetical order of the
    # parameter names (capitals first: Cl_pop, Tlag_pop, ..., a, b, beta_..., ka_pop, omega_...), not in
    # the order of the convergence columns: they are matched by name. A folder reading is named already.
    if (is.null(names(se)) && length(se) == nrow(out)) names(se) <- sort(out$parameter, method = "radix")
    out$se  <- if (!is.null(names(se))) as.numeric(se[out$parameter]) else NA_real_
    out$run <- factor(idx, levels = run_levels)
    out
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(parameter = factor(parameter, levels = param_levels),
                  ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se) |>
    dplyr::select(run, parameter, estimate, se, ymin, ymax)

  list(iterations = iterations, final = final)
}


#' Plot the results of a convergence assessment
#'
#' Two views of the same runs. `"iterations"` draws the SAEM trajectory of each
#' parameter, one line per run: when SAEM converges, the lines of all the runs
#' end on the same value. `"estimates"` shows the final estimate of each
#' parameter in each run: the points should line up.
#'
#' @param results The object returned by `lixoftConnectors::getAssessmentResults()`
#'   (one element per run), the path of an .rds file holding it, or the `Assessment` folder
#'   written by Monolix (see [read_assessment()]).
#' @param which `"iterations"` (default) or `"estimates"`.
#' @param ci If TRUE, the final estimates are drawn with their 95 percent
#'   interval (when standard errors are available). Only used for `"estimates"`.
#' @param base_size Base font size of the plot (ggplot2 theme).
#'
#' @return A ggplot object.
#' @export
plot_convergence_assessment <- function(results, which = c("iterations", "estimates"),
                                        ci = FALSE, base_size = 11) {
  which <- match.arg(which)
  tidy  <- tidy_convergence_assessment(results)

  if (which == "iterations") {
    if (dplyr::n_distinct(tidy$iterations$iteration) < 2) {
      stop("The assessment holds the final estimates only, no SAEM trajectories: they are read from ",
           "ChartsData/Saem/CvParam.txt of each run (compute the charts of the runs), or use the ",
           "output of getAssessmentResults(). Use which = \"estimates\".", call. = FALSE)
    }
    ggplot2::ggplot(tidy$iterations, ggplot2::aes(x = iteration, y = value, color = run)) +
      ggplot2::geom_line(alpha = 0.7, linewidth = 0.6) +
      ggplot2::facet_wrap(~ parameter, scales = "free_y") +
      ggplot2::labs(x = "SAEM iteration", y = NULL) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(legend.position = "none")
  } else {
    p <- ggplot2::ggplot(tidy$final, ggplot2::aes(x = run, y = estimate, color = run))
    if (ci) p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = 0.2, alpha = 0.8, na.rm = TRUE)
    p +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ parameter, scales = "free_y") +
      ggplot2::labs(x = "Run", y = "Final estimate") +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(legend.position = "none")
  }
}


# Accepts the results object or the path of an .rds file
#' @noRd
.read_assessment <- function(x) {
  if (is.character(x) && length(x) == 1) {
    if (dir.exists(x)) {
      x <- read_assessment(x)
    } else {
      if (!file.exists(x)) stop("File or folder not found: ", x, call. = FALSE)
      x <- readRDS(x)
    }
  }
  if (!is.list(x) || length(x) == 0) {
    stop("`results` must be the output of getAssessmentResults() (a list with one element per run) ",
         "or the path of an .rds file holding it.", call. = FALSE)
  }
  x
}

#' Read the results of a convergence assessment from the folder Monolix wrote
#'
#' A convergence assessment run from Monolix (GUI or `runAssessment()`) leaves an `Assessment`
#' folder with one sub-folder per run (`Assessment_001`, ...) and an `assessment.txt` that gathers
#' the final estimates. This reads that folder into the same structure as
#' `lixoftConnectors::getAssessmentResults()`, so the results can be plotted
#' ([plot_convergence_assessment()]) or included in [run_diagnostics()] without keeping the R object.
#'
#' What each run gives: its estimates and standard errors (`populationParameters.txt`), its
#' likelihood criteria when computed (`LogLikelihood/logLikelihood.txt`) and its SAEM trajectory
#' when the charts data of the run exist (`ChartsData/Saem/CvParam.txt`). Without the trajectories
#' only the final estimates are available (the `"iterations"` plot is then not possible). If the
#' run folders are missing, `assessment.txt` is used (final estimates only, no standard errors).
#'
#' @param assessment_dir The `Assessment` folder (e.g. `file.path(model_dir, "Assessment")`); the
#'   model folder is accepted as well, its `Assessment` sub-folder is then read. If an
#'   `AssessmentResults.rds` is next to the runs, in the folder or in the model folder (what
#'   `saveRDS(getAssessmentResults(), ...)` writes), it is read instead: it holds the SAEM trajectories,
#'   which the run folders do not when their charts were not computed.
#'
#' @return A list with one element per run, laid out like the output of
#'   `getAssessmentResults()` (`populationParameters$convergence`, `standardErrors`,
#'   `loglikelihood`). The attribute `"trajectories"` tells whether the SAEM iterations are there.
#' @examples
#' \dontrun{
#' res <- read_assessment("warfarin_final/Assessment")
#' plot_convergence_assessment(res, which = "estimates", ci = TRUE)
#' }
#' @export
read_assessment <- function(assessment_dir) {
  if (!dir.exists(assessment_dir)) stop("Assessment folder not found: ", assessment_dir, call. = FALSE)

  # The model folder is accepted too: its Assessment sub-folder is the one that is read
  folders <- unique(c(assessment_dir, if (dir.exists(file.path(assessment_dir, "Assessment"))) file.path(assessment_dir, "Assessment")))
  has_runs <- function(d) any(grepl("^Assessment_[0-9]+$", list.dirs(d, recursive = FALSE, full.names = FALSE)))
  run_folder <- folders[vapply(folders, has_runs, logical(1))][1]

  # An AssessmentResults.rds kept next to the runs holds the SAEM trajectories that the folder lacks
  for (d in folders) {
    rds <- file.path(d, "AssessmentResults.rds")
    if (file.exists(rds)) {
      message("Assessment read from ", rds, " (it holds the SAEM trajectories).")
      runs <- readRDS(rds)
      attr(runs, "trajectories") <- any(vapply(runs, function(r) nrow(r$populationParameters$convergence) > 1, logical(1)))
      return(runs)
    }
  }

  runs <- if (!is.na(run_folder)) {
    # Only the Assessment_NNN folders: the model's own populationParameters.txt is not a run
    pop_files <- list.files(run_folder, pattern = "^populationParameters\\.txt$", recursive = TRUE, full.names = TRUE)
    pop_files <- pop_files[grepl("^Assessment_[0-9]+$", basename(dirname(pop_files)))]
    run_no    <- as.integer(sub("^Assessment_", "", basename(dirname(pop_files))))
    lapply(pop_files[order(run_no, pop_files)], .read_assessment_run)
  } else {
    .read_assessment_summary(file.path(folders[length(folders)], "assessment.txt"))
  }
  attr(runs, "trajectories") <- any(vapply(runs, function(r) nrow(r$populationParameters$convergence) > 1, logical(1)))
  runs
}

# One run: the folder that holds its populationParameters.txt
#' @noRd
.read_assessment_run <- function(pop_file) {
  run_dir <- dirname(pop_file)
  pop <- utils::read.csv(pop_file, stringsAsFactors = FALSE, check.names = FALSE)
  se_col <- intersect(c("se_sa", "se_lin"), names(pop))
  se_col <- se_col[vapply(se_col, function(k) any(!is.na(pop[[k]])), logical(1))][1]

  # SAEM trajectory when the charts exist, otherwise the final estimates as a single row
  cv_file <- file.path(run_dir, "ChartsData", "Saem", "CvParam.txt")
  if (file.exists(cv_file)) {
    cv <- utils::read.csv(cv_file, stringsAsFactors = FALSE, check.names = FALSE)
    conv <- cv[, setdiff(names(cv), c("iteration", "phase")), drop = FALSE]
    n_expl <- sum(cv$phase == 1); n_smooth <- sum(cv$phase == 2)
  } else {
    conv <- as.data.frame(as.list(stats::setNames(pop$value, pop$parameter)), check.names = FALSE)
    conv <- cbind(convergenceIndicator = NA_real_, conv)
    n_expl <- n_smooth <- NA_real_
  }

  ll_file <- file.path(run_dir, "LogLikelihood", "logLikelihood.txt")
  ll <- if (file.exists(ll_file)) {
    x <- utils::read.csv(ll_file, stringsAsFactors = FALSE)
    stats::setNames(as.list(x[[2]]), x[[1]])
  } else list()
  ll$method <- if (file.exists(ll_file)) names(utils::read.csv(ll_file, nrows = 1))[2] else NA_character_

  list(
    populationParameters = list(nbexploratoryiterations = n_expl, nbsmoothingiterations = n_smooth,
                                convergence = conv),
    standardErrors = list(method = if (is.na(se_col)) NA_character_
                                    else c(se_sa = "stochasticApproximation", se_lin = "linearization")[[se_col]],
                          values = if (is.na(se_col)) stats::setNames(rep(NA_real_, nrow(pop)), pop$parameter)
                                   else stats::setNames(pop[[se_col]], pop$parameter)),
    loglikelihood = ll
  )
}

# Fallback: assessment.txt (Parameters, Run_1, Run_2, ...), final estimates only
#' @noRd
.read_assessment_summary <- function(file) {
  if (!file.exists(file)) stop("No run folder and no assessment.txt found in the assessment folder.", call. = FALSE)
  tab <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  tab <- tab[!tab[[1]] %in% "-2LL", , drop = FALSE]
  lapply(names(tab)[-1], function(run) {
    conv <- as.data.frame(as.list(stats::setNames(tab[[run]], tab[[1]])), check.names = FALSE)
    list(populationParameters = list(nbexploratoryiterations = NA_real_, nbsmoothingiterations = NA_real_,
                                     convergence = cbind(convergenceIndicator = NA_real_, conv)),
         standardErrors = list(method = NA_character_, values = stats::setNames(rep(NA_real_, nrow(tab)), tab[[1]])),
         loglikelihood = list())
  })
}
