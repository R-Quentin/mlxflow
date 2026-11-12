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
#'   (one element per run), or the path of an .rds file holding it.
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
    out$se  <- if (length(se) == nrow(out)) as.numeric(se) else NA_real_
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
#'   (one element per run), or the path of an .rds file holding it.
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
    if (!file.exists(x)) stop("File not found: ", x, call. = FALSE)
    x <- readRDS(x)
  }
  if (!is.list(x) || length(x) == 0) {
    stop("`results` must be the output of getAssessmentResults() (a list with one element per run) ",
         "or the path of an .rds file holding it.", call. = FALSE)
  }
  x
}
