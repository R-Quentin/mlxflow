# Functions: Estimation diagnostics  ===========================================
  ## A) SAEM convergence                                              ----------
#' read_convergence(): ChartsData/Saem/CvParam.txt -- one row per SAEM
#'
#' iteration, wide (one column per model parameter, plus convergenceIndicator
#' itself). Reshaped to long format (iteration, phase, parameter, value,
#' category), where `category` is a coarse classification of each parameter
#' from its name (fixed effect, covariate effect, random effect SD, IOV SD,
#' correlation, residual error, or "Convergence indicator" for
#' convergenceIndicator specifically), used to order/facet parameters
#' sensibly. A model with everything fixed still has convergenceIndicator to
#' plot, so that alone does not error -- only a fully degenerate/empty
#' CvParam.txt (neither parameters nor convergenceIndicator) does.
#'
#' @param model_dir Monolix results folder of the model (ChartsData must have been exported).
#'
#' @export
read_convergence <- function(model_dir) {
  
  f <- file.path(model_dir, "ChartsData", "Saem", "CvParam.txt")
  if (!file.exists(f)) stop("ChartsData/Saem/CvParam.txt not found in: ", model_dir)
  
  wide <- read.table(f, sep = ",", comment.char = "", header = TRUE,
                     stringsAsFactors = FALSE) |> as_tibble()
  param_cols <- setdiff(colnames(wide), c("iteration", "phase"))
  
  if (length(param_cols) == 0) {
    stop("ChartsData/Saem/CvParam.txt has no parameter or convergenceIndicator ",
         "column to plot for: ", model_dir)
  }
  
  categorize <- function(p) {
    dplyr::case_when(
      p == "convergenceIndicator" ~ "Convergence indicator",
      str_starts(p, "beta_")  ~ "Covariate effect",
      str_starts(p, "omega_") ~ "Random effect (SD)",
      str_starts(p, "gamma_") ~ "IOV (SD)",
      str_starts(p, "corr_")  ~ "Correlation",
      str_ends(p, "_pop")     ~ "Fixed effect",
      TRUE                    ~ "Residual error"
    )
  }
  category_levels <- c("Fixed effect", "Covariate effect", "Random effect (SD)",
                       "IOV (SD)", "Correlation", "Residual error", "Convergence indicator")
  # Facets follow the column order of CvParam.txt itself (convergenceIndicator
  # first, as in the file); `category` is still available for zoom=.
  ordered_params <- param_cols
  
  out <- wide |>
    pivot_longer(all_of(param_cols), names_to = "parameter", values_to = "value") |>
    mutate(category  = factor(categorize(parameter), levels = category_levels),
           parameter = factor(parameter, levels = ordered_params))
  
  attr(out, "model_name") <- basename(model_dir)
  out
}
#' SAEM convergence plot for one model
#'
#' Reads ChartsData/Saem/CvParam.txt (read_convergence()) and plots every
#' model parameter's trajectory across SAEM iterations, one panel per
#' parameter (free y-scale, since fixed effects, covariate effects, random
#' effect SDs and residual error parameters live on very different scales),
#' with a dashed vertical line at each burn-in / phase change.
#' Panels follow the column order of CvParam.txt. convergenceIndicator gets
#' its own panel drawn in a distinct colour (red), since it is the overall
#' SAEM objective rather than a model parameter.
#'
#' @param model_dir   model folder
#' @param zoom        optional character vector to only show parameters in
#'                    these categories (see read_convergence(): "Fixed
#'                    effect", "Covariate effect", "Random effect (SD)",
#'                    "IOV (SD)", "Correlation", "Residual error",
#'                    "Convergence indicator"). NULL (default) = all.
#'                    (Named zoom rather than category to avoid ever
#'                    colliding with the data's own `category` column --
#'                    see the note in plot_iparam_distribution() about the
#'                    same class of bug with `parameter`.)
#' @param ncol        passed to facet_wrap(); NULL lets ggplot choose.
#'
#' @return a ggplot object
#' @export
plot_convergence <- function(model_dir, zoom = NULL, ncol = NULL) {
  
  data       <- read_convergence(model_dir)
  model_name <- attr(data, "model_name")
  
  if (!is.null(zoom)) {
    data <- data |> dplyr::filter(as.character(.data$category) %in% zoom)
    if (nrow(data) == 0) {
      stop("zoom = ", paste(zoom, collapse = ", "), " matched no parameters.")
    }
    data$parameter <- droplevels(data$parameter)
  }
  
  data <- data |> mutate(line_colour = if_else(parameter == "convergenceIndicator", "#B2182B", "#2166AC"))
  
  phase_changes <- data |>
    dplyr::distinct(iteration, phase) |>
    dplyr::arrange(iteration) |>
    dplyr::filter(phase != dplyr::lag(phase)) |>
    dplyr::pull(iteration)
  
  p <- ggplot(data, aes(x = iteration, y = value, colour = line_colour)) +
    geom_line(linewidth = 0.4) +
    scale_colour_identity() +
    facet_wrap(~ parameter, scales = "free_y", ncol = ncol) +
    labs(x = "SAEM iteration", y = NULL,
         title = sprintf("SAEM convergence \u2014 %s", model_name),
         caption = "Dashed line(s): SAEM phase change. Convergence indicator shown in red.") +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold", size = 8),
          plot.title = element_text(face = "bold", size = 12),
          plot.caption = element_text(size = 8, color = "grey40"))
  
  if (length(phase_changes) > 0) {
    p <- p + geom_vline(xintercept = phase_changes, linetype = "dashed", colour = "grey60")
  }
  p
}

  ## B) Individual parameter distributions (shrinkage, visual)        ----------
#' Read the distributions of the individual parameters
#'
#' Reads `ChartsData/DistributionOfTheIndividualParameters/pdf.txt` (empirical
#' density of each individual parameter, per split) and, when present,
#' `visualGuides.txt` (the population reference curve of each parameter).
#'
#' @param model_dir Monolix results folder of the model (ChartsData must have
#'   been exported).
#'
#' @return A list with `density` (wide table of `<parameter>_abscissa` and
#'   `<parameter>_pdf` columns), `theoretical` (the reference curves, `NULL`
#'   if the file is absent) and `model_name`.
#' @export
read_iparam_pdf <- function(model_dir) {
  
  f <- file.path(model_dir, "ChartsData", "DistributionOfTheIndividualParameters", "pdf.txt")
  if (!file.exists(f)) {
    stop("ChartsData/DistributionOfTheIndividualParameters/pdf.txt not found in: ", model_dir)
  }
  
  density <- read.table(f, sep = ",", comment.char = "", header = TRUE,
                        stringsAsFactors = FALSE) |> as_tibble()
  if ("split" %in% colnames(density)) density <- bind_cols(density, parse_split(density$split))
  
  guides_file <- file.path(model_dir, "ChartsData", "DistributionOfTheIndividualParameters", "visualGuides.txt")
  theoretical <- if (file.exists(guides_file)) {
    read.table(guides_file, sep = ",", comment.char = "", header = TRUE,
               stringsAsFactors = FALSE) |> as_tibble()
  } else NULL
  
  list(density = density, theoretical = theoretical, model_name = basename(model_dir))
}
#' Individual parameter distribution(s) ("shrinkage", visual) for one model
#'
#' Reads ChartsData/DistributionOfTheIndividualParameters/pdf.txt
#' (read_iparam_pdf()) and plots the empirical distribution of one or more
#' individual model parameters' estimates across subjects (e.g. clearance,
#' absorption rate), against Monolix's population reference curve when
#' available for that parameter -- a narrow empirical distribution relative
#' to that reference is the visual sign of shrinkage.
#'
#' Defaults to every parameter found in the file, faceted together with
#' any Y_TYPE / split-derived stratification that also varies (route, etc).
#' Pass parameter = to restrict to one or a few. The "parameter" facet
#' itself is never subject to facet_max_levels (it IS the point of the
#' plot when parameter = NULL) -- only the Y_TYPE / split-derived columns
#' go through the usual resolve_facets() cardinality logic (Shared
#' helpers).
#'
#' The numeric shrinkage percentage is annotated directly inside each
#' parameter's panel (top-right corner) when show_shrinkage = TRUE (the
#' default) -- a caption would get crowded with more than one parameter
#' shown at once. Rather than adding a separate reader for shrinkage.txt,
#' this reuses read_summary(..., full = TRUE) (Imported helpers), whose
#' $indiv_params already carries the shrinkage for every parameter;
#' shrinkage_method picks which of the two methods it reports.
#' shrinkage.txt's only value not already in there (shrinkage_condDist)
#' isn't pulled in here.
#'
#' @param model_dir         model folder
#' @param parameter         character vector of parameter name(s) to show,
#'                          as they appear in the pdf.txt column prefix
#'                          (e.g. "ka2", "Cl"). NULL (default) = every
#'                          parameter found in the file. A wrong name lists
#'                          the available ones in the error.
#' @param show_shrinkage    TRUE (default) annotates each panel with the
#'                          parameter's shrinkage (%), from summary.txt via
#'                          read_summary(); silently skipped where
#'                          unavailable (missing summary.txt, or no random
#'                          effect on that parameter).
#' @param shrinkage_method  which shrinkage to annotate: "Conditional
#'                          Distribution" (default) or "Conditional Mode".
#' @param filter,facet_by,facet_max_levels,value_labels,x_trans,y_trans,xlab,ylab,title      same conventions as plot_npde(type = "density"), applying to the Y_TYPE / split-derived columns
#'                          only (not to `parameter`, see above).
#'
#' @return a ggplot object
#' @export
plot_iparam_distribution <- function(model_dir, parameter = NULL, show_shrinkage = TRUE,
                                     shrinkage_method = c("Conditional Distribution", "Conditional Mode"),
                                     filter = NULL, facet_by = "auto", facet_max_levels = 6,
                                     value_labels = NULL,
                                     x_trans = NULL, y_trans = NULL,
                                     xlab = NULL, ylab = NULL, title = NULL) {
  
  shrinkage_method <- match.arg(shrinkage_method)
  filter_quo <- rlang::enquo(filter)
  model_name <- basename(model_dir)
  
  rd          <- read_iparam_pdf(model_dir)
  data        <- relabel_columns(rd$density, value_labels)
  theoretical <- rd$theoretical
  
  available <- unique(str_remove(colnames(data)[str_ends(colnames(data), "_abscissa")], "_abscissa$"))
  # Copied into a plain variable (not reused as `parameter` from here on) so
  # that dplyr verbs below can freely use a bare `parameter` column without
  # ever needing .env$ to disambiguate it from this argument -- the same
  # class of bug fixed via .env$ in plot_convergence()'s zoom filter, sidestepped
  # here by simply not keeping two things named "parameter" in scope at once.
  requested_params <- if (is.null(parameter)) available else parameter
  if (!all(requested_params %in% available)) {
    stop("parameter = ", paste(setdiff(requested_params, available), collapse = ", "),
         " not found. Available parameters: ", paste(available, collapse = ", "))
  }
  
  # Reshape from wide ("<param>_abscissa" / "<param>_pdf" per parameter,
  # side by side) to long (one row per parameter, generic abscissa/pdf
  # columns), keeping only the requested parameter(s) and every non-
  # parameter column (split/Y_TYPE/DV/split-derived) untouched.
  param_value_cols <- unlist(lapply(available, function(p) paste0(p, c("_abscissa", "_pdf"))))
  other_cols       <- setdiff(colnames(data), param_value_cols)
  keep_cols        <- c(other_cols, unlist(lapply(requested_params, function(p) paste0(p, c("_abscissa", "_pdf")))))
  
  data_long <- data |>
    dplyr::select(all_of(keep_cols)) |>
    pivot_longer(cols = ends_with(c("_abscissa", "_pdf")),
                 names_to = c("parameter", ".value"),
                 names_pattern = "^(.*)_(abscissa|pdf)$") |>
    mutate(parameter = factor(parameter, levels = requested_params))
  
  fixed_cols       <- c("abscissa", "pdf", "parameter", "split", "DV")
  candidate_facets <- setdiff(colnames(data_long), fixed_cols)
  
  r <- resolve_facets(data_long, filter_quo, candidate_facets,
                      facet_by = facet_by, facet_max_levels = facet_max_levels)
  data_long  <- r$data
  facet_cols <- r$facet_cols
  
  multi_param <- length(requested_params) > 1
  if (multi_param) facet_cols <- base::union("parameter", facet_cols)
  
  data_long <- data_long |>
    group_by(across(any_of(facet_cols))) |>
    arrange(abscissa) |>
    mutate(.width = 0.8 * (lead(abscissa) - abscissa)) |>
    fill(.width, .direction = "down") |>
    mutate(.xmin = abscissa - .width / 2, .xmax = abscissa + .width / 2) |>
    ungroup()
  
  p <- ggplot(data_long, aes(xmin = .xmin, xmax = .xmax, ymin = 0, ymax = pdf)) +
    geom_rect(fill = "grey40", alpha = 0.5)
  
  has_reference <- FALSE
  if (!is.null(theoretical)) {
    theo_keep <- intersect(unlist(lapply(requested_params, function(p) paste0(p, c("_abscissa", "_pdf")))),
                           colnames(theoretical))
    if (length(theo_keep) > 0) {
      has_reference <- TRUE
      theo_long <- theoretical |>
        dplyr::select(all_of(theo_keep)) |>
        pivot_longer(everything(), names_to = c("parameter", ".value"),
                     names_pattern = "^(.*)_(abscissa|pdf)$") |>
        mutate(parameter = factor(parameter, levels = requested_params))
      p <- p + geom_line(data = theo_long, aes(x = abscissa, y = pdf),
                         colour = "black", linewidth = 0.6, inherit.aes = FALSE)
    }
  }
  
  p <- add_facets(p, facet_cols, scales = if (multi_param) "free" else "fixed")
  if (!is.null(x_trans)) p <- p + scale_x_continuous(transform = x_trans)
  if (!is.null(y_trans)) p <- p + scale_y_continuous(transform = y_trans)
  
  # Shrinkage, annotated per-panel (top-right corner) rather than in the
  # caption -- with several parameters shown at once the caption would get
  # crowded. Monolix's summary.txt only ever reports ONE shrinkage value
  # per parameter -- there is no per-route / per-split breakdown to pull
  # from (shrinkage.txt doesn't have one either, it just repeats the same
  # per-method values as summary.txt plus a shrinkage_condDist column not
  # used here) -- so when route/split panels are also shown, the same
  # annotation is repeated across them (ggplot2 replicates a layer's data
  # across facet dimensions it doesn't itself carry) and is labelled
  # "(global)" as a reminder that it isn't stratum-specific.
  if (show_shrinkage) {
    summary_file <- file.path(model_dir, "summary.txt")
    summary_obj  <- if (file.exists(summary_file)) {
      tryCatch(read_summary(summary_file, full = TRUE), error = function(e) NULL)
    } else NULL
    
    if (!is.null(summary_obj) && !is.null(summary_obj$indiv_params)) {
      has_other_facets <- length(setdiff(facet_cols, "parameter")) > 0
      shrink <- summary_obj$indiv_params |>
        dplyr::filter(parameter %in% requested_params, method == shrinkage_method, !is.na(shrinkage)) |>
        mutate(label = sprintf(if (has_other_facets) "Shrinkage: %.0f%% (global)" else "Shrinkage: %.0f%%", shrinkage),
               parameter = factor(parameter, levels = requested_params))
      if (nrow(shrink) > 0) {
        p <- p + geom_text(data = shrink, aes(label = label), x = Inf, y = Inf,
                           hjust = 1.1, vjust = 1.5, size = 2.6, fontface = "italic",
                           colour = "grey30", inherit.aes = FALSE)
      }
    }
  }
  
  default_title <- if (!multi_param) {
    sprintf("%s distribution \u2014 %s", requested_params, model_name)
  } else {
    sprintf("Individual parameter distributions \u2014 %s", model_name)
  }
  default_title <- title_with_filter(default_title, filter_quo)
  
  p + labs(x = if (!is.null(xlab)) xlab else if (!multi_param) requested_params else "Value",
           y = if (is.null(ylab)) "Density" else ylab,
           title = if (is.null(title)) default_title else title,
           caption = if (has_reference) {
             "Bars: empirical distribution of individual estimates. Line: population reference."
           } else {
             "Bars: empirical distribution of individual estimates. (No population reference curve found.)"
           }) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 12),
          plot.caption = element_text(size = 8, color = "grey40"))
}


  ## C) Covariate assessment (eta/parameter vs covariate correlation) ----------
#' read_covariate_assessment(): reads both Tests/*Covariates.txt files for
#'
#' one model. Neither is DV/Y_TYPE- or split-stratified.
#'   - random_effects: Tests/correlationRandomEffectsCovariates.txt -- tests
#'     EVERY covariate available in the dataset against each random effect
#'     (eta), whether or not that covariate is currently in the model. This
#'     is the covariate SEARCH assessment (Monolix's "statistical tests" /
#'     covariate assessment tab).
#'   - individual_parameters: Tests/correlationIndividualParametersCovariates.txt
#'     -- tests only the covariate effects ALREADY in the model (named by
#'     their beta_* coefficient), i.e. a post-hoc check that each declared
#'     effect remains supported by the individual estimates. Its p.value
#'     column sometimes reads "<2.2e-16" rather than a plain number --
#'     parsed here into a numeric p.value plus a p_bound logical flag
#'     (TRUE = "p is below this value", FALSE = "this is the exact p.value"),
#'     so p.value is always numeric in the returned tibble.
#'
#' @param model_dir Monolix results folder of the model (the Tests folder must exist).
#' @param covariates Optional character vector: keep only these covariates.
#'
#' @export
read_covariate_assessment <- function(model_dir, covariates = NULL) {
  
  re_file   <- file.path(model_dir, "Tests", "correlationRandomEffectsCovariates.txt")
  ip_file   <- file.path(model_dir, "Tests", "correlationIndividualParametersCovariates.txt")
  wald_file <- file.path(model_dir, "Tests", "fixedEffects.txt")
  
  cov_pattern <- if (!is.null(covariates)) paste(covariates, collapse = "|") else NULL
  filter_covariates <- function(data, col) {
    if (is.null(data) || is.null(cov_pattern)) return(data)
    filtered <- dplyr::filter(data, str_detect(.data[[col]], regex(cov_pattern, ignore_case = TRUE)))
    if (nrow(filtered) == 0) {
      stop("covariates = ", paste(covariates, collapse = ", "),
           " matched no rows (column '", col, "'). Available values: ",
           paste(unique(data[[col]]), collapse = ", "), call. = FALSE)
    }
    filtered
  }
  
  random_effects <- if (file.exists(re_file)) {
    tryCatch(
      read.table(re_file, sep = ",", comment.char = "", header = TRUE, stringsAsFactors = FALSE) |> as_tibble(),
      error = function(e) stop("Failed reading correlationRandomEffectsCovariates.txt: ", conditionMessage(e), call. = FALSE)
    )
  } else NULL
  random_effects <- filter_covariates(random_effects, "covariate")
  
  individual_parameters <- if (file.exists(ip_file)) {
    tryCatch(
      read.table(ip_file, sep = ",", comment.char = "", header = TRUE, stringsAsFactors = FALSE) |>
        as_tibble() |>
        mutate(p_bound = str_starts(p.value, "<"),
               p.value = as.numeric(str_remove(p.value, "^<"))),
      error = function(e) stop("Failed reading correlationIndividualParametersCovariates.txt: ", conditionMessage(e), call. = FALSE)
    )
  } else NULL
  individual_parameters <- filter_covariates(individual_parameters, "covariate")
  
  wald <- if (file.exists(wald_file)) {
    tryCatch({
      wald_raw <- read.table(wald_file, sep = ",", comment.char = "", header = TRUE,
                             stringsAsFactors = FALSE) |>
        as_tibble() |>
        rename(statistics = statistics.sa., p.value = p.value.sa.)
      
      summary_file <- file.path(model_dir, "summary.txt")
      pop_params   <- if (file.exists(summary_file)) {
        tryCatch(read_summary(summary_file, full = TRUE)$pop_params, error = function(e) NULL)
      } else NULL
      
      if (!is.null(pop_params) && "parameter" %in% colnames(pop_params) && "estimate" %in% colnames(pop_params)) {
        wald_raw |>
          left_join(pop_params |> select(parameter, value = estimate), by = "parameter") |>
          mutate(value_available = !is.na(value)) |>
          relocate(parameter, value, .before = statistics)
      } else {
        wald_raw |> mutate(value = NA_real_, value_available = FALSE) |> relocate(value, .before = statistics)
      }
    }, error = function(e) stop("Failed reading/joining fixedEffects.txt: ", conditionMessage(e), call. = FALSE))
  } else NULL
  wald <- filter_covariates(wald, "parameter")
  
  list(random_effects = random_effects, individual_parameters = individual_parameters,
       wald = wald, model_name = basename(model_dir))
}
#' Random effect x covariate correlation heatmap (covariate search)
#'
#' Reads Tests/correlationRandomEffectsCovariates.txt
#' (read_covariate_assessment()) and plots a heatmap of every (random
#' effect, covariate) pair Monolix tested: one tile per pair, filled by
#' -log10(p-value) on a white -> green scale, with the raw p-value printed
#' on top. Rows = random effects (one per model parameter with a random
#' effect); columns = every covariate Monolix tested against them,
#' including transformed variants (e.g. "tSEX" alongside "SEX") if the
#' model uses any. A greener tile means a stronger association between
#' that random effect and that covariate, i.e. a stronger argument for
#' ADDING the covariate to that parameter.
#'
#' The colour scale is capped at -log10(p) = max_log10p (3 by default,
#' i.e. p = 0.001) rather than stretched to each model's own smallest
#' p-value, so heatmaps from different models/runs use the same scale and
#' stay directly comparable; anything beyond the cap is shown at full
#' colour (legend: ">= 3"). -log10() is used instead of the raw p because
#' p-values span many orders of magnitude, which a linear colour scale
#' would squash into one pale colour, and because it makes "bigger = more
#' significant" match "greener = stronger".
#'
#' Text colour is chosen per tile from the tile's actual rendered lightness
#' (relative luminance of the interpolated fill colour), so it stays
#' readable wherever a p-value sits on the scale; significant cells
#' (p < alpha) are additionally shown in bold.
#'
#' @param model_dir           model folder
#' @param covariates          optional character vector to restrict the
#'                            heatmap to matching covariates (case-
#'                            insensitive substring match, e.g. "WT" also
#'                            matches "logtWT"). NULL (default) keeps all.
#' @param alpha               significance threshold shown in bold (0.05
#'                            by default)
#' @param max_log10p          upper bound of the colour scale, in -log10(p)
#'                            units (3 = p of 0.001)
#' @param fill_colour         colour of the strongest tiles ("#2E8B57")
#' @param light_text_below    tile luminance (0-255, weights .299/.587/.114)
#'                            below which the text turns white; 140 by default,
#'                            i.e. about 77% of the way up the colour scale
#'                            (p of roughly 0.005 with the default cap)
#' @param eta_labels,covariate_labels  optional named character vectors to
#'                            relabel eta / covariate names (via
#'                            apply_labels(), Shared helpers)
#'
#' @return a ggplot object
#' @export
plot_covariate_search <- function(model_dir, covariates = NULL, alpha = 0.05,
                                  max_log10p = 3, fill_colour = "#2E8B57", light_text_below = 140,
                                  eta_labels = NULL, covariate_labels = NULL) {

  ca   <- read_covariate_assessment(model_dir, covariates = covariates)
  data <- ca$random_effects
  if (is.null(data)) {
    stop("Tests/correlationRandomEffectsCovariates.txt not found in: ", model_dir)
  }

  if (!is.null(eta_labels))       data$eta       <- apply_labels(data$eta, eta_labels)
  if (!is.null(covariate_labels)) data$covariate <- apply_labels(data$covariate, covariate_labels)

  ramp <- grDevices::colorRamp(c("white", fill_colour))
  data <- data |>
    mutate(
      neg_log10_p = -log10(pmax(p.value, 1e-300)),
      significant = p.value < alpha,
      .scaled     = pmin(pmax(neg_log10_p / max_log10p, 0), 1),
      .rgb        = ramp(.scaled),
      .luminance  = 0.299 * .rgb[, 1] + 0.587 * .rgb[, 2] + 0.114 * .rgb[, 3],
      text_colour = if_else(.luminance < light_text_below, "white", "grey10")
    )

  brks <- sort(unique(c(seq(0, floor(max_log10p), by = 1), max_log10p)))
  lbls <- as.character(brks)
  lbls[length(lbls)] <- paste0("\u2265 ", max_log10p)

  ggplot(data, aes(x = covariate, y = eta, fill = neg_log10_p)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = signif(p.value, 2), colour = text_colour,
                 fontface = if_else(significant, "bold", "plain")), size = 2.8) +
    scale_colour_identity() +
    scale_fill_gradient(low = "white", high = fill_colour, name = "-log10(p)",
                        limits = c(0, max_log10p), oob = scales::squish,
                        breaks = brks, labels = lbls) +
    labs(x = NULL, y = NULL,
        title = sprintf("Random effect \u00d7 covariate correlation \u2014 %s", ca$model_name),
        caption = sprintf("Bold: p < %s. Numbers: raw p-values. Colour scale capped at -log10(p) = %s so plots stay comparable across models.",
                          alpha, max_log10p)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
         legend.position = "bottom",
         legend.key.width = grid::unit(1.5, "cm"),
         legend.key.height = grid::unit(0.4, "cm"),
         plot.title = element_text(face = "bold", size = 12),
         plot.caption = element_text(size = 8, color = "grey40"))
}


  ## D) Nicely formatted covariate assessment tables                  ----------
#' Nicely formatted covariate assessment table (search, confirmation, or Wald)
#'
#' Builds one of the three "Tests > Covariates" tables from
#' read_covariate_assessment() as a styled kableExtra table, reusing
#' auto_pack_rows() / auto_style_cells() / clean_cov_labels() (Imported
#' helpers). Helper columns used only for styling (bold/colour flags,
#' p_bound) are kept in the data passed to auto_style_cells() (which needs
#' them to find column positions) but dropped from `display`, the data
#' actually passed to kable(), so they never show up as extra columns.
#'
#' The two kinds of table read in opposite directions, which is what the
#' emphasis follows:
#'   - "search" tests candidate covariates NOT (necessarily) in the model:
#'     a SMALL p-value means the covariate is associated with the random
#'     effect, i.e. a candidate to ADD. Rows with p < alpha are bold, with
#'     a green scale by p-value tier (matching plot_covariate_search()).
#'   - "confirmation" and "wald" test covariate effects ALREADY in the
#'     model: a SMALL p-value means the effect is well supported (keep it);
#'     a LARGE p-value (p > alpha_remove) means it is not, i.e. a candidate
#'     to REMOVE. Those rows are bold and red; well-supported rows are left
#'     plain.
#'
#' @param model_dir     model folder
#' @param which         "search" (default; Tests/correlationRandomEffects
#'                      Covariates.txt), "confirmation" (correlation test
#'                      for effects already in the model;
#'                      Tests/correlationIndividualParametersCovariates.txt)
#'                      or "wald" (Wald test, stochastic approximation;
#'                      Tests/fixedEffects.txt, with the point estimate
#'                      joined in from summary.txt when available)
#' @param covariates    optional character vector to restrict the table to
#'                      matching covariates (case-insensitive substring
#'                      match, e.g. "WT" also matches "logtWT" and
#'                      "beta_Cl_logtWT"). NULL (default) keeps every row.
#' @param alpha         "search" only: p-value below which a row is bold
#'                      (0.10 by default)
#' @param alpha_remove  "confirmation"/"wald" only: p-value above which an
#'                      existing covariate effect is flagged as a removal
#'                      candidate (0.05 by default)
#'
#' @return a kableExtra table
#' @export
covariate_tests_table <- function(model_dir, which = c("search", "confirmation", "wald"),
                                  covariates = NULL, alpha = 0.10, alpha_remove = 0.05) {

  which <- match.arg(which)
  ca    <- read_covariate_assessment(model_dir, covariates = covariates)

  remove_red <- "firebrick"

  if (which == "search") {

    data <- ca$random_effects
    if (is.null(data)) stop("Tests/correlationRandomEffectsCovariates.txt not found in: ", model_dir)

    data <- data |>
      mutate(
        covariate  = clean_cov_labels(covariate),
        statistics = round(statistics, 2),
        bold_col   = as.integer(p.value < alpha),
        color_col  = case_when(
          p.value <= 0.01                    ~ "#1B5E3A",
          p.value >  0.01 & p.value <= 0.05   ~ "#2E8B57",
          p.value >  0.05 & p.value <= 0.10   ~ "#6FAF8A",
          TRUE                                ~ "grey30"
        ),
        p.value    = format(p.value, scientific = TRUE, digits = 2)
      ) |>
      select(eta, covariate, statistics, p.value, bold_col, color_col)

    display <- data |> select(-bold_col, -color_col)

    display |>
      kable(booktabs = TRUE, col.names = c("Eta", "Covariate", "Statistics", "p-value"),
           caption = sprintf("Random effects vs covariates \u2014 %s", ca$model_name)) |>
      kableExtra::footnote(general = c("Pearson's correlation test (continuous) or ANOVA (categorical).",
                           "Test whether covariates should be added to the model.",
                           sprintf("Bold green: p < %s, the covariate is a candidate to add (darker = stronger).", alpha))) |>
      kable_styling(full_width = FALSE, font_size = 11, bootstrap_options = c("striped", "hover")) |>
      auto_pack_rows(data, group_vars = "eta") |>
      auto_style_cells(data, target_col = "p.value", bold_col = "bold_col", color_col = "color_col")

  } else if (which == "confirmation") {

    data <- ca$individual_parameters
    if (is.null(data)) stop("Tests/correlationIndividualParametersCovariates.txt not found in: ", model_dir)
    test_label <- unique(data$test)

    data <- data |>
      mutate(
        covariate  = clean_cov_labels(covariate),
        statistics = round(statistics, 2),
        bold_col   = as.integer(p.value > alpha_remove),   # computed before p.value becomes text
        color_col  = if_else(p.value > alpha_remove, remove_red, "grey30"),
        p.value    = if_else(p_bound, paste0("<", format(p.value, scientific = TRUE, digits = 2)),
                             format(p.value, scientific = TRUE, digits = 2))
      ) |>
      select(parameter, covariate, statistics, p.value, bold_col, color_col)

    display <- data |> select(-bold_col, -color_col)

    display |>
      kable(booktabs = TRUE, col.names = c("Parameter", "Covariate", "Statistics", "p-value"),
           caption = sprintf("Individual parameters vs covariates \u2014 %s", ca$model_name)) |>
      kableExtra::footnote(general = c(sprintf("%s.", paste(test_label, collapse = ", ")),
                           "Test whether covariates should be removed from the model.",
                           sprintf("Bold red: p > %s, the covariate effect is not supported (candidate to remove).", alpha_remove))) |>
      kable_styling(full_width = FALSE, font_size = 11, bootstrap_options = c("striped", "hover")) |>
      auto_pack_rows(data, group_vars = "parameter") |>
      auto_style_cells(data, target_col = "p.value", bold_col = "bold_col", color_col = "color_col")

  } else {

    data <- ca$wald
    if (is.null(data)) stop("Tests/fixedEffects.txt not found in: ", model_dir)

    data <- data |>
      mutate(
        covariate  = clean_cov_labels(parameter),
        value      = round(value, 3),
        statistics = round(statistics, 2),
        bold_col   = as.integer(p.value > alpha_remove),   # computed before p.value becomes text
        color_col  = if_else(p.value > alpha_remove, remove_red, "grey30"),
        p.value    = format(p.value, scientific = TRUE, digits = 2)
      ) |>
      select(covariate, value, statistics, p.value, bold_col, color_col, value_available)

    display <- data |> select(-bold_col, -color_col, -value_available)

    value_note <- if (any(!data$value_available)) {
      "Value: point estimate from summary.txt (blank where unavailable)."
    } else {
      "Value: point estimate from summary.txt."
    }

    display |>
      kable(booktabs = TRUE, col.names = c("Covariate effect", "Value", "Statistics", "p-value"),
           caption = sprintf("Individual parameters vs covariates \u2014 %s", ca$model_name)) |>
      kableExtra::footnote(general = c("Wald test (stochastic approximation).",
                           "Test whether covariates should be removed from the model.",
                           sprintf("Bold red: p > %s, the covariate effect is not supported (candidate to remove).", alpha_remove),
                           value_note)) |>
      kable_styling(full_width = FALSE, font_size = 11, bootstrap_options = c("striped", "hover")) |>
      auto_style_cells(data, target_col = "p.value", bold_col = "bold_col", color_col = "color_col")
  }
}
