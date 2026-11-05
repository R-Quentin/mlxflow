# Shared helpers          ======================================================

# Generic building blocks reused across "Functions: Performances" (and
# usable by "Functions: Structure" too, should the need arise): relabeling
# a discrete column, parsing Monolix's "split" string, highlighting
# censored (BLQ) points, and the filter= / facet_by= scaffolding shared by
# plot_vpc() and plot_npde().

    ### A) apply_labels() -- recode + reorder a discrete column            -----
#' Recodes x using a named lookup vector (labels[old] = new), leaving any
#'
#' value absent from the lookup untouched -- safer than dplyr::recode()
#' here since a labels vector does not need to cover every level. When
#' labels is supplied, the result is also turned into a factor ordered by
#' the (unique) order of labels' values, so the order in which the vector
#' is written controls facet / table row order as well as the displayed
#' text (e.g. c("#OCC_Route:1" = "Oral", "#OCC_Route:2" = "IM") both
#' renames AND puts "Oral" before "IM"). Unmatched values, if any, are
#' appended after, in their order of appearance in the data. To reorder
#' without renaming, simply map values to themselves in the desired order
#' (e.g. y_type_labels = c(y1 = "y1", y2 = "y2")).
#'
#' @noRd
apply_labels <- function(x, labels) {
  if (is.null(labels)) return(x)
  recoded <- unname(labels[x])
  recoded <- if_else(is.na(recoded), x, recoded)
  ordered_levels <- unique(c(unname(labels), recoded))
  factor(recoded, levels = ordered_levels)
}
#' relabel_columns(): applies apply_labels() to every column named in
#'
#' value_labels (a named list of named character vectors, e.g.
#' list(OCC_Route = c("1" = "Oral", "2" = "IM"))), silently skipping any
#' name not present in `data`. Used by plot_vpc() / plot_npde() /
#' read_obspred() to relabel Y_TYPE and/or split-derived columns.
#'
#' @noRd
relabel_columns <- function(data, value_labels) {
  if (is.null(value_labels)) return(data)
  for (col in names(value_labels)) {
    if (col %in% colnames(data)) data[[col]] <- apply_labels(data[[col]], value_labels[[col]])
  }
  data
}

    ### B) parse_split() -- unpack Monolix's "split" stratification string -----
#' Splits a raw "split" column ("#Key1:Val1#Key2:Val2...", or "All" / NA
#'
#' when no stratification applies -- e.g. a Y_TYPE-only model like
#' Thoueille) into one column per key found ACROSS THE WHOLE VECTOR, so a
#' key present for only some rows still gets its own column (NA elsewhere)
#' instead of being silently dropped. Values are always returned as
#' character (Monolix split values are arbitrary strings, not guaranteed
#' numeric) -- R's `==` coerces automatically, so filter = OCC_Route == 2
#' works exactly like == "2". Returns a 0-column tibble (same row count)
#' when no "#key:value" pattern is found at all (e.g. split == "All").
#'
#' @noRd
parse_split <- function(split) {
  
  split   <- as.character(split)
  matches <- str_match_all(split, "#([^:#]+):([^#]+)")
  # str_match_all() returns a 0-row matrix for a non-matching string (e.g.
  # "All") but a 1-row all-NA matrix for an NA string -- drop that row so
  # an NA split value doesn't get counted as a spurious "NA"-named key.
  matches <- map(matches, function(m) m[!is.na(m[, 2]), , drop = FALSE])
  keys    <- unique(unlist(lapply(matches, function(m) m[, 2])))
  
  if (length(keys) == 0) return(tibble(.rows = length(split)))
  
  cols <- map(keys, function(k) {
    map_chr(matches, function(m) {
      hit <- m[m[, 2] == k, 3]
      if (length(hit) == 0) NA_character_ else hit[[1]]
    })
  })
  names(cols) <- keys
  as_tibble(cols)
}

    ### C) layer_censored_points() -- highlight BLQ observations           -----
#' A geom_point() layer highlighting BLQ (censored == 1) observations with
#'
#' a white fill, quantifiable ones (censored == 0) with `point_color`.
#' Expects `data` to already carry a `censored_label` factor column (levels
#' "Above"/"Below") -- callers compute it once before calling ggplot(),
#' since the layer itself only reads the plot's own data via aes(). Falls
#' back to plain colored points when highlight_censored = FALSE or no
#' censored_label column is available.
#'
#' @noRd
layer_censored_points <- function(data, highlight_censored = TRUE,
                                  point_color = "#2166AC", size = 1.6, alpha = 0.55) {
  if (!highlight_censored || !"censored_label" %in% colnames(data)) {
    return(geom_point(color = point_color, size = size, alpha = alpha))
  }
  list(
    geom_point(aes(fill = censored_label), shape = 21, colour = "grey30",
               size = size, alpha = alpha),
    scale_fill_manual(values = c(Above = point_color, Below = "white"),
                      breaks = c("Above", "Below"), name = "Quantification limit",
                      na.translate = FALSE)
  )
}

    ### D) filter= / facet_by= scaffolding (plot_vpc(), plot_npde())       -----
#' resolve_facets(): the filter= + facet_by= logic shared by plot_vpc() and
#'
#' plot_npde(). Applies `filter_quo` (a quosure captured by the caller via
#' rlang::enquo(filter)) to `data`, and STOPS with an informative error if
#' that leaves zero rows -- rather than silently building an empty,
#' confusing plot (this is what used to produce a cryptic "No shared levels
#' found..." warning from ggplot2's fill scale instead of a clear message:
#' it happened whenever filter= referenced a column's raw split code while
#' value_labels had already relabeled it, e.g. filter = OCC_Route == "IM"
#' with no value_labels supplied, when the raw values are still "1"/"2").
#'
#' It then decides which of `candidate_cols` (typically Y_TYPE plus
#' whichever columns parse_split() produced) to facet by:
#'   - facet_by = "auto" (default): every candidate that still varies after
#'     filtering, PROVIDED it has at most `facet_max_levels` distinct
#'     values. A candidate that varies but exceeds that cap is left out of
#'     the auto-selection and reported separately -- this is what a
#'     within-subject / repeated-measures dimension (e.g. an administration
#'     occasion index such as OCC_Admin, which can take a couple dozen
#'     values over a study) looks like, as opposed to a between-group
#'     stratification (e.g. an administration route, typically 2-3 values)
#'     that you *do* want side by side. There is no way to tell the two
#'     apart from the data alone, so the cardinality cap is the automatic
#'     default and facet_by= is the manual override for the cases it gets
#'     wrong either way (facet_max_levels= raises/lowers the cap itself if
#'     6 is not the right line for a given model).
#'   - facet_by = NULL: no faceting at all.
#'   - facet_by = a character vector: exactly those columns are used,
#'     regardless of facet_max_levels.
#' Any column that still varies and ends up not faceted -- whichever of the
#' three cases above -- triggers a warning, since strata would otherwise be
#' silently overplotted on the same panel.
#'
#' @return list(data = <filtered data>, facet_cols = <character vector>).
#'
#' @noRd
resolve_facets <- function(data, filter_quo, candidate_cols,
                           facet_by = "auto", facet_max_levels = 6) {
  
  if (!rlang::quo_is_null(filter_quo)) data <- dplyr::filter(data, !!filter_quo)
  if (nrow(data) == 0) {
    stop("filter = ", rlang::expr_text(rlang::quo_get_expr(filter_quo)),
         " matched no rows. If it references a column renamed via ",
         "value_labels, use the relabeled value in the filter -- ",
         "relabeling is applied before filtering.", call. = FALSE)
  }
  
  n_levels <- vapply(candidate_cols, function(col) n_distinct(data[[col]]), integer(1))
  varying  <- candidate_cols[n_levels > 1]
  high_card <- character(0)
  
  if (identical(facet_by, "auto")) {
    high_card <- varying[n_levels[varying] > facet_max_levels]
    chosen    <- setdiff(varying, high_card)
  } else if (is.null(facet_by)) {
    chosen <- character(0)
  } else {
    chosen <- facet_by
  }
  
  if (length(high_card) > 0) {
    warning("Column(s) not auto-faceted (more than facet_max_levels = ", facet_max_levels,
            " distinct values -- looks like a repeated-measures / occasion dimension ",
            "rather than a group to display side by side; pass facet_by = to include ",
            "one explicitly, filter= to restrict to one value, or raise facet_max_levels ",
            "if it should have been included): ", paste(high_card, collapse = ", "), call. = FALSE)
  }
  
  not_faceted <- setdiff(varying, base::union(chosen, high_card))
  if (length(not_faceted) > 0) {
    warning("Column(s) still vary but are not faceted -- strata will be overplotted on the same panel: ",
            paste(not_faceted, collapse = ", "), call. = FALSE)
  }
  
  list(data = data, facet_cols = chosen)
}
#' add_facets(): adds facet_wrap() over facet_cols (a character vector of
#'
#' column names, as returned by resolve_facets()) to a ggplot, or returns p
#' unchanged when facet_cols is empty.
#'
#' @noRd
add_facets <- function(p, facet_cols, scales = "fixed") {
  if (length(facet_cols) == 0) return(p)
  p + facet_wrap(as.formula(paste("~", paste(facet_cols, collapse = " + "))),
                 labeller = label_both, scales = scales)
}
#' title_with_filter(): appends the filter= expression (if any) to a
#'
#' default plot title, e.g. "NPDE vs time -- Han (OCC_Route == \\"IM\\")".
#'
#' @noRd
title_with_filter <- function(base, filter_quo) {
  if (rlang::quo_is_null(filter_quo)) return(base)
  paste0(base, " (", rlang::expr_text(rlang::quo_get_expr(filter_quo)), ")")
}

    ### E) ChartsData guard and DV resolution                              -----
#' Check that the ChartsData of a model have been exported
#'
#' Every plot-based diagnostic (convergence, distributions, observations vs
#' predictions, VPC, NPDE, covariate tests) reads from the `ChartsData` folder
#' of the model. This stops with an explanatory message when it is missing or
#' empty, instead of letting each diagnostic fail on its own.
#'
#' @param model_dir Monolix results folder of the model.
#'
#' @return `TRUE`, invisibly; an error otherwise.
#' @export
check_chartsdata <- function(model_dir) {
  charts_dir <- file.path(model_dir, "ChartsData")
  n_files    <- if (dir.exists(charts_dir)) length(list.files(charts_dir, recursive = TRUE)) else 0L
  if (n_files == 0) {
    stop("No ChartsData found in '", charts_dir, "' (the folder is missing or empty), ",
        "so the diagnostics cannot be built. Load the project with lixoftConnectors ",
        "(loadProject()) and run computeChartsData() first, then re-run this report.",
        call. = FALSE)
  }
  invisible(TRUE)
}
#' resolve_dv_layout(): decides the observed-variable name and the number of
#'
#' error-model levels used to read one ChartsData export family, i.e. files
#' named "<dv>_<suffix>.txt" (single error model) or "y1_<suffix>.txt",
#' "y2_<suffix>.txt", ... (one per level).
#'
#' The .mlxtran (or an explicit dv=) gives a first guess, but the files on
#' disk are the ground truth: Monolix names exports after the observation
#' NAME, which is not always what the .mlxtran parser (or the caller) ends
#' up with -- e.g. the dataset column may be called "dv" while the exports
#' are "concentration_*.txt". When the guess matches no file, the DV is
#' inferred from the exported file names (if there is exactly one), with a
#' message saying so; otherwise the error lists what is available.
#'
#' @noRd
resolve_dv_layout <- function(model_dir, charts_dir, suffix, dv = NULL) {

  file_pattern  <- sprintf("_%s\\.txt$", suffix)
  on_disk       <- sub(file_pattern, "", list.files(charts_dir, pattern = file_pattern))
  y_on_disk     <- on_disk[grepl("^y[0-9]+$", on_disk)]
  named_on_disk <- setdiff(on_disk, y_on_disk)

  dv_name  <- dv
  n_levels <- NULL
  source   <- "dv="

  if (is.null(dv)) {
    source   <- ".mlxtran"
    dv_table <- tryCatch(detect_dv_from_mlxtran(resolve_mlxtran_path(model_dir)),
                         error = function(e) NULL)
    if (!is.null(dv_table)) {
      if (nrow(dv_table) > 1) {
        warning(sprintf(
          "Several observed variables declared in .mlxtran (%s); using the first one ('%s'). Pass dv= explicitly to choose another.",
          paste(dv_table$dv_name, collapse = ", "), dv_table$dv_name[1]
        ), call. = FALSE)
      }
      dv_name  <- dv_table$dv_name[1]
      n_levels <- dv_table$n_levels[1]
    }
  }

  # (a) several error-model levels declared and y<k> files present
  if (!is.null(n_levels) && n_levels > 1 && length(y_on_disk) > 0) {
    return(list(dv_name = dv_name, n_levels = n_levels))
  }
  # (b) the guess matches a named export
  if (!is.null(dv_name) && dv_name %in% named_on_disk) {
    return(list(dv_name = dv_name, n_levels = 1))
  }
  # (c) y<k> exports only (explicit dv=: as many levels as files; .mlxtran: as declared)
  if (length(y_on_disk) > 0 && length(named_on_disk) == 0) {
    return(list(dv_name = if (is.null(dv_name)) "y" else dv_name,
                n_levels = if (is.null(dv) && !is.null(n_levels)) n_levels else length(y_on_disk)))
  }
  # (d) infer from the exported file names
  if (length(named_on_disk) == 1) {
    message("Observed-variable name ",
           if (is.null(dv_name)) "could not be read from the .mlxtran" else sprintf("'%s' (from %s) has no matching '%s_%s.txt'", dv_name, source, dv_name, suffix),
           " in ", charts_dir, " -- using '", named_on_disk, "' inferred from the exported file names.")
    return(list(dv_name = named_on_disk, n_levels = 1))
  }

  stop(if (length(named_on_disk) == 0) {
    sprintf("No '*_%s.txt' file found in %s", suffix, charts_dir)
  } else {
    sprintf("Could not match DV '%s' to the exported files in %s (found: %s). Pass dv = one of these.",
            if (is.null(dv_name)) "?" else dv_name, charts_dir, paste(named_on_disk, collapse = ", "))
  }, call. = FALSE)
}
#' Name of the observed variable of a model
#'
#' Found from the project (.mlxtran) and validated against the names of the
#' exported `ChartsData` files, which are the ground truth: a dataset column
#' called "dv" whose observation is named, say, "concentration" is handled.
#' Call it once per report and pass the result as `dv` to the diagnostic
#' functions.
#'
#' @param model_dir Monolix results folder of the model.
#'
#' @return The observed-variable name (character string).
#' @export
resolve_dv <- function(model_dir) {
  charts_dir <- file.path(model_dir, "ChartsData", "ObservationsVsPredictions")
  if (!dir.exists(charts_dir)) stop("ChartsData/ObservationsVsPredictions not found in: ", model_dir)
  resolve_dv_layout(model_dir, charts_dir, "obsVsPred", dv = NULL)$dv_name
}
