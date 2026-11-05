#' Run the model diagnostic report
#'
#' Knits an Rmarkdown report for one Monolix model. Two kinds of output,
#' controlled independently:
#' \itemize{
#'   \item the html **report**, always written: `<knit_prefix>_<template>.html`
#'     in `knit_dir`;
#'   \item the standalone **plots, tables and .rds**, written only if
#'     `save_output = TRUE`, in `output_dir`, named `<prefix>_<name>.<ext>`.
#' }
#'
#' @param model_dir Monolix results folder of the model (the one holding
#'   summary.txt; ChartsData must have been exported).
#' @param workflow The report template: the name of a template of the package
#'   (default `"model_checkup"`) or the path of your own .Rmd file. A template
#'   must declare the parameters `model_dir`, `save_output`, `output_dir` and `prefix`.
#' @param knit_dir Folder of the html report (default: the working directory).
#' @param knit_prefix Prefix of the report file name (default: `prefix` if
#'   given, otherwise the model folder name).
#' @param save_output If TRUE, also save plots and tables in `output_dir`.
#' @param output_dir Folder of the saved plots and tables (default: `model_dir`).
#' @param prefix Prefix of the saved plots and tables (default: the model folder name).
#'
#' @return The path of the html report, invisibly.
#' @export
run_diagnostics <- function(model_dir,
                            workflow    = "model_checkup",
                            knit_dir    = NULL,
                            knit_prefix = NULL,
                            save_output = FALSE,
                            output_dir  = "",
                            prefix      = "") {

  if (!dir.exists(model_dir)) stop("model_dir not found: ", model_dir, call. = FALSE)
  template <- .find_template(workflow)

  # Absolute paths: the report is knitted from another folder than the console
  model_dir  <- normalizePath(model_dir, winslash = "/")
  model_name <- basename(model_dir)
  if (!is.null(output_dir) && nzchar(output_dir)) {
    output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  }

  # Report name: custom prefix if given, otherwise the model name
  if (is.null(knit_prefix) || !nzchar(knit_prefix)) {
    knit_prefix <- if (nzchar(prefix)) prefix else model_name
  }

  .render_report(template,
                 output_file = paste0(knit_prefix, "_", template$name, ".html"),
                 knit_dir    = knit_dir,
                 params = list(model_dir   = model_dir,
                               save_output = save_output,
                               output_dir  = if (is.null(output_dir)) "" else output_dir,
                               prefix      = prefix))
}


#' Run the report comparing two models
#'
#' Knits an Rmarkdown report that puts two Monolix models side by side:
#' parameters and information criteria in common tables, then the predictive
#' performance (prediction-error metrics, observations vs predictions, VPC,
#' NPDE) as paired plots. Covariate effects and convergence are left to
#' [run_diagnostics()].
#'
#' The html report is always written; the standalone tables, plots and .rds
#' only if `save_output = TRUE`.
#'
#' @param model_dir_1,model_dir_2 Monolix results folders of the two models
#'   (the ones holding summary.txt). The plot-based sections need ChartsData
#'   for both models; a section is skipped, with a note, when it is missing.
#' @param label_1,label_2 Names of the two models in tables and plots
#'   (default: the folder names).
#' @param workflow The report template: the name of a template of the package
#'   (default `"model_comparison"`) or the path of your own .Rmd file. A
#'   template must declare the parameters `model_dir_1`, `model_dir_2`,
#'   `label_1`, `label_2`, `save_output`, `output_dir` and `prefix`.
#' @param knit_dir Folder of the html report (default: the working directory).
#' @param knit_prefix Prefix of the report file name (default: `prefix` if
#'   given, otherwise `<label_1>_vs_<label_2>`).
#' @param save_output If TRUE, also save tables and plots in `output_dir`.
#' @param output_dir Folder of the saved tables and plots (default:
#'   `<knit_dir>/<label_1>_vs_<label_2>_comparison`).
#' @param prefix Prefix of the saved file names (default: `<label_1>_vs_<label_2>`).
#'
#' @return The path of the html report, invisibly.
#' @export
run_comparison <- function(model_dir_1, model_dir_2,
                           label_1     = NULL,
                           label_2     = NULL,
                           workflow    = "model_comparison",
                           knit_dir    = NULL,
                           knit_prefix = NULL,
                           save_output = FALSE,
                           output_dir  = "",
                           prefix      = "") {

  for (d in c(model_dir_1, model_dir_2)) {
    if (!dir.exists(d)) stop("model_dir not found: ", d, call. = FALSE)
  }
  template <- .find_template(workflow)

  model_dir_1 <- normalizePath(model_dir_1, winslash = "/")
  model_dir_2 <- normalizePath(model_dir_2, winslash = "/")
  if (is.null(label_1)) label_1 <- basename(model_dir_1)
  if (is.null(label_2)) label_2 <- basename(model_dir_2)
  stem <- paste0(label_1, "_vs_", label_2)

  if (is.null(knit_prefix) || !nzchar(knit_prefix)) knit_prefix <- if (nzchar(prefix)) prefix else stem
  if (is.null(knit_dir) || !nzchar(knit_dir)) knit_dir <- getwd()
  if (is.null(output_dir) || !nzchar(output_dir)) output_dir <- file.path(knit_dir, paste0(stem, "_comparison"))
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  .render_report(template,
                 output_file = paste0(knit_prefix, "_", template$name, ".html"),
                 knit_dir    = knit_dir,
                 params = list(model_dir_1 = model_dir_1, model_dir_2 = model_dir_2,
                               label_1 = label_1, label_2 = label_2,
                               save_output = save_output, output_dir = output_dir,
                               prefix = if (nzchar(prefix)) prefix else stem))
}


# Template: an existing file is used as is, otherwise a template shipped with the package
#' @noRd
.find_template <- function(workflow) {
  name <- tools::file_path_sans_ext(basename(workflow))
  if (file.exists(workflow)) {
    return(list(path = normalizePath(workflow, winslash = "/"), name = name, packaged = FALSE))
  }
  path <- system.file("rmarkdown", "templates", name, "skeleton", "skeleton.Rmd", package = "mlxflow")
  if (!nzchar(path)) stop("Report template not found: ", workflow, call. = FALSE)
  list(path = path, name = name, packaged = TRUE)
}

# Knits the template with the given params. A packaged template lives in a read-only
# folder, so a temporary copy is knitted. The report is evaluated in an environment whose
# parent is the package namespace, so its chunks can call the internal functions too.
#' @noRd
.render_report <- function(template, output_file, knit_dir, params) {
  if (is.null(knit_dir) || !nzchar(knit_dir)) knit_dir <- getwd()
  dir.create(knit_dir, showWarnings = FALSE, recursive = TRUE)

  input <- template$path
  if (template$packaged) {
    work <- tempfile("mlxflow_report_"); dir.create(work)
    file.copy(list.files(dirname(input), full.names = TRUE), work, recursive = TRUE)
    input <- file.path(work, basename(input))
  }

  report <- rmarkdown::render(
    input         = input,
    output_file   = output_file,
    output_dir    = knit_dir,
    knit_root_dir = getwd(),                       # relative paths keep their console meaning
    envir         = new.env(parent = parent.env(environment())),
    params        = params
  )
  message("Report saved to: ", report)
  invisible(report)
}
