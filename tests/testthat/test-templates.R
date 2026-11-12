# End-to-end test of the two report templates. It needs the Monolix output of a real model: put the
# folder of a model with its ChartsData exported (e.g. the warfarin demo) in
# tests/testthat/fixtures/warfarin_final/ (summary.txt + ChartsData/). The test is skipped otherwise.

model_fixture <- function() normalizePath(test_path("fixtures", "warfarin_final"), mustWork = FALSE)

skip_if_no_model <- function() {
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  skip_if(utils::packageVersion("ggplot2") < "3.5.0", "ggplot2 >= 3.5.0 needed")
  skip_if_not(dir.exists(file.path(model_fixture(), "ChartsData")), "no model fixture with ChartsData")
}

visible_text <- function(html_file) {
  html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  gsub("\\s+", " ", gsub("<[^>]+>", " ", html))
}

test_that("the check-up report renders, leaves out the optional sections and prints no NULL", {
  skip_if_no_model()
  model <- model_fixture()                       # resolved before changing directory
  withr::local_dir(withr::local_tempdir())
  out  <- suppressWarnings(suppressMessages(run_diagnostics(model, knit_dir = ".")))
  text <- visible_text(out)
  expect_match(text, "Run configuration")
  expect_no_match(text, "## NULL", fixed = TRUE)
  expect_no_match(text, "Convergence assessment provided")
})

test_that("the check-up report adds the convergence assessment section when it is provided", {
  skip_if_no_model()
  model <- model_fixture()
  withr::local_dir(withr::local_tempdir())
  one_run <- function(shift) list(
    populationParameters = list(convergence = data.frame(
      convergenceIndicator = 1:20, ka_pop = 1.4 + shift + 0.2 * exp(-(1:20) / 8), V_pop = 8 + shift + 0 * (1:20))),
    standardErrors = list(values = c(0.4, 0.3)))
  out  <- suppressWarnings(suppressMessages(run_diagnostics(
    model, knit_dir = ".", assessment = list(one_run(0), one_run(0.02), one_run(-0.03)))))
  expect_match(visible_text(out), "Convergence assessment provided")
})

test_that("the comparison report renders with its tables and paired plots", {
  skip_if_no_model()
  model <- model_fixture()
  withr::local_dir(withr::local_tempdir())
  out <- suppressWarnings(suppressMessages(run_comparison(
    model, model, label_1 = "A", label_2 = "B", knit_dir = ".", save_output = TRUE)))
  text <- visible_text(out)
  expect_match(text, "Information criteria")
  expect_match(text, "Prediction-error metrics")
  expect_no_match(text, "## NULL", fixed = TRUE)
  expect_true(file.exists("A_vs_B_comparison/A_vs_B_metrics.html"))
  expect_true(file.exists("A_vs_B_comparison/A_vs_B_vpc.png"))
})
