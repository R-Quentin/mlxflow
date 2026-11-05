test_that("run_diagnostics names the report after the model and writes it where asked", {
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  fixture <- normalizePath(test_path("fixtures", "summary_warfarin_final.txt"))   # before changing directory
  dir <- withr::local_tempdir(); withr::local_dir(dir)
  dir.create("models/warfarin_final", recursive = TRUE)
  file.copy(fixture, "models/warfarin_final/summary.txt")
  writeLines(c("---", "title: stub", "output: html_document", "params:", "  model_dir: x", "  save_output: FALSE",
               "  output_dir: ''", "  prefix: ''", "---", "", "```{r}", "s <- read_summary(file.path(params$model_dir, 'summary.txt'))",
               "cat('covariates:', s$covariates)", "```"), "stub.Rmd")
  quiet <- function(...) suppressWarnings(suppressMessages(run_diagnostics(..., workflow = "stub.Rmd")))

  quiet("models/warfarin_final")
  expect_true(file.exists("warfarin_final_stub.html"))                        # model name, not a generic name
  quiet("models/warfarin_final", knit_dir = "reports", knit_prefix = "custom")
  expect_true(file.exists("reports/custom_stub.html"))
  quiet("models/warfarin_final", prefix = "PFX")
  expect_true(file.exists("PFX_stub.html"))
  expect_error(run_diagnostics("models/zzz", workflow = "stub.Rmd"), "model_dir not found")
  expect_error(run_diagnostics("models/warfarin_final", workflow = "no_such_template"), "template not found")
})

test_that("the packaged template exists and declares the expected parameters", {
  tpl <- system.file("rmarkdown", "templates", "model_checkup", "skeleton", "skeleton.Rmd", package = "mlxflow")
  expect_true(nzchar(tpl))
  yaml <- rmarkdown::yaml_front_matter(tpl)
  expect_true(all(c("model_dir", "save_output", "output_dir", "prefix") %in% names(yaml$params)))
})

test_that("run_comparison names the report after both models and saves outputs in its own folder", {
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  fixture <- normalizePath(test_path("fixtures", "summary_warfarin_final.txt"))
  dir <- withr::local_tempdir(); withr::local_dir(dir)
  for (m in c("m1", "m2")) {
    dir.create(file.path("models", m), recursive = TRUE)
    file.copy(fixture, file.path("models", m, "summary.txt"))
  }
  writeLines(c("---", "title: stub", "output: html_document", "params:", "  model_dir_1: x", "  model_dir_2: x",
               "  label_1: ''", "  label_2: ''", "  save_output: FALSE", "  output_dir: ''", "  prefix: ''", "---", "",
               "```{r}", "cat(params$label_1, params$label_2, basename(params$output_dir), params$prefix)", "```"), "cmp.Rmd")
  quiet <- function(...) suppressWarnings(suppressMessages(run_comparison("models/m1", "models/m2", workflow = "cmp.Rmd", ...)))

  quiet()
  expect_true(file.exists("m1_vs_m2_cmp.html"))
  quiet(label_1 = "Base", label_2 = "Final", knit_dir = "reports")
  expect_true(file.exists("reports/Base_vs_Final_cmp.html"))
  html <- paste(readLines("reports/Base_vs_Final_cmp.html"), collapse = " ")
  expect_match(html, "Base Final Base_vs_Final_comparison Base_vs_Final")      # labels, output folder, prefix
  expect_error(run_comparison("models/m1", "models/zzz", workflow = "cmp.Rmd"), "model_dir not found")
})

test_that("the packaged comparison template declares the expected parameters", {
  tpl <- system.file("rmarkdown", "templates", "model_comparison", "skeleton", "skeleton.Rmd", package = "mlxflow")
  expect_true(nzchar(tpl))
  yaml <- rmarkdown::yaml_front_matter(tpl)
  expect_true(all(c("model_dir_1", "model_dir_2", "label_1", "label_2", "save_output", "output_dir", "prefix") %in% names(yaml$params)))
})
