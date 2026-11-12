fixture <- function(x) testthat::test_path("fixtures", x)

# Builds a folder of models from the real warfarin summaries, in flat or nested layout
make_models_dir <- function(names, files, nested = FALSE, env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  for (i in seq_along(names)) {
    d <- if (nested) file.path(root, names[i], names[i]) else file.path(root, names[i])
    dir.create(d, recursive = TRUE)
    file.copy(fixture(files[i]), file.path(d, "summary.txt"))
  }
  root
}
files3 <- c("summary_warfarin_example.txt", "summary_ModelBuilding_007.txt", "summary_warfarin_final.txt")

test_that("read_models handles nested (ModelBuilding) and flat layouts, ids from folder numbers", {
  nested <- make_models_dir(sprintf("ModelBuilding_%03d", 1:3), files3, nested = TRUE)
  m <- suppressWarnings(suppressMessages(read_models(nested)))
  expect_equal(m$model_id, 1:3)
  expect_equal(m$model_name, sprintf("ModelBuilding_%03d", 1:3))
  expect_equal(m$n_cov, c(0, 1, 1))

  flat <- make_models_dir(sprintf("run%03d_I%d_x", 0:2, 0:2), files3)
  expect_equal(suppressWarnings(suppressMessages(read_models(flat)))$model_id, 0:2)
})

test_that("read_models falls back on alphabetical ids, with a warning", {
  dir <- make_models_dir(c("full", "no_cov", "final"), files3)
  warns <- testthat::capture_warnings(m <- suppressMessages(read_models(dir)))
  expect_true(any(grepl("alphabetical", warns)))
  expect_equal(m$model_name, c("final", "full", "no_cov"))
  expect_equal(m$model_id, 0:2)
})

test_that("read_models and compare_models fail clearly on an empty folder", {
  empty <- withr::local_tempdir()
  expect_error(read_models(empty), "No summary.txt found")
  expect_error(compare_models(empty), "No summary.txt found")
})

test_that("compare_models works on nested folders and on models without standard errors", {
  nested <- make_models_dir(sprintf("ModelBuilding_%03d", 1:3), files3, nested = TRUE)
  one  <- compare_models(nested, values = "estimate")
  expect_s3_class(one$kable, "knitr_kable")
  expect_true(all(sprintf("ModelBuilding_%03d", 1:3) %in% names(one$data)))
  multi <- compare_models(nested, values = c("estimate", "rse_pct"))
  expect_equal(ncol(multi$data), 2 + 3 * 2)
})

test_that("compare_two_models returns the data and a table, with optional highlighting", {
  r <- compare_two_models(fixture("summary_warfarin_example.txt"), fixture("summary_warfarin_final.txt"))
  expect_named(r, c("data", "kable"))
  expect_s3_class(r$kable, "knitr_kable")
  expect_true(all(c("estimate_1", "estimate_2", "delta_estimate", "delta_rel_estimate") %in% names(r$data)))
  r2 <- compare_two_models(fixture("summary_warfarin_example.txt"), fixture("summary_ModelBuilding_007.txt"),
                           values = c("estimate", "rse_pct"))
  expect_true(all(is.na(r2$data$rse_pct_2)))   # no standard errors in the second model

  # RSE above 30 percent in bold: Tlag_pop (22.9), ka_pop (28.6) stay plain, the 34.4 of
  # warfarin_final's ka_pop is bold
  hi <- compare_two_models(fixture("summary_warfarin_example.txt"), fixture("summary_warfarin_final.txt"),
                           values = "rse_pct", bold_above = c(rse_pct = 30))
  html <- as.character(hi$kable)
  expect_match(html, 'font-weight: bold;"> 34.4 </td>', fixed = TRUE)
  expect_no_match(html, 'font-weight: bold;"> 22.9 </td>', fixed = TRUE)
})

test_that("read_models warns when the log-likelihood methods differ", {
  # warfarin_example and warfarin_final: importance sampling; ModelBuilding_007: linearization
  mixed <- make_models_dir(sprintf("ModelBuilding_%03d", 1:3), files3, nested = TRUE)
  expect_warning(m <- suppressMessages(read_models(mixed)), "not computed by the same method")
  expect_equal(m$ll_method, c("is", "lin", "is"))
  same <- make_models_dir(sprintf("ModelBuilding_%03d", 1:2), files3[c(1, 3)], nested = TRUE)
  expect_no_warning(suppressMessages(read_models(same)))
})
