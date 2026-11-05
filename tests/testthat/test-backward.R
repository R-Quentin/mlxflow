# The Monolix side is replaced by a fake engine (helper-fake-engine.R): deterministic effects
# V_logtWt (25.68, 1 df), Cl_logtWt (4.93, 1 df), V_SEX (5.20, 2 df), Cl_SEX (12.00, 2 df)

run_backward <- function(out, ..., base = new_fake_project(), env = parent.frame()) {
  local_mocked_bindings(.lx = fake_lx, .env = env)
  suppressMessages(backward_elimination(base, alpha = 0.01, output_dir = out, ...))
}
history_ids <- function(res) paste(res$history$parameter, res$history$covariate, sep = "_")

test_that("exhaustive backward adopts the highest p-value and numbers runs strictly", {
  out <- withr::local_tempdir()
  res <- run_backward(out, fast_screening = FALSE)
  # iteration 1: V_SEX (5.20 on 2 df, p = 0.074) beats Cl_logtWt (4.93 on 1 df, p = 0.026)
  expect_equal(history_ids(res)[2:3], c("V_SEX", "Cl_logtWt"))
  folders <- list.files(out, pattern = "^run[0-9]{3}_")
  runs <- as.integer(substr(folders[!grepl("\\.mlxtran$", folders)], 4, 6))
  expect_equal(sort(runs), 0:9)                              # no duplicate, no gap
  expect_true(file.exists(file.path(out, "run000_I0_base.mlxtran")))
  log <- utils::read.table(file.path(out, "backward_log.csv"), sep = ";", header = TRUE)
  expect_equal(log$run, 0:9)
  expect_true(all(c("critical_ofv", "alpha", "ofv", "bicc") %in% names(log)))
  expect_match(res$stop_reason, "still significant")
})

test_that("fast screening reaches the same final model with fewer fits", {
  full <- run_backward(withr::local_tempdir(), fast_screening = FALSE)
  fast <- run_backward(withr::local_tempdir(), fast_screening = TRUE)
  expect_equal(history_ids(fast), history_ids(full))
  expect_lt(sum(fast$history$n_candidates_tested, na.rm = TRUE), sum(full$history$n_candidates_tested, na.rm = TRUE))
})

test_that("covariates are read correctly from folders whose summary has a Fisher section", {
  out <- withr::local_tempdir()
  run_backward(out, fast_screening = TRUE)                    # fast mode runs the SE step on parents
  models <- suppressMessages(read_models(out))
  true_n <- vapply(models$model_name, function(n)
    sum(grepl("^beta_.*:", readLines(file.path(out, n, "summary.txt")))), numeric(1))
  expect_equal(models$n_cov, unname(true_n))
  expect_true(any(vapply(models$model_name, function(n)
    any(grepl("FISHER", readLines(file.path(out, n, "summary.txt")))), logical(1))))
})

test_that("the flowchart builder retraces the path of the engine", {
  out <- withr::local_tempdir()
  res <- run_backward(out, fast_screening = FALSE)
  models <- suppressMessages(read_models(out)) |>
    dplyr::mutate(covariates = purrr::map(covariates, clean_cov_labels))
  b <- build_backward_data(models, alpha = 0.01)
  expect_equal(b$last_best_name, sub("\\.mlxtran$", "", res$final_project |> basename()))
})

test_that("an interrupted run resumes without duplicates and gives the same result", {
  ref <- run_backward(withr::local_tempdir(), fast_screening = FALSE)
  out <- withr::local_tempdir(); base <- new_fake_project()
  fk$eng$fail_at_run <- 7                                     # crash inside iteration 2
  expect_error(run_backward(out, fast_screening = FALSE, base = base), "simulated crash")
  state <- readRDS(file.path(out, "backward_state.rds"))
  expect_false(state$finished); expect_equal(state$iter, 1)
  fk$eng$fail_at_run <- Inf
  res <- run_backward(out, fast_screening = FALSE, base = base, resume = TRUE)
  log <- utils::read.table(file.path(out, "backward_log.csv"), sep = ";", header = TRUE)
  expect_false(anyDuplicated(log$run) > 0)
  expect_equal(res$history[, c("iter", "run", "parameter", "covariate", "ofv")],
               ref$history[, c("iter", "run", "parameter", "covariate", "ofv")])
})

test_that("existing outputs are protected, and overwrite only removes what the run created", {
  out <- withr::local_tempdir(); base <- new_fake_project()
  first <- run_backward(out, base = base)
  expect_error(run_backward(out, base = base), "already contains")                 # no silent mixing
  again <- run_backward(out, base = base, resume = TRUE)                           # finished: stored result
  expect_equal(again$stop_reason, first$stop_reason)
  # legacy folders and an unrelated file
  dir.create(file.path(out, "step01_V_logtWt")); writeLines("keep me", file.path(out, "my_notes.txt"))
  run_backward(out, base = base, overwrite = TRUE)
  expect_true(file.exists(file.path(out, "my_notes.txt")))
  expect_false(dir.exists(file.path(out, "step01_V_logtWt")))
})

test_that("save_all = FALSE keeps only the adopted path on disk, the log keeps everything", {
  out <- withr::local_tempdir()
  run_backward(out, fast_screening = FALSE, save_all = FALSE)
  kept <- list.files(out, pattern = "^run[0-9]{3}_.*\\.mlxtran$")
  log  <- utils::read.table(file.path(out, "backward_log.csv"), sep = ";", header = TRUE)
  expect_lt(length(kept), nrow(log))
  expect_equal(nrow(log), 10)
})

test_that("a failing diagnostic callback never stops the run", {
  out <- withr::local_tempdir()
  warns <- testthat::capture_warnings(
    res <- run_backward(out, diagnostic_fun = function(iter, project) stop("typo in my plot")))
  expect_true(any(grepl("diagnostic_fun failed", warns)))
  expect_true(nrow(res$history) > 1)
})
