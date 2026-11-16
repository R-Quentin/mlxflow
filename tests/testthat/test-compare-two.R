fixture <- function(x) testthat::test_path("fixtures", x)

test_that("compare_criteria gives the criteria side by side with the difference", {
  cc <- compare_criteria(fixture("summary_warfarin_example.txt"), fixture("summary_warfarin_final.txt"),
                         labels = c("Base", "Final"))
  expect_named(cc, c("criterion", "Base", "Final", "delta"))
  expect_equal(cc$delta[cc$criterion == "OFV"], -25.75)
  expect_equal(unname(attr(cc, "ll_method")), c("is", "is"))
})

test_that("compare_criteria warns when the likelihood methods differ", {
  expect_warning(compare_criteria(fixture("summary_warfarin_final.txt"), fixture("summary_ModelBuilding_007.txt")),
                 "different methods")
})

test_that("compare_obspred_metrics keeps the order of the first model and computes the change", {
  m1 <- tibble::tibble(DV = "conc", Y_TYPE = "conc", Prediction = c("Population predictions", "Individual predictions"),
                       N = 10, ipred_type = c(NA, "mean"), MPE = c(12, -1), RMSE = c(1.9, 0.5))
  m2 <- dplyr::mutate(m1, MPE = MPE + 1, RMSE = RMSE - 0.1)
  r <- compare_obspred_metrics(m1, m2, labels = c("A", "B"), format = "long")
  expect_equal(unique(r$Prediction), c("Population predictions", "Individual predictions"))   # not alphabetical
  expect_equal(unique(r$metric), c("N", "MPE", "RMSE"))
  expect_equal(r$delta[r$metric == "MPE"], c(1, 1))
  expect_equal(r$delta[r$metric == "RMSE"], c(-0.1, -0.1))
  # accepts the list returned by obspred_diagnostics()
  expect_equal(compare_obspred_metrics(list(metrics = m1), list(metrics = m2), c("A", "B"), format = "long"), r)
})

test_that("compare_obspred_metrics gives metrics in columns by default: models then delta, per block", {
  m1 <- tibble::tibble(DV = "conc", Y_TYPE = "conc", Prediction = c("Population predictions", "Individual predictions"),
                       N = 10, MPE = c(12, -1), RMSE = c(1.9, 0.5))
  m2 <- dplyr::mutate(m1, MPE = MPE + 1)
  r <- compare_obspred_metrics(m1, m2, labels = c("A", "B"))
  expect_equal(names(r), c("DV", "Y_TYPE", "Prediction", "Model", "N", "MPE", "RMSE"))
  expect_equal(r$Model, rep(c("A", "B", "delta"), 2))
  expect_equal(r$Prediction, rep(c("Population predictions", "Individual predictions"), each = 3))
  expect_equal(r$MPE[1:3], c(12, 13, 1))
})

test_that("the strata of a split are matched, not mistaken for a metric", {
  strata <- function(n) tibble::tibble(DV = "conc", Y_TYPE = "conc", split = c("assay 1", "assay 2"),
                                       Prediction = "Population predictions", N = n, MPE = c(40, 70), RMSE = 2)
  r <- compare_obspred_metrics(strata(c(296, 1893)), strata(c(296, 1893)) |> dplyr::mutate(MPE = c(45, 60)), c("A", "B"))
  expect_equal(nrow(r), 6)                                              # 2 strata x (A, B, delta)
  expect_false("split" %in% r$Model); expect_equal(r$split, rep(c("assay 1", "assay 2"), each = 3))
  expect_equal(r$MPE[r$Model == "delta"], c(5, -10))
  long <- compare_obspred_metrics(strata(c(296, 1893)), strata(c(296, 1893)), c("A", "B"), format = "long")
  expect_false("split" %in% long$metric)
  expect_false(anyDuplicated(long[c("split", "Prediction", "metric")]) > 0)
})

test_that("strata in one model only cannot be paired with a metrics table, and unmatched rows are reported", {
  with_split <- tibble::tibble(DV = "conc", Y_TYPE = "conc", split = c("a", "b"), Prediction = "P", N = c(3, 7), MPE = 1)
  pooled     <- tibble::tibble(DV = "conc", Y_TYPE = "conc", Prediction = "P", N = 10, MPE = 1)
  expect_error(compare_obspred_metrics(with_split, pooled), "Only one of the two models has strata")
  other <- dplyr::mutate(with_split, split = c("a", "z"))
  expect_warning(compare_obspred_metrics(with_split, other), "one model only.*b")
})

test_that("a single observed variable with different names in the two models is matched", {
  m1 <- tibble::tibble(DV = "CAB_conc", Y_TYPE = "CAB_conc", Prediction = "P", N = 10, MPE = 1)
  m2 <- dplyr::mutate(m1, DV = "y1", Y_TYPE = "y1", MPE = 3)
  expect_message(r <- compare_obspred_metrics(m1, m2, c("A", "B")), "matched as the same observed variable")
  expect_equal(r$MPE, c(1, 3, 2)); expect_equal(unique(r$Y_TYPE), "CAB_conc")
})

test_that("delta chooses what the comparison tables show", {
  a <- fixture("summary_warfarin_example.txt"); b <- fixture("summary_warfarin_final.txt")
  n_head <- function(k) lengths(regmatches(as.character(k$kable), gregexpr("<th[ >]", as.character(k$kable))))
  expect_equal(n_head(compare_two_models(a, b, delta = FALSE)), 3)           # parameter + the two models
  expect_equal(n_head(compare_two_models(a, b, delta = "abs")), 4)
  expect_equal(n_head(compare_two_models(a, b, delta = TRUE)), n_head(compare_two_models(a, b)))
  expect_true(all(c("delta_estimate", "delta_rel_estimate") %in% names(compare_two_models(a, b, delta = FALSE)$data)))
  expect_error(compare_two_models(a, b, delta = "x"), "`delta` must be")

  root <- withr::local_tempdir()
  for (i in 1:3) { d <- file.path(root, c("m1", "m2", "m3")[i]); dir.create(d)
    file.copy(c(a, b, b)[i], file.path(d, "summary.txt")) }
  expect_equal(names(compare_models(root)$data), c("type", "parameter", "m1", "m2", "m3"))    # unchanged by default
  r <- compare_models(root, delta = TRUE)
  expect_equal(names(r$data), c("type", "parameter", "m1", "m2", "delta_m2", "delta_rel_m2", "m3", "delta_m3", "delta_rel_m3"))
  expect_equal(r$data$delta_m2[r$data$parameter == "Tlag_pop"], 0.05)
  r2 <- compare_models(root, values = c("estimate", "rse_pct"), delta = "abs", reference = "m2")
  expect_false(any(grepl("^delta_m2", names(r2$data)))); expect_true("delta_m1_rse_pct" %in% names(r2$data))
  expect_error(compare_models(root, reference = "zz"), "`reference` must be")
})

test_that("the parameter comparison table shows estimate, RSE and shrinkage, highlighted by thresholds", {
  cmp <- compare_two_models(fixture("summary_warfarin_example.txt"), fixture("summary_warfarin_final.txt"),
                            values = c("estimate", "rse_pct", "shrinkage"))$data
  k <- kable_param_comparison(cmp, labels = c("Base", "Final"))
  html <- as.character(k)
  expect_match(html, "Shrinkage \\(%\\)\\*"); expect_match(html, "RSE \\(%\\)")
  expect_match(html, "Change = Final - Base")
  expect_match(html, "font-weight: bold;color: rgba\\(205, 155, 29, 1\\) !important;\"> 34.4", perl = TRUE)   # RSE of ka_pop above 30
  expect_error(kable_param_comparison(cmp[c("type", "parameter")]), "none of the compared values")
})

test_that("compare_models can bold the values above a threshold", {
  root <- withr::local_tempdir()
  for (i in 1:2) {
    d <- file.path(root, c("a", "b")[i]); dir.create(d)
    file.copy(fixture(c("summary_warfarin_example.txt", "summary_warfarin_final.txt")[i]), file.path(d, "summary.txt"))
  }
  r <- compare_models(root, values = "rse_pct", bold_above = c(rse_pct = 30))
  expect_match(as.character(r$kable), 'font-weight: bold;"> 34.4 </td>', fixed = TRUE)
  expect_no_match(as.character(r$kable), 'font-weight: bold;"> 22.9 </td>', fixed = TRUE)
})

test_that("compare_obspred_metrics accepts folders and summaries, and explains what it cannot use", {
  # a read_summary() object leads back to its model folder; without ChartsData the error says so
  s <- read_summary(fixture("summary_warfarin_final.txt"))
  expect_equal(basename(dirname(s$summary_path)), "fixtures")
  expect_error(compare_obspred_metrics(s, s), "ChartsData|not found|exported", ignore.case = TRUE)
  expect_error(compare_obspred_metrics(list(a = 1), list(a = 1)), "Cannot get the prediction-error metrics")
})

test_that("a convergence assessment is reshaped and plotted", {
  one_run <- function(shift) list(
    populationParameters = list(convergence = data.frame(
      convergenceIndicator = 1:5, ka_pop = 1 + shift + seq(0, 0.4, length.out = 5), V_pop = 8 + shift + 0 * 1:5)),
    standardErrors = list(values = c(0.2, 0.3)))
  res <- list(one_run(0), one_run(0.01), one_run(-0.02))
  tidy <- tidy_convergence_assessment(res)
  expect_equal(nrow(tidy$iterations), 3 * 5 * 3)               # runs x iterations x (indicator + 2 parameters)
  expect_equal(levels(tidy$iterations$parameter)[1], "convergenceIndicator")
  expect_equal(nrow(tidy$final), 3 * 2)
  expect_equal(tidy$final$estimate[1:2], c(1.4, 8))
  expect_equal(round(tidy$final$ymax[1], 3), round(1.4 + 1.96 * 0.2, 3))
  expect_s3_class(plot_convergence_assessment(res), "ggplot")
  expect_s3_class(plot_convergence_assessment(res, which = "estimates", ci = TRUE), "ggplot")
  rds <- withr::local_tempfile(fileext = ".rds"); saveRDS(res, rds)
  expect_equal(tidy_convergence_assessment(rds), tidy)
  expect_error(tidy_convergence_assessment(list(list(a = 1))), "getAssessmentResults")
})

test_that("backward_monitor prints one line per adopted model", {
  root <- withr::local_tempdir()
  for (i in 1:2) {
    d <- file.path(root, c("run000_I0_base", "run002_I1_x")[i]); dir.create(d)
    file.copy(fixture(c("summary_warfarin_example.txt", "summary_warfarin_final.txt")[i]), file.path(d, "summary.txt"))
    file.create(paste0(d, ".mlxtran"))
  }
  monitor <- backward_monitor(shrinkage_warn = 30)
  out <- testthat::capture_messages({
    monitor(0, file.path(root, "run000_I0_base.mlxtran"))
    monitor(1, file.path(root, "run002_I1_x.mlxtran"))
  })
  expect_length(out, 2)
  expect_match(out[1], "iteration 0.*0 relation.*OFV 656.85.*BICc 703.77.*shrinkage .*!")   # Tlag shrinkage 50 percent
  expect_match(out[2], "iteration 1.*1 relation.*OFV 631.10.*change -22\\.3")
})
