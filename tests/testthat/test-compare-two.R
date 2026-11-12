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
  r <- compare_obspred_metrics(m1, m2, labels = c("A", "B"))
  expect_equal(unique(r$Prediction), c("Population predictions", "Individual predictions"))   # not alphabetical
  expect_equal(unique(r$metric), c("N", "MPE", "RMSE"))
  expect_equal(r$delta[r$metric == "MPE"], c(1, 1))
  expect_equal(r$delta[r$metric == "RMSE"], c(-0.1, -0.1))
  # accepts the list returned by obspred_diagnostics()
  expect_equal(compare_obspred_metrics(list(metrics = m1), list(metrics = m2), c("A", "B")), r)
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
