fixture <- function(x) testthat::test_path("fixtures", x)

test_that("covariates come from the population block only (Fisher matrix rows are ignored)", {
  # warfarin_final has a Fisher matrix whose widest row label is beta_V_logtWt: it used to
  # add a second, glued "V_logtWt-0.0089253" covariate
  s <- read_summary(fixture("summary_warfarin_final.txt"), full = TRUE)
  expect_equal(s$covariates, "V_logtWt")
  expect_equal(read_summary(fixture("summary_warfarin_final.txt"))$covariates, "V_logtWt")
  expect_equal(read_summary(fixture("summary_warfarin_example.txt"))$covariates, character(0))
})

test_that("criteria are read for importance sampling and linearization summaries", {
  crit <- function(f, x) { s <- read_summary(fixture(f)); s$criteria$value[s$criteria$criterion == x] }
  expect_equal(crit("summary_warfarin_final.txt", "OFV"), 631.10)
  expect_equal(crit("summary_ModelBuilding_007.txt", "BICc"), 678.33)
})

test_that("parameter tables work with and without standard errors", {
  full <- build_param_table(read_summary(fixture("summary_warfarin_final.txt"), full = TRUE), add_criteria = TRUE)
  expect_true(all(c("estimate", "rse_pct") %in% names(full)))
  partial <- build_param_table(read_summary(fixture("summary_ModelBuilding_007.txt"), full = TRUE))
  expect_false("rse_pct" %in% names(partial))
  expect_true("beta_V_logtWt" %in% partial$parameter)
})

test_that("clean_cov_labels drops categorical suffixes only", {
  expect_equal(clean_cov_labels(c("Cl_tTabac_G_Current_smoker", "V_logtWt")), c("Cl_tTabac", "V_logtWt"))
})

test_that("critical_ofv matches the usual chi-square thresholds", {
  expect_equal(round(critical_ofv(0.05, 1), 2), 3.84)
  expect_equal(round(critical_ofv(0.01, c(1, 2)), 2), c(6.63, 9.21))
})

test_that("the likelihood method is read from the line above the OFV", {
  expect_equal(read_summary(fixture("summary_warfarin_final.txt"))$ll_method, "is")
  expect_equal(read_summary(fixture("summary_ModelBuilding_007.txt"))$ll_method, "lin")
})
