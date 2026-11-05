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
