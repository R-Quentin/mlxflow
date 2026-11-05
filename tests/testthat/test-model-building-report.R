fixture <- function(x) testthat::test_path("fixtures", x)

test_that("the Monolix model building report is parsed", {
  rep <- read_model_building_report(fixture("modelBuilding.txt"))
  expect_equal(nrow(rep$models), 16)
  expect_equal(rep$best_model_id, 7)
  expect_equal(rep$models$covariates[[7]], "V_logtWt")
  expect_setequal(rep$models$covariates[[16]], c("V_logtWt", "Cl_logtWt"))
  expect_equal(rep$models$ofv[7], 627.94)
  expect_equal(rep$models$n_rel[1], 0)
  expect_equal(read_model_building_report(dirname(fixture("modelBuilding.txt")))$best_model_id, 7)
})

test_that("get_monolix_best_model finds the model by number or, failing that, by covariate set", {
  rep <- read_model_building_report(fixture("modelBuilding.txt"))
  models <- tibble::tibble(model_id = c(1L, 7L, 16L), model_name = sprintf("ModelBuilding_%03d", c(1, 7, 16)),
                           covariates = list(character(0), "V_logtWt", c("V_logtWt", "Cl_logtWt")))
  expect_equal(get_monolix_best_model(rep, models)$model_name, "ModelBuilding_007")

  no_seven <- dplyr::filter(models, model_id != 7L) |>
    dplyr::bind_rows(tibble::tibble(model_id = 3L, model_name = "ModelBuilding_003", covariates = list("V_logtWt")))
  best <- get_monolix_best_model(rep, no_seven)
  expect_equal(best$model_name, "ModelBuilding_003")
  expect_match(best$matched_by, "covariate set")
})
