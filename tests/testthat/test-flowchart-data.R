mk <- function(id, covs, ofv, nbeta = length(covs)) {
  tibble::tibble(model_id = as.integer(id), model_name = sprintf("ModelBuilding_%03d", id), OFV = ofv, AIC = ofv + 10,
                 BICc = ofv + 50, BIC = ofv + 40, covariates = list(covs), n_cov = as.integer(nbeta))
}
scm_models <- function() {
  dplyr::bind_rows(mk(1, character(0), 700), mk(2, "Cl_logtWt", 690), mk(3, "V_logtWt", 680),
                   mk(4, c("V_SEX_G_B", "V_SEX_G_C"), 695, 2), mk(5, c("V_logtWt", "Cl_logtWt"), 672),
                   mk(6, c("V_logtWt", "V_SEX_G_B", "V_SEX_G_C"), 678, 3)) |>
    dplyr::mutate(covariates = purrr::map(covariates, clean_cov_labels),
                  n_rel = purrr::map_int(covariates, ~ length(unique(.x))), iter = n_rel)
}

test_that("alpha is mandatory in every builder", {
  m <- scm_models()
  expect_error(build_scm_data(m), "alpha")
  expect_error(build_cossac_data(m), "alpha")
  expect_error(build_backward_data(m), "alpha")
})

test_that("the SCM builder tests on the right degrees of freedom and follows alpha", {
  m <- scm_models()
  s05 <- build_scm_data(m, alpha = 0.05)
  v_sex <- dplyr::filter(s05$nodes, model_id == 4)           # 2 betas, delta OFV = -5
  expect_equal(v_sex$df, 2); expect_false(v_sex$significant) # 5 < qchisq(0.95, 2) = 5.99
  expect_equal(s05$alpha, 0.05)
  expect_equal(s05$last_best_name, "ModelBuilding_005")
  # a stricter alpha stops one step earlier: Cl_logtWt (delta -8, p = 0.0047) is no longer kept at 0.001
  expect_equal(build_scm_data(m, alpha = 0.001)$last_best_name, "ModelBuilding_003")
})

test_that("the backward builder adopts the highest p-value, like backward_elimination()", {
  m <- scm_models() |> dplyr::filter(model_id %in% c(5, 3, 2, 1))
  b <- build_backward_data(m, alpha = 0.01)
  expect_equal(b$nodes$iter[b$nodes$model_id == 5], 0)
  expect_equal(b$last_best_name, "ModelBuilding_005")        # nothing removable at 0.01
  expect_true(all(c("df", "p_value") %in% names(b$nodes)))
})

test_that("flowcharts are drawn, and a single node does not break the plot", {
  m <- scm_models()
  expect_s3_class(plot_flowchart(build_scm_data(m, alpha = 0.05), title = "scm"), "ggplot")
  one <- build_backward_data(m[1, ], alpha = 0.01)
  expect_equal(nrow(one$nodes), 1)
  expect_s3_class(plot_flowchart(one, title = "single"), "ggplot")
})
