test_that("p-values written as a bound (\"<1e-10\") are read as numbers and flagged", {
  d <- withr::local_tempdir(); dir.create(file.path(d, "Tests"))
  writeLines(c("eta,covariate,statistics,p-value,test", "eta_Cl,wt,5.2,<1e-10,Pearson",
               "eta_Cl,sex,0.8,0.37,Anova", "eta_V,wt,NaN,NA,Pearson"),
             file.path(d, "Tests", "correlationRandomEffectsCovariates.txt"))
  re <- read_covariate_assessment(d)$random_effects
  expect_type(re$p.value, "double")
  expect_equal(re$p.value[1:2], c(1e-10, 0.37)); expect_equal(re$p_bound, c(TRUE, FALSE, FALSE))
  expect_s3_class(suppressWarnings(plot_covariate_search(d)), "ggplot")          # was an error: log10() of text
  expect_match(as.character(covariate_tests_table(d, which = "search")), "&lt;1.0e-10|<1.0e-10")
})
