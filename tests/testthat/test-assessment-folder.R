# An Assessment folder as Monolix writes it: Assessment_NNN/Assessment_NNN/populationParameters.txt, ...
make_assessment <- function(root, n_runs = 3, trajectories = FALSE) {
  for (i in seq_len(n_runs)) {
    run <- sprintf("Assessment_%03d", i); d <- file.path(root, run, run); dir.create(d, recursive = TRUE)
    dir.create(file.path(d, ".Internals"))
    writeLines(c("parameter,value,CV,se_lin,rse_lin,P2.5_lin,P97.5_lin",
                 sprintf("V_pop,%g,,0.2,2.5,7,9", 8 + i / 100), sprintf("omega_V,%g,13,0.02,15,0.1,0.2", 0.13 + i / 1000)),
               file.path(d, "populationParameters.txt"))
    dir.create(file.path(d, "LogLikelihood")); writeLines(c("criteria,linearization", sprintf("OFV,%g", 600 + i)), file.path(d, "LogLikelihood", "logLikelihood.txt"))
    if (trajectories) {
      dir.create(file.path(d, "ChartsData", "Saem"), recursive = TRUE)
      writeLines(c("iteration,convergenceIndicator,phase,V_pop,omega_V", "1,100,1,5,1", "2,90,1,7,0.5", "3,80,2,8,0.2"),
                 file.path(d, "ChartsData", "Saem", "CvParam.txt"))
    }
  }
  writeLines(c("Parameters,Run_1,Run_2,Run_3", "V_pop,8.01,8.02,8.03", "omega_V,0.131,0.132,0.133", "-2LL,601,602,603"),
             file.path(root, "assessment.txt"))
  root
}

test_that("an Assessment folder is read like the output of getAssessmentResults()", {
  res <- read_assessment(make_assessment(withr::local_tempdir()))
  expect_length(res, 3)
  expect_false(attr(res, "trajectories"))
  tdy <- tidy_convergence_assessment(res)
  expect_equal(as.character(tdy$final$parameter[1:2]), c("V_pop", "omega_V"))
  expect_equal(tdy$final$estimate[tdy$final$parameter == "V_pop"], c(8.01, 8.02, 8.03))
  expect_equal(tdy$final$se[1], 0.2)                                       # standard errors by name
  expect_equal(res[[2]]$loglikelihood$OFV, 602)
  expect_s3_class(plot_convergence_assessment(res, which = "estimates", ci = TRUE), "ggplot")
  expect_error(plot_convergence_assessment(res), "no SAEM trajectories")
  expect_s3_class(plot_convergence_assessment(make_assessment(withr::local_tempdir()), which = "estimates"), "ggplot")   # the folder itself
})

test_that("SAEM trajectories are read when the charts of the runs exist, and assessment.txt is the fallback", {
  res <- read_assessment(make_assessment(withr::local_tempdir(), trajectories = TRUE))
  expect_true(attr(res, "trajectories"))
  expect_equal(nrow(res[[1]]$populationParameters$convergence), 3)
  expect_equal(res[[1]]$populationParameters$nbsmoothingiterations, 1)
  expect_s3_class(plot_convergence_assessment(res), "ggplot")
  only_txt <- withr::local_tempdir(); file.copy(file.path(make_assessment(withr::local_tempdir()), "assessment.txt"), only_txt)
  expect_length(read_assessment(only_txt), 3)
  expect_error(read_assessment(withr::local_tempdir()), "assessment.txt")
})

test_that("run_diagnostics uses the Assessment folder of the model when there is one", {
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  dir <- withr::local_tempdir(); withr::local_dir(dir)
  dir.create("mod"); make_assessment(file.path("mod", "Assessment"))
  writeLines(c("---", "title: stub", "output: html_document", "params:", "  model_dir: x", "  save_output: FALSE", "  output_dir: ''",
               "  prefix: ''", "  assessment: null", "---", "", "```{r}", "cat('ASSESSMENT', is.null(params$assessment), basename(params$assessment))", "```"), "t.Rmd")
  expect_message(suppressWarnings(run_diagnostics("mod", workflow = "t.Rmd")), "Convergence assessment read from")
  expect_match(paste(readLines("mod_report.html"), collapse = " "), "ASSESSMENT FALSE Assessment")
  # a template that does not declare the parameter keeps working
  writeLines(c("---", "title: stub", "output: html_document", "params:", "  model_dir: x", "  save_output: FALSE", "  output_dir: ''", "  prefix: ''", "---", "", "ok"), "t2.Rmd")
  expect_no_error(suppressWarnings(suppressMessages(run_diagnostics("mod", workflow = "t2.Rmd"))))
})
