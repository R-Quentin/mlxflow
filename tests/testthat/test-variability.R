# The Monolix side is the fake engine (helper-fake-engine.R): adding a random effect on OCC lowers the OFV by
# V 40, Cl 9 (both significant at 0.05, 1 df) and ka 0.7 (not significant)

run_variability <- function(out, ..., base = new_fake_project(), env = parent.frame()) {
  local_mocked_bindings(.lx = fake_lx, .env = env)
  suppressMessages(test_variability(base, output_dir = out, level = "OCC", ...))
}

test_that("a single pass tests every parameter against the base model and adopts nothing", {
  out <- withr::local_tempdir()
  res <- run_variability(out)
  expect_equal(res$runs$parameter, c(NA, "ka", "V", "Cl"))
  expect_equal(res$runs$run, 0:3)
  expect_equal(round(res$runs$delta_ofv[-1], 2), c(-0.7, -40, -9))
  expect_equal(res$runs$df[-1], c(1, 1, 1))
  expect_equal(res$runs$significant[-1], c(FALSE, TRUE, TRUE))
  expect_false(any(res$runs$adopted, na.rm = TRUE))
  expect_match(res$stop_reason, "single pass")
  expect_equal(basename(res$final_project), "run000_I0_base.mlxtran")
  expect_true(all(file.exists(file.path(out, c("run000_I0_base.mlxtran", "run002_I1_V_OCC.mlxtran", "variability_log.csv")))))
  log <- utils::read.table(file.path(out, "variability_log.csv"), sep = ";", header = TRUE)
  expect_equal(log$run, 0:3)
  expect_true(all(c("critical_ofv", "alpha", "p_value", "ofv", "bicc") %in% names(log)))
})

test_that("the forward selection adopts the lowest p-value at each iteration and stops when nothing is significant", {
  res <- run_variability(withr::local_tempdir(), forward = TRUE)
  expect_equal(res$history$parameter, c(NA, "V", "Cl"))
  expect_equal(res$runs$run, 0:6)                                        # base + 3 + 2 + 1, no gap
  expect_equal(res$runs$iter, c(0, 1, 1, 1, 2, 2, 3))
  expect_equal(res$runs$adopted[!is.na(res$runs$adopted)], c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_match(res$stop_reason, "iteration 3: no random effect is significant")
  expect_match(basename(res$final_project), "run005_I2_Cl_OCC")
  # the second iteration is tested against the model with V, not against the base
  expect_equal(round(res$runs$delta_ofv[res$runs$run == 4], 2), -0.7)
})

test_that("results never go to the folder of the starting project, and parameters can be chosen", {
  base <- new_fake_project(); out <- withr::local_tempdir()
  local_mocked_bindings(.lx = fake_lx)
  before <- readLines(base)
  res <- suppressMessages(test_variability(base, out, level = "OCC", parameters = c("V", "ka")))
  expect_equal(res$runs$parameter, c(NA, "V", "ka"))
  expect_equal(readLines(base), before)                                  # the original project is not saved over
  expect_error(suppressMessages(test_variability(base, out, level = "OCC")), "already contains")
  expect_no_error(suppressMessages(test_variability(base, out, level = "OCC", overwrite = TRUE)))
})

test_that("level and parameters are checked, and the default level is the first one that is not the individual level", {
  local_mocked_bindings(.lx = fake_lx)
  base <- new_fake_project()
  expect_error(suppressMessages(test_variability(base, withr::local_tempdir(), level = "XYZ")), "Unknown variability level")
  expect_error(suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", parameters = "zz")), "Unknown parameter")
  res <- test_variability(base, withr::local_tempdir(), parameters = "V") |> suppressMessages()
  expect_equal(res$runs$level[1], "OCC")
  # only the individual level: IIV must be asked for explicitly
  fk$eng$var <- list(id = c(ka = TRUE, V = TRUE, Cl = TRUE)); fk$saveProject(base)
  expect_error(suppressMessages(test_variability(base, withr::local_tempdir())), "use level = 'id'")
})

test_that("a parameter that already has the random effect is skipped", {
  local_mocked_bindings(.lx = fake_lx)
  base <- new_fake_project(); fk$eng$var$OCC[["V"]] <- TRUE; fk$saveProject(base)
  msgs <- testthat::capture_messages(res <- test_variability(base, withr::local_tempdir(), level = "OCC"))
  expect_match(paste(msgs, collapse = ""), "Skipping .*V")
  expect_equal(res$runs$parameter, c(NA, "ka", "Cl"))
})

test_that("an interrupted run resumes after the last estimated model and gives the same result", {
  out <- withr::local_tempdir(); base <- new_fake_project()
  local_mocked_bindings(.lx = fake_lx)
  fk$eng$fail_at_run <- 5; fk$eng$fail_interrupt <- TRUE              # base = 1, iteration 1 = 2-4, killed on run 4 of iteration 2
  interrupted <- FALSE
  suppressMessages(tryCatch(test_variability(base, out, level = "OCC", forward = TRUE), interrupt = function(e) interrupted <<- TRUE))
  expect_true(interrupted)
  expect_error(suppressMessages(test_variability(base, out, level = "OCC", forward = TRUE)), "resume = TRUE")
  expect_error(suppressMessages(test_variability(base, out, level = "OCC", forward = TRUE, alpha = 0.01, resume = TRUE)), "settings differ")

  fk$eng$fail_at_run <- Inf; fk$eng$fail_interrupt <- FALSE
  res <- suppressMessages(test_variability(base, out, level = "OCC", forward = TRUE, resume = TRUE))
  ref <- suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", forward = TRUE))
  expect_equal(res$runs[c("run", "iter", "parameter", "delta_ofv", "adopted")], ref$runs[c("run", "iter", "parameter", "delta_ofv", "adopted")])
  expect_equal(res$history$parameter, c(NA, "V", "Cl"))
  expect_equal(basename(res$final_project), basename(ref$final_project))
  log <- utils::read.table(file.path(out, "variability_log.csv"), sep = ";", header = TRUE)
  expect_false(anyDuplicated(log$run) > 0)                           # no model estimated or logged twice
  again <- suppressMessages(test_variability(base, out, resume = TRUE))   # finished: stored result
  expect_equal(again$runs, res$runs)
})

test_that("the charts data are computed once the project is saved", {
  calls <- character()
  local_mocked_bindings(.lx = function(fun, ...) { calls <<- c(calls, fun); fake_lx(fun, ...) })
  suppressMessages(test_variability(new_fake_project(), withr::local_tempdir(), level = "OCC", parameters = "V"))
  i <- which(calls == "computeChartsData")[1]
  expect_equal(calls[(i - 1):i], c("saveProject", "computeChartsData"))
})

test_that("a failing estimation is reported and does not stop the search", {
  out <- withr::local_tempdir(); base <- new_fake_project()
  local_mocked_bindings(.lx = fake_lx)
  fk$eng$fail_at_run <- 3                                                # base = 1, ka = 2, V crashes
  expect_warning(res <- suppressMessages(test_variability(base, out, level = "OCC", parameters = c("ka", "V"))), "failed")
  fk$eng$fail_at_run <- Inf
  expect_match(res$runs$status[res$runs$parameter %in% "V"], "error")
  expect_equal(res$runs$status[res$runs$parameter %in% "ka"], "ok")
})

test_that("a scenario runs for every model, and is checked before anything is fitted", {
  base <- new_fake_project(); local_mocked_bindings(.lx = fake_lx)
  sc <- list(tasks = c(populationParameterEstimation = TRUE, conditionalDistributionSampling = TRUE,
                       logLikelihoodEstimation = TRUE, standardErrorEstimation = FALSE), linearization = FALSE)
  res <- suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", parameters = c("V", "Cl"), scenario = sc))
  expect_equal(fk$eng$scenario_runs, 3)                                  # base + 2 candidates
  expect_equal(res$runs$run, 0:2)
  bad <- sc; bad$tasks[["logLikelihoodEstimation"]] <- FALSE
  expect_error(suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", scenario = bad)),
               "must enable: logLikelihoodEstimation")
  expect_error(check_scenario(list(a = 1)), "`tasks`")
})

test_that("load_project_quiet checks the file and wraps Monolix errors", {
  local_mocked_bindings(.lx = fake_lx)
  expect_error(load_project_quiet("no/such/file.mlxtran"), "Project file not found")
  expect_true(load_project_quiet(new_fake_project()))
  local_mocked_bindings(.lx = function(...) stop("boom"))
  expect_error(load_project_quiet(new_fake_project()), "Failed to load the Monolix project.*boom")
})

test_that("an OFV that is not a finite number is flagged, never adopted, and does not stop the search", {
  base <- new_fake_project(); local_mocked_bindings(.lx = fake_lx)
  old <- list(v = fk$eng$var_effects$V_OCC, ofv = fk$eng$base_ofv)         # the fake engine is shared by the tests: restore it
  withr::defer({ fk$eng$var_effects$V_OCC <- old$v; fk$eng$base_ofv <- old$ofv })
  fk$eng$var_effects$V_OCC <- Inf                                          # the likelihood of the V model breaks down
  res <- suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", forward = TRUE))
  v <- res$runs[res$runs$parameter %in% "V", ]
  expect_match(v$status[1], "^invalid"); expect_true(is.na(v$p_value[1])); expect_false(v$significant[1]); expect_false(v$adopted[1])
  expect_equal(res$history$parameter, c(NA, "Cl"))                         # Cl adopted instead, then nothing is left worth adopting
  expect_true(all(is.finite(res$history$ofv)))
  # a base model without a usable likelihood cannot be tested against
  fk$eng$var_effects$V_OCC <- 40; fk$eng$base_ofv <- Inf
  expect_error(suppressMessages(test_variability(new_fake_project(), withr::local_tempdir(), level = "OCC")), "not a finite number")
})

test_that("a significant effect whose variance is poorly estimated is not adopted when max_rse is given", {
  base <- new_fake_project(); local_mocked_bindings(.lx = fake_lx)
  sc <- list(tasks = c(populationParameterEstimation = TRUE, conditionalDistributionSampling = TRUE,
                       logLikelihoodEstimation = TRUE, standardErrorEstimation = TRUE), linearization = FALSE)
  expect_error(suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", max_rse = 100)), "standardErrorEstimation")
  res <- suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", forward = TRUE, scenario = sc, max_rse = 100))
  expect_equal(res$history$parameter, c(NA, "V"))                          # Cl is significant but its gamma has an RSE of 150 %
  cl <- res$runs[res$runs$parameter %in% "Cl", ]
  expect_true(all(cl$significant)); expect_false(any(cl$eligible)); expect_equal(cl$rse_new[1], 150)
  expect_match(res$stop_reason, "well estimated")
  no_filter <- suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", forward = TRUE, scenario = sc))
  expect_equal(no_filter$history$parameter, c(NA, "V", "Cl"))
})

test_that("the pre-selection drops the parameters that showed nothing in the first iteration", {
  base <- new_fake_project(); local_mocked_bindings(.lx = fake_lx)
  res  <- suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", forward = TRUE, prune_alpha = 0.05))
  expect_equal(res$runs$run, 0:4)                                          # base + 3 (iteration 1) + Cl: ka is not tested again
  expect_equal(res$runs$parameter[res$runs$iter == 2], "Cl")               # ka (p = 0.4) is not tested again
  expect_equal(res$history$parameter, c(NA, "V", "Cl"))
  expect_error(suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", prune_alpha = 0.05)), "forward = TRUE")
})

test_that("the flowchart data of a forward variability test follow the adopted models", {
  out <- withr::local_tempdir(); local_mocked_bindings(.lx = fake_lx)
  res <- suppressMessages(test_variability(new_fake_project(), out, level = "OCC", forward = TRUE))
  d <- build_variability_data(res)
  expect_equal(d$type, "variability"); expect_equal(d$alpha, 0.05)
  expect_equal(d$nodes$added_cov[d$nodes$iter == 0], "Base model")
  expect_equal(d$nodes$iter, c(0, 1, 1, 1, 2, 2, 3))
  expect_equal(d$nodes$added_cov[d$nodes$iter == 1], c("+ka", "+V", "+Cl"))
  expect_equal(d$nodes$is_best, c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))       # V then Cl adopted
  expect_equal(d$nodes$significant[d$nodes$iter == 1], c(FALSE, TRUE, TRUE))
  # arrows: base -> iteration 1, adopted V (run 2) -> iteration 2, adopted Cl (run 5) -> iteration 3
  expect_equal(as.data.frame(d$edges), data.frame(from = c(0, 0, 0, 2, 2, 5), to = c(1, 2, 3, 4, 5, 6)), ignore_attr = TRUE)
  expect_match(d$last_best_name, "run005_I2_Cl_OCC")
  expect_s3_class(plot_flowchart(d, title = "IOV"), "ggplot")
  # from the folder, also for a state that is not finished
  expect_equal(build_variability_data(out)$nodes, d$nodes)
  expect_equal(build_variability_data(res$runs, alpha = 0.05)$edges, d$edges)
  expect_error(build_variability_data(res$runs), "`alpha` is required")
  expect_error(build_variability_data(withr::local_tempdir()), "No variability_state.rds")
})

test_that("a single pass is drawn from the base model, with the most significant model highlighted", {
  local_mocked_bindings(.lx = fake_lx)
  res <- suppressMessages(test_variability(new_fake_project(), withr::local_tempdir(), level = "OCC"))
  d <- build_variability_data(res)
  expect_equal(d$edges$from, c(0, 0, 0)); expect_equal(d$nodes$is_best, c(FALSE, FALSE, TRUE, FALSE))   # V has the lowest p
  p <- plot_flowchart(d)
  expect_s3_class(p, "ggplot")
})

test_that("invalid, failed and imprecise models are labelled in the flowchart", {
  base <- new_fake_project(); local_mocked_bindings(.lx = fake_lx)
  old <- list(v = fk$eng$var_effects$V_OCC); withr::defer(fk$eng$var_effects$V_OCC <- old$v)
  fk$eng$var_effects$V_OCC <- Inf
  sc <- list(tasks = c(populationParameterEstimation = TRUE, conditionalDistributionSampling = TRUE,
                       logLikelihoodEstimation = TRUE, standardErrorEstimation = TRUE), linearization = FALSE)
  res <- suppressMessages(test_variability(base, withr::local_tempdir(), level = "OCC", forward = TRUE, scenario = sc, max_rse = 100))
  d <- build_variability_data(res)
  expect_true(any(grepl("\\(invalid\\)", d$nodes$added_cov)))
  expect_true(any(grepl("\\(RSE 150%\\)", d$nodes$added_cov)))
  expect_false(any(d$nodes$is_best & d$nodes$parameter %in% "V"))
  expect_s3_class(plot_flowchart(d), "ggplot")
})
