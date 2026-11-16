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
