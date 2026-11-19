# A fake lixoftConnectors: a deterministic synthetic "Monolix" good enough to exercise the
# backward driver (covariate toggling, result folders, project files, summary.txt).
# Tests replace mlxflow:::.lx by fake_lx(), so Monolix is not needed.
fk <- new.env()
local({
# Fake lixoftConnectors: a deterministic synthetic "Monolix" good enough to exercise
# the driver (covariate toggling, result folders, project files, summary.txt)
eng <- new.env()
eng$effects <- list(V_logtWt = c(25.68, 1), Cl_logtWt = c(4.93, 1), V_SEX = c(5.2, 2), Cl_SEX = c(12.0, 2))
eng$base_ofv <- 600; eng$n_runs <- 0; eng$fail_at_run <- Inf
eng$project <- NULL; eng$dir <- NULL; eng$cov <- NULL
# Variability (test_variability): OFV gain and one extra parameter (gamma) per random effect on OCC
eng$var_effects <- list(ka_OCC = 0.7, V_OCC = 40, Cl_OCC = 9)
eng$gamma_rse <- list(ka_OCC = 30, V_OCC = 20, Cl_OCC = 150)   # RSE (%) of the new gamma
eng$var <- NULL; eng$scenario <- NULL; eng$scenario_runs <- 0; eng$fail_interrupt <- FALSE
.default_var <- function() list(id = c(ka = TRUE, V = TRUE, Cl = TRUE), OCC = c(ka = FALSE, V = FALSE, Cl = FALSE))
.var_flags <- function() { out <- character()
  for (l in names(eng$var)) for (p in names(eng$var[[l]])) if (eng$var[[l]][[p]]) out <- c(out, paste(p, l, sep = "_")); out }

.rel_flags <- function() { out <- list(); for (p in names(eng$cov)) for (cv in names(eng$cov[[p]]))
  out[[paste(p, cv, sep = "_")]] <- eng$cov[[p]][[cv]]; out }
.beta_names <- function() { nm <- character()
  for (k in names(.rel_flags())) if (.rel_flags()[[k]]) {
    parts <- strsplit(k, "_")[[1]]; p <- parts[1]; cv <- parts[2]
    nm <- c(nm, if (cv == "SEX") paste0("beta_", p, "_SEX_G_", c("B", "C")) else paste0("beta_", k)) }
  nm }
.model_ofv <- function() { f <- .rel_flags(); o <- eng$base_ofv
  for (k in names(f)) if (!f[[k]]) o <- o + eng$effects[[k]][1]
  for (k in .var_flags()) if (k %in% names(eng$var_effects)) o <- o - eng$var_effects[[k]]; o }

initializeLixoftConnectors <- function(...) invisible(TRUE)
fake_new_project <- function(path) {
  eng$project <- path; eng$dir <- tools::file_path_sans_ext(path)
  eng$cov <- list(V = c(logtWt = TRUE, SEX = TRUE), Cl = c(logtWt = TRUE, SEX = TRUE))
  eng$var <- .default_var()
  saveProject(path) }
saveProject <- function(path = eng$project) {
  lines <- c(paste0("dir=", eng$dir),
             unlist(lapply(names(eng$cov), function(p) paste0(p, ":", names(eng$cov[[p]]), "=", eng$cov[[p]]))),
             unlist(lapply(names(eng$var), function(l) paste0("@var:", l, ":", names(eng$var[[l]]), "=", eng$var[[l]]))))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, path); eng$project <- path; invisible(TRUE) }
loadProject <- function(path) {
  stopifnot(file.exists(path)); l <- readLines(path)
  eng$dir <- sub("^dir=", "", l[1]); eng$project <- path
  cov <- list(); var <- list()
  for (x in l[-1]) { a <- strsplit(x, "[:=]")[[1]]
    if (a[1] == "@var") var[[a[2]]][a[3]] <- as.logical(a[4]) else cov[[a[1]]][a[2]] <- as.logical(a[3]) }
  eng$cov <- cov; eng$var <- if (length(var)) var else .default_var(); invisible(TRUE) }
setProjectSettings <- function(directory) { eng$dir <- directory; invisible(TRUE) }
getIndividualParameterModel <- function() list(covariateModel = eng$cov, name = c("ka", "V", "Cl"),
                                               variability = lapply(eng$var, as.list))
setIndividualParameterVariability <- function(...) { args <- list(...)
  for (l in names(args)) for (p in names(args[[l]])) eng$var[[l]][[p]] <- args[[l]][[p]]; invisible(TRUE) }
setInitialEstimatesToLastEstimates <- function(...) invisible(TRUE)
setScenario <- function(scenario) { eng$scenario <- scenario; invisible(TRUE) }
runScenario <- function() { eng$scenario_runs <- eng$scenario_runs + 1
  runPopulationParameterEstimation(); runConditionalDistributionSampling(); runLogLikelihoodEstimation(); invisible(TRUE) }
setCovariateModel <- function(...) { args <- list(...); for (p in names(args)) eng$cov[[p]] <- args[[p]]; invisible(TRUE) }
runPopulationParameterEstimation <- function() { eng$n_runs <- eng$n_runs + 1
  if (eng$n_runs >= eng$fail_at_run) {
    # an interruption (Ctrl-C, R crash) is not an error: nothing catches it, the run just stops
    stop(if (isTRUE(eng$fail_interrupt)) structure(class = c("interrupt", "condition"), list(message = "interrupted", call = NULL))
         else "simulated crash") }
  invisible(TRUE) }
runConditionalDistributionSampling <- function() invisible(TRUE)
runLogLikelihoodEstimation <- function(linearization = FALSE) {
  dir.create(eng$dir, showWarnings = FALSE, recursive = TRUE); o <- .model_ofv()
  nm <- sub("\\.mlxtran$", "", basename(eng$project))
  b <- .beta_names()
  writeLines(c("********************************************************************************",
    sprintf("*  %s.mlxtran", nm), "********************************************************************************", "",
    "ESTIMATION OF THE POPULATION PARAMETERS ________________________________________", "",
    "Fixed Effects ----------------------------", "V_pop      :     8.0", "Cl_pop     :     0.13",
    if (length(b)) sprintf("%-10s :     1.0", b), "",
    "Standard Deviation of the Random Effects -", "omega_V    :     0.2", "omega_Cl   :     0.3", "",
    "Error Model Parameters -------------------", "a          :     0.24", "",
    "Elapsed time (seconds):      1.5", "",
    "ESTIMATION OF THE LOG-LIKELIHOOD _______________________________________________", "",
    sprintf("-2 x log-likelihood                       (OFV) :         %.2f", o),
    sprintf("Akaike Information Criteria               (AIC) :         %.2f", o + 10),
    sprintf("Corrected Bayesian Information Criteria  (BICc) :         %.2f", o + 50 + 5.5 * length(b)),
    sprintf("Bayesian Information Criteria             (BIC) :         %.2f", o + 40), "",
    "DATASET INFORMATION", "Number of individuals: 32"), file.path(eng$dir, "summary.txt"))
  invisible(TRUE) }
getEstimatedLogLikelihood <- function() { o <- .model_ofv()
  list(importanceSampling = c(OFV = o, AIC = o + 10, BICc = o + 50 + 5.5 * length(.beta_names()), BIC = o + 40)) }
getEstimatedPopulationParameters <- function() {
  b <- .beta_names(); v <- .var_flags(); v <- v[grepl("_OCC$", v)]; g <- if (length(v)) paste0("gamma_", v) else character()
  c(setNames(rep(1, 5 + length(b) + length(g)), c("V_pop", "Cl_pop", "omega_V", "omega_Cl", "a", b, g))) }
runStandardErrorEstimation <- function(linearization = FALSE) {
  f <- file.path(eng$dir, "summary.txt"); b <- .beta_names()
  if (file.exists(f) && length(b)) {
    l <- readLines(f); i <- grep("ESTIMATION OF THE LOG-LIKELIHOOD", l)
    w <- max(nchar(b))
    rows <- c(sprintf("%-*s%s", w, b[1], "-0.0089253 -0.016995"), sprintf("%*s  0.1 0.2", w, "V_pop"))
    fis <- c("ESTIMATION OF THE FISHER INFORMATION MATRIX ____________________________________", "",
             "Correlation Matrix :", rows, "")
    writeLines(append(l, fis, after = i - 1), f) }
  invisible(TRUE) }
getEstimatedStandardErrors <- function() {
  b <- .beta_names(); se <- vapply(b, function(x) {
    parts <- strsplit(x, "_")[[1]]; k <- paste(parts[2], parts[3], sep = "_"); e <- eng$effects[[k]]
    1 / sqrt(e[1] / e[2]) }, numeric(1))
  g <- .var_flags(); g <- g[grepl("_OCC$", g)]
  list(stochasticApproximation = data.frame(parameter = c(b, if (length(g)) paste0("gamma_", g)), se = c(unname(se), rep(1, length(g))),
                                            rse = c(100 * unname(se), unname(unlist(eng$gamma_rse[g]))))) }
computeChartsData <- function() invisible(TRUE)

}, envir = fk)

fake_lx <- function(fn, ...) get(fn, envir = fk)(...)

# A fresh fake project (all four relations active) in a temporary folder
new_fake_project <- function(dir = withr::local_tempdir(.local_envir = parent.frame())) {
  path <- file.path(dir, "Full_model.mlxtran")
  fk$fake_new_project(path)
  fk$eng$n_runs <- 0; fk$eng$fail_at_run <- Inf; fk$eng$scenario_runs <- 0; fk$eng$fail_interrupt <- FALSE
  fk$eng$effects <- list(V_logtWt = c(25.68, 1), Cl_logtWt = c(4.93, 1), V_SEX = c(5.2, 2), Cl_SEX = c(12.0, 2))
  path
}
