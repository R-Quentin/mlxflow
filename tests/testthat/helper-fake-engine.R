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

.rel_flags <- function() { out <- list(); for (p in names(eng$cov)) for (cv in names(eng$cov[[p]]))
  out[[paste(p, cv, sep = "_")]] <- eng$cov[[p]][[cv]]; out }
.beta_names <- function() { nm <- character()
  for (k in names(.rel_flags())) if (.rel_flags()[[k]]) {
    parts <- strsplit(k, "_")[[1]]; p <- parts[1]; cv <- parts[2]
    nm <- c(nm, if (cv == "SEX") paste0("beta_", p, "_SEX_G_", c("B", "C")) else paste0("beta_", k)) }
  nm }
.model_ofv <- function() { f <- .rel_flags(); o <- eng$base_ofv
  for (k in names(f)) if (!f[[k]]) o <- o + eng$effects[[k]][1]; o }

initializeLixoftConnectors <- function(...) invisible(TRUE)
fake_new_project <- function(path) {
  eng$project <- path; eng$dir <- tools::file_path_sans_ext(path)
  eng$cov <- list(V = c(logtWt = TRUE, SEX = TRUE), Cl = c(logtWt = TRUE, SEX = TRUE))
  saveProject(path) }
saveProject <- function(path) {
  lines <- c(paste0("dir=", eng$dir),
             unlist(lapply(names(eng$cov), function(p) paste0(p, ":", names(eng$cov[[p]]), "=", eng$cov[[p]]))))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, path); eng$project <- path; invisible(TRUE) }
loadProject <- function(path) {
  stopifnot(file.exists(path)); l <- readLines(path)
  eng$dir <- sub("^dir=", "", l[1]); eng$project <- path
  cov <- list()
  for (x in l[-1]) { a <- strsplit(x, "[:=]")[[1]]; cov[[a[1]]][a[2]] <- as.logical(a[3]) }
  eng$cov <- cov; invisible(TRUE) }
setProjectSettings <- function(directory) { eng$dir <- directory; invisible(TRUE) }
getIndividualParameterModel <- function() list(covariateModel = eng$cov)
setCovariateModel <- function(...) { args <- list(...); for (p in names(args)) eng$cov[[p]] <- args[[p]]; invisible(TRUE) }
runPopulationParameterEstimation <- function() { eng$n_runs <- eng$n_runs + 1
  if (eng$n_runs >= eng$fail_at_run) stop("simulated crash"); invisible(TRUE) }
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
  b <- .beta_names(); c(setNames(rep(1, 5 + length(b)), c("V_pop", "Cl_pop", "omega_V", "omega_Cl", "a", b))) }
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
  list(stochasticApproximation = data.frame(parameter = b, se = unname(se), rse = 100 * unname(se))) }
computeChartsData <- function() invisible(TRUE)

}, envir = fk)

fake_lx <- function(fn, ...) get(fn, envir = fk)(...)

# A fresh fake project (all four relations active) in a temporary folder
new_fake_project <- function(dir = withr::local_tempdir(.local_envir = parent.frame())) {
  path <- file.path(dir, "Full_model.mlxtran")
  fk$fake_new_project(path)
  fk$eng$n_runs <- 0; fk$eng$fail_at_run <- Inf
  fk$eng$effects <- list(V_logtWt = c(25.68, 1), Cl_logtWt = c(4.93, 1), V_SEX = c(5.2, 2), Cl_SEX = c(12.0, 2))
  path
}
