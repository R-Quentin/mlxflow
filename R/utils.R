#' Critical OFV difference of a likelihood ratio test
#'
#' Chi-square quantile used to decide whether removing (or adding) a covariate
#' relation changes the -2 log-likelihood significantly. For one degree of
#' freedom, `critical_ofv(0.05, 1)` is 3.84 and `critical_ofv(0.01, 1)` is 6.63.
#'
#' @param alpha Significance level of the likelihood ratio test (e.g. 0.01).
#' @param df Degrees of freedom, i.e. the number of coefficients removed or
#'   added with the relation (1 for a continuous covariate, number of levels
#'   minus 1 for a categorical one).
#'
#' @return The OFV difference above which the change is significant (numeric,
#'   same length as `df`).
#' @examples
#' critical_ofv(0.01, 1)
#' critical_ofv(0.01, 1:3)
#' @export
critical_ofv <- function(alpha, df) qchisq(1 - alpha, df = df)


# Calls a lixoftConnectors function by name. All the connector calls of the package go
# through here: lixoftConnectors is not on CRAN (it ships with MonolixSuite), so it is only
# required when a function really talks to Monolix, and tests can replace this function.
.lx <- function(fn, ...) {
  if (!requireNamespace("lixoftConnectors", quietly = TRUE)) {
    stop("This function needs the 'lixoftConnectors' package, installed with MonolixSuite ",
         "(see the Lixoft documentation: install.packages(<MonolixSuite>/connectors/lixoftConnectors.tar.gz)).",
         call. = FALSE)
  }
  getExportedValue("lixoftConnectors", fn)(...)
}
