#' EmpiricalSimSurv: Empirical Copula-Based Simulation for Survival and Mixed Data
#'
#' Generates simulated clinical trial data matching target summary statistics
#' without requiring patient-level target data.  Supports continuous, binary,
#' ordinal, and survival endpoints in any combination.
#'
#' @section Main entry points:
#' \describe{
#'   \item{\code{\link{run_simulation}}}{Mixed data types (survival + any combination)}
#'   \item{\code{\link{run_survival_sim}}}{Survival-only pipeline}
#'   \item{\code{\link{impute_censored}}}{KM + parametric tail imputation}
#' }
#'
#' @references Ding Y, Liu Y, Qu Y (2025), Commun Stat Simul Comput.
#'
#' @importFrom MASS mvrnorm
#' @importFrom survival Surv survfit
#' @importFrom stats sd quantile cov pnorm qnorm rnorm runif approx AIC
#' @importFrom flexsurv flexsurvreg
#' @keywords internal
"_PACKAGE"
