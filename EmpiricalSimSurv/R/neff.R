#' Compute Peto's Effective Sample Size from a survfit Object
#'
#' @param fit_km  A \code{survfit} object (stratified or unstratified).
#' @return Data frame: time, n_risk, n_event, surv, var_s, neff, arm.
#' @export
calc_neff <- function(fit_km) {
  times   <- fit_km$time
  n_risk  <- fit_km$n.risk
  n_event <- fit_km$n.event
  surv    <- fit_km$surv

  denom   <- n_risk * (n_risk - n_event)
  gw_term <- ifelse(denom > 0, n_event / denom, 0)
  cum_gw  <- cumsum(gw_term)
  var_s   <- surv^2 * cum_gw

  neff <- ifelse(var_s > 0 & surv > 0 & surv < 1,
                 surv * (1 - surv) / var_s,
                 NA_real_)
  res <- data.frame(
    time    = times,
    n_risk  = n_risk,
    n_event = n_event,
    surv    = surv,
    var_s   = var_s,
    neff    = neff,
    stringsAsFactors = FALSE
  )

  if (!is.null(fit_km$strata)) {
    res$arm <- rep(names(fit_km$strata), fit_km$strata)
  } else {
    res$arm <- "Single Group"
  }

  res
}
