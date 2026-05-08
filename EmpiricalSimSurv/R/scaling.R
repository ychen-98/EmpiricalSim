#' Compute Scaling Parameters (tau, delta) for Survival Simulation
#'
#' Derives a log-linear mapping  \eqn{Y = exp(log(X) * tau + delta)} that
#' aligns two quantiles of the historical distribution with two target values.
#'
#' @param hist_data  Historical time data (numeric vector, can be imputed).
#' @param target_val_low,target_val_high  Target time values at the
#'   lower/upper scaling quantile levels.
#' @param q_scale_range  Two quantile levels, e.g. \code{c(0.1, 0.9)}.
#' @param obs_time,status  Optional observed times and event indicators
#'   for KM-based quantile extraction.
#' @param use_km  If TRUE and obs_time/status given, use KM percentiles
#'   for the historical arm instead of empirical quantiles.
#'
#' @return List with \code{tau}, \code{delta}, \code{hist_val_low},
#'   \code{hist_val_high}.
#' @export
compute_scaling <- function(hist_data, target_val_low, target_val_high,
                            q_scale_range = c(0.1, 0.9),
                            obs_time = NULL, status = NULL,
                            use_km = FALSE) {

  if (use_km && !is.null(obs_time) && !is.null(status)) {
    fit_km <- survival::survfit(survival::Surv(obs_time, status) ~ 1)
    hist_val_low  <- stats::quantile(fit_km, probs = q_scale_range[1])$quantile
    hist_val_high <- stats::quantile(fit_km, probs = q_scale_range[2])$quantile
    if (is.na(hist_val_low))  hist_val_low  <- min(obs_time[status == 1])
    if (is.na(hist_val_high)) hist_val_high <- max(obs_time)
  } else {
    hist_val_low  <- stats::quantile(hist_data, q_scale_range[1], names = FALSE)
    hist_val_high <- stats::quantile(hist_data, q_scale_range[2], names = FALSE)
  }

  tau   <- (log(target_val_high) - log(target_val_low)) /
    (log(hist_val_high)   - log(hist_val_low))
  delta <- log(target_val_low) - log(hist_val_low) * tau

  list(tau = tau, delta = delta,
       hist_val_low = hist_val_low, hist_val_high = hist_val_high)
}


#' Compute Simple Range-Based Scaling (tau, delta) for Non-Survival Data
#'
#' \eqn{Y = X * tau + delta} where tau and delta map the range of the
#' reference (global) data to the range of the target (local) data.
#'
#' Only the min and max of both the target and reference are needed.
#'
#' @param target_mins  Numeric vector of target minimums (one per column).
#' @param target_maxs  Numeric vector of target maximums (one per column).
#' @param dat_ref      Reference data frame (same column order).
#' @param col_names    Optional column names for output labelling.
#'
#' @return List with numeric vectors \code{tau} and \code{delta},
#'   one element per column.
#' @export
compute_scaling_range <- function(target_mins, target_maxs, dat_ref,
                                  col_names = NULL) {
  p <- ncol(dat_ref)
  stopifnot(length(target_mins) == p, length(target_maxs) == p)

  tau   <- numeric(p)
  delta <- numeric(p)

  for (i in seq_len(p)) {
    range_ref <- range(dat_ref[, i], na.rm = TRUE)

    denom <- diff(range_ref)
    if (abs(denom) < .Machine$double.eps) {
      tau[i]   <- 1
      delta[i] <- target_mins[i] - range_ref[1]
    } else {
      tau[i]   <- (target_maxs[i] - target_mins[i]) / denom
      delta[i] <- target_mins[i] - range_ref[1] * tau[i]
    }
  }

  nms <- if (!is.null(col_names)) col_names else colnames(dat_ref)
  names(tau)   <- nms
  names(delta) <- nms
  list(tau = tau, delta = delta)
}
