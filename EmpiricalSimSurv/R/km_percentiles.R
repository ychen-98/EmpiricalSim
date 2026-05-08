#' Extract KM-Based Percentiles and Auto-Assign Pipeline Roles
#'
#' Fits a KM curve to the supplied data and evaluates it at the requested
#' quantile levels.  Returns a table of (level, time, reachable) plus
#' automatic role assignment for the pipeline (scaling endpoints,
#' shape-matching endpoints, and median).
#'
#' @section Role naming convention:
#'   \code{val_*} = time values.
#'   \code{q_*}   = quantile levels (probabilities).
#'
#' @param obs_time  Numeric event/censor times.
#' @param status    Event indicator (1 = event).
#' @param q_probs   Quantile levels to evaluate.
#' @return List: percentiles (data.frame), stats, roles, neff_table.
#' @export
extract_km_percentiles <- function(obs_time, status,
                                   q_probs = c(0.10, 0.30, 0.50, 0.70, 0.90)) {

  fit_km <- survival::survfit(survival::Surv(obs_time, status) ~ 1)
  neff_table   <- calc_neff(fit_km)
  km_quantiles <- stats::quantile(fit_km, probs = q_probs)$quantile

  min_surv  <- min(fit_km$surv)
  reachable <- (1 - q_probs) >= min_surv

  res_df <- data.frame(q_level = q_probs,
                       time_val = as.numeric(km_quantiles),
                       reachable = reachable)

  # Auto-assign roles
  usable <- res_df[res_df$reachable & !is.na(res_df$time_val), ]
  n_u    <- nrow(usable)

  if (n_u >= 3) {
    i_SL <- 1;  i_SH <- n_u;  i_Med <- ceiling(n_u / 2)
    i_Low  <- max(1, i_Med - 1);  i_High <- min(n_u, i_Med + 1)
    if (n_u == 3) { i_Low <- 1; i_High <- 3 }

    roles <- list(
      val_SL = usable$time_val[i_SL],   q_SL = usable$q_level[i_SL],
      val_low = usable$time_val[i_Low],  q_low = usable$q_level[i_Low],
      val_med = usable$time_val[i_Med],  q_med = usable$q_level[i_Med],
      val_high = usable$time_val[i_High], q_high = usable$q_level[i_High],
      val_SH = usable$time_val[i_SH],   q_SH = usable$q_level[i_SH])
  } else {
    roles <- list(val_SL = NA, q_SL = NA, val_low = NA, q_low = NA,
                  val_med = NA, q_med = NA, val_high = NA, q_high = NA,
                  val_SH = NA, q_SH = NA)
    warning("Fewer than 3 usable KM percentiles; role assignment unavailable.")
  }

  list(percentiles = res_df,
       stats = list(min_surv_reached = min_surv,
                    max_q_available  = 1 - min_surv,
                    max_time_observed = max(obs_time),
                    pct_censored = mean(status == 0)),
       roles = roles, neff_table = neff_table)
}


#' Compute Safe Quantile Levels from a KM Curve (N_eff-Based)
#'
#' @param obs_time,status  Event/censor data.
#' @param n_quantiles  Number of grid points (default 5).
#' @param margin_frac  Fraction to trim from each end (default 0.1).
#' @param min_neff     N_eff threshold (default 10).
#' @param digits       Rounding precision.
#' @return Numeric vector of safe quantile levels.
#' @export
auto_quantile_levels <- function(obs_time, status, n_quantiles = 5,
                                 margin_frac = 0.1, min_neff = 10,
                                 digits = 2) {

  fit_km  <- survival::survfit(survival::Surv(obs_time, status) ~ 1)
  neff_df <- calc_neff(fit_km)

  event_mask <- fit_km$n.event > 0
  neff_ok    <- !is.na(neff_df$neff) & neff_df$neff >= min_neff
  usable     <- which(event_mask & neff_ok)

  min_surv <- if (length(usable) > 0) fit_km$surv[max(usable)] else min(fit_km$surv)
  max_q    <- 1 - min_surv

  margin <- max(round(max_q * margin_frac, digits), 0.005)
  lower  <- margin;  upper <- max_q - margin
  if (upper <= lower) { lower <- max_q * 0.05; upper <- max_q * 0.95 }

  q_probs <- unique(round(seq(lower, upper, length.out = n_quantiles), digits))
  q_probs
}


#' Assign Pipeline Roles from a Pre-Specified Percentile Table
#'
#' When the user supplies target percentiles directly (no patient-level data),
#' this function assigns the same role structure that
#' \code{extract_km_percentiles()} would produce: scaling endpoints (SL, SH),
#' shape-matching endpoints (low, high), and the median (med).
#'
#' @param target_percentiles  A data.frame with columns:
#'   \describe{
#'     \item{q_level}{Quantile levels, e.g. \code{c(0.10, 0.30, 0.50, 0.70, 0.90)}}
#'     \item{time_val}{Target survival times at those levels}
#'   }
#'   Must contain at least 3 rows.  Rows are assumed "usable" (reachable).
#'
#' @return A list with the same \code{roles} structure as
#'   \code{extract_km_percentiles()}: \code{val_SL}, \code{q_SL},
#'   \code{val_low}, \code{q_low}, \code{val_med}, \code{q_med},
#'   \code{val_high}, \code{q_high}, \code{val_SH}, \code{q_SH}.
#' @export
assign_roles_from_percentiles <- function(target_percentiles) {

  stopifnot(is.data.frame(target_percentiles))
  stopifnot(all(c("q_level", "time_val") %in% colnames(target_percentiles)))

  tp <- target_percentiles[order(target_percentiles$q_level), , drop = FALSE]
  tp <- tp[!is.na(tp$time_val), , drop = FALSE]
  n_u <- nrow(tp)

  if (n_u < 3) {
    stop("target_percentiles must have at least 3 rows with non-NA time_val.")
  }

  i_SL  <- 1
  i_SH  <- n_u
  i_Med <- ceiling(n_u / 2)
  i_Low  <- max(1, i_Med - 1)
  i_High <- min(n_u, i_Med + 1)
  if (n_u == 3) { i_Low <- 1; i_High <- 3 }

  list(
    val_SL  = tp$time_val[i_SL],   q_SL  = tp$q_level[i_SL],
    val_low = tp$time_val[i_Low],  q_low = tp$q_level[i_Low],
    val_med = tp$time_val[i_Med],  q_med = tp$q_level[i_Med],
    val_high = tp$time_val[i_High], q_high = tp$q_level[i_High],
    val_SH  = tp$time_val[i_SH],   q_SH  = tp$q_level[i_SH]
  )
}
