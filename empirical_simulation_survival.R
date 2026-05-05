################################################################################################
####  empirical_simulation.R
####  Unified Empirical Simulation Framework
####  Supports: continuous, binary, ordinal, survival data types
####  Use cases:
####    1. Survival data alone
####    2. Mixed data (survival + continuous/binary/ordinal)
####    3. Non-survival data only (continuous, binary, ordinal)
####
####  Primary interface: summary statistics (target_means/target_sds,
####    target_percentiles).  Raw target data (dat_target) accepted as
####    optional backup for automatic extraction.
####
####  R version >= 4.1.2 
####  Last update: 05/05/2026
####  Reference: Ding Y, Liu Y, Qu Y (2025), Commun Stat Simul Comput.
################################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(MASS)
  library(survival)
})

# =============================================================================
# 1.  CORE ENGINE — GRID SEARCH FOR OPTIMAL (alpha, beta)
# =============================================================================

#' Find Optimal Alpha/Beta via Two-Stage Adaptive Grid Search (Vectorized)
#'
#' Maps a historical reference distribution to match target summary statistics
#' using a power-distortion transformation  u^alpha * (1-u)^beta.
#'
#' @section Data types:
#'   continuous / binary / ordinal  -> match target mean & SD
#'   survival                       -> match target percentiles
#'
#' @section Naming convention (survival):
#'   q_low, q_med, q_high           = quantile LEVELS (e.g. 0.30, 0.50, 0.70)
#'   target_low, target_med, target_high = target VALUES at those levels
#'
#' @param a1,a2,a_step  Alpha search range and step size (Stage 1).
#' @param b1,b2,b_step  Beta search range and step size (Stage 1).
#' @param relax_index  Relaxation factor for the two-step selection filter.
#' @param search_index_alpha  Width of Stage 2 refinement window.
#' @param tau,delta  Linear (mean/SD types) or log-scale (survival) shift.
#' @param datref  Reference data — data.frame (first column used) or vector.
#' @param tarmean,tarsd  Target mean and SD (required for non-survival).
#' @param type  One of "continuous", "binary", "ordinal", "survival".
#' @param target_med,target_low,target_high  Target time values (survival).
#' @param q_low,q_med,q_high  Quantile levels (default 0.30, 0.50, 0.70).
#'
#' @return A list with Stage 1/2 grids, best parameters, diagnostics.

find_alpha_beta <- function(a1 = 0.2, a2 = 2.0, a_step = 0.01,
                            b1 = 0, b2 = 0.150, b_step = 0.001,
                            relax_index = 3,
                            search_index_alpha = 5,
                            tau, delta, datref,
                            tarmean = NULL, tarsd = NULL,
                            type,
                            target_med = NULL, target_low = NULL,
                            target_high = NULL,
                            q_low = 0.30, q_med = 0.50, q_high = 0.70) {
  
  # --- Input Validation ----------------------------------------------------
  valid_types <- c("continuous", "binary", "ordinal", "survival")
  if (!type %in% valid_types) {
    stop("type must be one of: ", paste0("'", valid_types, "'", collapse = ", "))
  }
  
  if (is.data.frame(datref) || is.list(datref)) {
    var_name <- names(datref)[1]
    datk     <- as.numeric(datref[[1]])
  } else {
    var_name <- deparse(substitute(datref))
    datk     <- as.numeric(datref)
  }
  
  n0 <- length(datk)
  if (n0 < 2) stop("datref must contain at least 2 numeric values.")
  
  if (type == "survival") {
    if (is.null(target_med) || is.null(target_low) || is.null(target_high))
      stop("For type='survival', target_med, target_low, and target_high are required.")
  } else {
    if (is.null(tarmean) || is.null(tarsd))
      stop("For type='", type, "', tarmean and tarsd are required.")
    if (!is.finite(tarsd) || tarsd <= 0) stop("tarsd must be a positive finite number.")
  }
  
  # Sort for continuous & survival; keep original for binary/ordinal
  if (type %in% c("continuous", "survival")) datk <- sort(datk)
  
  # Uniform grid — precompute logs once
  w       <- seq(0.0001, 1 - 0.0005, by = 0.0005)
  nw      <- length(w)
  log_w   <- log(w)
  log_1mw <- log(1 - w)
  
  # --- Precompute ordinal cumulative probs (if needed) ---------------------
  if (type == "ordinal") {
    tab     <- table(datk)
    ord_p   <- as.numeric(tab) / n0
    ord_cuts <- as.numeric(names(tab))
    ord_K   <- length(ord_cuts)
    ord_pcum <- cumsum(ord_p)
  }
  
  # =========================================================================
  # VECTORIZED GRID EVALUATOR
  # =========================================================================
  
  eval_grid <- function(alphas, betas) {
    ng <- length(alphas)
    
    X <- exp(outer(log_w, alphas) + outer(log_1mw, betas))
    
    if (type == "continuous") {
      IDX <- floor(n0 * X)
      IDX[IDX < 1]  <- 1
      IDX[IDX > n0] <- n0
      VALS <- matrix(datk[IDX], nrow = nw, ncol = ng)
      
      m_vec  <- colMeans(VALS)
      m2_vec <- colMeans(VALS^2)
      sd_vec <- sqrt(pmax(m2_vec - m_vec^2, 0))
      
      data.frame(mean.e = m_vec * tau + delta,
                 sd.e   = sd_vec * tau)
      
    } else if (type == "binary") {
      p_obs <- sum(datk) / n0
      VALS <- (X < p_obs) * 1.0
      m_vec  <- colMeans(VALS)
      m2_vec <- colMeans(VALS^2)
      sd_vec <- sqrt(pmax(m2_vec - m_vec^2, 0))
      
      data.frame(mean.e = m_vec * tau + delta,
                 sd.e   = sd_vec * tau)
      
    } else if (type == "ordinal") {
      VALS <- matrix(ord_cuts[ord_K], nrow = nw, ncol = ng)
      for (k in ord_K:1) {
        mask <- X < ord_pcum[k]
        VALS[mask] <- ord_cuts[k]
      }
      
      m_vec  <- colMeans(VALS)
      m2_vec <- colMeans(VALS^2)
      sd_vec <- sqrt(pmax(m2_vec - m_vec^2, 0))
      
      data.frame(mean.e = m_vec * tau + delta,
                 sd.e   = sd_vec * tau)
      
    } else if (type == "survival") {
      IDX <- floor(n0 * X)
      IDX[IDX < 1]  <- 1
      IDX[IDX > n0] <- n0
      VALS <- matrix(exp(log(datk[IDX]) * tau + delta), nrow = nw, ncol = ng)
      
      q_fun <- function(col) {
        quantile(col, probs = c(q_med, q_low, q_high), na.rm = TRUE, names = FALSE)
      }
      Q <- apply(VALS, 2, q_fun)
      
      data.frame(med_est  = Q[1, ],
                 low_est  = Q[2, ],
                 high_est = Q[3, ])
    }
  }
  
  # =========================================================================
  # STAGE 1 — Coarse Grid
  # =========================================================================
  al1_grid <- seq(a1, a2, by = a_step)
  bl1_grid <- if (b2 == 0) 0 else seq(b1, b2, by = b_step)
  G1 <- expand.grid(Alpha = al1_grid, Beta = bl1_grid, KEEP.OUT.ATTRS = FALSE)
  
  est1 <- eval_grid(G1$Alpha, G1$Beta)
  G1   <- cbind(G1, est1)
  
  if (type == "survival") {
    G1$diff_med <- G1$med_est - target_med
    best1_idx   <- which.min(abs(G1$diff_med))
    best1       <- G1[best1_idx, , drop = FALSE]
    alpha_star  <- best1$Alpha[1]
    error_1     <- abs(best1$diff_med[1])
    width_a     <- a_step
  } else {
    G1$diff_m   <- G1$mean.e - tarmean
    G1$diff_s_p <- G1$sd.e / tarsd
    best1_idx   <- which.min(abs(G1$diff_m))
    best1       <- G1[best1_idx, , drop = FALSE]
    alpha_star  <- best1$Alpha[1]
    width_a     <- abs(best1$diff_s_p[1] - 1)
    error_1     <- abs(best1$diff_m[1])
  }
  
  # Diagnostics
  msgs <- character(0)
  if (isTRUE(all.equal(alpha_star, a1)))
    msgs <- c(msgs, "Alpha hit lower bound a1; consider lowering a1.")
  if (isTRUE(all.equal(alpha_star, a2)))
    msgs <- c(msgs, "Alpha hit upper bound a2; consider increasing a2.")
  if (type != "survival" && is.finite(width_a) && width_a > 0.4)
    msgs <- c(msgs, "Substantial SD mismatch (sd_ratio > 1.4 or < 0.6); consider data transformation.")
  
  # =========================================================================
  # STAGE 2 — Refined Grid
  # =========================================================================
  if (type == "survival") {
    a1_ref <- max(a1, alpha_star - search_index_alpha * a_step)
    a2_ref <- min(a2, alpha_star + search_index_alpha * a_step)
    as_ref <- a_step / 10
    b1_ref <- b1;  b2_ref <- b2;  bs_ref <- if (b_step > 0) b_step else 0.001
  } else { 
    as_ref <- 0.001
    b1_ref <- 0;  b2_ref <- min(width_a / 10, 0.5);  bs_ref <- 0.001
    a1_ref <- max(a1, alpha_star - search_index_alpha * width_a)
    a2_ref <- min(a2, alpha_star + search_index_alpha * width_a)  
  }
  
  al2_grid <- seq(a1_ref, a2_ref, by = as_ref)
  bl2_grid <- seq(b1_ref, max(b1_ref, b2_ref), by = bs_ref)
  G2 <- expand.grid(Alpha = al2_grid, Beta = bl2_grid, KEEP.OUT.ATTRS = FALSE)
  
  chunk_size <- 10000L
  ng2 <- nrow(G2)
  
  if (ng2 <= chunk_size) {
    est2 <- eval_grid(G2$Alpha, G2$Beta)
  } else {
    chunks <- split(seq_len(ng2), ceiling(seq_len(ng2) / chunk_size))
    est2_list <- lapply(chunks, function(idx) {
      eval_grid(G2$Alpha[idx], G2$Beta[idx])
    })
    est2 <- do.call(rbind, est2_list)
  }
  
  G2 <- cbind(G2, est2)
  
  if (type == "survival") {
    G2$diff_med  <- G2$med_est  - target_med
    G2$diff_low  <- G2$low_est  - target_low
    G2$diff_high <- G2$high_est - target_high
    G2$diff_perc_total <- abs(G2$diff_low) + abs(G2$diff_high)
    
    G2_keep <- G2[abs(G2$diff_med) <= (relax_index * error_1), , drop = FALSE]
    if (nrow(G2_keep) == 0L) G2_keep <- G2[which.min(abs(G2$diff_med)), , drop = FALSE]
    best2 <- G2_keep[which.min(G2_keep$diff_perc_total), , drop = FALSE]
  } else {
    G2$diff_m   <- G2$mean.e - tarmean
    G2$diff_s_p <- G2$sd.e / tarsd
    G2$diff_s   <- G2$sd.e - tarsd
    
    G2_keep <- G2[abs(G2$diff_m) <= relax_index * error_1, , drop = FALSE]
    if (nrow(G2_keep) == 0L) G2_keep <- best1
    best2 <- G2_keep[which.min(abs(G2_keep$diff_s)), , drop = FALSE]
  }
  
  # =========================================================================
  # OUTPUT
  # =========================================================================
  if (type == "survival") {
    best_summary <- list(
      alpha = best2$Alpha[1], beta = best2$Beta[1],
      est_med = best2$med_est[1], est_low = best2$low_est[1],
      est_high = best2$high_est[1],
      target_med = target_med, target_low = target_low,
      target_high = target_high,
      q_low = q_low, q_med = q_med, q_high = q_high,
      diff_med = best2$diff_med[1],
      diff_perc_total = best2$diff_perc_total[1])
  } else {
    best_summary <- list(
      alpha = best2$Alpha[1], beta = best2$Beta[1],
      est_mean = best2$mean.e[1], est_sd = best2$sd.e[1],
      target_mean = tarmean, target_sd = tarsd,
      mean_diff = best2$diff_m[1], sd_ratio = best2$diff_s_p[1])
  }
  
  guide_txt <- NULL
  if (any(grepl("transform", msgs, ignore.case = TRUE))) {
    guide_txt <- paste0(
      "Transformation guide:\n",
      "- Right-skewed continuous: try log.\n",
      "- Proportion/percent: rescale to (0,1), apply logit.\n",
      "- Ordinal with rare levels: consider collapsing.")
  }
  
  list(
    historical_summary = c(mean = mean(datk), sd = stats::sd(datk)),
    target_summary = if (type == "survival") {
      c(target_med = target_med, target_low = target_low,
        target_high = target_high,
        q_low = q_low, q_med = q_med, q_high = q_high)
    } else { c(mean = tarmean, sd = tarsd) },
    input = c(tau = tau, delta = delta, var_name = var_name, type = type),
    stage1_grid = c(alpha = paste0("[", a1, ", ", a2, "], by=", a_step),
                    beta  = paste0("[", b1, ", ", b2, "], by=", b_step)),
    stage2_grid = c(alpha = paste0("[", round(a1_ref, 6), ", ", round(a2_ref, 6), "], by=", as_ref),
                    beta  = paste0("[", round(b1_ref, 6), ", ", round(b2_ref, 6), "], by=", bs_ref)),
    stage1_results = G1, stage1_best = best1,
    stage2_results = G2, stage2_filtered = G2_keep,
    best_raw = best2, best_summary = best_summary,
    messages = msgs, transformation_guide = guide_txt
  )
}

# =============================================================================
# 2.  SCALING — compute tau & delta
# =============================================================================

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
compute_scaling <- function(hist_data, target_val_low, target_val_high,
                            q_scale_range = c(0.1, 0.9),
                            obs_time = NULL, status = NULL,
                            use_km = FALSE) {
  
  if (use_km && !is.null(obs_time) && !is.null(status)) {
    fit_km <- survfit(Surv(obs_time, status) ~ 1)
    hist_val_low  <- quantile(fit_km, probs = q_scale_range[1])$quantile
    hist_val_high <- quantile(fit_km, probs = q_scale_range[2])$quantile
    if (is.na(hist_val_low))  hist_val_low  <- min(obs_time[status == 1])
    if (is.na(hist_val_high)) hist_val_high <- max(obs_time)
  } else {
    hist_val_low  <- quantile(hist_data, q_scale_range[1], names = FALSE)
    hist_val_high <- quantile(hist_data, q_scale_range[2], names = FALSE)
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
#' @param dat_target  Target data frame (columns = variables).
#' @param dat_ref     Reference data frame (same column order).
#'
#' @return List with numeric vectors \code{tau} and \code{delta},
#'   one element per column.
compute_scaling_range <- function(dat_target, dat_ref) {
  stopifnot(ncol(dat_target) == ncol(dat_ref))
  p <- ncol(dat_target)
  
  tau   <- numeric(p)
  delta <- numeric(p)
  
  for (i in seq_len(p)) {
    range_tar <- range(dat_target[, i], na.rm = TRUE)
    range_ref <- range(dat_ref[, i],    na.rm = TRUE)
    
    denom <- diff(range_ref)
    if (abs(denom) < .Machine$double.eps) {
      tau[i]   <- 1
      delta[i] <- range_tar[1] - range_ref[1]
    } else {
      tau[i]   <- diff(range_tar) / denom
      delta[i] <- range_tar[1] - range_ref[1] * tau[i]
    }
  }
  
  names(tau)   <- colnames(dat_target)
  names(delta) <- colnames(dat_target)
  list(tau = tau, delta = delta)
}


# =============================================================================
# 3.  N_eff — Peto's Effective Sample Size
# =============================================================================

#' Compute Peto's Effective Sample Size from a survfit Object
#'
#' @param fit_km  A \code{survfit} object (stratified or unstratified).
#' @return Data frame: time, n_risk, n_event, surv, var_s, neff, arm.

calc_neff <- function(fit_km) {
  times   <- fit_km$time
  n_risk  <- fit_km$n.risk
  n_event <- fit_km$n.event
  surv    <- fit_km$surv
  K       <- length(times)
  
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


# =============================================================================
# 4.  KM PERCENTILE EXTRACTION
# =============================================================================

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
extract_km_percentiles <- function(obs_time, status,
                                   q_probs = c(0.10, 0.30, 0.50, 0.70, 0.90)) {
  
  fit_km <- survfit(Surv(obs_time, status) ~ 1)
  neff_table   <- calc_neff(fit_km)
  km_quantiles <- quantile(fit_km, probs = q_probs)$quantile
  
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
auto_quantile_levels <- function(obs_time, status, n_quantiles = 5,
                                 margin_frac = 0.1, min_neff = 10,
                                 digits = 2) {
  
  fit_km  <- survfit(Surv(obs_time, status) ~ 1)
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


# =============================================================================
# 4b. ASSIGN ROLES FROM A USER-SUPPLIED PERCENTILE TABLE
# =============================================================================

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

assign_roles_from_percentiles <- function(target_percentiles) {
  
  stopifnot(is.data.frame(target_percentiles))
  stopifnot(all(c("q_level", "time_val") %in% colnames(target_percentiles)))
  
  # Sort by q_level and remove rows with NA time_val
  tp <- target_percentiles[order(target_percentiles$q_level), , drop = FALSE]
  tp <- tp[!is.na(tp$time_val), , drop = FALSE]
  n_u <- nrow(tp)
  
  if (n_u < 3) {
    stop("target_percentiles must have at least 3 rows with non-NA time_val.")
  }
  
  # Same logic as extract_km_percentiles role assignment
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


# =============================================================================
# 5.  IMPUTATION — KM + Parametric Fallback
# =============================================================================

#' Compute a tail-weighted AIC for a flexsurvreg fit.
tail_weighted_aic <- function(fit, obs_time, status,
                              tail_start, multiplier = 3) {
  if (is.null(fit)) return(Inf)
  if (tail_start <= 0 || multiplier <= 1) return(AIC(fit))
  
  dist_name <- fit$dlist$name
  pars <- fit$res[, "est"]
  
  dfn <- fit$dfns$d
  pfn <- fit$dfns$p
  
  n <- length(obs_time)
  ll_i <- numeric(n)
  
  for (j in seq_len(n)) {
    if (status[j] == 1) {
      val <- do.call(dfn, c(list(x = obs_time[j]), as.list(pars), list(log = TRUE)))
      ll_i[j] <- val
    } else {
      val <- do.call(pfn, c(list(q = obs_time[j]), as.list(pars), list(lower.tail = FALSE, log.p = TRUE)))
      ll_i[j] <- val
    }
  }
  
  weights <- rep(1, n)
  weights[obs_time >= tail_start] <- multiplier
  weighted_ll <- sum(weights * ll_i)
  
  k <- length(pars)
  -2 * weighted_ll + 2 * k
}


#' Impute Censored Observations Using KM Tail with Parametric Fallback
#'
#' @param obs_time   Numeric vector of observed times
#' @param status     Event indicator (1 = event, 0 = censored)
#' @param min_neff   Minimum Peto's N_eff required to trust the KM estimate.
#' @param seed       Random seed for reproducibility
#' @param verbose    If TRUE, print diagnostics
#' @param candidate_dists Character vector of flexsurv distribution names.
#' @param tail_weight_fraction Fraction of trunc_time for tail weighting.
#' @param tail_weight_multiplier Extra weight for tail observations.
#' @param blend_window Width of blending window around trunc_time.
#' @param M Number of multiply-imputed datasets to return.
#' @param s_param_floor Minimum S_param(c) for inverse-CDF sampling.
#'
#' @return If M == 1, a list with imputed data and diagnostics.
#'   If M > 1, a list with imputations and shared info.
impute_censored <- function(obs_time, status, seed = 123,
                            min_neff = 10,
                            verbose = TRUE,
                            candidate_dists = c("weibull",
                                                "lognormal",
                                                "llogis",
                                                "gengamma",
                                                "gompertz"),
                            tail_weight_fraction = 0.5,
                            tail_weight_multiplier = 3,
                            blend_window = 0,
                            M = 1,
                            s_param_floor = 1e-6) {
  
  require(survival)
  require(flexsurv)
  
  n <- length(obs_time)
  stopifnot(length(status) == n)
  stopifnot(all(status %in% c(0, 1)))
  stopifnot(all(obs_time > 0))
  
  censored_idx <- which(status == 0)
  
  if (length(censored_idx) == 0) {
    if (verbose) cat("No censored observations. Returning original data.\n")
    single <- list(
      imputed_time          = obs_time,
      method_used           = rep("event", n),
      km_truncation_time    = NA_real_,
      km_surv_at_truncation = NA_real_,
      neff_at_truncation    = NA_real_,
      parametric_dist       = NA_character_,
      parametric_fit        = NULL,
      aic_table             = NULL,
      neff_table            = NULL
    )
    if (M == 1) return(single)
    return(list(
      imputations = replicate(M, single, simplify = FALSE),
      shared      = single
    ))
  }
  
  # ================================================================
  # A. Fit KM curve, compute N_eff, determine truncation point
  # ================================================================
  fit_km  <- survfit(Surv(obs_time, status) ~ 1)
  neff_df <- calc_neff(fit_km)
  
  event_mask <- fit_km$n.event > 0
  neff_ok    <- neff_df$neff >= min_neff
  
  usable <- which(event_mask & neff_ok)
  
  if (length(usable) > 0) {
    trunc_idx   <- max(usable)
    trunc_time  <- fit_km$time[trunc_idx]
    trunc_surv  <- fit_km$surv[trunc_idx]
    trunc_neff  <- neff_df$neff[trunc_idx]
    
    keep        <- seq_len(trunc_idx)
    km_times    <- c(0, fit_km$time[keep])
    km_surv     <- c(1, fit_km$surv[keep])
    km_min_surv <- min(km_surv)
  } else {
    trunc_time  <- 0
    trunc_surv  <- 1
    trunc_neff  <- NA_real_
    km_times    <- c(0)
    km_surv     <- c(1)
    km_min_surv <- 1
  }
  
  if (verbose) {
    cat("--- KM Truncation (N_eff-based) ---\n")
    cat("  min_neff threshold:  ", min_neff, "\n")
    cat("  KM truncated at time:", trunc_time, "\n")
    cat("  S_KM at truncation:  ", round(trunc_surv, 4), "\n")
    cat("  N_eff at truncation: ", round(trunc_neff, 2), "\n")
    if (any(event_mask)) {
      cat("  Original last event time:",
          max(fit_km$time[fit_km$n.event > 0]), "\n\n")
    }
  }
  
  # ================================================================
  # B. Fit parametric candidates, select by tail-weighted AIC
  # ================================================================
  surv_obj <- Surv(obs_time, status)
  
  fits <- setNames(
    lapply(candidate_dists, function(d) {
      tryCatch(flexsurvreg(surv_obj ~ 1, dist = d),
               error = function(e) NULL)
    }),
    candidate_dists
  )
  
  fits <- Filter(Negate(is.null), fits)
  
  if (length(fits) == 0) {
    stop("All parametric fits failed. Check your data.")
  }
  
  tail_start <- trunc_time * tail_weight_fraction
  
  aic_values <- sapply(fits, function(f) {
    tryCatch(
      tail_weighted_aic(f, obs_time, status,
                        tail_start  = tail_start,
                        multiplier  = tail_weight_multiplier),
      error = function(e) Inf
    )
  })
  
  standard_aic <- sapply(fits, function(f) tryCatch(AIC(f), error = function(e) Inf))
  
  aic_table <- data.frame(
    distribution     = names(fits),
    AIC_standard     = standard_aic,
    AIC_tail_weighted = aic_values,
    stringsAsFactors = FALSE
  )
  aic_table <- aic_table[order(aic_table$AIC_tail_weighted), ]
  rownames(aic_table) <- NULL
  
  best_dist <- aic_table$distribution[1]
  param_fit <- fits[[best_dist]]
  
  if (verbose) {
    cat("--- Parametric Model Selection (Tail-Weighted AIC) ---\n")
    cat("  Tail weight region:  time >=", round(tail_start, 2), "\n")
    cat("  Tail weight multiplier:", tail_weight_multiplier, "\n")
    print(aic_table)
    cat("  Selected:", best_dist, "\n\n")
  }
  
  # ================================================================
  # C. Parametric helpers
  # ================================================================
  param_surv_at <- function(t) {
    s <- tryCatch(
      summary(param_fit, t = t, type = "survival")[[1]]$est,
      error = function(e) NA_real_
    )
    if (is.na(s) || !is.finite(s)) return(1e-12)
    max(s, 1e-12)
  }
  
  param_quantile_at <- function(p_surv) {
    q_prob <- 1 - p_surv
    q_prob <- max(1e-12, min(1 - 1e-12, q_prob))
    tryCatch(
      summary(param_fit, type = "quantile", quantiles = q_prob)[[1]]$est,
      error = function(e) NA_real_
    )
  }
  
  param_conditional_median <- function(c_time) {
    s_at_c <- param_surv_at(c_time)
    target_surv <- s_at_c / 2
    if (target_surv < 1e-12) return(c_time * 1.5)
    t_med <- param_quantile_at(target_surv)
    if (is.na(t_med) || !is.finite(t_med)) return(c_time * 1.5)
    t_med
  }
  
  # ================================================================
  # D. Single-imputation engine
  # ================================================================
  impute_once <- function(rng_seed) {
    set.seed(rng_seed)
    
    imputed_time <- obs_time
    method_used  <- rep("event", n)
    
    km_count    <- 0
    param_count <- 0
    blend_count <- 0
    
    blend_lo <- trunc_time - blend_window / 2
    blend_hi <- trunc_time
    
    for (i in censored_idx) {
      c_time <- obs_time[i]
      
      in_blend_zone <- (blend_window > 0) &&
        (c_time >= blend_lo) &&
        (c_time < blend_hi)
      use_km        <- (c_time < blend_lo) ||
        (blend_window == 0 && c_time < trunc_time)
      
      draw_km <- function() {
        s_km_at_c <- approx(km_times, km_surv, xout = c_time,
                            method = "constant", rule = 2, f = 0)$y
        u <- runif(1, 0, s_km_at_c)
        if (u < km_min_surv) {
          return(draw_param())
        }
        t_new <- approx(km_surv, km_times, xout = u,
                        method = "constant", rule = 2, ties = max)$y
        max(t_new, c_time + 0.001)
      }
      
      draw_param <- function() {
        s_at_c <- param_surv_at(c_time)
        if (s_at_c < s_param_floor) {
          t_new <- param_conditional_median(c_time)
          return(max(t_new, c_time + 0.001))
        }
        u     <- runif(1, 0, s_at_c)
        t_new <- param_quantile_at(u)
        if (is.na(t_new) || !is.finite(t_new)) {
          t_new <- param_conditional_median(c_time)
        }
        if (is.na(t_new) || !is.finite(t_new)) {
          t_new <- c_time * 1.5
        }
        max(t_new, c_time + 0.001)
      }
      
      if (use_km) {
        imputed_time[i] <- draw_km()
        method_used[i]  <- "km"
        km_count <- km_count + 1
        
      } else if (in_blend_zone) {
        alpha <- (blend_hi - c_time) / blend_window
        alpha <- max(0, min(1, alpha))
        
        t_km    <- draw_km()
        t_param <- draw_param()
        imputed_time[i] <- alpha * t_km + (1 - alpha) * t_param
        imputed_time[i] <- max(imputed_time[i], c_time + 0.001)
        method_used[i]  <- "blend"
        blend_count <- blend_count + 1
        
      } else {
        imputed_time[i] <- draw_param()
        method_used[i]  <- "parametric"
        param_count <- param_count + 1
      }
    }
    
    if (verbose && M == 1) {
      cat("--- Imputation Summary ---\n")
      cat("  Total censored:", length(censored_idx), "\n")
      cat("  KM-imputed:    ", km_count,
          sprintf("(censored before t=%.1f)", blend_lo), "\n")
      if (blend_window > 0) {
        cat("  Blended:       ", blend_count,
            sprintf("(censored in [%.1f, %.1f))", blend_lo, blend_hi), "\n")
      }
      cat("  Parametric:    ", param_count,
          paste0("(", best_dist, ")"),
          sprintf("(censored at or after t=%.1f)", trunc_time), "\n\n")
    }
    
    list(
      imputed_time          = imputed_time,
      method_used           = method_used,
      km_truncation_time    = trunc_time,
      km_surv_at_truncation = trunc_surv,
      neff_at_truncation    = trunc_neff,
      parametric_dist       = best_dist,
      parametric_fit        = param_fit,
      aic_table             = aic_table,
      neff_table            = neff_df,
      counts                = c(km = km_count, blend = blend_count,
                                parametric = param_count)
    )
  }
  
  # ================================================================
  # E. Run single or multiple imputations
  # ================================================================
  if (M == 1) {
    return(impute_once(seed))
  }
  
  seeds <- seed + seq_len(M) - 1
  imputations <- lapply(seeds, function(s) {
    if (verbose) cat(sprintf("=== Imputation %d (seed=%d) ===\n",
                             which(seeds == s), s))
    impute_once(s)
  })
  
  if (verbose) {
    all_methods <- do.call(rbind, lapply(imputations, function(x) x$counts))
    cat("--- Multiple Imputation Summary (M =", M, ") ---\n")
    cat("  KM-imputed (mean):  ", round(mean(all_methods[, "km"]), 1), "\n")
    if (blend_window > 0) {
      cat("  Blended (mean):     ", round(mean(all_methods[, "blend"]), 1), "\n")
    }
    cat("  Parametric (mean):  ", round(mean(all_methods[, "parametric"]), 1), "\n\n")
  }
  
  shared <- list(
    km_truncation_time    = trunc_time,
    km_surv_at_truncation = trunc_surv,
    neff_at_truncation    = trunc_neff,
    parametric_dist       = best_dist,
    parametric_fit        = param_fit,
    aic_table             = aic_table,
    neff_table            = neff_df
  )
  
  list(imputations = imputations, shared = shared)
}

# =============================================================================
# 6.  DATA GENERATION — multivariate correlated simulation
# =============================================================================

#' Generate Simulated Data via Empirical Distortion
#'
#' @param dat_ref     Reference data frame (columns = variables).
#' @param N           Number of rows to simulate.
#' @param alpha,beta  Named or positional numeric vectors (length = ncol).
#' @param tau,delta   Named or positional numeric vectors (length = ncol).
#' @param types       Character vector of column types.
#' @param cor_method  "rank" or "identity".
#'
#' @return Data frame with N rows, same column names as \code{dat_ref}.
generate_data <- function(dat_ref, N, alpha, beta, tau, delta,
                          types, cor_method = "rank") {
  
  p  <- ncol(dat_ref)
  n0 <- nrow(dat_ref)
  stopifnot(length(alpha) == p, length(beta) == p,
            length(tau) == p, length(delta) == p, length(types) == p)
  
  if (cor_method == "rank") {
    F_u     <- (apply(dat_ref, 2, rank) - 0.5) / n0
    Phi_inv <- qnorm(F_u)
    Sigma   <- cov(Phi_inv)
  } else {
    Sigma <- diag(p)
  }
  
  Z <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = Sigma, tol = 1e-6)
  U <- pnorm(Z)
  
  U_prime <- matrix(NA_real_, N, p)
  for (j in seq_len(p)) {
    U_prime[, j] <- U[, j]^alpha[j] * (1 - U[, j])^beta[j]
  }
  
  dat_sim <- matrix(NA_real_, N, p)
  
  for (j in seq_len(p)) {
    col_data <- dat_ref[, j]
    u_j      <- U_prime[, j]
    
    if (types[j] == "continuous") {
      sorted_vals <- sort(col_data)
      idx <- floor(n0 * u_j); idx[idx < 1] <- 1; idx[idx > n0] <- n0
      dat_sim[, j] <- sorted_vals[idx] * tau[j] + delta[j]
      
    } else if (types[j] == "binary") {
      p_obs <- mean(col_data)
      dat_sim[, j] <- as.numeric(u_j < p_obs)
      
    } else if (types[j] == "ordinal") {
      tab  <- table(col_data)
      probs <- as.numeric(tab) / sum(tab)
      cuts  <- as.numeric(names(tab))
      K     <- length(cuts)
      p_cum <- cumsum(probs)
      
      col_sim <- rep(cuts[K], N)
      for (k in K:1) col_sim[u_j < p_cum[k]] <- cuts[k]
      col_sim[u_j < p_cum[1]] <- cuts[1]
      dat_sim[, j] <- col_sim
      
    } else if (types[j] == "survival") {
      sorted_vals <- sort(col_data)
      idx <- floor(n0 * u_j); idx[idx < 1] <- 1; idx[idx > n0] <- n0
      dat_sim[, j] <- exp(log(sorted_vals[idx]) * tau[j] + delta[j])
      
    } else {
      stop("Unknown type '", types[j], "' for column ", j)
    }
  }
  
  dat_sim <- as.data.frame(dat_sim)
  colnames(dat_sim) <- colnames(dat_ref)
  dat_sim
}


# =============================================================================
# 7.  SURVIVAL-ONLY PIPELINE
# =============================================================================

#' Run Full Empirical Simulation Pipeline (Survival Data)
#'
#' End-to-end pipeline.  Accepts EITHER:
#'   (a) target_percentiles — a data.frame(q_level, time_val) of
#'       pre-specified percentile targets (primary interface), OR
#'   (b) target_time / target_status — raw patient-level data from
#'       which KM percentiles are extracted (backup / convenience).
#'
#' @param target_percentiles  A data.frame with columns \code{q_level}
#'   and \code{time_val}.  At least 3 rows required.  When supplied,
#'   \code{target_time}/\code{target_status} are ignored.
#' @param target_time,target_status  Raw target arm data (optional backup).
#' @param hist_time,hist_status      Historical/reference arm data.
#' @param N_sim           Number of simulated observations.
#' @param n_quantiles     Grid size for auto quantile levels (data path).
#' @param margin_frac     Margin trimming fraction (data path).
#' @param min_neff        N_eff threshold (data path).
#' @param q_probs_manual  If supplied with data path, overrides auto levels.
#' @param use_km_scaling  Use KM-based percentiles for scaling.
#' @param a1,a2,a_step,b1,b2,b_step  Grid search parameters.
#' @param relax_index,search_index_alpha  Two-stage tuning.
#' @param rho_matrix      Optional external correlation matrix.
#' @param verbose         Print diagnostics.
#'
#' @return List: sim_data, scales, shape_params, target_roles,
#'   hist_km, q_probs_used, etc.
run_survival_sim <- function(target_percentiles = NULL,
                             target_time = NULL, target_status = NULL,
                             hist_time, hist_status,
                             N_sim = 5000,
                             n_quantiles = 5,
                             margin_frac = 0.1,
                             min_neff = 10,
                             q_probs_manual = NULL,
                             use_km_scaling = FALSE,
                             a1 = 0.2, a2 = 2.0, a_step = 0.01,
                             b1 = 0, b2 = 0.200, b_step = 0.001,
                             relax_index = 5,
                             search_index_alpha = 7,
                             rho_matrix = NULL,
                             verbose = TRUE) {
  
  # ================================================================
  # PATH A: Pre-specified percentiles (summary statistics only)
  # ================================================================
  if (!is.null(target_percentiles)) {
    
    stopifnot(is.data.frame(target_percentiles))
    stopifnot(all(c("q_level", "time_val") %in% colnames(target_percentiles)))
    
    tp <- target_percentiles[order(target_percentiles$q_level), , drop = FALSE]
    tp <- tp[!is.na(tp$time_val), , drop = FALSE]
    
    if (nrow(tp) < 3)
      stop("target_percentiles must have at least 3 rows with non-NA time_val.")
    
    q_probs <- tp$q_level
    
    if (verbose) cat("--- PERCENTILE-BASED target (no patient data) ---\n",
                     "  Levels:", q_probs, "\n",
                     "  Times: ", tp$time_val, "\n\n")
    
    # Build roles directly
    r <- assign_roles_from_percentiles(tp)
    
    # Build a pseudo target_km result for output compatibility
    target_km <- list(
      percentiles = data.frame(q_level = tp$q_level,
                               time_val = tp$time_val,
                               reachable = TRUE),
      stats = list(min_surv_reached = NA,
                   max_q_available  = max(tp$q_level),
                   max_time_observed = NA,
                   pct_censored = NA),
      roles = r,
      neff_table = NULL
    )
    
  # ================================================================
  # PATH B: Raw target data (backup)
  # ================================================================
  } else if (!is.null(target_time) && !is.null(target_status)) {
    
    # 0. Quantile grid
    if (!is.null(q_probs_manual)) {
      stopifnot(is.numeric(q_probs_manual),
                all(q_probs_manual > 0 & q_probs_manual < 1))
      q_probs <- sort(q_probs_manual)
      if (verbose) cat("--- MANUAL quantile levels:", q_probs, "---\n\n")
    } else {
      q_probs <- auto_quantile_levels(target_time, target_status,
                                      n_quantiles = n_quantiles,
                                      margin_frac = margin_frac,
                                      min_neff = min_neff)
      if (verbose) cat("--- AUTO quantile levels (min_neff =", min_neff, "):",
                       q_probs, "---\n\n")
    }
    
    # 1. KM percentiles
    target_km <- extract_km_percentiles(target_time, target_status,
                                        q_probs = q_probs)
    r <- target_km$roles
    
    if (verbose) {
      cat("--- Target KM ---\n"); print(target_km$percentiles); cat("\n")
    }
    
  } else {
    stop("Provide either target_percentiles (data.frame with q_level, time_val) ",
         "or target_time + target_status.")
  }
  
  # ================================================================
  # Common path: Historical KM, scaling, shape matching, simulation
  # ================================================================
  
  # Historical KM
  hist_km <- extract_km_percentiles(hist_time, hist_status, q_probs = q_probs)
  
  if (verbose) {
    cat("--- Historical KM ---\n"); print(hist_km$percentiles); cat("\n")
  }
  
  # Check roles
  if (any(is.na(c(r$val_SL, r$val_med, r$val_SH))))
    stop("Not enough usable KM percentiles for role assignment. ",
         "Provide at least 3 percentile-time pairs.")
  
  if (verbose) {
    cat("Roles: scaling q_SL=", r$q_SL, "->", r$val_SL,
        ", q_SH=", r$q_SH, "->", r$val_SH, "\n")
    cat("       median  q_med=", r$q_med, "->", r$val_med, "\n")
    cat("       shape   q_low=", r$q_low, "->", r$val_low,
        ", q_high=", r$q_high, "->", r$val_high, "\n\n")
  }
  
  # 3. Scaling
  scales <- compute_scaling(
    hist_data = hist_time, target_val_low = r$val_SL, target_val_high = r$val_SH,
    q_scale_range = c(r$q_SL, r$q_SH),
    obs_time = hist_time, status = hist_status, use_km = use_km_scaling)
  
  if (verbose) cat("Scaling: tau =", scales$tau, ", delta =", scales$delta, "\n\n")
  
  # 4. Shape matching
  shape <- find_alpha_beta(
    a1 = a1, a2 = a2, a_step = a_step, b1 = b1, b2 = b2, b_step = b_step,
    relax_index = relax_index, search_index_alpha = search_index_alpha,
    tau = scales$tau, delta = scales$delta, datref = hist_time,
    type = "survival",
    target_med = r$val_med, target_low = r$val_low, target_high = r$val_high,
    q_low = r$q_low, q_med = r$q_med, q_high = r$q_high)
  
  if (verbose) { cat("--- Shape ---\n"); print(shape$best_summary); cat("\n") }
  
  # 5. Simulate
  n0 <- length(hist_time); sorted_hist <- sort(hist_time)
  if (is.null(rho_matrix)) { Z <- matrix(rnorm(N_sim), ncol = 1) }
  else { Z <- MASS::mvrnorm(N_sim, rep(0, ncol(rho_matrix)), rho_matrix) }
  U <- pnorm(Z)
  a_opt <- shape$best_summary$alpha;  b_opt <- shape$best_summary$beta
  U_prime <- U^a_opt * (1 - U)^b_opt
  idx <- floor(n0 * U_prime); idx[idx < 1] <- 1; idx[idx > n0] <- n0
  sim_vec <- pmax(0.001, exp(log(sorted_hist[idx]) * scales$tau + scales$delta))
  
  dat_sim <- data.frame(sim_time = as.vector(sim_vec), sim_status = 1L)
  
  if (verbose) {
    cat("--- Simulation Complete (N =", nrow(dat_sim), ") ---\n")
    print(round(quantile(dat_sim$sim_time, c(0.05, 0.25, 0.50, 0.75, 0.95)), 4))
    cat("\n")
  }
  
  list(sim_data = dat_sim, scales = scales, shape_params = shape,
       target_km = target_km, hist_km = hist_km,
       q_probs_used = q_probs,
       q_probs_mode = if (!is.null(target_percentiles)) "percentiles"
                      else if (!is.null(q_probs_manual)) "manual"
                      else "auto",
       min_neff_used = min_neff)
}


# =============================================================================
# 8.  HELPER — Extract Target Summaries from Data (backup utility)
# =============================================================================

#' Extract Target Summary Statistics from a Data Frame
#'
#' Convenience function to compute per-column mean and SD from a target
#' data frame, for use as \code{target_means} / \code{target_sds} in
#' \code{run_simulation()}.  Binary columns with zero SD are replaced
#' with the Bernoulli SD.
#'
#' @param dat_target  Target data frame.
#' @param types       Character vector of column types.
#' @return List with vectors \code{means} and \code{sds}.
extract_target_summaries <- function(dat_target, types) {
  p <- ncol(dat_target)
  means <- sds <- numeric(p)
  for (j in seq_len(p)) {
    means[j] <- mean(dat_target[, j], na.rm = TRUE)
    sds[j]   <- sd(dat_target[, j], na.rm = TRUE)
    if (types[j] == "binary" && (is.na(sds[j]) || sds[j] == 0)) {
      sds[j] <- sqrt(means[j] * (1 - means[j]))
    }
  }
  names(means) <- names(sds) <- colnames(dat_target)
  list(means = means, sds = sds)
}


# =============================================================================
# 9.  MIXED-TYPE PIPELINE — the main user-facing function
# =============================================================================

#' Run Full Empirical Simulation (Mixed Data Types)
#'
#' Orchestrates the complete simulation for a data set that may contain
#' any mix of continuous, binary, ordinal, and survival columns.
#'
#' \strong{Primary interface (summary statistics):}
#'   For non-survival columns, supply \code{target_means} and
#'   \code{target_sds} (numeric vectors, one element per non-survival
#'   column).  For survival columns, supply
#'   \code{surv_target_percentiles} (a data.frame with \code{q_level}
#'   and \code{time_val}).
#'
#' \strong{Backup interface (raw data):}
#'   Optionally supply \code{dat_target} — a data frame from which
#'   target means, SDs, and survival KM percentiles are extracted
#'   automatically via \code{extract_target_summaries()} and
#'   \code{extract_km_percentiles()}.
#'
#' @param dat_ref       Reference data frame (historical / global control).
#' @param types         Character vector of column types (length = ncol).
#' @param N_sim         Number of rows to simulate.
#' @param target_means  Numeric vector of target means for non-survival
#'   columns (length = number of non-survival columns).
#' @param target_sds    Numeric vector of target SDs for non-survival
#'   columns (length = number of non-survival columns).
#' @param dat_target    Optional backup: target data frame for automatic
#'   extraction of means, SDs, and survival data.
#' @param scaling_method  \code{"range"} (requires \code{dat_target}),
#'   \code{"summary"} (uses target_means to derive delta; tau = 1), or
#'   \code{"manual"} (supply \code{tau_manual}/\code{delta_manual}).
#' @param tau_manual,delta_manual  Optional manual tau/delta vectors.
#' @param surv_target_percentiles  Data.frame with \code{q_level} and
#'   \code{time_val} for survival columns (primary interface).
#' @param surv_target_time,surv_target_status  Raw target survival vectors
#'   (backup, used only when \code{surv_target_percentiles} is NULL).
#' @param surv_time_col,surv_status_col  Column names for survival
#'   time and status in \code{dat_ref} / \code{dat_target}.
#' @param surv_settings  List of survival pipeline settings passed to
#'   \code{run_survival_sim()}.
#' @param search_settings  List of grid search settings for non-survival
#'   columns passed to \code{find_alpha_beta()}.
#' @param cor_method    \code{"rank"} or \code{"identity"}.
#' @param verbose       Print diagnostics.
#'
#' @return List: sim_data (data.frame), alpha, beta, tau, delta,
#'   type_results (per-column search details), surv_pipeline (if applicable).
run_simulation <- function(dat_ref, types, N_sim = 5000,
                           target_means = NULL, target_sds = NULL,
                           dat_target = NULL,
                           scaling_method = "summary",
                           tau_manual = NULL, delta_manual = NULL,
                           surv_target_percentiles = NULL,
                           surv_target_time = NULL, surv_target_status = NULL,
                           surv_time_col = NULL, surv_status_col = NULL,
                           surv_settings = list(),
                           search_settings = list(),
                           cor_method = "rank",
                           verbose = TRUE) {
  
  p <- ncol(dat_ref)
  stopifnot(length(types) == p)
  col_names <- colnames(dat_ref)
  has_surv  <- any(types == "survival")
  surv_cols <- which(types == "survival")
  nonsurv_cols <- which(types != "survival")
  
  # =======================================================================
  # RESOLVE TARGET SUMMARIES for non-survival columns
  # =======================================================================
  if (length(nonsurv_cols) > 0) {
    if (is.null(target_means) || is.null(target_sds)) {
      if (is.null(dat_target)) {
        stop("For non-survival columns, provide either target_means/target_sds ",
             "or dat_target for automatic extraction.")
      }
      if (verbose) cat("Extracting target means/SDs from dat_target...\n")
      summ <- extract_target_summaries(dat_target[, nonsurv_cols, drop = FALSE],
                                       types[nonsurv_cols])
      target_means <- summ$means
      target_sds   <- summ$sds
    }
    
    # Validate lengths
    if (length(target_means) != length(nonsurv_cols) ||
        length(target_sds)   != length(nonsurv_cols)) {
      stop("target_means and target_sds must have length = ",
           length(nonsurv_cols), " (one per non-survival column).")
    }
    
    # Name them for clarity
    names(target_means) <- col_names[nonsurv_cols]
    names(target_sds)   <- col_names[nonsurv_cols]
  }
  
  # =======================================================================
  # TAU / DELTA for non-survival columns
  # =======================================================================
  if (length(nonsurv_cols) > 0) {
    if (scaling_method == "range") {
      if (is.null(dat_target)) {
        stop("scaling_method='range' requires dat_target. ",
             "Use 'summary' or 'manual' instead.")
      }
      sc <- compute_scaling_range(dat_target[, nonsurv_cols, drop = FALSE],
                                  dat_ref[, nonsurv_cols, drop = FALSE])
      tau_ns   <- sc$tau
      delta_ns <- sc$delta
      
    } else if (scaling_method == "summary") {
      # tau = 1, delta = shift to match target mean
      tau_ns   <- rep(1, length(nonsurv_cols))
      ref_means <- colMeans(dat_ref[, nonsurv_cols, drop = FALSE], na.rm = TRUE)
      delta_ns  <- target_means - ref_means
      names(tau_ns)   <- col_names[nonsurv_cols]
      names(delta_ns) <- col_names[nonsurv_cols]
      
    } else if (scaling_method == "manual") {
      stopifnot(!is.null(tau_manual), !is.null(delta_manual))
      tau_ns   <- tau_manual[nonsurv_cols]
      delta_ns <- delta_manual[nonsurv_cols]
      
    } else {
      stop("scaling_method must be 'range', 'summary', or 'manual'.")
    }
  } else {
    tau_ns   <- numeric(0)
    delta_ns <- numeric(0)
  }
  
  # =======================================================================
  # ALPHA / BETA for each non-survival column
  # =======================================================================
  alpha_vec <- numeric(p);  beta_vec <- numeric(p)
  tau_vec   <- numeric(p);  delta_vec <- numeric(p)
  type_results <- vector("list", p);  names(type_results) <- col_names
  
  search_defaults <- list(a1 = 0.2, a2 = 2.0, a_step = 0.01,
                          b1 = 0, b2 = 0.150, b_step = 0.001,
                          relax_index = 3,
                          search_index_alpha = 5)
  search_cfg <- modifyList(search_defaults, search_settings)
  
  for (idx in seq_along(nonsurv_cols)) {
    j <- nonsurv_cols[idx]
    col_j    <- dat_ref[, j]
    tar_mean <- target_means[idx]
    tar_sd   <- target_sds[idx]
    
    res_j <- find_alpha_beta(
      a1 = search_cfg$a1, a2 = search_cfg$a2, a_step = search_cfg$a_step,
      b1 = search_cfg$b1, b2 = search_cfg$b2, b_step = search_cfg$b_step,
      relax_index = search_cfg$relax_index,
      search_index_alpha = search_cfg$search_index_alpha,
      tau = tau_ns[idx], delta = delta_ns[idx],
      datref = col_j, tarmean = tar_mean, tarsd = tar_sd,
      type = types[j])
    
    alpha_vec[j] <- res_j$best_summary$alpha
    beta_vec[j]  <- res_j$best_summary$beta
    tau_vec[j]   <- tau_ns[idx]
    delta_vec[j] <- delta_ns[idx]
    type_results[[j]] <- res_j
    
    if (verbose) {
      cat(sprintf("  [%s] %-15s type=%-11s alpha=%.4f beta=%.4f | est_mean=%.3f (tar=%.3f) est_sd=%.3f (tar=%.3f)\n",
                  j, col_names[j], types[j],
                  res_j$best_summary$alpha, res_j$best_summary$beta,
                  res_j$best_summary$est_mean, tar_mean,
                  res_j$best_summary$est_sd, tar_sd))
    }
  }
  
  # =======================================================================
  # SURVIVAL columns
  # =======================================================================
  surv_pipeline <- NULL
  if (has_surv) {
    
    # --- Historical survival vectors ---
    if (!is.null(surv_time_col) && surv_time_col %in% col_names) {
      hist_time   <- dat_ref[[surv_time_col]]
      hist_status <- dat_ref[[surv_status_col]]
    } else {
      hist_time   <- dat_ref[, surv_cols[1]]
      hist_status <- if (length(surv_cols) >= 2) dat_ref[, surv_cols[2]] else rep(1, nrow(dat_ref))
    }
    
    # --- Determine target for survival: percentiles > raw vectors > dat_target ---
    if (!is.null(surv_target_percentiles)) {
      # Primary path: summary statistics
      surv_pipeline <- run_survival_sim(
        target_percentiles = surv_target_percentiles,
        hist_time = hist_time, hist_status = hist_status,
        N_sim = N_sim,
        use_km_scaling = surv_settings$use_km_scaling %||% FALSE,
        a1 = surv_settings$a1 %||% 0.2,
        a2 = surv_settings$a2 %||% 2.0,
        a_step = surv_settings$a_step %||% 0.01,
        b1 = surv_settings$b1 %||% 0,
        b2 = surv_settings$b2 %||% 0.200,
        b_step = surv_settings$b_step %||% 0.001,
        relax_index = surv_settings$relax_index %||% 5,
        search_index_alpha = surv_settings$search_index_alpha %||% 7,
        verbose = verbose)
      
    } else {
      # Backup path: raw target data
      if (!is.null(surv_target_time) && !is.null(surv_target_status)) {
        tgt_time   <- surv_target_time
        tgt_status <- surv_target_status
      } else if (!is.null(dat_target)) {
        if (!is.null(surv_time_col) && surv_time_col %in% colnames(dat_target)) {
          tgt_time   <- dat_target[[surv_time_col]]
          tgt_status <- dat_target[[surv_status_col]]
        } else {
          tgt_time   <- dat_target[, surv_cols[1]]
          tgt_status <- if (length(surv_cols) >= 2) {
            dat_target[, surv_cols[2]]
          } else {
            rep(1, nrow(dat_target))
          }
        }
      } else {
        stop("For survival columns, provide surv_target_percentiles, ",
             "surv_target_time/surv_target_status, or dat_target.")
      }
      
      surv_defaults <- list(n_quantiles = 5, margin_frac = 0.1, min_neff = 10,
                            use_km_scaling = FALSE, q_probs_manual = NULL,
                            a1 = 0.2, a2 = 2.0, a_step = 0.01,
                            b1 = 0, b2 = 0.200, b_step = 0.001,
                            relax_index = 5, search_index_alpha = 7)
      surv_cfg <- modifyList(surv_defaults, surv_settings)
      
      surv_pipeline <- run_survival_sim(
        target_time = tgt_time, target_status = tgt_status,
        hist_time = hist_time, hist_status = hist_status,
        N_sim = N_sim,
        n_quantiles = surv_cfg$n_quantiles, margin_frac = surv_cfg$margin_frac,
        min_neff = surv_cfg$min_neff, q_probs_manual = surv_cfg$q_probs_manual,
        use_km_scaling = surv_cfg$use_km_scaling,
        a1 = surv_cfg$a1, a2 = surv_cfg$a2, a_step = surv_cfg$a_step,
        b1 = surv_cfg$b1, b2 = surv_cfg$b2, b_step = surv_cfg$b_step,
        relax_index = surv_cfg$relax_index,
        search_index_alpha = surv_cfg$search_index_alpha,
        verbose = verbose)
    }
    
    # Fill in survival alpha/beta/tau/delta
    for (j in surv_cols) {
      alpha_vec[j] <- surv_pipeline$shape_params$best_summary$alpha
      beta_vec[j]  <- surv_pipeline$shape_params$best_summary$beta
      tau_vec[j]   <- surv_pipeline$scales$tau
      delta_vec[j] <- surv_pipeline$scales$delta
    }
  }
  
  # =======================================================================
  # GENERATE full multivariate data
  # =======================================================================
  names(alpha_vec) <- col_names;  names(beta_vec) <- col_names
  names(tau_vec)   <- col_names;  names(delta_vec) <- col_names
  
  sim_data <- generate_data(
    dat_ref = dat_ref, N = N_sim,
    alpha = alpha_vec, beta = beta_vec,
    tau = tau_vec, delta = delta_vec,
    types = types, cor_method = cor_method)
  
  if (verbose) {
    cat("\n--- Simulation Complete ---\n")
    cat("  N =", nrow(sim_data), " | columns =", ncol(sim_data), "\n\n")
  }
  
  list(sim_data = sim_data,
       alpha = alpha_vec, beta = beta_vec,
       tau = tau_vec, delta = delta_vec,
       types = types,
       target_means = if (length(nonsurv_cols) > 0) target_means else NULL,
       target_sds   = if (length(nonsurv_cols) > 0) target_sds   else NULL,
       type_results = type_results,
       surv_pipeline = surv_pipeline)
}
