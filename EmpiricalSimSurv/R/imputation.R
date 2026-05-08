#' Compute a Tail-Weighted AIC for a flexsurvreg Fit
#'
#' @param fit         A \code{flexsurvreg} object, or NULL.
#' @param obs_time    Numeric vector of observed times.
#' @param status      Event indicator (1 = event, 0 = censored).
#' @param tail_start  Time at which tail weighting begins.
#' @param multiplier  Weight multiplier for tail observations (default 3).
#' @return Numeric scalar — weighted AIC (or standard AIC if tail_start <= 0).
#' @export
tail_weighted_aic <- function(fit, obs_time, status,
                              tail_start, multiplier = 3) {
  if (is.null(fit)) return(Inf)
  if (tail_start <= 0 || multiplier <= 1) return(stats::AIC(fit))

  pars <- fit$res[, "est"]
  dfn  <- fit$dfns$d
  pfn  <- fit$dfns$p

  n    <- length(obs_time)
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
#' Uses a KM curve for early censoring and falls back to a parametric model
#' (selected by tail-weighted AIC) for late censoring. The truncation point
#' is determined by Peto's effective sample size (N_eff).
#'
#' @param obs_time   Numeric vector of observed times.
#' @param status     Event indicator (1 = event, 0 = censored).
#' @param seed       Random seed for reproducibility.
#' @param min_neff   Minimum Peto's N_eff required to trust the KM estimate.
#' @param verbose    If TRUE, print diagnostics.
#' @param candidate_dists Character vector of flexsurv distribution names.
#' @param tail_weight_fraction Fraction of trunc_time for tail weighting.
#' @param tail_weight_multiplier Extra weight for tail observations.
#' @param blend_window Width of blending window around trunc_time.
#' @param M Number of multiply-imputed datasets to return.
#' @param s_param_floor Minimum S_param(c) for inverse-CDF sampling.
#'
#' @return If M == 1, a list with imputed data and diagnostics.
#'   If M > 1, a list with \code{imputations} and \code{shared} info.
#' @export
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
  fit_km  <- survival::survfit(survival::Surv(obs_time, status) ~ 1)
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
  surv_obj <- survival::Surv(obs_time, status)

  fits <- setNames(
    lapply(candidate_dists, function(d) {
      tryCatch(flexsurv::flexsurvreg(surv_obj ~ 1, dist = d),
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

  standard_aic <- sapply(fits, function(f) tryCatch(stats::AIC(f), error = function(e) Inf))

  aic_table <- data.frame(
    distribution      = names(fits),
    AIC_standard      = standard_aic,
    AIC_tail_weighted = aic_values,
    stringsAsFactors  = FALSE
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
      use_km_flag <- (c_time < blend_lo) ||
        (blend_window == 0 && c_time < trunc_time)

      draw_km <- function() {
        s_km_at_c <- stats::approx(km_times, km_surv, xout = c_time,
                                   method = "constant", rule = 2, f = 0)$y
        u <- stats::runif(1, 0, s_km_at_c)
        if (u < km_min_surv) {
          return(draw_param())
        }
        t_new <- stats::approx(km_surv, km_times, xout = u,
                               method = "constant", rule = 2, ties = max)$y
        max(t_new, c_time + 0.001)
      }

      draw_param <- function() {
        s_at_c <- param_surv_at(c_time)
        if (s_at_c < s_param_floor) {
          t_new <- param_conditional_median(c_time)
          return(max(t_new, c_time + 0.001))
        }
        u     <- stats::runif(1, 0, s_at_c)
        t_new <- param_quantile_at(u)
        if (is.na(t_new) || !is.finite(t_new)) {
          t_new <- param_conditional_median(c_time)
        }
        if (is.na(t_new) || !is.finite(t_new)) {
          t_new <- c_time * 1.5
        }
        max(t_new, c_time + 0.001)
      }

      if (use_km_flag) {
        imputed_time[i] <- draw_km()
        method_used[i]  <- "km"
        km_count <- km_count + 1

      } else if (in_blend_zone) {
        alpha_blend <- (blend_hi - c_time) / blend_window
        alpha_blend <- max(0, min(1, alpha_blend))

        t_km    <- draw_km()
        t_param <- draw_param()
        imputed_time[i] <- alpha_blend * t_km + (1 - alpha_blend) * t_param
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
