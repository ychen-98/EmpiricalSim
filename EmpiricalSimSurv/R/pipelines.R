#' Run Full Empirical Simulation Pipeline (Survival Data)
#'
#' End-to-end pipeline for survival-only simulations.  Accepts EITHER:
#'   (a) \code{target_percentiles} — a data.frame(q_level, time_val) of
#'       pre-specified percentile targets (primary interface), OR
#'   (b) \code{target_time} / \code{target_status} — raw patient-level data
#'       from which KM percentiles are extracted (backup / convenience).
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
#' @return List: sim_data, scales, shape_params, target_km, hist_km,
#'   q_probs_used, q_probs_mode, min_neff_used.
#' @export
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

    r <- assign_roles_from_percentiles(tp)

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
  # Common path: historical KM, scaling, shape matching, simulation
  # ================================================================

  hist_km <- extract_km_percentiles(hist_time, hist_status, q_probs = q_probs)

  if (verbose) {
    cat("--- Historical KM ---\n"); print(hist_km$percentiles); cat("\n")
  }

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

  scales <- compute_scaling(
    hist_data = hist_time, target_val_low = r$val_SL, target_val_high = r$val_SH,
    q_scale_range = c(r$q_SL, r$q_SH),
    obs_time = hist_time, status = hist_status, use_km = use_km_scaling)

  if (verbose) cat("Scaling: tau =", scales$tau, ", delta =", scales$delta, "\n\n")

  shape <- find_alpha_beta(
    a1 = a1, a2 = a2, a_step = a_step, b1 = b1, b2 = b2, b_step = b_step,
    relax_index = relax_index, search_index_alpha = search_index_alpha,
    tau = scales$tau, delta = scales$delta, datref = hist_time,
    type = "survival",
    target_med = r$val_med, target_low = r$val_low, target_high = r$val_high,
    q_low = r$q_low, q_med = r$q_med, q_high = r$q_high)

  if (verbose) { cat("--- Shape ---\n"); print(shape$best_summary); cat("\n") }

  n0 <- length(hist_time); sorted_hist <- sort(hist_time)
  if (is.null(rho_matrix)) {
    Z <- matrix(stats::rnorm(N_sim), ncol = 1)
  } else {
    Z <- MASS::mvrnorm(N_sim, rep(0, ncol(rho_matrix)), rho_matrix)
  }
  U <- stats::pnorm(Z)
  a_opt <- shape$best_summary$alpha;  b_opt <- shape$best_summary$beta
  U_prime <- U^a_opt * (1 - U)^b_opt
  idx <- floor(n0 * U_prime); idx[idx < 1] <- 1; idx[idx > n0] <- n0
  sim_vec <- pmax(0.001, exp(log(sorted_hist[idx]) * scales$tau + scales$delta))

  dat_sim <- data.frame(sim_time = as.vector(sim_vec), sim_status = 1L)

  if (verbose) {
    cat("--- Simulation Complete (N =", nrow(dat_sim), ") ---\n")
    print(round(stats::quantile(dat_sim$sim_time, c(0.05, 0.25, 0.50, 0.75, 0.95)), 4))
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


#' Run Full Empirical Simulation (Mixed Data Types)
#'
#' Orchestrates the complete simulation for a data set that may contain
#' any mix of continuous, binary, ordinal, and survival columns.
#'
#' \strong{Primary interface (summary statistics):}
#'   For non-survival columns, supply \code{target_means},
#'   \code{target_sds}, \code{target_mins}, and \code{target_maxs}
#'   (numeric vectors, one element per non-survival column).
#'   For survival columns, supply \code{surv_target_percentiles}
#'   (a data.frame with \code{q_level} and \code{time_val}).
#'
#' \strong{Backup interface (raw data):}
#'   Optionally supply \code{dat_target} — a data frame from which
#'   target means, SDs, mins, maxs, and survival KM percentiles are
#'   extracted automatically via \code{extract_target_summaries()} and
#'   \code{extract_km_percentiles()}.
#'
#' @param dat_ref       Reference data frame (historical / global control).
#' @param types         Character vector of column types (length = ncol).
#' @param N_sim         Number of rows to simulate.
#' @param target_means  Numeric vector of target means for non-survival
#'   columns (length = number of non-survival columns).
#' @param target_sds    Numeric vector of target SDs for non-survival columns.
#' @param target_mins   Numeric vector of target minimums for non-survival
#'   columns (required for \code{scaling_method = "range"}).
#' @param target_maxs   Numeric vector of target maximums for non-survival
#'   columns (required for \code{scaling_method = "range"}).
#' @param dat_target    Optional backup: target data frame for automatic
#'   extraction of means, SDs, mins, maxs, and survival data.
#' @param scaling_method  \code{"range"} (default; uses target_mins/target_maxs), or
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
#' @export
run_simulation <- function(dat_ref, types, N_sim = 5000,
                           target_means = NULL, target_sds = NULL,
                           target_mins = NULL, target_maxs = NULL,
                           dat_target = NULL,
                           scaling_method = "range",
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
    needs_mean_sd  <- is.null(target_means) || is.null(target_sds)
    needs_min_max  <- (scaling_method == "range") &&
                      (is.null(target_mins) || is.null(target_maxs))

    if (needs_mean_sd || needs_min_max) {
      if (is.null(dat_target)) {
        if (needs_mean_sd) {
          stop("For non-survival columns, provide target_means/target_sds ",
               "or dat_target for automatic extraction.")
        }
        if (needs_min_max) {
          stop("scaling_method='range' requires target_mins/target_maxs ",
               "or dat_target for automatic extraction.")
        }
      }
      if (verbose) cat("Extracting target summaries from dat_target...\n")
      summ <- extract_target_summaries(dat_target[, nonsurv_cols, drop = FALSE],
                                       types[nonsurv_cols])
      if (needs_mean_sd) {
        target_means <- summ$means
        target_sds   <- summ$sds
      }
      if (needs_min_max) {
        target_mins <- summ$mins
        target_maxs <- summ$maxs
      }
    }

    if (length(target_means) != length(nonsurv_cols) ||
        length(target_sds)   != length(nonsurv_cols)) {
      stop("target_means and target_sds must have length = ",
           length(nonsurv_cols), " (one per non-survival column).")
    }
    if (scaling_method == "range") {
      if (length(target_mins) != length(nonsurv_cols) ||
          length(target_maxs) != length(nonsurv_cols)) {
        stop("target_mins and target_maxs must have length = ",
             length(nonsurv_cols), " (one per non-survival column).")
      }
    }

    names(target_means) <- col_names[nonsurv_cols]
    names(target_sds)   <- col_names[nonsurv_cols]
    if (!is.null(target_mins)) names(target_mins) <- col_names[nonsurv_cols]
    if (!is.null(target_maxs)) names(target_maxs) <- col_names[nonsurv_cols]
  }

  # =======================================================================
  # TAU / DELTA for non-survival columns
  # =======================================================================
  if (length(nonsurv_cols) > 0) {
    if (scaling_method == "range") {
      sc <- compute_scaling_range(target_mins, target_maxs,
                                  dat_ref[, nonsurv_cols, drop = FALSE],
                                  col_names = col_names[nonsurv_cols])
      tau_ns   <- sc$tau
      delta_ns <- sc$delta

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
  search_cfg <- utils::modifyList(search_defaults, search_settings)

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

    if (!is.null(surv_time_col) && surv_time_col %in% col_names) {
      hist_time   <- dat_ref[[surv_time_col]]
      hist_status <- dat_ref[[surv_status_col]]
    } else {
      hist_time   <- dat_ref[, surv_cols[1]]
      hist_status <- if (length(surv_cols) >= 2) dat_ref[, surv_cols[2]] else rep(1, nrow(dat_ref))
    }

    if (!is.null(surv_target_percentiles)) {
      surv_pipeline <- run_survival_sim(
        target_percentiles = surv_target_percentiles,
        hist_time = hist_time, hist_status = hist_status,
        N_sim = N_sim,
        use_km_scaling       = surv_settings$use_km_scaling       %||% FALSE,
        a1                   = surv_settings$a1                   %||% 0.2,
        a2                   = surv_settings$a2                   %||% 2.0,
        a_step               = surv_settings$a_step               %||% 0.01,
        b1                   = surv_settings$b1                   %||% 0,
        b2                   = surv_settings$b2                   %||% 0.200,
        b_step               = surv_settings$b_step               %||% 0.001,
        relax_index          = surv_settings$relax_index          %||% 5,
        search_index_alpha   = surv_settings$search_index_alpha   %||% 7,
        verbose = verbose)

    } else {
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
      surv_cfg <- utils::modifyList(surv_defaults, surv_settings)

      surv_pipeline <- run_survival_sim(
        target_time = tgt_time, target_status = tgt_status,
        hist_time = hist_time, hist_status = hist_status,
        N_sim = N_sim,
        n_quantiles        = surv_cfg$n_quantiles,
        margin_frac        = surv_cfg$margin_frac,
        min_neff           = surv_cfg$min_neff,
        q_probs_manual     = surv_cfg$q_probs_manual,
        use_km_scaling     = surv_cfg$use_km_scaling,
        a1                 = surv_cfg$a1,
        a2                 = surv_cfg$a2,
        a_step             = surv_cfg$a_step,
        b1                 = surv_cfg$b1,
        b2                 = surv_cfg$b2,
        b_step             = surv_cfg$b_step,
        relax_index        = surv_cfg$relax_index,
        search_index_alpha = surv_cfg$search_index_alpha,
        verbose = verbose)
    }

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
       target_mins  = if (length(nonsurv_cols) > 0) target_mins  else NULL,
       target_maxs  = if (length(nonsurv_cols) > 0) target_maxs  else NULL,
       type_results = type_results,
       surv_pipeline = surv_pipeline)
}
