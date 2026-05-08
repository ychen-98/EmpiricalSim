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
#' @return A list with Stage 1/2 grids, best parameters, and diagnostics.
#' @export
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
        stats::quantile(col, probs = c(q_med, q_low, q_high), na.rm = TRUE, names = FALSE)
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
