#' Generate Simulated Data via Empirical Distortion
#'
#' @param dat_ref     Reference data frame (columns = variables).
#' @param N           Number of rows to simulate.
#' @param alpha,beta  Named or positional numeric vectors (length = ncol).
#' @param tau,delta   Named or positional numeric vectors (length = ncol).
#' @param types       Character vector of column types.
#' @param cor_method  \code{"rank"} (preserve rank correlations) or
#'   \code{"identity"} (independent columns).
#'
#' @return Data frame with N rows, same column names as \code{dat_ref}.
#' @export
generate_data <- function(dat_ref, N, alpha, beta, tau, delta,
                          types, cor_method = "rank") {

  p  <- ncol(dat_ref)
  n0 <- nrow(dat_ref)
  stopifnot(length(alpha) == p, length(beta) == p,
            length(tau) == p, length(delta) == p, length(types) == p)

  if (cor_method == "rank") {
    F_u     <- (apply(dat_ref, 2, rank) - 0.5) / n0
    Phi_inv <- stats::qnorm(F_u)
    Sigma   <- stats::cov(Phi_inv)
  } else {
    Sigma <- diag(p)
  }

  Z <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = Sigma, tol = 1e-6)
  U <- stats::pnorm(Z)

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
      tab   <- table(col_data)
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
