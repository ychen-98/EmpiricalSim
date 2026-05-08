#' Extract Target Summary Statistics from a Data Frame
#'
#' Convenience function to compute per-column mean, SD, min, and max
#' from a target data frame, for use as \code{target_means},
#' \code{target_sds}, \code{target_mins}, \code{target_maxs} in
#' \code{run_simulation()}.  Binary columns with zero SD are replaced
#' with the Bernoulli SD.
#'
#' @param dat_target  Target data frame.
#' @param types       Character vector of column types.
#' @return List with vectors \code{means}, \code{sds}, \code{mins},
#'   \code{maxs}.
#' @export
extract_target_summaries <- function(dat_target, types) {
  p <- ncol(dat_target)
  means <- sds <- mins <- maxs <- numeric(p)
  for (j in seq_len(p)) {
    means[j] <- mean(dat_target[, j], na.rm = TRUE)
    sds[j]   <- stats::sd(dat_target[, j], na.rm = TRUE)
    mins[j]  <- min(dat_target[, j], na.rm = TRUE)
    maxs[j]  <- max(dat_target[, j], na.rm = TRUE)
    if (types[j] == "binary" && (is.na(sds[j]) || sds[j] == 0)) {
      sds[j] <- sqrt(means[j] * (1 - means[j]))
    }
  }
  names(means) <- names(sds) <- names(mins) <- names(maxs) <- colnames(dat_target)
  list(means = means, sds = sds, mins = mins, maxs = maxs)
}
