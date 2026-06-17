#' Impute Missing Values in a Reference / Target Data Frame (Non-Survival)
#'
#' Fills \code{NA}/\code{NaN} entries in the non-survival columns of a data
#' frame using simple, widely used, type-aware rules, so that downstream
#' functions (\code{\link{find_alpha_beta}}, \code{\link{generate_data}},
#' \code{\link{run_simulation}}) — which all require complete columns —
#' can run.  Survival columns are NOT imputed here; censored survival times
#' are handled separately by \code{\link{impute_censored}}.
#'
#' @section Imputation rules (the most common defaults):
#'   \describe{
#'     \item{continuous}{\code{method = "mean"} (default) or \code{"median"}.
#'       Median is more robust to skew/outliers.}
#'     \item{binary}{Mode (most frequent of the two levels). Ties resolve to
#'       the larger code.}
#'     \item{ordinal}{Mode (most frequent observed level). Ties resolve to
#'       the larger code. Rounded to an observed category so no impossible
#'       levels are introduced.}
#'   }
#'   Survival columns (and any column whose name is in \code{skip_cols}) are
#'   passed through untouched.
#'
#' @param dat          A data frame.
#' @param types        Character vector of column types, length = \code{ncol(dat)}.
#'   One of \code{"continuous"}, \code{"binary"}, \code{"ordinal"},
#'   \code{"survival"}.
#' @param continuous_method  \code{"mean"} (default) or \code{"median"} for
#'   continuous columns.
#' @param skip_cols    Optional character vector of column names to leave
#'   untouched (e.g. survival time/status columns).
#' @param max_missing_frac  Numeric in \code{(0, 1]}. If the missing fraction
#'   in any imputable column exceeds this, a warning is emitted (the column is
#'   still imputed). Default \code{0.3}.
#' @param verbose      If TRUE, print a per-column summary of what was imputed.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{data}{The data frame with imputable \code{NA}/\code{NaN} filled.}
#'     \item{report}{A data.frame: column, type, n_missing, frac_missing,
#'       method, fill_value, imputed (logical), note.}
#'     \item{any_remaining}{Logical — TRUE if any \code{NA}/\code{NaN} remain
#'       (e.g. in skipped/survival columns or columns that could not be
#'       imputed).}
#'   }
#' @examples
#' df <- data.frame(age = c(50, NA, 60, NaN, 45),
#'                  sex = c(1, 0, NA, 1, 0),
#'                  ecog = c(0, 2, 2, NA, 1))
#' out <- impute_missing(df, types = c("continuous", "binary", "ordinal"))
#' out$data
#' out$report
#' @export
impute_missing <- function(dat, types,
                           continuous_method = "mean",
                           skip_cols = NULL,
                           max_missing_frac = 0.3,
                           verbose = TRUE) {

  # --- Input validation ----------------------------------------------------
  if (!is.data.frame(dat)) stop("`dat` must be a data frame.")
  p <- ncol(dat)
  if (length(types) != p) {
    stop("`types` must have length = ncol(dat) (", p, ").")
  }
  continuous_method <- match.arg(continuous_method, c("mean", "median"))
  if (!is.numeric(max_missing_frac) ||
      max_missing_frac <= 0 || max_missing_frac > 1) {
    stop("`max_missing_frac` must be a single number in (0, 1].")
  }

  col_names <- colnames(dat)
  if (is.null(col_names)) col_names <- paste0("V", seq_len(p))
  n_row <- nrow(dat)

  # treat both NA and NaN as missing
  is_missing <- function(x) is.na(x)   # is.na() is TRUE for NaN as well

  # mode for discrete codes (ties -> larger code)
  discrete_mode <- function(v) {
    v <- v[!is_missing(v)]
    if (length(v) == 0) return(NA_real_)
    tb <- table(v)
    top <- names(tb)[tb == max(tb)]
    max(as.numeric(top))
  }

  report <- data.frame(
    column       = col_names,
    type         = types,
    n_missing    = NA_integer_,
    frac_missing = NA_real_,
    method       = NA_character_,
    fill_value   = NA_real_,
    imputed      = FALSE,
    note         = NA_character_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(p)) {
    col      <- dat[[j]]
    miss_idx <- which(is_missing(col))
    n_miss   <- length(miss_idx)
    report$n_missing[j]    <- n_miss
    report$frac_missing[j] <- if (n_row > 0) n_miss / n_row else 0

    # --- Skip survival and user-skipped columns ----------------------------
    if (types[j] == "survival" || col_names[j] %in% skip_cols) {
      report$method[j] <- "skipped"
      report$note[j]   <- if (types[j] == "survival") {
        "survival column - use impute_censored() for censored times"
      } else {
        "in skip_cols"
      }
      next
    }

    if (n_miss == 0) {
      report$method[j] <- "none"
      report$note[j]   <- "no missing values"
      next
    }

    # --- Warn on heavy missingness -----------------------------------------
    if (report$frac_missing[j] > max_missing_frac) {
      warning(sprintf(
        "Column '%s' is %.1f%% missing (> %.0f%% threshold); imputation may be unreliable.",
        col_names[j], 100 * report$frac_missing[j], 100 * max_missing_frac))
    }

    # --- All values missing -> cannot impute -------------------------------
    if (n_miss == n_row) {
      report$method[j] <- "failed"
      report$note[j]   <- "all values missing - cannot impute"
      warning(sprintf("Column '%s' is entirely missing; cannot impute. Left as-is.",
                      col_names[j]))
      next
    }

    # --- Type-aware fill ----------------------------------------------------
    if (types[j] == "continuous") {
      num <- as.numeric(col)
      fill <- if (continuous_method == "median") {
        stats::median(num, na.rm = TRUE)
      } else {
        mean(num, na.rm = TRUE)
      }
      meth <- continuous_method

    } else if (types[j] %in% c("binary", "ordinal")) {
      fill <- discrete_mode(col)
      meth <- "mode"

    } else {
      report$method[j] <- "failed"
      report$note[j]   <- paste0("unknown type '", types[j], "'")
      warning(sprintf("Column '%s' has unknown type '%s'; not imputed.",
                      col_names[j], types[j]))
      next
    }

    if (is.na(fill) || (is.numeric(fill) && !is.finite(fill))) {
      report$method[j] <- "failed"
      report$note[j]   <- "fill value non-finite (insufficient observed data)"
      warning(sprintf("Column '%s': could not compute a finite fill value; left as-is.",
                      col_names[j]))
      next
    }

    col[miss_idx]        <- fill
    dat[[j]]             <- col
    report$method[j]     <- meth
    report$fill_value[j] <- as.numeric(fill)
    report$imputed[j]    <- TRUE
    report$note[j]       <- sprintf("filled %d value(s)", n_miss)
  }

  any_remaining <- any(vapply(dat, function(x) any(is_missing(x)), logical(1)))

  # --- Summary -------------------------------------------------------------
  total_missing <- sum(report$n_missing, na.rm = TRUE)
  if (verbose) {
    cat("--- Missing-Data Imputation (non-survival) ---\n")
    cat("  Rows:", n_row, "| Columns:", p, "\n")
    cat("  Total missing cells:", total_missing, "\n")
    print(report[, c("column", "type", "n_missing", "frac_missing",
                     "method", "fill_value")], row.names = FALSE)
    if (any_remaining) {
      cat("  NOTE: missing values remain (skipped/survival/failed columns).\n")
    }
    cat("\n")
  }

  if (total_missing == 0 && verbose) {
    message("No NA/NaN detected; data returned unchanged.")
  }
  if (any_remaining) {
    warning("Imputation finished but some NA/NaN values remain. ",
            "Downstream functions (find_alpha_beta, generate_data, ",
            "run_simulation) require complete non-survival columns.")
  }

  list(data = dat, report = report, any_remaining = any_remaining)
}