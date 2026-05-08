test_that("calc_neff returns expected columns", {
  fit <- survival::survfit(survival::Surv(obs_time, status) ~ 1)
  res <- calc_neff(fit)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("time", "n_risk", "n_event", "surv", "neff", "arm") %in% names(res)))
  expect_true(all(res$neff[!is.na(res$neff)] > 0))
})

test_that("extract_km_percentiles returns correct structure", {
  res <- extract_km_percentiles(obs_time, status)
  expect_named(res, c("percentiles", "stats", "roles", "neff_table"))
  expect_s3_class(res$percentiles, "data.frame")
  expect_true(all(c("q_level", "time_val", "reachable") %in% names(res$percentiles)))
})

test_that("assign_roles_from_percentiles assigns 5 roles correctly", {
  tp <- data.frame(
    q_level  = c(0.10, 0.30, 0.50, 0.70, 0.90),
    time_val = c(3.5,  7.2,  12.0, 18.5, 32.0)
  )
  r <- assign_roles_from_percentiles(tp)
  expect_equal(r$q_SL,  0.10); expect_equal(r$val_SL,  3.5)
  expect_equal(r$q_med, 0.50); expect_equal(r$val_med, 12.0)
  expect_equal(r$q_SH,  0.90); expect_equal(r$val_SH,  32.0)
})

test_that("assign_roles_from_percentiles requires >= 3 rows", {
  tp2 <- data.frame(q_level = c(0.25, 0.75), time_val = c(6, 20))
  expect_error(assign_roles_from_percentiles(tp2), "at least 3")
})

test_that("compute_scaling_range maps reference range to target range", {
  ref <- data.frame(x = 1:10, y = 1:10)
  sc  <- compute_scaling_range(c(0, 5), c(10, 20), ref)
  expect_equal(sc$tau["x"],   (10 - 0) / (10 - 1))
  expect_equal(length(sc$tau), 2L)
})

test_that("find_alpha_beta returns valid output for continuous type", {
  set.seed(42)
  ref <- rnorm(100, 50, 10)
  res <- find_alpha_beta(tau = 1, delta = 0, datref = ref,
                         tarmean = 55, tarsd = 8, type = "continuous")
  expect_true(!is.null(res$best_summary$alpha))
  expect_true(res$best_summary$alpha > 0)
})

test_that("generate_data returns correct dimensions", {
  set.seed(7)
  ref <- data.frame(x = rnorm(50, 60, 10), y = rbinom(50, 1, 0.4))
  out <- generate_data(ref, N = 200,
                       alpha = c(1, 1), beta = c(0, 0),
                       tau   = c(1, 1), delta = c(0, 0),
                       types = c("continuous", "binary"))
  expect_equal(nrow(out), 200L)
  expect_equal(ncol(out), 2L)
  expect_true(all(out$y %in% c(0, 1)))
})

test_that("run_survival_sim (percentile path) returns sim_data", {
  tp <- data.frame(
    q_level  = c(0.25, 0.50, 0.75),
    time_val = c(6.0,  12.0, 20.0)
  )
  imp <- impute_censored(obs_time, status, seed = 1, verbose = FALSE)
  res <- run_survival_sim(
    target_percentiles = tp,
    hist_time   = imp$imputed_time,
    hist_status = rep(1L, length(imp$imputed_time)),
    N_sim = 500, verbose = FALSE
  )
  expect_s3_class(res$sim_data, "data.frame")
  expect_equal(nrow(res$sim_data), 500L)
  expect_true(all(res$sim_data$sim_time > 0))
})

test_that("run_simulation (non-survival) returns correct columns", {
  ref_ns <- hist_df[, c("age", "sex", "ecog")]
  res <- run_simulation(
    dat_ref      = ref_ns,
    types        = c("continuous", "binary", "ordinal"),
    N_sim        = 300,
    target_means = c(55, 0.50, 0.90),
    target_sds   = c(8,  0.50, 0.70),
    target_mins  = c(25, 0.0,  0.0),
    target_maxs  = c(85, 1.0,  2.0),
    verbose      = FALSE
  )
  expect_equal(nrow(res$sim_data), 300L)
  expect_named(res$sim_data, c("age", "sex", "ecog"))
  expect_true(all(res$sim_data$sex %in% c(0, 1)))
})

test_that("extract_target_summaries returns correct structure", {
  summ <- extract_target_summaries(hist_df[, c("age", "sex", "ecog")],
                                   types = c("continuous", "binary", "ordinal"))
  expect_named(summ, c("means", "sds", "mins", "maxs"))
  expect_length(summ$means, 3L)
})
