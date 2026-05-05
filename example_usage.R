################################################################################################
####  example_usage.R
####  Example usage for empirical_simulation_survival.R
####  Demonstrates all input paths:
####    1. Survival-only with pre-specified percentiles (no patient data)
####    2. Survival-only with raw target data (backup path)
####    3. Mixed-type with summary statistics (primary interface)
####    4. Mixed-type with dat_target backup
####    5. Non-survival only
####
####  Last update: 05/05/2026
################################################################################################

source("empirical_simulation_survival.R")

set.seed(42)

# =============================================================================
# GENERATE EXAMPLE HISTORICAL (REFERENCE) DATA
# =============================================================================
# In practice, this is your historical control arm or pooled reference dataset.

n_hist <- 300

hist_time   <- rweibull(n_hist, shape = 1.2, scale = 15)
hist_status <- rbinom(n_hist, 1, prob = 0.75)   # 25% censored
hist_age    <- rnorm(n_hist, mean = 60, sd = 10)
hist_sex    <- rbinom(n_hist, 1, prob = 0.45)    # binary: 1 = female
hist_ecog   <- sample(0:2, n_hist, replace = TRUE, prob = c(0.4, 0.4, 0.2))

hist_df <- data.frame(
  time = hist_time,
  age  = hist_age,
  sex  = hist_sex,
  ecog = hist_ecog
)

cat("========================================================\n")
cat("  Historical data summary\n")
cat("========================================================\n")
cat("  N =", n_hist, "\n")
cat("  Time:  mean =", round(mean(hist_time), 2),
    " sd =", round(sd(hist_time), 2), "\n")
cat("  Age:   mean =", round(mean(hist_age), 2),
    " sd =", round(sd(hist_age), 2), "\n")
cat("  Sex:   prop female =", round(mean(hist_sex), 3), "\n")
cat("  ECOG:  ", paste(names(table(hist_ecog)),
                       round(prop.table(table(hist_ecog)), 3),
                       sep = "=", collapse = ", "), "\n\n")


# =============================================================================
# EXAMPLE 1: Survival-only — pre-specified percentiles (PRIMARY PATH)
# =============================================================================
# This is the most common real-world scenario: you have published KM
# percentiles from a paper or protocol assumptions, but no patient data.

cat("========================================================\n")
cat("  EXAMPLE 1: Survival with target percentiles\n")
cat("========================================================\n\n")

# Target percentiles — e.g. from a published KM curve or protocol assumption
target_pctl <- data.frame(
  q_level  = c(0.10, 0.30, 0.50, 0.70, 0.90),
  time_val = c(3.5,  7.2,  12.0, 18.5, 32.0)
)

cat("Target percentiles:\n")
print(target_pctl)
cat("\n")

# Impute censored observations in historical data first
imp <- impute_censored(hist_time, hist_status, seed = 123, verbose = TRUE)

# Run survival simulation using percentiles (no target patient data needed)
result_1 <- run_survival_sim(
  target_percentiles = target_pctl,
  hist_time   = imp$imputed_time,
  hist_status = rep(1, length(imp$imputed_time)),  # fully imputed
  N_sim       = 5000,
  verbose     = TRUE
)

cat("Simulated quantiles vs targets:\n")
sim_q <- quantile(result_1$sim_data$sim_time, probs = target_pctl$q_level)
comparison_1 <- data.frame(
  q_level    = target_pctl$q_level,
  target     = target_pctl$time_val,
  simulated  = round(as.numeric(sim_q), 2)
)
print(comparison_1)
cat("\n\n")


# =============================================================================
# EXAMPLE 2: Survival-only — raw target data (BACKUP PATH)
# =============================================================================
# When you do have patient-level target data (e.g. from a completed trial).

cat("========================================================\n")
cat("  EXAMPLE 2: Survival with raw target data\n")
cat("========================================================\n\n")

# Simulate some "target arm" data for demonstration
n_target <- 200
target_time   <- rweibull(n_target, shape = 1.0, scale = 12)
target_status <- rbinom(n_target, 1, prob = 0.80)

result_2 <- run_survival_sim(
  target_time   = target_time,
  target_status = target_status,
  hist_time     = imp$imputed_time,
  hist_status   = rep(1, length(imp$imputed_time)),
  N_sim         = 5000,
  verbose       = TRUE
)

cat("Simulated quantiles:\n")
print(round(quantile(result_2$sim_data$sim_time,
                     c(0.10, 0.25, 0.50, 0.75, 0.90)), 2))
cat("\n\n")


# =============================================================================
# EXAMPLE 3: Mixed-type — summary statistics (PRIMARY PATH)
# =============================================================================
# Survival column uses percentiles; non-survival columns use mean/SD.
# No target patient data needed at all.

cat("========================================================\n")
cat("  EXAMPLE 3: Mixed-type with summary statistics\n")
cat("========================================================\n\n")

# Column types for hist_df: time=survival, age=continuous, sex=binary, ecog=ordinal
col_types <- c("survival", "continuous", "binary", "ordinal")

# Target summaries for non-survival columns (from literature, protocol, etc.)
# Order matches non-survival columns: age, sex, ecog
tgt_means <- c(age = 55.0, sex = 0.50, ecog = 0.90)
tgt_sds   <- c(age = 8.0,  sex = 0.50, ecog = 0.70)
tgt_mins  <- c(age = 25.0, sex = 0.0,  ecog = 0.0)
tgt_maxs  <- c(age = 85.0, sex = 1.0,  ecog = 2.0)

# Target survival percentiles (same as Example 1)
surv_pctl <- data.frame(
  q_level  = c(0.10, 0.30, 0.50, 0.70, 0.90),
  time_val = c(3.5,  7.2,  12.0, 18.5, 32.0)
)

cat("Target means:", tgt_means, "\n")
cat("Target SDs:  ", tgt_sds, "\n")
cat("Target mins: ", tgt_mins, "\n")
cat("Target maxs: ", tgt_maxs, "\n")
cat("Survival percentiles:\n")
print(surv_pctl)
cat("\n")

# Use imputed historical time in dat_ref
hist_df_imp <- hist_df
hist_df_imp$time <- imp$imputed_time

result_3 <- run_simulation(
  dat_ref       = hist_df_imp,
  types         = col_types,
  N_sim         = 5000,
  target_means  = tgt_means,
  target_sds    = tgt_sds,
  target_mins   = tgt_mins,
  target_maxs   = tgt_maxs,
  scaling_method = "range",
  surv_target_percentiles = surv_pctl,
  surv_time_col    = "time",
  surv_status_col  = NULL,            # fully imputed, not needed
  verbose       = TRUE
)

cat("--- Simulated data summary ---\n")
cat("  Time median: ", round(median(result_3$sim_data$time), 2), "\n")
cat("  Age  mean:   ", round(mean(result_3$sim_data$age), 2),
    " (target:", tgt_means["age"], ")\n")
cat("  Age  sd:     ", round(sd(result_3$sim_data$age), 2),
    " (target:", tgt_sds["age"], ")\n")
cat("  Sex  prop:   ", round(mean(result_3$sim_data$sex), 3),
    " (target:", tgt_means["sex"], ")\n")
cat("  ECOG mean:   ", round(mean(result_3$sim_data$ecog), 2),
    " (target:", tgt_means["ecog"], ")\n\n\n")


# =============================================================================
# EXAMPLE 4: Mixed-type — dat_target backup path
# =============================================================================
# When you have a full target dataset and want automatic extraction.

cat("========================================================\n")
cat("  EXAMPLE 4: Mixed-type with dat_target (backup)\n")
cat("========================================================\n\n")

# Create a "target" dataset for demonstration
dat_target_demo <- data.frame(
  time = rweibull(n_target, shape = 1.0, scale = 12),
  age  = rnorm(n_target, mean = 55, sd = 8),
  sex  = rbinom(n_target, 1, prob = 0.50),
  ecog = sample(0:2, n_target, replace = TRUE, prob = c(0.5, 0.35, 0.15))
)

result_4 <- run_simulation(
  dat_ref        = hist_df_imp,
  types          = col_types,
  N_sim          = 5000,
  dat_target     = dat_target_demo,        # auto-extracts means/SDs
  scaling_method = "summary",
  surv_target_time   = dat_target_demo$time,
  surv_target_status = rep(1, n_target),   # assume all events for demo
  surv_time_col      = "time",
  verbose        = TRUE
)

cat("--- Simulated vs target ---\n")
cat("  Age  sim mean:", round(mean(result_4$sim_data$age), 2),
    " target mean:", round(mean(dat_target_demo$age), 2), "\n")
cat("  Sex  sim prop:", round(mean(result_4$sim_data$sex), 3),
    " target prop:", round(mean(dat_target_demo$sex), 3), "\n\n\n")


# =============================================================================
# EXAMPLE 5: Non-survival only — continuous + binary + ordinal
# =============================================================================
# No survival column at all.

cat("========================================================\n")
cat("  EXAMPLE 5: Non-survival only\n")
cat("========================================================\n\n")

hist_nonsurv <- data.frame(
  age  = hist_age,
  sex  = hist_sex,
  ecog = hist_ecog
)

col_types_ns <- c("continuous", "binary", "ordinal")

result_5 <- run_simulation(
  dat_ref        = hist_nonsurv,
  types          = col_types_ns,
  N_sim          = 3000,
  target_means   = c(age = 58.0, sex = 0.55, ecog = 0.80),
  target_sds     = c(age = 9.0,  sex = 0.50, ecog = 0.65),
  target_mins    = c(age = 30.0, sex = 0.0,  ecog = 0.0),
  target_maxs    = c(age = 88.0, sex = 1.0,  ecog = 2.0),
  scaling_method = "range",
  verbose        = TRUE
)

cat("--- Simulated data ---\n")
cat("  Age  mean:", round(mean(result_5$sim_data$age), 2), " sd:",
    round(sd(result_5$sim_data$age), 2), "\n")
cat("  Sex  prop:", round(mean(result_5$sim_data$sex), 3), "\n")
cat("  ECOG mean:", round(mean(result_5$sim_data$ecog), 2), "\n\n")


# =============================================================================
# EXAMPLE 6: Using extract_target_summaries() as a helper
# =============================================================================
# Shows how to extract summaries from data for reuse across simulations.

cat("========================================================\n")
cat("  EXAMPLE 6: extract_target_summaries() helper\n")
cat("========================================================\n\n")

summ <- extract_target_summaries(dat_target_demo[, c("age", "sex", "ecog")],
                                 types = c("continuous", "binary", "ordinal"))
cat("Extracted from dat_target_demo:\n")
cat("  means:", summ$means, "\n")
cat("  sds:  ", summ$sds, "\n")
cat("  mins: ", summ$mins, "\n")
cat("  maxs: ", summ$maxs, "\n\n")

# Now use these in a new simulation
result_6 <- run_simulation(
  dat_ref        = hist_nonsurv,
  types          = col_types_ns,
  N_sim          = 3000,
  target_means   = summ$means,
  target_sds     = summ$sds,
  target_mins    = summ$mins,
  target_maxs    = summ$maxs,
  scaling_method = "range",
  verbose        = FALSE
)

cat("Simulation with extracted summaries:\n")
cat("  Age  mean:", round(mean(result_6$sim_data$age), 2),
    " (target:", round(summ$means["age"], 2), ")\n")
cat("  Sex  prop:", round(mean(result_6$sim_data$sex), 3),
    " (target:", round(summ$means["sex"], 3), ")\n")
cat("  ECOG mean:", round(mean(result_6$sim_data$ecog), 2),
    " (target:", round(summ$means["ecog"], 2), ")\n\n")


# =============================================================================
# EXAMPLE 7: assign_roles_from_percentiles() — inspect role assignment
# =============================================================================
# Shows how percentiles are mapped to pipeline roles.

cat("========================================================\n")
cat("  EXAMPLE 7: Role assignment from percentiles\n")
cat("========================================================\n\n")

# 5 percentiles -> SL=P10, low=P30, med=P50, high=P70, SH=P90
pctl_5 <- data.frame(
  q_level  = c(0.10, 0.30, 0.50, 0.70, 0.90),
  time_val = c(3.5,  7.2,  12.0, 18.5, 32.0)
)
roles_5 <- assign_roles_from_percentiles(pctl_5)
cat("5 percentiles:\n")
cat("  Scaling:  q_SL=", roles_5$q_SL, "->", roles_5$val_SL,
    "  q_SH=", roles_5$q_SH, "->", roles_5$val_SH, "\n")
cat("  Median:   q_med=", roles_5$q_med, "->", roles_5$val_med, "\n")
cat("  Shape:    q_low=", roles_5$q_low, "->", roles_5$val_low,
    "  q_high=", roles_5$q_high, "->", roles_5$val_high, "\n\n")

# 3 percentiles (minimum) -> SL=P25, med=P50, SH=P75
pctl_3 <- data.frame(
  q_level  = c(0.25, 0.50, 0.75),
  time_val = c(6.0,  12.0, 20.0)
)
roles_3 <- assign_roles_from_percentiles(pctl_3)
cat("3 percentiles (minimum):\n")
cat("  Scaling:  q_SL=", roles_3$q_SL, "->", roles_3$val_SL,
    "  q_SH=", roles_3$q_SH, "->", roles_3$val_SH, "\n")
cat("  Median:   q_med=", roles_3$q_med, "->", roles_3$val_med, "\n")
cat("  Shape:    q_low=", roles_3$q_low, "->", roles_3$val_low,
    "  q_high=", roles_3$q_high, "->", roles_3$val_high, "\n\n")

cat("========================================================\n")
cat("  All examples complete.\n")
cat("========================================================\n")
