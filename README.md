# Empirical Simulation Framework

A unified R framework for generating simulated clinical trial data that matches target summary statistics, without requiring patient-level target data. Supports continuous, binary, ordinal, and survival endpoints — individually or in any combination.

**Reference:** Ding Y, Liu Y, Qu Y (2025), *Commun Stat Simul Comput*.

---

## Overview

The framework maps a historical reference distribution to a target distribution using a power-distortion transformation on the probability integral. Rather than requiring individual patient data from the target population, you specify target characteristics as summary statistics:

- **Survival endpoints:** percentile–time pairs from a published Kaplan-Meier curve (e.g., "median OS = 12 months, P30 = 7.2 months")
- **Continuous/binary/ordinal endpoints:** target mean and standard deviation

A two-stage adaptive grid search finds the optimal distortion parameters (alpha, beta) that best reproduce these targets from the historical data.

## Requirements

- R >= 4.1.2
- Required packages: `dplyr`, `MASS`, `survival`
- Optional (for censoring imputation): `flexsurv`

## Files

| File | Description |
|------|-------------|
| `empirical_simulation_survival.R` | Core framework — all functions |
| `example_usage.R` | Seven worked examples covering all input paths |

## Quick Start

### Survival data with published percentiles (most common use case)

```r
source("empirical_simulation_survival.R")

# Target from a published KM curve — no patient data needed
target_pctl <- data.frame(
  q_level  = c(0.10, 0.30, 0.50, 0.70, 0.90),
  time_val = c(3.5,  7.2,  12.0, 18.5, 32.0)
)

# Historical reference arm (your available data)
# hist_time and hist_status are your reference vectors

# Step 1: Impute censored observations
imp <- impute_censored(hist_time, hist_status, seed = 123)

# Step 2: Simulate
result <- run_survival_sim(
  target_percentiles = target_pctl,
  hist_time   = imp$imputed_time,
  hist_status = rep(1, length(imp$imputed_time)),
  N_sim       = 5000
)

# Simulated data
head(result$sim_data)
```

### Mixed-type data with summary statistics

```r
# Historical reference data frame
# Columns: time (survival), age (continuous), sex (binary), ecog (ordinal)
col_types <- c("survival", "continuous", "binary", "ordinal")

# Target summaries — one value per non-survival column
tgt_means <- c(age = 55.0, sex = 0.50, ecog = 0.90)
tgt_sds   <- c(age = 8.0,  sex = 0.50, ecog = 0.70)
tgt_mins  <- c(age = 25.0, sex = 0.0,  ecog = 0.0)
tgt_maxs  <- c(age = 85.0, sex = 1.0,  ecog = 2.0)

# Target survival percentiles
surv_pctl <- data.frame(
  q_level  = c(0.10, 0.30, 0.50, 0.70, 0.90),
  time_val = c(3.5,  7.2,  12.0, 18.5, 32.0)
)

result <- run_simulation(
  dat_ref       = hist_df,
  types         = col_types,
  N_sim         = 5000,
  target_means  = tgt_means,
  target_sds    = tgt_sds,
  target_mins   = tgt_mins,
  target_maxs   = tgt_maxs,
  scaling_method = "range",
  surv_target_percentiles = surv_pctl,
  surv_time_col = "time"
)
```

### Non-survival data only

```r
result <- run_simulation(
  dat_ref      = hist_df,
  types        = c("continuous", "binary", "ordinal"),
  N_sim        = 3000,
  target_means = c(age = 58.0, sex = 0.55, ecog = 0.80),
  target_sds   = c(age = 9.0,  sex = 0.50, ecog = 0.65),
  target_mins  = c(age = 30.0, sex = 0.0,  ecog = 0.0),
  target_maxs  = c(age = 88.0, sex = 1.0,  ecog = 2.0)
)
```

## Input Paths

The framework supports two input paths. The primary path uses summary statistics only; the backup path accepts raw patient data for convenience.

### Primary: summary statistics (no target patient data)

| Endpoint type | What you provide |
|---------------|-----------------|
| Survival | `target_percentiles`: a data.frame with `q_level` and `time_val` columns (minimum 3 rows) |
| Continuous | `target_means`, `target_sds`, `target_mins`, `target_maxs` (one element per non-survival column) |
| Binary | `target_means` (proportion), `target_sds`, `target_mins`, `target_maxs` |
| Ordinal | `target_means`, `target_sds`, `target_mins`, `target_maxs` |

### Backup: raw target data

Supply `dat_target` (a data.frame) and/or `surv_target_time` + `surv_target_status`. The framework auto-extracts means, SDs, and KM percentiles. Use `extract_target_summaries()` to inspect what gets extracted.

## Function Reference

### Core pipeline

| Function | Purpose |
|----------|---------|
| `run_simulation()` | Main entry point — mixed data types, any combination |
| `run_survival_sim()` | Survival-only pipeline (called internally or standalone) |
| `find_alpha_beta()` | Two-stage grid search for optimal distortion parameters |
| `generate_data()` | Multivariate correlated data generation |

### Scaling

| Function | Purpose |
|----------|---------|
| `compute_scaling()` | Log-linear tau/delta for survival data |
| `compute_scaling_range()` | Range-based tau/delta for non-survival data |

### Survival utilities

| Function | Purpose |
|----------|---------|
| `extract_km_percentiles()` | KM fit → percentile table with auto role assignment |
| `assign_roles_from_percentiles()` | Role assignment from user-supplied percentile table |
| `auto_quantile_levels()` | N_eff-based safe quantile grid |
| `calc_neff()` | Peto's effective sample size from a survfit object |
| `impute_censored()` | KM + parametric tail imputation for censored data |

### Helpers

| Function | Purpose |
|----------|---------|
| `extract_target_summaries()` | Extract mean/SD from a target data.frame (backup utility) |

## Key Concepts

### Percentile role assignment

When you supply target percentiles, the framework assigns five roles used by the pipeline:

| Role | Purpose | Default with 5 percentiles (P10–P90) |
|------|---------|--------------------------------------|
| SL (scaling low) | Lower anchor for tau/delta | P10 |
| SH (scaling high) | Upper anchor for tau/delta | P90 |
| med (median) | Primary shape-matching target | P50 |
| low (shape low) | Secondary shape target | P30 |
| high (shape high) | Secondary shape target | P70 |

With 3 percentiles (minimum), SL/low share the lowest level and SH/high share the highest.

### Scaling methods for non-survival data

The `scaling_method` parameter in `run_simulation()` controls how tau and delta are computed:

| Method | Description | Required inputs |
|--------|-------------|----------------|
| `"range"` (default) | Maps reference min/max to target min/max | `target_mins`, `target_maxs` |
| `"manual"` | User-supplied `tau_manual` / `delta_manual` | `tau_manual`, `delta_manual` |

All methods also require `target_means` and `target_sds` for the alpha/beta shape-matching step. The scaling method only affects how tau and delta are derived.

### Censoring imputation

Before simulation, censored historical data should be imputed using `impute_censored()`. This function uses a KM curve for early censoring and falls back to a parametric model (selected by tail-weighted AIC) for late censoring. The truncation point is determined by Peto's effective sample size (N_eff), not a fixed sample-size threshold.

## Examples

`example_usage.R` contains seven self-contained examples:

| # | Scenario | Target input |
|---|----------|-------------|
| 1 | Survival, percentiles | `target_percentiles` — no patient data |
| 2 | Survival, raw data | `target_time` + `target_status` (backup) |
| 3 | Mixed-type, summaries | `target_means`/`target_sds` + `surv_target_percentiles` |
| 4 | Mixed-type, dat_target | `dat_target` auto-extraction (backup) |
| 5 | Non-survival only | `target_means`/`target_sds` |
| 6 | Helper utility | `extract_target_summaries()` for reuse |
| 7 | Role inspection | `assign_roles_from_percentiles()` output |

Run all examples:

```r
source("example_usage.R")
```
