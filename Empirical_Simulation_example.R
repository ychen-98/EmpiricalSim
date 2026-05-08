################################################################################################
####  empirical_simulation_example.R
####  Uses the EmpiricalSimSurv package on example_new.csv (Non-survival data only)
####
####  Author: Emma "Yao" Chen  
####  R version >= 4.1.2 
####  Last update: 05/05/2026
####  Reference: Ding Y, Liu Y, Qu Y (2025), Commun Stat Simul Comput.
################################################################################################



#---------------------------------------------------------------------#
#-------------------------- Library Calls ----------------------------#
#---------------------------------------------------------------------#
library(tidyverse)
library(haven)
library(reshape2)
library(table1)
library(EmpiricalSimSurv)

#---------------------------------------------------------------------#
#---------------------------- Data input -----------------------------#
#---------------------------------------------------------------------#
# Load the dataset with complete records
dat.c <- read.csv("data/example_new.csv", header = TRUE)

# Dataset by treatment arm
dat.c.ly <- dat.c %>% filter(TRT01P=="LY") %>% dplyr::select(-TRT01P)
dat.c.gl <- dat.c %>% filter(TRT01P=="GL") %>% dplyr::select(-TRT01P)

#---------------------------------------------------------------------------#
#------------------------ Summary table by treatment -----------------------#
#---------------------------------------------------------------------------#

#### Table 1 #####
table1(~ AGE+factor(SEX)+HBA1CBL+HBA1C_wk4+HBA1C_wk12+HBA1C_wk26+HBA1C_wk39+HBA1C_wk52+BFSGMGDL+FSG_wk4+FSG_wk12+FSG_wk26+FSG_wk39+FSG_wk52+BTRGMGDL+
         TRG_wk4+TRG_wk12+TRG_wk26+TRG_wk39+TRG_wk52+factor(NHypoe_BL)+NHypoe_wk12+NHypoe_wk26+NHypoe_wk52+THypoe_BL+THypoe_wk12+THypoe_wk26+THypoe_wk52| TRT01P, data=dat.c, overall = FALSE)

colnames(dat.c.gl)  # check this matches your type_vec positions


type_vec <- c(
  "continuous", # AGE
  "binary",     # SEX
  rep("continuous", 3),  # HBA1CBL, BFSGMGDL, BTRGMGDL
  "binary",    # NHypoe_BL
  "ordinal",
  rep("continuous", 15), # THypoe_BL -> HBA1C_wk4 -> THypoe_wk52
  rep("ordinal", 6)
)

# ── Step 1: Extract target summaries from LY arm ──────────────────────────
summ_ly <- extract_target_summaries(dat.c.ly, types = type_vec)

cat("Target means:\n");  print(round(summ_ly$means, 3))
# Target means:
#   AGE         SEX     HBA1CBL    BFSGMGDL    BTRGMGDL   NHypoe_BL   THypoe_BL   HBA1C_wk4  HBA1C_wk12  HBA1C_wk26 
# 66.303       0.647       7.408     133.448     148.418       0.003       0.077       7.126       6.685       6.593 
# HBA1C_wk39  HBA1C_wk52     FSG_wk4    FSG_wk12    FSG_wk26    FSG_wk39    FSG_wk52     TRG_wk4    TRG_wk12    TRG_wk26 
# 6.706       6.713     107.713     129.748     103.287     110.943     107.245     169.919     169.281     166.552 
# TRG_wk39    TRG_wk52 NHypoe_wk12 NHypoe_wk26 NHypoe_wk52 THypoe_wk12 THypoe_wk26 THypoe_wk52 
# 174.606     176.353       0.821       1.991       3.127       3.934       8.898      13.878 

cat("Target SDs:\n");    print(round(summ_ly$sds,   3))
# Target SDs:
#   AGE         SEX     HBA1CBL    BFSGMGDL    BTRGMGDL   NHypoe_BL   THypoe_BL   HBA1C_wk4  HBA1C_wk12  HBA1C_wk26 
# 8.718       0.478       0.820      41.327      96.623       0.058       0.245       0.932       0.642       0.821 
# HBA1C_wk39  HBA1C_wk52     FSG_wk4    FSG_wk12    FSG_wk26    FSG_wk39    FSG_wk52     TRG_wk4    TRG_wk12    TRG_wk26 
# 0.750       0.912      31.637      55.718      24.195      41.337      42.217      47.823      59.782      67.340 
# TRG_wk39    TRG_wk52 NHypoe_wk12 NHypoe_wk26 NHypoe_wk52 THypoe_wk12 THypoe_wk26 THypoe_wk52 
# 67.260      92.269       1.891       4.840       8.042       5.981      13.333      19.952 

# ── Step 2: Simulate ───────────────────────────────────────────────────────
set.seed(42)
result <- run_simulation(
  dat_ref        = dat.c.gl,
  types          = type_vec,
  N_sim          = 5000,    
  target_means   = summ_ly$means,
  target_sds     = summ_ly$sds,
  target_mins    = summ_ly$mins,
  target_maxs    = summ_ly$maxs,
  scaling_method = "range",
  verbose        = TRUE
)
# 
# [1] AGE             type=continuous  alpha=1.1729 beta=0.0030 | est_mean=66.303 (tar=66.303) est_sd=8.892 (tar=8.718)
# [2] SEX             type=binary      alpha=1.2063 beta=0.0000 | est_mean=0.646 (tar=0.647) est_sd=0.478 (tar=0.478)
# [3] HBA1CBL         type=continuous  alpha=0.7330 beta=0.0320 | est_mean=7.408 (tar=7.408) est_sd=0.667 (tar=0.820)
# [4] BFSGMGDL        type=continuous  alpha=1.0552 beta=0.0020 | est_mean=133.442 (tar=133.448) est_sd=42.425 (tar=41.327)
# [5] BTRGMGDL        type=continuous  alpha=1.8580 beta=0.0010 | est_mean=148.412 (tar=148.418) est_sd=89.074 (tar=96.623)
# [6] NHypoe_BL       type=binary      alpha=0.8679 beta=0.0000 | est_mean=0.004 (tar=0.003) est_sd=0.059 (tar=0.058)
# [7] THypoe_BL       type=ordinal     alpha=1.1590 beta=0.0010 | est_mean=0.077 (tar=0.077) est_sd=0.310 (tar=0.245)
# [8] HBA1C_wk4       type=continuous  alpha=1.5980 beta=0.0030 | est_mean=7.126 (tar=7.126) est_sd=0.951 (tar=0.932)
# [9] HBA1C_wk12      type=continuous  alpha=0.7860 beta=0.0000 | est_mean=6.685 (tar=6.685) est_sd=0.645 (tar=0.642)
# [10] HBA1C_wk26      type=continuous  alpha=1.0960 beta=0.0000 | est_mean=6.593 (tar=6.593) est_sd=0.784 (tar=0.821)
# [11] HBA1C_wk39      type=continuous  alpha=1.9928 beta=0.0030 | est_mean=6.706 (tar=6.706) est_sd=0.723 (tar=0.750)
# [12] HBA1C_wk52      type=continuous  alpha=1.1750 beta=0.0020 | est_mean=6.713 (tar=6.713) est_sd=0.913 (tar=0.912)
# [13] FSG_wk4         type=continuous  alpha=1.2350 beta=0.0050 | est_mean=107.713 (tar=107.713) est_sd=29.032 (tar=31.637)
# [14] FSG_wk12        type=continuous  alpha=1.7900 beta=0.1420 | est_mean=129.751 (tar=129.748) est_sd=39.927 (tar=55.718)
# [15] FSG_wk26        type=continuous  alpha=0.8630 beta=0.0010 | est_mean=103.284 (tar=103.287) est_sd=24.568 (tar=24.195)
# [16] FSG_wk39        type=continuous  alpha=1.5900 beta=0.0730 | est_mean=110.943 (tar=110.943) est_sd=34.279 (tar=41.337)
# [17] FSG_wk52        type=continuous  alpha=1.7000 beta=0.1250 | est_mean=107.248 (tar=107.245) est_sd=28.151 (tar=42.217)
# [18] TRG_wk4         type=continuous  alpha=0.3270 beta=0.0020 | est_mean=169.919 (tar=169.919) est_sd=34.449 (tar=47.823)
# [19] TRG_wk12        type=continuous  alpha=0.3610 beta=0.0330 | est_mean=169.279 (tar=169.281) est_sd=36.716 (tar=59.782)
# [20] TRG_wk26        type=continuous  alpha=0.5640 beta=0.0000 | est_mean=166.552 (tar=166.552) est_sd=60.475 (tar=67.340)
# [21] TRG_wk39        type=continuous  alpha=0.4450 beta=0.0060 | est_mean=174.619 (tar=174.606) est_sd=52.722 (tar=67.260)
# [22] TRG_wk52        type=continuous  alpha=0.5510 beta=0.0000 | est_mean=176.356 (tar=176.353) est_sd=89.855 (tar=92.269)
# [23] NHypoe_wk12     type=ordinal     alpha=1.9150 beta=0.0000 | est_mean=0.821 (tar=0.821) est_sd=1.913 (tar=1.891)
# [24] NHypoe_wk26     type=ordinal     alpha=1.9110 beta=0.0120 | est_mean=1.991 (tar=1.991) est_sd=3.877 (tar=4.840)
# [25] NHypoe_wk52     type=ordinal     alpha=1.8730 beta=0.0170 | est_mean=3.127 (tar=3.127) est_sd=5.663 (tar=8.042)
# [26] THypoe_wk12     type=ordinal     alpha=1.3210 beta=0.0020 | est_mean=3.934 (tar=3.934) est_sd=5.373 (tar=5.981)
# [27] THypoe_wk26     type=ordinal     alpha=1.4860 beta=0.0020 | est_mean=8.897 (tar=8.898) est_sd=11.867 (tar=13.333)
# [28] THypoe_wk52     type=ordinal     alpha=1.2020 beta=0.0210 | est_mean=13.878 (tar=13.878) est_sd=14.299 (tar=19.952)
# 
# --- Simulation Complete ---
#   N = 5000  | columns = 28 

sim_data <- result$sim_data

# ── Step 3: Diagnostic comparison (simulated vs LY target) ────────────────
cat("\n=== Mean comparison: Simulated vs LY target ===\n")
comp <- data.frame(
  column    = colnames(dat.c.ly),
  type      = type_vec,
  LY_mean   = round(colMeans(dat.c.ly,  na.rm = TRUE), 3),
  Sim_mean  = round(colMeans(sim_data,  na.rm = TRUE), 3),
  LY_sd     = round(apply(dat.c.ly, 2, sd, na.rm = TRUE), 3),
  Sim_sd    = round(apply(sim_data, 2, sd, na.rm = TRUE), 3)
)
print(comp)

# ── Step 4: table1 side-by-side ───────────────────────────────────────────
dat_compare <- dplyr::bind_rows(
  dat.c.ly  %>% dplyr::mutate(Source = "LY (Target)"),
  sim_data  %>% dplyr::mutate(Source = "Simulated from GL")
)

table1(~ AGE + factor(SEX) +
         HBA1CBL + HBA1C_wk4 + HBA1C_wk12 + HBA1C_wk26 + HBA1C_wk39 + HBA1C_wk52 +
         BFSGMGDL + FSG_wk4  + FSG_wk12  + FSG_wk26  + FSG_wk39  + FSG_wk52  +
         BTRGMGDL + TRG_wk4  + TRG_wk12  + TRG_wk26  + TRG_wk39  + TRG_wk52  +
         factor(NHypoe_BL) + NHypoe_wk12 + NHypoe_wk26 + NHypoe_wk52 +
         THypoe_BL + THypoe_wk12 + THypoe_wk26 + THypoe_wk52 | Source,
       data = dat_compare, overall = FALSE)
