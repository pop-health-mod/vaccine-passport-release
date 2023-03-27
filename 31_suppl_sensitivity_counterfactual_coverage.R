# Libraries & functions ----
library(fixest)
library(splines)
library(tidyverse)
library(data.table)
library(lubridate)
library(MetBrewer)

source("./src/utils_load_data.R")
source("./src/utils_its_bootstrap.R")
source("./src/utils_sensitivity.R")
source("./02_setup_regression_formulas.R")

# which dataset to use
which_dataset <- "passport"

for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  DO_CMA <- FALSE
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  cat(paste("========== Running models for", ifelse(DO_CMA, CMA, PROVINCE), "==========\n"))
  
  # add CMA name to the end of saved files
  CMA_suffix <- ifelse(DO_CMA, paste("_", CMA, sep = ""), "")
  
  path_out <- sprintf("./vaccine-passport-data/out/its-sensitivity-%s/ctfl-coverage", PROVINCE)
  
  # Load vaccination coverage data ----
  source("./04_setup_load_data.R", echo = TRUE)
  rm(data_coverage)
  
  # Set up variables for plot and regressions ----
  date_breaks <- seq(min(data_passport$date_wk_end), max(data_passport$date_wk_end), by = 14)
  
  # specify knot positions (use only pre-intervention period, exclude first day of passport impact)
  knot_pos <- quantile(unique(data_passport[date_wk_end < PASSPORT_START]$week), probs = c(.10, .50, .90))
  knot_dates <- min(data_passport$date_wk_end) + knot_pos*7
  print(knot_pos)
  print(knot_dates)
  
  ## Model 1: age only ----
  ## create the vaccine coverage x passport interaction variable
  data_passport[, vaxcov_passport := vaxcov_grp]
  
  fmla_age <- paste(fmla_age,
                    "vaxcov_passport:pass_anno", "vaxcov_passport:pass_anno:week_anno",
                    sep = " + ")
  
  fit_1 <- femlm(as.formula(fmla_age),
                 data = data_passport, family = "negbin",
                 offset =~log(n_unvax),
                 notes = F, mem.clean = T)
  
  # save outputs as csv and plot rate and coverage
  predict_ctfl_vaxcovgrp(data_passport, its_fit = fit_1,
                         var_select = c("week_anno", "vaxcov_grp", "vaxcov_passport"),
                         strat_var = "age")
  
  ## Model 2: SDOH only ----
  data_passport_da[, vaxcov_passport := vaxcov_grp]
  
  fmla_income <- paste(fmla_income,
                       "vaxcov_passport:pass_anno", "vaxcov_passport:pass_anno:week_anno",
                       sep = " + ")
  fmla_vismin <- paste(fmla_vismin,
                       "vaxcov_passport:pass_anno", "vaxcov_passport:pass_anno:week_anno",
                       sep = " + ")
  
  ## 2a: income ----
  fit_2a <- femlm(as.formula(fmla_income),
                  data = data_passport_da, family = "negbin",
                  offset =~log(n_unvax),
                  notes = F, mem.clean = T)
  
  # save outputs as csv and plot rate and coverage
  predict_ctfl_vaxcovgrp(data_passport_da, its_fit = fit_2a,
                         var_select = c("week_anno", "vaxcov_grp", "vaxcov_passport"),
                         strat_var = "quin_income")
  
  ## 2b: visible minority ----
  fit_2b <- femlm(as.formula(fmla_vismin),
                  data = data_passport_da, family = "negbin",
                  offset =~log(n_unvax),
                  notes = F, mem.clean = T)
  
  # save outputs as csv and plot rate and coverage
  predict_ctfl_vaxcovgrp(data_passport_da, its_fit = fit_2b,
                         var_select = c("week_anno", "vaxcov_grp", "vaxcov_passport"),
                         strat_var = "quin_vismin")
}
