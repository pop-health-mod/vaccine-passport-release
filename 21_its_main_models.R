# Libraries & functions ----
library(fixest)
library(splines)
library(tidyverse)
library(data.table)
library(lubridate)
library(MetBrewer)

source("./src/utils_load_data.R")
source("./src/utils_its_bootstrap.R")
source("./02_setup_regression_formulas.R")

# dataset to use
which_dataset <- "passport"

for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  for(DO_CMA in c(FALSE, TRUE)){
    # add CMA name to the end of saved files
    CMA_suffix <- ifelse(DO_CMA, paste("_", CMA, sep = ""), "")
    
    if(DO_CMA){
      path_out <- sprintf("../vaccine-passport-data/out/its-fit-%s-%s", PROVINCE, CMA)
    } else {
      path_out <- sprintf("../vaccine-passport-data/out/its-fit-%s", PROVINCE)
    }
    
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
    
    # ---- Model 1: age only ----
    ## fit model and predict rate and coverage
    start_t <- Sys.time()
    fit_1 <- femlm(as.formula(fmla_age), 
                   data = data_passport, family = "negbin",
                   offset =~log(n_unvax),
                   notes = F, mem.clean = T)
    end_t <- Sys.time()
    print(end_t - start_t)
    
    # save outputs
    pred_1 <- get_rate_and_coverage(data_passport, fit_1, c("week_anno", "vaxcov_grp"), "age")
    save_model_fit(pred_1$model_fit[order(age, date_wk_end)][order(-type)], "model1_age")
    save_model_coeff(fit_1, "model_1age")

    # Model 2: SDOH only ----
    ## 2a: income ----
    # note: ran without vaxcov_grp for aggregate dataset
    fit_2a <- femlm(as.formula(fmla_income),
                    data = data_passport_da, family = "negbin",
                    offset = ~log(n_unvax),
                    notes = F, mem.clean = T)
    
    # save outputs
    pred_2a <- get_rate_and_coverage(data_passport_da, fit_2a, c("week_anno", "vaxcov_grp"), "quin_income")
    save_model_fit(pred_2a$model_fit[order(quin_income, date_wk_end)][order(-type)], "model2a_income")
    save_model_coeff(fit_2a, "model_2a_income")
    
    ## 2b: visible minority ----
    fit_2b <- femlm(as.formula(fmla_vismin),
                    data = data_passport_da, family = "negbin",
                    offset = ~log(n_unvax),
                    notes = F, mem.clean = T)
    
    # save outputs
    pred_2b <- get_rate_and_coverage(data_passport_da, fit_2b, c("week_anno", "vaxcov_grp"), "quin_vismin")
    save_model_fit(pred_2b$model_fit[order(quin_vismin, date_wk_end)][order(-type)], "model2b_vismin")
    save_model_coeff(fit_2b, "model_2b_vismin")

    # Model 3: age x SDOH ----
    ## 3a: income ----
    start_t <- Sys.time()
    fit_3a <- femlm(as.formula(fmla_age.income),
                    data = data_passport, family = "negbin",
                    offset = ~log(n_unvax),
                    notes = F, mem.clean = T)
    end_t <- Sys.time()
    print(end_t - start_t)
    
    # save outputs
    pred_3a <- get_rate_and_coverage(data_passport, fit_3a, c("week_anno", "vaxcov_grp"), c("age", "quin_income"), pred_by_age = TRUE)
    save_model_fit(pred_3a$model_fit[order(age, quin_income, date_wk_end)][order(-type)], "model3a_age.income")
    save_model_coeff(fit_3a, "model_3a_income")

    ## 3b: visible minority ----
    start_t <- Sys.time()
    fit_3b <- femlm(as.formula(fmla_age.vismin),
                    data = data_passport, family = "negbin",
                    offset = ~log(n_unvax),
                    notes = F, mem.clean = T)
    end_t <- Sys.time()
    print(end_t - start_t)
    
    # save outputs as csv and plot rate and coverage
    pred_3b <- get_rate_and_coverage(data_passport, fit_3b, c("week_anno", "vaxcov_grp"), c("age", "quin_vismin"), pred_by_age = TRUE)
    save_model_fit(pred_3b$model_fit[order(age, quin_vismin, date_wk_end)][order(-type)], "model3b_age.vismin")
    save_model_coeff(fit_3b, "model_3b_vismin")
  }
}
