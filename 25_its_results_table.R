# Libraries ----
library(tidyverse)
library(lubridate)

source("./src/utils_its_bootstrap.R")
source("./src/utils_load_data.R")

# province and CMA for troubleshooting
# PROVINCE <- "qc"
# CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
# DO_CMA <- FALSE

# lists to store tables
plt_cover_ls <- vector("list", 4)
plt_impact_abs <- vector("list", 4)
plt_impact_rel <- vector("list", 4)

names(plt_cover_ls) <- c("qc", "on", "mtl", "tor")
names(plt_impact_abs) <- c("qc", "on", "mtl", "tor")
names(plt_impact_rel) <- c("qc", "on", "mtl", "tor")

for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  # input directories for data
  for(DO_CMA in c(FALSE, TRUE)){
    # Compute table components ----
    if(DO_CMA){
      path_its <- sprintf("../vaccine-passport-data/out/its-fit-%s-%s", PROVINCE, CMA)
      path_boot <- sprintf("../vaccine-passport-data/out/its-boot-%s-%s", PROVINCE, CMA)
      path_obs <- sprintf("../vaccine-passport-data/out/observed-%s-%s", PROVINCE, CMA)
      
      city_suffix <- paste("_", CMA, sep = "")
    } else {
      path_its <- sprintf("../vaccine-passport-data/out/its-fit-%s", PROVINCE)
      path_boot <- sprintf("../vaccine-passport-data/out/its-boot-%s", PROVINCE)
      path_obs <- sprintf("../vaccine-passport-data/out/observed-%s", PROVINCE)
      
      city_suffix <- ""
    }
    
    ## Pre-announcement coverage ----
    ## load observed coverage
    data_observed <- load_merge_data(path_obs, "observed")
    
    # keep only June 2021 onwards
    data_observed <- data_observed %>% filter(date_wk_end >= "2021-06-01")
    
    # compute coverage
    data_observed <- data_observed %>%
      mutate(cv_dose = (n_vacc_1dose / pop_rpdb) * 100)
    
    # TODO code for table coverage
    
    ## ITS impact estimates ----
    ## load ITS point estimates
    data_its <- load_merge_data(path_its, "its")
    
    ## bootstrap estimates
    data_boot <- load_merge_data(path_boot, "bootstrap", R = 1000)
    
    # compute rate and coverage, and get quantiles
    boot_ci <- compute_boot_ci(data_boot, c("strat_var", "strat_lvl", "original_lvl"),
                               numer = "rate_1dose", denom = "n_unvax")
    
    ## impact (point estimate)
    # compare to observed
    data_impact <- bind_rows(data_observed %>% mutate(type = "a_observed"),
                             data_its %>% filter(type == "No passport\n(counterfactual)"))
    
    impact_pt <- compute_impact(data_impact, MODEL_END,
                                var_group_by = c("strat_var", "strat_lvl", "original_lvl"))
    
    ## impact (confidence interval)
    impact_ci <- compute_impact_ci(data_boot, data_observed, MODEL_END,
                                   var_group_by = c("strat_var", "strat_lvl",
                                                    "original_lvl"))
    
    ## join pt and CI together
    impact_tbl <- left_join(impact_pt %>% select(type:vax_cover),
                            impact_ci %>% select(strat_var:impact_uci),
                            by = c("strat_var", "strat_lvl", "original_lvl", "date_wk_end"))
    
    if(DO_CMA){
      impact_tbl <- impact_tbl %>% 
        mutate(region = CMA, .before = 1)
      plt_impact_abs[[CMA]] <- impact_tbl
    } else {
      impact_tbl <- impact_tbl %>% 
        mutate(region = PROVINCE, .before = 1)
      plt_impact_abs[[PROVINCE]] <- impact_tbl
    }
    
    ## Relative impact ----
    # TODO add code to compute relative increase
    
    
    if(DO_CMA){
      impact_tbl <- impact_tbl %>% 
        mutate(region = CMA, .before = 1)
      plt_impact_rel[[CMA]] <- impact_tbl
    } else {
      impact_tbl <- impact_tbl %>% 
        mutate(region = PROVINCE, .before = 1)
      plt_impact_rel[[PROVINCE]] <- impact_tbl
    }
  }
}


