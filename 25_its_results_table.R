# Libraries ----
library(tidyverse)
library(lubridate)

source("./src/utils_its_bootstrap.R")
source("./src/utils_load_data.R")

# lists to store tables
tbl_cover_ls <- vector("list", 2)
tbl_abs_ls <- vector("list", 2)
tbl_rel_ls <- vector("list", 2)

names(tbl_cover_ls) <- c("qc", "on")
names(tbl_abs_ls) <- c("qc", "on")
names(tbl_rel_ls) <- c("qc", "on")

for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  # input directories for data
  for(DO_CMA in c(FALSE)){
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
    
    # subset to last pre-passport timepoint
    tbl_cover <- data_observed %>% 
      select(date_wk_start:strat_var, original_lvl, pop_rpdb:n_vacc_1dose, cv_dose) %>% 
      filter(date_wk_end == (PASSPORT_START - 7))
    
    if(DO_CMA){
      tbl_cover <- tbl_cover %>% 
        mutate(region = CMA, .before = 1)
      tbl_cover_ls[[CMA]] <- tbl_cover
    } else {
      tbl_cover <- tbl_cover %>% 
        mutate(region = PROVINCE, .before = 1)
      tbl_cover_ls[[PROVINCE]] <- tbl_cover
    }
    
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
    tbl_impact <- left_join(impact_pt %>% select(type:vax_cover),
                            impact_ci %>% select(strat_var:impact_uci),
                            by = c("strat_var", "strat_lvl", "original_lvl", "date_wk_end"))
    
    if(DO_CMA){
      tbl_impact <- tbl_impact %>% 
        mutate(region = CMA, .before = 1)
      tbl_abs_ls[[CMA]] <- tbl_impact
    } else {
      tbl_impact <- tbl_impact %>% 
        mutate(region = PROVINCE, .before = 1)
      tbl_abs_ls[[PROVINCE]] <- tbl_impact
    }
    
    ## Relative impact ----
    # rescale the data, i.e., substract the number of doses prior to the passport
    data_obs_rescale <- subtract_base_doses(data_observed, data_observed, PASSPORT_START)
    data_its_rescale <- subtract_base_doses(data_its, data_observed, PASSPORT_START)
    data_boot_rescale <- subtract_base_doses(data_boot, data_observed, PASSPORT_START)
    
    ## relative impact (point estimate)
    # compare to observed
    data_impact_rel <- bind_rows(data_obs_rescale %>% mutate(type = "a_observed"),
                                 data_its_rescale %>% filter(type == "No passport\n(counterfactual)"))
    impact_rel_pt <- compute_impact_rel(data_impact_rel, MODEL_END,
                                        var_group_by = c("strat_var", "strat_lvl", "original_lvl"),
                                        numer = "n_vacc_adjusted")
    
    ## impact (confidence interval)
    impact_rel_ci <- compute_impact_rel_ci(data_boot_rescale, data_obs_rescale, MODEL_END,
                                           var_group_by = c("strat_var", "strat_lvl", "original_lvl"),
                                           numer = "n_vacc_adjusted")
    
    ## join pt and CI together
    tbl_impact_rel <- left_join(impact_rel_pt, impact_rel_ci,
                                by = c("strat_var", "strat_lvl", "original_lvl", "date_wk_end"))
    tbl_impact_rel <- tbl_impact_rel %>% select(strat_var:impact, impact_lci, impact_uci, n_vacc_adjusted)
    
    if(DO_CMA){
      tbl_impact_rel <- tbl_impact_rel %>% 
        mutate(region = CMA, .before = 1)
      tbl_rel_ls[[CMA]] <- tbl_impact_rel
    } else {
      tbl_impact_rel <- tbl_impact_rel %>% 
        mutate(region = PROVINCE, .before = 1)
      tbl_rel_ls[[PROVINCE]] <- tbl_impact_rel
    }
  }
}

## Put tables together ----
# join tables from list
tbl_cover <- bind_rows(tbl_cover_ls)
tbl_impact_abs <- bind_rows(tbl_abs_ls)
tbl_impact_rel <- bind_rows(tbl_rel_ls)

# format table outputs
tbl_cover <- tbl_cover %>% mutate(cv_dose = round(cv_dose, 1))

tbl_impact_abs <- tbl_impact_abs %>% 
  group_by(region, strat_var, original_lvl) %>% 
  mutate(
    abs_impact = sprintf("%s (%s\u2013%s)",
                         format(round(impact, 1), nsmall = 1),
                         format(round(impact_lci, 1), nsmall = 1),
                         format(round(impact_uci, 1), nsmall = 1))
  ) %>% 
  ungroup()

tbl_impact_rel <- tbl_impact_rel %>% 
  group_by(region, strat_var, original_lvl) %>% 
  mutate(
    rel_impact = sprintf("%s (%s\u2013%s)",
                         format(round_fn(impact), nsmall = 1),
                         format(round_fn(impact_lci), nsmall = 1),
                         format(round_fn(impact_uci), nsmall = 1))
  ) %>% 
  ungroup()

# create master table
tbl_cover_impact <- full_join(
  select(tbl_cover, region, strat_var, original_lvl, pre_pass_cov = cv_dose),
  select(tbl_impact_abs, region, strat_var, original_lvl, abs_impact,
                         abs = impact, abs_lci = impact_lci, abs_uci = impact_uci),
  by = c("region", "strat_var", "original_lvl")
)

tbl_cover_impact <- full_join(
  tbl_cover_impact,
  select(tbl_impact_rel, region, strat_var, original_lvl, rel_impact,
                         rel = impact, rel_lci = impact_lci, rel_uci = impact_uci),
  by = c("region", "strat_var", "original_lvl")
)

tbl_out <- tbl_cover_impact %>% 
  select(region, strat_var, original_lvl, pre_pass_cov, abs_impact, rel_impact)

# save table
write.csv(tbl_out,
          sprintf("./out/manuscript-tables/table_S3_cover_impact_abs_rel.csv"),
          row.names = FALSE)
