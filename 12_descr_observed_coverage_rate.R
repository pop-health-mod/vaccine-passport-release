# Libraries and setup ----
library(tidyverse)
library(data.table)
library(lubridate)
library(MetBrewer)

source("./src/utils_load_data.R")

# which dataset to use
which_dataset <- "descriptive"

for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  for(DO_CMA in c(FALSE, TRUE)){
    # add CMA name to the end of saved files
    CMA_suffix <- ifelse(DO_CMA, paste("_", CMA, sep = ""), "")
    
    if(DO_CMA){
      out_path <- sprintf("../vaccine-passport-data/out/observed-%s-%s", PROVINCE, CMA)
    } else {
      out_path <- sprintf("../vaccine-passport-data/out/observed-%s", PROVINCE)
    }
    
    source("./04_setup_load_data.R", echo = TRUE)
    
    data_coverage_da <- data_passport_da
    rm(data_passport_da)
    
    # Age or SDOH stratified coverage ----
    # single stratification by age
    data_cov_age <- data_coverage[,
                                  lapply(.SD, sum),
                                  by = .(date_wk_start, date_wk_end, age),
                                  .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose")]
    write.csv(data_cov_age, sprintf("%s/cover_1_age%s.csv", out_path, CMA_suffix), row.names = F)
    
    ### SDOH, single panels
    # income
    data_cov_income <- data_coverage_da[,
                                        lapply(.SD, sum), 
                                        by = .(date_wk_start, date_wk_end, quin_income),
                                        .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose")][order(quin_income, date_wk_start)]
    write.csv(data_cov_income, sprintf("%s/cover_2a_income%s.csv", out_path, CMA_suffix), row.names = F)
    
    # visible minority
    data_cov_vismin <- data_coverage_da[,
                                        lapply(.SD, sum), 
                                        by = .(date_wk_start, date_wk_end, quin_vismin),
                                        .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose")][order(quin_vismin, date_wk_start)]
    write.csv(data_cov_vismin, sprintf("%s/cover_2b_vismin%s.csv", out_path, CMA_suffix), row.names = F)
    
    # Age x SDOH interaction coverage ----
    # income x age
    data_cov_age.income <- data_coverage[,
                                         lapply(.SD, sum), 
                                         by = .(date_wk_start, date_wk_end, quin_income, age),
                                         .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose")][order(quin_income, age, date_wk_start)]
    write.csv(data_cov_age.income, sprintf("%s/cover_3a_age.income%s.csv", out_path, CMA_suffix), row.names = F)
    
    # visible minority x age
    data_cov_age.vismin <- data_coverage[,
                                         lapply(.SD, sum), 
                                         by = .(date_wk_start, date_wk_end, quin_vismin, age),
                                         .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose")][order(quin_vismin, age, date_wk_start)]
    write.csv(data_cov_age.vismin, sprintf("%s/cover_3b_age.vismin%s.csv", out_path, CMA_suffix), row.names = F)
    
  }
}
