# Libraries and setup ----
library(tidyverse)
library(data.table)
library(lubridate)

source("./src/utils_load_data.R")

# which dataset to use
which_dataset <- "descriptive"

# list to store tables
tbl_1_prov <- vector("list", 2)
names(tbl_1_prov) <- c("qc", "on")
tbl_s1_cma <- vector("list", 2)
names(tbl_s1_cma) <- c("mtl", "tor")

for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  for(DO_CMA in c(FALSE, TRUE)){
    # add CMA name to the end of saved files
    CMA_suffix <- ifelse(DO_CMA, paste("_", CMA, sep = ""), "")
    
    # Load data  and keep select timepoints ----
    source("./04_setup_load_data.R", echo = TRUE)
    
    data_coverage <- data_coverage[
      date_wk_end %in% c(as.Date("2021-07-03"), PASSPORT_START - 7, MODEL_END)
    ]
    
    ## compute pop size and number of doses
    # for each age group
    data_cov_age <- data_coverage[, lapply(.SD, sum),
                                  by = .(age, date_wk_start, date_wk_end),
                                  .SDcols = c("pop_rpdb", "n_vacc_1dose")]
    
    # for total population
    data_cov_prov <- data_coverage[, lapply(.SD, sum),
                                   by = .(date_wk_start, date_wk_end),
                                   .SDcols = c("pop_rpdb", "n_vacc_1dose")]
    
    # join and set factor for reordering
    data_cov <- bind_rows(data_cov_age, data_cov_prov)
    rm(data_cov_age, data_cov_prov)
    
    data_cov$age[is.na(data_cov$age)] <- "total"
    data_cov[, age := factor(age, c("total", unique(data_coverage$age)))]
    
    # compute vaccination coverage
    data_cov[, cv_1dose := (n_vacc_1dose / pop_rpdb) * 100]
    
    # turn data.table into wide format
    data_cov[, date := paste(month.abb[month(date_wk_end)], day(date_wk_end), sep = "_")]
    
    tbl_1_tmp <- data_cov %>% 
      select(date, age, pop_rpdb, n_vacc_1dose, cv_1dose) %>% 
      pivot_wider(names_from = date, values_from = c(n_vacc_1dose, cv_1dose)) %>% 
      arrange(age)
    
    # get number of DAs
    nb_da <- data_coverage[date_wk_start == min(date_wk_start) & age == min(age)] %>% nrow()
    tbl_1_tmp <- tbl_1_tmp %>% 
      mutate(da_number = ifelse(age == "total", nb_da, NA), .after = 1)
    
    # Save current region into list ----
    if(DO_CMA){
      tbl_1_tmp <- tbl_1_tmp %>% 
        mutate(region = CMA, .before = 1)
      tbl_s1_cma[[CMA]] <- tbl_1_tmp
    } else {
      tbl_1_tmp <- tbl_1_tmp %>% 
        mutate(region = PROVINCE, .before = 1)
      tbl_1_prov[[PROVINCE]] <- tbl_1_tmp
    }
  }
}

# Join regions together and save ----
# Quebec and Ontario
tbl_1_final <- bind_rows(tbl_1_prov)
tbl_1_final <- select(tbl_1_final, region, age, da_number, pop_rpdb, starts_with("cv"), starts_with("n_Vacc"))
fwrite(tbl_1_final, "./out/manuscript-tables/table_1_popsize_vc_prov.csv")

# Montreal and Toronto
tbl_s1_final <- bind_rows(tbl_s1_cma)
tbl_s1_final <- select(tbl_s1_final, region, age, da_number, pop_rpdb, starts_with("cv"), starts_with("n_Vacc"))
fwrite(tbl_s1_final, "./out/manuscript-tables/table_S1_popsize_vc_cma.csv")

