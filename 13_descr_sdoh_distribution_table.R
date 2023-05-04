# Libraries & functions ----
library(tidyverse)
library(data.table)
library(lubridate)

source("./src/utils_load_data.R")
source("./src/utils_its_bootstrap.R")

# which dataset to use
which_dataset <- "passport"

# list to store tables
tbl_1prov <- vector("list", 2)
names(tbl_1prov) <- c("qc", "on")
tbl_2cma <- vector("list", 2)
names(tbl_2cma) <- c("mtl", "tor")

for(PROVINCE in c("qc", "on")){
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  for(DO_CMA in c(FALSE, TRUE)){
    # add CMA name to the end of saved files
    CMA_suffix <- ifelse(DO_CMA, paste("_", CMA, sep = ""), "")
    
    # Load vaccination coverage data ----
    source("./04_setup_load_data.R", echo = TRUE)
    
    data_inc <- data_coverage %>% 
      filter(date_wk_end == min(date_wk_end)) %>% 
      group_by(quin_income) %>% 
      summarize(nb_da = n()/6, nb_ppl = sum(pop_rpdb),
                across(ATIPPE, .fns = c(min = min, max = max, mean = mean, median = median)))
    data_vis <- data_coverage %>% 
      filter(date_wk_end == min(date_wk_end)) %>% 
      group_by(quin_vismin) %>% 
      summarize(nb_da = n()/6, nb_ppl = sum(pop_rpdb),
                across(Visible_minority_Overall, .fns = c(min = min, max = max, mean = mean, median = median)))
    
    # TODO rename column(s) + join income and visible minority into single table
    tbl_sdoh_tmp <- bind_rows(data_inc, data_vis)
    
    # Save current region into list ----
    if(DO_CMA){
      tbl_sdoh_tmp <- tbl_sdoh_tmp %>% 
        mutate(region = CMA, .before = 1)
      tbl_2cma[[CMA]] <- tbl_sdoh_tmp
    } else {
      tbl_sdoh_tmp <- tbl_sdoh_tmp %>% 
        mutate(region = PROVINCE, .before = 1)
      tbl_1prov[[PROVINCE]] <- tbl_sdoh_tmp
    }
  }
}

# TODO subset columns and reorder
tbl_1prov <- bind_rows(tbl_1prov)
tbl_2cma <- bind_rows(tbl_2cma)

write.csv(tbl_1prov, "./out/manuscript-tables/table_S1a_sdoh_distr_provs.csv", row.names = FALSE)
write.csv(tbl_2cma, "./out/manuscript-tables/table_S1b_sdoh_distr_CMAs.csv", row.names = FALSE)
