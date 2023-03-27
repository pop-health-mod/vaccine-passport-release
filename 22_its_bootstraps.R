# Libraries & functions ----
library(fixest)
library(splines)
library(tidyverse)
library(data.table)
library(lubridate)

source("./src/utils_load_data.R")
source("./src/utils_its_bootstrap.R")
source("./02_setup_regression_formulas.R")

# how many replicates to do and how to break them up
nb_replicates <- 500
split_rep <- 250
split_rep_intx <- 100

loop_main <- as.integer(nb_replicates / split_rep)
loop_intx <- as.integer(nb_replicates / split_rep_intx)

# which round of bootstraps is it
boot_rnd <- 1
boot_rnd_ttl <- 2

# which dataset to use
which_dataset <- "passport"

# specify whether to skip a region
do_region <- data.table(region = c("qc", "mtl", "on", "tor"),
                        done = c(F, F, F, F))

master_start_time <- Sys.time()

# run loop over each province and CMA
for(PROVINCE in c("qc", "on")){
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  for(DO_CMA in c(FALSE, TRUE)){
    cat(paste("========== Running models for", ifelse(DO_CMA, CMA, PROVINCE), "==========\n"))
    source("03_setup_policy_dates.R")
    
    # variable to know whether to skip current analysis
    if(DO_CMA & do_region[region == CMA]$done){
      cat("\tSkipped\n")
      next
    } else if(!DO_CMA & do_region[region == PROVINCE]$done){
      cat("\tSkipped\n")
      next
    }
    
    if(DO_CMA){
      path_out <- sprintf("../vaccine-passport-data/out/its-boot-part/%s-%s", PROVINCE, CMA)
    } else {
      path_out <- sprintf("../vaccine-passport-data/out/its-boot-part/%s", PROVINCE)
    }
    
    # Data loading and processing ----
    source("04_setup_load_data.R", echo = TRUE)
    rm(data_coverage)
    
    # specify knot positions (use only pre-intervention period)
    knot_pos <- quantile(unique(data_passport[date_wk_end < PASSPORT_START]$week), probs = c(.10, .50, .90))
    knot_dates <- min(data_passport$date_wk_end) + knot_pos*7
    print(knot_pos)
    print(knot_dates)
    
    ## Model 1: age only ----
    cat("1:\tAge model ==========\n")
    for(i in ((1:loop_main) + loop_main * (boot_rnd-1))){
      # track progress
      cat(sprintf("Run %s out of %s.\n", i, loop_main * boot_rnd))
      
      # run and save bootstrap
      boot_repl_age <- bootstrap_rate.cover(data_passport, resample_by = "CTuid", fmla = fmla_age,
                                            var_select = c("week_anno", "vaxcov_grp"), strat_var = "age",
                                            R = split_rep, track_time = TRUE)
      
      fwrite(boot_repl_age, sprintf("%s/bootstrap_model1_age.R%s-%s.csv",
                                    path_out, i * split_rep, nb_replicates * boot_rnd_ttl))
      rm(boot_repl_age)
    }
    
    ## Model 2: SDOH only ----
    ## 2a: income ----
    cat("\n2a:\tIncome model==========\n")
    for(i in ((1:loop_main) + loop_main * (boot_rnd-1))){
      # track progress
      cat(sprintf("Run %s out of %s.\n", i, loop_main * boot_rnd))
      # run and save bootstrap
      boot_repl_income <- bootstrap_rate.cover(data_passport_da, resample_by = "CTuid", fmla = fmla_income,
                                               var_select = c("week_anno", "vaxcov_grp"), strat_var = "quin_income",
                                               R = split_rep, track_time = TRUE)
      
      fwrite(boot_repl_income, sprintf("%s/bootstrap_model2a_income.R%s-%s.csv",
                                       path_out, i * split_rep, nb_replicates * boot_rnd_ttl))
      rm(boot_repl_income)
    }
    
    ## 2b: visible minority ----
    cat("\n2b:\tVisible minority==========\n")
    for(i in ((1:loop_main) + loop_main * (boot_rnd-1))){
      # track progress
      cat(sprintf("Run %s out of %s.\n", i, loop_main * boot_rnd))
      # run and save bootstrap
      boot_repl_vismin <- bootstrap_rate.cover(data_passport_da, resample_by = "CTuid", fmla = fmla_vismin,
                                               var_select = c("week_anno", "vaxcov_grp"), strat_var = "quin_vismin",
                                               R = split_rep, track_time = TRUE)
      
      fwrite(boot_repl_vismin, sprintf("%s/bootstrap_model2b_vismin.R%s-%s.csv",
                                       path_out, i * split_rep, nb_replicates * boot_rnd_ttl))
      rm(boot_repl_vismin)
    }
    
    # Model 3: age x SDOH ----
    ## 3a: income ----
    cat("\n3a:\tIncome x age model==========\n")
    for(i in ((1:loop_intx) + loop_i../vaccine-passport-data/x * (boot_rnd-1))){
      # track progress
      cat(sprintf("Run %s out of %s.\n", i, loop_intx * boot_rnd))
      # run and save bootstrap
      boot_repl_age.income <- bootstrap_rate.cover(data_passport, resample_by = "CTuid", fmla = fmla_age.income,
                                                   var_select = c("week_anno", "vaxcov_grp"), strat_var = c("age", "quin_income"),
                                                   R = split_rep_intx, track_time = TRUE,
                                                   pred_by_age = (!DO_CMA & PROVINCE == "on"))
      
      fwrite(boot_repl_age.income, sprintf("%s/bootstrap_model3a_age.income.R%s-%s.csv",
                                           path_out, i * split_rep_intx, nb_replicates * boot_rnd_ttl))
      rm(boot_repl_age.income)
    }
    
    ## 3b: visible minority ----
    cat("\n3b:\tVisible minority x age model==========\n")
    for(i in ((1:loop_intx) + loop_intx * (boot_rnd-1))){
      # track progress
      cat(sprintf("Run %s out of %s.\n", i, loop_intx * boot_rnd))
      # run and save bootstrap
      boot_repl_age.vismin <- bootstrap_rate.cover(data_passport, resample_by = "CTuid", fmla = fmla_age.vismin,
                                                   var_select = c("week_anno", "vaxcov_grp"), strat_var = c("age", "quin_vismin"),
                                                   R = split_rep_intx, track_time = TRUE,
                                                   pred_by_age = (!DO_CMA & PROVINCE == "on"))
      
      fwrite(boot_repl_age.vismin, sprintf("%s/bootstrap_model3b_age.vismin.R%s-%s.csv",
                                           path_out, i * split_rep_intx, nb_replicates * boot_rnd_ttl))
      rm(boot_repl_age.vismin)
    }
    
    cat("\n")
  }
}
master_end_time <- Sys.time()
total_run_time <- master_end_time - master_start_time
