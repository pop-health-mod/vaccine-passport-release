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

for(PROVINCE in c("qc", "on")){
  # policy dates
  CMA <- "none"
  DO_CMA <- FALSE
  source("./03_setup_policy_dates.R")
  
  out_path <- sprintf("./out/its-sensitivity-%s", PROVINCE)
  
  # 1. Change starting date of time-series ----
  cat("Sensitivity 1:\tChange starting date of time-series ==========================\n")
  
  # load dataset
  which_dataset <- "sensitivity_ts_start"
  source("./04_setup_load_data.R", echo = TRUE)
  rm(data_coverage, data_passport_da)
  
  # dates to test as an end date for the impact of the passport
  proposed_start_ts_dates <- as.character(as.Date("2021-07-03") + c(-7, 0, 7))
  
  df_model_dx <- tibble(start_ts_date = character(), knot_position.wk = character(),
                        knot_position.date = character(),
                        aic_unvax = numeric(), bic_unvax = numeric(),
                        ssr_null = numeric(), dev_fit = numeric())
  df_fit_rate.cover <- tibble()
  
  start_time <- Sys.time()
  for(start_date in proposed_start_ts_dates){
    ## Model loop ----
    ## track progress
    run_time <- Sys.time() - start_time
    print(round(run_time, 2))
    
    ## Temporary dataset ----
    ## create temporary dataset with the time-series starting at the specified end point
    start_date <- as.Date(start_date)
    data_tmp <- copy(data_passport)
    data_tmp <- data_tmp[date_wk_end >= start_date]
    
    # add vaccine coverage as of start of the time-series
    data_tmp <- get_vaxcoverage(data_tmp, time_coverage = min(data_tmp$date_wk_end))
    data_tmp[, vaxcov_grp := case_when(vaxcov_at_t >= .90 ~ "90_over",
                                       vaxcov_at_t >= .80 ~ "80_89.9_", vaxcov_at_t >= .70 ~ "70_79.9_",
                                       vaxcov_at_t >= .60 ~ "60_69.9_", vaxcov_at_t >= .50 ~ "50_59.9_",
                                       vaxcov_at_t <  .50 ~ "49.9_under",
                                       is.na(vaxcov_at_t) ~ "90_over")]
    
    data_tmp <- data_tmp %>% 
      select(codeDA:date_wk_end, vaxcov_grp, week, week_anno,
             pass_anno, pop_rpdb:rate_1dose)
    
    # adjust week indicator so that first day of time-series is 0th date
    data_tmp[, week := as.numeric(date_wk_end - as.Date(start_date))/7]
    
    # inspect indicators
    data_tmp[codeDA == min(codeDA) & age == min(age),
             .(date_wk_start, date_wk_end, week, week_anno, pass_anno)]
    
    ## Model fit and diagnostics ----
    # specify knot positions (use only pre-intervention period, exclude first day of passport impact)
    knot_pos <- quantile(unique(data_tmp[date_wk_end < PASSPORT_START]$week), probs = c(.10, .50, .90))
    knot_dates <- min(data_passport$date_wk_end) + knot_pos*7
    print(knot_pos)
    print(knot_dates)
    
    # model fit
    negb_fit_ts_start <- femlm(as.formula(fmla_age), 
                               data = data_tmp, family = "negbin",
                               offset =~log(n_unvax),
                               notes = F, mem.clean = T)
    
    # store model diagnostics
    df_model_dx <- df_model_dx %>% 
      add_row(start_ts_date = as.character(start_date),
              knot_position.wk = paste(knot_pos, collapse = "--"),
              knot_position.date = paste(knot_dates, collapse = "--"),
              
              aic_unvax = AIC(negb_fit_ts_start), bic_unvax = BIC(negb_fit_ts_start), 
              ssr_null = negb_fit_ts_start$ssr_null, dev_fit = deviance(negb_fit_ts_start))
    
    ## Estimate predicted and observed rates and coverage ----
    predictions <- get_predicted_coverage(data_tmp, negb_fit_ts_start, var_select = c("week_anno", "vaxcov_grp"),
                                          omit_observed = TRUE, return_rate = FALSE)
    
    # rename type, add marker of which parameters were used, and extract rate and coverage from list
    predictions[, type := factor(type,
                                 levels = c("Fitted", "Counterfactual"),
                                 labels = c("Passport", "No passport\n(counterfactual)"))]
    predictions[, `:=`(start_ts_date = start_date)]
    
    # summarize rates and coverage by age group
    fitted_vals <- predictions[, 
                               lapply(.SD, sum), 
                               by = .(type, start_ts_date, age, date_wk_end),
                               .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose", "rate_1dose", "rate_1dose_alt", "rate_pred_cumul")]
    
    # store in dataframe
    df_fit_rate.cover <- bind_rows(df_fit_rate.cover, fitted_vals)
  }
  run_time <- Sys.time() - start_time
  print(round(run_time, 2))
  rm(negb_fit_ts_start, data_tmp, start_date, start_time, run_time)
  
  # save outputs
  save_model_fit(df_fit_rate.cover, "start_of_timeseries")
  
  df_model_dx <- df_model_dx[, c("start_ts_date", "knot_position.wk", "knot_position.date",
                                 "aic_unvax", "bic_unvax", "ssr_null", "dev_fit")]
  write.csv(df_model_dx, sprintf("%s/start_of_timeseries_dx.csv", out_path), row.names = FALSE)
  
  # Load data for 2nd and 3rd sensitivity analyses ----
  which_dataset <- "passport"
  source("./04_setup_load_data.R", echo = TRUE)
  rm(data_coverage, data_passport_da)
  
  # 2. Change length of impact period ----
  cat("Sensitivity 2:\tChange length of impact period ===============================\n")
  
  # dates to test as an end date for the impact of the passport
  proposed_stop_dates <- as.character(PASSPORT_END + c(-7, 0, 7))
  
  df_model_dx <- tibble(date_stop = character(), knot_position.wk = character(),
                        knot_position.date = character(),
                        aic_unvax = numeric(), bic_unvax = numeric(),
                        ssr_null = numeric(), dev_fit = numeric())
  df_fit_rate.cover <- tibble()
  
  start_time <- Sys.time()
  for(stop_date in proposed_stop_dates){
    ## Model loop ----
    ## track progress
    run_time <- Sys.time() - start_time
    print(round(run_time, 2))
    
    ## Temporary dataset ----
    ## create temporary dataset with the vaccine passport impact going up to specified stop date
    stop_date <- as.Date(stop_date)
    data_tmp <- copy(data_passport)
    
    # impact of vaccine passport assumed to last from Aug 7th (exclusive) to post_imp_trend
    # if using a different trend after the impact of the passport "washes off", add an indicator
    data_tmp[, `:=`(pass_anno = case_when(date_wk_end <  PASSPORT_START ~ 0, 
                                          date_wk_end >= PASSPORT_START & date_wk_end <= stop_date ~ 1,
                                          date_wk_end >  stop_date ~ 0))]
    # centered time indicators
    data_tmp[, week_anno := week - min(data_passport[pass_anno == 1]$week)]
    
    # inspect indicators
    data_tmp[codeDA == min(codeDA) & age == min(age),
             .(date_wk_start, date_wk_end, week, week_anno, pass_anno)]
    
    data_tmp <- data_tmp %>% 
      select(codeDA:date_wk_end, vaxcov_grp, week, week_anno, pass_anno, pop_rpdb:rate_1dose)
    
    
    ## Model fit and diagnostics ----
    # specify knot positions (use only pre-intervention period, exclude first day of passport impact)
    knot_pos <- quantile(unique(data_tmp[date_wk_end < PASSPORT_START]$week), probs = c(.10, .50, .90))
    knot_dates <- min(data_passport$date_wk_end) + knot_pos*7
    print(knot_pos)
    print(knot_dates)
    
    # fit model
    negb_fit_imp_length <- femlm(as.formula(fmla_age), 
                                 data = data_tmp, family = "negbin",
                                 offset =~log(n_unvax),
                                 notes = F, mem.clean = T)
    
    # store model diagnostics
    df_model_dx <- df_model_dx %>% 
      add_row(date_stop = as.character(stop_date),
              knot_position.wk = paste(knot_pos, collapse = "--"),
              knot_position.date = paste(knot_dates, collapse = "--"),
              
              aic_unvax = AIC(negb_fit_imp_length), bic_unvax = BIC(negb_fit_imp_length),
              ssr_null = negb_fit_imp_length$ssr_null, dev_fit = deviance(negb_fit_imp_length))
    
    ## Estimate predicted and observed rates and coverage ----
    predictions <- get_predicted_coverage(data_tmp, negb_fit_imp_length, var_select = c("week_anno", "vaxcov_grp"),
                                          omit_observed = TRUE, return_rate = FALSE)
    
    # rename type, add marker of which parameters were used, and extract rate and coverage from list
    predictions[, type := factor(type,
                                 levels = c("Fitted", "Counterfactual"),
                                 labels = c("Passport", "No passport\n(counterfactual)"))]
    predictions[, `:=`(stop_date = stop_date)]
    
    # summarize rates and coverage by age group
    fitted_vals <- predictions[, 
                               lapply(.SD, sum), 
                               by = .(type, stop_date, age, date_wk_end),
                               .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose", "rate_1dose", "rate_1dose_alt", "rate_pred_cumul")]
    
    # store in dataframe
    df_fit_rate.cover <- bind_rows(df_fit_rate.cover, fitted_vals)
  }
  run_time <- Sys.time() - start_time
  print(round(run_time, 2))
  rm(negb_fit_imp_length, data_tmp, stop_date, start_time, run_time)
  
  ## save outputs
  # add length of impact and reorder
  df_model_dx$impact_length <- paste((ceiling(as.Date(df_model_dx$date_stop) - PASSPORT_START) / 7)+1, "weeks")
  df_model_dx <- df_model_dx[, c("impact_length", "date_stop", 
                                 "knot_position.wk", "knot_position.date",
                                 "aic_unvax", "bic_unvax", "ssr_null", "dev_fit")]
  
  write.csv(df_model_dx, sprintf("%s/length_of_impact_dx.csv", out_path), row.names = FALSE)
  
  df_fit_rate.cover
  save_model_fit(df_fit_rate.cover, "length_of_impact")
  
  # 3. Change model specification ----
  cat("Sensitivity 3:\tChange model specification ===================================\n")
  data_tmp <- copy(data_passport)
  
  # add indicator for either 
  #     the post-passport period (QC) or 
  #     squared time and the July trend (ON)
  if(PROVINCE == "qc"){
    data_tmp[, `:=`(post_trend = case_when(date_wk_end <= PASSPORT_END ~ 0,
                                           date_wk_end >  PASSPORT_END ~ 1))]
    data_tmp[, week_post := week - min(data_tmp[post_trend == 1]$week)]
  } else if(PROVINCE == "on"){
    data_tmp[, `:=`(july_trend = case_when(date_wk_end <= "2021-07-31" ~ 1, 
                                           date_wk_end >  "2021-07-31" ~ 0),
                    week.sq = (week - median(unique(data_tmp$week)))^2)]
  }
  
  # which variables to select
  if(PROVINCE == "qc"){
    model_spec_vars <- c("post_trend", "week_post")
  } else if(PROVINCE == "on"){
    model_spec_vars <- c("july_trend", "week.sq")
  }
  
  data_tmp <- data_tmp %>% 
    select(codeDA, SAC, CTname, age, date_wk_start, date_wk_end,
           week, week_anno, pass_anno, all_of(model_spec_vars),
           pop_rpdb:rate_1dose, vaxcov_grp, vaxcov_at_t)
  
  # examine time variables
  data_tmp[codeDA == min(codeDA) & age == min(age)] %>% 
    select(all_of(c("date_wk_start", "date_wk_end",
                    "week", "pass_anno", "week_anno", model_spec_vars))) %>% 
    print()
  
  df_model_dx <- tibble(model_specification = character(),
                        knot_position.wk = character(), knot_position.date = character(),
                        aic_unvax = numeric(), bic_unvax = numeric(),
                        ssr_null = numeric(), dev_fit = numeric())
  df_fit_rate.cover <- tibble()
  
  # specify knot positions (use only pre-intervention period)
  knot_pos <- quantile(unique(data_tmp[date_wk_end < PASSPORT_START]$week), probs = c(.10, .50, .90))
  knot_dates <- min(data_tmp$date_wk_end) + knot_pos*7
  print(knot_dates)
  
  # Model specifications to use ----
  # use a linear model for both provinces
  ls_fmla <- list(final_model = fmla_age,
                  linear_model = gsub("ns\\(week, knots = knot_pos\\)", "week", fmla_age),
                  alt_model = NA)
  if(PROVINCE == "qc"){
    # QC alternative model is a linear model with a dummy variable for the post-passport period
    ls_fmla$alt_model <- paste(ls_fmla$linear_model, 
                               "+ post_trend + week_post:post_trend", 
                               "+ age:post_trend + age:week_post:post_trend")
  } else if(PROVINCE == "on"){
    # ON alternative model is a linear model with a dummy variable for July + a squared time term
    ls_fmla$alt_model <- paste(ls_fmla$linear_model, 
                               "+ july_trend + week:july_trend", 
                               "+ age:july_trend + age:week:july_trend",
                               "+ week.sq + age:week.sq")
  }
  
  start_time <- Sys.time()
  for(i in 1:length(ls_fmla)){
    ## Model loop ----
    ## track progress
    run_time <- Sys.time() - start_time
    cat("Model being tested:", names(ls_fmla)[i], "\n", sep = "\t")
    print(round(run_time, 2))
    
    ## Model fit and diagnostics ----
    # fit model
    negb_fit <- femlm(as.formula(ls_fmla[[i]]),
                      data = data_tmp, family = "negbin",
                      offset =~log(n_unvax),
                      notes = F, mem.clean = T)
    
    # store model diagnostics
    df_model_dx <- df_model_dx %>% 
      add_row(model_specification = names(ls_fmla)[i],
              knot_position.wk = paste(knot_pos, collapse = "--"),
              knot_position.date = paste(knot_dates, collapse = "--"),
              
              aic_unvax = AIC(negb_fit), bic_unvax = BIC(negb_fit),
              ssr_null = negb_fit$ssr_null, dev_fit = deviance(negb_fit))
    
    ## Estimate predicted and observed rates and coverage ----
    predictions <- get_predicted_coverage(data_tmp, negb_fit, 
                                          var_select = c("week_anno", "vaxcov_grp", model_spec_vars),
                                          omit_observed = TRUE, return_rate = FALSE)
    
    # rename type, add marker of which parameters were used, and extract rate and coverage from list
    predictions[, type := factor(type,
                                 levels = c("Fitted", "Counterfactual"),
                                 labels = c("Passport", "No passport\n(counterfactual)"))]
    predictions[, `:=`(model_specification = names(ls_fmla)[i])]
    
    # summarize rates and coverage by age group
    fitted_vals <- predictions[, 
                               lapply(.SD, sum), 
                               by = .(type, model_specification, age, date_wk_end),
                               .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose", "rate_1dose", "rate_1dose_alt", "rate_pred_cumul")]
    
    # store in dataframe
    df_fit_rate.cover <- bind_rows(df_fit_rate.cover, fitted_vals)
  }
  run_time <- Sys.time() - start_time
  print(round(run_time, 2))
  rm(negb_fit, knot_pos, knot_dates, start_time, run_time)
  
  # save outputs
  save_model_fit(df_fit_rate.cover, "alt_model_specs")
  
  df_model_dx <- df_model_dx[, c("model_specification",
                                 "knot_position.wk", "knot_position.date",
                                 "aic_unvax", "bic_unvax", "ssr_null", "dev_fit")]
  
  write.csv(df_model_dx, sprintf("%s/alt_model_specs_dx.csv", out_path), row.names = FALSE)
}
