# ---- Fns to add indicator variables ----
# adds an indicator for each DA/age of the vaccine coverage observed at time_coverage
get_vaxcoverage <- function(data, time_coverage = "2021-07-31", strata = "age"){
  # filter to only specified date
  data_vaccine_cover <- data[date_wk_end == time_coverage]
  
  # compute coverage
  data_vaccine_cover[, vaxcov_at_t := n_vacc_1dose / pop_rpdb]
  var_select <- c("codeDA", strata, "vaxcov_at_t")
  data_vaccine_cover <- data_vaccine_cover[, ..var_select]
  
  # add coverage to original data
  data <- data[data_vaccine_cover, on = c("codeDA", strata)]
  
  return(data)
}

# ---- Analysis ----
bootstrap_rate.cover <- function(data, resample_by = "CTuid", fmla,
                                 var_select = c("vaxcov_grp"), 
                                 strat_var = "age",
                                 pred_by_age = FALSE,
                                 R = 1000,
                                 track_time = TRUE){
  if(track_time){
    start_t <- Sys.time()
  }
  ### Step 1: sample by unit (with replacement)
  # variables to keep
  var_select <- c("codeDA", "CTuid", strat_var, "vaxcov_grp",
                  "week", "week_anno", "pass_anno",
                  "pop_rpdb", "n_unvax", "rate_1dose")
  
  # get the list of unique resampling units
  resampling_units <- unique(data[[resample_by]])
  
  # create a list to store the aggregated results
  boot_rate.cover <- vector(mode = "list", length = R)
  for(i in 1:R){
    # resample units with replacement
    resampled <- sample(resampling_units, size = length(resampling_units), replace = TRUE)
    
    # create the new dataframe by retrieving each unit individually
    list_resample_data <- vector("list", length = length(resampled))
    for(j in 1:length(resampled)){
      list_resample_data[[j]] <- data[data[[resample_by]] == resampled[j],
                                      ..var_select]
    }
    resampled_data <- bind_rows(list_resample_data)
    rm(list_resample_data)
    
    ### Step 2: fit model to the resampled data
    fit <- femlm(as.formula(fmla), data = resampled_data,
                 family = "negbin", offset =~log(n_unvax),
                 notes = F, mem.clean = T)
    
    ### Step 3: predict using the original data
    # simulate predicted and observed rates and coverage
    predictions <- get_predicted_coverage(data, fit, omit_observed = TRUE, return_rate = FALSE, 
                                          var_select = var_select, strat_var = strat_var, pred_by_age = pred_by_age)
    
    # summarize rates and coverage by stratifying variable
    boot_rate.cover[[i]] <- predictions[, 
                                        lapply(.SD, sum), 
                                        by = c("type", strat_var, "date_wk_end"),
                                        .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose", "rate_1dose", "rate_1dose_alt", "rate_pred_cumul")]
    # store replicate number
    boot_rate.cover[[i]][, replicate_nb := i]
    
    # remove variables to avoid memory issues
    rm(fit, predictions)
    
    # print progress at every 10% of the replicates
    if((i/R*10) %in% (1:10) & track_time){
      print(sprintf("Done %s out of %s bootstrap replicates.", i, R))
      print(Sys.time() - start_t)
    }
  }
  boot_rate.cover <- bind_rows(boot_rate.cover)
  
  return(boot_rate.cover)
}

# estimates predicted rate for supplied data, using its_fit
get_predicted_rate <- function(data, its_fit, predict_fit = FALSE,
                               date_start_cov = NULL, var_select = NULL, strat_var = "age",
                               pred_by_age = FALSE){
  if(is.null(date_start_cov)){
    date_start_cov <- min(data$date_wk_end)
  }
  
  ### create dataset with observed, predicted counterfactual and predicted intervention
  # the variables to be selected MUST include the stratifying element in the formula
  vars_to_select <- c("codeDA", var_select, strat_var, "date_wk_end", "week", 
                      "pass_anno", "pop_rpdb", "n_unvax", "n_vacc_1dose")
  
  # 2) predicted counterfactual rate
  # get data and set indicator to 0
  data_pred_fits_2 <- data %>% 
    select(all_of(vars_to_select)) %>% 
    mutate(rate_1dose = NA, type = "Counterfactual", pass_anno = 0)
  
  # 3) predicted intervention (i.e., observed) rate
  if(predict_fit){
    data_pred_fits_3 <- data %>% 
      select(all_of(vars_to_select)) %>% 
      mutate(rate_1dose = NA, type = "Fitted")
  }
  
  if(predict_fit){
    data_pred_fits <- bind_rows(data_pred_fits_2, data_pred_fits_3)
    rm(data_pred_fits_2, data_pred_fits_3)
  } else {
    data_pred_fits <- data_pred_fits_2
    rm(data_pred_fits_2)
  }
  
  # for 2&3), predict rate and translate to coverage
  # predict rate, set to 0 pre-intervention and get cumulative intervention
  # predicting each age separately can be used to avoid memory issues (in interaction analyses)
  if(pred_by_age){
    age_vec <- unique(data_pred_fits$age)
    data_pred_fits$rate_1dose <- NA_real_
    
    for(cur_age in age_vec){
      data_pred_fits$rate_1dose[data_pred_fits$age == cur_age] <- predict(its_fit, data_pred_fits[age == cur_age], type = "response")
    }
  } else {
    data_pred_fits$rate_1dose <- predict(its_fit, data_pred_fits, type = "response")
  }
  
  return(data_pred_fits)
}

# estimates predicted coverage for supplied data, using its_fit
#     step 1: predict rates using get_predicted_rate()
#     step 2: apply fitted vaccination rates by computing cumulative sum of vaccination starting from date_start_cov
get_predicted_coverage <- function(data, its_fit, predict_fit = TRUE,
                                   date_start_cov = NULL, compute_cv = FALSE, 
                                   var_select = NULL, strat_var = "age",
                                   omit_observed = FALSE, return_rate = FALSE, pred_by_age = FALSE){
  if(is.null(date_start_cov)){
    date_start_cov <- min(data$date_wk_end)
  }
  
  ### create dataset with observed, predicted counterfactual and predicted intervention
  # the variables to be selected must include the stratifying element in the formula
  vars_to_select <- c("codeDA", var_select, strat_var, "date_wk_end", "week", 
                      "pass_anno", "pop_rpdb", "n_unvax", "n_vacc_1dose")
  # 1) observed
  data_pred_cover <- data %>% 
    select(all_of(vars_to_select)) %>% 
    mutate(type = "Observed")
  
  # get rate for 2) predicted counterfactual and 3) predicted intervention
  data_pred_fits <- get_predicted_rate(data, its_fit, predict_fit, date_start_cov, 
                                       var_select = var_select, strat_var = strat_var, pred_by_age = pred_by_age)
  if(return_rate){
    data_pred_rates <- data_pred_fits
  }
  
  # set all time points not to be predicted to a rate of 0
  data_pred_fits$rate_1dose_alt <- ifelse(data_pred_fits$date_wk_end <= date_start_cov, 0, data_pred_fits$rate_1dose)
  
  data_pred_fits <- data_pred_fits %>%
    group_by(across(c("type", "codeDA", all_of(strat_var)))) %>%
    mutate(rate_pred_cumul = cumsum(rate_1dose_alt)) %>%
    ungroup()
  
  # get vaccination coverage at starting point pre-intervention
  data_starting_coverage <- data %>% 
    filter(date_wk_end == date_start_cov) %>% select(codeDA, all_of(strat_var), n_1dose_starting = n_vacc_1dose)
  
  # retrieve cumulative vaccination rate and join to coverage pre-intervention
  data_pred_fits <- data_pred_fits %>%
    left_join(data_starting_coverage, by = c("codeDA", strat_var))
  
  # compute cumulative doses and predicted coverage
  data_pred_fits <- data_pred_fits %>% 
    mutate(n_vacc_1dose = ifelse(date_wk_end <= date_start_cov, n_vacc_1dose, n_1dose_starting + rate_pred_cumul))
  
  # add predicted data to observed
  data_pred_cover <- bind_rows(data_pred_cover, data_pred_fits)
  
  # remove observed if not required, e.g. if computing the predicted rate for different models
  if(omit_observed){
    data_pred_cover <- data_pred_cover[type != "Observed"]
  }
  
  if(compute_cv){
    data_pred_cover$cv_1dose <- data_pred_cover$n_vacc_1dose / data_pred_cover$pop_rpdb * 100
    
    # correct for "over-vaccination"
    data_pred_cover <- data_pred_cover %>% mutate(cv_1dose = ifelse(cv_1dose > 100, 100, cv_1dose))
  }
  
  if(return_rate){
    return(list(rate = data_pred_rates, cover = data_pred_cover))
  } else {
    return(data_pred_cover)
  }
}

get_rate_and_coverage <- function(data, its_fit, var_select = NULL, strat_var, 
                                  predict_fit = TRUE, denom = "n_unvax", pred_by_age = FALSE){
  # simulate predicted and observed rates and coverage
  predictions <- get_predicted_coverage(data, its_fit, predict_fit = predict_fit,
                                        omit_observed = TRUE, return_rate = FALSE,
                                        var_select = var_select, strat_var = strat_var, pred_by_age = pred_by_age)
  
  # rename type
  predictions[, type := factor(type,
                               levels = c("Fitted", "Counterfactual"),
                               labels = c("Passport", "No passport\n(counterfactual)"))]
  
  # summarize rates and coverage by stratifying variable
  fitted_vals <- predictions[, 
                             lapply(.SD, sum), 
                             by = c("type", strat_var, "date_wk_end"),
                             .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose", "rate_1dose", "rate_1dose_alt", "rate_pred_cumul")]
  
  # summarize for the observed data
  data_obs <- data[,
                   lapply(.SD, sum), 
                   by = c(strat_var, "date_wk_end"),
                   .SDcols = c("pop_rpdb", "n_unvax", "n_vacc_1dose", "rate_1dose")]
  
  # vaccine coverage
  data_obs$rate_plt <- data_obs[["rate_1dose"]] / data_obs[[denom]] * 10^5
  data_obs$cv_dose <- data_obs$n_vacc_1dose / data_obs$pop_rpdb * 100
  
  return(list(model_fit = fitted_vals, data_obs = data_obs))
}

# ---- Process result ----
# compute impact point estimates across selected grouping variables
# point estimate is estimated as coverage_fit - coverage_counterfactual
# code assumes only fit (or observed) and counterfactual are supplied, and that 
# the counterfactual is the second level in the factor
compute_impact <- function(data, date_select, var_group_by = "age",
                           numer = "n_vacc_1dose", denom = "pop_rpdb"){
  # limit to selected date for impact and order
  data <- data %>% filter(date_wk_end == date_select)
  data <- data %>% arrange(across(all_of(c(var_group_by, "type"))))
  
  # compute coverage
  data$vax_cover <- data[[numer]] / data[[denom]] * 100
  
  # compute gap between observed/fitted and counterfactual scenario
  data <- data %>% 
    group_by(across(all_of(var_group_by))) %>% 
    mutate(impact = vax_cover - lead(vax_cover)) %>% 
    ungroup()
  data <- data %>% filter(!is.na(impact))
  
  # select only necessary variables
  data <- data %>% select(type, all_of(var_group_by), date_wk_end, impact, 
                          vax_cover, all_of(c(numer, denom)))
  
  return(data)
}

compute_impact_ci <- function(data_bootstrap, data_observed_,
                              date_select, var_group_by = "age",
                              numer = "n_vacc_1dose", denom = "pop_rpdb"){
  # create a replicate of the observed data for each bootstrap
  indx_rep <- which(data_observed_$date_wk_end == MODEL_END)
  data_obs_dupl <- data_observed_[rep(indx_rep, max(data_bootstrap$replicate_nb)),] %>% 
    group_by(across(all_of(var_group_by))) %>% 
    mutate(type = "a_observed", replicate_nb = 1:n()) %>% 
    ungroup()
  
  # compute the impact in each replicate
  data_impact <- bind_rows(data_bootstrap %>% filter(type == "No passport\n(counterfactual)"),
                           data_obs_dupl)
  impact_ci <- compute_impact(data_impact, date_select = date_select, 
                              var_group_by = c("replicate_nb", var_group_by),
                              numer = numer, denom = denom)
  
  # compute 95% confidence interval
  impact_ci <- impact_ci %>% 
    group_by(across(all_of(var_group_by)), date_wk_end) %>% 
    summarize(impact_lci = quantile(impact, .025), impact_uci = quantile(impact, .975),
              vax_cover_lci = quantile(vax_cover, .025), vax_cover_uci = quantile(vax_cover, .975),
              .groups = "drop")
  
  # return dataset
  return(impact_ci)
}

# compute bootstrapped confidence intervals
compute_boot_ci <- function(data, group_var,
                            numer = "rate_1dose", denom = "n_unvax"){
  
  # instantiate new object if data is a data.table
  if(class(data)[1] == "data.table"){
    data <- copy(data)
  }
  
  # compute coverage
  data$cv_dose <- (data[["n_vacc_1dose"]] / data[["pop_rpdb"]]) * 100
  
  # compute rate
  data$rate_plt <- (data[[numer]] / data[[denom]]) * 10^5
  
  # get 95% confidence intervals
  data <- data %>% 
    group_by(type, across(all_of(group_var)), date_wk_end) %>% 
    summarize(rate_plt_lci = quantile(rate_plt, .025), rate_plt_uci = quantile(rate_plt, .975),
              cv_reg_lci = quantile(cv_dose, .025), cv_reg_uci = quantile(cv_dose, .975),
              .groups = "drop")
  
  return(data)
}

# ---- Helper/miscellaneous fns ----
# saves ITS outputs
save_model_fit <- function(predictions, model_name){
  if(DO_CMA){
    save_path <- sprintf("%s/predicted_%s_fitted_values_%s.csv", 
                         path_out, model_name, CMA)
  } else {
    save_path <- sprintf("%s/predicted_%s_fitted_values.csv", 
                         path_out, model_name)
  }
  write.csv(predictions, save_path, row.names = FALSE)
}

# store model coefficients
save_model_coeff <- function(fit, model_name){
  broom::tidy(fit) %>% 
      mutate(estimate_exp = exp(estimate), .after = estimate) %>% 
      fwrite(file = sprintf("%s/its-coeff/%s_%s%s.csv", path_out, model_name, PROVINCE, CMA_suffix))
}
