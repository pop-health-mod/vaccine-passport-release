# Libraries & functions ----
library(tidyverse)
library(data.table)
library(lubridate)
library(MetBrewer)

theme_set(theme_bw())

source("./src/utils_load_data.R")
source("./src/utils_its_bootstrap.R")
source("./src/utils_sensitivity.R")

source("./src/plot_main.R")
source("./src/plot_setup.R")

fig_path <- "./response-reviewers/fig"

## functions
# redefine compute_impact for this context
# Alternative model specification functions ----
compute_impact_sens <- function(data_sensitivity, data_obs, sens_var,
                                drop_main_model = TRUE){
  # drop main ITS from sensitivity (using fit from main model)
  if(drop_main_model){
    data_sensitivity <- data_sensitivity[!grepl("main model", data_sensitivity[[sens_var]]), ]
  }
  nb_lvls <- length(unique(data_sensitivity[[sens_var]]))
  
  ## impact point estimates, comparing to observed observed
  # create a replicate of the observed data (only last timepoint) for each sensitivity fit
  indx_rep <- which(data_obs$date_wk_end == MODEL_END)
  data_obs_dupl <- data_obs[rep(indx_rep, nb_lvls), ] %>% 
    group_by(age) %>% 
    mutate(type = "A_observed", obs_rep = (1:n()) + drop_main_model) %>% 
    ungroup()
  
  data_obs_dupl[[sens_var]] <- factor(data_obs_dupl$obs_rep, 
                                      1:(nb_lvls + drop_main_model),
                                      levels(data_sensitivity[[sens_var]]))
  
  data_impact <- bind_rows(
    data_obs_dupl,
    data_sensitivity %>% filter(type == "No passport\n(counterfactual)")
  )
  
  impact_pt_estim <- compute_impact(data_impact, MODEL_END, var_group_by = c(sens_var, "age"))
  
  ## join pt and CI together
  impact_pt_estim <- impact_pt_estim %>% select(type:vax_cover, pop_rpdb, n_vacc_1dose)
  
  # create label (grouping is to avoid creating a space in the non-negative estimates)
  # case_when is to do CI only for main model
  impact_pt_estim <- impact_pt_estim %>% 
    group_by(across(all_of(sens_var)), age) %>% 
    mutate(
      impact_lab = sprintf("%s %% pts.", format(round(impact, 1), nsmall = 1))
    ) %>% 
    ungroup()
  
  return(impact_pt_estim)
}

compute_impact_global <- function(data_sensitivity, data_obs, sens_var,
                                  drop_main_model = FALSE){
  # drop main ITS from sensitivity (using fit from main model)
  if(drop_main_model){
    data_sensitivity <- data_sensitivity[!grepl("main model", data_sensitivity[[sens_var]]), ]
  }
  nb_lvls <- length(unique(data_sensitivity[[sens_var]]))
  # browser()
  # compute province-wide impact
  data_obs <- data_obs %>% 
    group_by(date_wk_end) %>% 
    summarize(across(c(pop_rpdb, n_unvax, n_vacc_1dose), sum))
  
  ## impact point estimates, comparing to observed observed
  # create a replicate of the observed data (only last timepoint) for each sensitivity fit
  indx_rep <- which(data_obs$date_wk_end == MODEL_END)
  data_obs_dupl <- data_obs[rep(indx_rep, nb_lvls), ] %>% 
    mutate(type = "A_observed", obs_rep = (1:n()) + drop_main_model) %>% # only shift by 1 if main model already computed
    ungroup()
  
  data_obs_dupl[[sens_var]] <- factor(data_obs_dupl$obs_rep, 
                                      1:(nb_lvls + drop_main_model),
                                      levels(data_sensitivity[[sens_var]]))
  
  data_impact <- bind_rows(
    data_obs_dupl,
    data_sensitivity %>% filter(type == "No passport\n(counterfactual)")
  )
  
  impact_pt_estim <- compute_impact(data_impact, MODEL_END, var_group_by = c(sens_var))
  
  ## join pt and CI together
  impact_pt_estim <- impact_pt_estim %>% select(type:vax_cover, pop_rpdb, n_vacc_1dose)
  
  # create label (grouping is to avoid creating a space in the non-negative estimates)
  # case_when is to do CI only for main model
  impact_pt_estim <- impact_pt_estim %>% 
    group_by(across(all_of(sens_var))) %>% 
    mutate(
      impact_lab = sprintf("%s %% pts.", format(round(impact, 1), nsmall = 1))
    ) %>% 
    ungroup()
  
  return(impact_pt_estim)
}

# function to plot impact
plot_impact_sens <- function(data_impact, sens_var_name, plt_name, col_name = NULL){
  ## create variable labels
  data_impact <- data_impact %>%
    mutate(age_updated = case_when(
      # add dashes or '+' sign
      age != "60_" ~ gsub("_", "\u2013", age),
      age == "60_" ~ gsub("_", "+", age)
      )
    )
  
  ggplot(data_impact, aes(x = age_updated, col = get(sens_var_name))) +
    # null effect
    geom_hline(yintercept = 0, linetype = "dashed") +
    # impact
    geom_pointrange(aes(y = impact, ymin = impact_lci, ymax = impact_uci),
                    position = position_dodge(width = 0.25), fatten = 1.6) +
    
    scale_colour_viridis_d(end = 0.8) +
    # aesthetics
    # coord_cartesian(ylim = impact_limits) +
    coord_cartesian(ylim = c(-0.5, 3.5)) +
    # labels
    labs(x = "Age group", y = "Vaccine passport impact\n(in percentage points)",
         col = col_name, title = plt_name) +
    theme(
      # resize title tag
      plot.title = element_text(size = 7),
      
      # make axis and facet text legible
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      
      # make legend key bigger and text legible
      legend.key.width = unit(0.7, "line"),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.position = "bottom"
    )
}

which_dataset <- "passport"

for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  DO_CMA <- FALSE
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  # input directories for data
  path_sens <- sprintf("../vaccine-passport-data/out/its-sensitivity-%s", PROVINCE)
  path_its <- sprintf("../vaccine-passport-data/out/its-fit-%s", PROVINCE)
  path_boot <- sprintf("../vaccine-passport-data/out/its-boot-%s", PROVINCE)
  path_ctfl <- sprintf("../vaccine-passport-data/out/its-sensitivity-%s/ctfl-coverage", PROVINCE)
  path_obs <- sprintf("../vaccine-passport-data/out/observed-%s", PROVINCE)
  city_suffix <- ""
  
  # Load data ----
  ## Observed passport impact ----
  ## observed
  data_observed <- load_results(path_obs, "age", "observed")
  data_observed <- data_observed %>% filter(date_wk_end >= "2021-06-01")
  
  # compute coverage and rate
  data_observed <- data_observed %>% 
    group_by(age) %>% 
    mutate(rate_1dose = n_vacc_1dose - lag(n_vacc_1dose),
           rate_plt = (rate_1dose / n_unvax) * 10^5,
           cv_dose = (n_vacc_1dose / pop_rpdb) * 100) %>% 
    ungroup()
  
  ## ITS point estimate
  data_its <- load_results(path_its, "age", "its")
  
  ## bootstrap estimates
  data_boot <- load_results(path_boot, "age", "bootstrap", R = 1000)
  
  # change labels to match nb of lines in legend
  data_its$type <- factor(data_its$type, 
                          levels = c("Passport", "No passport\n(counterfactual)"),
                          labels = c("Passport\n", "No passport\n(counterfactual)"))
  data_boot$type <- factor(data_boot$type, 
                           levels = c("Passport", "No passport\n(counterfactual)"),
                           labels = c("Passport\n", "No passport\n(counterfactual)"))
  
  # get quantiles for rate and coverage
  boot_ci <- compute_boot_ci(data_boot, "age")
  
  ## impact (point estimate and CI)
  # compare to observed
  data_impact <- bind_rows(data_observed %>% mutate(type = "A_observed"),
                           data_its %>% filter(type == "No passport\n(counterfactual)"))
  
  impact_pt <- compute_impact(data_impact, MODEL_END, var_group_by = "age")
  impact_ci <- compute_impact_ci(data_boot, data_observed, MODEL_END, var_group_by = "age")
  
  ## join point estimates and CI together
  impact_its <- full_join(impact_pt %>% select(age:pop_rpdb),
                          impact_ci %>% select(age:impact_uci),
                          by = c("age", "date_wk_end"))
  
  # create text label for impact (with CI)
  # create label (grouping is to avoid creating a space in the non-negative estimates)
  impact_its <- impact_its %>% 
    group_by(age) %>% 
    mutate(impact_lab = sprintf("%s %% pts.\n(%s\u2013%s)", 
                                format(round(impact, 1), nsmall = 1), 
                                format(round(impact_lci, 1), nsmall = 1), 
                                format(round(impact_uci, 1), nsmall = 1))) %>% 
    ungroup()
  
  ## Sensitivity analyses ----
  ### 1. Alternative start dates ----
  # load data
  data_ts_start <- load_data_sensitivity(
    sprintf("%s/predicted_start_of_timeseries_fitted_values.csv", path_sens),
    original_levels = as.Date(c("2021-07-03", "2021-06-26", "2021-07-10")), 
    new_lvl_names = c("July 7th\n(main model)", "June 26th\n(-1 week)", "July 10th\n(+1 week)"),
    sens_var_original_name = "start_ts_date", 
    sens_var_new_name = "start_ts"
  )
  
  # compute impact
  impact_ts_start <- compute_impact_sens(data_ts_start, data_observed, "start_ts")
  
  impact_ts_start <- bind_rows(
    impact_its %>% mutate(start_ts = factor("July 7th\n(main model)")),
    impact_ts_start
  )
  impact_ts_start <- impact_ts_start %>% 
    select(start_ts, date_wk_end, age, vax_cover, impact, impact_lab,
           pop_rpdb, n_vacc_1dose, impact_lci, impact_uci)
  
  ### 2. Alternative lengths of passport impact ----
  # load data
  data_alt_length <- load_data_sensitivity(
    sprintf("%s/predicted_length_of_impact_fitted_values.csv", path_sens),
    original_levels = PASSPORT_END + c(0, -7, 7),
    new_lvl_names = c("6 weeks\n(main model)", "5 weeks", "7 weeks"),
    sens_var_original_name = "stop_date", 
    sens_var_new_name = "pass_length"
  )
  
  # compute impact
  impact_alt_length <- compute_impact_sens(data_alt_length, data_observed, "pass_length")
  
  impact_alt_length <- bind_rows(
    impact_its %>% mutate(pass_length = factor("6 weeks\n(main model)")),
    impact_alt_length
  )
  impact_alt_length <- impact_alt_length %>% 
    select(pass_length, date_wk_end, age, vax_cover, impact, impact_lab,
           pop_rpdb, n_vacc_1dose, impact_lci, impact_uci)
  
  ### 3. Alternative model specifications ----
  # load data
  data_altmod <- load_data_sensitivity(
    sprintf("%s/predicted_alt_model_specs_fitted_values.csv", path_sens),
    original_levels = c("final_model", "alt_model", "linear_model"),
    new_lvl_names = case_when(PROVINCE == "qc" ~ c("Natural spline\n(main model)", "Different slope and intercept\nafter impact period", "Log-linear\nmodel"),
                              PROVINCE == "on" ~ c("Natural spline\n(main model)", "Different slope and\nintercept for July", "Log-linear\nmodel")),
    sens_var_original_name = "model_specification", 
    sens_var_new_name = "model_spec"
  )
  
  # compute impact
  impact_altmod <- compute_impact_sens(data_altmod, data_observed, "model_spec")
  
  impact_altmod <- bind_rows(
    impact_its %>% mutate(model_spec = factor("Natural spline\n(main model)")),
    impact_altmod
  )
  impact_altmod <- impact_altmod %>% 
    select(model_spec, date_wk_end, age, vax_cover, impact, impact_lab,
           pop_rpdb, n_vacc_1dose, impact_lci, impact_uci)
  
  ## Plot ----
  # change start of timeseries
  plot_impact_sens(impact_ts_start, "start_ts", sprintf("%s) Change start of time series",
                                                        ifelse(PROVINCE == "qc", "A", "B")))
  ggsave(sprintf("%s/fig_R5_%s.1_sensitivity_ts_start.png", fig_path, PROVINCE),
         device = "png", dpi = 320,
         height = 5, width = 8, units = "cm")
  
  # change in length of impact period
  plot_impact_sens(impact_alt_length, "pass_length", sprintf("%s) Change length of the impact period",
                                                             ifelse(PROVINCE == "qc", "C", "D")))
  ggsave(sprintf("%s/fig_R5_%s.2_sensitivity_impact_length.png", fig_path, PROVINCE),
         device = "png", dpi = 320,
         height = 5, width = 8, units = "cm")
  
  # change in model specification for the vaccinations ~ time relationship
  plot_impact_sens(impact_altmod, "model_spec", sprintf("%s) Change model specification",
                                                        ifelse(PROVINCE == "qc", "E", "F")))
  ggsave(sprintf("%s/fig_R5_%s.3_sensitivity_model_spec.png", fig_path, PROVINCE),
         device = "png", dpi = 320,
         height = 5, width = 8, units = "cm")
}
