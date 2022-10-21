# Libraries ----
library(tidyverse)
library(lubridate)
library(grid)
library(gridExtra)
library(MetBrewer)

source("./src/utils_its_bootstrap.R")
source("./src/utils_sensitivity.R")
source("./src/utils_load_data.R")

source("./src/plot_main.R")
source("./src/plot_sensitivity.R")
source("./src/plot_setup.R")

## store plots in nested list, 
#   1st level = province
#   2nd level = sensitivity analysis
plt_cover_ls <- vector("list", 2)
plt_rate_ls <- vector("list", 2)

names(plt_cover_ls) <- c("qc", "on")
names(plt_rate_ls) <- c("qc", "on")

plt_ls <- vector("list", 3)
names(plt_ls) <- c("ts_start", "alt_length", "altmod")

for(prov in c("qc", "on")){
  plt_cover_ls[[prov]] <- plt_ls
  plt_rate_ls[[prov]] <- plt_ls
}
rm(prov, plt_ls)

# figure directory
fig_path <- sprintf("./fig")

for(PROVINCE in c("qc", "on")){
  # policy dates
  CMA <- "none"
  DO_CMA <- FALSE
  source("./03_setup_policy_dates.R")
  
  # store each province's in separate variables
  if(PROVINCE == "qc"){
    MODEL_END_QC <- MODEL_END
  } else if(PROVINCE == "on"){
    MODEL_END_ON <- MODEL_END
  }
  
  # directories
  path_sens <- sprintf("./out/its-sensitivity-%s", PROVINCE)
  path_obs <- sprintf("./out/observed-%s", PROVINCE)
  path_its <- sprintf("./out/its-fit-%s", PROVINCE)
  path_boot <- sprintf("./out/its-boot-%s", PROVINCE)
  
  # tibble for policy milestones labels
  policy_txt <- tibble(event = c("Announ.", "Implem."),
                       dates = POLICY,
                       start_ts = factor("main model"),
                       model_spec = factor("main model"),
                       pass_length = factor("main model"),
                       age = "12_17",
                       cv_pos = c(0, 0),
                       rate_pos = c(7950, 7000))
  
  # where to place announcement and implementation text
  if(PROVINCE == "qc"){
    policy_txt$cv_pos <- c(97.5, 93)
  } else {
    policy_txt$cv_pos <- c(66.5, 62)
  }
  
  # Load obserbed and main model results ----
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
  data_impact <- bind_rows(data_observed %>% mutate(type = "a_observed"),
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
  
  ## Age variable levels ----
  ## create variable labels
  variable_txt <- impact_its %>% 
    select(age)
  
  variable_txt <- variable_txt %>%
    mutate(var_text = case_when(
      # add dashes or '+' sign
      age != "60_" ~ gsub("_", "\u2013", age),
      age == "60_" ~ gsub("_", "+", age)
      )
    )
  
  # Sensitivity analyses ----
  # plot setup
  date_breaks <- seq(min(data_its$date_wk_end), max(data_its$date_wk_end), by = 28)
  
  ## 1. Alternative start dates ----
  # load data
  data_ts_start <- load_data_sensitivity(
    sprintf("%s/predicted_start_of_timeseries_fitted_values.csv", path_sens),
    data_its = data_its,
    original_levels = as.Date(c("2021-07-03", "2021-06-26", "2021-07-10")), 
    new_lvl_names = c("main model", "minus 1 week", "plus 1 week"),
    sens_var_original_name = "start_ts_date", 
    sens_var_new_name = "start_ts"
  )
  
  # compute impact
  impact_ts_start <- compute_impact_sens(data_ts_start, data_observed, "start_ts")
  
  impact_ts_start <- bind_rows(
    impact_its %>% mutate(start_ts = factor("main model")),
    impact_ts_start
  )
  impact_ts_start <- impact_ts_start %>% 
    select(start_ts, date_wk_end, age, vax_cover, impact, impact_lab,
           pop_rpdb, n_vacc_1dose, impact_lci, impact_uci)
  
  ## plot coverage
  plt_cover_ls[[PROVINCE]][["ts_start"]] <- plot_sensitivity_cover(
    data_ts_start, impact_ts_start, "start_ts",
    boot_ci, data_observed, variable_txt
  )
  
  ## 2. Alternative lengths of passport impact ----
  # load data
  data_alt_length <- load_data_sensitivity(
    sprintf("%s/predicted_length_of_impact_fitted_values.csv", path_sens),
    data_its = data_its,
    # original_levels = c("6 weeks", "5 weeks", "7 weeks"),
    original_levels = PASSPORT_END + c(-7, 0, 7),
    new_lvl_names = c("main model", "5 weeks", "7 weeks"),
    sens_var_original_name = "stop_date", 
    sens_var_new_name = "pass_length"
  )
  
  # compute impact
  impact_alt_length <- compute_impact_sens(data_alt_length, data_observed, "pass_length")
  
  impact_alt_length <- bind_rows(
    impact_its %>% mutate(pass_length = factor("main model")),
    impact_alt_length
  )
  impact_alt_length <- impact_alt_length %>% 
    select(pass_length, date_wk_end, age, vax_cover, impact, impact_lab,
           pop_rpdb, n_vacc_1dose, impact_lci, impact_uci)
  
  ## plot coverage
  plt_cover_ls[[PROVINCE]][["alt_length"]] <- plot_sensitivity_cover(
    data_alt_length, impact_alt_length, "pass_length",
    boot_ci, data_observed, variable_txt
  )
  
  ## plot rate
  plt_rate_ls[[PROVINCE]][["alt_length"]] <- plot_sensitivity_rate(
    data_alt_length, impact_alt_length, "pass_length",
    boot_ci, data_observed, variable_txt
  )
  
  ## 3. Alternative model specifications ----
  # load data
  data_altmod <- load_data_sensitivity(
    sprintf("%s/predicted_alt_model_specs_fitted_values.csv", path_sens),
    data_its = data_its,
    original_levels = c("final_model", "alt_model", "linear_model"),
    new_lvl_names = c("main model", "alt_model", "linear_model"),
    sens_var_original_name = "model_specification", 
    sens_var_new_name = "model_spec"
  )
  
  # compute impact
  impact_altmod <- compute_impact_sens(data_altmod, data_observed, "model_spec")
  
  impact_altmod <- bind_rows(
    impact_its %>% mutate(model_spec = factor("main model")),
    impact_altmod
  )
  impact_altmod <- impact_altmod %>% 
    select(model_spec, date_wk_end, age, vax_cover, impact, impact_lab,
           pop_rpdb, n_vacc_1dose, impact_lci, impact_uci)
  
  ## Plot coverage
  plt_cover_ls[[PROVINCE]][["altmod"]] <- plot_sensitivity_cover(
    data_altmod, impact_altmod, "model_spec",
    boot_ci, data_observed, variable_txt
  )
  
  ## plot rate
  plt_rate_ls[[PROVINCE]][["altmod"]] <- plot_sensitivity_rate(
    data_altmod, impact_altmod, "model_spec",
    boot_ci, data_observed, variable_txt
  )
}

# Save plots ----
plt_legend <- cowplot::get_legend(
  plt_rate_ls[["qc"]][["alt_length"]] + 
    theme(legend.position = "bottom")
)

## save alternative start dates (coverage only)
png(sprintf("%s/fig_S6_sensitivity_start_of_timeseries_cover.png", fig_path),
    width = 40, height = 47, units = "cm", res = 320,
    type = "cairo-png")
grid.arrange(plt_cover_ls[["qc"]][["ts_start"]],
             plt_cover_ls[["on"]][["ts_start"]],
             plt_legend,
             ncol = 1,  heights = unit(c(45/2, 45/2, 2), c("cm")))
draw_sensitivity_labs("Start of timeseries: July 3rd (main model)",
                      "Start of timeseries: June 26th",
                      "Start of timeseries: July 10th",
                      y_pos = 0.05)
dev.off()

## save alternative lengths of passport impact
# coverage
png(sprintf("%s/fig_S8_sensitivity_length_of_impact_cover.png", fig_path),
    width = 40, height = 47, units = "cm", res = 320,
    type = "cairo-png")
grid.arrange(plt_cover_ls[["qc"]][["alt_length"]],
             plt_cover_ls[["on"]][["alt_length"]],
             plt_legend,
             ncol = 1,  heights = unit(c(45/2, 45/2, 2), c("cm")))
draw_sensitivity_labs("The effect of the vaccine passport lasts 6 weeks (main model)",
                      "The effect of the vaccine passport lasts 5 weeks",
                      "The effect of the vaccine passport lasts 7 weeks",
                      y_pos = 0.05)
dev.off()

# rate
png(sprintf("%s/fig_S7_sensitivity_length_of_impact_rate.png", fig_path),
    width = 40, height = 47, units = "cm", res = 320,
    type = "cairo-png")
grid.arrange(plt_rate_ls[["qc"]][["alt_length"]],
             plt_rate_ls[["on"]][["alt_length"]],
             plt_legend,
             ncol = 1,  heights = unit(c(45/2, 45/2, 2), c("cm")))
draw_sensitivity_labs("The effect of the vaccine passport lasts 6 weeks (main model)",
                      "The effect of the vaccine passport lasts 5 weeks",
                      "The effect of the vaccine passport lasts 7 weeks",
                      y_pos = 0.072)
dev.off()

## save alternative model specifications
# coverage
png(sprintf("%s/fig_S10_sensitivity_alt_model_specs_cover.png", fig_path),
    width = 40, height = 47, units = "cm", res = 320,
    type = "cairo-png")
grid.arrange(plt_cover_ls[["qc"]][["altmod"]],
             plt_cover_ls[["on"]][["altmod"]],
             plt_legend,
             ncol = 1,  heights = unit(c(45/2, 45/2, 2), c("cm")))
draw_sensitivity_labs("Vaccination rate modeled with a natural cubic spline (main model)",
                      "Vaccination rate modeled with a different intercept and slope after the vaccine passport impact period",
                      "Vaccination rate modeled as log-linear",
                      lab_e = "Vaccination rate modeled with a quadratic term and a different intercept and slope for July",
                      y_pos = 0.05)
dev.off()

# rate
png(sprintf("%s/fig_S9_sensitivity_alt_model_specs_rate.png", fig_path),
    width = 40, height = 47, units = "cm", res = 320,
    type = "cairo-png")
grid.arrange(plt_rate_ls[["qc"]][["altmod"]],
             plt_rate_ls[["on"]][["altmod"]],
             plt_legend,
             ncol = 1,  heights = unit(c(45/2, 45/2, 2), c("cm")))
draw_sensitivity_labs("Vaccination rate modeled with a natural cubic spline (main model)",
                      "Vaccination rate modeled with a different intercept and slope after the vaccine passport impact period",
                      "Vaccination rate modeled as log-linear",
                      lab_e = "Vaccination rate modeled with a quadratic term and a different intercept and slope for July",
                      y_pos = 0.072)
dev.off()
