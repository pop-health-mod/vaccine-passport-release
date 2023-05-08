# Libraries & functions ----
library(tidyverse)
library(gridExtra)
library(grid)
library(data.table)
library(lubridate)
library(MetBrewer)

theme_set(theme_bw())

source("./src/utils_load_data.R")
source("./src/utils_its_bootstrap.R")

source("./src/plot_main.R")
source("./src/plot_sensitivity.R")
source("./src/plot_setup.R")

# output directory for figures
fig_path <- sprintf("./fig")

# done for QC only
PROVINCE <- "qc"
CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
DO_CMA <- F
which_dataset <- "passport"

source("./03_setup_policy_dates.R")

# Load data ----
path_sens <- sprintf("../vaccine-passport-data/out/its-sensitivity-%s", PROVINCE)
path_its <- sprintf("../vaccine-passport-data/out/its-fit-%s", PROVINCE)
path_boot <- sprintf("../vaccine-passport-data/out/its-boot-%s", PROVINCE)
path_obs <- sprintf("../vaccine-passport-data/out/observed-%s", PROVINCE)

## Main results ----
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

## Sensitivity (fit with log-cases variable) ----
data_sens <- read_csv(sprintf("%s/predicted_qc_fit_w_cases_fitted_values.csv", path_sens),
                      col_types = cols())

# reorder variables to output
data_sens <- data_sens %>% 
  select(type, age, date_wk_end:rate_pred_cumul)

# Plot results ----
date_breaks <- seq(min(data_its$date_wk_end), max(data_its$date_wk_end), by = 28)

# change variable labels
data_sens_cases <- bind_rows(
  data_its %>% mutate(model_spec = "main model"),
  data_sens %>% mutate(model_spec = "log(cases)")
)

# change labels to match nb of lines in legend
data_sens_cases$type <- factor(data_sens_cases$type,
                               levels = c("Passport", "No passport\r\n(counterfactual)", "No passport\n(counterfactual)"),
                               labels = c("Passport\n", "No passport\n(counterfactual)", "No passport\n(counterfactual)"))
data_boot_plt <- data_boot
data_boot_plt$type <- factor(data_boot_plt$type,
                             levels = c("Passport", "No passport\n(counterfactual)"),
                             labels = c("Passport\n", "No passport\n(counterfactual)"))

# get quantiles for rate and coverage
boot_ci <- compute_boot_ci(data_boot_plt, "age")

#### plot pars
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

## create variable labels
variable_txt <- data_sens_cases %>% 
  select(age) %>% 
  unique()

variable_txt <- variable_txt %>%
  mutate(
    var_text = case_when(
      # add dashes or '+' sign
      age != "60_" ~ gsub("_", "\u2013", age),
      age == "60_" ~ gsub("_", "+", age)
    )
  )

## plot rate
data_sens_cases$model_spec <- factor(data_sens_cases$model_spec,
                                     levels = c("main model", "log(cases)"))

p_sens_case <- plot_sensitivity_rate(
  data_sens_cases, "impact_logcases", "model_spec",
  boot_ci, data_observed, variable_txt
)

## Save plot ----
plt_legend <- cowplot::get_legend(
  p_sens_case + 
    theme(legend.position = "bottom")
)

# rate
# height of 47cm is for: 45/2 per province (3 plots each) + 2cm for labels
# 45/2 * 2/3
png("./fig/fig_S11_sensitivity_cases.png",
    width = 40, height = 17, units = "cm", res = 320,
    type = "cairo-png")
grid.arrange(p_sens_case, plt_legend,
             ncol = 1,  heights = unit(c(15, 2), c("cm")))
draw_sensitivity_labs_qc("Main model",
                         "Sensitivity (log-cases included in the model)",
                         x_pos = 0.072)
dev.off()

# Table results ----
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

# create text for table (with CI)
impact_its <- impact_its %>% 
  group_by(age) %>% 
  mutate(impact_txt = sprintf("%s (%s\u2013%s)", 
                              format(round(impact, 1), nsmall = 1), 
                              format(round(impact_lci, 1), nsmall = 1), 
                              format(round(impact_uci, 1), nsmall = 1))) %>% 
  ungroup()

impact_its

# sensitivity
data_impact_sens <- bind_rows(data_observed %>% mutate(type = "a_observed"),
                              data_sens %>% filter(type == "No passport\r\n(counterfactual)"))

impact_sens <- compute_impact(data_impact_sens, MODEL_END, var_group_by = "age")
impact_sens <- impact_sens %>% select(age:pop_rpdb)

# create text for table (with CI)
impact_sens <- impact_sens %>% 
  group_by(age) %>% 
  mutate(impact_txt = format(round(impact, 1), nsmall = 1)) %>% 
  ungroup()

impact_sens

impact_ttl <- bind_rows(
  impact_its %>% mutate(model_spec = "main model"),
  impact_sens %>% mutate(model_spec = "log(cases)")
) %>% 
  select(model_spec, age, date_wk_end, impact_txt)

write.csv(impact_ttl,
          "./out/manuscript-tables/table_S3_case_fit_impact_estimates.csv",
          row.names = FALSE)
