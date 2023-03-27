# Libraries and setup ----
library(tidyverse)
library(data.table)
library(lubridate)
library(MetBrewer)
library(gridExtra)

source("./src/utils_load_data.R")
source("./src/plot_main.R")
source("./src/plot_sensitivity.R")
source("./src/plot_setup.R")

# which dataset to use
which_dataset <- "descriptive"

# province and CMA for troubleshooting
PROVINCE <- "qc"
CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")


# list to store individual plots
list_plots <- vector("list", 2)
names(list_plots) <- c("qc", "on")

list_plots$qc <- vector("list", 3)
list_plots$on <- vector("list", 3)

for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  DO_CMA <- FALSE
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  # Load vaccination coverage data ----
  source("./04_setup_load_data.R", echo = TRUE)
  data_coverage <- data_coverage[date_wk_end >= "2021-07-03" & date_wk_end <= MODEL_END]
  data_passport_da <- data_passport_da[date_wk_end >= "2021-07-03" & date_wk_end <= MODEL_END]
  
  # Vaccination rate ----
  # set dates for x axis (single-panel plots)
  date_breaks <- seq.Date(min(data_coverage$date_wk_end), max(data_coverage$date_wk_end)+7, by = 14)
  
  ## Age-stratified distribution ----
  data_mean_age <- data_coverage[
    , .(rate_wkly = (sum(rate_1dose) / sum(n_unvax))*10^5,
        nb_rate_0 = sum(rate_1dose == 0)),
    by = .(age, date_wk_end)]
  
  # plot distribution, mean and variance
  ## create variable labels
  variable_txt <- data_coverage %>% 
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
  
  var_labels <- variable_txt$var_text
  names(var_labels) <- variable_txt$age
  
  list_plots[[PROVINCE]][[1]] <- plot_distr_rate(
    data_coverage, data_mean_age,
    "age", label_fn = labeller(age = var_labels)
  )
  
  ## Income-stratified distribution ----
  # statistics
  data_mean_income <- data_passport_da[
    , .(rate_wkly = (sum(rate_1dose) / sum(n_unvax))*10^5,
        nb_rate_0 = sum(rate_1dose == 0)),
    by = .(quin_income, date_wk_end)]
  
  # plot distribution, mean and variance
  ## create variable labels
  variable_txt <- data_passport_da %>% 
    select(quin_income) %>% 
    unique()
  
  variable_txt <- variable_txt %>%
    mutate(var_text = c("Lowest", "2nd lowest", "Middle",
                        "2nd highest", "Highest"))
  
  var_labels <- variable_txt$var_text
  names(var_labels) <- sort(variable_txt$quin_income)
  
  list_plots[[PROVINCE]][[2]] <- plot_distr_rate(
    data_passport_da, data_mean_income,
    "quin_income", label_fn = labeller(quin_income = var_labels)
  )
  
  ## Visible minority-stratified distribution ----
  # statistics
  data_mean_vismin <- data_passport_da[
    , .(rate_wkly = (sum(rate_1dose) / sum(n_unvax))*10^5,
        nb_rate_0 = sum(rate_1dose == 0)),
    by = .(quin_vismin, date_wk_end)]
  
  # plot distribution, mean and variance
  ## create variable labels
  variable_txt <- data_passport_da %>% 
    select(quin_vismin) %>% 
    unique()
  
  variable_txt <- variable_txt %>%
    mutate(var_text = c("Highest", "2nd highest", "Middle",
                        "2nd lowest", "Lowest"))
  
  var_labels <- variable_txt$var_text
  names(var_labels) <- sort(variable_txt$quin_vismin)
  
  list_plots[[PROVINCE]][[3]] <- plot_distr_rate(
    data_passport_da, data_mean_vismin,
    "quin_vismin", label_fn = labeller(quin_vismin = var_labels)
  )
}

rm(data_coverage, data_passport_da)

p_legend <- cowplot::get_legend(
  list_plots$qc[[1]] +
    theme(legend.position = "bottom",
          legend.background = element_blank(),
          legend.key = element_blank())
)

{
  # age
  png("./fig/fig_S11_distrib_rate_age.png",
      width = 17, height = 20, units = "cm", res = 600,
      type = "cairo-png")
  grid.arrange(
    list_plots$qc[[1]],
    list_plots$on[[1]],
    p_legend,
    
    heights = unit(c(9.8, 9.8, 0.4), c("cm")),
    ncol = 1
  )
  dev.off()
  
  # income
  png("./fig/fig_S12_distrib_rate_income.png",
      width = 17, height = 20, units = "cm", res = 600,
      type = "cairo-png")
  grid.arrange(
    list_plots$qc[[2]],
    list_plots$on[[2]],
    p_legend,

    heights = unit(c(9.8, 9.8, 0.4), c("cm")),
    ncol = 1
  )
  dev.off()

  # visible minority
  png("./fig/fig_S13_distrib_rate_vismin.png",
      width = 17, height = 20, units = "cm", res = 600,
      type = "cairo-png")
  grid.arrange(
    list_plots$qc[[3]],
    list_plots$on[[3]],
    p_legend,

    heights = unit(c(9.8, 9.8, 0.4), c("cm")),
    ncol = 1
  )
  dev.off()
}
