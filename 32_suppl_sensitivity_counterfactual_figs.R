# Libraries & functions ----
library(tidyverse)
library(lubridate)
library(gridExtra)
library(MetBrewer)

source("./src/utils_load_data.R")
source("./src/utils_its_bootstrap.R")
source("./src/utils_sensitivity.R")

source("./src/plot_sensitivity.R")
source("./src/plot_setup.R")

# output directory for figures
fig_path <- sprintf("./fig")

# province and CMA for troubleshooting
PROVINCE <- "qc"
CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")

# lists to store plots
plt_age_ls <- vector("list", 2)
names(plt_age_ls) <- c("qc", "on")

plt_income_ls <- plt_age_ls
plt_vismin_ls <- plt_age_ls

for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  DO_CMA <- FALSE
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  # input directories for data
  path_its <- sprintf("./out/its-fit-%s", PROVINCE)
  path_boot <- sprintf("./out/its-boot-%s", PROVINCE)
  path_ctfl <- sprintf("./out/its-sensitivity-%s/ctfl-coverage", PROVINCE)
  path_obs <- sprintf("./out/observed-%s", PROVINCE)
  city_suffix <- ""
  
  # Load data ----
  ## Observed passport impact ----
  ## observed
  data_observed <- load_merge_data(path_obs, "observed")
  
  # keep only June 2021 onwards
  data_observed <- data_observed %>% filter(date_wk_end >= "2021-06-01")
  
  # compute coverage and rate
  data_observed <- data_observed %>% 
    mutate(rate_plt = (rate_1dose / n_unvax) * 10^5,
           cv_dose = (n_vacc_1dose / pop_rpdb) * 100)
  
  ## ITS point estimates
  data_its <- load_merge_data(path_its, "its")
  
  ## bootstrap estimates
  data_boot <- load_merge_data(path_boot, "bootstrap", R = 1000)
  
  ## impact (point estimate)
  # compare to observed
  data_impact <- bind_rows(data_observed %>% mutate(type = "a_observed"),
                           data_its %>% filter(type == "No passport\n(counterfactual)"))
  
  impact_pt_obs <- compute_impact(data_impact, MODEL_END, 
                                  var_group_by = c("strat_var", "strat_lvl", "original_lvl"))
  
  ## impact (confidence interval)
  impact_ci <- compute_impact_ci(data_boot, data_observed, MODEL_END,
                                 var_group_by = c("strat_var", "strat_lvl",
                                                  "original_lvl"))
  
  ## Counterfactuals ----
  ## ITS point estimates
  # age
  data_ctfl_age <- load_ctfl_cov_results(path_ctfl, "age")
  impact_age <- compute_impact(data_ctfl_age, MODEL_END,
                               var_group_by = c("vaxcov_grp_set_to", "age"))
  # income
  data_ctfl_inc <- load_ctfl_cov_results(path_ctfl, "income")
  impact_income <- compute_impact(data_ctfl_inc, MODEL_END,
                                  var_group_by = c("vaxcov_grp_set_to", "quin_income"))
  
  # vismin
  data_ctfl_vis <- load_ctfl_cov_results(path_ctfl, "vismin")
  impact_vismin <- compute_impact(data_ctfl_vis, MODEL_END,
                                  var_group_by = c("vaxcov_grp_set_to", "quin_vismin"))
  
  # add the respective observed impacts to the data
  impact_age <- join_impact_obs_ctfl(impact_age, impact_pt_obs, "Age", "age")
  impact_income <- join_impact_obs_ctfl(impact_income, impact_pt_obs,
                                        "Income quintile", "quin_income")
  impact_vismin <- join_impact_obs_ctfl(impact_vismin, impact_pt_obs,
                                        "Visible minority quintile", "quin_vismin")
  
  # Plot estimates (individual) ----
  # parameters for plotting the absolute and relative (to unvaccinated before 
  # the announcement) impact of vaccine passport
  # change plot label based on location
  plt_title <- case_when(PROVINCE == "qc" ~ "A) Québec",
                         PROVINCE == "on" ~ "B) Ontario")
  
  ## Clean up labels ----
  lvls_vax <- c("Observed", "49.9_under", "50_59.9_", "60_69.9_", 
                "70_79.9_", "80_89.9_", "90_over")
  labs_vax <- c("Observed", "<50%", "50%\u2013<60%", "60%\u2013<70%", 
                "70%\u2013<80%", "80%\u2013<90%", ">90%")
  
  # age
  impact_age <- impact_age %>% 
    mutate(
      vaxcov_lab = factor(
        vaxcov_grp_set_to,
        levels = lvls_vax,
        labels = labs_vax
      ),
      age = case_when(
        age == "60_" ~ gsub("_", "+", age),
        T ~ gsub("_", "\u2013", age)
      )
    )
  impact_ci <- impact_ci %>% 
    mutate(
      age = case_when(
        strat_var != "Age" ~ NA_character_,
        original_lvl == "60_" ~ gsub("_", "+", original_lvl),
        T ~ gsub("_", "\u2013", original_lvl)
      )
    )
  
  # SDOH
  impact_income <- impact_income %>% 
    mutate(
      vaxcov_lab = factor(
        vaxcov_grp_set_to,
        levels = lvls_vax,
        labels = labs_vax
      )
    )
  impact_vismin <- impact_vismin %>% 
    mutate(
      vaxcov_lab = factor(
        vaxcov_grp_set_to,
        levels = lvls_vax,
        labels = labs_vax
      )
    )
  
  impact_ci <- impact_ci %>% 
    mutate(
      quin_income = case_when(
        strat_var != "Income quintile" ~ NA_character_,
        T ~ original_lvl
      ),
      quin_vismin = case_when(
        strat_var != "Visible minority quintile" ~ NA_character_,
        T ~ original_lvl
      )
    )
  
  ## Plot ----
  # list of quintile ranking names
  # corresponds to 1:5 for income
  # and to 5:1 for visible minority
  sdoh_labels <- c("Lowest", "2nd\nlowest", "Middle", "2nd\nhighest", "Highest")
  
  # age
  plt_age_ls[[PROVINCE]] <- plot_ctfl_cov_points(
    impact_age, impact_ci,
    "age", "Age", "Age", plt_title
  ) +
    scale_colour_viridis_d(end = 0.8, guide = guide_legend(reverse = TRUE))
  
  # income
  plt_income_ls[[PROVINCE]] <- plot_ctfl_cov_points(
    impact_income, impact_ci,
    "quin_income", "Income quintile", "Income\nquintile", plt_title
  ) +
    scale_colour_viridis_d(end = 0.8, labels = sdoh_labels)
  
  # visible minority
  plt_vismin_ls[[PROVINCE]] <- plot_ctfl_cov_points(
    impact_vismin, impact_ci,
    "quin_vismin", "Visible minority quintile", "Proportion\nracialized quintile", plt_title
  ) +
    scale_colour_viridis_d(end = 0.8, labels = rev(sdoh_labels))
}

# Plot final figures ----

## Age ----
# format plot
theme_age_common <- theme(
  legend.position = "none",
  # titles text
  plot.title = element_text(size = 8),
  legend.title = element_text(size = 8),
  axis.title.y = element_blank(),
  axis.title.x = element_text(size = 8),
  # smaller text
  axis.text.y = element_text(size = 6),
  axis.text.x = element_text(size = 6, hjust = 1, angle = 45),
  legend.text = element_text(size = 6)
)

plt_age_ls[["qc"]] <- plt_age_ls[["qc"]] + theme_age_common + theme(plot.background = element_blank())
plt_age_ls[["on"]] <- plt_age_ls[["on"]] + theme_age_common

# extract legend
p_legend_age <- cowplot::get_legend(
  plt_age_ls[["qc"]] +
    theme(legend.position = "right",
          legend.key.width = unit(0.4, "line"),
          legend.key.height = unit(0.8, "line"))
)

p_yaxis_age <- grid::textGrob(
  "Vaccine passport impact\n(in percentage points)",
  gp = grid::gpar(fontsize = 8, lineheight = 0.8), rot = 90
)

# plot together
png(sprintf("%s/fig_S3_ctfl_cover_age.png", fig_path), type = "cairo-png",
    width = 17, height = 6, units = "cm", res = 320)
grid.arrange(p_yaxis_age, plt_age_ls[["qc"]], plt_age_ls[["on"]], p_legend_age,
             nrow = 1, widths = unit(c(0.4, 14.5/2, 14.5/2, 1), c("cm")))
dev.off()

## Income and visible minority quintile ----
# format plot
theme_sdoh_common <- theme(
  legend.position = "none",
  # titles text
  plot.title = element_text(size = 8),
  legend.title = element_text(size = 8),
  axis.title.y = element_blank(),
  axis.title.x = element_text(size = 8),
  # smaller text
  axis.text.y = element_text(size = 6),
  axis.text.x = element_text(size = 6, hjust = 1, angle = 45),
  legend.text = element_text(size = 6),
  plot.background = element_blank()
)

# income
plt_income_ls[["qc"]] <- plt_income_ls[["qc"]] + theme_sdoh_common
plt_income_ls[["on"]] <- plt_income_ls[["on"]] + theme_sdoh_common

# visible minority
plt_vismin_ls[["qc"]] <- plt_vismin_ls[["qc"]] + theme_sdoh_common
plt_vismin_ls[["on"]] <- plt_vismin_ls[["on"]] + theme_sdoh_common

# extract legends
p_legend_income <- cowplot::get_legend(
  plt_income_ls[["qc"]] +
    theme(legend.position = "bottom",
          legend.key.width = unit(0.4, "line"),
          legend.key.height = unit(0.8, "line"))
)

p_legend_vismin <- cowplot::get_legend(
  plt_vismin_ls[["qc"]] +
    theme(legend.position = "bottom",
          legend.key.width = unit(0.4, "line"),
          legend.key.height = unit(0.8, "line"))
)

p_yaxis_sdoh <- grid::textGrob(
  "Vaccine passport impact (in percentage points)",
  gp = grid::gpar(fontsize = 8), rot = 90
)

# plot together
png(sprintf("%s/fig_S5_ctfl_cover_sdoh.png", fig_path), type = "cairo-png",
    width = 17, height = 6*2+1, units = "cm", res = 320)
grid.arrange(
  p_yaxis_sdoh, p_yaxis_sdoh,
  plt_income_ls[["qc"]], plt_vismin_ls[["qc"]],
  plt_income_ls[["on"]], plt_vismin_ls[["on"]],
  p_legend_income, p_legend_vismin,
  # layout
  layout_matrix = matrix(c(1, 3, 2, 4,
                           1, 5, 2, 6,
                           7, 7, 8, 8),
                         ncol = 4, byrow = T),
  # sizes
  widths = unit(c(0.5, 16.3/2, 0.3, 16.3/2), c("cm")),
  heights = unit(c(6, 6, 1), c("cm"))
)
dev.off()


