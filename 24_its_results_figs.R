# Libraries ----
library(tidyverse)
library(lubridate)
library(grid)
library(gridExtra)
library(MetBrewer)

source("./src/utils_its_bootstrap.R")
source("./src/utils_load_data.R")

source("./src/plot_main.R")
source("./src/plot_setup.R")

# output directory for figures
fig_path <- sprintf("./fig")

# PROVINCE <- "qc"
# CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
# DO_CMA <- FALSE

# lists to store plots
plt_cover_ls <- vector("list", 4)
plt_rate_ls <- vector("list", 4)

names(plt_cover_ls) <- c("qc", "on", "mtl", "tor")
names(plt_rate_ls) <- c("qc", "on", "mtl", "tor")

for(PROVINCE in c("qc", "on")){
  # Single-variable ITS results in single figure ----
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  # input directories for data
  for(DO_CMA in c(FALSE, TRUE)){
    if(DO_CMA){
      path_its <- sprintf("./out/its-fit-%s-%s", PROVINCE, CMA)
      path_boot <- sprintf("./out/its-boot-%s-%s", PROVINCE, CMA)
      path_obs <- sprintf("./out/observed-%s-%s", PROVINCE, CMA)
      
      city_suffix <- paste("_", CMA, sep = "")
    } else {
      path_its <- sprintf("./out/its-fit-%s", PROVINCE)
      path_boot <- sprintf("./out/its-boot-%s", PROVINCE)
      path_obs <- sprintf("./out/observed-%s", PROVINCE)
      
      city_suffix <- ""
    }
    
    theme_set(theme_bw())
    
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
    
    # compute rate and coverage, and get quantiles
    boot_ci <- compute_boot_ci(data_boot, c("strat_var", "strat_lvl", "original_lvl"),
                               numer = "rate_1dose", denom = "n_unvax")
    
    ## impact (point estimate)
    # compare to observed
    data_impact <- bind_rows(data_observed %>% mutate(type = "a_observed"),
                             data_its %>% filter(type == "No passport\n(counterfactual)"))
    
    impact_pt <- compute_impact(data_impact, MODEL_END,
                                var_group_by = c("strat_var", "strat_lvl", "original_lvl"))
    
    ## impact (confidence interval)
    impact_ci <- compute_impact_ci(data_boot, data_observed, MODEL_END,
                                   var_group_by = c("strat_var", "strat_lvl",
                                                    "original_lvl"))
    
    ## Tables for label text ----
    ## join pt and CI together
    impact_txt <- left_join(impact_pt %>% select(type:vax_cover),
                            impact_ci %>% select(strat_var:impact_uci),
                            by = c("strat_var", "strat_lvl", "original_lvl", "date_wk_end"))
    
    # create effect labels (grouping avoids space when CIs are <0)
    impact_txt <- impact_txt %>% 
      group_by(strat_var, strat_lvl, original_lvl) %>% 
      mutate(impact_lab = sprintf("%s %% pts.\n(%s\u2013%s)",
                                  format(round(impact, 1), nsmall = 1),
                                  format(round(impact_lci, 1), nsmall = 1),
                                  format(round(impact_uci, 1), nsmall = 1))) %>% 
      ungroup()
    
    ## create variable labels
    variable_txt <- impact_txt %>%
      select(strat_var:original_lvl)
    
    variable_txt <- variable_txt %>%
      mutate(var_text = case_when(
        # add dashes or '+' sign
        strat_var == "Age" & strat_lvl != 6 ~ gsub("_", "\u2013", original_lvl),
        strat_var == "Age" & strat_lvl == 6 ~ gsub("_", "+", original_lvl),
        # label lower ones (income is 1, vismin is 5)
        (grepl("Income", strat_var) & original_lvl == 1) |
          (grepl("Visible", strat_var) & original_lvl == 5) ~ "Lowest",
        (grepl("Income", strat_var) & original_lvl == 2) |
          (grepl("Visible", strat_var) & original_lvl == 4) ~ "2nd lowest",
        # label middle
        original_lvl == 3 ~ "Middle",
        # label higher ones (income is 5, vismin is 1)
        (grepl("Income", strat_var) & original_lvl == 5) |
          (grepl("Visible", strat_var) & original_lvl == 1) ~ "Highest",
        (grepl("Income", strat_var) & original_lvl == 4) |
          (grepl("Visible", strat_var) & original_lvl == 2) ~ "2nd highest"
        )
      )
    
    policy_txt <- tibble(event = c("Announ.", "Implem."),
                         dates = POLICY,
                         strat_var = variable_txt$strat_var[1],
                         strat_lvl = variable_txt$strat_lvl[1],
                         cv_pos = c(0, 0),
                         rate_pos = c(7950, 7000))
    
    # where to place announcement and implementation text
    if(PROVINCE == "qc"){
      policy_txt$cv_pos <- c(97.5, 93)
    } else {
      policy_txt$cv_pos <- c(66.5, 62)
    }

    rm(data_impact, data_boot)
    
    ## Plot coverage ----
    # setup
    date_breaks <- seq(min(data_its$date_wk_end), max(data_its$date_wk_end), by = 28)
    
    # also plot 18-29 under for Toronto (TEMPORAL FIX, RE-CHECK IF DENOMINATOR ISSUE RESOLVED)
    if(DO_CMA & CMA == "mtl"){
      grp_plot_under <- c("50_59", "60_")
    } else if(DO_CMA & CMA == "tor"){
      grp_plot_under <- c("18_29", "60_")
    } else {
      grp_plot_under <- c("60_")
    }
    
    # plot
    plt_cover <- plot_vaccine_cover(data_its,
                                    x_min = as.Date("2021-07-03")-5,
                                    x_max = max(data_its$date_wk_end)+5,
                                    y_min = 60, y_max = 100,
                                    colour_var = NULL,
                                    x_axis = FALSE) +
      # confidence interval
      geom_ribbon(data = boot_ci, aes(x = date_wk_end, ymin = cv_reg_lci, ymax = cv_reg_uci, fill = type),
                  col = NA, alpha = alpha_plt, inherit.aes = F) +
      # fitted data
      geom_line(aes(col = type)) +
      # observed data
      geom_point(data = data_observed, aes(shape = "Observed\ndata"), alpha = alpha_plt) +
      # impact (some age groups are plotted under line because of high coverage)
      geom_label(data = impact_txt %>% filter(!original_lvl %in% grp_plot_under),
                 aes(x = MODEL_END, y = vax_cover + 1, label = impact_lab),
                 hjust = 1, vjust = 0, size = 4.3, alpha = 0.8,
                 label.size = NA,
                 inherit.aes = FALSE) +
      geom_label(data = impact_txt %>% filter(original_lvl %in% grp_plot_under),
                 aes(x = MODEL_END, y = vax_cover - 10, label = impact_lab),
                 hjust = 1, vjust = 0, size = 4.3, alpha = 0.8,
                 label.size = NA,
                 inherit.aes = FALSE) +
      
      # variable level
      geom_label(data = variable_txt, aes(x = min(data_its$date_wk_end)-4, y = 99.7,
                                          label = var_text),
                 hjust = 0, vjust = 1, size = 5) +
      
      # policy announcements
      geom_label(data = policy_txt, aes(x = dates + 2, y = cv_pos, label = event),
                 hjust = 0, vjust = 0.5, size = 5,
                 label.size = NA, alpha = 0.8) +
      
      facet_grid(strat_var ~ strat_lvl) +
      
      # aesthetics
      scale_y_continuous(breaks = c(6:10 * 10)) +
      labs(title = "", x = "Date (end of epidemiological week)",
           shape = "", col = "Model fit", fill = "Model fit") +
      colour_model + fill_model +
      
      # format legend and text sizes put shape before colour and make legend key bigger
      guides_man_fig1.2 +
      theme_man_fig1.2
    
    ### remove empty panels
    gtable_cover <- format_ggtable_fig1.2(plt_cover)
    
    ## store
    plt_cover_ls[[ifelse(DO_CMA, CMA, PROVINCE)]] <- gtable_cover
    
    ## Plot rate ----
    # setup
    date_breaks <- seq(min(data_its$date_wk_end), max(data_its$date_wk_end), by = 28)
    
    # plot
    plt_rate <- plot_vaccine_rate(data_its,
                                  x_min = as.Date("2021-07-03")-5,
                                  x_max = max(data_its$date_wk_end)+5,
                                  y_max = 8400,
                                  colour_var = NULL,
                                  x_axis = FALSE) +
      # confidence interval
      geom_ribbon(data = boot_ci,
                  aes(x = date_wk_end, ymin = rate_plt_lci, ymax = rate_plt_uci, fill = type),
                  col = NA, alpha = alpha_plt, inherit.aes = F) +
      # fitted data
      geom_line(aes(col = type)) +
      # observed data
      geom_point(data = data_observed, aes(shape = "Observed\ndata"), alpha = alpha_plt) +
      
      # variable level
      geom_label(data = variable_txt, aes(x = min(data_its$date_wk_end)-4, y = 8350,
                                          label = var_text),
                 hjust = 0, vjust = 1, size = 5) +
      
      # policy announcements
      geom_label(data = policy_txt, aes(x = dates + 1, y = rate_pos, label = event),
                 hjust = 0, vjust = 0.5, size = 5,
                 label.size = NA, alpha = 0.8) +
      
      facet_grid(strat_var ~ strat_lvl) +
      
      # aesthetics
      scale_y_continuous(breaks = c(0:4 * 2000)) +
      labs(title = "", x = "Date (end of epidemiological week)",
           shape = "", col = "Model fit", fill = "Model fit") +
      colour_model + fill_model +
      
      # format legend and text sizes put shape before colour and make legend key bigger
      guides_man_fig1.2 +
      theme_man_fig1.2
    
    ### remove empty panels
    gtable_rate <- format_ggtable_fig1.2(plt_rate)
    
    ## store
    plt_rate_ls[[ifelse(DO_CMA, CMA, PROVINCE)]] <- gtable_rate
  }
}
rm(boot_ci, data_its, data_observed, gtable_cover, gtable_rate,
   impact_ci, impact_pt, impact_txt)

# set-up variables for names
title_age <- "age group"
title_inc <- "income quintile"
title_vis <- "proportion racialized quintile"

## Plot coverage ----
y_pos <- 0.05
## save
for(i in 1:2){
  # label for either province or CMA
  if(i == 1){
    title_qc <- "Québec"
    title_on <- "Ontario"
  } else {
    title_qc <- "Montréal"
    title_on <- "Toronto"
  }
  
  ## plot
  png(sprintf("%s/fig_%s2_cover.png", fig_path, ifelse(i == 1, "", "S")),
      width = 40, height = 45, units = "cm", res = 320,
      type = "cairo-png")
  # pdf(sprintf("%s/fig_%s2_cover.pdf", fig_path, ifelse(i == 1, "", "S")),
  #     width = 40/cm(1), height = 45/cm(1))
  grid.arrange(plt_cover_ls[[ifelse(i == 1, "qc", "mtl")]],
               plt_cover_ls[[ifelse(i == 1, "on", "tor")]])
  # Quebec labels
  grid.text(sprintf("A) %s, %s", title_qc, title_age), y_pos, 0.992, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("B) %s, %s", title_qc, title_inc), y_pos, 0.839, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("C) %s, %s", title_qc, title_vis), y_pos, 0.686, just = "left",
            gp = gpar(fontsize = 18))
  # Ontario labels
  grid.text(sprintf("D) %s, %s", title_on, title_age), y_pos, 0.492, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("E) %s, %s", title_on, title_inc), y_pos, 0.339, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("F) %s, %s", title_on, title_vis), y_pos, 0.186, just = "left",
            gp = gpar(fontsize = 18))
  dev.off()
}

## Plot rate ----
y_pos <- 0.072
## save
for(i in 1:2){
  # label for either province or CMA
  if(i == 1){
    title_qc <- "Québec"
    title_on <- "Ontario"
  } else {
    title_qc <- "Montréal"
    title_on <- "Toronto"
  }
  
  ## plot
  png(sprintf("%s/fig_%s1_rate.png", fig_path, ifelse(i == 1, "", "S")),
      width = 40, height = 45, units = "cm", res = 320,
      type = "cairo-png")
  # pdf(sprintf("%s/fig_%s1_rate.pdf", fig_path, ifelse(i == 1, "", "S")),
  #     width = 40/cm(1), height = 45/cm(1))
  grid.arrange(plt_rate_ls[[ifelse(i == 1, "qc", "mtl")]],
               plt_rate_ls[[ifelse(i == 1, "on", "tor")]])
  # Quebec labels
  grid.text(sprintf("A) %s, %s", title_qc, title_age), y_pos, 0.992, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("B) %s, %s", title_qc, title_inc), y_pos, 0.839, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("C) %s, %s", title_qc, title_vis), y_pos, 0.686, just = "left",
            gp = gpar(fontsize = 18))
  # Ontario labels
  grid.text(sprintf("D) %s, %s", title_on, title_age), y_pos, 0.492, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("E) %s, %s", title_on, title_inc), y_pos, 0.339, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("F) %s, %s", title_on, title_vis), y_pos, 0.186, just = "left",
            gp = gpar(fontsize = 18))
  dev.off()
}

# Heatmap for ITS interaction results ----
# lists to store plots
plt_heatmap_ls <- vector("list", 4)
names(plt_heatmap_ls) <- c("qc", "on", "mtl", "tor")

# create 2nd level of list storage
for(region in names(plt_heatmap_ls)){
  plt_heatmap_ls[[region]] <- vector("list", 2)
  names(plt_heatmap_ls[[region]]) <- c("income", "vismin")
}

## Generate individual heatmaps ----
for(PROVINCE in c("qc", "on")){
  # which CMA to use
  CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
  
  # passport dates
  source("./03_setup_policy_dates.R")
  
  # input directories for data
  for(DO_CMA in c(FALSE, TRUE)){
    if(DO_CMA){
      path_its <- sprintf("./out/its-fit-%s-%s", PROVINCE, CMA)
      path_boot <- sprintf("./out/its-boot-%s-%s", PROVINCE, CMA)
      path_obs <- sprintf("./out/observed-%s-%s", PROVINCE, CMA)
      
      city_suffix <- paste("_", CMA, sep = "")
    } else {
      path_its <- sprintf("./out/its-fit-%s", PROVINCE)
      path_boot <- sprintf("./out/its-boot-%s", PROVINCE)
      path_obs <- sprintf("./out/observed-%s", PROVINCE)
      
      city_suffix <- ""
    }
    
    theme_set(theme_minimal())
    
    ### Load data
    # observed results
    data_obs_income <- load_results(path_obs, "age.income", "observed")
    data_obs_vismin <- load_results(path_obs, "age.vismin", "observed")
    
    # ITS point estimates
    results_income <- load_results(path_its, "age.income", "its")
    results_vismin <- load_results(path_its, "age.vismin", "its")
    
    # bootstrap data
    data_boot_income <- load_results(path_boot, "age.income", "bootstrap", R = 1000)
    data_boot_vismin <- load_results(path_boot, "age.vismin", "bootstrap", R = 1000)
    
    ### compute impact of vaccine passport
    ## point estimates
    # income
    data_impact <- bind_rows(data_obs_income %>% mutate(type = "a_observed"),
                             results_income %>% filter(type == "No passport\n(counterfactual)"))
    impact_income <- compute_impact(data_impact, MODEL_END, c("age", "quin_income"))
    
    # visible minority
    data_impact <- bind_rows(data_obs_vismin %>% mutate(type = "a_observed"),
                             results_vismin %>% filter(type == "No passport\n(counterfactual)"))
    impact_vismin <- compute_impact(data_impact, MODEL_END, c("age", "quin_vismin"))
    rm(data_impact)
    
    ## confidence intervals
    impact_income_ci <- compute_impact_ci(data_boot_income, data_obs_income,
                                          MODEL_END, var_group_by = c("age", "quin_income"))
    impact_vismin_ci <- compute_impact_ci(data_boot_vismin, data_obs_vismin,
                                          MODEL_END, var_group_by = c("age", "quin_vismin"))
    
    ### Plot
    ## plot set-up
    heatmap_limits <- c(-1, 6)
    
    # change plot label based on location
    plt_title <- case_when(PROVINCE == "qc" & !DO_CMA ~ "A) Québec",
                           PROVINCE == "qc" & DO_CMA ~ "A) Montréal",
                           PROVINCE == "on" & !DO_CMA ~ "C) Ontario",
                           PROVINCE == "on" & DO_CMA ~ "C) Toronto")
    
    ## format labels 
    # age
    age_labs <- unique(impact_income$age)
    age_labs <- case_when(
      age_labs != "60_" ~ gsub("_", "\u2013", age_labs),
      age_labs == "60_" ~ gsub("_", "+", age_labs)
      )
    age_labs <- rev(age_labs)
    
    # SDOH
    inc_labs <- c("Lowest", "2nd\nlowest", "Middle", "2nd\nhighest", "Highest")
    vis_labs <- rev(inc_labs)
    
    ## plot income x age
    plt_heatmap_ls[[ifelse(DO_CMA, CMA, PROVINCE)]][["income"]] <-
      plot_heatmap_impact(impact_income,
                          x_name = "quin_income", x_title = "Income quintile",
                          plt_title = plt_title,
                          impact_limits = heatmap_limits,
                          show_ci = TRUE,
                          data_ci = impact_income_ci,
                          txt_size = 3.5) +
      scale_y_discrete(limits = rev, labels = age_labs) +
      scale_x_discrete(labels = inc_labs) +
      labs(y = "Age group")
    
    ## plot visible minority x age
    plt_title <- ifelse(PROVINCE == "qc",
                        gsub("A", "B", plt_title),
                        gsub("C", "D", plt_title))
    
    plt_heatmap_ls[[ifelse(DO_CMA, CMA, PROVINCE)]][["vismin"]] <-
      plot_heatmap_impact(impact_vismin,
                          x_name = "quin_vismin", x_title = "Proportion racialized quintile",
                          plt_title = plt_title,
                          impact_limits = heatmap_limits,
                          show_ci = TRUE,
                          data_ci = impact_vismin_ci,
                          txt_size = 3.5) +
      scale_y_discrete(limits = rev, labels = age_labs) +
      scale_x_discrete(labels = vis_labs) +
      labs(y = "Age group")
  }
}

## Plot together ----
# extract label
tmp <- ggplot_gtable(
  ggplot_build(
    plt_heatmap_ls[["qc"]][["income"]] +
      labs(fill = "Vaccine\npassport\nimpact\n(p.p.)") +
      theme(legend.position = "right")
  )
)
legend_pos <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
heatmap_legend <- tmp$grobs[[legend_pos]]

# plot
for(i in 1:2){
  plt <- grid.arrange(
    plt_heatmap_ls[[ifelse(i == 1, "qc", "mtl")]][["income"]],
    plt_heatmap_ls[[ifelse(i == 1, "qc", "mtl")]][["vismin"]],
    plt_heatmap_ls[[ifelse(i == 1, "on", "tor")]][["income"]],
    plt_heatmap_ls[[ifelse(i == 1, "on", "tor")]][["vismin"]],
    heatmap_legend,
    layout_matrix = matrix(c(1, 2, 5,
                             3, 4, 5),
                           nrow = 2, byrow = T),
    widths = c(20.5/2, 20.5/2, 2)
  )
  # add white background
  plt <- cowplot::ggdraw(plt) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  
  ggsave(sprintf("%s/fig_%s_interaction_age.sdoh.png", fig_path,
                 ifelse(i == 1, "3", "S4")),
         plot = plt,
         device = "png", units = "cm",
         width = 20.5+2, height = 24, dpi = 320)
  # ggsave(sprintf("%s/fig_%s_interaction_age.sdoh.pdf", fig_path,
  #                ifelse(i == 1, "3", "S4")),
  #        plot = plt,
  #        device = "pdf", units = "cm",
  #        width = 20.5+2, height = 24)
}
