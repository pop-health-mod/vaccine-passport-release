# Alternative model specification functions ----
plot_sensitivity_cover <- function(data_sens, impact_txt_sens, sens_var,
                                   data_ci, data_obs, var_labels){
  # add sensitivity variable to bootstrap so that it only appears for main model
  data_ci[[sens_var]] <- factor("main model")
  
  # need to store model end separately to plot Quebec's labels properly
  if(PROVINCE == "qc"){
    model_end_plt <- MODEL_END_QC
  } else if(PROVINCE == "on"){
    model_end_plt <- MODEL_END_ON
  }
  
  # plot
  plot_vaccine_cover(data_sens,
                     x_min = min(data_sens$date_wk_end)-5, 
                     x_max = max(data_sens$date_wk_end)+5, 
                     y_min = 60, y_max = 100, 
                     colour_var = NULL,
                     x_axis = FALSE) +
    # confidence interval
    geom_ribbon(data = data_ci,
                aes(x = date_wk_end, ymin = cv_reg_lci, ymax = cv_reg_uci, fill = type),
                col = NA, alpha = alpha_plt, inherit.aes = F) +
    # fitted data
    geom_line(aes(col = type)) +
    # observed data
    geom_point(data = data_obs, aes(shape = "Observed\ndata"), alpha = alpha_plt) +
    # impact (need to plot old age UNDER)
    geom_label(data = impact_txt_sens %>% filter(age != "60_"), 
               aes(x = model_end_plt, y = vax_cover + 1, label = impact_lab),
               hjust = 1, vjust = 0, size = 4.3, alpha = 0.8,
               label.size = NA,
               inherit.aes = FALSE) +
    geom_label(data = impact_txt_sens %>% filter(age == "60_"), 
               aes(x = model_end_plt, y = vax_cover - 10, label = impact_lab),
               hjust = 1, vjust = 0, size = 4.3, alpha = 0.8,
               label.size = NA,
               inherit.aes = FALSE) +
    
    # variable level
    geom_label(data = var_labels, aes(x = min(data_sens$date_wk_end)-4, y = 99.7, 
                                      label = var_text),
               hjust = 0, vjust = 1, size = 5) +
    # policy announcements
    geom_label(data = policy_txt, aes(x = dates + 2, y = cv_pos, label = event),
               hjust = 0, vjust = 0.5, size = 5,
               label.size = NA, alpha = 0.8) +
    
    facet_grid(c(sens_var, "age")) +
    
    # aesthetics
    scale_y_continuous(breaks = c(6:10 * 10)) +
    labs(title = "", x = "Date (end of epidemiological week)",
         shape = "", col = "Model fit", fill = "Model fit") +
    colour_model + fill_model +
    
    # format text sizes, remove legend
    guides_man_fig1.2 +
    theme_man_fig1.2 +
    theme(legend.position = "none")
}

plot_sensitivity_rate <- function(data_sens, impact_txt_sens, sens_var,
                                  data_ci, data_obs, var_labels){
  # add sensitivity variable to bootstrap so that it only appears for main model
  data_ci[[sens_var]] <- factor("main model")
  
  # plot
  plot_vaccine_rate(data_sens,
                    x_min = min(data_sens$date_wk_end)-5, 
                    x_max = max(data_sens$date_wk_end)+5, 
                    y_max = 8400,
                    colour_var = NULL,
                    x_axis = FALSE) +
    # confidence interval
    geom_ribbon(data = data_ci, 
                aes(x = date_wk_end, ymin = rate_plt_lci, ymax = rate_plt_uci, fill = type),
                col = NA, alpha = alpha_plt, inherit.aes = F) +
    # fitted data
    geom_line(aes(col = type)) +
    # observed data
    geom_point(data = data_obs, aes(shape = "Observed\ndata"), alpha = alpha_plt) +
    
    # variable level
    geom_label(data = var_labels, aes(x = min(data_sens$date_wk_end)-4, y = 8350, 
                                      label = var_text),
               hjust = 0, vjust = 1, size = 5) +
    # policy announcements
    geom_label(data = policy_txt, aes(x = dates + 1, y = rate_pos, label = event),
               hjust = 0, vjust = 0.5, size = 5,
               label.size = NA, alpha = 0.8) +
    
    facet_grid(c(sens_var, "age")) +
    
    # aesthetics
    scale_y_continuous(breaks = c(0:4 * 2000)) +
    labs(title = "", x = "Date (end of epidemiological week)",
         shape = "", col = "Model fit", fill = "Model fit") +
    colour_model + fill_model +
    
    # format text sizes, remove legend
    guides_man_fig1.2 +
    theme_man_fig1.2 +
    theme(legend.position = "none")
}

#' y_pos should be 0.05 for coverage and 0.072 for rate
draw_sensitivity_labs <- function(lab_a, lab_b, lab_c,
                                  lab_d = NULL, lab_e = NULL, lab_f = NULL,
                                  y_pos){
  # Quebec labels
  grid.text(sprintf("A) %s", lab_a), y_pos, 0.994, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("B) %s", lab_b), y_pos, 0.847, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("C) %s", lab_c), y_pos, 0.70, just = "left",
            gp = gpar(fontsize = 18))
  # Ontario labels
  grid.text(sprintf("D) %s", ifelse(is.null(lab_d), lab_a, lab_d)), y_pos, 0.515, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("E) %s", ifelse(is.null(lab_e), lab_b, lab_e)), y_pos, 0.368, just = "left",
            gp = gpar(fontsize = 18))
  grid.text(sprintf("F) %s", ifelse(is.null(lab_f), lab_c, lab_f)), y_pos, 0.221, just = "left",
            gp = gpar(fontsize = 18))
}

# Hold baseline coverage constant (figures S3 and S5) ----
plot_ctfl_cov_points <- function(data_impact, data_impact_ci,
                                 var_name, strat_var_name,
                                 var_title, plt_title,
                                 impact_limits = c(-1, 3)){
  # extract CIs
  data_impact_ci <- data_impact_ci %>% 
    filter(strat_var == strat_var_name) %>% 
    mutate(vaxcov_lab = "Observed")
  
  # plot
  ggplot(data_impact, aes(x = vaxcov_lab, col = get(var_name))) +
    # null effect
    geom_hline(yintercept = 0, linetype = "dashed") +
    # data
    geom_point(aes(y = impact), position = position_dodge(0.5), size = 0.85) +
    geom_linerange(data = data_impact_ci, aes(ymin = impact_lci, ymax = impact_uci),
                   position = position_dodge(0.5), linewidth = 0.35) +
    # aesthetics
    coord_cartesian(ylim = impact_limits) +
    # labels
    labs(x = "Baseline vaccine coverage set to", y = "Vaccine passport impact\n(in percentage points)",
         col = var_title, title = plt_title)
}
