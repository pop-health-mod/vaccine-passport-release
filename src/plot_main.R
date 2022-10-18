plot_vaccine_rate <- function(data, numer = "rate_1dose", denom = "n_unvax",
                              x_min = NULL, x_max = NULL, y_min = 0, y_max = 25000, 
                              colour_var = "age",
                              show_dates = 1:2, x_axis = TRUE, denom_factor = 5){
  # set default dates
  if(is.null(x_min)){
    x_min <- min(data$date_wk_end)
  }
  if(is.null(x_max)){
    x_max <- max(data$date_wk_end)
  }
  
  # compute rate
  data$rate_plt <- data[[numer]] / data[[denom]] * 10^denom_factor
  
  # plot
  if(is.null(colour_var)){
    p <- ggplot(data, aes_string(x = "date_wk_end", y = "rate_plt"))
  } else {
    p <- ggplot(data, aes_string(x = "date_wk_end", y = "rate_plt", col = colour_var))
  }
  p <- p +
    geom_vline(xintercept = POLICY[show_dates], linetype = "dashed") +
    
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = F) +
    
    labs(y = "Vaccination rate\n(vaccinations/100k unvaccinated people)") +
    
    theme(panel.grid.minor.x = element_blank())
  
  # display x axis and use month abbreviation (default is yes)
  if(x_axis){
    p <- p + 
      scale_x_date(breaks = date_breaks, date_labels = "%d %b") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p <- p +
      scale_x_date(breaks = date_breaks, date_labels = "%m-%d")
  }
  p
}

plot_vaccine_cover <- function(data, numer = "n_vacc_1dose", denom = "pop_rpdb",
                               x_min = NULL, x_max = NULL, y_min = 60, y_max = 100,
                               colour_var = "age",
                               show_dates = 1:2, x_axis = TRUE){
  # set default dates
  if(is.null(x_min)){
    x_min <- min(data$date_wk_end)
  }
  if(is.null(x_max)){
    x_max <- max(data$date_wk_end)
  }
  
  # compute coverage
  data$cv_dose <- data[[numer]] / data[[denom]] * 100
  
  # plot
  if(is.null(colour_var)){
    p <- ggplot(data, aes_string(x = "date_wk_end", y = "cv_dose"))
  } else {
    p <- ggplot(data, aes_string(x = "date_wk_end", y = "cv_dose", col = colour_var))
  }
  p <- p +
    geom_vline(xintercept = POLICY[show_dates], linetype = "dashed") +
    
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = F) +
    
    labs(y = "Vaccine coverage (1 dose, %)") +
    
    theme(panel.grid.minor.x = element_blank())
  
  # display x axis and use month abbreviation (default is yes)
  if(x_axis){
    p <- p + 
      scale_x_date(breaks = date_breaks, date_labels = "%d %b") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p <- p +
      scale_x_date(breaks = date_breaks, date_labels = "%m-%d")
  }
  p
}

# displays impact point estimates across age (y axis) and selected decile (x axis)
plot_heatmap_impact <- function(data, x_name, x_title,
                                y_name = "age", y_title = "Age",
                                plt_title,
                                impact_limits = c(-1, 6),
                                show_ci = FALSE, data_ci = NULL,
                                txt_size = 4){
  p <- ggplot(data, aes_string(x = x_name, y = y_name, fill = "impact")) +
    geom_tile(color = "black") +
    scale_fill_gradient2(limit = impact_limits, low = "#65156EFF", high = "red", mid = "white",
                         midpoint = 0) +
    
    labs(subtitle = plt_title, x = x_title, y = y_title) +
    coord_fixed() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      title = element_text(size = 12),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 11)
    )
  
  # display with confidence interval or just point estimate
  if(show_ci){
    # make table combining point estimate and confidence interval limits
    data_lab <- data %>% 
      select(age, all_of(x_name), impact) %>% 
      left_join(data_ci %>% select(age, all_of(x_name), impact_lci, impact_uci),
                by = c("age", x_name))
    
    # create label
    data_lab <- data_lab %>% 
      group_by(age, across(all_of(x_name))) %>% 
      mutate(impact_txt = sprintf("%s\n(%s\u2013%s)", 
                                  format(round(impact, 1), nsmall = 1), 
                                  format(round(impact_lci, 1), nsmall = 1), 
                                  format(round(impact_uci, 1), nsmall = 1))
             ) %>% 
      ungroup()
    
    p <- p +
      geom_text(data = data_lab, aes(label = impact_txt), color = "black", size = txt_size)
  } else {
    p <- p +
      geom_text(aes(label = format(round(impact, 1), nsmall = 1)), color = "black", size = txt_size)
  }
  return(p)
}

# takes in the plot for Figures 1/2 and
#     1) removes empty panels
#     2) moves x-axis next to top panel
format_ggtable_fig1.2 <- function(p){
  g <- ggplotGrob(p)
  
  # get the grobs that must be removed
  rm_grobs <- g$layout$name %in% c("panel-2-6", "panel-3-6")
  
  # remove grobs
  g$grobs[rm_grobs] <- NULL
  g$layout <- g$layout[!rm_grobs, ]
  
  # move axis closer to panel
  g$layout[g$layout$name == "axis-b-6", c("t", "b")] = c(9.5, 9.5)
  
  return(g)
}
