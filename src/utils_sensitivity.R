# Alternative model specification functions ----
compute_impact_sens <- function(data_sensitivity, data_obs, sens_var){
  # drop main ITS from sensitivity (using fit from main model)
  data_sensitivity <- data_sensitivity[data_sensitivity[[sens_var]] != "main model", ]
  
  ## impact point estimates, comparing to observed observed
  # create a replicate of the observed data (only last timepoint) for each sensitivity fit
  indx_rep <- which(data_obs$date_wk_end == MODEL_END)
  data_obs_dupl <- data_obs[rep(indx_rep, 2), ] %>% 
    group_by(age) %>% 
    mutate(type = "a_observed", obs_rep = 2:(n()+1)) %>% 
    ungroup()
  
  data_obs_dupl[[sens_var]] <- factor(data_obs_dupl$obs_rep, 
                                      1:3,
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
