# Alternative model specification functions ----
compute_impact_sens <- function(data_sensitivity, data_obs, sens_var,
                                drop_main_model = TRUE){
  # drop main ITS from sensitivity (using fit from main model)
  if(drop_main_model){
    data_sensitivity <- data_sensitivity[data_sensitivity[[sens_var]] != "main model", ]
  }
  
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

# Hold baseline coverage constant (figures S3 and S5) ----
predict_ctfl_vaxcovgrp <- function(data_da, its_fit,
                                   var_select = c("week_anno", "vaxcov_grp"),
                                   strat_var){
  # predict rate and coverage setting vaxcov_passport to all possible values
  for(cur_vaxgrp in unique(data_da$vaxcov_grp)){
    # figure and output names
    if(strat_var == "age"){
      model_name <- sprintf("1_age_set_vaxc_%s", cur_vaxgrp)
    } else if(strat_var == "quin_income"){
      model_name <- sprintf("2a_inc_set_vaxc_%s", cur_vaxgrp)
    } else if(strat_var == "quin_vismin"){
      model_name <- sprintf("2b_vis_set_vaxc_%s", cur_vaxgrp)
    }
    
    # create dataset in which all DAs have same coverage group
    # Q: what would the results have been if we remove effect of baseline coverage
    #    from the passport's impact
    data_same_cv <- copy(data_da)
    data_same_cv <- data_same_cv[, vaxcov_passport := cur_vaxgrp]
    
    ## predict rate + coverage using new dataset
    prediction <- get_rate_and_coverage(data_same_cv, its_fit, var_select = var_select,
                                        strat_var = strat_var,
                                        predict_fit = TRUE)
    
    # save outputs as csv and plot rate and coverage
    if(strat_var == "age"){
      save_model_fit(prediction$model_fit[order(age, date_wk_end)][order(-type)], 
                     model_name = model_name)
    } else if(strat_var == "quin_income"){
      save_model_fit(prediction$model_fit[order(quin_income, date_wk_end)][order(-type)], 
                     model_name = model_name)
    } else if(strat_var == "quin_vismin"){
      save_model_fit(prediction$model_fit[order(quin_vismin, date_wk_end)][order(-type)], 
                     model_name = model_name)
    }
  }
}

join_impact_obs_ctfl <- function(impact_ctfl, impact_obs,
                                 strat_var_name, var_name){
  # keep only observed impact for target variable
  impact_obs <- impact_obs %>% 
    filter(strat_var == strat_var_name) %>% 
    mutate(vaxcov_grp_set_to = "Observed") %>% 
    select(type, vaxcov_grp_set_to, original_lvl, date_wk_end,
           impact, vax_cover, n_vacc_1dose, pop_rpdb)
  names(impact_obs)[names(impact_obs) == "original_lvl"] <- var_name
  
  # add observed impacts to counterfactual data table
  impact_ctfl <- bind_rows(impact_obs, impact_ctfl)
  
  # output
  impact_ctfl
}
