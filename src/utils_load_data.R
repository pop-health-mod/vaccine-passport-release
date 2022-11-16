#' load primary datasets
load_vaccination_data <- function(pr, dataset = NULL){
  if(pr == "qc"){
    # data from June 2021 to December 2022 with age up to 60+
    # either the whole dataset (else statement) or a sample of 10% of the data (~1,300 DAs; if statement)
    if(dataset == "passport_sample"){
      data <- fread("../vaccine-passport-data/vaccination_qc_passport_sample.csv")
    } else {
      data <- fread("../vaccine-passport-data/vaccination_qc_passport.csv")
    }
    
  } else if(pr == "on"){
    # data from June 2021 to December 2022 with age up to 60+
    data <- fread("../vaccine-passport-data/vaccination_on_passport.csv")
  } else {
    # throw error if province not available
    stop("Province must be 'qc' for Quebec or 'on' for Ontario")
  }
  
  # format
  data[, `:=`(date_wk_start = as.Date(date_wk_start),
              date_wk_end = as.Date(date_wk_end),
              codeDA = as.character(codeDA))]
  
  # remove province-wide "DA"
  data <- data[codeDA != "0" & codeDA != "99999999"]

  return(data)
}

# Load results (main analyses) ----
#' function to load 
#'        observed data
#'        ITS results
#'        bootstrap results
#'        counterfactual predicted vaccine coverage holding vaxcov constant
load_results <- function(directory, variable, result_type, data_obs, 
                         vaxcov_grp = NULL, R = 500){
  ### create file path
  model_nb <- case_when(variable == "age" ~ 1,
                        variable %in% c("income", "vismin") ~ 2,
                        variable %in% c("age.income", "age.vismin") ~ 3)
  model_letter <- case_when(variable == "age" ~ "",
                            grepl("income", variable) ~ "a",
                            grepl("vismin", variable) ~ "b")
  
  # add city suffix if necessary, just for observed and ITS
  city_suffix <- ifelse(DO_CMA, paste("_", CMA, sep = ""), "")
  
  if(result_type == "observed"){
    file_path <- sprintf("cover_%s%s_%s%s.csv", 
                         model_nb, model_letter, variable, city_suffix)
  } else if(result_type == "its"){
    file_path <- sprintf("predicted_model%s%s_%s_fitted_values%s.csv", 
                         model_nb, model_letter, variable, city_suffix)
  } else if(result_type == "bootstrap"){
    file_path <- sprintf("bootstrap_model%s%s_%s.R%s.csv", 
                         model_nb, model_letter, variable, R)
  } else if(result_type == "counterfactual_vaxcov"){
    file_path <- sprintf("predicted_%s%s_%s_set_vaxc_%s_fitted_values%s.csv",
                         model_nb, model_letter, substr(variable, 1, 3), 
                         vaxcov_grp, city_suffix)
  } else {
    stop("Provided 'result_type' is not valid, please verify.")
  }
  file_path <- paste(directory, file_path, sep = "/")
  
  ### read in results
  results_out <- read_csv(file_path, col_types = cols())
  
  ### fix variable types
  # model fits
  # rename if bootstrap, otherwise fix the string
  if(result_type == "bootstrap"){
    results_out$type <- factor(results_out$type, 
                               levels = c("Fitted", "Counterfactual"),
                               labels = c("Passport", "No passport\n(counterfactual)"))
  } else if(result_type %in% c("its", "counterfactual_vaxcov")){
    results_out$type <- factor(results_out$type, 
                               levels = c("Passport", "No passport\r\n(counterfactual)", "No passport\n(counterfactual)"),
                               labels = c("Passport", "No passport\n(counterfactual)", "No passport\n(counterfactual)"))
  }
  
  # quintile variables
  if(grepl("income", variable)){
    results_out$quin_income <- factor(results_out$quin_income)
  } else if(grepl("vismin", variable)){
    results_out$quin_vismin <- factor(results_out$quin_vismin)
  }

  return(results_out)
}

#' load all 3 single-variable models (age, income, vismin)
#' into a single data frame and to format it for plotting
load_merge_data <- function(directory, result_type, R = 500,
                            strat_vars = c("age", "quin_income", "quin_vismin")){
  ls_results <- vector("list", length = length(strat_vars))
  
  for(i in 1:length(strat_vars)){
    # load observed data
    cur_var <- strat_vars[i]
    ls_results[[i]] <- load_results(directory = directory, 
                                    gsub("quin_", "", cur_var), 
                                    result_type = result_type,
                                    R = R)
    
    # add variable name to data and rename
    ls_results[[i]]$strat_var <- cur_var
    if(cur_var == "age"){
      ls_results[[i]]$strat_lvl <- factor(
        ls_results[[i]][[cur_var]],
        levels = c("12_17", "18_29", "30_39", "40_49", "50_59", "60_"),
        labels = 1:6
      )
    } else {
      ls_results[[i]]$strat_lvl <- factor(ls_results[[i]][[cur_var]], 1:6)
    }
    names(ls_results[[i]])[names(ls_results[[i]]) == cur_var] <- "original_lvl"
    
    # compute rate for observed data
    if(result_type == "observed"){
      ls_results[[i]] <- ls_results[[i]] %>% 
        group_by(strat_lvl) %>% 
        mutate(rate_1dose = n_vacc_1dose - lag(n_vacc_1dose)) %>% 
        ungroup()
    }
  }
  
  # bind data, and reorder
  data_results <- bind_rows(ls_results)
  
  data_results <- data_results %>% 
    select(contains("replicate"), contains("type"), starts_with("date_wk"), 
           strat_var, strat_lvl, original_lvl,
           pop_rpdb, n_unvax, n_vacc_1dose, rate_1dose)
  
  # rename variables
  data_results$strat_var <- factor(data_results$strat_var,
                                   levels = c("age", "quin_income", "quin_vismin"),
                                   labels = c("Age", "Income quintile", "Visible minority quintile"))
  
  # rename type for plotting (so that colour legends match in size)
  if(result_type != "observed"){
    data_results$type <- factor(data_results$type, 
                                levels = c("Passport", "No passport\n(counterfactual)"),
                                labels = c("Passport\n", "No passport\n(counterfactual)"))
  }
  
  return(data_results)
}

join_bootstraps <- function(file_path, R_partition, R_total){
  # check number of partitions
  nb_part <- as.integer(R_total / R_partition)
  
  # list to store all results
  ls_boot <- vector("list", nb_part)
  
  # loop through all bootstrap files
  for(i in 1:nb_part){
    boot_results <- fread(sprintf(file_path, i * R_partition, R_total))
    
    # fix replicate number if necessary
    if(i > 1){
      boot_results$replicate_nb <- boot_results$replicate_nb + (R_partition * (i-1))
    }
    
    # store in list
    ls_boot[[i]] <- boot_results
  }
  
  # consolidate into a single data.table
  boot_results_compiled <- bind_rows(ls_boot)
  
  return(boot_results_compiled)
}

# Load results (supplementary analyses) ----
# original_levels must always have the main model first
load_data_sensitivity <- function(model_fit_path, data_its,
                                  original_levels, new_lvl_names,
                                  sens_var_original_name, sens_var_new_name){
  # load data
  data_sens <- read_csv(model_fit_path, col_types = cols())
  data_sens$type <- factor(data_sens$type, 
                           levels = c("Passport", "No passport\r\n(counterfactual)", "No passport\n(counterfactual)"),
                           labels = c("Passport", "No passport\n(counterfactual)", "No passport\n(counterfactual)"))
  
  data_sens <- data_sens[data_sens[[sens_var_original_name]] %in% original_levels, ]
  
  # reorder so that first model shown is main model
  data_sens[[sens_var_new_name]] <- factor(
    data_sens[[sens_var_original_name]],
    levels = as.character(original_levels),
    labels = new_lvl_names
  )
  
  # change labels to match nb of lines for both fits
  data_sens$type <- factor(data_sens$type, 
                           levels = c("Passport", "No passport\n(counterfactual)"),
                           labels = c("Passport\n", "No passport\n(counterfactual)"))
  
  # replace the chosen ITS model in sensitivity data with the main ITS data
  data_its[[sens_var_new_name]] <- factor(new_lvl_names[1])
  
  data_sens <- bind_rows(data_its,
                         data_sens[data_sens[[sens_var_new_name]] != new_lvl_names[1], ])
  
  # reorder variables to output
  data_sens <- data_sens %>% 
    select(all_of(c(sens_var_new_name, sens_var_original_name)), type, age,
           date_wk_end:rate_pred_cumul)
  
  return(data_sens)
}

#' function to load all counterfactual vaxcov results for a single variable
#' (age, income, vismin) into a single data frame and to format it for plotting
load_ctfl_cov_results <- function(directory, variable,
                                  R = 1000){
  vaxcov_grp_all <- c("49.9_under", "50_59.9_", "60_69.9_",
                      "70_79.9_", "80_89.9_", "90_over")
  
  ls_results <- vector("list", length = length(vaxcov_grp_all))
  
  for(i in 1:length(vaxcov_grp_all)){
    # load main ITS data
    ls_results[[i]] <- load_results(directory = directory, 
                                    variable = variable, result_type = "counterfactual_vaxcov",
                                    vaxcov_grp = vaxcov_grp_all[i],
                                    R = R)
    
    # add baseline variable to data
    ls_results[[i]]$vaxcov_grp_set_to <- vaxcov_grp_all[i]
  }
  
  # bind data, and reorder
  data_results <- bind_rows(ls_results)
  
  data_results <- data_results %>% 
    select(vaxcov_grp_set_to, type, starts_with("age"), starts_with("quin"), starts_with("date_wk"), 
           pop_rpdb, n_unvax, n_vacc_1dose, rate_1dose)
  
  return(data_results)
}

