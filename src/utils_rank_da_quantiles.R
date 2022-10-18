# note: script assumes 'data' variable is of class data.table
# divide a variable into k quantiles (generalized)
#' data         dataset to be ranked
#' strat_var    name of variable on which to do ranking
#' rename_var   name to assign to new quantile variable (if not provided, will be quantiles)
#' k            number of quantiles to divide the data into
#' reverse      whether to rank by increasing (TRUE) or decreasing (FALSE)
#' by_cma.ca    whether to rank by CMA
#' balance_pop  whether to balance the population when ranking
#' denom        population count to use for balancing
get_quants_var <- function(data, strat_var, rename_var = NULL, k = 10, reverse = F,
                           by_cma.ca = FALSE, balance_pop = FALSE, denom = "pop_rpdb"){
  # get names of the other census identifiers if available in data, as well as denominators
  addtl_census_id <- names(data)[names(data) %in% c("CSDuid", "CSDname", "CTname", "SAC")]
  
  # restrict dataset to one row per DA
  if(!balance_pop){
    denom <- NULL
  }
  var_select <- c("codeDA", denom, addtl_census_id, strat_var)
  data_unique <- data[, ..var_select]
  data_unique <- unique(data_unique)
  
  # get CMA names, make into string vector to avoid issues with using them as list names
  if(by_cma.ca){
    data_unique[, SAC := paste("sac_", SAC, sep = "")]
    cma_list <- sort(unique(data_unique$SAC))
  }
  
  # balance population, i.e., ensure that each quantile has the same number of people in it
  # otherwise just take a single row per DA
  if(balance_pop){
    duplicate_indx <- rep(1:nrow(data_unique), data_unique[[denom]])
    data_strat_var <- data_unique[duplicate_indx, ]
  } else {
    data_strat_var <- data_unique
  }
  
  # get the quantile boundaries for the census variable, 
  #     divide into k groups and compute the (k-1) boundary points for the quantiles
  #     infinities added as end-pts of the vector for simplicity in the for-loop
  # if ranking by CMA instead of province-wide, need to get the quantile boundaries for each CMA
  if(by_cma.ca){
    ls_cma_quant_boundaries <- vector("list", length = length(cma_list))
    names(ls_cma_quant_boundaries) <- cma_list
    for(cur_sac in cma_list){
      ls_cma_quant_boundaries[[cur_sac]] <- quantile(data_strat_var[SAC == cur_sac][[strat_var]], probs = (1/k) * (1:(k-1)))
      ls_cma_quant_boundaries[[cur_sac]] <- c(-Inf, ls_cma_quant_boundaries[[cur_sac]], Inf)
    }
  } else {
    quant_boundaries <- quantile(data_strat_var[[strat_var]], probs = (1/k) * (1:(k-1)))
    quant_boundaries <- c(-Inf, quant_boundaries, Inf)
  }
  
  ### assign quantiles
  #   get the indices of all the DAs within the bounds of the i-th quantile, and assign
  #   i as the quantile for those DAs
  # NB: if there is not enough heterogeneity, it will not be possible to balance the groups
  data_unique$quantiles <- NA_integer_
  if(!by_cma.ca){
    for(i in 1:k){
      cur_indx <- (data_unique[[strat_var]] >  quant_boundaries[i] & 
                     data_unique[[strat_var]] <= quant_boundaries[i+1])
      data_unique$quantiles[cur_indx] <- i
    }
  } else {
    for(cur_sac in cma_list){
      for(i in 1:k){
        cur_indx <- data_unique$SAC == cur_sac
        cur_indx <- (cur_indx & 
                       data_unique[[strat_var]] >  ls_cma_quant_boundaries[[cur_sac]][i] & 
                       data_unique[[strat_var]] <= ls_cma_quant_boundaries[[cur_sac]][i+1])
        data_unique$quantiles[cur_indx] <- i
      }
    }
  }
  
  # reverse order of quantiles if specified
  if(reverse){
    data_unique$quantiles <- ((k+1) - data_unique$quantiles)
  }
  
  data_unique <- data_unique[, c("codeDA", "quantiles")]
  data <- left_join(data, data_unique, by = "codeDA")
  
  # rename if quantile name specified
  if(!is.null(rename_var)){
    setnames(data, "quantiles", rename_var)
  }
  
  return(data)
}
