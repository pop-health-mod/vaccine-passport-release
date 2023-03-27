# Join bootstrap datasets together -----
library(tidyverse)
library(data.table)

source("./src/utils_load_data.R")

# create list of all directories to extract from
list_directory <- c("its-boot-part/qc", "its-boot-part/qc-mtl",
                    "its-boot-part/on", "its-boot-part/on-tor")
list_model_names <- c("bootstrap_model1_age",
                      "bootstrap_model2a_income", "bootstrap_model2b_vismin",
                      "bootstrap_model3a_age.income", "bootstrap_model3b_age.vismin")
list_paths <- c()

# create full paths
for(i in 1:length(list_directory)){
  for(j in 1:length(list_model_names)){
    # create base of the file path, add wildcards, then store
    str_path <- sprintf("../vaccine-passport-data/out/%s/%s.R", list_directory[i], list_model_names[j])
    str_path <- paste(str_path, "%s-%s.csv", sep = "")
    
    list_paths <- append(list_paths, str_path)
  }
}
rm(list_directory, list_model_names, i, j)

# load results
for(path_data in list_paths){
  print(path_data)
  # create path for final bootstrapped dataset
  new_path <- gsub("%s-%s", "10", path_data)
  new_path <- gsub("/its-boot-part/", "/its-boot-", new_path)
  
  # age and SDOH models are partitioned into 250 reps
  # interaction models are partitioned into 100 reps
  nb_part <- ifelse(grepl("model1_age|model2", path_data), 250, 100)
  
  # load partitioned data and write full dataset
  boot_current <- join_bootstraps(path_data, R_partition = 250, R_total = 1000)
  fwrite(boot_current, new_path)
}
