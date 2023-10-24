library(tidyverse)

data_censuspop <- tibble()
for(cur_prov in c("QC", "ON")){
  ## load population (from https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=9810002001&pickMembers%5B0%5D=1.26&pickMembers%5B1%5D=2.1)
  census_pop <- read_csv(sprintf("./response-reviewers/data-public-raw/census_pop_prov.%s_2021.csv", cur_prov))
  
  census_pop <- census_pop %>% 
    select(age_old = `Age (in single years), average age and median age (128)`,
           gender = `Gender (3a)`, pop = VALUE)
  
  unique(census_pop$age_old)
  
  # convert age to number and keep only single years and totals (drop gender stratification)
  census_pop <- census_pop[census_pop$gender == "Total - Gender", ]
  census_pop <- census_pop[, c("age_old", "pop")]
  
  # compute size of each age group
  census_pop <- census_pop %>% 
    mutate(age = case_when(as.numeric(age_old) %in% 12:17 ~ "12_17", # 12 to 17 and 18 to 29 need to sum up individual age groups
                           as.numeric(age_old) %in% 18:19 ~ "18_29", 
                           age_old == "20 to 24 years" ~ "18_29", age_old == "25 to 29 years" ~ "18_29",
                           age_old == "30 to 34 years" ~ "30_39", age_old == "35 to 39 years" ~ "30_39",
                           # age groups 30+ can be created from the 5-year age groups
                           age_old == "30 to 34 years" ~ "30_39", age_old == "35 to 39 years" ~ "30_39",
                           age_old == "40 to 44 years" ~ "40_49", age_old == "45 to 49 years" ~ "40_49",
                           age_old == "50 to 54 years" ~ "50_59", age_old == "55 to 59 years" ~ "50_59",
                           age_old == "60 to 64 years" ~ "60_",
                           age_old == "65 years and over" ~ "60_",
                           age_old == "Total - Age" ~ "total_pop",
                           T ~ NA_character_))
  table(census_pop$age, census_pop$age_old)
  
  census_pop <- census_pop %>% group_by(age) %>% summarize(pop = sum(pop), .groups = "drop")
  
  census_pop <- census_pop %>% mutate(prov = cur_prov, .before = 1)
  
  data_censuspop <- bind_rows(data_censuspop, census_pop)
}

# rename province, remove extra, and save
data_censuspop <- data_censuspop[!is.na(data_censuspop$age), ]

write.csv(data_censuspop,
          "./response-reviewers/data-public/census_pop.csv", row.names = FALSE)
