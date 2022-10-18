# Libraries ----
library(tidyverse)
library(data.table)

source("./src/utils_rank_da_quantiles.R")
PROVINCE <- "qc"
denom_used <- "pop_adjusted"

# Load DA and census data -----
## load DA identifiers
da_file <- ifelse(PROVINCE == "qc", 
                  "../vaccine-passport-data/list_DA_pop_qc.csv",
                  "../vaccine-passport-data/list_DA_pop_on.csv")

data_da <- fread(file = da_file)
data_da$codeDA <- as.character(data_da$codeDA)
data_da <- data_da[codeDA != "0" & codeDA != "99999999"]

# compute DA-level population size
data_da <- data_da[,
                   lapply(.SD, sum), 
                   by = .(codeDA),
                   .SDcols = c("pop_rpdb", "pop_adjusted", "n_vacc_end_ts")]

## load census data
census <- read_csv(sprintf("../vaccine-passport-data/Census_SDOH_DA_%s.csv", PROVINCE))
census <- census[, c("DAuid", "Population", "Visible_minority_Overall", "ATIPPE")]
census$codeDA <- as.character(census$DAuid)

# 100% of DAs linked to census
data_da <- left_join(data_da, census, by = "codeDA")
data_da <- data_da[!is.na(data_da$ATIPPE) & !is.na(data_da$Visible_minority_Overall)]

## add census identifiers
census_id <- read.csv("../vaccine-passport-data/DA_identifiers_13.17.22.csv", colClasses = c("CTname" = "character"))
census_id <- census_id %>% mutate(across(everything(), as.character))
data_da <- left_join(data_da, census_id[, c("DAuid", "CTname", "SAC")], by = c("codeDA" = "DAuid"))

# inspect NAs in Statistical Area Classification and Census Tract
sum(is.na(data_da$SAC))
sum(is.na(data_da$CTname))

# impute CTname and SAC as 999 for ON and QC
data_da[is.na(SAC), SAC := "999"]
if(PROVINCE == "qc"){
  data_da[is.na(CTname), CTname := "9924.00"]
} else if(PROVINCE == "on"){
  data_da[is.na(CTname), CTname := "9935.00"]
}

# Add quintiles to dataset ----
# set DAs with NA population to 0
if(PROVINCE == "on"){
  data_da[is.na(pop_rpdb), pop_rpdb := 0]
}

# income quintiles
data_da <- get_quants_var(data_da, "ATIPPE", "quin_income", k = 5, by_cma.ca = TRUE,
                          balance_pop = TRUE, denom = denom_used)
table(data_da$quin_income)

# visible minority quintiles
data_da <- get_quants_var(data_da, "Visible_minority_Overall", "quin_vismin", k = 5, reverse = TRUE,
                          balance_pop = TRUE, denom = denom_used)
table(data_da$quin_vismin)

# Save dataset (province-wide) ----
data_ranking <- data_da[, .(codeDA, SAC, CTname, pop_rpdb, pop_adjusted, ATIPPE, Visible_minority_Overall, 
                            quin_income, quin_vismin)]
data_ranking <- data_ranking[order(codeDA)]

data_ranking <- data_ranking[, Visible_minority_Overall := round(Visible_minority_Overall, 5)]
write.csv(data_ranking, sprintf("../vaccine-passport-data/da_ranked_%s_prov_balanced.csv", PROVINCE), row.names = FALSE)

# Re-rank by CMA ----
# remove ranked columns
data_da[, c("quin_income", "quin_vismin") := NULL]

# keep only CMA of interest
cma_code <- ifelse(PROVINCE == "qc", 462, 535)
data_da <- data_da[SAC == cma_code]

## ranking
# income quintiles
data_da <- get_quants_var(data_da, "ATIPPE", "quin_income", k = 5, balance_pop = TRUE, denom = denom_used)
table(data_da$quin_income)

# visible minority quintiles
data_da <- get_quants_var(data_da, "Visible_minority_Overall", "quin_vismin", k = 5, reverse = TRUE,
                          balance_pop = TRUE, denom = denom_used)
table(data_da$quin_vismin)

# Save dataset (CMA) ----
data_ranking <- data_da[, .(codeDA, SAC, CTname, pop_rpdb, pop_adjusted, ATIPPE, Visible_minority_Overall, 
                            quin_income, quin_vismin)]
data_ranking <- data_ranking[order(codeDA)]

cma_name <- ifelse(PROVINCE == "qc", "mtl", "tor")
data_ranking <- data_ranking[, Visible_minority_Overall := round(Visible_minority_Overall, 5)]
write.csv(data_ranking, sprintf("../vaccine-passport-data/da_ranked_%s_%s_balanced.csv", PROVINCE, cma_name), row.names = FALSE)
