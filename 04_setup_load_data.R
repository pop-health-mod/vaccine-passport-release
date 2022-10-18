# ---- Load vaccination coverage data ----
data_coverage <- load_vaccination_data(PROVINCE, which_dataset)

## adjust population
# assign the number of doses as of the end of the timeseries if doses > number of people
data_pop_adjust <- data_coverage[date_wk_end == MODEL_END, 
                                 lapply(.SD, sum), 
                                 by = .(codeDA, age, date_wk_start, date_wk_end),
                                 .SDcols = c("pop_rpdb", "n_vacc_1dose")]
data_pop_adjust[, vacc_diff := as.double(abs(pop_rpdb - n_vacc_1dose)), by = .(codeDA, age)]
data_pop_adjust <- data_pop_adjust[n_vacc_1dose > pop_rpdb]
data_pop_adjust[, 
                .(nb_da = .N, pop_ttl = sum(pop_rpdb), pop_mdn = median(as.double(pop_rpdb)),
                  vacc_diff_min = min(vacc_diff), vacc_diff_max = max(vacc_diff),
                  vacc_diff_median = median(vacc_diff), vacc_diff_.90 = quantile(vacc_diff, .90)), 
                by = .(age)][order(age)]

# use adjusted population as population variable 
data_coverage <- merge(data_coverage, data_pop_adjust[, .(codeDA, age, pop_adjusted = n_vacc_1dose)], 
                       by = c("codeDA", "age"), all = T)
data_coverage[!is.na(pop_adjusted), pop_rpdb := pop_adjusted]

# ---- Add vaccination rate to data (nb vaccinated) ----
## vaccination rate by age and DA
if(PROVINCE == "qc"){
  data_coverage[, `:=`(rate_1dose = n_vacc_1dose - lag(n_vacc_1dose)),
                by = .(codeDA, age)]
}
data_coverage[, `:=`(n_unvax = lag(pop_rpdb - n_vacc_1dose)),
              by = .(codeDA, age)]

data_coverage <- data_coverage %>% 
  select(codeDA, age, date_wk_start, date_wk_end, weekCDC, pop_rpdb, n_unvax, n_vacc_1dose, rate_1dose)

summary(data_coverage[date_wk_end < MODEL_END]$n_unvax)

## aggregate at DA level (SDOH analyses)
data_passport_da <- data_coverage[, 
                                  lapply(.SD, sum), 
                                  by = .(codeDA, date_wk_start, date_wk_end, weekCDC),
                                  .SDcols = c("pop_rpdb", "n_vacc_1dose")]
data_passport_da[, `:=`(rate_1dose = n_vacc_1dose - lag(n_vacc_1dose),
                        n_unvax = lag(pop_rpdb - n_vacc_1dose)),
                 by = .(codeDA)]
data_passport_da <- data_passport_da %>% 
  select(codeDA, date_wk_start, date_wk_end, weekCDC, pop_rpdb, n_unvax, n_vacc_1dose, rate_1dose)

summary(data_passport_da[date_wk_end < MODEL_END]$n_unvax)

# ---- Add census data ----
census <- fread(
  sprintf("../vaccine-passport-data/da_ranked_%s_%s_balanced.csv", 
          PROVINCE, ifelse(DO_CMA, CMA, "prov"))
)
census$codeDA <- as.character(census$codeDA)
census[, pop_rpdb := NULL]

# turn into factors
census[, quin_income := factor(quin_income)]
census[, quin_vismin := factor(quin_vismin)]

data_coverage <- merge(data_coverage, census, by = "codeDA", all.x = TRUE)
data_passport_da <- merge(data_passport_da, census, by = "codeDA", all.x = TRUE)

data_coverage <- data_coverage[!is.na(quin_vismin) & !is.na(quin_income)]
data_passport_da <- data_passport_da[!is.na(quin_vismin) & !is.na(quin_income)]

# ---* Restrict to CMA if running analysis at the CMA level ----
# Montreal    SAC == 462
# Toronto     SAC == 535
if(DO_CMA & CMA == "mtl"){
  data_coverage <- data_coverage[SAC == 462]
  data_passport_da <- data_passport_da[SAC == 462]
} else if(DO_CMA & CMA == "tor"){
  data_coverage <- data_coverage[SAC == 535]
  data_passport_da <- data_passport_da[SAC == 535]
}

# ---* Add CT unique identifier ----
data_coverage[, CTuid := paste(SAC, CTname, sep = "")]
data_passport_da[, CTuid := paste(SAC, CTname, sep = "")]

# ---* Check SDOH quintiles ----
# income deciles and quintiles
table(data_coverage[date_wk_end == min(date_wk_end) & age == min(age)]$quin_income)

# visible minority deciles and quintiles
table(data_coverage[date_wk_end == min(date_wk_end) & age == min(age)]$quin_vismin)

# ---- Subset data to date when trend became linear and add indicator ----
if(which_dataset %in% c("passport", "passport_sample")){
  data_passport <- data_coverage[date_wk_end >= "2021-07-03" & date_wk_end <= MODEL_END]
  unique(data_passport$date_wk_end)
  
  # create dummy date with the first data point as the 0-th date
  data_passport[, week := (weekCDC - min(weekCDC))]
  
  # indicator for during impact
  data_passport[, `:=`(pass_anno = case_when(date_wk_end <  PASSPORT_START ~ 0, 
                                             date_wk_end >= PASSPORT_START & date_wk_end <= PASSPORT_END ~ 1,
                                             date_wk_end >  PASSPORT_END ~ 0))]
  
  # centered time indicators
  data_passport[, week_anno := week - min(data_passport[pass_anno == 1]$week)]
  
  data_passport <- data_passport %>% 
    select(codeDA, SAC, CTname, CTuid, age, date_wk_start, date_wk_end, 
           week, week_anno, pass_anno,
           pop_rpdb:rate_1dose, ATIPPE, Visible_minority_Overall,
           ends_with("_income"), ends_with("_vismin"))
  
  # inspect indicators
  data_passport[codeDA == min(codeDA) & age == min(age),
                .(date_wk_start, date_wk_end, week, week_anno, pass_anno)]
  
  ## add time and dummy indicators to aggregate dataset
  data_passport_da <- data_passport_da[date_wk_end >= "2021-07-03" & date_wk_end <= MODEL_END]
  
  data_passport_da <- left_join(data_passport_da,
                                data_passport[codeDA == min(codeDA) & age == "12_17", 
                                              .(date_wk_end, week, week_anno, pass_anno)],
                                by = c("date_wk_end"))
  
  data_passport_da <- data_passport_da %>% 
    select(codeDA, SAC, CTname, CTuid, date_wk_start, date_wk_end, 
           week, week_anno, pass_anno, 
           pop_rpdb:rate_1dose, ATIPPE, Visible_minority_Overall,
           ends_with("_income"), ends_with("_vismin"))
  
  # add vaccine coverage on time-series start date
  # all rows with population = 0 have NA as vaccination coverage, impute 90_over
  data_passport <- get_vaxcoverage(data_passport, time_coverage = min(data_passport$date_wk_end))
  data_passport[, vaxcov_grp := case_when(vaxcov_at_t >= .90 ~ "90_over",
                                          vaxcov_at_t >= .80 ~ "80_89.9_", vaxcov_at_t >= .70 ~ "70_79.9_",
                                          vaxcov_at_t >= .60 ~ "60_69.9_", vaxcov_at_t >= .50 ~ "50_59.9_",
                                          vaxcov_at_t <  .50 ~ "49.9_under",
                                          is.na(vaxcov_at_t) ~ "90_over")]
  
  data_passport_da <- get_vaxcoverage(data_passport_da, time_coverage = min(data_passport_da$date_wk_end), strata = NULL)
  data_passport_da[, vaxcov_grp := case_when(vaxcov_at_t >= .90 ~ "90_over",
                                             vaxcov_at_t >= .80 ~ "80_89.9_", vaxcov_at_t >= .70 ~ "70_79.9_",
                                             vaxcov_at_t >= .60 ~ "60_69.9_", vaxcov_at_t >= .50 ~ "50_59.9_",
                                             vaxcov_at_t <  .50 ~ "49.9_under",
                                             is.na(vaxcov_at_t) ~ "90_over")]
}
