# Libraries & functions ----
library(fixest)
library(splines)
library(tidyverse)
library(data.table)
library(lubridate)

theme_set(theme_bw())

source("./src/utils_load_data.R")
source("./src/utils_its_bootstrap.R")

source("./02_setup_regression_formulas.R")

# done for QC only
PROVINCE <- "qc"
CMA <- ifelse(PROVINCE == "qc", "mtl", "tor")
DO_CMA <- F
which_dataset <- "passport"

# passport dates
source("./03_setup_policy_dates.R")

# file paths
path_out <- "../vaccine-passport-data/out/its-sensitivity-qc"

# Load vaccination coverage data ----
source("./04_setup_load_data.R", echo = TRUE)
rm(data_coverage)

# Set up variables for regressions ----
# specify knot positions (use only pre-intervention period, exclude first day of passport impact)
knot_pos <- quantile(unique(data_passport[date_wk_end < PASSPORT_START]$week), probs = c(.10, .50, .90))
knot_dates <- min(data_passport$date_wk_end) + knot_pos*7
print(knot_pos)
print(knot_dates)

# Load case data
# add DA code to the passport data
da_rss <- read_csv("../vaccine-passport-data/data-covid-incidence/da_rss_link.csv")
da_rss <- da_rss[, c("codeDA", "rss_code")]
da_rss$codeDA <- as.character(da_rss$codeDA)

# data_passport <- left_join(data_passport, da_rss, by = "codeDA")
data_passport <- merge(data_passport, da_rss, by = c("codeDA"), all.x = T, all.y = F)

# inspect
unique(data_passport$rss_code)
da_rss_uniq <- data_passport[date_wk_start == min(date_wk_start) & age == min(age),
                             .(codeDA, rss_code)]
table(da_rss_uniq$rss_code, useNA = "ifany")
# da_rss_uniq %>% filter(is.na(rss_code)) %>% View()
da_rss_uniq <- fill(da_rss_uniq, rss_code)

# merge with simple fill
data_passport <- data_passport[, rss_code := NULL]
data_passport <- merge(data_passport, da_rss_uniq, by = c("codeDA"), all.x = T, all.y = F)

## add cases
data_case <- read_csv("../vaccine-passport-data/data-covid-incidence/ts_by_rss.csv")
# data_case <- data_case %>% 
#   mutate(case_lag1 = lag(case), case_lag2 = lag(case, 2))

data_passport <- merge(data_passport,
                       # data_case[, c("rss", "date_wk_end", "case", "case_lag1", "case_lag2")],
                       data_case[, c("rss", "date_wk_end", "case")],
                       by.x = c("rss_code", "date_wk_end"),
                       by.y = c("rss", "date_wk_end"),
                       all.x = T, all.y = F)

# check for missing data
sum(is.na(data_passport$case))
# sum(is.na(data_passport$case_lag1))
# sum(is.na(data_passport$case_lag2))

# ---- Model 1: age only ----
## to avoid 0's in fit
data_passport$case <- data_passport$case + 0.0001

## fit model and predict rate and coverage
fmla_age <- paste(fmla_age, "+ log(case)")
fmla_age

fit_case <- femlm(as.formula(fmla_age),
                  data = data_passport,
                  family = "negbin",
                  offset =~log(n_unvax),
                  notes = F, mem.clean = T)

summary(fit_case)

pred_case <- get_rate_and_coverage(data_passport, fit_case, c("week_anno", "vaxcov_grp", "case"), "age")
save_model_fit(pred_case$model_fit[order(age, date_wk_end)][order(-type)], "qc_fit_w_cases")
