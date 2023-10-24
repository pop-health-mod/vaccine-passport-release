library(tidyverse)
library(gridExtra)
library(lubridate)

theme_set(theme_bw())

# Process rate data ----
### load Quebec
data_qc <- read_csv("./response-reviewers/data-public-raw/vax-age-qc-COVID19_Qc_Vaccination_CatAge.csv")

# keep only daily rates for 1st and 2nd doses
data_qc <- data_qc[, grepl("date|Numero1_jour$|Numero2_jour$", names(data_qc))]
names(data_qc)

# split column names into (1) age group, (2) number of doses and (3) daily rate
data_qc <- data_qc %>% 
  pivot_longer(cols = Age_0_4_ans_DOSE_Numero1_jour:Age_85_110_ans_DOSE_Numero2_jour, names_prefix = "Age_",
               names_to = c("age_old", "dose"), values_to = "rate_daily",
               names_sep = "_ans_DOSE_Numero") %>% 
  mutate(dose = gsub("_jour", "", dose))

### load Ontario (NOTE: all data supplied is in cumulative)
data_on <- read_csv("./response-reviewers/data-public-raw/vax-age-on-All Covid-19 vaccine trends data.csv")

# remove columns with percent fully vaccinated
data_on <- data_on[, c("Geographic area", "Date", "Epi week",
                       grep("At least 1 dose count", names(data_on), value = T))]

# keep only province-wide data
data_on <- subset(data_on, `Geographic area` == "Ontario")

# create date variable
data_on$date <- as.Date(data_on$Date, "%B %d, %Y")

names(data_on)
names(data_on) <- case_when(
  names(data_on) == "At least 1 dose count" ~ "ttl_pop",
  grepl("At least 1 dose", names(data_on)) ~ gsub("At least 1 dose count\\: | ", "", names(data_on)),
  T ~ names(data_on)
)

# turn into one row per age group
data_on <- data_on %>% 
  pivot_longer(cols = ttl_pop:`80+`, names_to = "age_old", values_to = "vax_cumul")

### Ignore
# keep only 1st and 2nd dose (99 is flag for 'fully vax', i.e. 2 doses)
# after 2021-12-01, two-dose counting is moved from "2 doses" to "fully vaccinated"
###

### ensure age names match and join datasets
unique(data_qc$age_old)
unique(data_on$age_old)

# QC's age groups all fit within ON's groups, so regroup for QC and edit ON's names
data_qc <- data_qc %>% 
  mutate(age = case_when(age_old == "12_17" ~ "12_17",
                         age_old %in% c("18_24", "25_29") ~ "18_29",
                         age_old %in% c("30_34", "35_39") ~ "30_39",
                         age_old %in% c("40_44", "45_49") ~ "40_49",
                         age_old %in% c("50_54", "55_59") ~ "50_59",
                         age_old %in% c("60_64", "65_69") ~ "60_69",
                         age_old %in% c("70_74", "75_79") ~ "70_79",
                         age_old %in% c("80_84", "85_110") ~ "80+",
                         TRUE ~ NA_character_))

data_qc <- data_qc[!is.na(data_qc$age), ]

data_qc <- data_qc %>% 
  group_by(date, age, dose) %>% 
  summarize(rate_daily = sum(rate_daily), .groups = "drop")

# process age groups for Ontario
data_on$age <- gsub("-", "_", data_on$age_old)
data_on$age <- gsub("yrs", "", data_on$age)

data_on <- data_on[!(data_on$age %in% c("ttl_pop", "under5", "5_11")), ]

# check names match  
unique(data_qc$age) %>% sort()
unique(data_on$age) %>% sort()

# compute rate
data_on <- data_on %>% 
  group_by(age) %>% 
  mutate(rate_daily = vax_cumul - lag(vax_cumul)) %>% 
  ungroup()

data_on <- data_on[, c("date", "age", "rate_daily")]

# create joint dataset of vaccinations and reduce number of age groups
data_vax <- bind_rows(
  mutate(data_qc, prov = "QC"),
  mutate(data_on, prov = "ON")
) %>% 
  select(prov, date:rate_daily)

# keep only first dose
data_vax <- data_vax %>% filter(dose == 1 | prov == "ON") %>% select(-dose)

# group ages over 60 and group by week
data_vax <- data_vax %>% 
  mutate(age = ifelse(age %in% c("60_69", "70_79", "80+"), "60+", age))

data_vax <- data_vax %>% 
  mutate(epiyr = epiyear(date), epiwk = epiweek(date)) %>%
  
  group_by(prov, age, epiyr, epiwk) %>%
  
  summarize(
    date_wk_start = min(date), date_wk_end = max(date),
    rate_wkly = sum(rate_daily),
    .groups = "drop"
  )

## Save processed data ----
write.csv(data_vax, "./response-reviewers/data-public/vax_by_age.csv", row.names = FALSE)

## Add census population counts ----
census_pop <- read_csv("./response-reviewers/data-public/census_pop.csv")
census_pop <- subset(census_pop, age != "total_pop")

data_vax <- left_join(data_vax, census_pop, by = c("prov", "age"))
rm(census_pop)

# compute unvaccinated
data_vax <- data_vax %>% arrange(prov, age, date_wk_end)

data_vax <- data_vax %>% 
  group_by(prov, age) %>% 
  mutate(rate_wkly = ifelse(is.na(rate_wkly), 0, rate_wkly),
         vacc_ttl = cumsum(rate_wkly),
         pop_unvax = ifelse(vacc_ttl == 0, pop, pop - lag(vacc_ttl))) %>% 
  ungroup()

# Plot vaccinations ----
data_vax <- data_vax %>% mutate(pr_name = factor(prov,
                                                 c("QC", "ON"),
                                                 c("Québec", "Ontario"))
                                )

data_vax <- filter(data_vax, date_wk_end <= "2022-08-01")

# format age labels
data_vax <- data_vax %>% mutate(age = gsub("_", "\u2013", age))

## function
plot_vax <- function(data, rate_var = "rate_wkly", study_period = FALSE){
  p <- ggplot(data, aes(x = date_wk_end, y = get(rate_var))) 
  
  if(study_period){
    p <- p +
      # background for study period
      annotate(geom = "rect", alpha = .2,
               xmin = as.Date("2021-07-03"), xmax = as.Date("2021-11-13"),
               ymin = -5e4, ymax = max(data[[rate_var]], na.rm = T) + 5e4) +
      annotate(geom = "text", label = "Study\nperiod", hjust = 0, vjust = 0.6,
               size = 2.5,
               x = as.Date("2021-07-03")+10,
               y = max(data[[rate_var]] * .9, na.rm = T))
  }
  
  p <- p +
    # vaccination data
    geom_line(aes(col = age), linewidth = 0.3) +
    
    facet_wrap(~pr_name) +
    scale_x_date(date_labels = "%m\n%Y") +
    scale_colour_viridis_d(end = 0.9, direction = -1) +
    
    labs(x = "Date (month)", y = "First doses administered (weekly)",
         col = "Age group")
  
  p + theme(
    # resize title
    plot.title = element_text(size = 7),
    
    # make axis and facet text legible
    axis.title = element_text(size = 8),
    axis.text.y = element_text(size = 6.5),
    # axis.text.x = element_text(size = 6.5, angle = 25, hjust = 1),
    axis.text.x = element_text(size = 6.5),
    strip.text = element_text(size = 6.5, margin = margin(0.08, 0, 0.08, 0, "cm")),
    
    # make legend key bigger and text legible
    # legend.key.width = unit(0.7, "line"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6.5),
    legend.position = "none",
    
    plot.background = element_blank()
  )
}

## entire plot
data_vax <- mutate(data_vax, rate_unvax = rate_wkly / pop_unvax * 1e5)
data_vax <- mutate(data_vax, rate_pop = rate_wkly / pop * 1e5)

# start up to end of summer 2022
p1 <- plot_vax(data_vax, "rate_pop", study_period = TRUE) +
  coord_cartesian(xlim = c(min(data_vax$date_wk_end), as.Date("2022-06-15")),
                  ylim = c(0, max(data_vax$rate_pop, na.rm = T))) +
  scale_x_date(date_breaks = "2 months", date_labels = "%m\n%Y") +
  labs(title = "A) Vaccination from January 2021 to July 2022",
       y = "Weekly vaccination rate (1st-dose vaccinations/100k people)")

# our study period
p2 <- plot_vax(data_vax, "rate_pop") +
  coord_cartesian(xlim = as.Date(c("2021-06-01", "2021-11-13")), ylim = c(0, 3000)) +
  geom_vline(
    data = data.frame(
      pr_name = factor(rep(c("Québec", "Ontario"), each = 2),
                       levels = c("Québec", "Ontario")),
      date = as.Date(c("2021-08-05", "2021-09-15",
                       "2021-09-01", "2021-09-22"))
    ),
    aes(xintercept = date), linetype = "dashed"
  ) +
  labs(title = "B) Vaccination from June to November 2021 (study period)",
       y = "Weekly vaccination rate (1st-dose vaccinations/100k people)")

# time when passports are rescinded
#' QC: March 14th 2022, https://www.msss.gouv.qc.ca/ministere/salle-de-presse/communique-3437/
#'     announced on Feb 15th
#' ON: March 1st 2022, https://news.ontario.ca/en/release/1001600/ontario-moving-to-next-phase-of-reopening-on-february-17
#'     announced on Feb 14th
df_pass_end <- data.frame(
  pr_name = factor(c("Québec", "Québec", "Ontario", "Ontario"),
                   levels = c("Québec", "Ontario")),
  date_wk_end = as.Date(c("2022-02-15", "2022-03-14",
                          "2022-02-14", "2022-03-01"))
)
 

# ggplot(subset(data_vax, date_wk_end >= "2022-01-01" & date_wk_end <= "2022-04-30")) +
p3 <- plot_vax(data_vax, "rate_pop") +
  geom_vline(data = df_pass_end, aes(xintercept = date_wk_end), linetype = "dashed") +
  
  coord_cartesian(xlim = as.Date(c("2021-12-01", "2022-04-30")), ylim = c(0, 1500)) +
  labs(title = "C) Vaccination from December 2021 to April 2022",
       y = "Weekly vaccination rate (1st-dose vaccinations/100k people)")

# extract colour label
tmp <- ggplot_gtable(
  ggplot_build(
    p1 + theme(legend.position = "right")
  )
)
legend_pos <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
rate_legend <- tmp$grobs[[legend_pos]]

# extract x axis
x_axis <- cowplot::get_plot_component(p1, "xlab-b")
p1 <- p1 + theme(axis.title.x = element_blank()) + labs(y = "")
p2 <- p2 + theme(axis.title.x = element_blank())
p3 <- p3 + theme(axis.title.x = element_blank()) + labs(y = "")

# plot
# grid.arrange(p1, p2, p3, ncol = 1)

png("./response-reviewers/fig/fig_R2_vax_rate_by_period.png",
    width = 16, height = 12.5, units = "cm", res = 320,
    type = "cairo-png")
grid.arrange(p1, p2, p3, x_axis, rate_legend,
             ncol = 2,
             layout_matrix = matrix(c(1, 5,
                                      2, 5,
                                      3, 5,
                                      4, 5),
                                    ncol = 2, byrow = T),
             heights = unit(c(4, 4, 4, 0.5), "cm"),
             widths = unit(c(14, 2), "cm"))
dev.off()
