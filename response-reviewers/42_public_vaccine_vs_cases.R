library(tidyverse)
library(gridExtra)
library(lubridate)

theme_set(theme_bw())

data_vax <- read_csv("./response-reviewers/data-public/vax_by_age.csv")
data_vax <- rename(data_vax, rate_vax_wkly = rate_wkly)

# dates for passport
passport_mtx <- matrix(c("2021-08-05", "2021-09-15",
                         "2021-09-01", "2021-09-22"),
                       dimnames = list(c("QC", "ON"), c("start", "end")),
                       byrow = T, nrow = 2)

# Load and process case data ----
## quebec
data_case_qc <- read_csv("./response-reviewers/data-public-raw/covid-case-qc-graph_1-1_page_par_region.csv")
names(data_case_qc) <- c("date", "case_new", "avg_wkly", "notes")
data_case_qc$date <- as.Date(data_case_qc$date)
data_case_qc$prov <- "QC"
data_case_qc$notes <- NULL

## ontario
data_case_on <- read_csv("./response-reviewers/data-public-raw/covid-case-on-All case trends data.csv")
data_case_on <- data_case_on[, c(1, 2, 6, 7)]
names(data_case_on) <- c("date", "phu", "case_new", "avg_wkly")
data_case_on$date <- as.Date(data_case_on$date, "%B %d, %Y")
data_case_on$prov <- "ON"
data_case_on$avg_wkly <- as.double(data_case_on$avg_wkly)

data_case <- bind_rows(data_case_qc, data_case_on)
rm(data_case_qc, data_case_on)

# compute weekly totals
data_case <- data_case %>% 
  mutate(epiyr = epiyear(date), epiwk = epiweek(date)) %>%
  
  group_by(prov, epiyr, epiwk) %>%
  
  summarize(
    date_wk_start = min(date), date_wk_end = max(date),
    rate_case_wkly = sum(case_new),
    .groups = "drop"
  )

## Add population counts ----
census_pop <- read_csv("./response-reviewers/data-public/census_pop.csv")
census_pop$age[census_pop$age == "60_"] <- "60+"

data_vax <- left_join(data_vax,
                      subset(census_pop, age != "total_pop"),
                      by = c("prov", "age"))

data_case <- left_join(data_case,
                       subset(census_pop, age == "total_pop")[, c("prov", "pop")],
                       by = c("prov"))

rm(census_pop)

# group all 12+ vaccinations
data_vax_prov <- data_vax %>% 
  group_by(prov, epiyr, epiwk, date_wk_start, date_wk_end) %>% 
  summarize(across(rate_vax_wkly:pop, sum), .groups = "drop")

# Plot cases along vaccinations ----
plot_case_and_vax <- function(data_cases, data_vaccinations){
  ggplot(mapping = aes(x = date_wk_end)) +
    geom_line(data = data_cases, aes(y = rate_case_wkly, col = "Cases")) +
    geom_line(data = data_vaccinations, aes(y = rate_vax_wkly, col = "Vaccina-\ntions")) +
    
    scale_x_date(date_labels = "%m\n%Y") +
    scale_colour_manual(values = c("red", "black")) +
    
    facet_wrap(~pr_name, nrow = 1) +
    
    theme(
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
    ) +
    labs(x = "Date (month)", col = NULL)
}

## Plot ----
#' overall weak association between cases and vaccination
#'      During vaccine passport "impact period", increase in vaccinations
#'      occurs in QC *before* increase in cases this likely signals that it
#'      is the passport, not case increases, that is leading to people
#'      getting vaccinated. In Ontario they occur concurrently,
#'      meaning that it is unlikely it was the case loads that lead to increases.

# Compute per population rates
data_case <- data_case %>% 
  mutate(rate_raw = rate_case_wkly,
         rate_case_wkly = rate_raw / pop * 1e5)
data_vax <- data_vax %>% 
  mutate(rate_raw = rate_vax_wkly,
         rate_vax_wkly = rate_raw / pop * 1e5)
data_vax_prov <- data_vax_prov %>% 
  mutate(rate_raw = rate_vax_wkly,
         rate_vax_wkly = rate_raw / pop * 1e5)

data_case <- data_case %>% mutate(
  pr_name = factor(prov, c("QC", "ON"), c("Québec", "Ontario"))
)
data_vax_prov <- data_vax_prov %>% mutate(
  pr_name = factor(prov, c("QC", "ON"), c("Québec", "Ontario"))
)

# start up to end of summer 2022
p1 <- plot_case_and_vax(
  mutate(data_case, rate_case_wkly = rate_case_wkly * 10),
  data_vax_prov) +
  coord_cartesian(xlim = as.Date(c("2020-03-01", "2022-06-15")),
                  ylim = c(0, 14000)) +
  scale_y_continuous(
    # name = "Weekly vaccination rate (1st-dose vaccinations/100k people)",
    name = "",
    breaks = c(0, 4, 8, 12) * 1000,
    sec.axis = sec_axis(trans =~. * 10^-1,
                        # name = "Weekly confirmed COVID-19 cases per 100,000 people",
                        name = "",
                        breaks = c(0, 4, 8, 12) * 100)
  ) +
  labs(title = "A) Vaccination & cases from March 2020 to July 2022")

## our study period
p2 <- plot_case_and_vax(
  mutate(data_case, rate_case_wkly = rate_case_wkly * 10),
  data_vax_prov) +
  # passport dates
  geom_vline(
    data = data.frame(
      pr_name = factor(rep(c("Québec", "Ontario"), each = 2),
                       levels = c("Québec", "Ontario")),
      date = as.Date(c("2021-08-05", "2021-09-15",
                       "2021-09-01", "2021-09-22"))
    ),
    aes(xintercept = date), linetype = "dashed"
  ) +
  # aesthetics
  coord_cartesian(xlim = as.Date(c("2021-06-01", "2021-11-13")),
                  ylim = c(0, 3000)) +
  scale_y_continuous(
    name = "Weekly vaccination rate (1st-dose vaccinations/100k people)",
    breaks = 0:3 * 1000,
    sec.axis = sec_axis(trans =~. * 10^-1, name = "Weekly confirmed COVID-19 cases per 100,000 people",
                        breaks = 0:3 * 100)
  ) +
  labs(title = "B) Vaccination & cases from June to November 2021 (study period)")

## December 2021 to April 2022
p3 <- plot_case_and_vax(
  mutate(data_case, rate_case_wkly = rate_case_wkly),
  data_vax_prov) +
  # aesthetics
  coord_cartesian(xlim = as.Date(c("2021-12-01", "2022-04-30")),
                  ylim = c(0, 1000)) +
  scale_y_continuous(
    # name = "Weekly vaccination rate (1st-dose vaccinations/100k people)",
    name = "",
    breaks = seq(0, 1000, 250),
    sec.axis = sec_axis(trans =~.,
                        name = ""
                        # name = "Weekly confirmed COVID-19 cases per 100,000 people"
                        )
  ) +
  labs(title = "C) Vaccination from December 2021 to April 2022")

## Save plot ----
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
png("./response-reviewers/fig/fig_R1_weekly_vax_and_cases.png",
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

