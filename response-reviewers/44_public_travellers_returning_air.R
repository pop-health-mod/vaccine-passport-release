library(tidyverse)
library(lubridate)

theme_set(theme_bw())

# time when passports are rescinded
#' QC: March 14th 2022, https://www.msss.gouv.qc.ca/ministere/salle-de-presse/communique-3437/
#'     announced on Feb 15th
#' ON: March 1st 2022, https://news.ontario.ca/en/release/1001600/ontario-moving-to-next-phase-of-reopening-on-february-17
#'     announced on Feb 14th
df_pass_end <- data.frame(
  pr_name = factor(rep(c("Québec", "Ontario"), each = 4),
                   levels = c("Québec", "Ontario")),
  date = as.Date(c("2021-08-05", "2021-09-15", "2022-02-15", "2022-03-14",
                   "2021-09-01", "2021-09-22", "2022-02-14", "2022-03-01"))
)

# Format data ----
data_arrivals <- read_csv(
  "./response-reviewers/data-public-raw/statcan-travel-return-2410005301-eng.csv",
  skip = 9
)

# keep only columns for all travellers
# and first column
names(data_arrivals)
data_arrivals <- data_arrivals[, grepl("[A-Z]|^\\.\\.\\.1$", names(data_arrivals))]

# keep only rows with Canadian residents returning via air
data_arrivals <- data_arrivals[grep("Canadian residents.+air", data_arrivals$`...1`), ]

# transform to one row per month
data_arrivals <- data_arrivals %>% 
  select(-`...1`) %>% 
  pivot_longer(`January 2018`:`December 2022`, names_to = "date_m", values_to = "nb_travel")

# format number of travellers and date
data_arrivals <- data_arrivals %>% 
  mutate(nb_travel = as.integer(gsub(",", "", nb_travel)),
         date_m = as.Date(paste(date_m, "01"), "%B %Y %d"))

data_arrivals <- data_arrivals %>% 
  group_by(date_m) %>% 
  summarize(nb_travel = sum(nb_travel))

write.csv(data_arrivals, "./response-reviewers/data-public/air_travel_by_month.csv", row.names = F)

## make fake date to plot all years in single
data_arrivals <- data_arrivals %>% 
  mutate(
    yr = year(date_m),
    date_same = as.Date(sprintf("2005-%s-01", month(date_m)))
  )

# plot on same x axis
p_arrival <- ggplot(mapping = aes(x = date_same, y = nb_travel / 1e5, col = as.character(yr))) +
  geom_line(data = subset(data_arrivals, yr != 2021), linewidth = 0.5) +
  geom_line(data = subset(data_arrivals, yr == 2021), linewidth = 0.8) +
  
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_y_continuous(breaks = 0:5 * 5) +
  coord_cartesian(ylim = c(-1, max(data_arrivals$nb_travel) * 1.05 / 1e5), expand = F) +
  labs(x = "Date (month)", y = "Number of returning Canadian\ncitizens (in 100,000s)",
       col = NULL) +
  theme(
    # make axis and facet text legible
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6.5),
    
    # make legend key bigger and text legible
    # legend.key.width = unit(0.7, "line"),
    legend.text = element_text(size = 6.5),
    panel.grid.minor = element_blank()
  )

# plot
png("./response-reviewers/fig/fig_R4_air_arrivals.png",
    width = 12, height = 5, units = "cm", res = 320,
    type = "cairo-png")
p_arrival
dev.off()

# plot continuously
ggplot(data_arrivals, aes(x = date_m, y = nb_travel)) +
  # geom_point() +
  geom_line() +
  geom_vline(data = df_pass_end, aes(xintercept = date, col = pr_name)) +
  
  scale_x_date(date_labels = "%b %Y") +
  
  labs(x = "Date (month)", y = "Number of returning Canadian citizens",
       col = "Province")
