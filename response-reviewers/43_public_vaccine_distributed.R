library(tidyverse)

theme_set(theme_bw())

# load data
dat_admin <- read_csv("./response-reviewers/data-public-raw/phac-vaccination-administration.csv")
dat_distr <- read_csv("./response-reviewers/data-public-raw/phac-vaccination-distribution.csv")

head(dat_admin)
head(dat_distr)

# keep only total numbers
dat_admin <- select(dat_admin, pruid, prename, date = as_of_date, nb_dose = numtotal_all_administered)
dat_distr <- select(dat_distr, pruid, prename, date = report_date, nb_dose = numtotal_all_distributed)

## fix typo, at row 4579 it says 2021 when it should be 2022
dat_admin[4570:4580, ]
dat_admin$date[dat_admin$prename == "Quebec" & dat_admin$nb_dose == 15136037] <- as.Date("2022-01-02")

# bind together and reshape
dat_admin <- mutate(dat_admin, type = "Administered vaccine doses")
dat_distr <- mutate(dat_distr, type = "Distributed vaccine doses")

data_vax <- bind_rows(dat_admin, dat_distr)

data_vax <- filter(data_vax, prename %in% c("Quebec", "Ontario"))
data_vax <- mutate(data_vax,
                   prename = factor(prename,
                                    c("Quebec", "Ontario"),
                                    c("Québec", "Ontario")),
                   type = factor(type,
                                 c("Distributed vaccine doses", "Administered vaccine doses")))
data_vax <- arrange(data_vax, prename, date)

# keep only until start of summer 2022 (or December 2021?)
data_vax <- filter(data_vax, date <= "2022-06-20")
# data_vax <- filter(data_vax, date <= "2021-12-31")

# plot
p <- ggplot(data_vax, aes(x = date, y = nb_dose / 10^6)) +
  # background for study period
  annotate(geom = "rect", alpha = .2,
           xmin = as.Date("2021-07-03"), xmax = as.Date("2021-11-13"),
           ymin = -10, ymax = max(data_vax$nb_dose) + 40) +
  annotate(geom = "text", label = "Study\nperiod", hjust = 0, vjust = 0.8,
           size = 2.5,
           x = as.Date("2021-07-03")+10,
           y = 36) +
  
  # plot data
  geom_line(aes(col = type)) +
  scale_colour_manual(values = c("black", "red")) + ## review
  
  facet_wrap(~prename) +
  scale_x_date(date_labels = "%m-%Y", date_breaks = "2 months") +
  coord_cartesian(ylim = c(0, max(data_vax$nb_dose)  / 10^6)) +
  
  labs(x = "Date", y = "Number of SARS-CoV-2 vaccine doses\n(millions)",
       col = NULL) +
       # caption = "Data from PHAC's Canadian report on COVID-19 vaccine doses and administered, respectively (as of 2023-02-22)") +
  theme(
    # make axis and facet text legible
    axis.title = element_text(size = 8),
    axis.text.y = element_text(size = 6.5),
    axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1),
    strip.text = element_text(size = 6.5, margin = margin(0.08, 0, 0.08, 0, "cm")),
    
    # make legend key bigger and text legible
    legend.text = element_text(size = 6.5),
    legend.position = "top",
    # remove space between plot and legend
    legend.margin = margin(0, 0, 0, 0, "cm"),
    legend.box.margin = margin(0, 0, 0, 0, "cm"),
    legend.box.spacing = unit(0, "cm"),
    # remove minor grid line
    panel.grid.minor.x = element_blank()
  )
p
# Public Health Agency of Canada. Canadian report on COVID-19 vaccine doses distributed. Ottawa: Public Health Agency of Canada; February 16, 2023. https://health-infobase.canada.ca/covid-19/vaccine-distribution/
# Public Health Agency of Canada. Canadian report on COVID-19 vaccine doses administered. Ottawa: Public Health Agency of Canada; February 22, 2023. https://health-infobase.canada.ca/covid-19/vaccine-administration/

png("./response-reviewers/fig/fig_R3_distributed_vaccines_cumul.png",
    width = 16, height = 7.5, units = "cm", res = 320,
    type = "cairo-png")
p
dev.off()

# check gap as of announcement and end of study period
data_vax %>% 
  filter(
    pruid == 24 & date %in% c(as.Date("2021-08-05"), as.Date("2021-09-18") + 5*7+1) | # no administered data on Oct 23 2021
    pruid == 35 & date %in% c(as.Date("2021-09-01"), as.Date("2021-10-09") + 5*7)
  ) %>% 
  unique() %>% 
  group_by(prename, date) %>% 
  mutate(dose_gap = nb_dose - lag(nb_dose))

