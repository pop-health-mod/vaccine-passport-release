# dates for announcement and implementation of vaccine passport
#       and for start and end of impact period
POLICY <- case_when(PROVINCE == "qc" ~ as.Date(c("2021-08-05", "2021-09-15")),
                    PROVINCE == "on" ~ as.Date(c("2021-09-01", "2021-09-22")))
# first day passport assumed to have an impact
PASSPORT_START <- case_when(PROVINCE == "qc" ~ as.Date("2021-08-14"),
                            PROVINCE == "on" ~ as.Date("2021-09-04"))
# last day passport assumed to have an impact
PASSPORT_END <- case_when(PROVINCE == "qc" ~ as.Date("2021-09-18"),
                          PROVINCE == "on" ~ as.Date("2021-10-09"))
MODEL_END <- as.Date(PASSPORT_END + 5*7)
