#' ---
#' title: "contrasts of time series Bavaria"
#' author: "Nadja Sauter"
#' date: "2020-10-24"
#' ---

# three  time series (daily counts till 01.05.2020)
  # a) Officialy reported cases
  # b) Disease onsets
  # c) Back Projections


### Data
# loading necessary data
load("./data/timeseries.RData")

## Description:
## a) Officialy reported cases: "dat"
  # Data on reported COVID-19 cases in Bavaria are official case reporting data collected
  # based on the German Infection Protection Act (IfSG) and provided
  # by the Bavarian Health and Food Safety authority.
## b) Disease onsets : "onset_bp"
  # data calculated in disease_onset_bav.R
## c) Back Projections: "backprojection"
  # data calculated in backproj_bav.R


### Libraries
library(tidyverse)
library(ggplot2)


### Tieme series
# transform data to right form for plotting
timeseries2 <- merge(onset_bp, backprojection, by = "date", all.x = TRUE, suffix = c(".onset", ".backprojection"))
timeseries3 <- merge(timeseries2, dat, by= "date", all.x= TRUE)
timeseries3_gather <- timeseries3 %>% gather(Timeseries,numbers,-c(date, t.onset, t.backprojection ))

# Plotting time series
timeseries3_gather %>%
  ggplot(aes(x=t.onset, y=numbers, colour=Timeseries))  +
  geom_col(data = timeseries3, aes(x=t.onset, y=n_reported),
           col = "lightgrey", fill = "lightgrey", alpha = 0.1, show.legend = TRUE) +
  geom_line(data = timeseries3_gather)+
  scale_color_manual(name = "Data:", values = c("#018571", "steelblue", "lightgrey"),  labels = c("Back projections", "Disease onsets", "Officially reported cases") )+
  theme_bw()+
  scale_x_continuous(breaks = seq(max(timeseries3$t.onset), min(timeseries3$t.onset), by = -7),
                     labels = strftime(seq(max(timeseries3$date), min(timeseries3$date), by = "-1 weeks"),
                                       format = "%d.%m.")) +
  scale_y_continuous(expand=expansion(mult = c(0, .1))) +
  labs(x = "Date", y = "number of cases")  +
  theme(
    axis.text=element_text(size = rel(1.5)),
    axis.title=element_text(size = rel(2.2)),
    legend.text = element_text(size = rel(1.5)),
    legend.title =element_text(size = rel(1.5)),
    legend.position = "bottom")
  ggsave("./results/timeseries_bavaria.png", width = 9, height = 5)
