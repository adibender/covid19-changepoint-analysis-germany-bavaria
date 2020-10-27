#' ---
#' title: "contrasts of time series Germany"
#' author: "Nadja Sauter"
#' date: "2020-10-24"
#' ---

# three  timee series (daily number till 01.05.2020)
# a) Officialy reported cases
# b) Disease onsets
# c) Back Projections


### Data
# loading necessary data
load("./data/timeseries.RData")

## Description:
## a) Officialy reported cases: "meldedaten_bereinigt"
# Data on reported COVID-19 cases from Project CoronaMaps (https://corona.stat.uni-muenchen.de/maps/),
# with data from RKI
## b) Disease onsets : "onset"
# data calculated in disease_onset_ger.R
## c) Back Projections: "backprojection"
# data calculated in backproj_ger.R


### Libraries
library(tidyverse)
library(ggplot2)


### Tieme series
# transform data to right form for plotting
onset = onset %>% mutate(t=1:n())
timeseries2 <- merge(onset, backprojection, by = "date", all.x = TRUE, suffix = c(".onset", ".backprojection"))
timeseries3 <- merge(timeseries2, meldedaten_bereinigt, by = "date", all.x = TRUE)
gather_timesereis3 <- gather(timeseries3, Timeseries, numbers, c(-date, -t))


# Plotting time series
gather_timesereis3 %>%
  ggplot(aes(x=t, y=numbers, colour=Timeseries))  +
  geom_col(data = timeseries3, aes(x=t, y=n_reported),
           col = "lightgrey", fill = "lightgrey", alpha = 0.1, show.legend = TRUE) +
  geom_line()+
  scale_color_manual(name = "Data:", values = c("#018571", "steelblue", "lightgrey"),  labels = c("Back projections", "Disease onsets", "Officially reported cases") )+
  theme_bw()+
  scale_x_continuous(breaks = seq(max(timeseries3$t), min(timeseries3$t), by = -7),
                     labels = strftime(seq(max(timeseries3$date), min(timeseries3$date), by = "-1 weeks"),
                                       format = "%d.%m.")) +
  scale_y_continuous(expand=expansion(mult = c(0, .1))) +
  labs(x = "Date", y = "number of cases")   +
  theme(
    axis.text=element_text(size = rel(1.5)),
    axis.title=element_text(size = rel(2.2)),
    legend.text = element_text(size = rel(1.5)),
    legend.title =element_text(size = rel(1.5)),
    legend.position = "bottom")
  ggsave("./graphics/germany_timeseries.png", width = 9, height = 5)
