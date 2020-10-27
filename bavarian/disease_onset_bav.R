#' ---
#' title: "Segmented regression of daily COVID-19 infections using back-projected epidemic curve"
#' author: "Michael Höhle, Felix Günther, Helmut Küchenhoff, Andreas Bender"
#' date: "2020-10-23" new version with data unit May,1st
#' ---

#' # Segmented regression of daily COVID-19 infections using back-projected epidemic curve
#'
#' Authors: Michael Höhle (<hoehle@math.su.se>), Felix Günther, Helmut Küchenhoff, Andreas Bender<p>
#' ---

#+ load packages, message=FALSE, warning=FALSE
# Load packages
library(tidyverse)
library(segmented)
library(surveillance)
library(future)
library(future.apply)

#' ## Introduction
#'
#' This document contains a reproducible analysis to infer change-points for the back-projection
#' of the curve of daily disease onsets of reported German COVID-19 cases.
#'
#' ### RKI data
#'
#' For Bavaria we have access to the raw data and provide a thorough description and sensitivity
#' analysis: https://www.stablab.stat.uni-muenchen.de/_assets/docs/nowcasting_covid19_bavaria.pdf
#'
#' About 26% of the total cases in Germany originate from Bavaria, but in order to compare
#' results with Dehning et al. (2020) we use the processed RKI data covering entire Germany.
#'
#' The RKI provides a processed curve of daily new disease onsets of reported German COVID-19 cases -
#' the procedure is described in detail in an der Heiden and Hamouda (2020) and is
#' based on two two steps to reconstruct the epidemic curve. The approach is comparable to
#' our approach for Bavaria. More details can be found at https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Projekte_RKI/Nowcasting.html
#'
#' The steps are:
#'
#' 1. Impute missing date of disease onset using a missing-at-random assumption within 10 age-sex strata and reporting week.
#'
#' 2. Correct for right-truncation near the current "now" of data (in our case 2020-05-15)
#'    using a so called Nowcast procedure up to 30 days into the past.
#'
#' Note: With a now of "2020-05-17" the nowcast part is not relevant for the considered period
#' from 2020-02-23 and until 2020-04-21. Hence, only step 1 of the procedure - the missing data
#' imputation step - is of relevance. Approximately 50-60% of the reported cases have
#' a reported day of onset. However, only <10% have an explicit indication that they were symptomless
#' at the time of the report. Hence, the MAR assumption within the selected strata appears sufficiently
#' valid. This is further investigated in our sensitivity analyses for Bavaria.
#' https://www.stablab.stat.uni-muenchen.de/_assets/docs/nowcasting_covid19_bavaria.pdf
#'


now=lubridate::ymd("2020-09-24")
#' Pick date and estimated number of disease onsets that day as the only two columns
covid = read_tsv("./data/nowcasting_results_2020-09-24.csv") %>%
  rename(Date=date, Est_Onsets=nowcast_med) %>%
  dplyr::select(Date, Est_Onsets)


# Stop 21 days before now, in order to avoid any nowcasting effects.
covid_ts <- covid %>% filter(Date <= now - 21) %>% arrange(Date)

#' ### Segmented regression model for number of infections
#' We estimate a segmented regression model with three or four breakpoints based on a
#' quasi Poisson model for the number of new infections per day with expectation
#' $$
#' E( Y_t) = \exp \left( \beta_0 + \beta_1 t + \sum_{k=1}^K \gamma_k (t-CP_k)_+ \right),
#' $$
#' where $E(Y_t)$ is the expected number of new cases at time $t$, and $K$ is the number of break points. $x_+ = \max(x,0)$ is the positive part of $x$. The break points are used to partition the curve of new infections $Y_t$ into $K+1$ phases. These are characterized by different growth parameters.
#'
#' The breakpoints are found based on discrete optimization of the deviance for all
#' potential combinations of breakpoints. Those breakpoints are put as starting
#' values in the `segmented::segmented` function which can be used to obtain CIs for the
#' dates of the breakpoints.
#'
source("breakpoint_fun_onset.R")
# Get necessary data (date and  case numbers)

dat_bp =  covid_ts %>%
  select(date=Date, backpro=Est_Onsets) %>%
  filter(date<=lubridate::ymd("2020-05-01"))

# Estimate segmented regression with 3 breakpoints based on discrete optimization,
# can take several minutes, we load saved results
res_3_disc = estimate_bp_disc_optim(data = dat_bp, bp = 3)
save(res_3_disc, file = "./results/res_disc_opt_3bp_onset.RData")
load("./results/res_disc_opt_3bp_onset.RData")

# Estimate segmented regression with 4 Breakpoints based on discrete optimization,
# can take several minutes, we load saved results
res_4_disc = estimate_bp_disc_optim(data = dat_bp, bp = 4)
save(res_4_disc, file = "./results/res_disc_opt_4bp_onset.RData")
load("./results/res_disc_opt_4bp_onset.RData")

# Estimate segmented regression with 3 or 4 Breakpoints based on 'segmented' package with startpoints from
# discrete optimization
res_3_seg = estimate_bp_segmented(
  data           = dat_bp,
  bp             = 3,
  start_bp       = res_3_disc$bps,
  segmented_seed = 0520)
res_4_seg = estimate_bp_segmented(
  data           = dat_bp,
  bp             = 4,
  start_bp       = res_4_disc$bps-1,
  segmented_seed = 0520)

# Compare deviance and overdispersion of model with 3 and 4 breakpoints,
# model with 4 BPs appears to fits data better
round(
  c(
    dev_3_bp = summary(res_3_seg$segmented_model)$deviance,
    overdisp_3_bp = summary(res_3_seg$segmented_model)$dispersion,
    dev_4_bp = summary(res_4_seg$segmented_model)$deviance,
    overdisp_4_bp = summary(res_4_seg$segmented_model)$dispersion),
  2)

############################## Results #########################################
# Plot of estimated segmented regression model
pdf("results/res-bavaria-disease-onset.pdf", width = 9, height = 6)
res_4_seg$plot
dev.off()

# More info on estimated breakpoints
# Breakpoints + 95%-CIs
knitr::kable(res_4_seg$breakpoints, "html") %>%
  readr::write_file("results/disease-onset-change-points.html")
# Daily multiplicative change in infection numbers
knitr::kable(res_4_seg$coef, "html", digits = 3) %>%
  readr::write_file("results/disease-onset-factors.html")
