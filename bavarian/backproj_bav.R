#' ---
#' title: "Segmented regression of daily COVID-19 infections using back-projected epidemic curve bavarian Data"
#' author: "Michael Höhle, Felix Günther, Helmut Küchenhoff, Andreas Bender"
#' date: "2020-10-23" new version with data unit May,1st
#' ---

#' # Segmented regression of daily COVID-19 infections using back-projected epidemic curve
#'
#' Authors: Michael Höhle (<hoehle@math.su.se>), Felix Günther and Helmut Küchenhoff, Andreas Bender<p>
#' ---

# Load packages
library(tidyverse)
library(segmented)
library(surveillance)
library(future)
library(future.apply)

#' ## Introduction
#'
#' This document contains a reproducible analysis to infer change-points for the back-projection
#' of the curve of daily disease onsets of reported Bavaian COVID-19 cases.
#'
#'
#'
#' For Bavaria we have access to the raw data and provide a thorough description
#' https://www.stablab.stat.uni-muenchen.de/_assets/docs/nowcasting_covid19_bavaria.pdf
#'
#'
#'
#'
#'
#'
#' https://www.stablab.stat.uni-muenchen.de/_assets/docs/nowcasting_covid19_bavaria.pdf
#'

#'
now=lubridate::ymd("2020-09-24")
#' Pick date and estimated number of disease onsets that day as the only two columns
covid = read_tsv("./data/nowcasting_results_2020-09-24.csv") %>%
  rename(Date=date, Est_Onsets=nowcast_med) %>%
  dplyr::select(Date, Est_Onsets)


 # Stop 21 days before now, in order to avoid any nowcasting effects.
covid_ts <- covid %>% filter(Date <= now - 21) %>% arrange(Date)


#' ## Backprojection of the epidemic curve
#'
#' Non-parametric back-projection as in Becker et al. (1991). The exposure time is
#' the relevant time scale to assess interventions.
#'
#' We take a literature based approach to deduce an incubation time distribution.
#'
#' Lauer et al. (2020) - log normal distribution - same as in Dehning et al. (2020)
#' Source: [Lauer et al. (2020)](https://www.ncbi.nlm.nih.gov/pubmed/32150748).
quantiles_incu <- data.frame(q=c(0.5, 0.975), value=c(5.1, 11.5))
# Fit distribution
f_target2 <- function(theta) {
  qlnorm(quantiles_incu$q[c(1,2)], theta[1], theta[2]) - quantiles_incu$value[c(1,2)]
}
incu_lnorm <- nleqslv::nleqslv(c(1,1), f_target2)$x
# Compare observed and fitted value
data.frame(
  q        = quantiles_incu$q,
  observed = quantiles_incu$value,
  fitted   = qlnorm( quantiles_incu$q, incu_lnorm[1], incu_lnorm[2]))

inc_pdf = data.frame(t = seq(0,15, length=1000)) %>%
  mutate(pdf = dlnorm(t, incu_lnorm[1], incu_lnorm[2]))

# Discretize incubation period distribution
cdf <- plnorm(seq(0,14,by=1), incu_lnorm[1], incu_lnorm[2])
pmf <- structure(c(0,diff(cdf)), names=seq(0,14,by=1))
# Normalize the discrete incubation period distribution
pmf <- pmf/sum(pmf)
df <- data.frame(days=as.numeric(names(pmf)), pmf=pmf)

#' The backprojection can be done using the function `surveillance::backprojNP`
#' The backprojected curve shows the number of infections per day and can be
#' compared to interventions similar to Werber et al. (2013) (https://doi.org/10.1093/aje/kwt069)

# Extract data from nowcast
sts_symp <- sts(
  epoch       = covid_ts$Date,
  observed    = matrix(covid_ts$Est_Onsets, ncol = 1, nrow = nrow(covid_ts)),
  epochAsDate = TRUE)

# Perform back projection with smoothing to adjust for weekday effects (k=6)
bp <- backprojNP(sts_symp, incu.pmf=pmf, control=list(k=6, eq3a.method="C"))
# Reduce to relevant subset of the data making it comparable with Dehning et al. (2020)
bp <- bp[epoch(bp) <= as.Date("2020-05-01") & epoch(bp) >= as.Date("2020-02-27"),]

bpdf <- bp %>%
  as.data.frame() %>%
  mutate(epoch_numeric=as.numeric(epoch), t = epoch - max(epoch) + 1)


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
source("breakpoint_fun.R")
# Get necessary data (date and backprojected case numbers)
dat_bp = bpdf %>%
  select(date=epoch, backpro = upperbound) %>%
  filter(date>=lubridate::ymd("2020-02-27"))
# Estimate segmented regression with 3 breakpoints based on discrete optimization,
# can take several minutes, we load saved results
res_3_disc = estimate_bp_disc_optim(data = dat_bp, bp = 3)
save(res_3_disc, file = "./results/res_disc_opt_3bp.RData")
load("./results/res_disc_opt_3bp.RData")

# Estimate segmented regression with 4 Breakpoints based on discrete optimization,
# can take several minutes, we load saved results
res_4_disc = estimate_bp_disc_optim(data = dat_bp, bp = 4)
save(res_4_disc, file = "./results/res_disc_opt_4bp.RData")
load("./results/res_disc_opt_4bp.RData")

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
  start_bp       = res_4_disc$bps,
  segmented_seed = 0520)

# Compare deviance and overdispersion of model with 3 and 4 breakpoints,
# model with 4 BPs appears to fits data better
round(
  c(
    dev_3_bp      = summary(res_3_seg$segmented_model)$deviance,
    overdisp_3_bp = summary(res_3_seg$segmented_model)$dispersion,
    dev_4_bp      = summary(res_4_seg$segmented_model)$deviance,
    overdisp_4_bp = summary(res_4_seg$segmented_model)$dispersion),
  2)


############################## Results #########################################
# Plot of estimated segmented regression model
### Use model with 3 change points

pdf("results/res-bavaria-backprojections.pdf", width = 9, height = 6)
res_3_seg$plot
dev.off()

# More info on estimated breakpoints
# Breakpoints + 95%-CIs
knitr::kable(res_3_seg$breakpoints, "html") %>%
  readr::write_file("results/back-projection-change-points.html")
# Daily multiplicative change in infection numbers
knitr::kable(res_3_seg$coef, "html", digits = 3) %>%
  readr::write_file("results/back-projection-factors.html")
