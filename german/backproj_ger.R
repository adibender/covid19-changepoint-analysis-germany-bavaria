#' ---
#' title: "Segmented regression of daily COVID-19 infections using back-projected epidemic curve: German data "
#' author: "Michael Höhle, Felix Günther, Helmut Küchenhoff, Andreas Bender"
#' date: "2020-10-23" new version with data unit May,1st
#' ---

#' # Segmented regression of daily COVID-19 infections using back-projected epidemic curve
#'
#' Authors: Michael Höhle (<hoehle@math.su.se>), Felix Günther, Helmut Küchenhoff, Andreas Bender<p>
#' ---


# Load packages
library(readxl)
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

#'  Load public RKI data containing the epidemic curve as Excel file available from
#' https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Projekte_RKI/Nowcasting_Zahlen.xlsx?__blob=publicationFile
#'
now <- as.Date("2020-09-25")
daten_file <- str_c("./data/Nowcasting_Zahlen-", as.character(now), ".xlsx")

#' Pick date and estimated number of disease onsets that day as the only two columns
covid <- readxl::read_xlsx(path = daten_file, sheet = "Nowcast_R") %>%
  rename(Date = "Datum des Erkrankungsbeginns",
         Est_Onsets = "Punktschätzer der Anzahl Neuerkrankungen (ohne Glättung)") %>%
  mutate(Date = as.Date(Date)) %>%
  dplyr::select(Date, Est_Onsets)

#' Manually read values by digitizing the plot in Fig. 2 of an der Heiden & Hamouda (2020)
#' in the range of 2020-02-19 and until 2020-03-09, using the [Engauge Digitizer](http://markummitchell.github.io/engauge-digitizer/) software.
#' We do this because the above Excel file only contains the time series from March 02 on.
covid_extra <- read_csv("./data/rki_nowcast.csv", locale=locale(decimal_mark=",")) %>%
  mutate(Date = round(x-5.5) + as.Date("2020-02-24") ) %>%
  rename(Est_Onsets_MA3 = Curve1) %>%
  select(Date, Est_Onsets_MA3)

# Derive raw values from filtered values (moving average of 3 days)
covid_extra$Est_Onsets <- covid_extra$Est_Onsets_MA3
for (i in 3:nrow(covid_extra)) {
  covid_extra$Est_Onsets[i] <- 3*covid_extra$Est_Onsets_MA3[i] - covid_extra$Est_Onsets[i-1] - covid_extra$Est_Onsets[i-2]
}

# Check if digitizing worked by comparing MA(3) smooth with the obtained values.
covid <- covid %>%
  mutate(Est_Onsets_MA3 = 1/3*(Est_Onsets + lag(Est_Onsets) + lag(Est_Onsets, n=2)))
# Merge the two time series
covid_ts <- full_join(covid, covid_extra, by="Date") %>% arrange(Date)
# Check that reconstruction is ok
covid_ts %>% filter(Date >= as.Date("2020-03-01") & Date <= as.Date("2020-03-09"))

# Use true value if possible, otherwise use the reconstructed value
covid_ts <- covid_ts %>%
  mutate(Est_Onsets = if_else(is.na(Est_Onsets.x), Est_Onsets.y, Est_Onsets.x))

# Stop at may,1st.

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
# Plot incubation period distribution
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
bp <- bp[epoch(bp) <= as.Date("2020-05-01") & epoch(bp) >= as.Date("2020-02-24"),]

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
dat_bp = bpdf %>% select(date=epoch, backpro = upperbound)
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
    dev_3_bp = summary(res_3_seg$segmented_model)$deviance,
    overdisp_3_bp = summary(res_3_seg$segmented_model)$dispersion,
    dev_4_bp = summary(res_4_seg$segmented_model)$deviance,
    overdisp_4_bp = summary(res_4_seg$segmented_model)$dispersion),
  2)


############################## Results #########################################
# Plot of estimated segmented regression model
pdf("results/res-ger-backprojections.pdf", width = 9, height = 6)
res_4_seg$plot
dev.off()

# More info on estimated breakpoints
# Breakpoints + 95%-CIs
knitr::kable(res_4_seg$breakpoints, "html") %>%
  readr::write_file("results/back-projection-change-points.html")
# Daily multiplicative change in infection numbers
knitr::kable(res_4_seg$coef, "html", digits = 3) %>%
  readr::write_file("results/back-projection-factors.html")
