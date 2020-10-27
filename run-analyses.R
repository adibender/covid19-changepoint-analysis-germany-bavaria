# This is the main file for the analysis that reruns all the code necessary to
# obtain the analyses presented in
# "Analysis of the early Covid-19 epidemiccurve in Germany by regression modelswith change points"


# set multicore options
library(parallel)
plan("multicore")
options(mc.cores = 6)

# Analyses for Bavarian Data
setwd("bavarian")
dir.create("results") # results will be saved here
dir.create("graphics")
source("graphic_bav.R",echo=TRUE)
source("backproj_bav.R", echo = TRUE)
source("disease_onset_bav.R", echo = TRUE)


# Analyses for Germany
setwd("../german/")
dir.create("results") # results will be saved here
dir.create("graphics")
source("graphic_ger.R",echo=TRUE)
source("backproj_ger.R", echo = TRUE)
source("disease_onset_ger.R", echo = TRUE)
