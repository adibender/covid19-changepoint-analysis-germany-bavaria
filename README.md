# Analysis of the early Covid-19 epidemic curve in Germany by regression models with change points

**The code in this repository refers to an initial analysis provided in the preprint. The revised code is available here: https://zenodo.org/record/4449816**


This repository contains code and data needed in order to reproduce results presented in "Analysis of the early Covid-19 epidemic curve in Germany by regression models with change points".

It contains two folders one with data and code for Bavaria, one for Germany:

1. **`bavarian`**:
  - **`data`**: Data that contains the raw data for the analysis + additional data for Figure 1
  - `graphic_bav.R` : graphic for comparison of the three curves (reported, disease onset, back-projection)
  - `backproj_bav.R`: estimates the back-projection and conducts change-point analysis for the infections
  - `disease_onset_bav.R`: conducts change-point analysis for the (raw) disease onset data (without back-projection)
  - `breakpoint_fun.R` and `breakpoint_fun_onset.R`: These script contain the functions that calculate the breakpoints, which are sourced in the `backproj_bav.R` and `disease_onset_bav.R`


2. **`german`**: Equivalent to contents of **`bavarian`** folder except using data for whole Germany.


The whole analysis can be reproduced in one go by running

```r
source("run-analyses.R")
```
