## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4)
library(gcemod)

## ----data---------------------------------------------------------------------
data(hn)
Ind <- data.frame(event     = as.integer(hn$status == 1),   # recurrence
                  competing = as.integer(hn$status == 2))   # death w/o recurrence
Cov <- hn[, c("age", "smoker", "t_cat", "n_cat", "p16")]

## ----cox----------------------------------------------------------------------
fit <- gcecox(hn$time, Ind, Cov, M = 5, t = 5)
summary(fit)

## ----fg, eval = FALSE---------------------------------------------------------
# fit_fg <- gcefg(hn$time, Ind, Cov, M = 5, t = 5)
# summary(fit_fg)

## ----rs-----------------------------------------------------------------------
head(data.frame(riskscore = fit$riskscore,
                omegaplus = fit$omegaplus, omega = fit$omega))

## ----cut----------------------------------------------------------------------
cuts <- gce_cutpoints(fit, groups = 2, method = "optimal")
cuts$summary
cuts$cutpoints

## ----alligator, fig.height=4--------------------------------------------------
gce_alligator(fit, groups = cuts)

## ----calib--------------------------------------------------------------------
gce_calibration(fit, which = "omegaplus", groups = 5)

## ----perf---------------------------------------------------------------------
gce_performance(fit)

## ----prostate, eval = FALSE---------------------------------------------------
# data(prostate)
# Ind_p <- data.frame(event     = as.integer(prostate$status == 1),
#                     competing = as.integer(prostate$status == 2))
# Cov_p <- prostate[, c("age", "psa", "gleason", "t2b", "comorbidity")]
# 
# fit_p <- gcecox(prostate$time, Ind_p, Cov_p, M = 5, t = 5)
# summary(fit_p)
# gce_alligator(fit_p, groups = 2)

