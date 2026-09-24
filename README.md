# gcemod

Generalized competing event (GCE) modeling of **omega+** --- the ratio of the
hazard for an event of interest to the hazard for a competing event --- on the
cause-specific (Cox) and subdistribution (Fine-Gray) scales, with confidence
intervals and p-values from the Lunn-McNeil stacked-data approach.

## What it does

1. **Estimate covariate effects on omega+** for Cox and Fine-Gray models, with
   95% CIs, Wald p-values, and an omnibus test from the Lunn-McNeil stacked-data
   approach (`gcecox()`, `gcefg()`, `lunnmcneil()`).
2. **Build a risk score** from the model linear predictor (`$riskscore`).
3. **Find cutpoints** that maximize the omega+ difference between groups
   (`gce_cutpoints(method = "optimal")`).
4. **Alligator plots** --- cumulative incidence of the event of interest and the
   competing event within risk groups (`gce_alligator()`).
5. **Calibration plots** --- predicted vs. observed omega+ by rank, with a
   user-specified number of groups (`gce_calibration()`).
6. **Compare effects across events** --- primary, competing, total (composite),
   and relative (omega+) hazard ratios in one table
   (`gce_tidy(fit, effects = TRUE)`).
7. **Score new subjects** relative to the average subject, and expose the
   covariate means/SDs behind the risk score (`gce_riskscore()`, `gce_scaling()`).
8. **Cumulative-incidence-ratio GCE model** --- covariate effects on
   rho(t) = F1(t)/F2(t), the odds that a realized event is the event of interest,
   by pseudo-observation regression (`gcecif()`). Returns a `gcemod` object, so it
   works with `gce_tidy()`, `gce_scaling()`, and `gce_riskscore()` like the
   `gcecox()` / `gcefg()` fits.

Per-subject omega = omega+/(1+omega+) is kept as a descriptive quantity
(`$omega`) so model-predicted values can be compared with observed ones.
Discrimination is available via `gce_performance()` (cause-specific C-index).

## Install and check (from source)

```r
# one time: build the bundled example datasets
source("data-raw/make_sample.R")      # writes data/hn.rda and data/prostate.rda

# regenerate man/ and NAMESPACE from the roxygen headers
devtools::document()

devtools::check()
```

## Quick start

```r
library(gcemod)
data(hn)   # head-and-neck cohort: status 1 = recurrence, 2 = death w/o recurrence
Ind <- data.frame(event     = as.integer(hn$status == 1),
                  competing = as.integer(hn$status == 2))
Cov <- hn[, c("age", "smoker", "t_cat", "n_cat", "p16")]

fit  <- gcecox(hn$time, Ind, Cov, M = 5, t = 5)   # or gcefg(...) / gcecif(...)
summary(fit)

cif  <- gcecif(hn$time, Ind, Cov, t = 5)   # cumulative-incidence-ratio scale
gce_tidy(cif)                              # same tidy interface as gcecox/gcefg

gce_tidy(fit, effects = TRUE)   # primary / competing / total / relative HRs
cuts <- gce_cutpoints(fit, groups = 3, method = "optimal")
gce_alligator(fit, groups = cuts)
gce_calibration(fit, which = "omegaplus", groups = 5)
```

See `vignette("gcemod")` for the full workflow.

## References

Lunn M, McNeil D (1995) Applying Cox regression to competing risks.
*Biometrics* 51:524-32.

Carmona R, et al. (2014) Validated competing event model for the stage I-II
endometrial cancer population. *Int J Radiat Oncol Biol Phys* 89:888-98.

Carmona R, et al. (2016) Improved method to stratify elderly patients with
cancer at risk for competing events. *J Clin Oncol* 34:1270-77.

Mell LK, et al. (2019) Nomogram to predict the benefit of intensive treatment
for locoregionally advanced head and neck cancer. *Clin Cancer Res* 25:7078-7088.

Zakeri K, et al. (2020) Predictive classifier for intensive treatment of head
and neck cancer. *Cancer* 126:5263-5273.

Mell LK, et al. (2024) Effects of androgen deprivation therapy on prostate
cancer outcomes according to competing event risk. *Eur Urol* 85:373-381.
