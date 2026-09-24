# gcemod 0.3.0

First CRAN release.

* **`gce_tidy(effects = TRUE)` covariate-effects table.** `gce_tidy()` has an
  `effects` argument. With `effects = TRUE` it returns a comparison table that,
  for each covariate, juxtaposes the hazard ratio for the primary event, the
  competing event, the total (composite) event, and the relative hazard ratio
  (the omega+ ratio) from the GCE model, each with a confidence interval. The
  individual/total columns are estimated on the model's own scale --- cause-
  specific Cox for `gcecox()`, Fine-Gray subdistribution for `gcefg()` --- so the
  reader can see whether a covariate acts on the event of interest, the competing
  event, both, or their balance. `gcecox()` and `gcefg()` retain the model design
  matrix on the fitted object (`$design`) to support this.
* **`gcecif()` cumulative-incidence-ratio GCE model.** A third GCE estimator
  alongside `gcecox()` and `gcefg()`, on the cumulative-incidence scale:
  covariate effects on the ratio of cumulative incidences rho(t) = F1(t)/F2(t)
  at a landmark time `t` --- the odds that a subject's realized event by time t
  is the event of interest, a cumulative-incidence-scale analogue of omega+. Uses
  jackknife pseudo-observation regression of the Aalen-Johansen estimator with a
  robust variance; `method = "logit"` (default) fits the bounded share via a
  logit-link GEE, `method = "logratio"` fits log(F1/F2) via an identity-link GEE.
  Returns a `gcemod` object, so it works with `gce_tidy()`, `gce_scaling()`, and
  `gce_riskscore()` like the other fits. The full methodological development is
  described in a companion manuscript (in preparation).
* **Reported coefficients on the natural scale by default; risk score always
  normalized.** `gcecox()`, `gcefg()`, and `lunnmcneil()` report omega+ ratios
  per natural covariate unit by default; `standardize = TRUE` reports them per
  1 SD instead (a reparameterization, so p-values and the omnibus test are
  identical). The GCE risk score is *always* the normalized, mean-centered linear
  predictor, so the average patient has risk score 0 and cutpoints / plots /
  calibration are unaffected by the reporting choice.
* **New `gce_scaling()` and `gce_riskscore()`.** `gce_scaling()` returns each
  covariate's mean, SD, and per-SD omega+ ratio. `gce_riskscore(fit, newdata)`
  scores new subjects with the stored means/SDs and coefficients, returning the
  risk score and `omegaplus_rel` = the subject's omega+ relative to the average
  subject in the fitting data (with absolute omega+/omega for cause-specific
  models). `$scaling`, `$center`, `$scale`, `$coef_raw`, and `$coef_std` are
  stored on the fitted object.
* **Cutpoint selection on the log-omega+ scale (no tail chasing).** The optimal
  cutpoint criterion maximizes between-group separation in *log* omega+ (the
  risk score / linear predictor) rather than on the exponential omega+ scale.
  Because that is the additive scale of the joint omega+ contrast, the search is
  no longer dominated by a few extreme high-risk subjects and returns balanced,
  well-separated groups. `criterion = "cif"` (joint Gray statistic across both
  events) is available as an alternative.
* **Reported test is now the omega+ group contrast.** `gce_cutpoints()$test`
  reports the Lunn-McNeil test that the ratio of cause-specific (or
  subdistribution) hazards differs between groups, with the univariable omega+
  ratio estimate and CI -- the GCE test that separate per-event Gray/log-rank
  tests miss. The per-event Gray tests are also returned (`$individual_gray`) so
  the difference is visible, along with a selection-adjusted permutation p-value
  (`$p.perm`, via `nperm > 0`).
* **Alligator plot (`gce_alligator()`).** Draws the cumulative incidence of the
  event of interest and the competing event on a single panel, with risk groups
  distinguished by color and event types by line type (event of interest solid,
  competing event dashed). By default it forms 2 groups at the optimal omega+
  cutpoint (`method = "optimal"`, log-scale; `criterion` selectable) and
  annotates the cutpoint value and group sizes. The x-axis is truncated at the
  readout time `t` with 1-unit ticks; the y-axis auto-caps at the smallest of
  25/50/75/100% covering the peak incidence, with 5%/10% ticks (override with
  `ymax`), so low-incidence plots are not compressed; and numbers at risk are
  shown in an aligned panel below the plot (via `patchwork`). Groups of 2-3 are
  labeled Low/High or Low/Intermediate/High by default.
* **Two bundled example datasets (`hn`, `prostate`).** Two competing-risks
  example datasets show opposite GCE regimes: `hn` (head and neck, n = 1000) is
  recurrence-dominant (omega+ > 1), while `prostate` (n = 1000) is
  competing-death-dominant (omega+ < 1). Both code `status` as 0 = censored,
  1 = event of interest, 2 = competing, with `time` in years and sequential
  subject IDs. They are provided for illustration only and do not represent real
  patients. Categorical covariates (stage, site, Gleason, performance status)
  are factors. All examples and the vignette use these datasets.

# gcemod 0.2.0

Refocused the package on modeling omega+ (event of interest vs. competing event).

## Breaking changes

* **Removed the omega-ratio regression.** Covariate effects on omega (event of
  interest vs. all events) are no longer modeled.
* Per-subject omega = omega+/(1+omega+) is **kept** as a descriptive quantity
  (`$omega`) for comparing model-predicted with observed values.
* `gce_tidy()` no longer takes `which` (omega+ only). `gce_cutpoints()` now
  targets omega+ separation rather than a log-rank split.

## New features

* **Fine-Gray Lunn-McNeil.** `lunnmcneil()` gains `type = "finegray"`: a stacked,
  risk-set-weighted (`survival::finegray`) model giving omega+ CIs and p-values
  on the subdistribution scale from one joint model, replacing the bootstrap in
  `gcefg()`.
* `gce_cutpoints(method = "optimal")` finds cutpoints maximizing the between-group
  omega+ separation (greedy maximum between-group sum of squares).
* `gce_alligator()` --- cumulative incidence of the event of interest and the
  competing event within risk groups, with Gray's test.
* `gce_calibration()` --- predicted vs. observed omega+ (or omega) by rank group,
  user-specified number of groups.

## Robustness

* Constant covariates are dropped up front with a warning instead of causing
  a downstream error.

# gcemod 0.1.0

First release. `gcemod` is a successor to `gcerisk`.

## Changes relative to `gcerisk`

* **Variance for the omega ratio.** Where `gcerisk` estimates the variance of
  `log(omega ratio)` by sampling independent normal deviates for the
  cause-specific and all-cause coefficients and taking the variance of their
  difference, `gcemod` uses a case-resampling (nonparametric) bootstrap that
  refits the models on resampled subjects. (This omega-ratio regression was
  removed in 0.2.0.)

* **Hypothesis testing added (Lunn-McNeil).** The omega-plus ratio
  (event of interest vs. competing event) is estimated from a single
  stacked/augmented Cox model in the style of Lunn & McNeil (1995). The
  event-type interaction coefficient is exactly `log(omega+ ratio)` and carries
  a model-based (robust) standard error, a per-covariate Wald test, and an
  omnibus multivariate Wald test of no differential covariate effect across
  event types. Because omega is a monotone transform of omega-plus
  (`omega = omega+/(1 + omega+)`), this test applies to both.

* **Cumulative hazard at time `t`.** `gcerisk` estimates cumulative hazards by
  regressing `-log(S(t))` on time with `lm()` and extrapolating, which forces a
  linear cumulative hazard. `gcemod` reads the Nelson-Aalen cumulative hazard
  (cause-specific model) or the subdistribution cumulative hazard derived from
  the Fine-Gray CIF (`-log(1 - F(t))`) directly at the requested time via step
  functions.

* **Edge cases.** Formula construction now uses `reformulate()` instead of manual
  string pasting, which failed for a single covariate (the `1:(p-1)` loop). NA
  handling (complete-case with a warning), single-covariate models, and
  non-convergent bootstrap replicates (skipped with a warning) are handled.

## New features

* `lunnmcneil()` exposes the stacked-model estimator and tests directly.
* `gce_cutpoints()` derives risk-stratification groups from a GCE risk score
  (quantile groups, or an optimal two-group split by log-rank scan).
* `gce_performance()` reports the cause-specific concordance (C-index) and a
  grouped calibration table for the event of interest.
* `print()`, `summary()`, and `gce_tidy()` methods for compact, report-ready
  output.
