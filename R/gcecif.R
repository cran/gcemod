#' Generalized competing event model on the cumulative-incidence-ratio scale
#'
#' @description
#' Fits a generalized competing event (GCE) model on the \emph{cumulative
#' incidence} scale, estimating covariate effects on the ratio of cumulative
#' incidences \eqn{\rho(t) = F_1(t) / F_2(t)} at a landmark time \code{t}, where
#' \eqn{F_1} and \eqn{F_2} are the cumulative incidence functions of the event of
#' interest and the competing event. \eqn{\rho(t)} is the odds that a subject's
#' \emph{realized} event by time \code{t} is the event of interest rather than the
#' competing event, so \code{gcecif()} is a cumulative-incidence-scale companion
#' to the cause-specific (\code{\link{gcecox}}) and subdistribution
#' (\code{\link{gcefg}}) omega+ models and returns a compatible \code{"gcemod"}
#' object.
#'
#' Effects are estimated by a \emph{stacked pseudo-observation regression}: each
#' subject contributes one jackknife pseudo-observation for \eqn{F_1(t)} and one
#' for \eqn{F_2(t)} (both mean-type functionals of the Aalen-Johansen estimator,
#' for which conditional unbiasedness of pseudo-observations is established),
#' and a single log-link GEE is fit to the stacked pairs with an event-type
#' indicator. The main effects model \eqn{\log F_2(t \mid z)}; the event-type
#' interaction coefficients estimate \eqn{\log \rho(t \mid z) =
#' \log F_1(t \mid z) - \log F_2(t \mid z)} directly. Because each subject
#' appears twice, a cluster-robust (subject) sandwich variance is used, which
#' captures the dependence between a subject's two pseudo-observations induced
#' by the shared all-cause survivor. This is the cumulative-incidence-scale
#' analogue of the stacked Lunn-McNeil construction used on the two hazard
#' scales. A GCE risk score is formed from the interaction linear predictor, so
#' \code{exp(riskscore)} is a subject's \eqn{\rho} relative to the average
#' subject (see \code{\link{gce_riskscore}}, \code{\link{gce_scaling}}).
#'
#' @details
#' Because \eqn{F_k(t) = \int_0^t S(u^-)\,h_k(u)\,du} with \eqn{S} the all-cause
#' survivor, \eqn{\rho(t)} is a survival-weighted time-average of the
#' instantaneous cause-specific hazard ratio \eqn{\omega^{+}(u)}; under a
#' time-constant \eqn{\omega^{+}} it equals \eqn{\omega^{+}} exactly. The ratio
#' of cumulative incidences is a non-collapsible, marginal-over-time measure, so
#' \code{exp(coef)} is interpreted as a conditional association with the realized
#' event mix, not a mechanism-level per-cause effect.
#'
#' The estimating function is \emph{linear} in the two per-cause (mean-type)
#' pseudo-observations, so its conditional unbiasedness follows from the
#' standard pseudo-observation theory (Graw et al. 2009; Overgaard et al. 2017);
#' pseudo-observations of the nonlinear ratio itself are not used. The working
#' model comprises a log-linear model for \eqn{F_2(t \mid z)} (the nuisance main
#' effects, reported as \code{$coef_f2}) and a log-linear model for
#' \eqn{\rho(t \mid z)} (the target interactions, reported as \code{$coef}).
#'
#' The estimand is defined for a genuine competing-risks setting, in which both
#' \eqn{F_1(t)} and \eqn{F_2(t)} are appreciable at the landmark \code{t}. If
#' the competing event is effectively absent, or \code{t} precedes appreciable
#' competing-event incidence, \eqn{\log \rho} is weakly identified and the fit
#' may fail or return very large standard errors; this should be read as a
#' signal that the competing-risks premise does not hold.
#'
#' @param Time numeric vector of follow-up times.
#' @param Ind a data frame of event indicators: column 1 = event of interest,
#'   column 2 = competing event (\code{0}/\code{1}).
#' @param Cov a data frame of covariates (factors are expanded automatically).
#' @param M number of bins for the \eqn{\rho} risk plot (default \code{5}).
#' @param t landmark time at which \eqn{\rho} is evaluated.
#' @param conf.level confidence level (default \code{0.95}).
#'
#' @return An object of class \code{"gcemod"} (subclass \code{"gcecif"}) with:
#'   \code{coef} (log \eqn{\rho} ratios), \code{result} (ratio table with CIs and
#'   p-values), \code{omnibus} (Wald test that all \eqn{\rho} ratios equal 1),
#'   \code{riskscore} (interaction linear predictor, mean-centered),
#'   \code{scaling} (covariate means/SDs and per-SD effects), \code{F1},
#'   \code{F2}, \code{rho0} (marginal \eqn{F_1(t)/F_2(t)}), \code{rho}
#'   (per-subject fitted \eqn{\rho}), \code{pseudo} (a two-column data frame of
#'   per-cause pseudo-observations), \code{coef_f2} (the log-linear
#'   \eqn{F_2(t \mid z)} nuisance model), \code{converged}, and the design
#'   matrix, time, and cause coding used by \code{\link{gce_tidy}}.
#'
#' @references
#' Andersen PK, Klein JP, Rosthoj S (2003). Generalized linear models for
#' correlated pseudo-observations, with applications to multi-state models.
#' \emph{Biometrika} 90:15-27.
#'
#' Graw F, Gerds TA, Schumacher M (2009). On pseudo-values for regression
#' analysis in competing risks models. \emph{Lifetime Data Analysis} 15:241-255.
#'
#' Overgaard M, Parner ET, Pedersen J (2017). Asymptotic theory of generalized
#' estimating equations based on jack-knife pseudo-observations. \emph{The
#' Annals of Statistics} 45:1988-2015.
#'
#' Overgaard M, Parner ET, Pedersen J (2018). Estimating the variance in a
#' pseudo-observation scheme with competing risks. \emph{Scandinavian Journal
#' of Statistics} 45:923-940.
#'
#' @examples
#' data(prostate)
#' Ind <- data.frame(event = as.integer(prostate$status == 1),
#'                   competing = as.integer(prostate$status == 2))
#' Cov <- prostate[, c("age", "psa", "gleason", "t2b", "comorbidity")]
#' fit <- gcecif(prostate$time, Ind, Cov, t = 10)
#' summary(fit)
#' gce_tidy(fit)
#' @seealso \code{\link{gcecox}}, \code{\link{gcefg}}, \code{\link{gce_tidy}},
#'   \code{\link{gce_riskscore}}
#' @export
gcecif <- function(Time, Ind, Cov, M = 5, t, conf.level = 0.95) {
  cl <- match.call()
  chk <- .gce_check(Time, Ind, Cov, M, t)
  Time <- chk$Time; Ind <- chk$Ind; Cov <- chk$Cov
  n <- length(Time)

  CA <- as.integer(Ind[[1]]); CM <- as.integer(Ind[[2]])
  cause <- integer(n); cause[CM == 1] <- 2L; cause[CA == 1] <- 1L

  ut <- sort(unique(Time[cause > 0 & Time <= t]))
  if (length(ut) < 2L) stop("Too few events at or before 't'.", call. = FALSE)
  Y  <- vapply(ut, function(u) sum(Time >= u), numeric(1))
  d1 <- vapply(ut, function(u) sum(Time == u & cause == 1L), numeric(1))
  d2 <- vapply(ut, function(u) sum(Time == u & cause == 2L), numeric(1))

  .cif <- function(Y, d1, d2) {
    haz   <- ifelse(Y > 0, (d1 + d2) / Y, 0)
    Sprev <- c(1, cumprod(1 - haz)[-length(haz)])
    c(sum(Sprev * ifelse(Y > 0, d1 / Y, 0)),
      sum(Sprev * ifelse(Y > 0, d2 / Y, 0)))
  }
  f0 <- .cif(Y, d1, d2); F1 <- f0[1]; F2 <- f0[2]
  if (F1 <= 0 || F2 <= 0)
    stop("A cumulative incidence is zero at 't'; choose a larger 't' or check ",
         "that a genuine competing risk is present.", call. = FALSE)

  ## leave-one-out CIFs
  atrisk <- outer(Time, ut, ">=")
  idx    <- match(Time, ut)
  F1i <- F2i <- numeric(n)
  for (i in seq_len(n)) {
    Yi <- Y - atrisk[i, ]
    d1i <- d1; d2i <- d2
    if (!is.na(idx[i])) {
      if (cause[i] == 1L) d1i[idx[i]] <- d1i[idx[i]] - 1
      if (cause[i] == 2L) d2i[idx[i]] <- d2i[idx[i]] - 1
    }
    fi <- .cif(Yi, d1i, d2i); F1i[i] <- fi[1]; F2i[i] <- fi[2]
  }

  ## per-cause (mean-type) jackknife pseudo-observations
  ps1 <- n * F1 - (n - 1) * F1i
  ps2 <- n * F2 - (n - 1) * F2i

  ## stacked log-link GEE: main effects = log F2(t|z) model, event-type
  ## interactions = log rho(t|z) model, cluster-robust on subject
  X  <- stats::model.matrix(~ ., data = Cov)   # intercept + covariates
  p0 <- ncol(X)
  start <- c(log(F2), rep(0, p0 - 1L), log(F1) - log(F2), rep(0, p0 - 1L))
  sf <- .gce_stacked(X, ps1, ps2, start)
  theta <- sf$theta
  names(theta) <- c(colnames(X), paste0(colnames(X), ":rho"))
  main_idx <- seq_len(p0)
  int_idx  <- p0 + seq_len(p0)

  beta_int <- theta[int_idx]                          # (log rho at z = 0, effects)
  V_int    <- sf$V[int_idx, int_idx, drop = FALSE]
  se_int   <- sqrt(diag(V_int))
  rho_pred <- exp(as.numeric(X %*% beta_int))         # fitted rho(t | z_i)

  ## intercept vs covariate coefficients (within the interaction block)
  keepc <- colnames(X) != "(Intercept)"
  cn    <- colnames(X)[keepc]
  bcov  <- beta_int[keepc]; secov <- se_int[keepc]
  Vcov  <- V_int[keepc, keepc, drop = FALSE]
  names(bcov) <- names(secov) <- cn
  zc <- stats::qnorm(1 - (1 - conf.level) / 2)

  result <- cbind(`exp(coef) (rho ratio)` = round(exp(bcov), 5),
                  `lower .95` = round(exp(bcov - zc * secov), 5),
                  `upper .95` = round(exp(bcov + zc * secov), 5),
                  `P-value`   = round(2 * stats::pnorm(-abs(bcov / secov)), 5))
  rownames(result) <- cn

  ## risk score: mean-centered interaction linear predictor on the covariates
  Xc  <- X[, keepc, drop = FALSE]
  ctr <- colMeans(Xc); sds <- apply(Xc, 2, stats::sd)
  coef_raw <- bcov; coef_std <- coef_raw * sds
  names(ctr) <- names(sds) <- names(coef_raw) <- cn
  riskscore <- as.numeric(sweep(Xc, 2, ctr, "-") %*% coef_raw)
  normCER   <- as.numeric(scale(riskscore))

  ## omnibus Wald test (all covariate rho-ratios = 1)
  omnibus <- c(statistic = NA_real_, df = NA_real_, p.value = NA_real_)
  if (all(is.finite(bcov)) && all(is.finite(Vcov))) {
    stat <- tryCatch(as.numeric(t(bcov) %*% solve(Vcov, bcov)),
                     error = function(e) NA_real_)
    if (is.finite(stat)) {
      df <- length(bcov)
      omnibus <- c(statistic = stat, df = df,
                   p.value = stats::pchisq(stat, df, lower.tail = FALSE))
    }
  }

  scaling <- data.frame(term = cn, mean = as.numeric(ctr), sd = as.numeric(sds),
                        coef_perSD = as.numeric(coef_std),
                        rho_ratio_perSD = exp(as.numeric(coef_std)),
                        row.names = NULL, stringsAsFactors = FALSE)

  ## rho vs risk-score diagnostic (M bins)
  qb    <- seq(0, 1, length.out = M + 1)
  rs_q  <- stats::quantile(normCER,  probs = qb[-1], na.rm = TRUE)
  rho_q <- stats::quantile(rho_pred, probs = qb[-1], na.rm = TRUE)
  ratioplot <- .gce_scatter(rs_q, rho_q, "Risk score", expression(rho))

  structure(list(
    call = cl,
    method = "cumulative-incidence ratio (stacked pseudo-observation)",
    coef = coef_raw, result = result,
    coef_raw = coef_raw, coef_std = coef_std,
    center = ctr, scale = sds, scaling = scaling,
    omnibus = omnibus, ratioplot = ratioplot,
    riskscore = riskscore, normCER = normCER,
    F1 = F1, F2 = F2, rho0 = F1 / F2, rho = rho_pred,
    pseudo = data.frame(F1 = ps1, F2 = ps2),
    coef_f2 = stats::setNames(theta[main_idx], colnames(X)),
    vcov = sf$V, converged = sf$converged,
    time = Time, cause = cause, design = Xc,
    covnames = cn, covariate.names = names(Cov),
    n = n, n_event = sum(cause == 1L), n_competing = sum(cause == 2L),
    t = t, estimator = "stacked", conf.level = conf.level
  ), class = c("gcecif", "gcemod"))
}

## Stacked log-link GEE for the two per-cause pseudo-observation vectors.
## Gauss-Newton on the working-independence estimating equation
##   U(theta) = sum_r mu_r w_r (y_r - mu_r) = 0,   mu_r = exp(w_r' theta),
## with backtracking on the residual sum of squares, and a cluster-robust
## (subject) sandwich variance: each subject's two rows are grouped before the
## outer product, capturing the within-subject pseudo-observation dependence.
.gce_stacked <- function(X, ps1, ps2, start) {
  n <- nrow(X); p0 <- ncol(X); p2 <- 2L * p0
  W  <- rbind(cbind(X, X),                        # cause-1 rows: main + interaction
              cbind(X, matrix(0, n, p0)))         # cause-2 rows: main only
  y  <- c(ps1, ps2)
  id <- rep.int(seq_len(n), 2L)
  rss <- function(th) {
    m <- exp(as.numeric(W %*% th)); s <- sum((y - m)^2)
    if (is.finite(s)) s else Inf
  }
  theta <- start; f_cur <- rss(theta); step_max <- Inf
  for (it in 1:200) {
    mu <- exp(as.numeric(W %*% theta))
    U  <- as.numeric(crossprod(W, mu * (y - mu)))
    A  <- crossprod(W * mu)
    step <- tryCatch(
      as.numeric(solve(A + diag(1e-8 * mean(diag(A)), p2), U)),
      error = function(e) NULL)
    if (is.null(step))
      stop("Stacked pseudo-observation fit failed (singular information); ",
           "check that both events carry appreciable incidence at 't'.",
           call. = FALSE)
    h <- 1
    repeat {
      f_new <- rss(theta + h * step)
      if (f_new <= f_cur + 1e-12) break
      h <- h / 2
      if (h < 1e-10) break
    }
    theta <- theta + h * step; f_cur <- rss(theta)
    step_max <- max(abs(h * step))
    if (step_max < 1e-9) break
  }
  mu <- exp(as.numeric(W %*% theta))
  Gi <- rowsum(W * (mu * (y - mu)), id)           # per-subject score contributions
  A  <- crossprod(W * mu)
  Ai <- solve(A + diag(1e-12 * mean(diag(A)), p2))
  V  <- Ai %*% crossprod(Gi) %*% Ai
  list(theta = theta, V = V, converged = step_max < 1e-9, iter = it)
}
