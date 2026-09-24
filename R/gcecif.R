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
#' Effects are estimated by jackknife pseudo-observation regression. Two
#' equivalent parameterizations are provided (both estimate the effect on
#' \eqn{\log \rho(t)}): \code{method = "logit"} (default) regresses pseudo-values
#' of the bounded share \eqn{\pi_1 = F_1/(F_1+F_2)} by a logit-link GEE, and
#' \code{method = "logratio"} regresses pseudo-values of \eqn{\log \rho} directly
#' by an identity-link GEE (ordinary least squares). Both use a robust (sandwich)
#' variance. A GCE risk score is formed from the model linear predictor, so
#' \code{exp(riskscore)} is a subject's \eqn{\rho} relative to the average subject
#' (see \code{\link{gce_riskscore}}, \code{\link{gce_scaling}}).
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
#' The estimand is defined for a genuine competing-risks setting, in which both
#' \eqn{F_1(t)} and \eqn{F_2(t)} are appreciable at the landmark \code{t}. If the
#' competing event is effectively absent, or \code{t} precedes appreciable
#' competing-event incidence, \eqn{\log \rho} is weakly identified; the optional
#' \code{ridge} penalty stabilizes such fits.
#'
#' @param Time numeric vector of follow-up times.
#' @param Ind a data frame of event indicators: column 1 = event of interest,
#'   column 2 = competing event (\code{0}/\code{1}).
#' @param Cov a data frame of covariates (factors are expanded automatically).
#' @param M number of bins for the \eqn{\rho} risk plot (default \code{5}).
#' @param t landmark time at which \eqn{\rho} is evaluated.
#' @param method estimation scale: \code{"logit"} (default; logit-link GEE on the
#'   share) or \code{"logratio"} (identity-link GEE on \eqn{\log\rho}).
#' @param ridge non-negative ridge penalty on the covariate coefficients (not the
#'   intercept), used only with \code{method = "logit"} to stabilize weakly
#'   identified fits. Default \code{0}.
#' @param conf.level confidence level (default \code{0.95}).
#'
#' @return An object of class \code{"gcemod"} (subclass \code{"gcecif"}) with:
#'   \code{coef} (log \eqn{\rho} ratios), \code{result} (ratio table with CIs and
#'   p-values), \code{omnibus} (Wald test that all \eqn{\rho} ratios equal 1),
#'   \code{riskscore} (linear predictor), \code{scaling} (covariate means/SDs and
#'   per-SD effects), \code{F1}, \code{F2}, \code{rho0} (baseline \eqn{\rho} at the
#'   mean covariate profile), \code{rho} (per-subject fitted \eqn{\rho}),
#'   \code{pseudo} (per-subject pseudo-values), and the design matrix, time, and
#'   cause coding used by \code{\link{gce_tidy}}.
#'
#' @references
#' Andersen PK, Klein JP, Rosthoj S (2003). Generalized linear models for
#' correlated pseudo-observations, with applications to multi-state models.
#' \emph{Biometrika} 90:15-27.
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
gcecif <- function(Time, Ind, Cov, M = 5, t,
                   method = c("logit", "logratio"),
                   ridge = 0, conf.level = 0.95) {
  cl <- match.call()
  method <- match.arg(method)
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

  X   <- stats::model.matrix(~ ., data = Cov)   # intercept + covariates
  eps <- 1e-9
  if (method == "logratio") {
    theta <- log(F1) - log(F2)
    ps    <- n * theta - (n - 1) * (log(F1i + eps) - log(F2i + eps))
    fitr  <- .gce_ols(X, ps)
    rho_pred <- exp(as.numeric(X %*% fitr$beta))
  } else {
    pi0 <- F1 / (F1 + F2)
    ps  <- n * pi0 - (n - 1) * (F1i / (F1i + F2i + eps))
    fitr <- .gce_logit(X, ps, ridge = ridge)
    pi_pred <- 1 / (1 + exp(-as.numeric(X %*% fitr$beta)))
    rho_pred <- pi_pred / (1 - pi_pred)
  }
  beta <- fitr$beta; se <- fitr$se; V <- fitr$V
  names(beta) <- names(se) <- colnames(X)

  ## intercept vs covariate coefficients
  keepc <- colnames(X) != "(Intercept)"
  cn    <- colnames(X)[keepc]
  bcov  <- beta[keepc]; secov <- se[keepc]
  Vcov  <- V[keepc, keepc, drop = FALSE]
  zc    <- stats::qnorm(1 - (1 - conf.level) / 2)

  result <- cbind(`exp(coef) (rho ratio)` = round(exp(bcov), 5),
                  `lower .95` = round(exp(bcov - zc * secov), 5),
                  `upper .95` = round(exp(bcov + zc * secov), 5),
                  `P-value`   = round(2 * stats::pnorm(-abs(bcov / secov)), 5))
  rownames(result) <- cn

  ## risk score: mean-centered linear predictor on the covariates
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
    call = cl, method = sprintf("cumulative-incidence ratio (%s)", method),
    coef = coef_raw, result = result,
    coef_raw = coef_raw, coef_std = coef_std,
    center = ctr, scale = sds, scaling = scaling,
    omnibus = omnibus, ratioplot = ratioplot,
    riskscore = riskscore, normCER = normCER,
    F1 = F1, F2 = F2, rho0 = F1 / F2, rho = rho_pred, pseudo = ps,
    time = Time, cause = cause, design = Xc,
    covnames = cn, covariate.names = names(Cov),
    n = n, n_event = sum(cause == 1L), n_competing = sum(cause == 2L),
    t = t, estimator = method, ridge = ridge, conf.level = conf.level
  ), class = c("gcecif", "gcemod"))
}

## identity-link GEE (OLS) with robust HC0 variance
.gce_ols <- function(X, y) {
  B <- solve(crossprod(X)); beta <- as.numeric(B %*% crossprod(X, y))
  r <- as.numeric(y - X %*% beta)
  V <- B %*% (t(X) %*% (X * r^2)) %*% B
  list(beta = beta, se = sqrt(diag(V)), V = V)
}

## logit-link GEE via IRLS (optional ridge on non-intercept), robust variance
.gce_logit <- function(X, y, ridge = 0) {
  p <- ncol(X); beta <- rep(0, p)
  Pen <- diag(c(0, rep(ridge, p - 1)), p, p)
  for (it in 1:100) {
    eta <- as.numeric(X %*% beta); mu <- 1 / (1 + exp(-eta))
    W   <- pmax(mu * (1 - mu), 1e-6); z <- eta + (y - mu) / W
    XtW <- t(X * W)
    bnew <- solve(XtW %*% X + Pen, XtW %*% z + Pen %*% beta)
    if (max(abs(bnew - beta)) < 1e-8) { beta <- as.numeric(bnew); break }
    beta <- as.numeric(bnew)
  }
  eta <- as.numeric(X %*% beta); mu <- 1 / (1 + exp(-eta))
  W <- pmax(mu * (1 - mu), 1e-6); r <- y - mu
  D <- solve(t(X * W) %*% X + Pen)
  V <- D %*% (t(X) %*% (X * r^2)) %*% D
  list(beta = beta, se = sqrt(diag(V)), V = V)
}
