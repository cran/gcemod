#' Generalized competing event model via cause-specific Cox regression
#'
#' @description
#' Fits a generalized competing event (GCE) model on the cause-specific hazard
#' scale and estimates covariate effects on \eqn{\omega^{+}} (event of interest
#' vs. competing event) using the Lunn-McNeil stacked model (see
#' \code{\link{lunnmcneil}}), with 95\% confidence intervals, Wald p-values, and
#' an omnibus test. A GCE risk score is formed from the model linear predictor,
#' and per-subject \eqn{\omega = \omega^{+}/(1+\omega^{+})} is returned for
#' comparison with observed values.
#'
#' @param Time numeric vector of follow-up times.
#' @param Ind a data frame of event indicators. Column 1 = event of interest,
#'   column 2 = competing event (\code{0}/\code{1}). A third column, if present,
#'   is ignored.
#' @param Cov a data frame of covariates.
#' @param M number of bins for the \eqn{\omega^{+}} risk plot (default \code{5}).
#' @param t time point at which per-subject \eqn{\omega} / \eqn{\omega^{+}} are
#'   evaluated.
#' @param conf.level confidence level (default \code{0.95}).
#' @param cluster.se logical; cluster-robust SEs in the Lunn-McNeil model
#'   (default \code{TRUE}).
#' @param standardize logical; controls only the \emph{reported} omega+ ratio
#'   scale (default \code{FALSE} = per natural covariate unit; \code{TRUE} =
#'   per 1 SD). The risk score is \emph{always} the normalized (mean-centered)
#'   linear predictor regardless of this setting, so cutpoints, plots, and the
#'   predicted-vs-observed omega calibration are unaffected. Covariate means/SDs
#'   are available in \code{$scaling} (see \code{\link{gce_scaling}}).
#'
#' @return An object of class \code{"gcemod"} (subclass \code{"gcecox"}) with:
#'   \code{coef} (log \eqn{\omega^{+}} ratios), \code{result} (ratio table with
#'   CIs and p-values), \code{omnibus} (Lunn-McNeil global test),
#'   \code{riskscore} (linear predictor), \code{omega} / \code{omegaplus}
#'   (per-subject values at \code{t}), \code{omegaplot} and \code{omegatimeplot}
#'   (ggplot objects), and the fitted Lunn-McNeil model.
#'
#' @examples
#' data(hn)
#' Ind <- data.frame(event = as.integer(hn$status == 1),
#'                   competing = as.integer(hn$status == 2))
#' Cov <- hn[, c("age", "smoker", "t_cat", "n_cat", "p16")]
#' fit <- gcecox(hn$time, Ind, Cov, M = 5, t = 5)
#' summary(fit)
#' @seealso \code{\link{gcefg}}, \code{\link{lunnmcneil}}, \code{\link{gce_cutpoints}}
#' @export
gcecox <- function(Time, Ind, Cov, M = 5, t,
                   conf.level = 0.95, cluster.se = TRUE, standardize = FALSE) {
  cl <- match.call()
  chk <- .gce_check(Time, Ind, Cov, M, t)
  Time <- chk$Time; Ind <- chk$Ind; Cov <- chk$Cov
  n <- length(Time)

  CA <- as.integer(Ind[[1]])
  CM <- as.integer(Ind[[2]])

  # Mutually exclusive cause coding (event of interest wins ties).
  cause <- integer(n)
  cause[CM == 1] <- 2L
  cause[CA == 1] <- 1L

  ## ---- omega+ via Lunn-McNeil (cause-specific) ----
  lm.fit <- lunnmcneil(Time, cause, Cov, event = 1, competing = 2,
                       type = "coxph", cluster.se = cluster.se,
                       conf.level = conf.level, standardize = standardize)
  coefop <- lm.fit$omegaplus$log_ratio
  names(coefop) <- lm.fit$omegaplus$term
  cn <- names(coefop)

  ## ---- risk score: always the normalized (mean-centered) linear predictor ----
  ## Uses the RAW coefficients on mean-centered covariates, equivalently the
  ## per-SD coefficients on standardized covariates. Independent of the reported
  ## coefficient scale, and centered so the average patient has riskscore = 0.
  Xd <- stats::model.matrix(~ ., data = Cov)   # '~ .' tolerates non-syntactic names
  Xd <- Xd[, colnames(Xd) != "(Intercept)", drop = FALSE]
  Xd <- Xd[, cn, drop = FALSE]
  ctr <- lm.fit$center[cn]; sds <- lm.fit$scale[cn]
  coef_raw <- lm.fit$coef_raw[cn]
  coef_std <- coef_raw * sds
  riskscore <- as.numeric(sweep(Xd, 2, ctr, "-") %*% coef_raw)
  normCER <- as.numeric(scale(riskscore))

  ## ---- baseline cause-specific cumulative hazards at t ----
  md <- data.frame(Time = Time, CA = CA, CM = CM)
  sf.ca <- survfit(Surv(Time, CA) ~ 1, data = md)
  sf.cm <- survfit(Surv(Time, CM) ~ 1, data = md)
  Hca_t <- .cumhaz_at(sf.ca, t)
  Hcm_t <- .cumhaz_at(sf.cm, t)
  w0plus <- Hca_t / Hcm_t

  omegaplus_sub <- w0plus * exp(riskscore)        # per-subject omega+
  omega_sub     <- omegaplus_sub / (1 + omegaplus_sub)  # per-subject omega

  ## ---- omega+ vs risk-score plot (M bins) ----
  qb <- seq(0, 1, length.out = M + 1)
  rs_q <- stats::quantile(normCER, probs = qb[-1], na.rm = TRUE)
  op_q <- stats::quantile(omegaplus_sub, probs = qb[-1], na.rm = TRUE)
  z1p <- .gce_scatter(rs_q, op_q, "Risk score", expression(omega^"+"))

  ## ---- omega vs time plot (population, direct cumulative hazards) ----
  grid <- seq(t / 20, t, length.out = 20)
  Hca_g <- .cumhaz_at(sf.ca, grid)
  Hcm_g <- .cumhaz_at(sf.cm, grid)
  omega_t <- Hca_g / (Hca_g + Hcm_g)
  z3p <- .gce_scatter(grid, omega_t, "Time", expression(omega),
                      ylim = c(0, min(1, 2 * max(omega_t, na.rm = TRUE))))

  ## ---- result table ----
  op <- lm.fit$omegaplus
  result <- cbind(`exp(coef) (omega+ ratio)` = round(op$omegaplus_ratio, 5),
                  `lower .95` = round(op$conf.low, 5),
                  `upper .95` = round(op$conf.high, 5),
                  `P-value`   = round(op$p.value, 5))
  rownames(result) <- op$term

  scaling <- data.frame(term = cn, mean = as.numeric(ctr), sd = as.numeric(sds),
                        coef_perSD = as.numeric(coef_std),
                        omegaplus_ratio_perSD = exp(as.numeric(coef_std)),
                        row.names = NULL, stringsAsFactors = FALSE)

  structure(list(
    call = cl, method = "Cox (cause-specific)",
    coef = coefop, result = result,
    coef_raw = coef_raw, coef_std = coef_std,
    center = ctr, scale = sds, scaling = scaling, w0plus = w0plus,
    omnibus = lm.fit$omnibus, lunnmcneil = lm.fit,
    omegaplot = z1p, omegatimeplot = z3p,
    riskscore = riskscore, normCER = normCER,
    omega = omega_sub, omegaplus = omegaplus_sub,
    time = Time, cause = cause, design = Xd,
    covnames = cn, covariate.names = names(Cov),
    n = n, n_event = sum(cause == 1), n_competing = sum(cause == 2),
    t = t, conf.level = conf.level
  ), class = c("gcecox", "gcemod"))
}
