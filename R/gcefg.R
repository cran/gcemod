# Internal: subdistribution cumulative hazard from a crr prediction matrix.
# pc is the matrix returned by predict.crr (column 1 = time, remaining columns
# = CIF for each covariate row). Returns a length(times) x nsubject matrix of
# H(t) = -log(1 - F(t)), read directly at the requested times (no extrapolation).
.crr_H <- function(pc, times) {
  tt  <- pc[, 1]
  cif <- pc[, -1, drop = FALSE]
  idx <- findInterval(times, tt)
  idx[idx < 1L] <- 1L
  cift <- cif[idx, , drop = FALSE]
  cift <- pmin(pmax(cift, 0), 1 - 1e-12)
  -log(1 - cift)
}

#' Generalized competing event model via Fine-Gray regression
#'
#' @description
#' Fits a generalized competing event (GCE) model on the subdistribution hazard
#' scale and estimates covariate effects on \eqn{\omega^{+}} using a
#' Lunn-McNeil-style stacked Fine-Gray model (see \code{\link{lunnmcneil}} with
#' \code{type = "finegray"}), giving 95\% confidence intervals, Wald p-values,
#' and an omnibus test from one joint model. A GCE risk score is formed from the
#' model linear predictor, and per-subject \eqn{\omega^{+}} / \eqn{\omega} at
#' time \code{t} are returned for comparison with observed values.
#'
#' @param Time numeric vector of follow-up times.
#' @param Ind a data frame of event indicators. Column 1 = event of interest,
#'   column 2 = competing event (\code{0}/\code{1}).
#' @param Cov a data frame of covariates (factors are expanded to indicators).
#' @param M number of bins for the \eqn{\omega^{+}} risk plot (default \code{5}).
#' @param t time point at which per-subject \eqn{\omega^{+}} / \eqn{\omega} are
#'   evaluated.
#' @param conf.level confidence level (default \code{0.95}).
#' @param standardize logical; controls only the \emph{reported} omega+ ratio
#'   scale (default \code{FALSE} = per natural covariate unit; \code{TRUE} =
#'   per 1 SD). The risk score is always the normalized (mean-centered) linear
#'   predictor. Covariate means/SDs are in \code{$scaling}
#'   (see \code{\link{gce_scaling}}).
#'
#' @return An object of class \code{"gcemod"} (subclass \code{"gcefg"}) with the
#'   same components as \code{\link{gcecox}} on the subdistribution scale.
#'
#' @examples
#' data(hn)
#' Ind <- data.frame(event = as.integer(hn$status == 1),
#'                   competing = as.integer(hn$status == 2))
#' Cov <- hn[, c("age", "smoker", "t_cat", "n_cat", "p16")]
#' fit <- gcefg(hn$time, Ind, Cov, M = 5, t = 5)
#' summary(fit)
#' @seealso \code{\link{gcecox}}, \code{\link{lunnmcneil}}
#' @export
gcefg <- function(Time, Ind, Cov, M = 5, t, conf.level = 0.95, standardize = FALSE) {
  cl <- match.call()
  chk <- .gce_check(Time, Ind, Cov, M, t)
  Time <- chk$Time; Ind <- chk$Ind; Cov <- chk$Cov
  n <- length(Time)

  CA <- as.integer(Ind[[1]])
  CM <- as.integer(Ind[[2]])
  cause <- integer(n)
  cause[CM == 1] <- 2L
  cause[CA == 1] <- 1L

  ## ---- omega+ via Lunn-McNeil (subdistribution / Fine-Gray) ----
  lm.fit <- lunnmcneil(Time, cause, Cov, event = 1, competing = 2,
                       type = "finegray", conf.level = conf.level,
                       standardize = standardize)
  coefop <- lm.fit$omegaplus$log_ratio
  names(coefop) <- lm.fit$omegaplus$term
  cn <- names(coefop)

  ## ---- design matrix and normalized risk score (always mean-centered) ----
  Xd <- stats::model.matrix(~ ., data = Cov)   # '~ .' tolerates non-syntactic names
  Xd <- Xd[, colnames(Xd) != "(Intercept)", drop = FALSE]
  Xd <- Xd[, cn, drop = FALSE]
  ctr <- lm.fit$center[cn]; sds <- lm.fit$scale[cn]
  coef_raw <- lm.fit$coef_raw[cn]
  coef_std <- coef_raw * sds
  riskscore <- as.numeric(sweep(Xd, 2, ctr, "-") %*% coef_raw)
  normCER <- as.numeric(scale(riskscore))

  ## ---- per-subject subdistribution cumulative hazards at t (via crr) ----
  fit1 <- suppressWarnings(cmprsk::crr(Time, cause, Xd, failcode = 1, cencode = 0))
  fit2 <- suppressWarnings(cmprsk::crr(Time, cause, Xd, failcode = 2, cencode = 0))
  pc1 <- stats::predict(fit1, as.matrix(Xd))
  pc2 <- stats::predict(fit2, as.matrix(Xd))
  Hca_t <- as.numeric(.crr_H(pc1, t))
  Hcm_t <- as.numeric(.crr_H(pc2, t))
  omegaplus_sub <- Hca_t / Hcm_t
  omega_sub     <- Hca_t / (Hca_t + Hcm_t)

  ## ---- omega+ vs risk-score plot (M bins) ----
  qb <- seq(0, 1, length.out = M + 1)
  rs_q <- stats::quantile(normCER, probs = qb[-1], na.rm = TRUE)
  op_q <- stats::quantile(omegaplus_sub, probs = qb[-1], na.rm = TRUE)
  z1p <- .gce_scatter(rs_q, op_q, "Risk score", expression(omega^"+"))

  ## ---- omega vs time (population average, direct CIF) ----
  grid <- seq(t / 20, t, length.out = 20)
  Hca_g <- .crr_H(pc1, grid)
  Hcm_g <- .crr_H(pc2, grid)
  omega_t <- rowMeans(Hca_g) / (rowMeans(Hca_g) + rowMeans(Hcm_g))
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
    call = cl, method = "Fine-Gray (subdistribution)",
    coef = coefop, result = result,
    coef_raw = coef_raw, coef_std = coef_std,
    center = ctr, scale = sds, scaling = scaling,
    omnibus = lm.fit$omnibus, lunnmcneil = lm.fit,
    omegaplot = z1p, omegatimeplot = z3p,
    riskscore = riskscore, normCER = normCER,
    omega = omega_sub, omegaplus = omegaplus_sub,
    fit1 = fit1, fit2 = fit2,
    time = Time, cause = cause, design = Xd,
    covnames = cn, covariate.names = names(Cov),
    n = n, n_event = sum(cause == 1), n_competing = sum(cause == 2),
    t = t, conf.level = conf.level
  ), class = c("gcefg", "gcemod"))
}
