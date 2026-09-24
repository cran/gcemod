#' Lunn-McNeil stacked-model estimation of omega+ ratios
#'
#' @description
#' Estimates covariate effects on \eqn{\omega^{+}} (the ratio of the hazard for
#' the event of interest to the hazard for the competing event) from a single
#' stacked ("augmented") model, following Lunn and McNeil (1995). Each subject
#' contributes one record per event type; the model is stratified by event type
#' and the covariate-by-event-type interaction estimates
#' \eqn{\log(\omega^{+}\text{ ratio})} with model-based standard errors,
#' per-covariate Wald tests, and an omnibus multivariate Wald test.
#'
#' @details
#' Two hazard scales are supported:
#' \describe{
#'   \item{\code{type = "coxph"}}{Cause-specific hazards. The augmented data set
#'     duplicates each subject into an "event" stratum and a "competing" stratum;
#'     the interaction coefficient is exactly the difference in cause-specific
#'     log-hazards, i.e. \eqn{\log(\omega^{+})}. This is the original Lunn-McNeil
#'     construction and the interaction MLE equals the difference of two separate
#'     cause-specific Cox fits.}
#'   \item{\code{type = "finegray"}}{Subdistribution hazards. Risk-set-weighted
#'     data are generated for each event type with
#'     \code{\link[survival]{finegray}}, stacked, and fit with a single weighted,
#'     event-type-stratified \code{\link[survival]{coxph}}. The interaction
#'     estimates the difference in subdistribution log-hazards. Because subjects
#'     recur across the weighted risk sets, cluster-robust (subject) standard
#'     errors are used. This is an extension of Lunn-McNeil to the
#'     subdistribution scale.}
#' }
#' Subject records are duplicated in both cases, so cluster-robust standard
#' errors (clustered on subject) are the default.
#'
#' @param Time numeric vector of follow-up times.
#' @param cause event indicator per subject: \code{0} = censored, and the codes
#'   in \code{event} and \code{competing} for the two event types.
#' @param Cov a data frame of covariates (factors are expanded to indicators).
#' @param event code marking the event of interest (default \code{1}).
#' @param competing code marking the competing event (default \code{2}).
#' @param type \code{"coxph"} (cause-specific, default) or \code{"finegray"}
#'   (subdistribution).
#' @param cluster.se logical; cluster-robust (subject) SEs (default \code{TRUE}).
#' @param conf.level confidence level (default \code{0.95}).
#' @param standardize logical; controls only the \emph{reported} coefficient
#'   scale (default \code{FALSE}). \code{FALSE} reports omega+ ratios per natural
#'   covariate unit; \code{TRUE} reports them per 1 SD. This is a
#'   reparameterization, so Wald tests, p-values, and the omnibus test are
#'   identical either way. The covariate means/SDs (\code{$center}, \code{$scale})
#'   are always returned.
#' @param ... passed to \code{\link[survival]{coxph}}.
#'
#' @return An object of class \code{"gce_lunnmcneil"}: a list with the fitted
#'   model (\code{$fit}), the omega+ coefficient table (\code{$omegaplus}), the
#'   reference (competing-event) coefficients (\code{$beta_competing}), the
#'   omnibus test (\code{$omnibus}), the design column names (\code{$terms}),
#'   and the hazard scale (\code{$type}).
#'
#' @references Lunn M, McNeil D (1995) Applying Cox regression to competing
#'   risks. \emph{Biometrics} 51:524-32.
#'
#' @examples
#' data(hn)
#' Cov <- hn[, c("age", "smoker", "t_cat", "n_cat", "p16")]
#' lunnmcneil(hn$time, hn$status, Cov, type = "coxph")
#' lunnmcneil(hn$time, hn$status, Cov, type = "finegray")
#' @export
lunnmcneil <- function(Time, cause, Cov, event = 1, competing = 2,
                       type = c("coxph", "finegray"),
                       cluster.se = TRUE, conf.level = 0.95,
                       standardize = FALSE, ...) {
  type <- match.arg(type)
  Cov <- as.data.frame(Cov)
  if (length(Time) != nrow(Cov) || length(cause) != nrow(Cov)) {
    stop("'Time', 'cause', and 'Cov' must have the same number of observations.",
         call. = FALSE)
  }
  keep <- stats::complete.cases(Time, cause, Cov)
  if (!all(keep)) {
    warning(sprintf("Dropping %d row(s) with missing values.", sum(!keep)),
            call. = FALSE)
    Time <- Time[keep]; cause <- cause[keep]; Cov <- Cov[keep, , drop = FALSE]
  }
  bad <- !(cause %in% c(0, event, competing))
  if (any(bad)) {
    stop("'cause' contains codes other than 0, 'event', and 'competing'. ",
         "Recode multi-cause data to a single competing category first.",
         call. = FALSE)
  }

  # Numeric design matrix (factors expanded); drop intercept and constants.
  # '~ .' (rather than reformulate(names(Cov))) so arbitrary, non-syntactic
  # column names (spaces, punctuation, leading digits) are handled.
  X <- stats::model.matrix(~ ., data = Cov)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  X <- .gce_drop_constant(X)
  cn <- colnames(X)
  p <- ncol(X)
  n <- nrow(X)
  if (p < 1L) stop("No usable (non-constant) covariates remain.", call. = FALSE)

  # Covariate means/SDs. Always computed and returned: they define the
  # normalized risk score (centered at the mean covariate profile) and let new
  # subjects be scored. The model is fit on the RAW covariate scale; the per-SD
  # report is obtained by a reparameterization (theta_sd = theta_raw * SD),
  # which leaves the Wald tests, p-values, and omnibus test unchanged.
  center <- colMeans(X)
  scale  <- apply(X, 2, stats::sd)
  names(center) <- names(scale) <- cn

  if (type == "coxph") {
    fit_info <- .lm_coxph(Time, cause, X, cn, event, competing, cluster.se, ...)
  } else {
    fit_info <- .lm_finegray(Time, cause, X, cn, event, competing, ...)
  }
  fit <- fit_info$fit
  idx <- fit_info$idx        # positions of the interaction (omega+) coefficients

  cf <- coef(fit)
  V  <- stats::vcov(fit)
  theta_raw <- cf[idx]
  se_raw    <- sqrt(diag(V)[idx])
  names(theta_raw) <- names(se_raw) <- cn

  # Reported scale: raw per-unit (default) or per-SD (standardize = TRUE).
  if (isTRUE(standardize)) {
    theta <- theta_raw * scale
    se    <- se_raw * scale
  } else {
    theta <- theta_raw
    se    <- se_raw
  }

  zcrit <- stats::qnorm(1 - (1 - conf.level) / 2)
  z <- theta_raw / se_raw                 # invariant to the reporting scale
  pval <- 2 * stats::pnorm(-abs(z))

  omegaplus <- data.frame(
    term            = cn,
    log_ratio       = as.numeric(theta),
    se              = as.numeric(se),
    omegaplus_ratio = exp(as.numeric(theta)),
    conf.low        = exp(theta - zcrit * se),
    conf.high       = exp(theta + zcrit * se),
    z               = as.numeric(z),
    p.value         = as.numeric(pval),
    row.names = NULL, stringsAsFactors = FALSE
  )

  Vsub <- V[idx, idx, drop = FALSE]
  omni <- tryCatch({
    stat <- as.numeric(t(theta_raw) %*% solve(Vsub, theta_raw))
    c(statistic = stat, df = p,
      p.value = stats::pchisq(stat, df = p, lower.tail = FALSE))
  }, error = function(e) {
    warning("Omnibus test could not be computed (singular covariance).",
            call. = FALSE)
    c(statistic = NA_real_, df = p, p.value = NA_real_)
  })

  beta_competing <- cf[fit_info$main_idx]
  names(beta_competing) <- cn

  structure(list(
    call = match.call(), type = type, fit = fit,
    omegaplus = omegaplus, coef_raw = theta_raw, beta_competing = beta_competing,
    omnibus = omni, terms = cn, conf.level = conf.level,
    standardize = isTRUE(standardize), center = center, scale = scale,
    n = n, n_event = sum(cause == event), n_competing = sum(cause == competing)
  ), class = "gce_lunnmcneil")
}

# ---- cause-specific (Cox) Lunn-McNeil augmentation ----
# The model is fit on synthetic, syntactic column names (x1, x2, ...) so the
# formula is robust to arbitrary user covariate names (spaces, punctuation,
# etc.). Coefficients are extracted by position in lunnmcneil() and relabelled
# with the original names, so the synthetic names never surface to the user.
.lm_coxph <- function(Time, cause, X, cn, event, competing, cluster.se, ...) {
  p <- ncol(X); n <- nrow(X)
  vn  <- paste0("x", seq_len(p))            # main-effect names
  vni <- paste0("x", seq_len(p), "_ev")     # interaction names
  Xstack <- rbind(X, X); colnames(Xstack) <- vn
  d1     <- rep(c(1, 0), each = n)          # event-stratum indicator
  Xint   <- Xstack * d1; colnames(Xint) <- vni

  base_df <- data.frame(
    .time   = rep(Time, 2),
    .status = c(as.integer(cause == event), as.integer(cause == competing)),
    .etype  = factor(rep(c("event", "competing"), each = n),
                     levels = c("competing", "event")),
    .id     = rep(seq_len(n), 2),
    check.names = FALSE
  )
  aug <- cbind(base_df, Xstack, Xint)

  termlabels <- c(vn, vni, "strata(.etype)")
  f <- stats::reformulate(termlabels, response = "Surv(.time, .status)")
  fit <- if (isTRUE(cluster.se)) {
    coxph(f, data = aug, cluster = aug$.id, robust = TRUE, ...)
  } else {
    coxph(f, data = aug, ...)
  }
  list(fit = fit, main_idx = seq_len(p), idx = (p + 1L):(2L * p))
}

# ---- subdistribution (Fine-Gray) Lunn-McNeil-style augmentation ----
.lm_finegray <- function(Time, cause, X, cn, event, competing, ...) {
  p <- ncol(X); n <- nrow(X)
  vn  <- paste0("x", seq_len(p))            # synthetic, syntactic names
  vni <- paste0("x", seq_len(p), "_ev")
  colnames(X) <- vn
  fstat <- factor(cause, levels = c(0, event, competing),
                  labels = c("censor", "event", "competing"))
  df <- data.frame(.time = Time, .fstat = fstat, .id = seq_len(n),
                   X, check.names = FALSE)
  # Carry .id through the RHS so finegray() retains it for clustering (it keeps
  # only variables named in the formula).
  covform <- stats::reformulate(c(vn, ".id"),
                                response = "Surv(.time, .fstat)")

  fg_e <- survival::finegray(covform, data = df, etype = "event")
  fg_c <- survival::finegray(covform, data = df, etype = "competing")

  Xe <- as.matrix(fg_e[, vn, drop = FALSE])
  Xc <- as.matrix(fg_c[, vn, drop = FALSE])
  Xstack <- rbind(Xe, Xc); colnames(Xstack) <- vn
  d1     <- c(rep(1, nrow(fg_e)), rep(0, nrow(fg_c)))
  Xint   <- Xstack * d1; colnames(Xint) <- vni

  base_df <- data.frame(
    .start  = c(fg_e$fgstart, fg_c$fgstart),
    .stop   = c(fg_e$fgstop,  fg_c$fgstop),
    .status = c(fg_e$fgstatus, fg_c$fgstatus),
    .wt     = c(fg_e$fgwt,    fg_c$fgwt),
    .etype  = factor(c(rep("event", nrow(fg_e)), rep("competing", nrow(fg_c))),
                     levels = c("competing", "event")),
    .id     = c(fg_e$.id, fg_c$.id),
    check.names = FALSE
  )
  aug <- cbind(base_df, Xstack, Xint)

  termlabels <- c(vn, vni, "strata(.etype)")
  f <- stats::reformulate(termlabels, response = "Surv(.start, .stop, .status)")
  fit <- coxph(f, data = aug, weights = aug$.wt, cluster = aug$.id,
               robust = TRUE, ...)
  list(fit = fit, main_idx = seq_len(p), idx = (p + 1L):(2L * p))
}

#' @export
print.gce_lunnmcneil <- function(x, ...) {
  scl <- if (x$type == "finegray") "subdistribution (Fine-Gray)" else "cause-specific (Cox)"
  cat(sprintf("Lunn-McNeil stacked model [%s]\n", scl))
  cat(sprintf("  n = %d;  events of interest = %d;  competing events = %d\n",
              x$n, x$n_event, x$n_competing))
  cat("\nomega+ ratios (event of interest vs. competing event):\n")
  tab <- x$omegaplus
  out <- data.frame(
    `omega+ ratio` = round(tab$omegaplus_ratio, 4),
    `lower .95`    = round(tab$conf.low, 4),
    `upper .95`    = round(tab$conf.high, 4),
    `p-value`      = .fmt_p(tab$p.value),
    check.names = FALSE, row.names = tab$term)
  print(out)
  cat(sprintf("\nOmnibus Wald test (all omega+ ratios = 1): chi-sq = %.3f, df = %d, p = %s\n",
              x$omnibus["statistic"], as.integer(x$omnibus["df"]),
              .fmt_p(x$omnibus["p.value"])))
  invisible(x)
}
