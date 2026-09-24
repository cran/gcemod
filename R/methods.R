#' @export
print.gcemod <- function(x, ...) {
  cat(sprintf("Generalized competing event model [%s]\n", x$method))
  cat(sprintf("  n = %d;  events of interest = %d;  competing events = %d;  t = %g\n",
              x$n, x$n_event, x$n_competing, x$t))
  cat("\nomega+ ratios (event of interest vs. competing event):\n")
  print(as.data.frame(x$result))
  if (!is.null(x$omnibus) && !is.na(x$omnibus["statistic"])) {
    cat(sprintf(
      "\nLunn-McNeil omnibus test (all omega+ ratios = 1): chi-sq = %.3f, df = %d, p = %s\n",
      x$omnibus["statistic"], as.integer(x$omnibus["df"]),
      .fmt_p(x$omnibus["p.value"])))
  }
  invisible(x)
}

#' @export
summary.gcemod <- function(object, ...) {
  cat(sprintf("Generalized competing event model [%s]\n", object$method))
  cat("Call: ", paste(deparse(object$call), collapse = "\n"), "\n", sep = "")
  cat(sprintf("  n = %d;  events of interest = %d;  competing events = %d;  t = %g\n",
              object$n, object$n_event, object$n_competing, object$t))
  cat("\n--- omega+ ratio (event of interest vs. competing event) ---\n")
  print(as.data.frame(object$result))
  if (!is.null(object$omnibus) && !is.na(object$omnibus["statistic"])) {
    cat(sprintf(
      "\nLunn-McNeil omnibus test: chi-sq = %.3f, df = %d, p = %s\n",
      object$omnibus["statistic"], as.integer(object$omnibus["df"]),
      .fmt_p(object$omnibus["p.value"])))
  }
  invisible(list(result = object$result, omnibus = object$omnibus))
}

#' @export
print.gcecif <- function(x, ...) {
  cat(sprintf("Generalized competing event model [%s]\n", x$method))
  cat(sprintf("  n = %d;  events of interest = %d;  competing events = %d;  t = %g\n",
              x$n, x$n_event, x$n_competing, x$t))
  cat(sprintf("  F1(t) = %.4f;  F2(t) = %.4f;  baseline rho = F1/F2 = %.4f\n",
              x$F1, x$F2, x$rho0))
  cat("\nrho ratios (odds the realized event is the event of interest):\n")
  print(as.data.frame(x$result))
  if (!is.null(x$omnibus) && !is.na(x$omnibus["statistic"])) {
    cat(sprintf(
      "\nOmnibus Wald test (all rho ratios = 1): chi-sq = %.3f, df = %d, p = %s\n",
      x$omnibus["statistic"], as.integer(x$omnibus["df"]),
      .fmt_p(x$omnibus["p.value"])))
  }
  invisible(x)
}

#' @export
summary.gcecif <- function(object, ...) {
  cat(sprintf("Generalized competing event model [%s]\n", object$method))
  cat("Call: ", paste(deparse(object$call), collapse = "\n"), "\n", sep = "")
  cat(sprintf("  n = %d;  events of interest = %d;  competing events = %d;  t = %g\n",
              object$n, object$n_event, object$n_competing, object$t))
  cat(sprintf("  F1(t) = %.4f;  F2(t) = %.4f;  baseline rho = F1/F2 = %.4f\n",
              object$F1, object$F2, object$rho0))
  cat("\n--- rho ratio (cumulative incidence of event of interest vs. competing) ---\n")
  print(as.data.frame(object$result))
  if (!is.null(object$omnibus) && !is.na(object$omnibus["statistic"])) {
    cat(sprintf(
      "\nOmnibus Wald test (all rho ratios = 1): chi-sq = %.3f, df = %d, p = %s\n",
      object$omnibus["statistic"], as.integer(object$omnibus["df"]),
      .fmt_p(object$omnibus["p.value"])))
  }
  invisible(list(result = object$result, omnibus = object$omnibus))
}

#' Tidy a GCE model into a data frame of estimates
#'
#' With \code{effects = FALSE} (the default) the omega+ ratio estimates are
#' returned. With \code{effects = TRUE} a covariate-effects comparison table is
#' returned that juxtaposes, for each covariate, the hazard ratio for the primary
#' (event-of-interest) hazard, the competing-event hazard, the total (composite,
#' any-event) hazard, and the relative hazard ratio (the omega+ ratio) from the
#' GCE model, each with a confidence interval. This lets the investigator see
#' whether a covariate acts mainly on the event of interest, the competing event,
#' both, or on their balance.
#'
#' The individual- and total-event columns are estimated on the same hazard scale
#' as the fitted model: cause-specific Cox models for a \code{\link{gcecox}} fit,
#' and Fine-Gray subdistribution models for a \code{\link{gcefg}} fit. The total
#' (composite) event has no competing event, so its subdistribution and
#' cause-specific hazards coincide and a Cox model on the composite endpoint is
#' used in both cases.
#'
#' @param object a \code{"gcemod"} object from \code{\link{gcecox}},
#'   \code{\link{gcefg}}, or \code{\link{gcecif}}. For a \code{\link{gcecif}} fit
#'   the \code{relative} column is the cumulative-incidence ratio (rho ratio) and
#'   the primary/competing/total columns are cause-specific hazard ratios.
#' @param effects logical; if \code{TRUE}, return the primary/competing/total/
#'   relative effects comparison table instead of the omega+ ratio table.
#'   Requires an object fit with \pkg{gcemod} >= 0.3.0 (which stores the design
#'   matrix). Default \code{FALSE}.
#' @param ... unused.
#' @return If \code{effects = FALSE}, a data frame with columns \code{term},
#'   \code{log_ratio}, \code{ratio}, \code{conf.low}, \code{conf.high},
#'   \code{p.value}. If \code{effects = TRUE}, a data frame with columns
#'   \code{term}, \code{primary}, \code{competing}, \code{total}, and
#'   \code{relative}, each a formatted "HR (low, high)" string.
#' @examples
#' data(hn)
#' Ind <- data.frame(event = as.integer(hn$status == 1),
#'                   competing = as.integer(hn$status == 2))
#' fit <- gcecox(hn$time, Ind, hn[, c("age", "t_cat", "p16")], M = 5, t = 5)
#' gce_tidy(fit)
#' gce_tidy(fit, effects = TRUE)
#' @export
gce_tidy <- function(object, effects = FALSE, ...) {
  if (!inherits(object, "gcemod")) stop("'object' must be a gcemod object.", call. = FALSE)
  tab  <- object$result
  base <- data.frame(
    term      = rownames(tab),
    log_ratio = as.numeric(object$coef),
    ratio     = as.numeric(tab[, 1]),
    conf.low  = as.numeric(tab[, 2]),
    conf.high = as.numeric(tab[, 3]),
    p.value   = as.numeric(tab[, 4]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(effects)) return(base)

  X <- object$design
  if (is.null(X))
    stop("This object does not carry a design matrix; refit with gcemod >= 0.3.0 ",
         "to use 'effects = TRUE'.", call. = FALSE)
  time  <- object$time
  cause <- object$cause
  cn    <- object$covnames
  z     <- stats::qnorm(1 - (1 - object$conf.level) / 2)
  fmt   <- function(b, se) sprintf("%.2f (%.2f, %.2f)",
                                   exp(b), exp(b - z * se), exp(b + z * se))

  ## coefficients come back in the column order of X (== cn), so align by position
  cox_hr <- function(ev) {
    f <- survival::coxph(survival::Surv(time, ev) ~ X)
    fmt(stats::coef(f), sqrt(diag(stats::vcov(f))))
  }
  crr_hr <- function(fit) fmt(fit$coef, sqrt(diag(fit$var)))

  if (inherits(object, "gcefg")) {          # subdistribution scale
    primary   <- crr_hr(object$fit1)
    competing <- crr_hr(object$fit2)
  } else {                                   # cause-specific (Cox) scale
    primary   <- cox_hr(as.integer(cause == 1L))
    competing <- cox_hr(as.integer(cause == 2L))
  }
  total    <- cox_hr(as.integer(cause > 0L)) # composite endpoint: no competing event
  relative <- sprintf("%.2f (%.2f, %.2f)", base$ratio, base$conf.low, base$conf.high)

  data.frame(
    term      = cn,
    primary   = primary,
    competing = competing,
    total     = total,
    relative  = relative,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
