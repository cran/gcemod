# Internal helpers for gcemod. Not exported.

#' Validate and clean GCE inputs
#'
#' Coerces \code{Cov} to a data frame, checks argument types, and drops rows with
#' missing values (complete-case) with a warning. Returns cleaned objects.
#'
#' @param Time numeric vector of follow-up times.
#' @param Ind data frame of event indicators (structure checked by the caller).
#' @param Cov covariate data frame (or object coercible to one).
#' @param M,t positive scalars.
#' @return A list with cleaned \code{Time}, \code{Ind}, \code{Cov}, and the row
#'   index \code{keep} that was retained.
#' @keywords internal
#' @noRd
.gce_check <- function(Time, Ind, Cov, M, t) {
  if (missing(t) || is.null(t)) stop("'t' (evaluation time) must be supplied.", call. = FALSE)
  Cov <- as.data.frame(Cov)
  Ind <- as.data.frame(Ind)

  if (!is.numeric(Time)) stop("'Time' must be numeric.", call. = FALSE)
  if (nrow(Cov) != length(Time)) stop("'Cov' and 'Time' imply different sample sizes.", call. = FALSE)
  if (nrow(Ind) != length(Time)) stop("'Ind' and 'Time' imply different sample sizes.", call. = FALSE)
  if (ncol(Cov) < 1L) stop("'Cov' must contain at least one covariate.", call. = FALSE)
  for (nm in c("M", "t")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || v <= 0) {
      stop(sprintf("'%s' must be a single positive number.", nm), call. = FALSE)
    }
  }
  if (any(Time < 0, na.rm = TRUE)) stop("'Time' must be non-negative.", call. = FALSE)

  # Guard against invalid column names that would break formula building.
  if (is.null(names(Cov)) || any(names(Cov) == "" ) || anyDuplicated(names(Cov))) {
    stop("'Cov' must have unique, non-empty column names.", call. = FALSE)
  }

  keep <- stats::complete.cases(Time, Ind, Cov)
  if (!all(keep)) {
    warning(sprintf("Dropping %d row(s) with missing values (complete-case analysis).",
                    sum(!keep)), call. = FALSE)
  }
  Cov <- Cov[keep, , drop = FALSE]

  # Drop constant covariates up front (a constant column is collinear with the
  # intercept and makes the Cox / Fine-Gray models unidentifiable).
  const <- vapply(Cov, function(x) length(unique(x[!is.na(x)])) <= 1L, logical(1))
  if (any(const)) {
    warning(sprintf("Dropping constant covariate(s): %s.",
                    paste(names(Cov)[const], collapse = ", ")), call. = FALSE)
    Cov <- Cov[, !const, drop = FALSE]
  }
  if (ncol(Cov) < 1L) stop("No usable (non-constant) covariates remain.", call. = FALSE)

  list(Time = Time[keep],
       Ind  = Ind[keep, , drop = FALSE],
       Cov  = Cov,
       keep = keep)
}

#' Drop constant columns from a numeric design matrix
#' @keywords internal
#' @noRd
.gce_drop_constant <- function(X) {
  keep <- apply(X, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
  X[, keep, drop = FALSE]
}

#' Build a survival formula safely for any number of covariates
#'
#' @param response a character string giving the left-hand side, e.g.
#'   \code{"Surv(time, status)"}.
#' @param terms character vector of right-hand-side terms.
#' @return a \code{formula}.
#' @keywords internal
#' @noRd
.gce_formula <- function(response, terms) {
  # reformulate handles the single-term case that manual paste()+loop did not.
  # A string response such as "Surv(time, status)" is inserted verbatim on the LHS.
  stats::reformulate(termlabels = terms, response = response)
}

#' Step-function value of a cumulative hazard at requested times
#'
#' Reads the Nelson-Aalen cumulative hazard carried by a \code{survfit} object
#' and evaluates it at \code{times} using a right-continuous step function. This
#' replaces the \code{lm(-log(S) ~ time)} extrapolation used by \pkg{gcerisk}.
#'
#' @param sf a \code{survfit} object (single curve).
#' @param times numeric vector of evaluation times.
#' @return numeric vector of cumulative hazards, one per element of \code{times}.
#' @keywords internal
#' @noRd
.cumhaz_at <- function(sf, times) {
  ch <- sf$cumhaz
  tt <- sf$time
  if (is.null(ch)) {
    # Fall back to -log(S) if cumhaz is unavailable.
    ch <- -log(sf$surv)
  }
  ok <- is.finite(ch)
  tt <- tt[ok]; ch <- ch[ok]
  if (length(tt) == 0L) return(rep(NA_real_, length(times)))
  f <- stats::stepfun(tt, c(0, ch))
  f(times)
}

#' Aalen-Johansen cumulative incidence of both causes at a landmark time
#'
#' Returns \code{c(F1, F2)}, the cumulative incidence of the event of interest
#' (cause 1) and the competing event (cause 2) at time \code{t}, from the
#' Aalen-Johansen estimator (product-limit all-cause survivor times the
#' cause-specific increments). Used for the observed cumulative-incidence ratio
#' in \code{\link{gce_calibration}} of a \code{\link{gcecif}} fit.
#'
#' @param time numeric follow-up times.
#' @param cause integer cause coding (0 censored, 1 event of interest, 2 competing).
#' @param t landmark time.
#' @return numeric length-2 vector \code{c(F1, F2)}.
#' @keywords internal
#' @noRd
.aj_cif_at <- function(time, cause, t) {
  ut <- sort(unique(time[cause > 0 & time <= t]))
  if (length(ut) < 1L) return(c(NA_real_, NA_real_))
  Y  <- vapply(ut, function(u) sum(time >= u), numeric(1))
  d1 <- vapply(ut, function(u) sum(time == u & cause == 1L), numeric(1))
  d2 <- vapply(ut, function(u) sum(time == u & cause == 2L), numeric(1))
  haz   <- ifelse(Y > 0, (d1 + d2) / Y, 0)
  Sprev <- c(1, cumprod(1 - haz)[-length(haz)])
  c(sum(Sprev * ifelse(Y > 0, d1 / Y, 0)),
    sum(Sprev * ifelse(Y > 0, d2 / Y, 0)))
}

#' Format a numeric p-value for printing
#' @keywords internal
#' @noRd
.fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, "<1e-04", formatC(p, digits = 4, format = "f")))
}
