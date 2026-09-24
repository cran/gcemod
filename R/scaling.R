#' Covariate scaling (means and SDs) used by the GCE risk score
#'
#' @description
#' Returns each covariate's mean and standard deviation (from the data used to
#' fit the model), together with the per-SD log-omega+ coefficient and omega+
#' ratio. These are the quantities needed to reproduce the normalized risk score
#' for new subjects (see \code{\link{gce_riskscore}}) and to interpret how far a
#' subject sits from the average patient.
#'
#' @param object a \code{"gcemod"} object from \code{\link{gcecox}},
#'   \code{\link{gcefg}}, or \code{\link{gcecif}}.
#' @return A data frame with columns \code{term}, \code{mean}, \code{sd},
#'   \code{coef_perSD} (log ratio per 1 SD) and the per-SD ratio
#'   (\code{omegaplus_ratio_perSD} for omega+ models,
#'   \code{rho_ratio_perSD} for a \code{\link{gcecif}} fit).
#' @examples
#' data(prostate)
#' Ind <- data.frame(as.integer(prostate$status == 1), as.integer(prostate$status == 2))
#' fit <- gcecox(prostate$time, Ind, prostate[, c("age", "psa", "gleason")], M = 5, t = 10)
#' gce_scaling(fit)
#' @seealso \code{\link{gce_riskscore}}
#' @export
gce_scaling <- function(object) {
  if (!inherits(object, "gcemod")) stop("'object' must be a gcemod object.", call. = FALSE)
  object$scaling
}

#' Normalized risk score and relative omega+ for new subjects
#'
#' @description
#' Computes the GCE risk score for new subjects using the covariate means/SDs
#' and coefficients stored in a fitted model. The score is the normalized
#' (mean-centered) linear predictor, so \code{exp(riskscore)} is the subject's
#' omega+ relative to the \emph{average} subject in the fitting data (a value of
#' 1 = same as average, 2 = twice the average omega+). For cause-specific
#' (\code{\link{gcecox}}) models the absolute predicted omega+ and omega at the
#' model's time \code{t} are also returned.
#'
#' @param object a \code{"gcemod"} object from \code{\link{gcecox}},
#'   \code{\link{gcefg}}, or \code{\link{gcecif}}.
#' @param newdata a data frame containing the covariates used to fit
#'   \code{object} (same column names). Factor levels not seen in the fitting
#'   data are treated as the reference level.
#' @return A data frame with \code{riskscore} and \code{exp(riskscore)} relative
#'   to the average subject: \code{omegaplus_rel} for omega+ models (with
#'   absolute \code{omegaplus} and \code{omega} at time \code{t} for
#'   cause-specific fits), or \code{rho_rel} for a \code{\link{gcecif}} fit.
#' @examples
#' data(prostate)
#' Ind <- data.frame(as.integer(prostate$status == 1), as.integer(prostate$status == 2))
#' fit <- gcecox(prostate$time, Ind, prostate[, c("age", "psa", "gleason")], M = 5, t = 10)
#' gce_riskscore(fit, newdata = prostate[1:5, ])
#' @seealso \code{\link{gce_scaling}}
#' @export
gce_riskscore <- function(object, newdata) {
  if (!inherits(object, "gcemod")) stop("'object' must be a gcemod object.", call. = FALSE)
  newdata <- as.data.frame(newdata)
  covn <- object$covariate.names
  miss <- setdiff(covn, names(newdata))
  if (length(miss)) stop(sprintf("'newdata' is missing covariate(s): %s.",
                                 paste(miss, collapse = ", ")), call. = FALSE)

  # Build the design matrix on newdata and align to the fitted design columns.
  # Subset to the covariates first, then '~ .', so non-syntactic column names
  # (and extra columns in newdata) are handled.
  Xf <- stats::model.matrix(~ ., data = newdata[, covn, drop = FALSE])
  Xf <- Xf[, colnames(Xf) != "(Intercept)", drop = FALSE]
  cn <- object$covnames
  Xnew <- matrix(0, nrow = nrow(Xf), ncol = length(cn), dimnames = list(NULL, cn))
  common <- intersect(colnames(Xf), cn)
  if (length(common)) Xnew[, common] <- Xf[, common]

  riskscore <- as.numeric(sweep(Xnew, 2, object$center[cn], "-") %*% object$coef_raw[cn])
  out <- data.frame(riskscore = riskscore,
                    omegaplus_rel = exp(riskscore),
                    stringsAsFactors = FALSE)
  if (!is.null(object$w0plus)) {
    op <- object$w0plus * exp(riskscore)
    out$omegaplus <- op
    out$omega <- op / (1 + op)
  }
  ## gcecif() works on the cumulative-incidence ratio: relabel exp(riskscore) as
  ## rho relative to the average subject (no absolute omega+ at time t).
  if (inherits(object, "gcecif"))
    names(out)[names(out) == "omegaplus_rel"] <- "rho_rel"
  out
}
