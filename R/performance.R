#' Discrimination (C-index) for a GCE model
#'
#' @description
#' Computes the cause-specific concordance (C-index) of the GCE risk score for
#' the event of interest. For calibration (predicted vs. observed
#' \eqn{\omega^{+}}), see \code{\link{gce_calibration}}.
#'
#' @param object a \code{"gcemod"} object.
#' @param conf.level confidence level for the C-index interval (default
#'   \code{0.95}).
#'
#' @details
#' The C-index is computed with \code{\link[survival]{concordance}} using the
#' event of interest as the outcome and the GCE risk score as the predictor.
#' Higher risk scores denote higher risk; the statistic is oriented so that
#' values above 0.5 indicate discrimination in the expected direction, and the
#' \code{"flipped"} attribute of the result records whether the orientation was
#' reversed.
#'
#' @return An object of class \code{"gce_performance"}: a named numeric vector
#'   with elements \code{estimate}, \code{se}, \code{lower}, and \code{upper},
#'   carrying attributes \code{conf.level} and \code{flipped}. Printing the
#'   object displays the C-index with its confidence interval.
#'
#' @examples
#' data(hn)
#' Ind <- data.frame(event = as.integer(hn$status == 1),
#'                   competing = as.integer(hn$status == 2))
#' fit <- gcecox(hn$time, Ind, hn[, c("age", "t_cat", "p16")], M = 5, t = 5)
#' gce_performance(fit)
#' @seealso \code{\link{gce_calibration}}
#' @export
gce_performance <- function(object, conf.level = 0.95) {
  if (!inherits(object, "gcemod")) stop("'object' must be a gcemod object.", call. = FALSE)

  time <- object$time
  ev   <- as.integer(object$cause == 1)
  rs   <- object$riskscore

  cc <- survival::concordance(Surv(time, ev) ~ rs)
  Cval <- unname(cc$concordance)
  v    <- unname(cc$var)
  flipped <- FALSE
  if (is.finite(Cval) && Cval < 0.5) { Cval <- 1 - Cval; flipped <- TRUE }
  zc <- stats::qnorm(1 - (1 - conf.level) / 2)
  structure(c(estimate = Cval, se = sqrt(v),
              lower = Cval - zc * sqrt(v), upper = Cval + zc * sqrt(v)),
            conf.level = conf.level, flipped = flipped,
            class = "gce_performance")
}

#' @export
print.gce_performance <- function(x, digits = 3, ...) {
  cl <- attr(x, "conf.level"); if (is.null(cl)) cl <- 0.95
  cat(sprintf("C-index (event of interest): %.*f (%g%% CI %.*f-%.*f)%s\n",
              digits, x["estimate"], 100 * cl, digits, x["lower"],
              digits, x["upper"],
              if (isTRUE(attr(x, "flipped")))
                "  [risk-score orientation flipped]" else ""))
  invisible(x)
}
