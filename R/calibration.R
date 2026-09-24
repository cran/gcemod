#' Calibration of predicted vs. observed omega+ by risk rank
#'
#' @description
#' Groups subjects into equal-sized rank groups on the GCE risk score and
#' compares the model-predicted quantity with an observed estimate in each group.
#' For omega+ models (\code{\link{gcecox}}, \code{\link{gcefg}}) the predicted
#' \eqn{\omega^{+}} (or \eqn{\omega}) is compared with observed cause-specific
#' Nelson-Aalen cumulative hazards at time \code{t}:
#' \eqn{\omega^{+}_{obs} = H_1(t)/H_2(t)} and
#' \eqn{\omega_{obs} = H_1(t)/(H_1(t)+H_2(t))}. For a \code{\link{gcecif}} fit the
#' predicted cumulative-incidence ratio \eqn{\rho} (or share \eqn{\pi_1}) is
#' compared with the observed Aalen-Johansen cumulative incidence within the
#' group: \eqn{\rho_{obs} = F_1(t)/F_2(t)} and
#' \eqn{\pi_{1,obs} = F_1(t)/(F_1(t)+F_2(t))}. A calibration plot (predicted vs.
#' observed with the identity line) and the underlying table are returned.
#'
#' @param object a \code{"gcemod"} object from \code{\link{gcecox}},
#'   \code{\link{gcefg}}, or \code{\link{gcecif}}.
#' @param which quantity to calibrate. For omega+ models \code{"omegaplus"}
#'   (default) or \code{"omega"}; for a \code{\link{gcecif}} fit \code{"rho"}
#'   (default) or \code{"share"}. If \code{NULL} the model-appropriate default is
#'   used.
#' @param groups number of rank groups (default \code{5}).
#'
#' @return A \code{ggplot} object (predicted vs. observed). The calibration
#'   table is attached as attribute \code{"data"}.
#'
#' @examples
#' data(hn)
#' Ind <- data.frame(event = as.integer(hn$status == 1),
#'                   competing = as.integer(hn$status == 2))
#' fit <- gcecox(hn$time, Ind, hn[, c("age", "t_cat", "p16")], M = 5, t = 5)
#' gce_calibration(fit, which = "omegaplus", groups = 5)
#' @export
gce_calibration <- function(object, which = NULL, groups = 5) {
  if (!inherits(object, "gcemod")) stop("'object' must be a gcemod object.", call. = FALSE)
  if (groups < 2L) stop("'groups' must be >= 2.", call. = FALSE)
  is_cif  <- inherits(object, "gcecif")
  choices <- if (is_cif) c("rho", "share") else c("omegaplus", "omega")
  which   <- if (is.null(which)) choices[1] else match.arg(which, choices)

  score <- object$riskscore
  time  <- object$time
  cause <- object$cause
  t     <- object$t
  pred_all <- if (is_cif) {
    if (which == "rho") object$rho else object$rho / (1 + object$rho)
  } else {
    if (which == "omegaplus") object$omegaplus else object$omega
  }

  br <- stats::quantile(score, probs = seq(0, 1, length.out = groups + 1),
                        na.rm = TRUE)
  br[1] <- -Inf; br[length(br)] <- Inf
  grp <- cut(score, breaks = br, labels = FALSE, include.lowest = TRUE)

  tab <- do.call(rbind, lapply(sort(unique(grp)), function(g) {
    sel <- grp == g
    if (is_cif) {
      f <- .aj_cif_at(time[sel], cause[sel], t)          # Aalen-Johansen CIFs
      obs <- if (which == "rho") f[1] / f[2] else f[1] / (f[1] + f[2])
    } else {
      sfa <- survfit(Surv(time[sel], as.integer(cause[sel] == 1)) ~ 1)
      sfc <- survfit(Surv(time[sel], as.integer(cause[sel] == 2)) ~ 1)
      H1 <- .cumhaz_at(sfa, t); H2 <- .cumhaz_at(sfc, t)
      obs <- if (which == "omegaplus") H1 / H2 else H1 / (H1 + H2)
    }
    data.frame(group = g, n = sum(sel),
               predicted = mean(pred_all[sel], na.rm = TRUE),
               observed = as.numeric(obs), stringsAsFactors = FALSE)
  }))
  rownames(tab) <- NULL

  ylab <- switch(which,
                 omegaplus = expression(observed~omega^"+"),
                 omega     = expression(observed~omega),
                 rho       = expression(observed~rho),
                 share     = expression(observed~pi[1]))
  xlab <- switch(which,
                 omegaplus = expression(predicted~omega^"+"),
                 omega     = expression(predicted~omega),
                 rho       = expression(predicted~rho),
                 share     = expression(predicted~pi[1]))
  rng <- range(c(tab$predicted, tab$observed), na.rm = TRUE, finite = TRUE)
  p <- ggplot2::ggplot(
        tab, ggplot2::aes(x = .data[["predicted"]], y = .data[["observed"]])) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(alpha = 0.4) +
    ggplot2::coord_cartesian(xlim = rng, ylim = rng) +
    ggplot2::labs(x = xlab, y = ylab,
                  title = sprintf("Calibration by risk rank (%d groups)", groups)) +
    ggplot2::theme_bw()

  attr(p, "data") <- tab
  p
}
