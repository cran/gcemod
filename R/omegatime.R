#' GCE contrasts over time: a visual check of relative-hazard constancy
#'
#' @description
#' Plots a nonparametric estimate of a chosen GCE contrast over follow-up, as a
#' visual assessment of the proportional relative-hazards (PRH) assumption and
#' of the time range over which it holds approximately. All estimates are built
#' from the Nelson-Aalen cause-specific cumulative hazards
#' \eqn{\Lambda_k(t)} and the Aalen-Johansen cumulative incidences
#' \eqn{F_k(t)}:
#' \describe{
#'   \item{\code{"omega"} (default)}{\eqn{\Lambda_1(t)/\{\Lambda_1(t)+\Lambda_2(t)\}},
#'     the bounded cumulative hazard share, plotted on \eqn{[0,1]}.}
#'   \item{\code{"omegaplus"}}{\eqn{\omega^{+}(t)=\Lambda_1(t)/\Lambda_2(t)},
#'     the cause-specific relative hazard, on a log scale.}
#'   \item{\code{"omegatilde"}}{\eqn{\tilde\Lambda_1(t)/\{\tilde\Lambda_1(t)+\tilde\Lambda_2(t)\}},
#'     the subdistribution analogue of \code{"omega"}, with
#'     \eqn{\tilde\Lambda_k(t) = -\log\{1-F_k(t)\}}, on \eqn{[0,1]}.}
#'   \item{\code{"omegatildeplus"}}{\eqn{\tilde\omega^{+}(t)=\tilde\Lambda_1(t)/\tilde\Lambda_2(t)},
#'     the subdistribution relative hazard, on a log scale.}
#'   \item{\code{"rho"}}{\eqn{\rho(t)=F_1(t)/F_2(t)}, the cumulative-incidence
#'     ratio, on a log scale.}
#'   \item{\code{"pi"}}{\eqn{\pi(t)=F_1(t)/\{F_1(t)+F_2(t)\} = \rho(t)/\{1+\rho(t)\}},
#'     the bounded share form of \eqn{\rho} (so \eqn{\mathrm{logit}\,\pi = \log\rho}),
#'     i.e. the probability that a realized event by \eqn{t} is the event of
#'     interest, on \eqn{[0,1]}.}
#' }
#'
#' If the cause-specific relative hazard is constant in time then
#' \eqn{\Lambda_1(t)=\omega^{+}\Lambda_2(t)} and, because
#' \eqn{F_k(t)=\int_0^t S(u^-)\lambda_k(u)\,du}, also
#' \eqn{F_1(t)=\omega^{+}F_2(t)}: the \code{"omega"}, \code{"omegaplus"},
#' \code{"rho"}, and \code{"pi"} displays are then exactly flat, and drift in any of them
#' delineates where a single time-constant summary is and is not adequate. The
#' subdistribution displays are flat only approximately (exactly so when both
#' incidences are small, since \eqn{\tilde\Lambda_k \approx F_k}).
#'
#' With \code{groups}, one curve is drawn per risk group. Approximately
#' parallel group curves on a log-scaled contrast indicate a time-constant
#' group \emph{contrast}---the modeling assumption of the stacked
#' estimators---even where the individual curves drift.
#'
#' @param object a \code{"gcemod"} object from \code{\link{gcecox}},
#'   \code{\link{gcefg}}, or \code{\link{gcecif}}.
#' @param which the contrast to plot: \code{"omega"} (default),
#'   \code{"omegaplus"}, \code{"omegatilde"}, \code{"omegatildeplus"},
#'   \code{"rho"}, or \code{"pi"}. See Description.
#' @param groups optional grouping: a \code{\link{gce_cutpoints}} result, a
#'   grouping factor of length n, or a single integer number of risk groups
#'   (derived with \code{\link{gce_cutpoints}}). \code{NULL} (default) draws
#'   one curve for the whole sample.
#' @param method grouping rule used when \code{groups} is an integer:
#'   \code{"optimal"} (default) or \code{"quantile"}.
#' @param xmax optional upper limit of the time axis; the curves are truncated
#'   at \code{xmax}. Defaults to \code{NULL} (full follow-up). Late follow-up
#'   is often carried by few subjects at risk, so truncating (e.g. at or
#'   somewhat beyond the readout time \code{t}) can give a fairer picture of
#'   the range that matters for the fitted summary.
#'
#' @return A \code{ggplot} object. The plotted step-function data (group, time,
#'   value) are attached as attribute \code{"data"}. The model's readout time
#'   \code{t} is marked with a dashed vertical line.
#'
#' @examples
#' data(hn)
#' Ind <- data.frame(event = as.integer(hn$status == 1),
#'                   competing = as.integer(hn$status == 2))
#' fit <- gcecox(hn$time, Ind, hn[, c("age", "t_cat", "p16")], M = 5, t = 5)
#' gce_omegatime(fit)
#' gce_omegatime(fit, which = "rho")
#' gce_omegatime(fit, which = "omegaplus", groups = 2)
#' @seealso \code{\link{gce_cutpoints}}, \code{\link{gce_alligator}}
#' @export
gce_omegatime <- function(object,
                          which = c("omega", "omegaplus", "omegatilde",
                                    "omegatildeplus", "rho", "pi"),
                          groups = NULL, method = "optimal", xmax = NULL) {
  which <- match.arg(which)
  if (!inherits(object, "gcemod"))
    stop("'object' must be a gcemod object.", call. = FALSE)
  time  <- object$time
  cause <- object$cause
  grp <- if (is.null(groups)) factor(rep.int("All", length(time)))
         else factor(.resolve_groups(object, groups, method = method))
  is_share <- which %in% c("omega", "omegatilde", "pi")

  one <- function(g) {
    sel <- grp == g
    tt <- time[sel]; cc <- cause[sel]
    ut <- sort(unique(tt[cc > 0]))
    if (length(ut) < 2L) return(NULL)
    Y  <- vapply(ut, function(u) sum(tt >= u), numeric(1))
    d1 <- vapply(ut, function(u) sum(tt == u & cc == 1L), numeric(1))
    d2 <- vapply(ut, function(u) sum(tt == u & cc == 2L), numeric(1))
    v <- switch(which,
      omega     = { L1 <- cumsum(d1 / Y); L2 <- cumsum(d2 / Y)
                    L1 / (L1 + L2) },
      omegaplus = { L1 <- cumsum(d1 / Y); L2 <- cumsum(d2 / Y)
                    L1 / L2 },
      {  # CIF-based contrasts: Aalen-Johansen
        S  <- cumprod(1 - (d1 + d2) / Y)
        Sm <- c(1, S[-length(S)])          # S(u-)
        F1 <- cumsum(Sm * d1 / Y)
        F2 <- cumsum(Sm * d2 / Y)
        switch(which,
          rho            = F1 / F2,
          pi             = F1 / (F1 + F2),
          omegatilde     = { Lt1 <- -log(pmax(1 - F1, 1e-12))
                             Lt2 <- -log(pmax(1 - F2, 1e-12))
                             Lt1 / (Lt1 + Lt2) },
          omegatildeplus = { Lt1 <- -log(pmax(1 - F1, 1e-12))
                             Lt2 <- -log(pmax(1 - F2, 1e-12))
                             Lt1 / Lt2 })
      })
    keep <- is.finite(v) & v > 0
    if (is_share) keep <- is.finite(v)
    keep <- keep & !is.na(v)
    if (!any(keep)) return(NULL)
    data.frame(group = as.character(g), time = ut[keep], value = v[keep],
               stringsAsFactors = FALSE)
  }
  dat <- do.call(rbind, lapply(levels(grp), one))
  if (!is.null(xmax)) {
    if (!is.numeric(xmax) || length(xmax) != 1L || xmax <= 0)
      stop("'xmax' must be a single positive number.", call. = FALSE)
    if (!is.null(dat)) dat <- dat[dat$time <= xmax, , drop = FALSE]
  }
  if (is.null(dat) || !nrow(dat))
    stop("Not enough events to plot the requested contrast over time.",
         call. = FALSE)
  dat$group <- factor(dat$group, levels = levels(grp))

  ylab <- switch(which,
    omega          = expression(omega(t)),
    omegaplus      = expression(omega^"+"*(t)),
    omegatilde     = expression(tilde(omega)(t)),
    omegatildeplus = expression(tilde(omega)^"+"*(t)),
    rho            = expression(rho(t)),
    pi             = expression(pi(t)))
  p <- if (nlevels(grp) > 1L) {
    ggplot2::ggplot(dat, ggplot2::aes(x = .data[["time"]], y = .data[["value"]],
                                      color = .data[["group"]]))
  } else {
    ggplot2::ggplot(dat, ggplot2::aes(x = .data[["time"]], y = .data[["value"]]))
  }
  p <- p + ggplot2::geom_step(linewidth = 0.7) +
    ggplot2::labs(x = "Time", y = ylab, color = "Risk group") +
    ggplot2::theme_bw()
  p <- if (is_share) {
    p + ggplot2::scale_y_continuous(limits = c(0, 1))
  } else {
    p + ggplot2::scale_y_log10()
  }
  if (!is.null(object$t) && is.finite(object$t)) {
    p <- p + ggplot2::geom_vline(xintercept = object$t, linetype = "dashed",
                                 alpha = 0.5)
  }
  attr(p, "data") <- dat
  p
}
