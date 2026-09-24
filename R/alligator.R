#' Alligator plot: cumulative incidence by risk group and event type
#'
#' @description
#' Plots the cause-specific cumulative incidence functions (CIFs) of the event
#' of interest and the competing event over time on a single panel, with risk
#' groups distinguished by color and event types by line type (event of
#' interest solid, competing event dashed). In higher-risk groups the
#' event-of-interest CIF rises above the competing-event CIF and the curves
#' separate like an alligator's open jaws.
#'
#' The x-axis is truncated at the model's readout time \code{t}, with tick marks
#' at 1-unit (e.g. 1-year) intervals. The y-axis auto-scales its cap to the
#' smallest of 25/50/75/100\% that covers the peak cumulative incidence (so
#' low-incidence settings are not compressed), with ticks at 10\% (or 5\% when
#' capped at 25\%). Numbers at risk for each group are shown in an aligned panel
#' below the plot. Gray's test across groups is attached.
#'
#' @param object a \code{"gcemod"} object from \code{\link{gcecox}},
#'   \code{\link{gcefg}}, or \code{\link{gcecif}}.
#' @param groups a \code{\link{gce_cutpoints}} result, a grouping factor of
#'   length n, or a single integer giving the number of risk groups
#'   (default \code{2}). When an integer is given, the groups are derived with
#'   \code{\link{gce_cutpoints}} using \code{method}.
#' @param method grouping rule used when \code{groups} is an integer:
#'   \code{"optimal"} (default) or \code{"quantile"} (equal-sized groups).
#'   Ignored when \code{groups} is a cutpoints object or a factor.
#' @param criterion optimal-search criterion when \code{groups} is an integer and
#'   \code{method = "optimal"}: \code{"omegaplus"} (default, log-\eqn{\omega^{+}}
#'   separation; not tail-sensitive) or \code{"cif"} (joint Gray statistic on the
#'   observed cumulative incidence of both events). See \code{\link{gce_cutpoints}}.
#' @param xlab,ylab axis labels.
#' @param group.lab legend title for the risk groups (default "Risk group").
#' @param risk.table logical; show numbers at risk in an aligned panel below the
#'   plot (default \code{TRUE}).
#' @param ymax optional numeric upper cap for the y-axis (0-1). If \code{NULL}
#'   (default) the cap auto-scales to 25/50/75/100\% based on the peak incidence.
#'
#' @return A \code{ggplot}/\code{patchwork} object. The Gray's-test result, the
#'   plotted CIF data, the numbers-at-risk table, and the cutpoints are attached
#'   as attributes \code{"grays_test"}, \code{"data"}, \code{"n.risk"}, and
#'   \code{"cutpoints"}.
#'
#' @examples
#' data(hn)
#' Ind <- data.frame(event = as.integer(hn$status == 1),
#'                   competing = as.integer(hn$status == 2))
#' fit <- gcecox(hn$time, Ind, hn[, c("age", "t_cat", "p16")], M = 5, t = 5)
#' gce_alligator(fit, groups = 2)
#' @export
gce_alligator <- function(object, groups = 2,
                          method = c("optimal", "quantile"),
                          criterion = c("omegaplus", "cif"),
                          xlab = "Time (years)", ylab = "Cumulative incidence",
                          group.lab = "Risk group", risk.table = TRUE,
                          ymax = NULL) {
  if (!inherits(object, "gcemod")) stop("'object' must be a gcemod object.", call. = FALSE)
  method <- match.arg(method)
  criterion <- match.arg(criterion)

  # Resolve groups, retaining the cutpoint(s) when available for annotation.
  if (is.list(groups) && !is.null(groups$groups)) {
    grp <- groups$groups; cutpoints <- groups$cutpoints
  } else if (length(groups) == length(object$riskscore)) {
    grp <- as.factor(groups); cutpoints <- NULL
  } else if (length(groups) == 1L && is.numeric(groups)) {
    cc <- gce_cutpoints(object, groups = groups, method = method, criterion = criterion)
    grp <- cc$groups; cutpoints <- cc$cutpoints
  } else {
    stop("'groups' must be a gce_cutpoints() result, a grouping factor, or a single integer.",
         call. = FALSE)
  }
  glev <- levels(grp)
  tmax <- object$t
  xbrk <- seq(0, floor(tmax + 1e-9), by = 1)

  ci <- cmprsk::cuminc(ftime = object$time, fstatus = object$cause,
                       group = as.integer(grp), cencode = 0)
  gtest <- ci[["Tests"]]
  els <- ci[setdiff(names(ci), "Tests")]

  # CIF curves, truncated at tmax and carried forward to the right edge.
  df <- do.call(rbind, lapply(names(els), function(nm) {
    parts <- strsplit(nm, " ", fixed = TRUE)[[1]]
    gi <- as.integer(parts[1]); cj <- parts[length(parts)]
    if (!cj %in% c("1", "2")) return(NULL)
    e <- els[[nm]]
    keep <- e$time <= tmax
    tt <- e$time[keep]; est <- e$est[keep]
    if (length(tt) == 0L) { tt <- 0; est <- 0 }
    tt <- c(tt, tmax); est <- c(est, est[length(est)])   # carry forward to tmax
    data.frame(
      group = factor(glev[gi], levels = glev),
      event = factor(if (cj == "1") "Event of interest" else "Competing event",
                     levels = c("Event of interest", "Competing event")),
      time  = tt, cif = est, stringsAsFactors = FALSE)
  }))

  # Numbers at risk (subjects with follow-up >= each time), per group.
  ng0 <- length(glev)
  nrisk <- data.frame(
    xt    = rep(xbrk, times = ng0),
    group = factor(rep(glev, each = length(xbrk)), levels = glev),
    stringsAsFactors = FALSE)
  nrisk$n <- mapply(function(x, g) sum(object$time[grp == g] >= x, na.rm = TRUE),
                    nrisk$xt, as.character(nrisk$group))

  # Adaptive y-axis cap: smallest of 25/50/75/100% that covers the peak CIF.
  peak <- suppressWarnings(max(df$cif, na.rm = TRUE))
  if (!is.finite(peak)) peak <- 1
  cap <- if (!is.null(ymax)) ymax else {
    br <- c(0.25, 0.5, 0.75, 1.0)
    br[which(br >= min(peak * 1.02, 1))][1]
  }
  if (length(cap) == 0L || is.na(cap)) cap <- 1
  ystep <- if (cap <= 0.25) 0.05 else 0.10
  yb <- seq(0, cap, by = ystep)

  # Subtitle: cutpoint(s) on the risk-score scale and group sizes.
  gsz <- as.integer(table(grp))
  sizes_txt <- paste(paste0(glev, ": n=", gsz), collapse = ",  ")
  sub <- if (!is.null(cutpoints) && length(cutpoints) > 0) {
    paste0("Risk-score cutpoint", if (length(cutpoints) > 1) "s" else "", " = ",
           paste(round(cutpoints, 3), collapse = ", "),
           "   |   group sizes: ", sizes_txt)
  } else {
    paste0("Group sizes: ", sizes_txt)
  }

  p_main <- ggplot2::ggplot(
        df, ggplot2::aes(x = .data[["time"]], y = .data[["cif"]],
                         colour = .data[["group"]],
                         linetype = .data[["event"]])) +
    ggplot2::geom_step() +
    ggplot2::scale_linetype_manual(
      values = c("Event of interest" = "solid", "Competing event" = "dashed")) +
    ggplot2::scale_x_continuous(breaks = xbrk) +
    ggplot2::scale_y_continuous(breaks = yb, labels = paste0(round(yb * 100), "%")) +
    ggplot2::coord_cartesian(xlim = c(0, tmax), ylim = c(0, cap)) +
    ggplot2::labs(y = ylab, colour = group.lab, linetype = NULL,
                  title = "Cumulative incidence by risk group", subtitle = sub) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")

  if (!isTRUE(risk.table)) {
    p <- p_main + ggplot2::labs(x = xlab)
    attr(p, "grays_test") <- gtest
    attr(p, "data") <- df
    attr(p, "n.risk") <- nrisk
    attr(p, "cutpoints") <- cutpoints
    return(p)
  }

  # Numbers-at-risk table as an aligned panel below the plot.
  p_main <- p_main +
    ggplot2::labs(x = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  nrisk$row <- factor(as.character(nrisk$group), levels = rev(glev))
  p_tbl <- ggplot2::ggplot() +
    ggplot2::geom_text(
      data = nrisk,
      mapping = ggplot2::aes(x = .data[["xt"]], y = .data[["row"]],
                             label = .data[["n"]], colour = .data[["group"]]),
      size = 3, show.legend = FALSE) +
    ggplot2::scale_x_continuous(breaks = xbrk) +
    ggplot2::coord_cartesian(xlim = c(0, tmax)) +
    ggplot2::labs(x = xlab, y = NULL, title = "No. at risk") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 10, face = "italic"),
                   legend.position = "none")

  ng <- length(glev)
  p <- patchwork::wrap_plots(p_main, p_tbl, ncol = 1,
                             heights = c(1, 0.08 * ng + 0.10))

  attr(p, "grays_test") <- gtest
  attr(p, "data") <- df
  attr(p, "n.risk") <- nrisk
  attr(p, "cutpoints") <- cutpoints
  p
}
