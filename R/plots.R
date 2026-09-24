# Internal ggplot helpers. Not exported.

#' Simple scatter plot used for omega / omega+ displays
#' @keywords internal
#' @noRd
.gce_scatter <- function(x, y, xlab, ylab, ylim = NULL) {
  df <- data.frame(x = x, y = y)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(alpha = 0.4) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw()
  if (!is.null(ylim)) p <- p + ggplot2::ylim(ylim[1], ylim[2])
  p
}
