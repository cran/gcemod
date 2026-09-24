# Between-group sum of squares of a value vector for a grouping (ANOVA numerator).
# Measured on the risk-score / log-omega+ scale it is not tail-sensitive.
.between_ss <- function(y, grp) {
  gm <- mean(y, na.rm = TRUE)
  parts <- tapply(y, grp, function(v) length(v) * (mean(v, na.rm = TRUE) - gm)^2)
  sum(parts, na.rm = TRUE)
}

# Joint Gray's test statistic across BOTH competing events (event of interest +
# competing event), summed. Returns a list(stat, df, p.value).
.cif_joint <- function(time, cause, grp) {
  if (length(unique(grp[!is.na(grp)])) < 2L)
    return(list(stat = NA_real_, df = NA_real_, p.value = NA_real_))
  ci <- tryCatch(cmprsk::cuminc(ftime = time, fstatus = cause,
                                group = grp, cencode = 0),
                 error = function(e) NULL)
  if (is.null(ci) || is.null(ci[["Tests"]]))
    return(list(stat = NA_real_, df = NA_real_, p.value = NA_real_))
  tst <- ci[["Tests"]]
  rows <- rownames(tst) %in% c("1", "2")     # the two competing events
  stat <- sum(tst[rows, "stat"], na.rm = TRUE)
  df   <- sum(tst[rows, "df"],   na.rm = TRUE)
  list(stat = stat, df = df,
       p.value = stats::pchisq(stat, df = df, lower.tail = FALSE))
}

.assign_groups <- function(score, cuts) {
  cut(score, breaks = c(-Inf, sort(cuts), Inf), labels = FALSE,
      include.lowest = TRUE)
}

# Univariable Lunn-McNeil omega+ contrast between risk groups: tests whether the
# ratio of cause-specific (or subdistribution) hazards differs across groups and
# returns the omega+ ratio estimate. This is the GCE test that a per-event Gray
# or log-rank test misses.
.omegaplus_group_test <- function(time, cause, grp, type) {
  gf <- data.frame(risk_group = factor(grp))
  tryCatch(lunnmcneil(time, cause, gf, type = type,
                      standardize = FALSE, cluster.se = TRUE),
           error = function(e) NULL)
}

# Per-event Gray's test p-values (event of interest and competing), to show that
# individual-event tests can be non-significant when the omega+ contrast is not.
.individual_gray <- function(time, cause, grp) {
  ci <- tryCatch(cmprsk::cuminc(ftime = time, fstatus = cause,
                                group = grp, cencode = 0),
                 error = function(e) NULL)
  if (is.null(ci) || is.null(ci[["Tests"]])) return(NULL)
  tst <- ci[["Tests"]]
  grab <- function(cc, col) { i <- match(cc, rownames(tst))
    if (is.na(i)) NA_real_ else as.numeric(tst[i, col]) }
  data.frame(event = c("event of interest", "competing event"),
             gray.stat = c(grab("1", "stat"), grab("2", "stat")),
             gray.p    = c(grab("1", "pv"),   grab("2", "pv")),
             stringsAsFactors = FALSE)
}

# Greedy forward search for cutpoints maximizing a criterion function crit(grp).
.search_cuts <- function(score, crit, groups, min.prop) {
  n <- length(score)
  cand <- unique(stats::quantile(score, probs = seq(min.prop, 1 - min.prop, by = 0.01),
                                 na.rm = TRUE))
  cuts <- numeric(0); best_val <- NA_real_
  for (k in seq_len(groups - 1L)) {
    best <- NA_real_; bv <- -Inf
    for (cc in setdiff(cand, cuts)) {
      trial <- sort(c(cuts, cc))
      grp <- .assign_groups(score, trial)
      props <- table(factor(grp, levels = seq_len(length(trial) + 1L))) / n
      if (any(props < min.prop)) next
      v <- crit(grp)
      if (is.finite(v) && v > bv) { bv <- v; best <- cc }
    }
    if (is.na(best)) return(NULL)
    cuts <- sort(c(cuts, best)); best_val <- bv
  }
  list(cuts = cuts, value = best_val)
}

# Resolve a grouping argument into a factor of length n.
.resolve_groups <- function(object, groups, method = "quantile", criterion = "omegaplus") {
  if (is.list(groups) && !is.null(groups$groups)) return(groups$groups)
  if (length(groups) == length(object$riskscore)) return(as.factor(groups))
  if (length(groups) == 1L && is.numeric(groups)) {
    return(gce_cutpoints(object, groups = groups, method = method,
                         criterion = criterion)$groups)
  }
  stop("'groups' must be a gce_cutpoints() result, a grouping factor, or a single integer.",
       call. = FALSE)
}

#' Risk-score cutpoints for GCE risk stratification
#'
#' @description
#' Splits subjects into ordered risk groups on the GCE risk score. With
#' \code{method = "optimal"} the cutpoint(s) maximize a separation criterion by a
#' greedy forward search; with \code{method = "quantile"} equal-sized groups are
#' used.
#'
#' The default criterion, \code{"omegaplus"}, maximizes the between-group
#' separation in \strong{log} \eqn{\omega^{+}} (equivalently, the risk score /
#' linear predictor). This is the additive scale on which the joint
#' event-of-interest-vs-competing \eqn{\omega^{+}} contrast is defined, and --
#' unlike separation measured on the exponential \eqn{\omega^{+}} scale -- it is
#' \strong{not} sensitive to a few extreme high-risk subjects, so the search does
#' not chase high-risk tails and tends to return balanced, well-separated groups.
#' The alternative \code{"cif"} maximizes a joint Gray's test statistic (summed
#' across both competing events) on the observed cumulative incidence.
#'
#' @param object a \code{"gcemod"} object from \code{\link{gcecox}},
#'   \code{\link{gcefg}}, or \code{\link{gcecif}}. For a \code{\link{gcecif}} fit
#'   the per-group summary reports mean \eqn{\rho} and share instead of
#'   \eqn{\omega^{+}} / \eqn{\omega}; the \code{"omegaplus"} criterion still means
#'   separation on the risk-score (linear-predictor) scale.
#' @param groups number of groups (default \code{2}).
#' @param method \code{"optimal"} (default) or \code{"quantile"}.
#' @param criterion optimal-search criterion: \code{"omegaplus"} (default,
#'   log-\eqn{\omega^{+}} separation) or \code{"cif"} (joint Gray statistic on
#'   observed cumulative incidence of both events).
#' @param min.prop minimum proportion of subjects per group (default \code{0.1}).
#' @param nperm number of permutations for the selection-adjusted p-value
#'   (default \code{0} = not computed). When \code{nperm > 0}, the
#'   maximally-selected joint Gray statistic is permuted (risk score vs. outcome)
#'   to give a p-value accounting for the data-driven choice of cutpoint.
#' @param labels optional character labels for the groups (low to high risk).
#'
#' @return A list with: \code{groups} (ordered factor), \code{cutpoints} (on the
#'   risk-score scale), \code{summary} (per-group n, events, mean
#'   \eqn{\omega^{+}} / \eqn{\omega}), and \code{test}. The \code{test} element
#'   reports the Lunn-McNeil \eqn{\omega^{+}} contrast between groups -- a test
#'   that the ratio of cause-specific (or subdistribution) hazards differs across
#'   groups -- with the univariable \eqn{\omega^{+}} ratio estimate and CI
#'   (\code{$omegaplus}), the joint test (\code{$statistic}, \code{$df},
#'   \code{$p.value}), the per-event Gray tests for comparison
#'   (\code{$individual_gray}), and a selection-adjusted \code{$p.perm} when
#'   \code{nperm > 0}.
#'
#' @examples
#' data(hn)
#' Ind <- data.frame(event = as.integer(hn$status == 1),
#'                   competing = as.integer(hn$status == 2))
#' fit <- gcecox(hn$time, Ind, hn[, c("age", "t_cat", "p16")], M = 5, t = 5)
#' gce_cutpoints(fit, groups = 2)
#' \donttest{
#' # with a selection-adjusted permutation p-value (slower)
#' gce_cutpoints(fit, groups = 2, nperm = 200)
#' }
#' @export
gce_cutpoints <- function(object, groups = 2,
                          method = c("optimal", "quantile"),
                          criterion = c("omegaplus", "cif"),
                          min.prop = 0.1, nperm = 0, labels = NULL) {
  method <- match.arg(method)
  criterion <- match.arg(criterion)
  if (!inherits(object, "gcemod")) stop("'object' must be a gcemod object.", call. = FALSE)
  if (groups < 2L) stop("'groups' must be >= 2.", call. = FALSE)

  score <- object$riskscore                 # log-omega+ scale (linear predictor)
  time  <- object$time
  cause <- object$cause

  if (method == "quantile") {
    br <- stats::quantile(score, probs = seq(0, 1, length.out = groups + 1))
    cutpoints <- unname(br[-c(1, length(br))])
  } else {
    crit <- if (criterion == "omegaplus") {
      function(grp) .between_ss(score, grp)             # separation in log omega+
    } else {
      function(grp) .cif_joint(time, cause, grp)$stat   # joint Gray statistic
    }
    res <- .search_cuts(score, crit, groups, min.prop)
    if (is.null(res)) stop("Could not find a valid split; try fewer groups or a smaller min.prop.",
                           call. = FALSE)
    cutpoints <- res$cuts
  }

  grp <- factor(.assign_groups(score, cutpoints), levels = seq_len(groups))
  # Default risk-direction labels for 2-3 groups (group 1 = lowest risk score).
  if (is.null(labels)) {
    labels <- switch(as.character(groups),
                     "2" = c("Low", "High"),
                     "3" = c("Low", "Intermediate", "High"),
                     NULL)
  }
  if (!is.null(labels)) {
    if (length(labels) != nlevels(grp))
      stop("length(labels) must equal the number of groups.", call. = FALSE)
    levels(grp) <- labels
  }

  ev <- as.integer(cause == 1)
  cm <- as.integer(cause == 2)
  is_cif <- inherits(object, "gcecif")
  summ <- do.call(rbind, lapply(levels(grp), function(g) {
    sel <- grp == g
    base <- data.frame(group = g, n = sum(sel),
                       events_of_interest = sum(ev[sel]),
                       competing_events   = sum(cm[sel]),
                       stringsAsFactors = FALSE)
    if (is_cif) {
      rho <- object$rho
      base$mean_rho   <- mean(rho[sel], na.rm = TRUE)
      base$mean_share <- mean((rho / (1 + rho))[sel], na.rm = TRUE)
    } else {
      base$mean_omegaplus <- mean(object$omegaplus[sel], na.rm = TRUE)
      base$mean_omega     <- mean(object$omega[sel], na.rm = TRUE)
    }
    base
  }))
  rownames(summ) <- NULL

  # Reported test: Lunn-McNeil omega+ contrast between groups (the difference in
  # relative cause-specific / subdistribution hazards), plus the univariable
  # omega+ ratio estimate; and the per-event Gray tests for comparison.
  ltype <- if (!is.null(object$lunnmcneil)) object$lunnmcneil$type else "coxph"
  scale_lab <- if (identical(ltype, "finegray")) "relative subdistribution hazards"
               else "relative cause-specific hazards"
  lm_g <- .omegaplus_group_test(time, cause, as.integer(grp), ltype)
  test <- list(
    method = paste0("Lunn-McNeil omega+ contrast (", scale_lab, ")"),
    omegaplus = if (!is.null(lm_g)) lm_g$omegaplus else NULL,   # ratio(s), CI, p
    statistic = if (!is.null(lm_g)) unname(lm_g$omnibus["statistic"]) else NA_real_,
    df        = if (!is.null(lm_g)) unname(lm_g$omnibus["df"]) else NA_real_,
    p.value   = if (!is.null(lm_g)) unname(lm_g$omnibus["p.value"]) else NA_real_,
    individual_gray = .individual_gray(time, cause, as.integer(grp)),
    p.perm = NA_real_)

  # Selection-adjusted permutation p-value for the omega+ group contrast.
  if (nperm > 0 && method == "optimal" && !is.na(test$statistic)) {
    crit <- if (criterion == "omegaplus") function(g) .between_ss(score, g)
            else function(g) .cif_joint(time, cause, g)$stat
    group_stat <- function(sc) {
      r <- .search_cuts(sc, crit, groups, min.prop)
      if (is.null(r)) return(NA_real_)
      g <- .assign_groups(sc, r$cuts)
      lg <- .omegaplus_group_test(time, cause, g, ltype)
      if (is.null(lg)) NA_real_ else unname(lg$omnibus["statistic"])
    }
    Tperm <- vapply(seq_len(nperm), function(b) group_stat(sample(score)), numeric(1))
    ok <- sum(!is.na(Tperm))
    test$p.perm <- (1 + sum(Tperm >= test$statistic, na.rm = TRUE)) / (ok + 1)
  }

  list(groups = grp, cutpoints = cutpoints, method = method,
       criterion = if (method == "optimal") criterion else NA,
       summary = summ, test = test)
}
