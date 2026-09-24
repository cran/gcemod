test_that("gcecox returns a well-formed omega+ gcemod object", {
  d <- make_test_data()
  Ind <- data.frame(d$CA, d$CM)
  fit <- gcecox(d$Time, Ind, d$Cov, M = 4, t = 5)

  expect_s3_class(fit, "gcemod")
  expect_s3_class(fit, "gcecox")
  expect_length(fit$coef, ncol(d$Cov))
  expect_equal(ncol(fit$result), 4L)                 # includes P-value column
  expect_true(all(fit$result[, "P-value"] >= 0 & fit$result[, "P-value"] <= 1))
  expect_s3_class(fit$omegaplot, "ggplot")
  expect_true(all(fit$omega > 0 & fit$omega < 1))
  expect_null(fit$result1)                            # omega regression removed
})

test_that("gcecox omega+ coef matches the Lunn-McNeil engine", {
  d <- make_test_data()
  Ind <- data.frame(d$CA, d$CM)
  fit <- gcecox(d$Time, Ind, d$Cov, M = 4, t = 5)
  lmfit <- lunnmcneil(d$Time, d$cause, d$Cov, type = "coxph")
  expect_equal(unname(fit$coef), lmfit$omegaplus$log_ratio, tolerance = 1e-6)
})

test_that("constant covariate is dropped with a warning, not an error", {
  d <- make_test_data()
  Cov <- cbind(d$Cov, const = 0)
  Ind <- data.frame(d$CA, d$CM)
  expect_warning(fit <- gcecox(d$Time, Ind, Cov, M = 4, t = 5), "constant")
  expect_false("const" %in% fit$covnames)
})

test_that("cutpoints, alligator, calibration, performance work on a gcecox fit", {
  d <- make_test_data()
  Ind <- data.frame(d$CA, d$CM)
  fit <- gcecox(d$Time, Ind, d$Cov, M = 4, t = 5)

  cut <- gce_cutpoints(fit, groups = 3, method = "optimal")
  expect_equal(nlevels(cut$groups), 3L)
  expect_equal(sum(cut$summary$n), fit$n)

  al <- gce_alligator(fit, groups = cut)
  expect_true(inherits(al, "ggplot") || inherits(al, "patchwork"))
  expect_s3_class(gce_calibration(fit, which = "omegaplus", groups = 4), "ggplot")

  td <- gce_tidy(fit)
  expect_equal(nrow(td), ncol(d$Cov))
  expect_true(all(c("term", "ratio", "conf.low", "conf.high", "p.value") %in% names(td)))

  ci <- gce_performance(fit)
  expect_gte(ci["estimate"], 0.5)
})

test_that("gce_tidy(effects = TRUE) returns a primary/competing/total/relative table", {
  d <- make_test_data()
  Ind <- data.frame(d$CA, d$CM)
  fit <- gcecox(d$Time, Ind, d$Cov, M = 4, t = 5)

  expect_false(is.null(fit$design))                    # design retained for refits
  eff <- gce_tidy(fit, effects = TRUE)
  expect_equal(nrow(eff), ncol(d$Cov))
  expect_identical(names(eff),
                   c("term", "primary", "competing", "total", "relative"))
  expect_true(all(grepl("^[0-9.]+ \\([0-9.]+, [0-9.]+\\)$", eff$primary)))
  expect_identical(eff$term, fit$covnames)
})

test_that("non-syntactic covariate names are handled end to end", {
  d <- make_test_data()
  Cov <- d$Cov
  names(Cov)[names(Cov) == "age"] <- "age at dx"   # space => non-syntactic name
  Ind <- data.frame(d$CA, d$CM)

  fit <- gcecox(d$Time, Ind, Cov, M = 4, t = 5)     # must not error
  expect_s3_class(fit, "gcemod")
  expect_length(fit$coef, ncol(Cov))

  # the Lunn-McNeil engine agrees regardless of the odd name
  lmfit <- lunnmcneil(d$Time, d$cause, Cov, type = "coxph")
  expect_equal(unname(fit$coef), lmfit$omegaplus$log_ratio, tolerance = 1e-6)

  # scoring new data carrying the same odd name works
  rs <- gce_riskscore(fit, newdata = Cov[1:5, , drop = FALSE])
  expect_equal(nrow(rs), 5L)
  expect_true(all(is.finite(rs$riskscore)))
})
