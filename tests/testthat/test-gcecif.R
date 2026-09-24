test_that("gcecif returns a well-formed gcemod object", {
  d <- make_test_data(n = 300)
  Ind <- data.frame(d$CA, d$CM)
  t0  <- stats::median(d$Time)
  fit <- gcecif(d$Time, Ind, d$Cov, t = t0)

  expect_s3_class(fit, "gcemod")
  expect_s3_class(fit, "gcecif")
  expect_true(fit$F1 > 0 && fit$F2 > 0)
  expect_equal(fit$rho0, fit$F1 / fit$F2, tolerance = 1e-10)
  expect_length(fit$pseudo, length(d$Time))
  expect_length(fit$riskscore, length(d$Time))

  # result table: one row per covariate term, intercept excluded
  expect_true(all(c("exp(coef) (rho ratio)", "lower .95", "upper .95", "P-value")
                  %in% colnames(fit$result)))
  expect_false("(Intercept)" %in% rownames(fit$result))
})

test_that("gcecif integrates with gce_tidy and gce_riskscore", {
  d <- make_test_data(n = 300)
  Ind <- data.frame(d$CA, d$CM)
  fit <- gcecif(d$Time, Ind, d$Cov, t = stats::median(d$Time))

  td <- gce_tidy(fit)
  expect_true(all(c("term", "log_ratio", "ratio", "conf.low", "conf.high",
                    "p.value") %in% names(td)))
  expect_equal(nrow(td), nrow(fit$result))

  rs <- gce_riskscore(fit, newdata = d$Cov[1:5, , drop = FALSE])
  expect_true("rho_rel" %in% names(rs))
  expect_equal(nrow(rs), 5L)
})

test_that("rho0 approximates omega+ under a constant cause-specific ratio", {
  skip_on_cran()
  set.seed(11)
  n <- 6000; wplus <- 1.5; h2 <- 0.15; h1 <- wplus * h2; hc <- 0.05
  t1 <- rexp(n, h1); t2 <- rexp(n, h2); tc <- rexp(n, hc)
  Time <- pmin(t1, t2, tc)
  cause <- ifelse(t1 <= t2 & t1 <= tc, 1L, ifelse(t2 < t1 & t2 <= tc, 2L, 0L))
  Ind <- data.frame(event = as.integer(cause == 1), competing = as.integer(cause == 2))
  Cov <- data.frame(x = rnorm(n))                 # null covariate
  fit <- gcecif(Time, Ind, Cov, t = 5)
  expect_equal(fit$rho0, wplus, tolerance = 0.1)  # rho(t) tracks omega+ = 1.5
})
