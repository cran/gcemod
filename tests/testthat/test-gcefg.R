test_that("gcefg runs and returns an omega+ gcemod object", {
  skip_on_cran()
  d <- make_test_data(n = 200)
  Ind <- data.frame(d$CA, d$CM)
  fit <- gcefg(d$Time, Ind, d$Cov, M = 4, t = 5)

  expect_s3_class(fit, "gcemod")
  expect_s3_class(fit, "gcefg")
  expect_length(fit$coef, ncol(d$Cov))
  expect_equal(ncol(fit$result), 4L)
  expect_s3_class(fit$omegatimeplot, "ggplot")
  expect_true(all(is.finite(fit$omega)))
})

test_that("gcefg se = 'boot' returns bootstrap standard errors", {
  skip_on_cran()
  set.seed(11)
  d <- make_test_data(n = 150)
  Ind <- data.frame(d$CA, d$CM)
  fit <- suppressWarnings(gcefg(d$Time, Ind, d$Cov, M = 4, t = 5,
                                se = "boot", B = 15))
  expect_identical(fit$lunnmcneil$se_method, "boot")
  expect_identical(fit$lunnmcneil$B, 15L)
  expect_true(all(is.finite(fit$result[, "lower .95"])))
  expect_true(all(fit$result[, "lower .95"] < fit$result[, "upper .95"]))
})

test_that("gce_omegatime plots omega(t) on a gcefg fit", {
  skip_on_cran()
  d <- make_test_data(n = 200)
  Ind <- data.frame(d$CA, d$CM)
  fit <- gcefg(d$Time, Ind, d$Cov, M = 4, t = 5)

  p <- gce_omegatime(fit)
  expect_s3_class(p, "ggplot")
  dat <- attr(p, "data")
  expect_true(is.data.frame(dat) && nrow(dat) > 0)
  expect_true(all(dat$value >= 0 & dat$value <= 1))

  p2 <- gce_omegatime(fit, which = "omegaplus", groups = 2)
  expect_s3_class(p2, "ggplot")
  expect_equal(nlevels(attr(p2, "data")$group), 2L)

  for (w in c("omegatilde", "omegatildeplus", "rho", "pi")) {
    pw <- gce_omegatime(fit, which = w)
    expect_s3_class(pw, "ggplot")
    expect_true(all(is.finite(attr(pw, "data")$value)))
  }
})

test_that("gce_calibration and gce_alligator work on a gcefg fit", {
  skip_on_cran()
  d <- make_test_data(n = 200)
  Ind <- data.frame(d$CA, d$CM)
  fit <- gcefg(d$Time, Ind, d$Cov, M = 4, t = 5)
  expect_s3_class(gce_calibration(fit, groups = 4), "ggplot")
  al <- gce_alligator(fit, groups = 3)
  expect_true(inherits(al, "ggplot") || inherits(al, "patchwork"))
})

test_that("gce_tidy(effects = TRUE) uses Fine-Gray fits on a gcefg object", {
  skip_on_cran()
  d <- make_test_data(n = 200)
  Ind <- data.frame(d$CA, d$CM)
  fit <- gcefg(d$Time, Ind, d$Cov, M = 4, t = 5)

  eff <- gce_tidy(fit, effects = TRUE)
  expect_identical(names(eff),
                   c("term", "primary", "competing", "total", "relative"))
  expect_equal(nrow(eff), ncol(d$Cov))
  # primary column is drawn from the stored Fine-Gray subdistribution fit (fit1)
  hr1 <- sprintf("%.2f", exp(as.numeric(fit$fit1$coef)))
  expect_true(all(mapply(startsWith, eff$primary, hr1)))
})
