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
