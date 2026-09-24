test_that("Cox Lunn-McNeil interaction equals difference of separate cause-specific fits", {
  d <- make_test_data()
  # Identity holds on the raw covariate scale, so disable standardization here.
  lmfit <- lunnmcneil(d$Time, d$cause, d$Cov, type = "coxph",
                      cluster.se = FALSE, standardize = FALSE)

  b_event <- coef(survival::coxph(
    survival::Surv(d$Time, as.integer(d$cause == 1)) ~ smoke + age + bmi, data = d$Cov))
  b_comp  <- coef(survival::coxph(
    survival::Surv(d$Time, as.integer(d$cause == 2)) ~ smoke + age + bmi, data = d$Cov))
  expected <- b_event - b_comp

  got <- lmfit$omegaplus$log_ratio
  names(got) <- lmfit$omegaplus$term
  expect_equal(unname(got[names(expected)]), unname(expected), tolerance = 1e-4)
})

test_that("Cox Lunn-McNeil returns a valid omnibus test", {
  d <- make_test_data()
  lmfit <- lunnmcneil(d$Time, d$cause, d$Cov, type = "coxph")
  expect_s3_class(lmfit, "gce_lunnmcneil")
  expect_equal(as.integer(lmfit$omnibus["df"]), ncol(d$Cov))
  expect_gte(lmfit$omnibus["p.value"], 0)
  expect_lte(lmfit$omnibus["p.value"], 1)
})

test_that("Fine-Gray Lunn-McNeil runs and returns omega+ estimates", {
  skip_on_cran()
  d <- make_test_data(n = 200)
  lmfit <- lunnmcneil(d$Time, d$cause, d$Cov, type = "finegray")
  expect_s3_class(lmfit, "gce_lunnmcneil")
  expect_equal(nrow(lmfit$omegaplus), ncol(d$Cov))
  expect_true(all(is.finite(lmfit$omegaplus$log_ratio)))
  expect_true(all(lmfit$omegaplus$p.value >= 0 & lmfit$omegaplus$p.value <= 1))
})

test_that("single-covariate models do not error", {
  d <- make_test_data()
  expect_s3_class(lunnmcneil(d$Time, d$cause, d$Cov[, "age", drop = FALSE], type = "coxph"),
                  "gce_lunnmcneil")
})

test_that("invalid cause codes are rejected", {
  d <- make_test_data()
  bad <- d$cause; bad[1] <- 5L
  expect_error(lunnmcneil(d$Time, bad, d$Cov), "codes other than")
})
