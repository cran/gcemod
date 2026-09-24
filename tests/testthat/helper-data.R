# Self-contained synthetic competing-risks data for tests (independent of the
# bundled Sample data, so tests run even before data-raw/make_sample.R is run).
make_test_data <- function(n = 250, seed = 42) {
  set.seed(seed)
  age  <- rnorm(n, 60, 10)
  smoke <- rbinom(n, 1, 0.4)
  bmi  <- rnorm(n, 27, 4)
  lp1  <- -0.4 + 0.4 * smoke + 0.02 * (age - 60)
  lp2  <- -0.7 + 0.03 * (age - 60) + 0.03 * (bmi - 27)
  t1   <- rexp(n, exp(lp1) / 30)
  t2   <- rexp(n, exp(lp2) / 40)
  tc   <- runif(n, 4, 100)
  te   <- pmin(t1, t2)
  obs  <- pmin(te, tc)
  cause <- ifelse(tc < te, 0L, ifelse(t1 <= t2, 1L, 2L))
  list(
    Time  = obs / 12,
    cause = cause,
    Cov   = data.frame(smoke = smoke, age = age, bmi = bmi),
    CA    = as.integer(cause == 1),
    CM    = as.integer(cause == 2),
    All   = as.integer(cause != 0)
  )
}
