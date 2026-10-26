test_that("vg_mle gives the values of the report", {
  set.seed(42)
  t <- c(0, cumsum(runif(1000, 0.5, 3.5)))
  g <- rgamma(1000, shape = diff(t) / 2, scale = 2)
  x <- c(0, cumsum(0.5 * g + 1.5 * sqrt(g) * rnorm(1000)))
  f <- vg_mle(x, t)
  expect_named(f, c("theta", "sigma", "nu"))
  expect_equal(unname(f), c(0.5223, 1.5422, 1.9162), tolerance = 1e-4)
})

test_that("meixner_mle gives the values of the report", {
  skip_if_not_installed("Runuran")
  set.seed(42)
  t <- c(0, cumsum(runif(1000, 0.5, 3.5)))
  x <- c(0, cumsum(vapply(diff(t), function(s) Runuran::ur(Runuran::pinvd.new(
    Runuran::udmeixner(alpha = 1, beta = 0.3, delta = 1.5 * s, mu = 0.2 * s)), 1), 0)))
  f <- meixner_mle(x, t)
  expect_named(f, c("a", "b", "d", "m"))
  expect_equal(unname(f), c(0.968, 0.387, 1.721, 0.092), tolerance = 1e-3)
})

test_that("invalid data are refused", {
  expect_error(vg_mle(c(0, 1, 1, 2), c(0, 1, 2, 3)), "zero increment")
  expect_error(meixner_mle(c(0, 1, 2), c(0, 1, 2)), "four increments")
})

test_that("the VG log-likelihood tends to the gamma one when sigma tends to 0", {
  dx <- c(0.0023, 0.41, 1.3, 0.07)
  dt <- c(0.5, 1, 2.5, 0.4)
  theta <- 1.2
  nu <- 0.3
  lg <- sum(dgamma(dx, shape = dt / nu, scale = theta * nu, log = TRUE))
  for (s in c(1e-6, 1e-12, 1e-50, 1e-100)) expect_equal(vg_loglik(dx, dt, theta, s, nu), lg, tolerance = 1e-8)
})
