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

test_that("the inverse gamma likelihood has total mass 1", {
  skip_on_cran()
  z <- c(0, 1, 2.5)
  distinct <- integrate(function(u1) sapply(u1, function(v1) integrate(function(u2) sapply(u2, function(v2)
    exp(igam_loglik_red(z, c(0, v1, v2)))), v1, v1 + 25, rel.tol = 1e-7)$value), 1e-8, 20, rel.tol = 1e-7)$value
  plateau <- integrate(function(u1) sapply(u1, function(v) exp(igam_loglik_red(z, c(0, v, v)))),
                       1e-8, 20, rel.tol = 1e-9)$value
  expect_equal(distinct + plateau, 1, tolerance = 1e-4)
})

test_that("igam_mle gives the values of the original code", {
  skip_on_cran()
  t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
  x <- c(0, 2.212327259, 2.81893462, 3.413725025, 4.801409076, 6.228550224,
         7.714096072, 7.949339032, 8.524511739, 9.238423227, 9.359099294)
  f <- igam_mle(x, t)
  expect_named(f, c("theta", "alpha", "beta"))
  expect_equal(unname(f), c(1.048168561, 2.136460775, 2.094853401), tolerance = 1e-6)
})

test_that("invalid data are refused", {
  expect_error(vg_mle(c(0, 1, 1, 2), c(0, 1, 2, 3)), "zero increment")
  expect_error(meixner_mle(c(0, 1, 2), c(0, 1, 2)), "four increments")
  expect_error(igam_mle(c(1, 2, 3), c(0, 1, 2)), "start at")
  expect_error(igam_mle(c(0, 2, 1), c(0, 1, 2)), "non-decreasing")
})

test_that("the VG log-likelihood tends to the gamma one when sigma tends to 0", {
  dx <- c(0.0023, 0.41, 1.3, 0.07)
  dt <- c(0.5, 1, 2.5, 0.4)
  theta <- 1.2
  nu <- 0.3
  lg <- sum(dgamma(dx, shape = dt / nu, scale = theta * nu, log = TRUE))
  for (s in c(1e-6, 1e-12, 1e-50, 1e-100)) expect_equal(vg_loglik(dx, dt, theta, s, nu), lg, tolerance = 1e-8)
})
