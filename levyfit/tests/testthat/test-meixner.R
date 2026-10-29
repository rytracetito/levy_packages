test_that("the complex trigamma matches trigamma and the derivative of the complex digamma", {
  x <- c(0.05, 0.3, 1, 2.7, 15, 80)
  expect_equal(Re(complex_trigamma(complex(real = x))), trigamma(x), tolerance = 1e-12)
  z <- complex(real = c(0.05, 0.4, 1.5, 3, 12), imaginary = c(-3, 0.7, -25, 4, 60))
  h <- 1e-5
  num <- (complex_digamma(z + h) - complex_digamma(z - h)) / (2 * h)
  expect_lt(max(Mod(complex_trigamma(z) - num) / Mod(num)), 1e-7)
})

test_that("the explicit Hessian of the Meixner log-likelihood matches the numerical one", {
  skip_if_not_installed("Runuran")
  set.seed(5)
  dt <- runif(200, 0.3, 2)
  dx <- vapply(dt, function(s) Runuran::ur(Runuran::pinvd.new(
    Runuran::udmeixner(alpha = 1, beta = 0.3, delta = 1.5 * s, mu = 0.2 * s)), 1), 0)
  f <- function(p) meixner_loglik(dx, dt, p[1], p[2], p[3], p[4])
  for (p in list(c(1, 0.3, 1.5, 0.2), c(0.7, -0.5, 2.5, 0.4))) {
    ha <- meixner_hessian(dx, dt, p[1], p[2], p[3], p[4])
    expect_equal(unname(ha), optimHess(p, f), tolerance = 1e-4)
  }
})

test_that("the moment estimates of the Meixner process are consistent", {
  skip_if_not_installed("Runuran")
  set.seed(6)
  gen <- Runuran::pinvd.new(Runuran::udmeixner(alpha = 1, beta = 0.3, delta = 1.5, mu = 0.2))
  dx <- Runuran::ur(gen, 50000)
  m <- meixner_moments(dx, rep(1, 50000))
  expect_equal(unname(m), c(1, 0.3, 1.5, 0.2), tolerance = 0.15)
  expect_null(meixner_moments(c(1, 2, 3, 4, 5), rep(1, 5)))
})

test_that("the Meixner estimate is the maximum reached from the moment estimates", {
  skip_if_not_installed("Runuran")
  set.seed(1)
  gen <- Runuran::pinvd.new(Runuran::udmeixner(alpha = 1, beta = 0.3, delta = 1.5, mu = 0.2))
  dx <- Runuran::ur(gen, 500)
  x <- c(0, cumsum(dx))
  t <- 0:500
  e <- meixner_mle(x, t)
  s <- meixner_score(dx, rep(1, 500), e[["a"]], e[["b"]], e[["d"]], e[["m"]])
  expect_lt(max(abs(s)), 1e-3)
  h <- meixner_hessian(dx, rep(1, 500), e[["a"]], e[["b"]], e[["d"]], e[["m"]])
  expect_lt(max(eigen(h, symmetric = TRUE, only.values = TRUE)$values), 0)
  ll <- function(p) meixner_loglik(dx, rep(1, 500), p[1], p[2], p[3], p[4])
  expect_gt(ll(e), ll(c(1, 0.3, 1.5, 0.2)))
  ci <- meixner_mle(x, t, ci = TRUE)
  expect_true(all(ci$lower < ci$estimate & ci$estimate < ci$upper))
})
