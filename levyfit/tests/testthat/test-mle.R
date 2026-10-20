plan <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)

test_that("wiener_mle gives the values of the report", {
  set.seed(42)
  x <- c(0, cumsum(0.5 * diff(plan) + 1.5 * sqrt(diff(plan)) * rnorm(10)))
  f <- wiener_mle(x, plan)
  expect_named(f, c("mu", "sigma2", "s2"))
  expect_equal(f[["mu"]], 1.064698, tolerance = 1e-6)
  expect_equal(f[["sigma2"]], 1.170411^2, tolerance = 1e-6)
  expect_equal(f[["s2"]], 1.170411^2 * 10 / 9, tolerance = 1e-6)
})

test_that("gamma_mle gives the values of the report", {
  set.seed(42)
  x <- c(0, cumsum(rgamma(10, shape = 2 * diff(plan), rate = 3)))
  f <- gamma_mle(x, plan)
  expect_named(f, c("a", "b"))
  expect_equal(f[["a"]], 1.781033, tolerance = 1e-6)
  expect_equal(f[["b"]], 3.147942, tolerance = 1e-6)
})

test_that("ig_mle gives the values of the report", {
  rig <- function(n, mean, shape) {
    nu <- rnorm(n)^2
    x <- mean + mean^2 * nu / (2 * shape) - mean / (2 * shape) * sqrt(4 * mean * shape * nu + mean^2 * nu^2)
    ifelse(runif(n) <= mean / (mean + x), x, mean^2 / x)
  }
  set.seed(42)
  t <- c(0, cumsum(runif(1000, 0.5, 3.5)))
  x <- c(0, cumsum(rig(1000, 0.5 * diff(t), 2 * diff(t)^2)))
  f <- ig_mle(x, t)
  expect_named(f, c("mu", "eta", "eta_unbiased"))
  expect_equal(f[["mu"]], 0.5082, tolerance = 1e-4)
  expect_equal(f[["eta"]], 2.0042, tolerance = 1e-4)
  expect_equal(f[["eta_unbiased"]], f[["eta"]] * 997 / 1000)
  expect_true(is.na(ig_mle(c(0, 1, 3, 4), c(0, 1, 2, 3))[["eta_unbiased"]]))
})

test_that("invalid data are refused", {
  expect_warning(wiener_mle(c(0, 1), c(0, 1)), "does not exist")
  expect_error(gamma_mle(c(0, 1, 1, 2), c(0, 1, 2, 3)), "strictly positive")
  expect_warning(gamma_mle(c(0, 1, 2, 3), c(0, 1, 2, 3)), "does not exist")
  expect_error(ig_mle(c(0, 1, 2), c(0, 1, 1)), "strictly increasing")
})
