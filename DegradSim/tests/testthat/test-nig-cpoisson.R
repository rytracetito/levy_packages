tg <- seq(0, 10, by = 0.5)
last_values <- function(x) x[nrow(x), -1]

test_that("sim_nig has the moments of the NIG distribution at t = 10", {
  set.seed(21)
  al <- 2; be <- 0.5; de <- 1; mu <- 0.2
  ga <- sqrt(al^2 - be^2)
  v <- last_values(sim_nig(tg, al, be, de, mu, n_traj = 4000))
  expect_equal(mean(v), (mu + de * be / ga) * 10, tolerance = 0.03)
  expect_equal(var(v), de * al^2 * 10 / ga^3, tolerance = 0.06)
})

test_that("the NIG increments follow the NIG density", {
  set.seed(22)
  al <- 2; be <- 0.5; de <- 1; mu <- 0.2; h <- 0.7
  dnig <- function(y) {
    q <- sqrt((de * h)^2 + (y - mu * h)^2)
    al * de * h / (pi * q) * besselK(al * q, 1, expon.scaled = TRUE) *
      exp(de * h * sqrt(al^2 - be^2) + be * (y - mu * h) - al * q)
  }
  expect_equal(integrate(dnig, -Inf, Inf)$value, 1, tolerance = 1e-6)
  x <-sim_nig(c(0, h), al, be, de, mu, n_traj = 3000)[2, -1]
  pnig <- function(q) vapply(q, function(v) integrate(dnig, -Inf, v)$value, 0)
  expect_gt(ks.test(x, pnig)$p.value, 0.01)
})

test_that("sim_cpoisson has the moments of the compound Poisson process", {
  set.seed(23)
  v <- last_values(sim_cpoisson(tg, lambda = 2, rjump = function(n) rgamma(n, 2, 3), n_traj = 4000))
  EJ <- 2 / 3
  EJ2 <- 2 / 9 + EJ^2
  expect_equal(mean(v), 2 * 10 * EJ, tolerance = 0.03)
  expect_equal(var(v), 2 * 10 * EJ2, tolerance = 0.06)
})

test_that("sim_cpoisson gives zero increments with probability exp(-lambda h)", {
  set.seed(24)
  x <- sim_cpoisson(seq(0, 200, by = 1), lambda = 1, n_traj = 10)
  dx <- as.vector(apply(unclass(x)[, -1], 2, diff))
  expect_equal(mean(dx == 0), exp(-1), tolerance = 0.05)
  expect_true(all(dx >= 0))
})

test_that("jumps equal to 1 give the Poisson process", {
  set.seed(25)
  x <- sim_cpoisson(c(0, 3), lambda = 2, rjump = function(n) rep(1, n), n_traj = 4000)
  v <- x[2, -1]
  expect_true(all(v == round(v)))
  expect_equal(mean(v), 6, tolerance = 0.03)
  expect_equal(var(v), 6, tolerance = 0.06)
})

test_that("the jump law is stored and printed", {
  set.seed(26)
  x <- sim_cpoisson(tg, lambda = 2, rjump = function(n) rgamma(n, 2, 3))
  expect_output(print(x), "jumps: rjump =")
  expect_equal(attr(x, "parameters"), c(lambda = 2))
  y <- sim_cpoisson(list(c(0, 1, 2), c(0, 2, 5, 7)), lambda = 2)
  expect_false(is.matrix(y))
  expect_s3_class(plot_traj(y), "ggplot")
})

test_that("invalid NIG and compound Poisson parameters are refused", {
  expect_error(sim_nig(tg, 1, 2, 1, 0), "beta")
  expect_error(sim_nig(tg, -1, 0, 1, 0), "alpha")
  expect_error(sim_nig(tg, 2, 0, 0, 0), "delta")
  expect_error(sim_cpoisson(tg, lambda = 0), "lambda")
  expect_error(sim_cpoisson(tg, lambda = 1, rjump = 2), "rjump")
  expect_error(sim_cpoisson(tg, lambda = 5, rjump = function(n) 1), "rjump")
})
