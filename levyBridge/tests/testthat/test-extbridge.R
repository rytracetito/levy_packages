test_that("ext_grid: time step, number of decimal places and grid", {
  g <- ext_grid(c(0, 2, 5, 11))
  expect_equal(g$h, 1)
  expect_equal(g$grid, 0:11)

  g <- ext_grid(c(0, 1.2, 2, 3.8, 5.5))
  expect_equal(g$h, 0.1)
  expect_length(g$grid, 56)
  expect_true(all(c(0, 1.2, 2, 3.8, 5.5) %in% g$grid))

  g <- ext_grid(c(0, 0.3, 0.7, 1.1))
  expect_true(0.3 %in% g$grid && 0.7 %in% g$grid)

  expect_equal(ext_grid(c(0, 0.35, 1.2, 1.85, 2.5))$h, 0.05)
  expect_equal(ext_grid(c(0, 2, 4, 6))$grid, c(0, 2, 4, 6))
})

test_that("ext_grid: rounding, warnings and errors", {
  expect_warning(g <- ext_grid(c(0, 1.25, 2.47, 3.86), digits = 1), "rounded to 1 decimal")
  expect_equal(g$t, c(0, 1.2, 2.5, 3.9))
  expect_equal(g$h, 0.1)
  expect_silent(ext_grid(c(0, 1.25, 2.47, 3.86), digits = 3))

  set.seed(1)
  expect_warning(g <- ext_grid(c(0, runif(3, 0, 5))), "rounded to 2 decimal")
  expect_true(all(g$t %in% g$grid))

  expect_error(suppressWarnings(ext_grid(c(0, 1.01, 1.04), digits = 1)), "equal after rounding")
  expect_error(ext_grid(c(0, 1, 2), digits = -1), "non-negative integer")
  expect_warning(ext_grid(c(0, 1.234, 400)), "time steps")
})

test_that("ExtBridge keeps all the observations, for every process", {
  t <- c(0, 1.2, 2, 3.8, 5.5, 7.1)
  x <- c(0, 1.1, 2.3, 3.0, 5.2, 6.4)
  set.seed(2)
  res <- list(gamma_bridge(x, t, a = 2, method = "ExtBridge"),
              wiener_bridge(x, t, sigma = 1.5, method = "ExtBridge"),
              ig_bridge(x, t, eta = 3, method = "ExtBridge"),
              vg_bridge(x, t, theta = 0.5, sigma = 1.5, nu = 2, method = "ExtBridge"),
              meixner_bridge(x, t, a = 1, d = 1.5, m = 0.2, method = "ExtBridge"))
  for (r in res) {
    expect_equal(attr(r, "method"), "ExtBridge")
    expect_equal(attr(r, "teq"), 0.1)
    expect_equal(r[[1]]$t, (0:71) / 10)
    expect_identical(r[[1]]$X[match(t, r[[1]]$t)], x)
    expect_false(anyNA(r[[1]]$X))
  }
  reg <- wiener_bridge(c(0, 1, 0.5, 3), c(0, 2, 4, 6), sigma = 1, method = "ExtBridge")
  expect_identical(reg[[1]]$X, c(0, 1, 0.5, 3))
})

test_that("ExtBridge with rounded times and with several trajectories", {
  set.seed(3)
  expect_warning(r <- gamma_bridge(c(0, 1, 2.2, 3), c(0, 1.25, 2.47, 3.86), a = 2, method = "ExtBridge", digits = 1),
                 "rounded")
  expect_equal(attr(r, "observed")[[1]]$t, c(0, 1.2, 2.5, 3.9))
  expect_equal(r[[1]]$X[match(c(1.2, 2.5, 3.9), r[[1]]$t)], c(1, 2.2, 3))

  liste <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.8, 2.1, 4.3)),
                data.frame(t = c(0, 0.5, 1.25, 2), X = c(0, 0.4, 1.1, 1.6)))
  r <- gamma_bridge(liste, a = 2, method = "ExtBridge")
  expect_equal(attr(r, "teq"), c(1, 0.25))
  expect_equal(sapply(r, nrow), c(7, 9))
  d <- increments(r)
  expect_equal(as.vector(tapply(d$dx, d$id, sum)), c(4.3, 1.6))
})

test_that("ExtBridge: increasing processes stay increasing", {
  set.seed(4)
  t <- c(0, 0.6, 1, 1.9)
  x <- c(0, 0.5, 0.7, 1.6)
  for (k in 1:20) {
    expect_true(all(diff(gamma_bridge(x, t, a = 0.5, method = "ExtBridge")[[1]]$X) >= 0))
    expect_true(all(diff(ig_bridge(x, t, eta = 0.5, method = "ExtBridge")[[1]]$X) >= 0))
  }
})

test_that("ExtBridge: the increments on the grid follow the law of X_h", {
  skip_on_cran()
  plan2 <- c(0, 0.6, 1, 1.9, 2.75)
  h2 <- 0.05
  one_increment <- function(rsim, bridge, pos, R = 1000) {
    vapply(seq_len(R), function(k) {
      x <- c(0, cumsum(rsim(diff(plan2))))
      diff(bridge(x, plan2)[[1]]$X)[pos]
    }, numeric(length(pos)))
  }
  set.seed(5)
  d <- one_increment(function(dt) rgamma(length(dt), 2 * dt, 3),
                     function(x, t) gamma_bridge(x, t, a = 2, method = "ExtBridge"), c(1, 13, 30))
  for (j in 1:3) expect_gt(suppressWarnings(ks.test(d[j, ], "pgamma", shape = 2 * h2, rate = 3))$p.value, 0.001)

  set.seed(6)
  d <- one_increment(function(dt) rnorm(length(dt), 0.5 * dt, 1.5 * sqrt(dt)),
                     function(x, t) wiener_bridge(x, t, sigma = 1.5, method = "ExtBridge"), c(1, 13, 30))
  for (j in 1:3) expect_gt(ks.test(d[j, ], "pnorm", mean = 0.5 * h2, sd = 1.5 * sqrt(h2))$p.value, 0.001)

  skip_if_not_installed("statmod")
  set.seed(7)
  d <- one_increment(function(dt) statmod::rinvgauss(length(dt), mean = 0.5 * dt, shape = 3 * dt^2),
                     function(x, t) ig_bridge(x, t, eta = 3, method = "ExtBridge"), c(1, 13, 30))
  for (j in 1:3) expect_gt(ks.test(d[j, ], statmod::pinvgauss, mean = 0.5 * h2, shape = 3 * h2^2)$p.value, 0.001)
})
