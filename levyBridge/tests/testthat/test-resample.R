t_obs <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
x_obs <- c(0, 1.31, 2.05, 3.62, 5.71, 7.15, 10.62, 12.83, 15.02, 18.94, 22.31)

test_that("the grid is regular, with the end points and the observations kept", {
  set.seed(1)
  r <- gamma_resampling(data.frame(t_obs, x_obs), theta = c(a = 2, b = 3))
  expect_s3_class(r, "levybridge")
  expect_equal(r[[1]]$t, seq(0, 22.5, length.out = 11))
  expect_equal(attr(r, "teq"), 2.25)
  expect_equal(r[[1]]$X[c(1, 11)], x_obs[c(1, 11)])
  expect_equal(attr(r, "observed")[[1]], data.frame(t = t_obs, X = x_obs))

  t2 <- c(0, 2, 3, 6, 7, 10)
  x2 <- c(0, 1.4, 2.0, 4.1, 4.6, 6.8)
  for (r in list(gamma_resampling(data.frame(t2, x2), theta = c(a = 2, b = 3)), wiener_resampling(data.frame(t2, x2), theta = c(mu = 0.5, sigma = 1)), ig_resampling(data.frame(t2, x2), theta = c(mu = 1, eta = 3)),
                 vg_resampling(data.frame(t2, x2), theta = c(theta = 0.5, sigma = 1.5, nu = 2)), meixner_resampling(data.frame(t2, x2), theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2)))) {
    expect_equal(r[[1]]$X[c(2, 4)], c(1.4, 4.1))
  }
})

test_that("a grid point equal to an inspection time up to rounding keeps the observed value", {
  t <- c(0, 0.3, 0.6, 0.9 + 1e-15, 1.2, 1.5, 1.8)
  x <- c(0, 0.2, 0.5, 0.61, 1.1, 1.2, 1.8)
  set.seed(1)
  r <- wiener_resampling(data.frame(t, x), theta = c(mu = 0.5, sigma = 1.5))
  expect_identical(r[[1]]$t[4], t[4])
  expect_identical(r[[1]]$X[4], x[4])
})

test_that("the resampled values stay between the known values for increasing processes", {
  set.seed(2)
  for (k in 1:20) {
    y <- gamma_resampling(data.frame(t_obs, x_obs), theta = c(a = 0.3, b = 3))[[1]]$X
    expect_true(all(diff(y) >= 0))
    y <- ig_resampling(data.frame(t_obs, x_obs), theta = c(mu = 1, eta = 0.5))[[1]]$X
    expect_true(all(diff(y) >= 0))
  }
})

test_that("equal observed values give a constant bridge, a decrease is refused", {
  expect_equal(gamma_resampling(data.frame(c(0, 1, 2.5, 3), c(0, 1, 1, 2)), theta = c(a = 2, b = 3))[[1]]$X[3], 1)
  expect_equal(ig_resampling(data.frame(c(0, 1, 2.5, 3), c(0, 1, 1, 2)), theta = c(mu = 1, eta = 3))[[1]]$X[3], 1)
  expect_error(gamma_resampling(data.frame(c(0, 1, 2, 3), c(0, 2, 1.5, 3)), theta = c(a = 2, b = 3)), "must not decrease")
  expect_error(ig_resampling(data.frame(c(0, 1, 2, 3), c(0, 2, 1.5, 3)), theta = c(mu = 1, eta = 3)), "must not decrease")
  expect_silent(wiener_resampling(data.frame(c(0, 1, 2, 3), c(0, 2, 1.5, 3)), theta = c(mu = 0.5, sigma = 1)))
})

test_that("invalid data, parameters and methods are refused with a clear message", {
  expect_error(gamma_resampling(data.frame(t_obs + 1, x_obs), theta = c(a = 2, b = 3)), "first inspection time must be 0")
  expect_error(gamma_resampling(data.frame(t_obs, replace(x_obs, 3, NA)), theta = c(a = 2, b = 3)), "missing")
  expect_error(gamma_resampling(data.frame(rev(t_obs), x_obs), theta = c(a = 2, b = 3)), "first inspection time|strictly increasing")
  expect_error(gamma_resampling(data.frame(0, 0), theta = c(a = 2, b = 3)), "at least two")
  expect_error(gamma_resampling(x_obs, theta = c(a = 2, b = 3)), "x must be a data frame")
  expect_error(gamma_resampling(data.frame(t = t_obs), theta = c(a = 2, b = 3)), "column 2")
  expect_error(gamma_resampling(data.frame(t_obs, x_obs), theta = c(a = c(1, 2), b = 3)), "must contain the parameter")
  expect_error(gamma_resampling(data.frame(t_obs, x_obs), theta = c(a = NA, b = 3)), "finite number")
  expect_error(gamma_resampling(data.frame(t_obs, x_obs), theta = c(a = -1, b = 3)), "positive")
  expect_error(vg_resampling(data.frame(t_obs, x_obs), theta = c(theta = "a", sigma = 1, nu = 1)), "named numeric vector")
  expect_silent(vg_resampling(data.frame(t_obs, x_obs), theta = c(theta = -1, sigma = 1, nu = 1)))
  expect_error(gamma_resampling(data.frame(t_obs, x_obs), theta = c(a = 2, b = 3), method = "other"))
})

test_that("several trajectories: one grid and one time step per trajectory", {
  liste <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.8, 2.1, 4.3)),
                data.frame(t = c(0, 2, 4, 5, 9), X = c(0, 1.2, 2.9, 3.1, 6.2)))
  set.seed(3)
  r <- gamma_resampling(liste, theta = c(a = 2, b = 3))
  expect_length(r, 2)
  expect_equal(attr(r, "teq"), c(2, 2.25))
  expect_equal(r[[2]]$t, seq(0, 9, by = 2.25))
  d <- incr(r)
  expect_equal(names(d), c("id", "teq", "dx"))
  expect_equal(as.vector(table(d$id)), c(3, 4))
  expect_equal(sum(d$dx[d$id == 2]), 6.2)
  expect_equal(nrow(gamma_resampling(liste[[1]], theta = c(a = 2, b = 3))[[1]]), 4)
  expect_error(gamma_resampling(list(data.frame(t = c(1, 2), X = c(0, 1))), theta = c(a = 2, b = 3)), "trajectory 1")
})

test_that("print, increments and plot_bridge work", {
  set.seed(4)
  r <- meixner_resampling(data.frame(t_obs, x_obs), theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2))
  expect_output(print(r), "Meixner bridge")
  expect_s3_class(plot_bridge(r), "ggplot")
  r2 <- wiener_resampling(list(data.frame(t = 0:3, X = c(0, 1, 0.5, 2)), data.frame(t = 0:2, X = c(0, 1, 3))), theta = c(mu = 0.5, sigma = 1))
  expect_s3_class(plot_bridge(r2, trajectory = 1:2), "ggplot")
  expect_error(plot_bridge(r2, trajectory = 3), "between 1 and 2")
})
