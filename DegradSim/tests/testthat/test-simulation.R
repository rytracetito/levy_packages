tg <- seq(0, 10, by = 0.5)
A15 <- function(t) 0.5 * t^1.5
last_values <- function(x) x[nrow(x), -1]

test_that("sim_gamma with A(t) = 2 t^1.5 reproduces the values of the report", {
  set.seed(1)
  v <- last_values(sim_gamma(tg, b = 3, A = function(t) 2 * t^1.5, n_traj = 5000))
  expect_equal(mean(v), 21.061, tolerance = 1e-4)
  expect_equal(var(v), 6.980, tolerance = 1e-3)
})

test_that("the moments at t = 10 match the theory, homogeneous processes", {
  set.seed(2)
  v <- last_values(sim_gamma(tg, a = 2, b = 3, n_traj = 4000))
  expect_equal(mean(v), 2 * 10 / 3, tolerance = 0.02)
  expect_equal(var(v), 2 * 10 / 9, tolerance = 0.06)
  v <- last_values(sim_wiener(tg, mu = 0.5, sigma = 1.5, n_traj = 4000))
  expect_equal(mean(v), 5, tolerance = 0.03)
  expect_equal(var(v), 2.25 * 10, tolerance = 0.06)
  v <- last_values(sim_ig(tg, mu = 0.5, eta = 2, n_traj = 4000))
  expect_equal(mean(v), 5, tolerance = 0.02)
  expect_equal(var(v), 0.5^3 * 10 / 2, tolerance = 0.06)
  v <- last_values(sim_vg(tg, theta = 0.5, sigma = 1.2, nu = 0.8, n_traj = 4000))
  expect_equal(mean(v), 5, tolerance = 0.04)
  expect_equal(var(v), (1.44 + 0.25 * 0.8) * 10, tolerance = 0.06)
  v <- last_values(sim_meixner(tg, a = 1, b = 0.3, d = 1.5, m = 0.2, n_traj = 2000))
  expect_equal(mean(v), (0.2 + 1.5 * tan(0.15)) * 10, tolerance = 0.03)
  expect_equal(var(v), 1.5 / (2 * cos(0.15)^2) * 10, tolerance = 0.08)
})

test_that("the moments at t = 10 match the theory, non-homogeneous processes", {
  set.seed(3)
  v <- last_values(sim_gamma(tg, b = 3, A = A15, n_traj = 4000))
  expect_equal(mean(v), A15(10) / 3, tolerance = 0.02)
  expect_equal(var(v), A15(10) / 9, tolerance = 0.06)
  v <- last_values(sim_wiener(tg, sigma = 1.5, A = A15, n_traj = 4000))
  expect_equal(mean(v), A15(10), tolerance = 0.02)
  expect_equal(var(v), 2.25 * A15(10), tolerance = 0.06)
  v <- last_values(sim_ig(tg, eta = 2, A = A15, n_traj = 4000))
  expect_equal(mean(v), A15(10), tolerance = 0.02)
  expect_equal(var(v), A15(10) / 2, tolerance = 0.06)
})

test_that("the increments of the non-homogeneous processes have the right distribution", {
  set.seed(8)
  x <- sim_gamma(c(0, 2, 5), b = 3, A = A15, n_traj = 3000)
  expect_gt(ks.test(x[3, -1] - x[2, -1], "pgamma", shape = A15(5) - A15(2), rate = 3)$p.value, 0.01)
  x <- sim_wiener(c(0, 2, 5), sigma = 1.5, A = A15, n_traj = 3000)
  dA <- A15(5) - A15(2)
  expect_gt(ks.test(x[3, -1] - x[2, -1], "pnorm", dA, 1.5 * sqrt(dA))$p.value, 0.01)
})

test_that("A(t) = a t, mu t gives back the homogeneous processes", {
  for (f in list(
    list(function() sim_gamma(tg, a = 2, b = 3), function() sim_gamma(tg, b = 3, A = function(t) 2 * t)),
    list(function() sim_wiener(tg, mu = 0.5, sigma = 1.5),
         function() sim_wiener(tg, sigma = 1.5 / sqrt(0.5), A = function(t) 0.5 * t)),
    list(function() sim_ig(tg, mu = 0.5, eta = 2), function() sim_ig(tg, eta = 8, A = function(t) 0.5 * t)))) {
    set.seed(9)
    y1 <- unclass(f[[1]]())[, 2]
    set.seed(9)
    y2 <- unclass(f[[2]]())[, 2]
    expect_equal(y1, y2)
  }
})

test_that("a constant time change gives zero increments", {
  set.seed(10)
  A <- function(t) pmin(t, 3)
  for (x in list(sim_gamma(tg, b = 3, A = A), sim_wiener(tg, sigma = 1, A = A), sim_ig(tg, eta = 2, A = A))) {
    expect_true(all(diff(x[tg >= 3, "X1"]) == 0))
  }
})

test_that("the two VG methods give the same distribution", {
  set.seed(6)
  v1 <- last_values(sim_vg(tg, theta = 0.5, sigma = 1.2, nu = 0.8, n_traj = 4000))
  v2 <- last_values(sim_vg(tg, theta = 0.5, sigma = 1.2, nu = 0.8, method = "difference", n_traj = 4000))
  expect_equal(mean(v2), 5, tolerance = 0.04)
  expect_equal(var(v2), (1.44 + 0.25 * 0.8) * 10, tolerance = 0.06)
  expect_gt(ks.test(v1, v2)$p.value, 0.01)
  expect_error(sim_vg(tg, 0.5, 1.2, 0.8, method = "autre"), "method")
})

test_that("the processes have the expected shape", {
  set.seed(4)
  expect_true(all(diff(sim_gamma(tg, 2, 3)[, "X1"]) > 0))
  expect_true(all(diff(sim_ig(tg, 0.5, 2)[, "X1"]) > 0))
  expect_true(any(diff(sim_wiener(tg, 0.5, 1.5)[, "X1"]) < 0))
})

test_that("several trajectories and inspection plans are handled", {
  set.seed(5)
  x <- sim_gamma(tg, 2, 3, n_traj = 4)
  expect_s3_class(x, "DegradSim")
  expect_true(is.matrix(x))
  expect_equal(colnames(x), c("t", "X1", "X2", "X3", "X4"))
  expect_equal(x[, "t"], tg)
  expect_equal(unname(x[1, -1]), rep(0, 4))
  expect_equal(colnames(sim_wiener(tg, 0.5, 1.5)), c("t", "X1"))
  same <- sim_ig(inspection_times("increasing", n = 10, n_traj = 3, tmax = 50), 0.5, 2)
  expect_true(is.matrix(same))
  expect_equal(dim(same), c(11, 4))
  y <- sim_wiener(list(c(0, 1, 2), c(0, 2, 5, 7)), 0.5, 1.5)
  expect_false(is.matrix(y))
  expect_equal(lengths(lapply(y, `[[`, "t")), c(3, 4))
  expect_equal(y[[2]][, 2], y[[2]]$X)
  z <- sim_ig(inspection_times("increasing", n = 10, n_traj = 3, tmax = c(50, 80, 100)), 0.5, 2)
  expect_false(is.matrix(z))
  expect_length(z, 3)
  d <- traj_df(z)
  expect_named(d, c("id", "t", "X"))
  expect_equal(nrow(d), 33)
  expect_s3_class(plot_traj(z), "ggplot")
  expect_s3_class(plot_traj(x), "ggplot")
  expect_equal(nrow(traj_df(x)), 4 * length(tg))
  expect_output(print(x), "matrix")
  expect_output(print(y), "list")
})

test_that("the model and the time change are stored and printed", {
  set.seed(11)
  x <- sim_gamma(tg, b = 3, A = function(t) 2 * t^1.5)
  expect_equal(attr(x, "model"), "non-homogeneous gamma")
  expect_equal(attr(x, "parameters"), c(b = 3))
  expect_output(print(x), "time change: A =")
  expect_equal(plot_traj(x)$labels$title, "non-homogeneous gamma process")
  expect_null(attr(sim_gamma(tg, 2, 3), "A"))
  expect_equal(attr(sim_wiener(tg, sigma = 1, A = A15), "model"), "time-scaled Wiener")
  expect_equal(attr(sim_ig(tg, eta = 1, A = A15), "model"), "non-homogeneous inverse Gaussian")
})

test_that("invalid inspection times are refused", {
  expect_error(sim_gamma(c(1, 2, 3), 2, 3), "start at 0")
  expect_error(sim_gamma(c(0, 2, 1), 2, 3), "strictly increasing")
  expect_error(sim_gamma(c(0, 1, NA), 2, 3), "finite")
  expect_error(sim_gamma(0, 2, 3), "at least two")
  expect_error(sim_gamma(tg, 2, 3, n_traj = 1.5), "n_traj")
})

test_that("invalid parameters are refused", {
  expect_error(sim_gamma(tg, -1, 3), "a must be")
  expect_error(sim_gamma(tg, b = 3), "a is missing")
  expect_error(sim_wiener(tg, sigma = 1), "mu is missing")
  expect_error(sim_ig(tg, eta = 1), "mu is missing")
  expect_error(sim_gamma(tg, 2, -3), "b must be")
  expect_error(sim_wiener(tg, 0.5, 0), "sigma")
  expect_error(sim_ig(tg, -0.5, 2), "mu must be")
  expect_error(sim_meixner(tg, 1, 4, 1.5, 0.2), "pi")
  expect_error(sim_vg(tg, 0.5, 1.2, -1), "nu")
})

test_that("a and A (or mu and A) cannot be given together", {
  expect_error(sim_gamma(tg, a = 2, b = 3, A = A15), "not both")
  expect_error(sim_wiener(tg, mu = 1, sigma = 1, A = A15), "not both")
  expect_error(sim_ig(tg, mu = 1, eta = 1, A = A15), "not both")
})

test_that("invalid time changes are refused", {
  expect_error(sim_gamma(tg, b = 3, A = 2), "R function")
  expect_error(sim_gamma(tg, b = 3, A = function(t) t + 1), "A\\(0\\) = 0")
  expect_error(sim_wiener(tg, sigma = 1, A = function(t) -t), "non-decreasing")
  expect_error(sim_ig(tg, eta = 2, A = function(t) sin(t)), "non-decreasing")
  expect_error(sim_gamma(tg, b = 3, A = function(t) t[1]), "vectorized")
  expect_error(sim_gamma(tg, b = 3, A = function(t) ifelse(t > 9, Inf, t)), "finite")
})

test_that("traj_df gives the time steps and increments", {
  set.seed(7)
  x <- sim_gamma(c(0, 1, 3, 6), 2, 3, n_traj = 2)
  d <- traj_df(x, increments = TRUE)
  expect_named(d, c("id", "t", "X", "dt", "dX"))
  expect_true(all(is.na(d$dt[d$t == 0])) && all(is.na(d$dX[d$t == 0])))
  expect_equal(d$dt[d$id == 1], c(NA, 1, 2, 3))
  expect_equal(d$dX[d$id == 2][-1], unname(diff(x[, "X2"])))
  y <- sim_wiener(list(c(0, 1, 2), c(0, 2, 5, 7)), 0.5, 1.5)
  expect_equal(traj_df(y, increments = TRUE)$dt, c(NA, 1, 1, NA, 2, 3, 2))
})
