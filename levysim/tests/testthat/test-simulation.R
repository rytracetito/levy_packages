tg <- seq(0, 10, by = 0.5)
L15 <- function(t) t^1.5
last_values <- function(x) x[nrow(x), -1]

test_that("sim_gamma reproduces the values of the report", {
  set.seed(1)
  v <- last_values(sim_gamma(tg, a = 2, b = 3, A = L15, n_traj = 5000))
  expect_equal(mean(v), 21.061, tolerance = 1e-4)
  expect_equal(var(v), 6.980, tolerance = 1e-3)
})

test_that("the moments at t = 10 match the theory", {
  set.seed(2)
  v <- last_values(sim_wiener(tg, mu = 0.5, sigma = 1.5, Lambda = L15, tau = function(t) t^1.2, n_traj = 4000))
  expect_equal(mean(v), 0.5 * L15(10), tolerance = 0.02)
  expect_equal(var(v), 2.25 * 10^1.2, tolerance = 0.06)
  v <- last_values(sim_ig(tg, mu = 0.5, eta = 2, Lambda = L15, n_traj = 4000))
  expect_equal(mean(v), 0.5 * L15(10), tolerance = 0.02)
  expect_equal(var(v), 0.5^3 * L15(10) / 2, tolerance = 0.06)
  v <- last_values(sim_vg(tg, theta = 0.5, sigma = 1.2, nu = 0.8, Lambda = L15, n_traj = 4000))
  expect_equal(mean(v), 0.5 * L15(10), tolerance = 0.03)
  expect_equal(var(v), (1.44 + 0.25 * 0.8) * L15(10), tolerance = 0.06)
  v <- last_values(sim_meixner(tg, a = 1, b = 0.3, d = 1.5, m = 0.2, Lambda = L15, n_traj = 2000))
  expect_equal(mean(v), (0.2 + 1.5 * tan(0.15)) * L15(10), tolerance = 0.03)
  expect_equal(var(v), 1.5 / (2 * cos(0.15)^2) * L15(10), tolerance = 0.08)
})

test_that("the two VG methods give the same distribution", {
  set.seed(6)
  v1 <- last_values(sim_vg(tg, theta = 0.5, sigma = 1.2, nu = 0.8, Lambda = L15, n_traj = 4000))
  v2 <- last_values(sim_vg(tg, theta = 0.5, sigma = 1.2, nu = 0.8, Lambda = L15,
                           method = "difference", n_traj = 4000))
  expect_equal(mean(v2), 0.5 * L15(10), tolerance = 0.03)
  expect_equal(var(v2), (1.44 + 0.25 * 0.8) * L15(10), tolerance = 0.06)
  expect_gt(ks.test(v1, v2)$p.value, 0.01)
  expect_error(sim_vg(tg, 0.5, 1.2, 0.8, method = "autre"), "method")
})

test_that("sim_igam follows the marginal distribution P(W(t) <= w) = Q(eta(w), t / theta)", {
  skip_on_cran()
  set.seed(3)
  t <- c(0, 1.2, 3.8)
  W <- sim_igam(t, theta = 1, alpha = 1, beta = 1.5, n_traj = 500)[, -1]
  u <- pgamma(t[3], shape = W[3, ]^1.5, lower.tail = FALSE)
  expect_gt(ks.test(u, "punif")$p.value, 0.01)
})

test_that("the processes have the expected shape", {
  set.seed(4)
  expect_true(all(diff(sim_gamma(tg, 2, 3)[, "X1"]) > 0))
  expect_true(all(diff(sim_ig(tg, 0.5, 2)[, "X1"]) > 0))
  expect_true(all(diff(sim_igam(c(0, 1, 2, 3), 1, 1, 1.5)[, "X1"]) >= 0))
  expect_true(any(diff(sim_wiener(tg, 0.5, 1.5)[, "X1"]) < 0))
})

test_that("several trajectories and inspection plans are handled", {
  set.seed(5)
  x <- sim_gamma(tg, 2, 3, n_traj = 4)
  expect_s3_class(x, "levysim")
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

test_that("invalid inputs are refused", {
  expect_error(sim_gamma(c(1, 2, 3), 2, 3), "start at 0")
  expect_error(sim_gamma(c(0, 2, 1), 2, 3), "strictly increasing")
  expect_error(sim_gamma(tg, -1, 3), "positive")
  expect_error(sim_wiener(tg, 0.5, 1.5, tau = function(t) -t), "tau")
  expect_error(sim_ig(tg, 0.5, 2, Lambda = function(t) 0 * t), "Lambda")
  expect_error(sim_meixner(tg, 1, 4, 1.5, 0.2), "pi")
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
