test_that("the drawn increments are kept in attr(x, \"dx\")", {
  set.seed(1)
  t <- inspection_times("regular", n = 20, tmax = 10)$time
  sims <- list(sim_gamma(t, a = 2, b = 3, n_traj = 2), sim_ig(t, mu = 1, eta = 3, n_traj = 2),
               sim_wiener(t, mu = 1, sigma = 1.5, n_traj = 2), sim_vg(t, theta = 0.5, sigma = 1.5, nu = 2, n_traj = 2),
               sim_nig(t, alpha = 2, beta = 0.5, delta = 1, mu = 0, n_traj = 2), sim_cpoisson(t, lambda = 2, n_traj = 2),
               sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2, n_traj = 2),
               sim_gamma(list(c(0, 1, 3), c(0, 2, 4, 6)), a = 2, b = 3))
  for (x in sims) {
    dx <- attr(x, "dx")
    expect_length(dx, if (is.matrix(x)) ncol(x) - 1 else length(x))
    for (k in seq_along(dx)) {
      lev <- if (is.matrix(x)) x[, k + 1] else x[[k]]$X
      expect_equal(cumsum(dx[[k]]), lev[-1])
    }
  }
  d <- traj_df(sims[[1]], increments = TRUE)
  expect_identical(d$dX[d$id == 2][-1], attr(sims[[1]], "dx")[[2]])
})

test_that("a tiny increment stays positive in attr(x, \"dx\") although the levels are equal", {
  set.seed(2)
  t <- c(0, 100, 100.05)
  x <- sim_gamma(t, a = 1, b = 1, n_traj = 200)
  lost <- sapply(2:201, function(k) x[3, k] == x[2, k])
  expect_true(any(lost))
  expect_true(all(sapply(attr(x, "dx"), function(v) v[2]) > 0 | !lost))
})

test_that("the method Bridge gives no attribute dx", {
  set.seed(3)
  x <- sim_gamma(c(0, 1, 2, 3), a = 2, b = 3, method = "Bridge")
  expect_null(attr(x, "dx"))
})
