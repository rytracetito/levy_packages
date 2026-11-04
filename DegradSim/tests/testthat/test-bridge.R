hb <- 5
tb <- seq(0, 100, by = hb)
sims <- list(
  gamma = list(f = function(...) sim_gamma(..., a = 2, b = 3), p = function(q) pgamma(q, 2 * hb, 3)),
  wiener = list(f = function(...) sim_wiener(..., mu = 1.5, sigma = 2), p = function(q) pnorm(q, 1.5 * hb, 2 * sqrt(hb))),
  ig = list(f = function(...) sim_ig(..., mu = 0.5, eta = 2), p = NULL),
  vg = list(f = function(...) sim_vg(..., theta = 0.5, sigma = 1.5, nu = 2), p = NULL),
  meixner = list(f = function(...) sim_meixner(..., a = 1, b = 0.3, d = 1.5, m = 0.2), p = NULL)
)

test_that("Bridge returns the same forms of result as RWA", {
  for (s in sims) {
    x <- s$f(tb, method = "Bridge", n_traj = 3)
    expect_s3_class(x, "DegradSim")
    expect_equal(dim(x), c(length(tb), 4L))
    expect_equal(unname(unclass(x)[1, -1]), rep(0, 3))
    expect_equal(attr(x, "model"), attr(s$f(tb, n_traj = 1), "model"))
    l <- s$f(list(c(0, 1, 3), c(0, 2, 4, 6)), method = "Bridge")
    expect_equal(vapply(l, nrow, 0L), c(3L, 4L))
    plan <- inspection_times("decreasing", n = 5:6, n_traj = 2, tmax = 50)
    expect_equal(vapply(s$f(plan, method = "Bridge"), nrow, 0L), c(6L, 7L))
    expect_equal(dim(s$f(c(0, 1), method = "Bridge")), c(2L, 2L))
  }
})

test_that("Bridge gives the distribution of the increments", {
  set.seed(21)
  x <- unclass(sims$gamma$f(tb, method = "Bridge", n_traj = 300))[, -1]
  expect_gt(suppressWarnings(ks.test(as.vector(diff(x)), sims$gamma$p))$p.value, 0.001)
  x <- unclass(sims$wiener$f(tb, method = "Bridge", n_traj = 300))[, -1]
  expect_gt(suppressWarnings(ks.test(as.vector(diff(x)), sims$wiener$p))$p.value, 0.001)
  for (n in c("ig", "vg", "meixner")) {
    a <- as.vector(diff(unclass(sims[[n]]$f(tb, method = "Bridge", n_traj = 150))[, -1]))
    b <- as.vector(diff(unclass(sims[[n]]$f(tb, n_traj = 150))[, -1]))
    expect_gt(suppressWarnings(ks.test(a, b))$p.value, 0.001)
  }
})

test_that("Bridge with a time-scaling function", {
  A <- function(t) 0.02 * t^1.5
  set.seed(22)
  dA <- diff(A(tb))
  x <- unclass(sim_gamma(tb, b = 3, A = A, method = "Bridge", n_traj = 300))[, -1]
  j <- length(dA)
  expect_gt(suppressWarnings(ks.test(diff(x)[j, ], "pgamma", dA[j], 3))$p.value, 0.001)
  x <- unclass(sim_wiener(tb, sigma = 2, A = A, method = "Bridge", n_traj = 300))[, -1]
  expect_gt(ks.test(diff(x)[j, ], "pnorm", dA[j], 2 * sqrt(dA[j]))$p.value, 0.001)
  expect_equal(attr(sim_ig(tb, eta = 2, A = A, method = "Bridge"), "model"), "time-scaled inverse Gaussian")
})

test_that("Bridge keeps the level where A is constant", {
  A <- function(t) pmin(t, 2)
  t0 <- c(0, 1, 2, 3, 4, 5)
  for (x in list(sim_gamma(t0, b = 3, A = A, method = "Bridge"),
                 sim_wiener(t0, sigma = 2, A = A, method = "Bridge"),
                 sim_ig(t0, eta = 2, A = A, method = "Bridge"))) {
    v <- unclass(x)[, "X1"]
    expect_equal(unname(diff(v)[3:5]), c(0, 0, 0))
  }
  expect_error(sim_gamma(c(0, 1, 2), b = 3, A = function(t) -t, method = "Bridge"), "non-decreasing")
})

test_that("Bridge is reproducible with a seed", {
  for (s in sims) {
    set.seed(23); a <- s$f(tb, method = "Bridge")
    set.seed(23); b <- s$f(tb, method = "Bridge")
    expect_identical(a, b)
  }
})
