plan <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
h <- 2.25

reconstruct <- function(rsim, bridge, N = 400) {
  as.vector(vapply(seq_len(N), function(k) {
    x <- c(0, cumsum(rsim(diff(plan))))
    incr(bridge(x, plan))$dx
  }, numeric(length(plan) - 1)))
}

test_that("gamma bridge: increments follow Gamma(a teq, b)", {
  set.seed(1)
  d <- reconstruct(function(dt) rgamma(length(dt), 2 * dt, 3), function(x, t) gamma_resampling(data.frame(t, x), theta = c(a = 2, b = 3)))
  expect_gt(ks.test(d, "pgamma", shape = 2 * h, rate = 3)$p.value, 0.001)
})

test_that("Brownian bridge: increments follow N(mu teq, sigma^2 teq)", {
  set.seed(2)
  d <- reconstruct(function(dt) 0.5 * dt + 1.5 * sqrt(dt) * rnorm(length(dt)), function(x, t) wiener_resampling(data.frame(t, x), theta = c(mu = 0.5, sigma = 1.5)))
  expect_gt(ks.test(d, "pnorm", mean = 0.5 * h, sd = 1.5 * sqrt(h))$p.value, 0.001)
})

test_that("inverse Gaussian bridge: one step and increments follow the exact distributions", {
  skip_if_not_installed("statmod")
  set.seed(3)
  d <- reconstruct(function(dt) statmod::rinvgauss(length(dt), mean = dt, shape = 3 * dt^2), function(x, t) ig_resampling(data.frame(t, x), theta = c(mu = 1, eta = 3)))
  expect_gt(ks.test(d, statmod::pinvgauss, mean = h, shape = 3 * h^2)$p.value, 0.001)

  pbridge <- function(v, t1, t2, D, eta) {
    T <- t1 + t2
    s <- sqrt(v * D * (D - v))
    pnorm(sqrt(eta) * (T * v - t1 * D) / s) +
      (1 - 2 * t1 / T) * exp(2 * eta * t1 * t2 / D) * pnorm(sqrt(eta) * ((t1 - t2) * v - t1 * D) / s)
  }
  v <- replicate(3000, ig_resampling(data.frame(c(0, 2, 3), c(0, 1, 1.5)), theta = c(mu = 1, eta = 2))[[1]]$X[2])
  expect_gt(ks.test(v, pbridge, t1 = 1.5, t2 = 0.5, D = 1, eta = 2)$p.value, 0.001)
})

test_that("variance-gamma bridge: increments follow the distribution of X_teq", {
  set.seed(4)
  rvg <- function(dt) { g <- rgamma(length(dt), shape = dt / 2, scale = 2); 0.5 * g + 1.5 * sqrt(g) * rnorm(length(dt)) }
  d <- reconstruct(rvg, function(x, t) vg_resampling(data.frame(t, x), theta = c(theta = 0.5, sigma = 1.5, nu = 2)))
  expect_gt(suppressWarnings(ks.test(d, rvg(rep(h, 2e5)))$p.value), 0.001)
})

test_that("Meixner bridge: exact sampler and increments", {
  set.seed(5)
  for (cfg in list(c(0.5, 1, -2), c(0.01, 0.02, 0.3), c(1.5, 0.4, 1.8))) {
    a1 <- cfg[1]; a2 <- cfg[2]; z <- cfg[3]
    lZ <- log(2 * pi) + lbeta(2 * a1, 2 * a2) + lh(a1 + a2, z)
    dq <- function(u) exp(lh(a1, u) + lh(a2, z - u) - lZ)
    y <- replicate(2000, rmeixner_bridge_y(a1, a2, z))
    pts <- sort(unique(c(y, 0, z, -a1, a1, z - a2, z + a2)))
    seg <- vapply(seq_len(length(pts) - 1), function(k) integrate(dq, pts[k], pts[k + 1], rel.tol = 1e-10)$value, 0)
    F <- integrate(dq, -Inf, pts[1])$value + c(0, cumsum(seg))
    expect_equal(F[length(F)] + integrate(dq, pts[length(pts)], Inf)$value, 1, tolerance = 1e-6)
    expect_gt(ks.test(y, approxfun(pts, F))$p.value, 0.001)
  }

  skip_if_not_installed("Runuran")
  gens <- lapply(unique(diff(plan)), function(s) Runuran::pinvd.new(Runuran::udmeixner(1, 0.3, 1.5 * s, 0.2 * s)))
  names(gens) <- format(unique(diff(plan)))
  rmx <- function(dt) vapply(format(dt), function(k) Runuran::ur(gens[[k]], 1), 0)
  d <- reconstruct(rmx, function(x, t) meixner_resampling(data.frame(t, x), theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2)), N = 200)
  g <- Runuran::pinvd.new(Runuran::udmeixner(1, 0.3, 1.5 * h, 0.2 * h))
  expect_gt(ks.test(d, function(q) Runuran::up(g, q))$p.value, 0.001)
})
