liste_t <- list(c(0, 3, 8), c(0, 2, 5, 9, 12), c(0, 4, 7, 11, 15, 18, 22, 27, 33, 39, 45),
                c(0, 2, 3, 5, 8, 9, 14, 15, 17, 20, 22, 25, 27, 30, 33, 46, 58))
sim_liste <- function(a, b) {
  lapply(liste_t, function(tk) data.frame(t = tk, X = c(0, cumsum(rgamma(length(tk) - 1, a * diff(tk), b)))))
}

test_that("MultiBridge with one trajectory is the method Bridge", {
  t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
  x <- c(0, cumsum(rgamma(10, 2 * diff(t), 3)))
  for (opt in 1:2) {
    set.seed(1); r1 <- gamma_bridge(x, t, a = 2)
    set.seed(1); r2 <- gamma_bridge(x, t, a = 2, method = "MultiBridge", option = opt, b = 3)
    expect_equal(r2[[1]], r1[[1]])
    expect_equal(attr(r2, "teq"), attr(r1, "teq"))
  }
})

test_that("MultiBridge: common step and grids of the two options", {
  set.seed(2)
  liste <- sim_liste(2, 3)
  h <- (8 + 12 + 45 + 58) / (2 + 4 + 10 + 16)
  r1 <- gamma_bridge(liste, a = 2, method = "MultiBridge")
  expect_equal(attr(r1, "teq"), rep(h, 4))
  expect_equal(sapply(r1, nrow), floor(c(8, 12, 45, 58) / h) + 1)
  expect_equal(attr(r1, "method"), "MultiBridge, option 1")
  r2 <- gamma_bridge(liste, a = 2, method = "MultiBridge", option = 2, b = 3)
  expect_equal(sapply(r2, nrow), rep(floor(58 / h) + 1, 4))
  expect_true(all(sapply(r2, function(d) max(d$t)) <= 58))
  # levels: increasing, and equal to the observations at the grid points that are inspection times
  for (k in 1:4) {
    expect_true(all(diff(r2[[k]]$X) >= 0))
    i <- match(round(liste[[k]]$t, 9), round(r1[[k]]$t, 9))
    expect_equal(r1[[k]]$X[i[!is.na(i)]], liste[[k]]$X[!is.na(i)])
  }
  # increments: sum of the increments of a trajectory = level at its last grid point
  d <- increments(r1)
  expect_equal(as.vector(tapply(d$dx, d$id, sum)), sapply(r1, function(p) p$X[nrow(p)]))
})

test_that("MultiBridge: errors and warnings", {
  set.seed(3)
  liste <- sim_liste(2, 3)
  expect_error(gamma_bridge(liste, a = 2, method = "MultiBridge", option = 2), "all the parameters of the process")
  expect_error(gamma_bridge(liste, a = 2, method = "MultiBridge", option = 3), "option must be 1 or 2")
  expect_error(wiener_bridge(liste, sigma = 1, method = "MultiBridge", option = 2), "all the parameters of the process")
  court <- list(data.frame(t = c(0, 0.5), X = c(0, 0.2)), data.frame(t = c(0, 10, 20), X = c(0, 4, 9)))
  expect_warning(r <- gamma_bridge(court, a = 2, method = "MultiBridge"), "shorter than the time step")
  expect_equal(nrow(r[[1]]), 1)
})

test_that("MultiBridge: the pooled reconstructed increments follow the law of X_h", {
  skip_on_cran()
  set.seed(4)
  h <- (8 + 12 + 45 + 58) / 32
  for (opt in 1:2) {
    dx <- unlist(replicate(300, increments(gamma_bridge(sim_liste(2, 3), a = 2, method = "MultiBridge",
                                                        option = opt, b = 3))$dx, simplify = FALSE))
    expect_gt(ks.test(dx, "pgamma", shape = 2 * h, rate = 3)$p.value, 0.001)
  }
  dx <- unlist(replicate(300, {
    w <- lapply(liste_t, function(tk) data.frame(t = tk, X = c(0, cumsum(rnorm(length(tk) - 1, 0.5 * diff(tk), 1.5 * sqrt(diff(tk)))))))
    increments(wiener_bridge(w, sigma = 1.5, method = "MultiBridge", option = 2, mu = 0.5))$dx
  }, simplify = FALSE))
  expect_gt(ks.test(dx, "pnorm", 0.5 * h, 1.5 * sqrt(h))$p.value, 0.001)
})
