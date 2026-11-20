test_that("the exact increments of attr(x, \"dx\") are used by the bridge", {
  x <- data.frame(t = c(0, 1, 2), X = c(0, 50, 50))
  # the second increment is lost in the levels: the bridge sees 0, the attribute gives 1e-15
  expect_equal(gamma_resampling(x, c(a = 1, b = 1), method = "ExtBridge", digits = 0)[[1]]$X, c(0, 50, 50))
  xs <- list(x)
  attr(xs, "dx") <- list(c(50, 1e-15))
  set.seed(1)
  r <- gamma_resampling(xs, c(a = 1, b = 1), method = "MultiBridge")
  expect_true(all(attr(r, "dx")[[1]] > 0))
  set.seed(1)
  r0 <- gamma_resampling(list(x), c(a = 1, b = 1), method = "MultiBridge")
  expect_equal(attr(r0, "dx")[[1]][2], 0)
})

test_that("an attribute dx that does not match the levels is ignored", {
  x <- list(data.frame(t = c(0, 1, 3), X = c(0, 2, 5)))
  attr(x, "dx") <- list(c(1, 1))
  set.seed(2); a <- gamma_resampling(x, c(a = 2, b = 3))
  attr(x, "dx") <- NULL
  set.seed(2); b <- gamma_resampling(x, c(a = 2, b = 3))
  expect_equal(a[[1]]$X, b[[1]]$X)
})

test_that("the results do not change when the exact increments equal the differences of the levels", {
  set.seed(3)
  x <- DegradSim::sim_gamma(list(c(0, 1.3, 2, 4.4, 6), c(0, 2.5, 3, 7)), a = 2, b = 3)
  y <- x
  attr(y, "dx") <- NULL
  set.seed(4); a <- gamma_resampling(x, c(a = 2, b = 3), method = "MultiBridge")
  set.seed(4); b <- gamma_resampling(y, c(a = 2, b = 3), method = "MultiBridge")
  expect_equal(a[[1]]$X, b[[1]]$X)
  expect_equal(attr(a, "dx"), attr(b, "dx"))
})
