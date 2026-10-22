test_that("the catalogue and the names are checked", {
  expect_equal(names(tests_disponibles()), c("HGP", "HWP", "InvGaus", "VG", "Meixner"))
  expect_equal(tests_disponibles("hgp"), tests_disponibles("HGP"))
  expect_error(tests_disponibles("poisson"), "unknown model")
  expect_error(test_adequation(rexp(10), "HGP", test = "XX"), "unknown test")
  expect_error(test_adequation(rexp(10), "HGP", alpha = 1), "alpha")
  expect_error(test_adequation(rexp(10), "HGP", B = 0), "B must")
  expect_error(test_adequation(c(1, 2), "HGP"), "at least 3")
})

test_that("test_adequation returns one row per test", {
  set.seed(1)
  r <- test_adequation(rnorm(40), "HWP")
  expect_equal(r$test, tests_disponibles("HWP"))
  expect_true(all(r$decision %in% c("reject H0", "do not reject H0")))
  expect_true(all(r$p.value > 0 & r$p.value <= 1))
  r <- test_adequation(rgamma(40, 2), "HGP", test = c("KS", "VGE", "HME1"), B = 30)
  expect_true(all(!is.na(r$p.value) & !is.na(r$critical.value)))
  expect_equal(r$decision == "reject H0", r$p.value <= 0.05)
  r <- test_adequation(rnorm(40), "HWP", test = "SW")
  expect_true(is.na(r$critical.value))
})

test_that("a negative increment under an increasing model is a certain rejection", {
  r <- test_adequation(c(1, 2, -0.5, 3), "InvGaus", test = c("KS", "AD"))
  expect_equal(r$p.value, c(0, 0))
  expect_true(all(r$decision == "reject H0"))
})

test_that("a failing test does not stop the others", {
  r <- test_adequation(c(1, 1, 1, 1, 1), "HWP", test = c("SW", "AD"))
  expect_true(all(r$decision == "not computed"))
  expect_true(all(nzchar(r$message)))
})

test_that("gof.levy reads the three data formats", {
  set.seed(2)
  t <- c(0, cumsum(round(runif(20, 0.5, 3.5), 1)))
  x <- c(0, cumsum(rnorm(20, 0.5 * diff(t), 1.5 * sqrt(diff(t)))))
  r1 <- gof.levy(data.frame(t, x), "HWP", test = "SW")
  expect_s3_class(r1, "gof.levy")
  expect_false(r1$regular)
  expect_equal(nrow(r1$results), 1)
  expect_equal(r1$steps, t[21] / 20)
  expect_equal(r1$bridge[[1]]$X[c(1, 21)], x[c(1, 21)])

  m <- cbind(t, x, x + 1)
  r2 <- gof.levy(m, "HWP", test = "SW")
  expect_equal(r2$results$group, c("all", "1", "2"))

  liste <- list(data.frame(t, x), data.frame(t = c(0, 1, 2.5, 4, 7), X = c(0, 1, 0.5, 2, 3)))
  r3 <- gof.levy(liste, "HWP", test = "SW")
  expect_equal(r3$results$group, c("1", "2"))
})

test_that("regular inspections are tested without bridge, and the origin is free", {
  set.seed(3)
  t <- seq(10, 70, by = 2)
  x <- cumsum(c(5, rgamma(30, 4, 3)))
  r <- gof.levy(data.frame(t, x), "HGP", test = "VGE", B = 30)
  expect_true(r$regular)
  expect_null(r$bridge)
  expect_null(r$parameters)
  expect_equal(r$increments$dx, diff(x))
  expect_equal(r$steps, 2)
})

test_that("a decreasing observation rejects the gamma model without bridge", {
  set.seed(4)
  t <- c(0, cumsum(runif(10, 0.5, 2)))
  x <- c(0, cumsum(rgamma(10, 2, 3)))
  x[6] <- x[5] - 0.1
  r <- gof.levy(data.frame(t, x), "HGP", test = c("KS", "AD"))
  expect_true(r$out_of_support)
  expect_equal(r$results$p.value, c(0, 0))
  expect_equal(r$results$group, c("all", "all"))
})

test_that("the method ExtBridge keeps all the observations", {
  set.seed(5)
  t <- c(0, cumsum(round(runif(15, 0.5, 3.5), 1)))
  x <- c(0, cumsum(rgamma(15, 10 * diff(t), 5)))
  r <- gof.levy(data.frame(t, x), "HGP", test = "VGE", method = "ExtBridge", B = 20)
  expect_equal(r$method, "ExtBridge")
  expect_equal(r$steps, 0.1)
  expect_equal(r$bridge[[1]]$X[match(round(t, 6), round(r$bridge[[1]]$t, 6))], x)
})

test_that("bad data are refused", {
  expect_error(gof.levy(1:10, "HGP"), "data must be")
  expect_error(gof.levy(data.frame(t = c(0, 2, 1, 3, 4), x = 1:5), "HGP"), "strictly increasing")
  expect_error(gof.levy(data.frame(t = 0:2, x = 0:2), "HGP"), "at least 4 increments")
  expect_error(gof.levy(data.frame(t = 0:5, x = c(0, 1, NA, 3, 4, 5)), "HGP"), "missing")
})

test_that("ExtBridge uses a common step and tests all the trajectories together", {
  set.seed(6)
  liste <- lapply(1:3, function(k) {
    tk <- c(0, cumsum(round(runif(12, 0.5, 3.5), 1)))
    data.frame(t = tk, X = c(0, cumsum(rgamma(12, 2 * diff(tk), 3))))
  })
  r <- gof.levy(liste, "HGP", test = "VGE", method = "ExtBridge", B = 20)
  expect_equal(r$steps, rep(0.1, 3))
  expect_equal(r$results$group, c("all", "1", "2", "3"))
  expect_true(all(r$increments$dx > 0))
  r2 <- gof.levy(liste, "HGP", test = "VGE", B = 20)
  expect_equal(r2$results$group, c("1", "2", "3"))
})

test_that("MultiBridge gives one common grid and a test on all the trajectories", {
  set.seed(7)
  liste <- lapply(c(10, 15, 25), function(n) {
    tk <- c(0, cumsum(round(runif(n, 0.5, 3.5), 1)))
    data.frame(t = tk, X = c(0, cumsum(rgamma(n, 2 * diff(tk), 3))))
  })
  h <- sum(sapply(liste, function(d) max(d$t))) / (10 + 15 + 25)
  r <- gof.levy(liste, "HGP", test = "VGE", method = "MultiBridge", B = 20)
  expect_equal(r$steps, rep(h, 3))
  expect_equal(r$results$group, c("all", "1", "2", "3"))
  expect_equal(r$method, "MultiBridge, option 1")
  r2 <- gof.levy(liste, "HGP", test = "VGE", method = "MultiBridge", option = 2, B = 20)
  expect_equal(length(unique(table(r2$increments$id))), 1)
  expect_error(gof.levy(liste, "HGP", method = "MultiBridge", option = 3), "option must be 1 or 2")
})

test_that("KS is accepted for the Wiener model as the Lilliefors test", {
  set.seed(9)
  x <- rnorm(40)
  r1 <- test_adequation(x, "HWP", test = c("KS", "AD"))
  r2 <- test_adequation(x, "HWP", test = c("Lilliefors", "AD"))
  expect_equal(r1, r2)
  expect_equal(r1$test, c("Lilliefors", "AD"))
  expect_equal(test_adequation(x, "HWP", test = c("ks", "Lilliefors"))$test, "Lilliefors")
  expect_error(test_adequation(x, "HWP", test = "Watson"), "unknown test")
})

test_that("the empty message column is not printed but stays in the data frame", {
  set.seed(10)
  r <- test_adequation(rnorm(40), "HWP", test = c("SW", "AD"))
  expect_true("message" %in% names(r))
  expect_false(any(grepl("message", capture.output(print(r)))))
  r2 <- test_adequation(c(1, 1, 1, 1, 1), "HWP", test = "SW")
  expect_true(any(grepl("message", capture.output(print(r2)))))
})
