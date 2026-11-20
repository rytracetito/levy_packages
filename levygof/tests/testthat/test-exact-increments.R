test_that("exact increments avoid a false rejection for sure", {
  x <- list(data.frame(t = c(0, 1, 2, 3, 4, 5), X = c(0, 2, 3.5, 50, 50, 52)))
  attr(x, "dx") <- list(c(2, 1.5, 46.5, 1e-15, 2))
  r <- gof.levy(x, "HGP", test = "KS", B = 20)
  expect_false(r$out_of_support)
  attr(x, "dx") <- NULL
  expect_true(gof.levy(x, "HGP", test = "KS", B = 20)$out_of_support)
})

test_that("the results do not change when the exact increments equal the differences of the levels", {
  set.seed(1)
  x <- DegradSim::sim_gamma(DegradSim::inspection_times("decreasing", n = 30, tmax = 60, digits = 1)$time, a = 2, b = 3)
  y <- x
  attr(y, "dx") <- NULL
  set.seed(2); a <- gof.levy(x, "HGP", test = c("KS", "AD"), B = 50)
  set.seed(2); b <- gof.levy(y, "HGP", test = c("KS", "AD"), B = 50)
  expect_equal(a$parameters, b$parameters)
  expect_equal(a$results$p.value, b$results$p.value)
})
