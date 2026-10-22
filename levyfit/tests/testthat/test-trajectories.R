plan <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
cut_in_two <- function(x) {
  list(data.frame(t = plan[1:6], X = x[1:6]), data.frame(t = plan[6:11], X = x[6:11]))
}

test_that("a trajectory cut in two gives the same estimates as the whole trajectory", {
  set.seed(3)
  xg <- c(0, cumsum(rgamma(10, shape = 2 * diff(plan), rate = 3)))
  xw <- c(0, cumsum(0.5 * diff(plan) + 1.5 * sqrt(diff(plan)) * rnorm(10)))
  expect_equal(gamma_mle(cut_in_two(xg)), gamma_mle(xg, plan))
  expect_equal(wiener_mle(cut_in_two(xw)), wiener_mle(xw, plan))
  expect_equal(ig_mle(cut_in_two(xg)), ig_mle(xg, plan))
})

test_that("vg_mle pools the increments of several trajectories", {
  set.seed(3)
  g <- rgamma(10, shape = diff(plan) / 2, scale = 2)
  xv <- c(0, cumsum(0.5 * g + 1.5 * sqrt(g) * rnorm(10)))
  expect_equal(vg_mle(cut_in_two(xv)), vg_mle(xv, plan))
})

test_that("gamma_mle on 30 units matches the Monte-Carlo of the report", {
  set.seed(1)
  est <- replicate(300, gamma_mle(lapply(1:30, function(i)
    data.frame(t = plan, X = c(0, cumsum(rgamma(10, shape = 2 * diff(plan), rate = 3)))))))
  expect_equal(rowMeans(est), c(a = 2.017, b = 3.028), tolerance = 0.02)
})

test_that("the data can be given in several forms", {
  set.seed(2)
  x <- c(0, cumsum(rgamma(10, shape = 2 * diff(plan), rate = 3)))
  expect_equal(gamma_mle(data.frame(time = plan, level = x)), gamma_mle(x, plan))
  expect_equal(gamma_mle(list(data.frame(t = plan, X = x))), gamma_mle(x, plan))
})

test_that("igam_mle accepts several trajectories", {
  skip_on_cran()
  x <- c(0, 2.212327259, 2.81893462, 3.413725025, 4.801409076, 6.228550224,
         7.714096072, 7.949339032, 8.524511739, 9.238423227, 9.359099294)
  expect_equal(igam_mle(list(data.frame(t = plan, X = x))), igam_mle(x, plan))
})

test_that("invalid lists of trajectories are refused", {
  expect_error(gamma_mle(list(data.frame(t = 0, X = 0), data.frame(t = 0:2, X = 0:2))),
               "at least two observations")
  expect_error(wiener_mle(list(data.frame(t = c(0, 2, 1), X = c(0, 1, 2)))), "strictly increasing")
  expect_error(igam_mle(list(data.frame(t = c(1, 2), X = c(0, 1)))), "start at")
})

test_that("increments_data gives the increments directly, without loss of precision", {
  set.seed(9)
  dx <- rgamma(300, shape = 0.02, rate = 3)
  expect_true(min(dx) < 1e-30)
  x <- c(1000, 1000 + cumsum(dx))
  expect_true(any(diff(x) == 0))
  e <- gamma_mle(increments_data(dx, dt = 0.1))
  expect_true(all(is.finite(e)))
  expect_equal(wiener_mle(increments_data(c(1, 2, 0.5, 3), c(1, 2, 1, 2))),
               wiener_mle(c(0, 1, 3, 3.5, 6.5), c(0, 1, 3, 4, 6)))
  expect_error(increments_data(c(1, NA)), "dx")
  expect_error(increments_data(1:3, dt = c(1, 2)), "dt")
  expect_error(igam_mle(increments_data(1:3)), "levels")
})
