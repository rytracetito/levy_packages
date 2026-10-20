test_that("inspection_times reproduces the original examples", {
  expect_equal(inspection_times("increasing", n = 10, digits = 0)$time,
               c(0, 29, 46, 58, 67, 75, 81, 87, 92, 96, 100))
  expect_equal(inspection_times(function(x) x^2, n = 10)$time, (0:10)^2)
  expect_equal(inspection_times(sqrt, n = 10, digits = 0)$time,
               c(0, 32, 45, 55, 63, 71, 77, 84, 89, 95, 100))
})

test_that("decreasing is the mirror image of increasing", {
  inc <- inspection_times("increasing", n = 20)$time
  dec <- inspection_times("decreasing", n = 20)$time
  expect_equal(dec, 100 - rev(inc))
})

test_that("several trajectories are returned with their id", {
  ins <- inspection_times("dec-inc", n = c(10, 20), n_traj = 2, tmax = c(50, 100))
  expect_equal(as.vector(table(ins$id)), c(11, 21))
  expect_equal(as.vector(tapply(ins$time, ins$id, max)), c(50, 100))
})

test_that("invalid schedules are refused", {
  expect_error(inspection_times("foo"), "unknown pattern")
  expect_error(inspection_times(function(x) x * (10 - x)), "not increasing")
  expect_error(inspection_times("decreasing", n = 50, digits = 0), "rounding")
})

test_that("a sets the contrast of the alternating patterns", {
  saut <- function(t) { g <- diff(t); max(g[-1] / g[-length(g)], g[-length(g)] / g[-1]) }
  for (p in c("inc-dec", "dec-inc", "inc-dec-inc", "dec-inc-dec")) {
    expect_equal(inspection_times(p, n = 20)$time, inspection_times(p, n = 20, a = 0.5)$time)
    expect_equal(inspection_times(p, n = 20, a = 0)$time, seq(0, 100, by = 5))
    g <- diff(inspection_times(p, n = 200, a = 0.5)$time)
    expect_equal(max(g) / min(g), 3, tolerance = 0.01)
    expect_lt(saut(inspection_times(p, n = 20)$time), 1.35)
    expect_silent(inspection_times(p, n = 50, digits = 0))
  }
  expect_equal(inspection_times("increasing", n = 20, a = 0.9)$time, inspection_times("increasing", n = 20)$time)
})

test_that("the alternating patterns have the expected phases", {
  phases <- function(t) { g <- diff(t); r <- rle(sign(diff(g))[abs(diff(g)) > 1e-9]); paste(ifelse(r$values > 0, "+", "-"), collapse = "") }
  for (n in c(10, 20, 50)) {
    expect_equal(phases(inspection_times("inc-dec", n = n)$time), "-+")
    expect_equal(phases(inspection_times("dec-inc", n = n)$time), "+-")
    expect_equal(phases(inspection_times("inc-dec-inc", n = n)$time), "-+-")
    expect_equal(phases(inspection_times("dec-inc-dec", n = n)$time), "+-+")
  }
})

test_that("invalid values of a are refused", {
  for (s in list(-0.1, 1, 2, NA, c(0.2, 0.3), "a")) expect_error(inspection_times("inc-dec", a = s), "a must be")
})
