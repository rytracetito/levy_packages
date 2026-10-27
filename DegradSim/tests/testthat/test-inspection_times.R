test_that("inspection_times reproduces the original examples", {
  expect_equal(inspection_times("decreasing", n = 10, digits = 0)$time,
               c(0, 29, 46, 58, 67, 75, 81, 87, 92, 96, 100))
  expect_equal(inspection_times(function(x) x^2, n = 10)$time, (0:10)^2)
  expect_equal(inspection_times(sqrt, n = 10, digits = 0)$time,
               c(0, 32, 45, 55, 63, 71, 77, 84, 89, 95, 100))
})

test_that("increasing is the mirror image of decreasing", {
  dec <- inspection_times("decreasing", n = 20)$time
  inc <- inspection_times("increasing", n = 20)$time
  expect_equal(inc, 100 - rev(dec))
})

test_that("the pattern decreasing gives a ratio of about 11 between the extreme gaps", {
  g <- diff(inspection_times("decreasing", n = 1000)$time)
  expect_equal(g[1] / g[1000], 11, tolerance = 0.01)
})

test_that("several trajectories are returned with their id", {
  ins <- inspection_times("inc-dec", n = c(10, 20), n_traj = 2, tmax = c(50, 100))
  expect_equal(as.vector(table(ins$id)), c(11, 21))
  expect_equal(as.vector(tapply(ins$time, ins$id, max)), c(50, 100))
})

test_that("invalid schedules are refused", {
  expect_error(inspection_times("foo"), "unknown pattern")
  expect_error(inspection_times("inc-dec-inc"), "unknown pattern")
  expect_error(inspection_times("dec-inc-dec"), "unknown pattern")
  expect_error(inspection_times(function(x) x * (1 - x)), "strictly increasing")
  expect_error(inspection_times(n = 2.5), "n must be")
  expect_error(inspection_times(n = 0), "n must be")
  expect_error(inspection_times(n = c(5, 6), n_traj = 3), "n must be")
  expect_error(inspection_times(tmax = -1), "tmax must be")
  expect_error(inspection_times(tmax = c(5, 6), n_traj = 3), "tmax must be")
  expect_error(inspection_times(n_traj = 0), "n_traj")
  expect_error(inspection_times(digits = -1), "digits")
  expect_error(inspection_times(digits = 0.5), "digits")
})

test_that("times made equal by the rounding are removed with a warning", {
  expect_warning(p <- inspection_times("increasing", n = 50, digits = 0), "inspections after time 0 are kept")
  expect_true(all(diff(p$time) > 0))
  expect_equal(p$time[1], 0)
  expect_equal(max(p$time), 100)
  expect_lt(nrow(p), 51)
  expect_silent(inspection_times("regular", n = 10, digits = 0))
  expect_error(inspection_times("regular", n = 10, tmax = 0.1, digits = 0), "too strong")
})

test_that("a sets the contrast of the alternating patterns", {
  saut <- function(t) { g <- diff(t); max(g[-1] / g[-length(g)], g[-length(g)] / g[-1]) }
  for (p in c("dec-inc", "inc-dec")) {
    expect_equal(inspection_times(p, n = 20)$time, inspection_times(p, n = 20, a = 0.5)$time)
    expect_equal(inspection_times(p, n = 20, a = 0)$time, seq(0, 100, by = 5))
    g <- diff(inspection_times(p, n = 200, a = 0.5)$time)
    expect_equal(max(g) / min(g), 3, tolerance = 0.01)
    expect_lt(saut(inspection_times(p, n = 20)$time), 1.35)
    expect_silent(inspection_times(p, n = 50, digits = 0))
  }
  expect_equal(inspection_times("decreasing", n = 20, a = 0.9)$time, inspection_times("decreasing", n = 20)$time)
})

test_that("the alternating patterns have the expected phases", {
  phases <- function(t) { g <- diff(t); r <- rle(sign(diff(g))[abs(diff(g)) > 1e-9]); paste(ifelse(r$values > 0, "+", "-"), collapse = "") }
  for (n in c(10, 20, 50)) {
    expect_equal(phases(inspection_times("dec-inc", n = n)$time), "-+")
    expect_equal(phases(inspection_times("inc-dec", n = n)$time), "+-")
  }
})

test_that("invalid values of a are refused", {
  for (s in list(-0.1, 1, 2, NA, c(0.2, 0.3), "a")) expect_error(inspection_times("dec-inc", a = s), "a must be")
})

test_that("the grid [0, 1] gives the same times as the former grid [0, 10]", {
  old <- function(f, n, tmax = 100) {
    fx <- f(seq(0, 10, length.out = n + 1))
    tmax * (fx - fx[1]) / (fx[n + 1] - fx[1])
  }
  a <- 0.5
  expect_equal(inspection_times("decreasing", n = 30)$time, old(function(x) log(x + 1), 30))
  expect_equal(inspection_times("increasing", n = 30)$time, old(function(x) -log(11 - x), 30))
  expect_equal(inspection_times("dec-inc", n = 30)$time,
               old(function(x) x + a * 10 / (2 * pi) * sin(2 * pi * x / 10), 30))
  expect_equal(inspection_times("inc-dec", n = 30)$time,
               old(function(x) x - a * 10 / (2 * pi) * sin(2 * pi * x / 10), 30))
})
