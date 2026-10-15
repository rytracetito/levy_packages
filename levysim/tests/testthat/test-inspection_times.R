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
