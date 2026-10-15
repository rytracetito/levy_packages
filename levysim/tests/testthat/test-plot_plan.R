test_that("plot_plan draws one point per inspection and one row per trajectory", {
  p <- plot_plan(seq(0, 10, by = 0.5))
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 21)
  expect_equal(unique(p$data$id), 1)

  plan <- inspection_times("increasing", n = c(5, 10, 15), n_traj = 3, tmax = 10)
  p <- plot_plan(plan, "Three plans")
  expect_equal(nrow(p$data), 3 + 5 + 10 + 15)
  expect_equal(sort(unique(p$data$id)), 1:3)
  expect_equal(p$labels$title, "Three plans")

  p <- plot_plan(list(c(0, 1, 3), c(0, 2, 4, 6)))
  expect_equal(as.vector(table(p$data$id)), c(3, 4))
})

test_that("plot_plan rejects invalid times", {
  expect_error(plot_plan(c(1, 2, 3)), "start at 0")
  expect_error(plot_plan(c(0, 2, 1)), "strictly increasing")
})
