#' Plot inspection times
#'
#' Draws the inspection times of one or several trajectories with `ggplot2`:
#' one point per inspection, and one row per trajectory when there are several.
#' It shows at a glance where the inspections are frequent and where they are
#' rare.
#'
#' @param t Inspection times, in one of the three forms accepted by the `sim_`
#'   functions:
#'   * a numeric vector, for one trajectory, for instance
#'     `inspection_times("regular", n = 20, tmax = 10)$time`;
#'   * a list of numeric vectors, one per trajectory, for instance
#'     `list(c(0, 1, 3), c(0, 2, 4, 6))`;
#'   * the result of [inspection_times()], one plan per trajectory.
#'
#'   The times of each trajectory start at 0 and are strictly increasing.
#' @param title Title of the plot, a character string; by default
#'   `"Inspection times"`.
#'
#' @return A `ggplot` object. It is displayed when printed, and it can be
#'   modified with the usual `ggplot2` functions, added with `+`: for instance
#'   `+ ggplot2::labs(x = "Time (hours)")` to change the axis label.
#'
#' @seealso [inspection_times()], [plot_traj()]
#'
#' @examples
#' # One trajectory with 20 regular inspections up to 10 (step 0.5)
#' plot_plan(inspection_times("regular", n = 20, tmax = 10)$time)
#'
#' # Inspections more and more frequent over time, with a custom title
#' t <- inspection_times("decreasing", n = 30, tmax = 100)
#' plot_plan(t, "Decreasing gaps between inspections")
#'
#' # Five trajectories, each with its own times: one row per trajectory
#' plan <- inspection_times("inc-dec", n = c(10, 15, 20, 25, 30), n_traj = 5, tmax = 100)
#' plot_plan(plan)
#'
#' # Different inspection times given as a list
#' plot_plan(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)))
#'
#' # Change the axis label
#' plot_plan(t) + ggplot2::labs(x = "Time (hours)")
#' @export
plot_plan <- function(t, title = "Inspection times") {
  plans <- as_plans(t, 1)
  d <- data.frame(id = rep(seq_along(plans), lengths(plans)), t = unlist(plans))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = t, y = id)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = title, x = "Time", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank())
  if (length(plans) == 1) {
    p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       panel.grid.major.y = ggplot2::element_blank())
  } else {
    p + ggplot2::scale_y_reverse(breaks = seq_along(plans)) +
      ggplot2::labs(y = "Trajectory")
  }
}
