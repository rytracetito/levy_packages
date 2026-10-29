#' Plot observed and resampled trajectories
#'
#' Draws, for one or several trajectories, the observed values at the
#' inspection times and the values resampled on the regular grid by a
#' resampling function, with `ggplot2`.
#'
#' @param x An object of class `"levybridge"`, returned by a resampling
#'   function of the package, for instance [gamma_resampling()].
#' @param trajectory Indices of the trajectories to plot in `x`, for instance
#'   `1` or `1:3` (default: the first one). With several trajectories, one
#'   panel is drawn per trajectory.
#' @param title Title of the plot, a character string; by default the name of
#'   the bridge, for instance `"gamma bridge"`.
#'
#' @return A `ggplot` object. It is displayed when printed, and it can be
#'   modified with the usual `ggplot2` functions, added with `+`: for instance
#'   `+ ggplot2::labs(x = "Time (hours)")`.
#'
#' @seealso [levyBridge-package]
#'
#' @examples
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3)
#' r <- gamma_resampling(x, theta = c(a = 2, b = 3))
#'
#' # Observed levels (red) and resampled levels on the regular grid (blue)
#' plot_bridge(r)
#'
#' # Custom title and axis labels
#' plot_bridge(r, title = "Gamma bridge, a = 2") +
#'   ggplot2::labs(x = "Time (hours)", y = "Degradation (%)")
#'
#' # Several trajectories: one panel per trajectory
#' plan <- DegradSim::inspection_times("decreasing", n = 10, n_traj = 3,
#'                                     tmax = c(20, 30, 40), digits = 1)
#' xs <- DegradSim::sim_gamma(plan, a = 2, b = 3)
#' plot_bridge(gamma_resampling(xs, theta = c(a = 2, b = 3)), trajectory = 1:3)
#' @export
plot_bridge <- function(x, trajectory = 1, title = NULL) {
  if (!inherits(x, "levybridge")) stop("x must be the result of a _resampling function")
  if (!all(trajectory %in% seq_along(x))) stop("trajectory must be between 1 and ", length(x))
  if (is.null(title)) title <- paste(attr(x, "model"), "bridge")
  obs <- attr(x, "observed")
  d <- do.call(rbind, lapply(trajectory, function(k) rbind(
    data.frame(id = k, t = obs[[k]]$t, X = obs[[k]]$X, type = "Observed"),
    data.frame(id = k, t = x[[k]]$t, X = x[[k]]$X, type = "Resampled (regular grid)"))))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = t, y = X, colour = type, group = type)) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(size = type)) +
    ggplot2::scale_colour_manual(values = c("Observed" = "#C2410C", "Resampled (regular grid)" = "#2563EB")) +
    ggplot2::scale_size_manual(values = c("Observed" = 2.6, "Resampled (regular grid)" = 1.6)) +
    ggplot2::labs(title = title, x = "Time", y = "Degradation", colour = NULL, size = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  if (length(trajectory) > 1) {
    p <- p + ggplot2::facet_wrap(~ id, scales = "free",
                                 labeller = ggplot2::as_labeller(function(v) paste("Trajectory", v)))
  }
  p
}

utils::globalVariables(c("t", "X", "type", "id"))
