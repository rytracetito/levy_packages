as_plans <- function(t, n_traj) {
  if (is.data.frame(t)) t <- unname(split(t[[2]], t[[1]]))
  if (!is.list(t)) t <- rep(list(t), n_traj)
  for (tt in t) {
    if (!is.numeric(tt) || length(tt) < 2) stop("each trajectory needs at least two inspection times")
    if (tt[1] != 0) stop("inspection times must start at 0")
    if (any(diff(tt) <= 0)) stop("inspection times must be strictly increasing")
  }
  t
}

simulate_paths <- function(t, n_traj, one_path, model, parameters) {
  plans <- as_plans(t, n_traj)
  paths <- lapply(plans, function(tt) data.frame(t = tt, X = one_path(tt)))
  structure(paths, class = "levysim", model = model, parameters = parameters)
}

#' Simulated data as a data frame
#'
#' Gathers the trajectories simulated by a `sim_` function into a single data
#' frame, one row per inspection, ready to be saved, summarized or used by
#' other functions.
#'
#' @param x Result of a simulation function: [sim_gamma()], [sim_wiener()],
#'   [sim_ig()], [sim_vg()], [sim_meixner()] or [sim_igam()].
#'
#' @return A data frame with one row per inspection and three columns:
#'   * `id`: number of the trajectory (1, 2, ...);
#'   * `t`: inspection time;
#'   * `X`: simulated degradation level.
#'
#' To work on a single trajectory, use either `x[[k]]` on the simulation
#' result or `subset(traj_df(x), id == k)`.
#'
#' @seealso [plot_traj()]
#'
#' @examples
#' set.seed(1)
#' x <- sim_gamma(seq(0, 10, by = 1), a = 2, b = 3, n_traj = 3)
#' d <- traj_df(x)
#' head(d)
#'
#' # Degradation level of each trajectory at the last inspection
#' d[d$t == 10, ]
#'
#' # Save the simulated data in a csv file
#' f <- tempfile(fileext = ".csv")
#' write.csv(d, f, row.names = FALSE)
#' @export
traj_df <- function(x) {
  do.call(rbind, lapply(seq_along(x), function(i) data.frame(id = i, t = x[[i]]$t, X = x[[i]]$X)))
}

#' Plot simulated trajectories
#'
#' Draws the trajectories simulated by a `sim_` function with `ggplot2`: one
#' grey line per trajectory, with a point at each inspection.
#'
#' @param x Result of a simulation function: [sim_gamma()], [sim_wiener()],
#'   [sim_ig()], [sim_vg()], [sim_meixner()] or [sim_igam()].
#' @param title Title of the plot, a character string; by default the name of
#'   the process, for instance `"gamma process"`.
#'
#' @return A `ggplot` object. It is displayed when printed, and it can be
#'   modified with the usual `ggplot2` functions, added with `+`: for instance
#'   `+ ggplot2::geom_hline(yintercept = 10)` to draw a failure threshold, or
#'   `+ ggplot2::labs(x = "Time (hours)")` to change an axis label.
#'
#' @seealso [traj_df()]
#'
#' @examples
#' set.seed(1)
#' x <- sim_wiener(seq(0, 10, by = 0.5), mu = 0.5, sigma = 1.5, n_traj = 5)
#' plot_traj(x)
#'
#' # Custom title
#' plot_traj(x, "Five trajectories of a Wiener process")
#'
#' # Add a failure threshold and change the axis labels
#' plot_traj(x) +
#'   ggplot2::geom_hline(yintercept = 5, linetype = "dashed", colour = "red") +
#'   ggplot2::labs(x = "Time (hours)", y = "Degradation (%)")
#' @export
plot_traj <- function(x, title = NULL) {
  if (is.null(title)) title <- paste(attr(x, "model"), "process")
  d <- traj_df(x)
  ggplot2::ggplot(d, ggplot2::aes(x = t, y = X, group = id)) +
    ggplot2::geom_line(colour = "grey50") +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::labs(title = title, x = "Time", y = "Degradation") +
    ggplot2::theme_minimal()
}

#' @export
print.levysim <- function(x, ...) {
  par <- attr(x, "parameters")
  cat(length(x), "simulated trajectory(ies) of the", attr(x, "model"), "process\n")
  cat("parameters:", paste(names(par), signif(par, 4), sep = " = ", collapse = ", "), "\n\n")
  print(utils::head(traj_df(x)))
  invisible(x)
}
