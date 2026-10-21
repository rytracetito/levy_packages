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
  same <- all(vapply(plans, function(tt) identical(tt, plans[[1]]), TRUE))
  if (same) {
    X <- vapply(plans, one_path, numeric(length(plans[[1]])))
    paths <- cbind(t = plans[[1]], X)
    colnames(paths) <- c("t", paste0("X", seq_along(plans)))
    class(paths) <- c("levysim", "matrix", "array")
  } else {
    paths <- lapply(plans, function(tt) data.frame(t = tt, X = one_path(tt)))
    class(paths) <- "levysim"
  }
  attr(paths, "model") <- model
  attr(paths, "parameters") <- parameters
  paths
}

as_path_list <- function(x) {
  if (is.matrix(x)) {
    return(lapply(seq_len(ncol(x) - 1) + 1, function(k) data.frame(t = x[, 1], X = x[, k])))
  }
  lapply(seq_along(x), function(k) x[[k]])
}

#' Simulated data as a data frame
#'
#' Gathers the trajectories simulated by a `sim_` function into a single data
#' frame, one row per inspection, ready to be saved, summarized or used by
#' other functions. Optionally adds the time steps and the increments between
#' two consecutive inspections.
#'
#' @param x Result of a simulation function of the package, for instance
#'   [sim_gamma()].
#' @param increments `FALSE` (default) for the levels only; `TRUE` to add the
#'   columns `dt` and `dX`.
#'
#' @return A data frame with one row per inspection and the columns:
#'   * `id`: number of the trajectory (1, 2, ...);
#'   * `t`: inspection time;
#'   * `X`: simulated degradation level at time `t`;
#'   * with `increments = TRUE`, also `dt`, the time elapsed since the previous
#'     inspection, and `dX`, the increment of the level since the previous
#'     inspection (`NA` at the first inspection, `t = 0`, of each trajectory).
#'
#' To work on a single trajectory, use `subset(traj_df(x), id == k)`.
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
#' # Time steps and increments between two inspections
#' head(traj_df(x, increments = TRUE))
#'
#' # Trajectories with different inspection times
#' y <- sim_gamma(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8)), a = 2, b = 3)
#' traj_df(y, increments = TRUE)
#'
#' # Save the simulated data in a csv file
#' f <- tempfile(fileext = ".csv")
#' write.csv(d, f, row.names = FALSE)
#' @export
traj_df <- function(x, increments = FALSE) {
  paths <- as_path_list(x)
  do.call(rbind, lapply(seq_along(paths), function(i) {
    d <- data.frame(id = i, t = paths[[i]]$t, X = paths[[i]]$X)
    if (increments) {
      d$dt <- c(NA, diff(d$t))
      d$dX <- c(NA, diff(d$X))
    }
    d
  }))
}

#' Plot simulated trajectories
#'
#' Draws the trajectories simulated by a `sim_` function with `ggplot2`: one
#' grey line per trajectory, with a point at each inspection.
#'
#' @param x Result of a simulation function of the package, for instance
#'   [sim_gamma()].
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
  n <- if (is.matrix(x)) ncol(x) - 1 else length(x)
  cat(n, "simulated trajectory(ies) of the", attr(x, "model"), "process")
  cat(if (is.matrix(x)) ", same inspection times (matrix)\n" else ", own inspection times (list)\n")
  cat("parameters:", paste(names(par), signif(par, 4), sep = " = ", collapse = ", "), "\n\n")
  if (is.matrix(x)) print(utils::head(unclass(x)[, seq_len(min(ncol(x), 6)), drop = FALSE])) else print(utils::head(traj_df(x)))
  invisible(x)
}
