#' Simulation of the gamma process
#'
#' Simulates exactly one or several trajectories of a gamma process, homogeneous
#' or not, at given inspection times.
#'
#' @section Model:
#' Let \eqn{A} be a non-decreasing function with \eqn{A(0) = 0}, \eqn{a > 0} and
#' \eqn{b > 0}. The gamma process satisfies \eqn{X(0) = 0}, has independent
#' increments, and for \eqn{0 \le s < t},
#' \deqn{X(t) - X(s) \sim Ga\big(a (A(t) - A(s)),\ b\big),}
#' where \eqn{Ga(k, b)} is the gamma distribution with shape \eqn{k} and
#' **rate** \eqn{b}, of density \eqn{b^k x^{k-1} e^{-bx} / \Gamma(k)}. Then
#' \eqn{E[X(t)] = a A(t) / b} and \eqn{Var(X(t)) = a A(t) / b^2}. The process is
#' increasing and purely made of jumps. With \eqn{A(t) = t} (default), it is
#' the homogeneous gamma process, a Lévy process. A classical non-homogeneous
#' choice is
#' \eqn{A(t) = t^\beta}.
#'
#' @section Exact simulation:
#' 1. \eqn{\Delta A_i = A(t_i) - A(t_{i-1})}, \eqn{i = 1, \dots, n}.
#' 2. Draw independently \eqn{\Delta X_i \sim Ga(a \Delta A_i, b)}.
#' 3. \eqn{X(t_0) = 0} and \eqn{X(t_i) = \sum_{j \le i} \Delta X_j}.
#'
#' @param t Inspection times. The times of each trajectory must start at 0 and
#'   be strictly increasing. Three forms are accepted:
#'   * a numeric vector, for instance `seq(0, 10, by = 0.5)`: all the
#'     trajectories share these times, and their number is `n_traj`;
#'   * a list of vectors, one per trajectory, for instance
#'     `list(c(0, 1, 3), c(0, 2, 4, 6))`: each trajectory has its own times
#'     (`n_traj` is then ignored);
#'   * the data frame returned by [inspection_times()], for instance
#'     `inspection_times("increasing", n = 20, n_traj = 5, tmax = 10)`.
#' @param a Shape parameter, \eqn{a > 0}.
#' @param b Rate parameter, \eqn{b > 0} (the scale is \eqn{1 / b}).
#' @param A Shape function of time: an R function of `t`, non-decreasing with
#'   \eqn{A(0) = 0}, for instance `function(t) t^1.5`. The default
#'   `function(t) t` gives the homogeneous process.
#' @param n_traj Number of trajectories to simulate when `t` is a vector
#'   (default 1).
#'
#' @return An object of class `levysim`: a list with one data frame per
#'   trajectory, with columns `t` (inspection times) and `X` (simulated
#'   degradation levels, with `X = 0` at `t = 0`). The \eqn{k}-th trajectory is
#'   `x[[k]]`. Printing the object shows the model, the parameters and the
#'   first rows. Use [plot_traj()] to plot the trajectories and [traj_df()] to
#'   gather them in a single data frame.
#'
#' @references
#' van Noortwijk, J. M. (2009). A survey of the application of gamma processes
#' in maintenance. *Reliability Engineering & System Safety*, 94(1), 2-21.
#' \doi{10.1016/j.ress.2007.03.019}
#'
#' @seealso [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()]
#'
#' @examples
#' # One trajectory inspected every 0.5 time unit up to 10
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_gamma(t, a = 2, b = 3)
#' x
#' plot_traj(x)
#'
#' # Ten trajectories with the same inspection times
#' x <- sim_gamma(t, a = 2, b = 3, n_traj = 10)
#' plot_traj(x)
#'
#' # Irregular inspections, more and more frequent over time
#' plan <- inspection_times("increasing", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_gamma(plan, a = 2, b = 3)
#' plot_traj(x)
#'
#' # Different inspection times for each trajectory
#' x <- sim_gamma(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), a = 2, b = 3)
#' x[[2]]
#'
#' # Non-homogeneous process with A(t) = t^1.5
#' x <- sim_gamma(t, a = 2, b = 3, A = function(t) t^1.5, n_traj = 10)
#' plot_traj(x, "Non-homogeneous gamma process")
#'
#' # All the simulated data in a single data frame
#' head(traj_df(x))
#' @export
sim_gamma <- function(t, a, b, A = function(t) t, n_traj = 1) {
  if (a <= 0 || b <= 0) stop("a and b must be positive")
  one_path <- function(t) {
    dA <- diff(A(t))
    if (any(dA < 0)) stop("A must be non-decreasing")
    c(0, cumsum(stats::rgamma(length(dA), shape = a * dA, rate = b)))
  }
  simulate_paths(t, n_traj, one_path, "gamma", c(a = a, b = b))
}
