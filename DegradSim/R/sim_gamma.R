#' Simulation of the gamma process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' a homogeneous gamma process (a Lévy process) or of a non-homogeneous gamma
#' process, obtained by a time change \eqn{A}.
#'
#' @section Model:
#' \eqn{\Gamma(k, b)} denotes the gamma distribution with shape \eqn{k > 0} and
#' **rate** \eqn{b > 0} (scale \eqn{1 / b}), of density
#' \eqn{b^k x^{k-1} e^{-bx} / \Gamma(k)}, \eqn{x > 0}.
#'
#' **Homogeneous gamma process** (`A = NULL`, default). Let \eqn{a > 0} and
#' \eqn{b > 0}. The process satisfies \eqn{X_0 = 0}, has independent
#' increments, and for \eqn{0 \le s < t},
#' \deqn{X_t - X_s \sim \Gamma\big(a (t - s),\ b\big).}
#' Then \eqn{E[X_t] = a t / b} and \eqn{Var(X_t) = a t / b^2}. It is a Lévy
#' process, increasing and purely made of jumps.
#'
#' **Non-homogeneous gamma process** (`A` given). With the time change
#' \eqn{A} defined in [DegradSim-package], \eqn{X_t = Y_{A(t)}}, where \eqn{Y}
#' is the homogeneous gamma process with \eqn{a = 1} (identifiability, see
#' [DegradSim-package]) and rate \eqn{b}:
#' \deqn{X_t - X_s \sim \Gamma\big(A(t) - A(s),\ b\big),}
#' \eqn{E[X_t] = A(t) / b} and \eqn{Var(X_t) = A(t) / b^2}. This is the
#' gamma process with shape function \eqn{A} of Kahle et al. (2016,
#' chapter 2) and of Cha and Mercier (2019, Example 2). A classical choice is
#' \eqn{A(t) = c\, t^\beta}.
#'
#' @section Exact simulation:
#' Random walk approximation, described in [DegradSim-package]: the
#' increments \eqn{\Delta X_i \sim \Gamma(\Delta_i, b)}, with
#' \eqn{\Delta_i = a (t_i - t_{i-1})} or \eqn{\Delta_i = A(t_i) - A(t_{i-1})},
#' are drawn independently with [stats::rgamma()], then cumulated.
#'
#' @param t Inspection times. The times of each trajectory must start at 0 and
#'   be strictly increasing. Three forms are accepted:
#'   * a numeric vector, for instance `seq(0, 10, by = 0.5)`: all the
#'     trajectories share these times, and their number is `n_traj`;
#'   * a list of vectors, one per trajectory, for instance
#'     `list(c(0, 1, 3), c(0, 2, 4, 6))`: each trajectory has its own times
#'     (`n_traj` is then ignored);
#'   * the data frame returned by [inspection_times()], which gives one
#'     plan per trajectory (`n_traj` is then ignored).
#' @param a Shape parameter of the homogeneous process, \eqn{a > 0}: the
#'   increment over a time \eqn{\Delta t} has shape \eqn{a \Delta t}. Not
#'   given when `A` is given (the shape is then \eqn{A(t) - A(s)}).
#' @param b Rate parameter, \eqn{b > 0} (the scale is \eqn{1 / b}).
#' @param A `NULL` (default) for the homogeneous process, or the time change
#'   \eqn{A} of the non-homogeneous process (defined in [DegradSim-package]):
#'   a vectorized R function of `t`, for instance `function(t) 2 * t^1.5`.
#' @param n_traj Number of trajectories to simulate when `t` is a vector
#'   (default 1), a positive integer.
#'
#' @return An object of class `DegradSim`, whose form depends on the
#'   inspection times:
#'   * **same times for all the trajectories** (`t` a vector, or a plan of
#'     [inspection_times()] with the same times for all the trajectories): a
#'     matrix with columns `t` (inspection times), `X1`, ..., `Xn` (simulated
#'     degradation levels of the trajectories 1 to \eqn{n}, equal to 0 at
#'     `t = 0`). The levels of the \eqn{k}-th trajectory are `x[, k + 1]` or
#'     `x[, "Xk"]`;
#'   * **different times** (`t` a list, or a plan of [inspection_times()] with
#'     different times): a list with one data frame per trajectory, with
#'     columns `t` and `X`. The \eqn{k}-th trajectory is the data frame
#'     `x[[k]]`, its times are `x[[k]]$t` (or `x[[k]][, 1]`) and its levels
#'     `x[[k]]$X` (or `x[[k]][, 2]`).
#'
#'   The model, the parameters and, for a non-homogeneous process, the time
#'   change are stored as attributes. Printing the object shows them with the
#'   first simulated values. Use [traj_df()] to gather the data in a single
#'   data frame (with the increments if needed) and [plot_traj()] to plot the
#'   trajectories.
#'
#' @references
#' Kahle, W., Mercier, S. and Paroissin, C. (2016). *Degradation Processes in
#' Reliability*. ISTE and Wiley, London. ISBN 978-1-84821-888-8.
#'
#' Cha, J. H. and Mercier, S. (2019). Transformed Lévy processes as
#' state-dependent wear models. *Advances in Applied Probability*, 51(2),
#' 468-486. \doi{10.1017/apr.2019.21}
#'
#' van Noortwijk, J. M. (2009). A survey of the application of gamma processes
#' in maintenance. *Reliability Engineering & System Safety*, 94(1), 2-21.
#' \doi{10.1016/j.ress.2007.03.019}
#'
#' @seealso [sim_ig()] and [sim_wiener()] for the other processes with a time
#'   change, [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()]
#'
#' @examples
#' # 1. Homogeneous process: one trajectory inspected every 0.5 time unit
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_gamma(t, a = 2, b = 3)
#' x
#' plot_traj(x)
#'
#' # 2. Ten trajectories with the same inspection times
#' x <- sim_gamma(t, a = 2, b = 3, n_traj = 10)
#' plot_traj(x)
#'
#' # Same times for all the trajectories: a matrix with the column t,
#' # then one column per trajectory (X1, X2, ...)
#' head(x[, 1:4])
#' x[, "X2"]
#' x[, 3]
#'
#' # 3. Irregular inspections, more and more frequent over time
#' plan <- inspection_times("increasing", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_gamma(plan, a = 2, b = 3)
#' plot_traj(x)
#'
#' # 4. Different inspection times for each trajectory: a list of data frames
#' x <- sim_gamma(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), a = 2, b = 3)
#' x[[2]]
#' x[[2]]$X
#' x[[2]][, 2]
#'
#' # 5. Non-homogeneous process with the time change A(t) = 2 t^1.5
#' #    (Cha and Mercier, 2019, Example 2): X_t - X_s ~ Gamma(A(t) - A(s), b)
#' x <- sim_gamma(t, b = 3, A = function(t) 2 * t^1.5, n_traj = 10)
#' x
#' plot_traj(x, "Non-homogeneous gamma process")
#'
#' # 6. A(t) = a t gives back the homogeneous process with shape a
#' set.seed(2)
#' y1 <- sim_gamma(t, a = 2, b = 3)
#' set.seed(2)
#' y2 <- sim_gamma(t, b = 3, A = function(t) 2 * t)
#' all.equal(unclass(y1)[, 2], unclass(y2)[, 2])
#'
#' # 7. All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#' @export
sim_gamma <- function(t, a, b, A = NULL, n_traj = 1) {
  if (missing(b) || !is.numeric(b) || length(b) != 1 || !(b > 0)) stop("b must be a positive number")
  if (is.null(A)) {
    if (missing(a)) stop("a is missing: give a (homogeneous process) or A (non-homogeneous process)")
    if (!is.numeric(a) || length(a) != 1 || !(a > 0)) stop("a must be a positive number")
    shape <- function(t) a * diff(t)
    model <- "gamma"
    parameters <- c(a = a, b = b)
  } else {
    if (!missing(a)) stop("give a or A, not both: with A, the shape of X_t - X_s is A(t) - A(s) (a = 1)")
    shape <- function(t) time_change(A, t)
    model <- "non-homogeneous gamma"
    parameters <- c(b = b)
  }
  one_path <- function(t) {
    k <- shape(t)
    c(0, cumsum(stats::rgamma(length(k), shape = k, rate = b)))
  }
  simulate_paths(t, n_traj, one_path, model, parameters, A)
}
