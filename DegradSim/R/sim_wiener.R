#' Simulation of the Wiener process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' a homogeneous Wiener process with drift (a Lévy process) or of a
#' time-scaled Wiener process, obtained by a time change \eqn{A}.
#'
#' @section Model:
#' \eqn{B} denotes a standard Brownian motion and \eqn{N(m, v)} the normal
#' distribution with mean \eqn{m} and variance \eqn{v}.
#'
#' **Homogeneous Wiener process** (`A = NULL`, default). Let \eqn{\mu} be
#' real and \eqn{\sigma > 0}. The process is
#' \deqn{X(t) = \mu t + \sigma B(t), \qquad t \ge 0,}
#' where \eqn{\mu} is the drift parameter (expected to be positive for a
#' degradation) and \eqn{\sigma^2} the variance (or volatility) parameter. It
#' has independent increments \eqn{X(t) - X(s) \sim N(\mu (t - s), \sigma^2 (t - s))},
#' \eqn{E[X(t)] = \mu t} and \eqn{Var(X(t)) = \sigma^2 t}. It is a Lévy
#' process with continuous, non-monotone trajectories.
#'
#' **Time-scaled Wiener process** (`A` given). Let \eqn{A} be a time change: a
#' non-decreasing function, right-continuous with left limits, with
#' \eqn{A(0) = 0} and \eqn{A(t) \to +\infty} as \eqn{t \to +\infty}. The
#' process is
#' \deqn{X(t) = A(t) + \sigma B(A(t)), \qquad t \ge 0,}
#' that is \eqn{X(t) = Y(A(t))} with \eqn{Y(u) = u + \sigma B(u)} (Cha and
#' Mercier, 2019, Example 1, with \eqn{\mu = 1}). Its increments are
#' independent, \eqn{X(t) - X(s) \sim N(A(t) - A(s), \sigma^2 (A(t) - A(s)))},
#' \eqn{E[X(t)] = A(t)} and \eqn{Var(X(t)) = \sigma^2 A(t)}.
#'
#' **Why \eqn{\mu = 1} in the time-scaled case.** With a drift \eqn{\mu}, the
#' processes \eqn{\mu A(t) + \sigma B(A(t))} and
#' \eqn{(\mu / c) (c A(t)) + (\sigma / \sqrt{c}) B(c A(t))} have the same
#' distribution for every \eqn{c > 0}: \eqn{\mu}, \eqn{\sigma} and \eqn{A}
#' could not be identified separately. The drift is therefore contained in
#' \eqn{A}. With \eqn{A(t) = \mu t} (\eqn{\mu > 0}), the time-scaled process
#' is the homogeneous process with drift \eqn{\mu} and variance parameter
#' \eqn{\sigma^2 \mu}: the homogeneous process with parameters
#' \eqn{(\mu, \sigma)} is obtained with \eqn{A(t) = \mu t} and
#' \eqn{\sigma / \sqrt{\mu}}.
#'
#' @section Exact simulation at the inspection times (random walk approximation):
#' The increments over disjoint intervals are independent normal variables,
#' so the values at the inspection times \eqn{0 = t_0 < t_1 < \dots < t_n} are
#' simulated exactly by the random walk approximation (RWA) of Kahle et al.
#' (2016, section 1.3.1, algorithm 1.7):
#' 1. \eqn{\Delta X_i = \mu \Delta t_i + \sigma \sqrt{\Delta t_i}\, Z_i}
#'    (homogeneous case) or
#'    \eqn{\Delta X_i = \Delta A_i + \sigma \sqrt{\Delta A_i}\, Z_i}
#'    (time-scaled case), with \eqn{\Delta t_i = t_i - t_{i-1}},
#'    \eqn{\Delta A_i = A(t_i) - A(t_{i-1})} and \eqn{Z_i} independent
#'    \eqn{N(0, 1)} drawn with [stats::rnorm()];
#' 2. \eqn{X(t_0) = 0} and \eqn{X(t_i) = \sum_{j \le i} \Delta X_j}.
#'
#' The values at the inspection times have exactly the distribution of the
#' process; only the path between two inspection times is not simulated.
#'
#' @inheritSection sim_gamma Conditions on the inspection times and on A
#'
#' @inheritParams sim_gamma
#' @param mu Drift of the homogeneous process, real: mean degradation per unit
#'   of time. Not given when `A` is given (the mean of \eqn{X(t)} is then
#'   \eqn{A(t)}).
#' @param sigma Diffusion coefficient, \eqn{\sigma > 0}: standard deviation
#'   (not variance) of \eqn{X(t) - X(s)} for \eqn{t - s = 1} (homogeneous
#'   case) or \eqn{A(t) - A(s) = 1} (time-scaled case).
#' @param A `NULL` (default) for the homogeneous process, or the time change
#'   for the time-scaled process: a vectorized R function of `t`,
#'   non-decreasing with \eqn{A(0) = 0}, for instance
#'   `function(t) 0.5 * t^1.5`.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Kahle, W., Mercier, S. and Paroissin, C. (2016). *Degradation Processes in
#' Reliability*. ISTE and Wiley, London. ISBN 978-1-84821-888-8.
#'
#' Cha, J. H. and Mercier, S. (2019). Transformed Lévy processes as
#' state-dependent wear models. *Advances in Applied Probability*, 51(2),
#' 468-486. \doi{10.1017/apr.2019.21}
#'
#' Whitmore, G. A. (1995). Estimating degradation by a Wiener diffusion process
#' subject to measurement error. *Lifetime Data Analysis*, 1, 307-319.
#'
#' @seealso [sim_gamma()] and [sim_ig()] for the other processes with a time
#'   change, [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()]
#'
#' @examples
#' # 1. Homogeneous process: one trajectory inspected every 0.5 time unit
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_wiener(t, mu = 0.5, sigma = 1.5)
#' x
#' plot_traj(x)
#'
#' # 2. Ten trajectories with the same inspection times
#' x <- sim_wiener(t, mu = 0.5, sigma = 1.5, n_traj = 10)
#' plot_traj(x)
#'
#' # Same times for all the trajectories: a matrix with the column t,
#' # then one column per trajectory (X1, X2, ...)
#' head(x[, 1:4])
#' x[, "X2"]
#'
#' # 3. Irregular inspections, frequent at the beginning and at the end
#' plan <- inspection_times("dec-inc", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_wiener(plan, mu = 0.5, sigma = 1.5)
#' plot_traj(x)
#'
#' # 4. Different inspection times for each trajectory: a list of data frames
#' x <- sim_wiener(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), mu = 0.5, sigma = 1.5)
#' x[[2]]
#' x[[2]]$X
#'
#' # 5. Time-scaled Wiener process X(t) = A(t) + sigma B(A(t)) with
#' #    A(t) = 0.5 t^1.5 (Cha and Mercier, 2019, Example 1)
#' x <- sim_wiener(t, sigma = 1.5, A = function(t) 0.5 * t^1.5, n_traj = 10)
#' x
#' plot_traj(x, "Time-scaled Wiener process")
#'
#' # 6. The homogeneous process (mu, sigma) is A(t) = mu t with sigma / sqrt(mu)
#' set.seed(2)
#' y1 <- sim_wiener(t, mu = 0.5, sigma = 1.5)
#' set.seed(2)
#' y2 <- sim_wiener(t, sigma = 1.5 / sqrt(0.5), A = function(t) 0.5 * t)
#' all.equal(unclass(y1)[, 2], unclass(y2)[, 2])
#'
#' # 7. All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#' @export
sim_wiener <- function(t, mu, sigma, A = NULL, n_traj = 1) {
  if (missing(sigma) || !is.numeric(sigma) || length(sigma) != 1 || !(sigma > 0)) stop("sigma must be a positive number")
  if (is.null(A)) {
    if (missing(mu)) stop("mu is missing: give mu (homogeneous process) or A (time-scaled process)")
    if (!is.numeric(mu) || length(mu) != 1 || !is.finite(mu)) stop("mu must be a real number")
    clock <- function(t) diff(t)
    drift <- mu
    model <- "Wiener"
    parameters <- c(mu = mu, sigma = sigma)
  } else {
    if (!missing(mu)) stop("give mu or A, not both: with A, the mean of X(t) is A(t) (mu = 1)")
    clock <- function(t) time_change(A, t)
    drift <- 1
    model <- "time-scaled Wiener"
    parameters <- c(sigma = sigma)
  }
  one_path <- function(t) {
    dA <- clock(t)
    c(0, cumsum(drift * dA + sigma * sqrt(dA) * stats::rnorm(length(dA))))
  }
  simulate_paths(t, n_traj, one_path, model, parameters, A)
}
