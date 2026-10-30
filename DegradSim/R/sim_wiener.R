#' Simulation of the Wiener process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' a Wiener process with drift (a Levy process) or of a time-scaled Wiener
#' process, obtained with a time-scaling function \eqn{A}.
#'
#' @section Model:
#' \eqn{B} denotes a standard Brownian motion and \eqn{\mathcal{N}(m, v)} the normal
#' distribution with mean \eqn{m} and variance \eqn{v}.
#'
#' **Wiener process** (`A = NULL`, default). Let \eqn{\mu} be
#' real and \eqn{\sigma > 0}. The process is
#' \deqn{X_t = \mu t + \sigma B_t, \qquad t \ge 0,}
#' where \eqn{\mu} is the drift parameter (expected to be positive for a
#' degradation) and \eqn{\sigma^2} the variance (or volatility) parameter. It
#' has independent increments \eqn{X_t - X_s \sim \mathcal{N}(\mu (t - s), \sigma^2 (t - s))},
#' \eqn{E[X_t] = \mu t} and \eqn{Var(X_t) = \sigma^2 t}. It is a Levy
#' process with continuous, non-monotone trajectories.
#'
#' **Time-scaled Wiener process** (`A` given). With the time-scaling function \eqn{A}
#' defined in [DegradSim-package],
#' \deqn{X_t = A(t) + \sigma B_{A(t)}, \qquad t \ge 0,}
#' that is \eqn{X_t = Y_{A(t)}} with \eqn{Y_u = u + \sigma B_u}: the drift
#' is set to \eqn{\mu = 1} (identifiability, see [DegradSim-package]), as in
#' Cha and Mercier (2019, Example 1). Its increments are independent,
#' \eqn{X_t - X_s \sim \mathcal{N}(A(t) - A(s), \sigma^2 (A(t) - A(s)))},
#' \eqn{E[X_t] = A(t)} and \eqn{Var(X_t) = \sigma^2 A(t)}.
#'
#' @section Exact simulation:
#' Random walk approximation, described in [DegradSim-package] (Kahle et
#' al., 2016, algorithm 1.7): the increments
#' \eqn{\Delta X_i = \mu \Delta t_i + \sigma \sqrt{\Delta t_i}\, Z_i} or
#' \eqn{\Delta X_i = \Delta A_i + \sigma \sqrt{\Delta A_i}\, Z_i}, with
#' \eqn{Z_i} independent \eqn{\mathcal{N}(0, 1)} drawn with [stats::rnorm()], are
#' cumulated.
#'
#' @inheritParams sim_gamma
#' @param mu Drift of the Wiener process, real: mean degradation per unit
#'   of time. Not given when `A` is given (the mean of \eqn{X_t} is then
#'   \eqn{A(t)}).
#' @param sigma Diffusion coefficient, \eqn{\sigma > 0}: standard deviation
#'   (not variance) of \eqn{X_t - X_s} for \eqn{t - s = 1} (Wiener process)
#'   or \eqn{A(t) - A(s) = 1} (time-scaled Wiener process).
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
#' @seealso
#' * [sim_gamma()], [sim_ig()]: the other time-scaled processes;
#' * [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()].
#'
#' @examples
#' # 1. One trajectory with 55 regular inspections up to 10 (step 10/55)
#' t <- inspection_times("regular", n = 55, tmax = 10)$time
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
#'
#' # 3. Irregular inspections, frequent at the beginning and at the end
#' plan <- inspection_times("inc-dec", n = 50, n_traj = 5, tmax = 10)
#' x <- sim_wiener(plan, mu = 0.5, sigma = 1.5)
#' plot_traj(x)
#'
#' # 4. Different inspection times for each trajectory: a list of data frames
#' x <- sim_wiener(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), mu = 0.5, sigma = 1.5)
#' x[[2]]
#' x[[2]]$X
#'
#' # 5. Time-scaled Wiener process X_t = A(t) + sigma B_{A(t)} with
#' #    A(t) = 0.5 t^1.5 (Cha and Mercier, 2019, Example 1)
#' x <- sim_wiener(t, sigma = 1.5, A = function(t) 0.5 * t^1.5, n_traj = 10)
#' x
#' plot_traj(x, "Time-scaled Wiener process")
#'
#' # 6. With A(t) = mu t, the time-scaled Wiener process with sigma is the
#' #    Wiener process with drift mu and diffusion coefficient sigma sqrt(mu)
#' set.seed(2)
#' y1 <- sim_wiener(t, sigma = 1.5, A = function(t) 0.5 * t)
#' set.seed(2)
#' y2 <- sim_wiener(t, mu = 0.5, sigma = 1.5 * sqrt(0.5))
#' all.equal(unclass(y1)[, 2], unclass(y2)[, 2])
#'
#' # 7. All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#'
#' # 8. Histogram of the increments, with their density in red. Constant
#' #    step h = 0.5: the increments are i.i.d. N(mu h, sigma^2 h)
#' t8 <- inspection_times("regular", n = 50, tmax = 25)$time
#' x <- sim_wiener(t8, mu = 0.5, sigma = 1.5, n_traj = 20)
#' dx <- traj_df(x, increments = TRUE)$dX
#' dx <- dx[!is.na(dx)]
#' ggplot2::ggplot(data.frame(dx), ggplot2::aes(dx)) +
#'   ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
#'                           bins = 30, boundary = 0, fill = "grey80", colour = "grey40") +
#'   ggplot2::stat_function(fun = function(y) dnorm(y, 0.5 * 0.5, 1.5 * sqrt(0.5)),
#'                          colour = "red") +
#'   ggplot2::labs(x = "dX", y = "Density", title = "Increments of the Wiener process")
#' @export
sim_wiener <- function(t, mu, sigma, A = NULL, method = "RWA", type = NULL, n_traj = 1) {
  check_sim_method(method, type)
  if (missing(sigma) || !is.numeric(sigma) || length(sigma) != 1 || !(sigma > 0)) stop("sigma must be a positive number")
  if (is.null(A)) {
    if (missing(mu)) stop("mu is missing: give mu (Wiener process) or A (time-scaled Wiener process)")
    if (!is.numeric(mu) || length(mu) != 1 || !is.finite(mu)) stop("mu must be a real number")
    clock <- function(t) diff(t)
    drift <- mu
    model <- "Wiener"
    parameters <- c(mu = mu, sigma = sigma)
  } else {
    if (!missing(mu)) stop("give mu or A, not both: with A, the mean of X_t is A(t) (mu = 1)")
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
