#' Simulation of the Wiener process
#'
#' Simulates exactly one or several trajectories of a Wiener process with drift,
#' homogeneous or with time transformations, at given inspection times.
#'
#' @section Model:
#' Let \eqn{\mu} be real, \eqn{\sigma > 0}, \eqn{\Lambda} a continuous function
#' with \eqn{\Lambda(0) = 0}, \eqn{\tau} a continuous strictly increasing
#' function with \eqn{\tau(0) = 0}, and \eqn{B} a standard Brownian motion. The
#' process is
#' \deqn{X(t) = \mu \Lambda(t) + \sigma B(\tau(t)),}
#' with independent increments
#' \eqn{X(t) - X(s) \sim N(\mu \Delta\Lambda, \sigma^2 \Delta\tau)}. Hence
#' \eqn{E[X(t)] = \mu \Lambda(t)} and \eqn{Var(X(t)) = \sigma^2 \tau(t)}. The
#' trajectories are continuous and not monotone. With
#' \eqn{\Lambda(t) = \tau(t) = t} (default), it is the homogeneous Wiener
#' process \eqn{\mu t + \sigma B(t)}, a Lévy process; \eqn{\Lambda = \tau} gives
#' the usual Wiener process with transformed time.
#'
#' @section Exact simulation:
#' \eqn{\Delta X_i = \mu \Delta\Lambda_i + \sigma \sqrt{\Delta\tau_i}\, Z_i},
#' with \eqn{Z_i} independent \eqn{N(0, 1)}, then cumulative sum.
#'
#' @inheritParams sim_gamma
#' @param mu Drift, real: mean degradation per unit of \eqn{\Lambda}.
#' @param sigma Diffusion coefficient, \eqn{\sigma > 0}: standard deviation
#'   (not variance) per square root of unit of \eqn{\tau}.
#' @param Lambda Function of time for the drift: an R function of `t`,
#'   continuous with \eqn{\Lambda(0) = 0}, for instance `function(t) t^1.5`.
#'   The default `function(t) t` gives a linear mean.
#' @param tau Function of time for the diffusion: an R function of `t`,
#'   strictly increasing with \eqn{\tau(0) = 0}. The default `function(t) t`
#'   gives the usual Brownian motion.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Whitmore, G. A. (1995). Estimating degradation by a Wiener diffusion process
#' subject to measurement error. *Lifetime Data Analysis*, 1, 307-319.
#'
#' @seealso [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()]
#'
#' @examples
#' # One trajectory inspected every 0.5 time unit up to 10
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_wiener(t, mu = 0.5, sigma = 1.5)
#' x
#' plot_traj(x)
#'
#' # Ten trajectories with the same inspection times
#' x <- sim_wiener(t, mu = 0.5, sigma = 1.5, n_traj = 10)
#' plot_traj(x)
#'
#' # Same times for all the trajectories: a matrix with the column t,
#' # then one column per trajectory (X1, X2, ...)
#' head(x[, 1:4])
#' x[, "X2"]
#'
#' # Irregular inspections, frequent at the start and at the end
#' plan <- inspection_times("dec-inc", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_wiener(plan, mu = 0.5, sigma = 1.5)
#' plot_traj(x)
#'
#' # Different inspection times for each trajectory: a list of data frames
#' x <- sim_wiener(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), mu = 0.5, sigma = 1.5)
#' x[[2]]
#'
#' # Transformed time: Lambda(t) = t^1.5 and tau(t) = t^1.2
#' x <- sim_wiener(t, mu = 0.5, sigma = 1.5, Lambda = function(t) t^1.5,
#'                 tau = function(t) t^1.2, n_traj = 10)
#' plot_traj(x, "Wiener process with transformed time")
#'
#' # All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#' @export
sim_wiener <- function(t, mu, sigma, Lambda = function(t) t, tau = function(t) t, n_traj = 1) {
  if (sigma <= 0) stop("sigma must be positive")
  one_path <- function(t) {
    dL <- diff(Lambda(t))
    dtau <- diff(tau(t))
    if (any(dtau <= 0)) stop("tau must be strictly increasing")
    c(0, cumsum(mu * dL + sigma * sqrt(dtau) * stats::rnorm(length(dL))))
  }
  simulate_paths(t, n_traj, one_path, "Wiener", c(mu = mu, sigma = sigma))
}
