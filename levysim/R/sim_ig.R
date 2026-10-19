#' Simulation of the inverse Gaussian process
#'
#' Simulates exactly one or several trajectories of an inverse Gaussian (IG)
#' process, homogeneous or not, at given inspection times.
#'
#' @section Model:
#' \eqn{IG(m, \lambda)} denotes the inverse Gaussian distribution with mean
#' \eqn{m} and shape \eqn{\lambda}, of density
#' \eqn{\sqrt{\lambda / (2 \pi x^3)} \exp(-\lambda (x - m)^2 / (2 m^2 x))},
#' \eqn{x > 0}. Let \eqn{\Lambda} be continuous, strictly increasing, with
#' \eqn{\Lambda(0) = 0}, \eqn{\mu > 0} and \eqn{\eta > 0}. The IG process has
#' \eqn{X(0) = 0}, independent increments and
#' \deqn{X(t) - X(s) \sim IG(\mu \Delta\Lambda, \eta \Delta\Lambda^2).}
#' Hence \eqn{E[X(t)] = \mu \Lambda(t)} and
#' \eqn{Var(X(t)) = \mu^3 \Lambda(t) / \eta}. The process is increasing and
#' purely made of jumps. With \eqn{\Lambda(t) = t} (default), it is the
#' homogeneous IG process, a Lévy process, with the parametrization of Wang and
#' Xu (2010).
#'
#' @section Exact simulation:
#' Each increment \eqn{IG(m, \lambda)}, with \eqn{m = \mu \Delta\Lambda_i} and
#' \eqn{\lambda = \eta \Delta\Lambda_i^2}, is drawn by the method of Michael,
#' Schucany and Haas (1976), based on
#' \eqn{\lambda (X - m)^2 / (m^2 X) \sim \chi^2_1}:
#' 1. \eqn{\nu = Z^2}, \eqn{Z \sim N(0, 1)};
#' 2. \eqn{x_1 = m + m^2 \nu / (2\lambda) - m \sqrt{4 m \lambda \nu + m^2 \nu^2} / (2\lambda)};
#' 3. return \eqn{x_1} with probability \eqn{m / (m + x_1)}, otherwise
#'    \eqn{m^2 / x_1}.
#'
#' The increments are independent; the trajectory is their cumulative sum.
#'
#' @inheritParams sim_gamma
#' @param mu Mean degradation per unit of \eqn{\Lambda}, \eqn{\mu > 0}.
#' @param eta Shape parameter, \eqn{\eta > 0}: the larger \eqn{\eta}, the less
#'   variable the increments (\eqn{Var(X(t)) = \mu^3 \Lambda(t) / \eta}).
#' @param Lambda Function of time: an R function of `t`, continuous and
#'   strictly increasing with \eqn{\Lambda(0) = 0}, for instance
#'   `function(t) t^1.5`. The default `function(t) t` gives the homogeneous
#'   process.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
#' variates using transformations with multiple roots. *The American
#' Statistician*, 30(2), 88-90. \doi{10.1080/00031305.1976.10479147}
#'
#' Wang, X. and Xu, D. (2010). An inverse Gaussian process model for
#' degradation data. *Technometrics*, 52(2), 188-197.
#' \doi{10.1198/TECH.2009.08197}
#'
#' @seealso [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()]
#'
#' @examples
#' # One trajectory inspected every 0.5 time unit up to 10
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_ig(t, mu = 0.5, eta = 2)
#' x
#' plot_traj(x)
#'
#' # Ten trajectories with the same inspection times
#' x <- sim_ig(t, mu = 0.5, eta = 2, n_traj = 10)
#' plot_traj(x)
#'
#' # Same times for all the trajectories: a matrix with the column t,
#' # then one column per trajectory (X1, X2, ...)
#' head(x[, 1:4])
#' x[, "X2"]
#'
#' # Irregular inspections, less and less frequent over time
#' plan <- inspection_times("decreasing", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_ig(plan, mu = 0.5, eta = 2)
#' plot_traj(x)
#'
#' # Different inspection times for each trajectory: a list of data frames
#' x <- sim_ig(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), mu = 0.5, eta = 2)
#' x[[2]]
#'
#' # Non-homogeneous process with Lambda(t) = t^1.5
#' x <- sim_ig(t, mu = 0.5, eta = 2, Lambda = function(t) t^1.5, n_traj = 10)
#' plot_traj(x, "Non-homogeneous inverse Gaussian process")
#'
#' # All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#' @export
sim_ig <- function(t, mu, eta, Lambda = function(t) t, n_traj = 1) {
  if (mu <= 0 || eta <= 0) stop("mu and eta must be positive")
  one_path <- function(t) {
    dL <- diff(Lambda(t))
    if (any(dL <= 0)) stop("Lambda must be strictly increasing")
    c(0, cumsum(rig(length(dL), mu * dL, eta * dL^2)))
  }
  simulate_paths(t, n_traj, one_path, "inverse Gaussian", c(mu = mu, eta = eta))
}

rig <- function(n, m, lambda) {
  nu <- stats::rnorm(n)^2
  x1 <- m + m^2 * nu / (2 * lambda) - m / (2 * lambda) * sqrt(4 * m * lambda * nu + m^2 * nu^2)
  ifelse(stats::runif(n) <= m / (m + x1), x1, m^2 / x1)
}
