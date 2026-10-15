#' Simulation of the Meixner process
#'
#' Simulates exactly one or several trajectories of a Meixner process,
#' homogeneous or not, at given inspection times.
#'
#' @section Model:
#' The \eqn{Meixner(a, b, d, m)} distribution, \eqn{a > 0},
#' \eqn{-\pi < b < \pi}, \eqn{d > 0}, \eqn{m} real, has density
#' \deqn{f(x) = \frac{(2 \cos(b/2))^{2d}}{2 a \pi\, \Gamma(2d)}
#'   \exp\Big(\frac{b (x - m)}{a}\Big) \Big|\Gamma\Big(d + i \frac{x - m}{a}\Big)\Big|^2.}
#' Let \eqn{\Lambda} be continuous, strictly increasing, with
#' \eqn{\Lambda(0) = 0}. The Meixner process has \eqn{X(0) = 0}, independent
#' increments and
#' \deqn{X(t) - X(s) \sim Meixner(a, b, d \Delta\Lambda, m \Delta\Lambda).}
#' Hence \eqn{E[X(t)] = (m + a d \tan(b/2)) \Lambda(t)} and
#' \eqn{Var(X(t)) = a^2 d \Lambda(t) / (2 \cos^2(b/2))}. The trajectories are
#' not monotone, purely made of jumps, of infinite variation; \eqn{b} sets the
#' asymmetry and \eqn{a} is a scale parameter. With \eqn{\Lambda(t) = t}
#' (default), it is the homogeneous Meixner process, a Lévy process.
#'
#' @section Exact simulation:
#' There is no Meixner generator in `stats`. Each increment is drawn from
#' \eqn{Meixner(a, b, d \Delta\Lambda_i, m \Delta\Lambda_i)} by numerical
#' inversion of its distribution function (method PINV of the package
#' `Runuran`, from the exact density, error in \eqn{u} of order
#' \eqn{10^{-10}}), then the increments are cumulated. One generator is built
#' for each distinct length \eqn{\Delta\Lambda_i} and reused for all the
#' trajectories.
#'
#' @inheritParams sim_gamma
#' @param a Scale parameter, \eqn{a > 0}.
#' @param b Asymmetry parameter, \eqn{-\pi < b < \pi}: \eqn{b = 0} gives
#'   symmetric increments, \eqn{b > 0} increments skewed to the right.
#' @param d Shape parameter per unit of \eqn{\Lambda}, \eqn{d > 0}.
#' @param m Location parameter per unit of \eqn{\Lambda}, real.
#' @param Lambda Function of time: an R function of `t`, continuous and
#'   strictly increasing with \eqn{\Lambda(0) = 0}, for instance
#'   `function(t) t^1.5`. The default `function(t) t` gives the homogeneous
#'   process.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Schoutens, W. (2003). *Lévy Processes in Finance: Pricing Financial
#' Derivatives*. Wiley.
#'
#' Schoutens, W. and Teugels, J. L. (1998). Lévy processes, polynomials and
#' martingales. *Communications in Statistics. Stochastic Models*, 14(1-2),
#' 335-349.
#'
#' @seealso [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()]
#'
#' @examples
#' # One trajectory inspected every 0.5 time unit up to 10
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' x
#' plot_traj(x)
#'
#' # Ten trajectories with the same inspection times
#' x <- sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2, n_traj = 10)
#' plot_traj(x)
#'
#' # Irregular inspections, frequent in the middle of the period
#' plan <- inspection_times("inc-dec", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_meixner(plan, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' plot_traj(x)
#'
#' # Different inspection times for each trajectory
#' x <- sim_meixner(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), a = 1, b = 0.3, d = 1.5, m = 0.2)
#' x[[2]]
#'
#' # Non-homogeneous process with Lambda(t) = t^1.5
#' x <- sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2, Lambda = function(t) t^1.5, n_traj = 10)
#' plot_traj(x, "Non-homogeneous Meixner process")
#'
#' # All the simulated data in a single data frame
#' head(traj_df(x))
#' @export
sim_meixner <- function(t, a, b, d, m, Lambda = function(t) t, n_traj = 1) {
  if (a <= 0 || d <= 0) stop("a and d must be positive")
  if (abs(b) >= pi) stop("b must be in (-pi, pi)")
  generators <- new.env()
  draw <- function(s) {
    key <- format(s, digits = 15)
    if (is.null(generators[[key]])) {
      loi <- Runuran::udmeixner(alpha = a, beta = b, delta = d * s, mu = m * s)
      generators[[key]] <- Runuran::pinvd.new(loi)
    }
    Runuran::ur(generators[[key]], 1)
  }
  one_path <- function(t) {
    dL <- diff(Lambda(t))
    if (any(dL <= 0)) stop("Lambda must be strictly increasing")
    c(0, cumsum(vapply(dL, draw, 0)))
  }
  simulate_paths(t, n_traj, one_path, "Meixner", c(a = a, b = b, d = d, m = m))
}
