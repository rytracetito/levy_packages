#' Simulation of the normal inverse Gaussian process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' a normal inverse Gaussian (NIG) process, a Levy process.
#'
#' @section Model:
#' The \eqn{NIG(\alpha, \beta, \delta, \mu)} distribution, with
#' \eqn{0 \le |\beta| < \alpha}, \eqn{\delta > 0} and \eqn{\mu} real, has
#' density
#' \deqn{f(x) = \frac{\alpha \delta\, K_1\big(\alpha \sqrt{\delta^2 + (x - \mu)^2}\big)}
#'   {\pi \sqrt{\delta^2 + (x - \mu)^2}} \exp\big(\delta \gamma + \beta (x - \mu)\big),
#'   \qquad x \in \mathbb{R},}
#' where \eqn{\gamma = \sqrt{\alpha^2 - \beta^2}} and \eqn{K_1} is the modified
#' Bessel function of the second kind (Barndorff-Nielsen, 1997). It is
#' infinitely divisible: for independent
#' \eqn{Y_1 \sim NIG(\alpha, \beta, \delta_1, \mu_1)} and
#' \eqn{Y_2 \sim NIG(\alpha, \beta, \delta_2, \mu_2)},
#' \deqn{Y_1 + Y_2 \sim NIG(\alpha, \beta, \delta_1 + \delta_2, \mu_1 + \mu_2).}
#'
#' The NIG process satisfies \eqn{X_0 = 0}, has independent increments, and for
#' \eqn{0 \le s < t},
#' \deqn{X_t - X_s \sim NIG\big(\alpha, \beta, \delta (t - s), \mu (t - s)\big).}
#' Hence \eqn{E[X_t] = (\mu + \delta \beta / \gamma)\, t} and
#' \eqn{Var(X_t) = \delta \alpha^2 t / \gamma^3}. It is a Levy process whose
#' trajectories are not monotone, purely made of jumps, of infinite
#' variation: \eqn{\alpha} sets the tails, \eqn{\beta} the asymmetry,
#' \eqn{\delta} the scale and \eqn{\mu} the location. It is a subordinated
#' Brownian motion (Hoyle, 2010, section 2.10):
#' \eqn{X_t = \mu t + \beta Z_t + W_{Z_t}}, where \eqn{W} is a standard
#' Brownian motion and \eqn{Z} an independent inverse Gaussian process with
#' \eqn{Z_t - Z_s \sim IG(\delta (t - s) / \gamma, \delta^2 (t - s)^2)}
#' (mean, shape), as the variance-gamma process is a Brownian motion
#' subordinated by a gamma process.
#'
#' @section Exact simulation:
#' Random walk approximation, described in [DegradSim-package]: with
#' \eqn{\Delta t_i = t_i - t_{i-1}}, each increment is drawn by subordination,
#' then the increments are cumulated:
#' 1. draw \eqn{\Delta Z_i \sim IG(\delta \Delta t_i / \gamma, (\delta \Delta t_i)^2)}
#'    (mean, shape) by the method of Michael, Schucany and Haas (1976), as in
#'    [sim_ig()];
#' 2. compute \eqn{\Delta X_i = \mu \Delta t_i + \beta \Delta Z_i + \sqrt{\Delta Z_i}\, \xi_i},
#'    with \eqn{\xi_i} independent \eqn{\mathcal{N}(0, 1)}.
#'
#' Conditionally on \eqn{\Delta Z_i}, the increment is
#' \eqn{\mathcal{N}(\mu \Delta t_i + \beta \Delta Z_i, \Delta Z_i)}, and its
#' mixture over the IG distribution of \eqn{\Delta Z_i} is the
#' \eqn{NIG(\alpha, \beta, \delta \Delta t_i, \mu \Delta t_i)} distribution.
#'
#' @inheritParams sim_gamma
#' @param alpha Tail parameter, \eqn{\alpha > 0}: the larger \eqn{\alpha}, the
#'   lighter the tails.
#' @param beta Asymmetry parameter, \eqn{|\beta| < \alpha}: \eqn{\beta = 0}
#'   gives symmetric increments, \eqn{\beta > 0} increments skewed to the right.
#' @param delta Scale parameter per unit of time, \eqn{\delta > 0}.
#' @param mu Location parameter per unit of time, real.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Barndorff-Nielsen, O. E. (1997). Normal inverse Gaussian distributions and
#' stochastic volatility modelling. *Scandinavian Journal of Statistics*,
#' 24(1), 1-13. \doi{10.1111/1467-9469.00045}
#'
#' Hoyle, A. E. V. (2010). *Information-Based Models for Finance and
#' Insurance*. PhD thesis, Imperial College London. arXiv:1010.0829.
#'
#' Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
#' variates using transformations with multiple roots. *The American
#' Statistician*, 30(2), 88-90.\cr
#' \doi{10.1080/00031305.1976.10479147}
#'
#' @seealso [sim_vg()], [plot_traj()], [traj_df()], [inspection_times()],
#'   [plot_plan()]
#'
#' @examples
#' # 1. One trajectory with 200 regular inspections up to 10 (step 0.05)
#' t <- inspection_times("regular", n = 200, tmax = 10)$time
#' set.seed(1)
#' x <- sim_nig(t, alpha = 2, beta = 0.5, delta = 1, mu = 0.2)
#' x
#' plot_traj(x)
#'
#' # 2. Ten trajectories with the same inspection times
#' x <- sim_nig(t, alpha = 2, beta = 0.5, delta = 1, mu = 0.2, n_traj = 10)
#' plot_traj(x)
#'
#' # Same times for all the trajectories: a matrix with the column t,
#' # then one column per trajectory (X1, X2, ...)
#' head(x[, 1:4])
#' x[, "X2"]
#'
#' # 3. Irregular inspections, more and more frequent over time
#' plan <- inspection_times("decreasing", n = 50, n_traj = 5, tmax = 10)
#' x <- sim_nig(plan, alpha = 2, beta = 0.5, delta = 1, mu = 0.2)
#' plot_traj(x)
#'
#' # 4. Different inspection times for each trajectory: a list of data frames
#' x <- sim_nig(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)),
#'              alpha = 2, beta = 0.5, delta = 1, mu = 0.2)
#' x[[2]]
#' x[[2]]$X
#'
#' # 5. Heavier tails (small alpha) and symmetric increments (beta = 0)
#' x <- sim_nig(t, alpha = 0.5, beta = 0, delta = 1, mu = 0.2, n_traj = 10)
#' plot_traj(x, "NIG process with alpha = 0.5")
#'
#' # 6. All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#'
#' # 7. Histogram of the increments, with the density of the increments in red.
#' #    Constant step h = 1: the increments are i.i.d.
#' #    with distribution NIG(alpha, beta, delta h, mu h)
#' t7 <- inspection_times("regular", n = 50, tmax = 50)$time
#' x <- sim_nig(t7, alpha = 2, beta = 0.5, delta = 1, mu = 0.2, n_traj = 20)
#' dx <- traj_df(x, increments = TRUE)$dX
#' dx <- dx[!is.na(dx)]
#' dnig <- function(y, alpha = 2, beta = 0.5, delta = 1, mu = 0.2) {
#'   gamma <- sqrt(alpha^2 - beta^2)
#'   q <- sqrt(delta^2 + (y - mu)^2)
#'   alpha * delta / (pi * q) * besselK(alpha * q, 1, expon.scaled = TRUE) *
#'     exp(delta * gamma + beta * (y - mu) - alpha * q)
#' }
#' ggplot2::ggplot(data.frame(dx), ggplot2::aes(dx)) +
#'   ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
#'                           bins = 30, boundary = 0, fill = "grey80", colour = "grey40") +
#'   ggplot2::stat_function(fun = dnig, colour = "red") +
#'   ggplot2::labs(x = "dX", y = "Density", title = "Increments of the NIG process")
#' @export
sim_nig <- function(t, alpha, beta, delta, mu, n_traj = 1) {
  if (!is.numeric(alpha) || length(alpha) != 1 || !(alpha > 0)) stop("alpha must be a positive number")
  if (!is.numeric(beta) || length(beta) != 1 || !(abs(beta) < alpha)) stop("beta must satisfy |beta| < alpha")
  if (!is.numeric(delta) || length(delta) != 1 || !(delta > 0)) stop("delta must be a positive number")
  if (!is.numeric(mu) || length(mu) != 1 || !is.finite(mu)) stop("mu must be a real number")
  gamma <- sqrt(alpha^2 - beta^2)
  one_path <- function(t) {
    dt <- diff(t)
    dZ <- rig(length(dt), delta * dt / gamma, (delta * dt)^2)
    c(0, cumsum(mu * dt + beta * dZ + sqrt(dZ) * stats::rnorm(length(dt))))
  }
  simulate_paths(t, n_traj, one_path, "NIG", c(alpha = alpha, beta = beta, delta = delta, mu = mu))
}
