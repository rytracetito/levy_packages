#' Brownian bridge
#'
#' Draws the value at time \eqn{t} of a Wiener process with drift known at two
#' times \eqn{t_A < t < t_B}, from the distribution of the Brownian bridge.
#' This function is used by [wiener_resampling()] at each point of the grid.
#'
#' @section Model:
#' The Wiener process with drift \eqn{\mu} and volatility \eqn{\sigma > 0} is
#' \eqn{X_t = \mu t + \sigma B_t}, with \eqn{B} a standard Brownian motion:
#' \eqn{X_0 = 0}, independent increments and
#' \eqn{X_t - X_s \sim \mathcal{N}(\mu (t - s), \sigma^2 (t - s))}. Its
#' trajectories are not monotone. The vector of its parameters is
#' \eqn{\Theta = (\mu, \sigma)}.
#'
#' @section Distribution of the bridge:
#' Given \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B},
#' \deqn{X_t \sim \mathcal{N}\Big(x_A + (x_B - x_A) \frac{t - t_A}{t_B - t_A},\ \sigma^2 \frac{(t - t_A)(t_B - t)}{t_B - t_A}\Big),}
#' which is the formula of the section Bridge of a Levy process of
#' [levyBridge-package] with Gaussian densities. This distribution does not
#' depend on the drift \eqn{\mu}.
#'
#' @inheritParams gamma_bridge
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(mu = , sigma = )`, with the volatility \eqn{\sigma}
#'   (not the variance). The function `wiener_mle()` of the package `levyfit`
#'   returns `mu` and `sigma2`: use
#'   `c(mu = est[["mu"]], sigma = sqrt(est[["sigma2"]]))`. The bridge uses
#'   only `sigma`.
#'
#' @inherit gamma_bridge return
#'
#' @seealso [wiener_resampling()] to resample trajectories on a regular grid;
#'   [gamma_bridge()], [ig_bridge()], [vg_bridge()], [meixner_bridge()].
#'
#' @examples
#' # 1. One draw of X_2, given X_1 = 3 and X_4 = 6 (mu = 0.5, sigma = 1.5)
#' set.seed(1)
#' wiener_bridge(t_A = 1, x_A = 3, t_B = 4, x_B = 6, t = 2, theta = c(mu = 0.5, sigma = 1.5))
#'
#' # 2. Many draws: X_2 ~ N(3 + 3 / 3, 1.5^2 * 1 * 2 / 3)
#' y <- wiener_bridge(1, 3, 4, 6, 2, theta = c(mu = 0.5, sigma = 1.5), n = 10000)
#' c(mean(y), var(y))
#' c(4, 1.5^2 * 2 / 3)
#'
#' # 3. A decrease is possible: the process is not monotone
#' wiener_bridge(0, 2, 2, 1, 0.5, theta = c(mu = 0.5, sigma = 1.5), n = 3)
#' @examplesIf requireNamespace("levyfit", quietly = TRUE)
#' # 4. Parameters estimated with levyfit (variance sigma2 turned into sigma)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' x <- DegradSim::sim_wiener(t, mu = 0.5, sigma = 1.5)
#' est <- levyfit::wiener_mle(x)
#' theta <- c(mu = est[["mu"]], sigma = sqrt(est[["sigma2"]]))
#' wiener_bridge(1.2, x[2, 2], 2, x[3, 2], 1.5, theta = theta, n = 3)
#' @export
wiener_bridge <- function(t_A, x_A, t_B, x_B, t, theta, n = 1) {
  p <- check_theta(theta, c("mu", "sigma"), real = "mu")
  bridge_draws(t_A, x_A, t_B, x_B, t, n, wiener_draw(p[["sigma"]]), "Wiener")
}

#' Wiener process: resampling on a regular grid
#'
#' Resamples one or several trajectories of a Wiener process with drift onto
#' a regular time grid: each grid point that is not an inspection time is
#' drawn with [wiener_bridge()] between the known levels around it. The grids
#' and the resampling algorithm are described in [levyBridge-package].
#'
#' @inheritParams gamma_resampling
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(mu = , sigma = )`, with the volatility \eqn{\sigma}
#'   (not the variance, see [wiener_bridge()]). The bridge uses only `sigma`;
#'   `mu` is used by the option 2 of the method `"MultiBridge"`, to extend the
#'   trajectories.
#'
#' @inherit gamma_resampling return
#'
#' @seealso [wiener_bridge()] for the distribution of the bridge,
#'   [levyBridge-package] for the grids and the algorithm, [plot_bridge()].
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (mu = 0.5, sigma = 1.5)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- DegradSim::sim_wiener(t, mu = 0.5, sigma = 1.5)
#' r <- wiener_resampling(x, theta = c(mu = 0.5, sigma = 1.5))
#' r
#' r[[1]]
#' plot_bridge(r)
#'
#' # 2. Method "ExtBridge": the grid 0, 0.1, ..., 22.5 contains all the times
#' r2 <- wiener_resampling(x, theta = c(mu = 0.5, sigma = 1.5), method = "ExtBridge")
#' attr(r2, "teq")
#' plot_bridge(r2)
#'
#' # 3. Three trajectories with their own inspection times, one common step
#' plan <- DegradSim::inspection_times("decreasing", n = 10, n_traj = 3,
#'                                     tmax = c(20, 30, 40), digits = 1)
#' xs <- DegradSim::sim_wiener(plan, mu = 0.5, sigma = 1.5)
#' r3 <- wiener_resampling(xs, theta = c(mu = 0.5, sigma = 1.5), method = "MultiBridge")
#' attr(r3, "teq")
#' plot_bridge(r3, trajectory = 1:3)
#' @export
wiener_resampling <- function(x, theta, method = "Bridge", digits = NULL, common_step = FALSE, option = 1) {
  p <- check_theta(theta, c("mu", "sigma"), real = "mu")
  sim <- function(s) stats::rnorm(1, p[["mu"]] * s, p[["sigma"]] * sqrt(s))
  bridge_paths(x, method, wiener_draw(p[["sigma"]]), "Wiener", p, digits = digits, common_step = common_step,
               option = option, sim = sim)
}

# Pont brownien : accroissements à gauche et à droite de t
wiener_draw <- function(sigma) {
  function(t_inf, delta, t, t_sup) {
    tau <- t_sup - t_inf
    e <- stats::rnorm(1, 0, sigma * sqrt((t - t_inf) * (t_sup - t) / tau))
    c(delta * (t - t_inf) / tau + e, delta * (t_sup - t) / tau - e)
  }
}
