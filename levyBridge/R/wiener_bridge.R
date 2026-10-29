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
#' @param theta Vector \eqn{\Theta = (\mu, \sigma)} of the parameters of the process,
#'   a numeric vector named `mu` and `sigma`, with the volatility \eqn{\sigma}
#'   (not the variance). The function `wiener_mle()` of the package `levyfit`
#'   returns `mu` and `sigma2`: use
#'   `c(mu = est[["mu"]], sigma = sqrt(est[["sigma2"]]))`. The bridge uses
#'   only `sigma`.
#'
#' @inherit gamma_bridge return
#'
#' @seealso
#' * [wiener_resampling()]: resampling of trajectories on a regular grid;
#' * [gamma_bridge()], [ig_bridge()], [vg_bridge()], [meixner_bridge()]: the bridges of the other processes.
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
#'
#' # Errors: t outside (t_A, t_B), theta without all the parameters
#' try(wiener_bridge(1, 3, 4, 6, 5, theta = c(mu = 0.5, sigma = 1.5)))
#' try(wiener_bridge(1, 3, 4, 6, 2, theta = c(sigma = 1.5)))
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
#' @param theta Vector \eqn{\Theta = (\mu, \sigma)} of the parameters of the process,
#'   a numeric vector named `mu` and `sigma`, with the volatility \eqn{\sigma}
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
#' theta <- c(mu = 0.5, sigma = 1.5)
#' r <- wiener_resampling(x, theta)
#' r
#'
#' # The resampled trajectory, on the regular grid 0, 2.25, 4.5, ..., 22.5
#' r[[1]]
#' attr(r, "teq")
#' plot_bridge(r)
#'
#' # 2. One trajectory given as a data frame (times, then levels): the grid is
#' #    0, 2, 4, 6, 8, 10, and the inspection times 2 and 6 keep their levels
#' x2 <- data.frame(t = c(0, 2, 3, 6, 7, 10), X = c(0, 1.4, 0.9, 3.1, 2.6, 4.8))
#' wiener_resampling(x2, theta)[[1]]
#'
#' # 3. Method "ExtBridge": the times have 1 decimal place, so the grid
#' #    0, 0.1, ..., 22.5 contains all of them
#' r3 <- wiener_resampling(x, theta, method = "ExtBridge")
#' attr(r3, "teq")
#' all(t %in% r3[[1]]$t)
#' plot_bridge(r3)
#'
#' # digits = 0 truncates the times to integers, with a warning
#' attr(wiener_resampling(x, theta, method = "ExtBridge", digits = 0), "teq")
#'
#' # 4. Five trajectories with the same inspection times (a matrix)
#' xm <- DegradSim::sim_wiener(t, mu = 0.5, sigma = 1.5, n_traj = 5)
#' head(xm)
#' r4 <- wiener_resampling(xm, theta)
#' attr(r4, "teq")
#' plot_bridge(r4, trajectory = 1:5)
#'
#' # 5. Three trajectories, each with its own inspection times (a list)
#' plan <- DegradSim::inspection_times("decreasing", n = 10, n_traj = 3,
#'                                     tmax = c(20, 30, 40), digits = 1)
#' xs <- DegradSim::sim_wiener(plan, mu = 0.5, sigma = 1.5)
#'
#' # One grid and one step per trajectory
#' attr(wiener_resampling(xs, theta), "teq")
#' attr(wiener_resampling(xs, theta, method = "ExtBridge"), "teq")
#'
#' # The same step for all the trajectories
#' attr(wiener_resampling(xs, theta, method = "ExtBridge", common_step = TRUE), "teq")
#' r5 <- wiener_resampling(xs, theta, method = "MultiBridge")
#' attr(r5, "teq")
#' sapply(r5, nrow)
#' plot_bridge(r5, trajectory = 1:3)
#'
#' # Option 2: all the trajectories up to 40; beyond its last inspection, each
#' #    one is extended by increments of the process (this uses mu)
#' r6 <- wiener_resampling(xs, theta, method = "MultiBridge", option = 2)
#' sapply(r6, nrow)
#'
#' # 6. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- wiener_resampling(x, theta)[[1]]$X
#' set.seed(2); y2 <- wiener_resampling(x, theta)[[1]]$X
#' identical(y1, y2)
#' @examplesIf requireNamespace("levyfit", quietly = TRUE)
#' # 7. With the parameters estimated by levyfit
#' est <- levyfit::wiener_mle(xs)
#' theta_hat <- c(mu = est[["mu"]], sigma = sqrt(est[["sigma2"]]))
#' theta_hat
#' r7 <- wiener_resampling(xs, theta_hat)
#' attr(r7, "parameters")
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
