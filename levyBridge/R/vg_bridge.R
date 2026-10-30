#' Variance-gamma bridge
#'
#' Draws the value at time \eqn{t} of a variance-gamma process known at two
#' times \eqn{t_A < t < t_B}, from the distribution of the variance-gamma
#' bridge. This function is used by [vg_resampling()] at each point of the
#' grid.
#'
#' @section Model:
#' The variance-gamma process with parameters \eqn{\theta}, \eqn{\sigma > 0}
#' and \eqn{\nu > 0} is \eqn{X_s = \theta \gamma_s + \sigma W_{\gamma_s}}
#' (Madan, Carr and Chang, 1998), where \eqn{\gamma} is a gamma process with
#' \eqn{\gamma_s \sim \Gamma(s / \nu)} with scale \eqn{\nu} (mean \eqn{s},
#' variance \eqn{\nu s}) and \eqn{W} an independent standard Brownian
#' motion. Its trajectories are not monotone. The vector of its parameters is
#' \eqn{\Theta = (\theta, \sigma, \nu)}: the drift \eqn{\theta} is one of the
#' parameters, not the whole vector \eqn{\Theta}.
#'
#' @section Distribution of the bridge:
#' Given \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B},
#' \deqn{X_t = x_A + (x_B - x_A)\, Y + \sqrt{S\, Y (1 - Y)}\ \xi,}
#' where the random variables \eqn{S}, \eqn{Y} and \eqn{\xi} are independent:
#' * the variable \eqn{S} follows the distribution
#'   \eqn{\mathrm{GIG}(p, \kappa^2, (x_B - x_A)^2)} in the notation of Peña and
#'   Jauch (2025), with density proportional to
#'   \deqn{s^{p - 1} \exp\big(-(\kappa^2 s + (x_B - x_A)^2 / s) / 2\big), \qquad s > 0,}
#'   where \eqn{p = (t_B - t_A) / \nu - 1/2} and
#'   \eqn{\kappa^2 = (\theta^2 + 2 \sigma^2 / \nu) / \sigma^4};
#' * \eqn{Y \sim \mathrm{Beta}((t - t_A) / \nu, (t_B - t) / \nu)};
#' * \eqn{\xi \sim \mathcal{N}(0, 1)}.
#'
#' This representation is given by Hoyle (2010, section 2.7.2): equations
#' (2.83) to (2.87) for the standard variance-gamma process, whose parameter
#' \eqn{m} is \eqn{1 / \nu} here, and equation (2.77), which writes the bridge
#' of any variance-gamma process as a scaled bridge of the standard one.
#'
#' When \eqn{x_B = x_A}, the last parameter of the GIG distribution is 0: the
#' variable \eqn{S} then has the gamma distribution with shape \eqn{p} and
#' rate \eqn{\kappa^2 / 2} if \eqn{p > 0}; if \eqn{p \le 0}, this case has
#' probability zero and the package takes \eqn{S = 0}.
#'
#' The variable \eqn{S / \sigma^2} is the increment
#' \eqn{\gamma_{t_B} - \gamma_{t_A}} of the gamma process given the
#' observations, and \eqn{Y} is the fraction of this increment reached at
#' time \eqn{t}, as in the gamma bridge; given \eqn{\gamma}, the Brownian
#' motion \eqn{W} gives the Gaussian term. The bridge uses the three
#' parameters. The GIG variable is drawn with [GIGrvg::rgig()] (Hörmann and
#' Leydold, 2014).
#'
#' @inheritParams gamma_bridge
#' @param theta Vector \eqn{\Theta = (\theta, \sigma, \nu)} of the parameters of the process,
#'   a numeric vector named `theta`, `sigma` and `nu` (drift, volatility of the
#'   subordinated Brownian motion and variance rate of the gamma process),
#'   for instance the estimate returned by `vg_mle()` of the package
#'   `levyfit`.
#'
#' @inherit gamma_bridge return
#'
#' @references
#' Hörmann, W. and Leydold, J. (2014). Generating generalized inverse
#' Gaussian random variates. *Statistics and Computing*, 24, 547-557.
#' \doi{10.1007/s11222-013-9387-3}
#'
#' Hoyle, E. (2010). *Information-Based Models for Finance and Insurance*.
#' PhD thesis, Imperial College London. arXiv:1010.0829.
#'
#' Peña, V. and Jauch, M. (2025). Properties of the generalized inverse
#' Gaussian with applications to Monte Carlo simulation and distribution
#' function evaluation. arXiv:2401.00749.
#'
#' Madan, D. B., Carr, P. P. and Chang, E. C. (1998). The variance gamma
#' process and option pricing. *European Finance Review*, 2(1), 79-105.
#' \doi{10.1023/A:1009703431535}
#'
#' @seealso
#' * [vg_resampling()]: resampling of trajectories on a regular grid;
#' * [gamma_bridge()], [wiener_bridge()], [ig_bridge()], [meixner_bridge()]: the bridges of the other processes.
#'
#' @examples
#' # 1. One draw of X_2, given X_1 = 3 and X_4 = 6 (theta = 0.5, sigma = 1.5, nu = 2)
#' set.seed(1)
#' vg_bridge(t_A = 1, x_A = 3, t_B = 4, x_B = 6, t = 2,
#'           theta = c(theta = 0.5, sigma = 1.5, nu = 2))
#'
#' # 2. Many draws: distribution of X_2 given the two levels
#' y <- vg_bridge(1, 3, 4, 6, 2, theta = c(theta = 0.5, sigma = 1.5, nu = 2), n = 10000)
#' hist(y, breaks = 50, freq = FALSE)
#'
#' # Errors: t outside (t_A, t_B), theta without all the parameters
#' try(vg_bridge(1, 3, 4, 6, 5, theta = c(theta = 0.5, sigma = 1.5, nu = 2)))
#' try(vg_bridge(1, 3, 4, 6, 2, theta = c(sigma = 1.5, nu = 2)))
#' # 3. Parameters estimated with levyfit on simulated trajectories
#' set.seed(1)
#' plan <- DegradSim::inspection_times("decreasing", n = 20, n_traj = 10, tmax = 40)
#' x <- DegradSim::sim_vg(plan, theta = 0.5, sigma = 1.5, nu = 2)
#' est <- levyfit::vg_mle(x)
#' est
#' vg_bridge(0, 0, 2, 1, 1, theta = est, n = 3)
#' @export
vg_bridge <- function(t_A, x_A, t_B, x_B, t, theta, n = 1) {
  p <- check_theta(theta, c("theta", "sigma", "nu"), real = "theta")
  bridge_draws(t_A, x_A, t_B, x_B, t, n, vg_draw(p), "variance-gamma")
}

#' Variance-gamma process: resampling on a regular grid
#'
#' Resamples one or several trajectories of a variance-gamma process onto a
#' regular time grid: each grid point that is not an inspection time is drawn
#' with [vg_bridge()] between the known levels around it. The grids and the
#' resampling algorithm are described in [levyBridge-package].
#'
#' @inheritParams gamma_resampling
#' @param theta Vector \eqn{\Theta = (\theta, \sigma, \nu)} of the parameters of the process,
#'   a numeric vector named `theta`, `sigma` and `nu`, for instance the estimate
#'   returned by `vg_mle()` of the package `levyfit`. The bridge and the
#'   option 2 of the method `"MultiBridge"` use the three parameters.
#'
#' @inherit gamma_resampling return
#'
#' @seealso [vg_bridge()] for the distribution of the bridge,
#'   [levyBridge-package] for the grids and the algorithm, [plot_bridge()].
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (theta = 0.5, sigma = 1.5, nu = 2)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- DegradSim::sim_vg(t, theta = 0.5, sigma = 1.5, nu = 2)
#' theta <- c(theta = 0.5, sigma = 1.5, nu = 2)
#' r <- vg_resampling(x, theta)
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
#' vg_resampling(x2, theta)[[1]]
#'
#' # 3. Method "ExtBridge": the times have 1 decimal place, so the grid
#' #    0, 0.1, ..., 22.5 contains all of them
#' r3 <- vg_resampling(x, theta, method = "ExtBridge")
#' attr(r3, "teq")
#' all(t %in% r3[[1]]$t)
#' plot_bridge(r3)
#'
#' # digits = 0 truncates the times to integers, with a warning
#' attr(vg_resampling(x, theta, method = "ExtBridge", digits = 0), "teq")
#'
#' # 4. Five trajectories with the same inspection times (a matrix)
#' xm <- DegradSim::sim_vg(t, theta = 0.5, sigma = 1.5, nu = 2, n_traj = 5)
#' head(xm)
#' r4 <- vg_resampling(xm, theta)
#' attr(r4, "teq")
#' plot_bridge(r4, trajectory = 1:5)
#'
#' # 5. Three trajectories, each with its own inspection times (a list)
#' plan <- DegradSim::inspection_times("decreasing", n = 10, n_traj = 3,
#'                                     tmax = c(20, 30, 40), digits = 1)
#' xs <- DegradSim::sim_vg(plan, theta = 0.5, sigma = 1.5, nu = 2)
#'
#' # One grid and one step per trajectory
#' attr(vg_resampling(xs, theta), "teq")
#' attr(vg_resampling(xs, theta, method = "ExtBridge"), "teq")
#'
#' # The same step for all the trajectories
#' attr(vg_resampling(xs, theta, method = "ExtBridge", common_step = TRUE), "teq")
#' r5 <- vg_resampling(xs, theta, method = "MultiBridge")
#' attr(r5, "teq")
#' sapply(r5, nrow)
#' plot_bridge(r5, trajectory = 1:3)
#'
#' # Option 2: all the trajectories up to 40; beyond its last inspection, each
#' #    one is extended by increments of the process
#' r6 <- vg_resampling(xs, theta, method = "MultiBridge", option = 2)
#' sapply(r6, nrow)
#'
#' # 6. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- vg_resampling(x, theta)[[1]]$X
#' set.seed(2); y2 <- vg_resampling(x, theta)[[1]]$X
#' identical(y1, y2)
#' # 7. With the parameters estimated by levyfit
#' set.seed(1)
#' plan_e <- DegradSim::inspection_times("decreasing", n = 20, n_traj = 10, tmax = 40)
#' xe <- DegradSim::sim_vg(plan_e, theta = 0.5, sigma = 1.5, nu = 2)
#' theta_hat <- levyfit::vg_mle(xe)
#' theta_hat
#' r7 <- vg_resampling(xe, theta_hat)
#' attr(r7, "parameters")
#' @export
vg_resampling <- function(x, theta, method = "Bridge", digits = NULL, common_step = FALSE, option = 1) {
  p <- check_theta(theta, c("theta", "sigma", "nu"), real = "theta")
  sim <- function(s) {
    g <- p[["nu"]] * exp(log_rgamma(s / p[["nu"]]))
    p[["theta"]] * g + p[["sigma"]] * sqrt(g) * stats::rnorm(1)
  }
  bridge_paths(x, method, vg_draw(p), "variance-gamma", p, digits = digits, common_step = common_step,
               option = option, sim = sim)
}

# Pont VG : GIG pour l'accroissement du temps gamma, bêta pour sa fraction en t, puis pont brownien
vg_draw <- function(p) {
  nu <- p[["nu"]]
  psi <- (p[["theta"]]^2 + 2 * p[["sigma"]]^2 / nu) / p[["sigma"]]^4
  function(t_inf, delta, t, t_sup) {
    lambda <- (t_sup - t_inf) / nu - 0.5
    s <- if (delta == 0 && lambda <= 0) 0 else GIGrvg::rgig(1, lambda = lambda, chi = delta^2, psi = psi)
    y <- rbeta_pair((t - t_inf) / nu, (t_sup - t) / nu)
    e <- sqrt(s * y[1] * y[2]) * stats::rnorm(1)
    c(y[1] * delta + e, y[2] * delta - e)
  }
}
