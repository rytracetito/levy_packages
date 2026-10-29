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
#' * \eqn{S \sim \mathrm{GIG}(p, \kappa^2, (x_B - x_A)^2)}, with density
#'   proportional to
#'   \eqn{s^{p - 1} \exp(-(\kappa^2 s + (x_B - x_A)^2 / s) / 2)},
#'   \eqn{p = (t_B - t_A) / \nu - 1/2} and
#'   \eqn{\kappa^2 = (\theta^2 + 2 \sigma^2 / \nu) / \sigma^4};
#' * \eqn{Y \sim \mathrm{Beta}((t - t_A) / \nu, (t_B - t) / \nu)};
#' * \eqn{\xi \sim \mathcal{N}(0, 1)}.
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
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(theta = , sigma = , nu = )` (drift, volatility of the
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
#' Madan, D. B., Carr, P. P. and Chang, E. C. (1998). The variance gamma
#' process and option pricing. *European Finance Review*, 2(1), 79-105.
#' \doi{10.1023/A:1009703431535}
#'
#' @seealso [vg_resampling()] to resample trajectories on a regular grid;
#'   [gamma_bridge()], [wiener_bridge()], [ig_bridge()], [meixner_bridge()].
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
#' @examplesIf requireNamespace("levyfit", quietly = TRUE)
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
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(theta = , sigma = , nu = )`, for instance the estimate
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
#' r <- vg_resampling(x, theta = c(theta = 0.5, sigma = 1.5, nu = 2))
#' r
#' r[[1]]
#' plot_bridge(r)
#'
#' # 2. Method "ExtBridge": the grid 0, 0.1, ..., 22.5 contains all the times
#' r2 <- vg_resampling(x, theta = c(theta = 0.5, sigma = 1.5, nu = 2), method = "ExtBridge")
#' attr(r2, "teq")
#' plot_bridge(r2)
#'
#' # 3. Three trajectories with their own inspection times, one common step
#' plan <- DegradSim::inspection_times("decreasing", n = 10, n_traj = 3,
#'                                     tmax = c(20, 30, 40), digits = 1)
#' xs <- DegradSim::sim_vg(plan, theta = 0.5, sigma = 1.5, nu = 2)
#' r3 <- vg_resampling(xs, theta = c(theta = 0.5, sigma = 1.5, nu = 2), method = "MultiBridge")
#' attr(r3, "teq")
#' plot_bridge(r3, trajectory = 1:3)
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
