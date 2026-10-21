#' Brownian bridge: resampling on a regular grid
#'
#' Resamples one or several trajectories of a Wiener process with drift,
#' observed at irregular inspection times, onto a regular time grid, by
#' exact simulation of the Brownian bridge between the known values.
#'
#' @section Model:
#' The Wiener process with drift \eqn{\mu} and volatility \eqn{\sigma > 0} is
#' \eqn{X_t = \mu t + \sigma B_t}, with \eqn{B} a standard Brownian motion:
#' \eqn{X_0 = 0}, independent increments and
#' \eqn{X_t - X_s \sim N(\mu (t - s), \sigma^2 (t - s))}. Its trajectories are
#' not monotone.
#'
#' @section Law of the bridge:
#' Let \eqn{t_A < t < t_B}, \eqn{\tau_1 = t - t_A}, \eqn{\tau_2 = t_B - t},
#' \eqn{T = \tau_1 + \tau_2} and \eqn{\Delta = x_B - x_A}. Given
#' \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B},
#' \deqn{X_t \sim N\Big(x_A + \Delta \frac{\tau_1}{T},\ \sigma^2 \frac{\tau_1 \tau_2}{T}\Big).}
#' The law of the bridge does not depend on the drift \eqn{\mu}: only
#' \eqn{\sigma} is needed.
#'
#' @inheritSection gamma_bridge Resampling algorithm (method "Bridge")
#' @inheritSection gamma_bridge Method "ExtBridge"
#'
#' @inheritParams gamma_bridge
#' @param sigma Volatility \eqn{\sigma > 0} (standard deviation per square
#'   root of time unit). The package `levyfit` returns the variance:
#'   use `sqrt(wiener_mle(x)[["sigma2"]])`.
#'
#' @inherit gamma_bridge return
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' @seealso
#' * [increments()]: the reconstructed increments;
#' * [plot_bridge()]: observed and resampled values;
#' * [gamma_bridge()]: bridge of the gamma process;
#' * [ig_bridge()]: bridge of the inverse Gaussian process;
#' * [vg_bridge()]: bridge of the variance-gamma process;
#' * [meixner_bridge()]: bridge of the Meixner process.
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (mu = 0.5, sigma = 1.5)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- c(0, cumsum(0.5 * diff(t) + 1.5 * sqrt(diff(t)) * rnorm(10)))
#' set.seed(1)
#' r <- wiener_bridge(x, t, sigma = 1.5)
#' r
#'
#' # The resampled trajectory, on the regular grid 0, 2.25, 4.5, ..., 22.5
#' r[[1]]
#'
#' # Time step of the grid: t_N / N = 22.5 / 10
#' attr(r, "teq")
#'
#' # 2. Reconstructed increments, the input of the goodness-of-fit tests
#' increments(r)
#'
#' # 3. Observed and resampled values
#' plot_bridge(r)
#'
#' # 4. Grid points equal to inspection times keep the observed values:
#' #    here the grid is 0, 2, 4, 6, 8, 10 and 2 and 6 are inspection times
#' t2 <- c(0, 2, 3, 6, 7, 10)
#' x2 <- c(0, 1.4, 0.9, 3.1, 2.6, 4.8)
#' wiener_bridge(x2, t2, sigma = 1.5)[[1]]
#'
#' # 5. One trajectory given as a data frame (times, then levels)
#' wiener_bridge(data.frame(t = t2, X = x2), sigma = 1.5)
#'
#' # 6. Several trajectories, each with its own inspection times:
#' #    one regular grid and one time step per trajectory
#' liste <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.4, 2.3, 1.9)),
#'               data.frame(t = c(0, 2, 4, 5, 9), X = c(0, 1.6, 0.8, 3.1, 4.2)))
#' r2 <- wiener_bridge(liste, sigma = 1.5)
#' r2
#' attr(r2, "teq")
#' increments(r2)
#' plot_bridge(r2, trajectory = 1:2)
#'
#' # 7. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- wiener_bridge(x, t, sigma = 1.5)[[1]]$X
#' set.seed(2); y2 <- wiener_bridge(x, t, sigma = 1.5)[[1]]$X
#' identical(y1, y2)
#'
#' # 8. Errors: the first time must be 0, and "ExtBridge" is not available yet
#' try(wiener_bridge(x, t + 1, sigma = 1.5))
#' try(wiener_bridge(x, t, sigma = 1.5, method = "ExtBridge"))
#' @examplesIf requireNamespace("levysim", quietly = TRUE) && requireNamespace("levyfit", quietly = TRUE)
#' # 9. Ten trajectories simulated with levysim, inspected more and more often
#' plan <- levysim::inspection_times("increasing", n = 20, n_traj = 10, tmax = 40)
#' sims <- levysim::sim_wiener(plan, mu = 0.5, sigma = 1.5)
#'
#' # Parameters estimated with levyfit
#' sigma_hat <- sqrt(levyfit::wiener_mle(sims)[["sigma2"]])
#' r <- wiener_bridge(sims, sigma = sigma_hat)
#' r
#' plot_bridge(r, trajectory = 1:4)
#'
#' # Reconstructed increments of all the trajectories
#' d <- increments(r)
#' head(d)
#' # Under the model, the increments follow N(mu * teq, sigma^2 * teq)
#' ks.test(d$dx, "pnorm", mean = 0.5 * d$teq[1], sd = 1.5 * sqrt(d$teq[1]))
#' @export
wiener_bridge <- function(x, t = NULL, sigma, method = "Bridge") {
  check_par(sigma = sigma)
  draw <- function(t_inf, x_inf, t, t_sup, x_sup) {
    m <- x_inf + (x_sup - x_inf) * (t - t_inf) / (t_sup - t_inf)
    v <- sigma^2 * (t - t_inf) * (t_sup - t) / (t_sup - t_inf)
    stats::rnorm(1, m, sqrt(v))
  }
  bridge_paths(x, t, method, draw, "Wiener", c(sigma = unname(sigma)))
}
