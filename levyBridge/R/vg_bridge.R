#' Variance-gamma bridge: resampling on a regular grid
#'
#' Resamples one or several trajectories of a variance-gamma process,
#' observed at irregular inspection times, onto a regular time grid, by
#' exact simulation of the variance-gamma bridge between the known values.
#'
#' @section Model:
#' The variance-gamma process with parameters \eqn{\theta}, \eqn{\sigma > 0}
#' and \eqn{\nu > 0} is \eqn{X_s = \theta \gamma_s + \sigma W_{\gamma_s}}, where
#' \eqn{\gamma} is a gamma process with \eqn{\gamma_s \sim \Gamma(s / \nu)}
#' with scale \eqn{\nu} (mean \eqn{s}, variance \eqn{\nu s}) and \eqn{W} an
#' independent standard Brownian motion. Its trajectories are not monotone.
#'
#' @section Law of the bridge:
#' Let \eqn{t_A < t < t_B}, \eqn{\tau_1 = t - t_A}, \eqn{\tau_2 = t_B - t},
#' \eqn{T = \tau_1 + \tau_2}, \eqn{\Delta = x_B - x_A},
#' \deqn{p = \frac{T}{\nu} - \frac{1}{2}, \qquad \kappa^2 = \frac{\theta^2 + 2 \sigma^2 / \nu}{\sigma^4}.}
#' Given \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B},
#' \deqn{X_t = x_A + Y \Delta + \sqrt{S\, Y (1 - Y)}\ \xi,}
#' with three independent variables
#' \eqn{S \sim \mathrm{GIG}(p, \kappa^2, \Delta^2)} (density proportional to
#' \eqn{s^{p - 1} \exp(-(\kappa^2 s + \Delta^2 / s) / 2)}),
#' \eqn{Y \sim \mathrm{Beta}(\tau_1 / \nu, \tau_2 / \nu)} and
#' \eqn{\xi \sim N(0, 1)}. \eqn{S / \sigma^2} is the increment of the clock
#' \eqn{\gamma} over \eqn{[t_A, t_B]}, \eqn{Y} the fraction of it reached at
#' \eqn{t}, and, given the clock, \eqn{X} is a Brownian bridge. The law of the
#' bridge depends on \eqn{(\theta, \sigma, \nu)} only through \eqn{\nu} and
#' \eqn{\kappa}. The GIG variable is drawn with [GIGrvg::rgig()]. The case
#' \eqn{\Delta = 0} with \eqn{p \le 0}, which has probability zero, gives
#' \eqn{S = 0}.
#'
#' @inheritSection gamma_bridge Resampling algorithm (method "Bridge")
#' @inheritSection gamma_bridge Method "ExtBridge"
#'
#' @inheritParams gamma_bridge
#' @param theta Drift of the subordinated Brownian motion, real.
#' @param sigma Volatility of the subordinated Brownian motion,
#'   \eqn{\sigma > 0}.
#' @param nu Variance rate of the gamma clock, \eqn{\nu > 0}.
#'
#' @inherit gamma_bridge return
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' Madan, D. B., Carr, P. P. and Chang, E. C. (1998). The variance gamma
#' process and option pricing. *European Finance Review*, 2(1), 79-105.
#' \doi{10.1023/A:1009703431535}
#'
#' Hörmann, W. and Leydold, J. (2014). Generating generalized inverse
#' Gaussian random variates. *Statistics and Computing*, 24, 547-557.
#' \doi{10.1007/s11222-013-9387-3}
#'
#' @seealso
#' * [increments()]: the reconstructed increments;
#' * [plot_bridge()]: observed and resampled values;
#' * [gamma_bridge()]: bridge of the gamma process;
#' * [wiener_bridge()]: bridge of the Wiener process;
#' * [ig_bridge()]: bridge of the inverse Gaussian process;
#' * [meixner_bridge()]: bridge of the Meixner process.
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (theta = 0.5, sigma = 1.5, nu = 2)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' x <- c(0, 0.9, 0.4, 2.8, 3.1, 5.6, 7.2, 6.5, 9.8, 11.3, 12.1)
#' set.seed(1)
#' r <- vg_bridge(x, t, theta = 0.5, sigma = 1.5, nu = 2)
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
#' vg_bridge(x2, t2, theta = 0.5, sigma = 1.5, nu = 2)[[1]]
#'
#' # 5. One trajectory given as a data frame (times, then levels)
#' vg_bridge(data.frame(t = t2, X = x2), theta = 0.5, sigma = 1.5, nu = 2)
#'
#' # 6. Several trajectories, each with its own inspection times:
#' #    one regular grid and one time step per trajectory
#' liste <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.4, 2.3, 1.9)),
#'               data.frame(t = c(0, 2, 4, 5, 9), X = c(0, 1.6, 0.8, 3.1, 4.2)))
#' r2 <- vg_bridge(liste, theta = 0.5, sigma = 1.5, nu = 2)
#' r2
#' attr(r2, "teq")
#' increments(r2)
#' plot_bridge(r2, trajectory = 1:2)
#'
#' # 7. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- vg_bridge(x, t, theta = 0.5, sigma = 1.5, nu = 2)[[1]]$X
#' set.seed(2); y2 <- vg_bridge(x, t, theta = 0.5, sigma = 1.5, nu = 2)[[1]]$X
#' identical(y1, y2)
#'
#' # 8. Method "ExtBridge": the times have 1 decimal place, the grid
#' #    0, 0.1, ..., 22.5 contains all of them
#' r3 <- vg_bridge(x, t, theta = 0.5, sigma = 1.5, nu = 2, method = "ExtBridge")
#' attr(r3, "teq")
#' all(t %in% r3[[1]]$t)
#' plot_bridge(r3)
#'
#' # Error: the first time must be 0
#' try(vg_bridge(x, t + 1, theta = 0.5, sigma = 1.5, nu = 2))
#' @examplesIf requireNamespace("levysim", quietly = TRUE) && requireNamespace("levyfit", quietly = TRUE)
#' # 9. Ten trajectories simulated with levysim, inspected more and more often
#' plan <- levysim::inspection_times("increasing", n = 20, n_traj = 10, tmax = 40)
#' sims <- levysim::sim_vg(plan, theta = 0.5, sigma = 1.5, nu = 2)
#'
#' # Parameters estimated with levyfit
#' # With the true parameters (they can be estimated with levyfit::vg_mle())
#' r <- vg_bridge(sims, theta = 0.5, sigma = 1.5, nu = 2)
#' r
#' plot_bridge(r, trajectory = 1:4)
#'
#' # Reconstructed increments of all the trajectories
#' d <- increments(r)
#' head(d)
#' # Under the model, the increments have mean theta * teq and
#' # variance (sigma^2 + nu * theta^2) * teq
#' c(mean(d$dx), 0.5 * d$teq[1], var(d$dx), (1.5^2 + 2 * 0.5^2) * d$teq[1])
#' @export
vg_bridge <- function(x, t = NULL, theta, sigma, nu, method = "Bridge", digits = NULL, common_step = FALSE) {
  check_par(theta = theta, positive = FALSE)
  check_par(sigma = sigma, nu = nu)
  psi <- (theta^2 + 2 * sigma^2 / nu) / sigma^4
  draw <- function(t_inf, delta, t, t_sup) {
    p <- (t_sup - t_inf) / nu - 0.5
    s <- if (delta == 0 && p <= 0) 0 else GIGrvg::rgig(1, lambda = p, chi = delta^2, psi = psi)
    y <- rbeta_pair((t - t_inf) / nu, (t_sup - t) / nu)
    e <- sqrt(s * y[1] * y[2]) * stats::rnorm(1)
    c(y[1] * delta + e, y[2] * delta - e)
  }
  bridge_paths(x, t, method, draw, "variance-gamma", c(theta = unname(theta), sigma = unname(sigma), nu = unname(nu)),
               digits = digits, common_step = common_step)
}
