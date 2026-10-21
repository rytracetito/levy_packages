#' Meixner bridge: resampling on a regular grid
#'
#' Resamples one or several trajectories of a Meixner process, observed at
#' irregular inspection times, onto a regular time grid, by exact simulation
#' of the Meixner bridge between the known values.
#'
#' @section Model:
#' The \eqn{Meixner(a, b, d, m)} distribution, \eqn{a > 0},
#' \eqn{-\pi < b < \pi}, \eqn{d > 0}, \eqn{m} real, has density
#' \deqn{f(x) = \frac{(2 \cos(b/2))^{2d}}{2 a \pi\, \Gamma(2d)}
#'   \exp\Big(\frac{b (x - m)}{a}\Big) \Big|\Gamma\Big(d + i \frac{x - m}{a}\Big)\Big|^2.}
#' The Meixner process has \eqn{X_0 = 0}, independent increments and
#' \eqn{X_t - X_s \sim Meixner(a, b, d (t - s), m (t - s))}. Its trajectories
#' are not monotone.
#'
#' @section Law of the bridge:
#' Let \eqn{t_A < t < t_B}, \eqn{h_\alpha(y) = |\Gamma(\alpha + i y)|^2},
#' \deqn{\alpha_1 = d (t - t_A), \qquad \alpha_2 = d (t_B - t), \qquad
#'   z = \frac{(x_B - x_A) - m (t_B - t_A)}{a}.}
#' Given \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B},
#' \eqn{X_t = x_A + m (t - t_A) + a Y}, where \eqn{Y} has the density
#' \deqn{q(y) = \frac{h_{\alpha_1}(y)\, h_{\alpha_2}(z - y)}{2 \pi B(2 \alpha_1, 2 \alpha_2)\, h_{\alpha_1 + \alpha_2}(z)}.}
#' The law of the bridge does not depend on \eqn{b}: only \eqn{a}, \eqn{d} and
#' \eqn{m} are needed.
#'
#' @section Simulation of the bridge:
#' \eqn{Y} is drawn exactly by acceptance-rejection with a Cauchy envelope
#' (algorithm 6 of Tito and Paroissin). For \eqn{\alpha > 0} and \eqn{y} real,
#' \eqn{h_\alpha(y) \le \Gamma(\alpha)^2} and
#' \eqn{h_\alpha(y) \le \Gamma(\alpha)^2 \alpha^2 / (\alpha^2 + y^2)}. With
#' \eqn{(y_0, s) = (0, \alpha_1)} if \eqn{\alpha_1 \le \alpha_2} and
#' \eqn{(y_0, s) = (z, \alpha_2)} otherwise, this gives, for all \eqn{y},
#' \deqn{h_{\alpha_1}(y)\, h_{\alpha_2}(z - y) \le \Gamma(\alpha_1)^2 \Gamma(\alpha_2)^2 \frac{s^2}{s^2 + (y - y_0)^2},}
#' a multiple of the Cauchy density of centre \eqn{y_0} and scale \eqn{s}.
#' The algorithm repeats, with \eqn{U} and \eqn{V} independent uniform
#' variables on \eqn{(0, 1)}:
#' 1. \eqn{Y = y_0 + s \tan(\pi (V - 1/2))}, a draw of the Cauchy law;
#' 2. accept \eqn{Y} if
#'    \deqn{U\, \Gamma(\alpha_1)^2 \Gamma(\alpha_2)^2 \frac{s^2}{s^2 + (Y - y_0)^2} \le h_{\alpha_1}(Y)\, h_{\alpha_2}(z - Y),}
#'    otherwise draw again.
#'
#' The number of draws follows a geometric law with parameter
#' \deqn{P_{acc} = \frac{2 B(2 \alpha_1, 2 \alpha_2)\, h_{\alpha_1 + \alpha_2}(z)}{s\, \Gamma(\alpha_1)^2 \Gamma(\alpha_2)^2}.}
#' \eqn{P_{acc}} decreases quickly when \eqn{|z|} increases (for
#' \eqn{\alpha_1 = \alpha_2 = 1}: 1/3 if \eqn{z = 0}, 0.04 if \eqn{z = 2},
#' \eqn{5 \cdot 10^{-4}} if \eqn{z = 4}): an increment far from its mean
#' \eqn{m (t_B - t_A)} can need many draws. The function \eqn{\log h_\alpha}
#' is computed with [gsl::lngamma_complex()].
#'
#' @inheritSection gamma_bridge Resampling algorithm (method "Bridge")
#' @inheritSection gamma_bridge Method "ExtBridge"
#'
#' @inheritParams gamma_bridge
#' @param a Scale parameter, \eqn{a > 0}.
#' @param d Shape parameter per unit of time, \eqn{d > 0}.
#' @param m Location parameter per unit of time, real.
#'
#' @inherit gamma_bridge return
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' Schoutens, W. and Teugels, J. L. (1998). Lévy processes, polynomials and
#' martingales. *Communications in Statistics. Stochastic Models*, 14(1-2),
#' 335-349.
#'
#' Devroye, L. (1986). *Non-Uniform Random Variate Generation*. Springer.
#'
#' Tito, A. R. and Paroissin, C. (in preparation).
#'
#' @seealso
#' * [increments()]: the reconstructed increments;
#' * [plot_bridge()]: observed and resampled values;
#' * [gamma_bridge()]: bridge of the gamma process;
#' * [wiener_bridge()]: bridge of the Wiener process;
#' * [ig_bridge()]: bridge of the inverse Gaussian process;
#' * [vg_bridge()]: bridge of the variance-gamma process.
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (a = 1, b = 0.3, d = 1.5, m = 0.2)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' x <- c(0, 0.9, 0.4, 2.8, 3.1, 5.6, 7.2, 6.5, 9.8, 11.3, 12.1)
#' set.seed(1)
#' r <- meixner_bridge(x, t, a = 1, d = 1.5, m = 0.2)
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
#' meixner_bridge(x2, t2, a = 1, d = 1.5, m = 0.2)[[1]]
#'
#' # 5. One trajectory given as a data frame (times, then levels)
#' meixner_bridge(data.frame(t = t2, X = x2), a = 1, d = 1.5, m = 0.2)
#'
#' # 6. Several trajectories, each with its own inspection times:
#' #    one regular grid and one time step per trajectory
#' liste <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.4, 2.3, 1.9)),
#'               data.frame(t = c(0, 2, 4, 5, 9), X = c(0, 1.6, 0.8, 3.1, 4.2)))
#' r2 <- meixner_bridge(liste, a = 1, d = 1.5, m = 0.2)
#' r2
#' attr(r2, "teq")
#' increments(r2)
#' plot_bridge(r2, trajectory = 1:2)
#'
#' # 7. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- meixner_bridge(x, t, a = 1, d = 1.5, m = 0.2)[[1]]$X
#' set.seed(2); y2 <- meixner_bridge(x, t, a = 1, d = 1.5, m = 0.2)[[1]]$X
#' identical(y1, y2)
#'
#' # 8. Errors: the first time must be 0, and "ExtBridge" is not available yet
#' try(meixner_bridge(x, t + 1, a = 1, d = 1.5, m = 0.2))
#' try(meixner_bridge(x, t, a = 1, d = 1.5, m = 0.2, method = "ExtBridge"))
#' @examplesIf requireNamespace("levysim", quietly = TRUE)
#' # 9. Ten trajectories simulated with levysim, inspected more and more often
#' plan <- levysim::inspection_times("increasing", n = 20, n_traj = 10, tmax = 40)
#' sims <- levysim::sim_meixner(plan, a = 1, b = 0.3, d = 1.5, m = 0.2)
#'
#' # With the true parameters (b is not needed by the bridge)
#' r <- meixner_bridge(sims, a = 1, d = 1.5, m = 0.2)
#' r
#' plot_bridge(r, trajectory = 1:4)
#'
#' # Reconstructed increments of all the trajectories
#' d <- increments(r)
#' head(d)
#' # Under the model, the increments have mean (m + a d tan(b/2)) teq
#' # and variance a^2 d teq / (2 cos(b/2)^2)
#' c(mean(d$dx), (0.2 + 1.5 * tan(0.15)) * d$teq[1])
#' c(var(d$dx), 1.5 * d$teq[1] / (2 * cos(0.15)^2))
#' @export
meixner_bridge <- function(x, t = NULL, a, d, m, method = "Bridge") {
  check_par(a = a, d = d)
  check_par(m = m, positive = FALSE)
  draw <- function(t_inf, x_inf, t, t_sup, x_sup) {
    z <- ((x_sup - x_inf) - m * (t_sup - t_inf)) / a
    x_inf + m * (t - t_inf) + a * rmeixner_bridge_y(d * (t - t_inf), d * (t_sup - t), z)
  }
  bridge_paths(x, t, method, draw, "Meixner", c(a = unname(a), d = unname(d), m = unname(m)))
}

lh <- function(alpha, y) {
  n <- max(length(alpha), length(y))
  2 * Re(gsl::lngamma_complex(rep_len(alpha, n), rep_len(y, n)))
}

rmeixner_bridge_y <- function(alpha1, alpha2, z) {
  if (alpha1 <= alpha2) {
    y0 <- 0
    s <- alpha1
  } else {
    y0 <- z
    s <- alpha2
  }
  lC <- 2 * lgamma(alpha1) + 2 * lgamma(alpha2)
  repeat {
    y <- y0 + s * tan(pi * (stats::runif(1) - 0.5))
    u <- stats::runif(1)
    if (log(u) + lC + 2 * log(s) - log(s^2 + (y - y0)^2) <= lh(alpha1, y) + lh(alpha2, z - y)) return(y)
  }
}
