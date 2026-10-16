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
#' \eqn{Y} is drawn exactly by acceptance-rejection. The envelope uses that
#' \eqn{h_\alpha(y)} decreases in \eqn{|y|}: on \eqn{[lo, hi]}, with
#' \eqn{lo = \min(0, z) - \ell}, \eqn{hi \ge \max(0, z) + \ell} and
#' \eqn{\ell = 3 + 2 \sqrt{\alpha_1 + \alpha_2}}, it is constant on cells of
#' width 0.1, refined geometrically around the peaks at 0 and \eqn{z}; on the
#' two tails, it is a Cauchy density, since for \eqn{|y| \ge |y_0|},
#' \eqn{h_\alpha(y) \le h_\alpha(y_0) (\alpha^2 + y_0^2) / (\alpha^2 + y^2)}.
#' The cells only change the speed, not the exactness. The function
#' \eqn{\log h_\alpha} is computed with [gsl::lngamma_complex()].
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
#' @seealso [increments()], [plot_bridge()], [gamma_bridge()],
#'   [wiener_bridge()], [ig_bridge()], [vg_bridge()]
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
#' c(mean(d$dx), (0.2 + 1.5 * tan(0.15)) * d$teq[1], var(d$dx), 1.5 * d$teq[1] / (2 * cos(0.15)^2))
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
  ell <- 3 + 2 * sqrt(alpha1 + alpha2)
  w0 <- 0.1
  lo <- min(0, z) - ell
  hi <- lo + w0 * ceiling((abs(z) + 2 * ell) / w0)
  refine <- function(centre, alpha) {
    if (alpha < w0) centre + c(0, outer(c(-1, 1), alpha * 2^(0:floor(log2(w0 / alpha))))) else centre
  }
  e <- c(seq(lo, hi, by = w0), refine(0, alpha1), refine(z, alpha2))
  e <- sort(unique(c(lo, e[e > lo & e < hi], hi)))
  K <- length(e) - 1
  w <- diff(e)

  d1 <- pmax(0, e[-(K + 1)], -e[-1])
  d2 <- pmax(0, e[-(K + 1)] - z, z - e[-1])
  lM <- lh(alpha1, d1) + lh(alpha2, d2)

  lCL <- lh(alpha1, lo) + lh(alpha2, z - lo) + log(alpha1^2 + lo^2)
  thL <- atan(alpha1 / abs(lo))
  lCR <- lh(alpha2, hi - z) + lh(alpha1, hi) + log(alpha2^2 + (hi - z)^2)
  thR <- atan(alpha2 / (hi - z))

  lmass <- c(lCL + log(thL / alpha1), lM + log(w), lCR + log(thR / alpha2))
  cum <- cumsum(exp(lmass - max(lmass)))

  repeat {
    k <- findInterval(stats::runif(1) * cum[K + 2], cum) + 1
    u <- stats::runif(1)
    if (k == 1) {
      y <- -alpha1 / tan(u * thL)
      le <- lCL - log(alpha1^2 + y^2)
    } else if (k == K + 2) {
      y <- z + alpha2 / tan(u * thR)
      le <- lCR - log(alpha2^2 + (y - z)^2)
    } else {
      y <- e[k - 1] + w[k - 1] * u
      le <- lM[k - 1]
    }
    if (log(stats::runif(1)) <= lh(alpha1, y) + lh(alpha2, z - y) - le) return(y)
  }
}
