#' Meixner bridge
#'
#' Draws the value at time \eqn{t} of a Meixner process known at two times
#' \eqn{t_A < t < t_B}, from the distribution of the Meixner bridge. This
#' function is used by [meixner_resampling()] at each point of the grid.
#'
#' @section Model:
#' The \eqn{Meixner(a, b, d, m)} distribution, \eqn{a > 0},
#' \eqn{b \in ]-\pi, \pi[}, \eqn{d > 0}, \eqn{m} real, has density
#' \deqn{f(x) = \frac{(2 \cos(b/2))^{2d}}{2 a \pi\, \Gamma(2d)}
#'   \exp\Big(\frac{b (x - m)}{a}\Big) \Big|\Gamma\Big(d + i \frac{x - m}{a}\Big)\Big|^2}
#' (Schoutens and Teugels, 1998). The Meixner process has \eqn{X_0 = 0},
#' independent increments and
#' \eqn{X_t - X_s \sim Meixner(a, b, d (t - s), m (t - s))}. Its
#' trajectories are not monotone. The vector of its parameters is
#' \eqn{\Theta = (a, b, d, m)}.
#'
#' @section Distribution of the bridge:
#' Write \eqn{h_\alpha(y) = |\Gamma(\alpha + i y)|^2} and
#' \deqn{\alpha_1 = d (t - t_A), \qquad \alpha_2 = d (t_B - t), \qquad
#'   z = \frac{(x_B - x_A) - m (t_B - t_A)}{a}.}
#' Given \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B},
#' \eqn{X_t = x_A + m (t - t_A) + a Y}, where the random variable \eqn{Y} has
#' the density
#' \deqn{q(y) = \frac{h_{\alpha_1}(y)\, h_{\alpha_2}(z - y)}{2 \pi B(2 \alpha_1, 2 \alpha_2)\, h_{\alpha_1 + \alpha_2}(z)},}
#' which is the formula of the section Bridge of a Levy process of
#' [levyBridge-package] with Meixner densities. This distribution does not
#' depend on \eqn{b}.
#'
#' @section Simulation of the bridge:
#' The random variable \eqn{Y} is drawn exactly by acceptance-rejection
#' (Devroye, 1986, section II.3). For \eqn{\alpha > 0} and \eqn{y} real, the
#' product formula
#' \eqn{h_\alpha(y) = \Gamma(\alpha)^2 \prod_{n \ge 0} (1 + y^2 / (\alpha + n)^2)^{-1}}
#' shows that \eqn{h_\alpha(y)} decreases with \eqn{|y|} and that
#' \eqn{h_\alpha(y) \le \Gamma(\alpha)^2 \alpha^2 / (\alpha^2 + y^2)}. With
#' \eqn{(y_0, s) = (0, \alpha_1)} if \eqn{\alpha_1 \le \alpha_2} and
#' \eqn{(y_0, s) = (z, \alpha_2)} otherwise, this gives, for all \eqn{y},
#' \deqn{h_{\alpha_1}(y)\, h_{\alpha_2}(z - y) \le \Gamma(\alpha_1)^2 \Gamma(\alpha_2)^2 \frac{s^2}{s^2 + (y - y_0)^2},}
#' a bound proportional to the Cauchy density with location \eqn{y_0} and
#' scale \eqn{s}. The probability of acceptance with this bound is
#' \deqn{P_{acc} = \frac{2 B(2 \alpha_1, 2 \alpha_2)\, h_{\alpha_1 + \alpha_2}(z)}{s\, \Gamma(\alpha_1)^2 \Gamma(\alpha_2)^2}.}
#' When this probability is small, a tighter bound is used, which depends only
#' on \eqn{\alpha_1}, \eqn{\alpha_2} and \eqn{z}: the draw stays exact. The
#' function \eqn{\log h_\alpha} is computed with [gsl::lngamma_complex()].
#'
#' @inheritParams gamma_bridge
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(a = , b = , d = , m = )` (scale, skewness, shape and
#'   location per unit of time), for instance the estimate returned by
#'   `meixner_mle()` of the package `levyfit`. The bridge uses only `a`, `d`
#'   and `m`.
#'
#' @inherit gamma_bridge return
#'
#' @references
#' Devroye, L. (1986). *Non-Uniform Random Variate Generation*. Springer.
#'
#' Schoutens, W. and Teugels, J. L. (1998). Lévy processes, polynomials and
#' martingales. *Communications in Statistics. Stochastic Models*, 14(1-2),
#' 335-349.
#'
#' @seealso [meixner_resampling()] to resample trajectories on a regular grid;
#'   [gamma_bridge()], [wiener_bridge()], [ig_bridge()], [vg_bridge()].
#'
#' @examples
#' # 1. One draw of X_2, given X_1 = 3 and X_4 = 6 (a = 1, b = 0.3, d = 1.5, m = 0.2)
#' set.seed(1)
#' meixner_bridge(t_A = 1, x_A = 3, t_B = 4, x_B = 6, t = 2,
#'                theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2))
#'
#' # 2. Many draws: distribution of X_2 given the two levels
#' y <- meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2), n = 5000)
#' hist(y, breaks = 50, freq = FALSE)
#'
#' # The skewness b does not change the bridge
#' set.seed(2); y1 <- meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2), n = 5)
#' set.seed(2); y2 <- meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, b = -2, d = 1.5, m = 0.2), n = 5)
#' identical(y1, y2)
#'
#' # Errors: t outside (t_A, t_B), theta without all the parameters
#' try(meixner_bridge(1, 3, 4, 6, 5, theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2)))
#' try(meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, d = 1.5, m = 0.2)))
#' try(meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, b = 4, d = 1.5, m = 0.2)))
#' @examplesIf requireNamespace("levyfit", quietly = TRUE)
#' # 3. Parameters estimated with levyfit on simulated trajectories
#' set.seed(1)
#' plan <- DegradSim::inspection_times("decreasing", n = 20, n_traj = 10, tmax = 40)
#' x <- DegradSim::sim_meixner(plan, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' est <- levyfit::meixner_mle(x)
#' est
#' meixner_bridge(0, 0, 2, 1, 1, theta = est, n = 3)
#' @export
meixner_bridge <- function(t_A, x_A, t_B, x_B, t, theta, n = 1) {
  p <- check_meixner(theta)
  bridge_draws(t_A, x_A, t_B, x_B, t, n, meixner_draw(p), "Meixner")
}

#' Meixner process: resampling on a regular grid
#'
#' Resamples one or several trajectories of a Meixner process onto a regular
#' time grid: each grid point that is not an inspection time is drawn with
#' [meixner_bridge()] between the known levels around it. The grids and the
#' resampling algorithm are described in [levyBridge-package].
#'
#' @inheritParams gamma_resampling
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(a = , b = , d = , m = )`, for instance the estimate
#'   returned by `meixner_mle()` of the package `levyfit`. The bridge uses
#'   only `a`, `d` and `m`; `b` is used by the option 2 of the method
#'   `"MultiBridge"`, to extend the trajectories.
#'
#' @inherit gamma_resampling return
#'
#' @seealso [meixner_bridge()] for the distribution of the bridge,
#'   [levyBridge-package] for the grids and the algorithm, [plot_bridge()].
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (a = 1, b = 0.3, d = 1.5, m = 0.2)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- DegradSim::sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' theta <- c(a = 1, b = 0.3, d = 1.5, m = 0.2)
#' r <- meixner_resampling(x, theta)
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
#' meixner_resampling(x2, theta)[[1]]
#'
#' # 3. Method "ExtBridge": the times have 1 decimal place, so the grid
#' #    0, 0.1, ..., 22.5 contains all of them
#' r3 <- meixner_resampling(x, theta, method = "ExtBridge")
#' attr(r3, "teq")
#' all(t %in% r3[[1]]$t)
#' plot_bridge(r3)
#'
#' # digits = 0 truncates the times to integers, with a warning
#' attr(meixner_resampling(x, theta, method = "ExtBridge", digits = 0), "teq")
#'
#' # 4. Five trajectories with the same inspection times (a matrix)
#' xm <- DegradSim::sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2, n_traj = 5)
#' head(xm)
#' r4 <- meixner_resampling(xm, theta)
#' attr(r4, "teq")
#' plot_bridge(r4, trajectory = 1:5)
#'
#' # 5. Three trajectories, each with its own inspection times (a list)
#' plan <- DegradSim::inspection_times("decreasing", n = 10, n_traj = 3,
#'                                     tmax = c(20, 30, 40), digits = 1)
#' xs <- DegradSim::sim_meixner(plan, a = 1, b = 0.3, d = 1.5, m = 0.2)
#'
#' # One grid and one step per trajectory
#' attr(meixner_resampling(xs, theta), "teq")
#' attr(meixner_resampling(xs, theta, method = "ExtBridge"), "teq")
#'
#' # The same step for all the trajectories
#' attr(meixner_resampling(xs, theta, method = "ExtBridge", common_step = TRUE), "teq")
#' r5 <- meixner_resampling(xs, theta, method = "MultiBridge")
#' attr(r5, "teq")
#' sapply(r5, nrow)
#' plot_bridge(r5, trajectory = 1:3)
#'
#' # Option 2: all the trajectories up to 40; beyond its last inspection, each
#' #    one is extended by increments of the process (this uses b)
#' r6 <- meixner_resampling(xs, theta, method = "MultiBridge", option = 2)
#' sapply(r6, nrow)
#'
#' # 6. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- meixner_resampling(x, theta)[[1]]$X
#' set.seed(2); y2 <- meixner_resampling(x, theta)[[1]]$X
#' identical(y1, y2)
#' @examplesIf requireNamespace("levyfit", quietly = TRUE)
#' # 7. With the parameters estimated by levyfit
#' set.seed(1)
#' plan_e <- DegradSim::inspection_times("decreasing", n = 20, n_traj = 10, tmax = 40)
#' xe <- DegradSim::sim_meixner(plan_e, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' theta_hat <- levyfit::meixner_mle(xe)
#' theta_hat
#' r7 <- meixner_resampling(xe, theta_hat)
#' attr(r7, "parameters")
#' @export
meixner_resampling <- function(x, theta, method = "Bridge", digits = NULL, common_step = FALSE, option = 1) {
  p <- check_meixner(theta)
  gen <- new.env()
  sim <- function(s) {
    key <- format(s, digits = 15)
    if (is.null(gen[[key]])) {
      gen[[key]] <- Runuran::pinvd.new(Runuran::udmeixner(alpha = p[["a"]], beta = p[["b"]], delta = p[["d"]] * s,
                                                          mu = p[["m"]] * s))
    }
    Runuran::ur(gen[[key]], 1)
  }
  bridge_paths(x, method, meixner_draw(p), "Meixner", p, digits = digits, common_step = common_step,
               option = option, sim = sim)
}

check_meixner <- function(theta) {
  p <- check_theta(theta, c("a", "b", "d", "m"), real = c("b", "m"))
  if (abs(p[["b"]]) >= pi) stop("b must be in (-pi, pi)", call. = FALSE)
  p
}

# Pont de Meixner : X_t - x_A = m (t - t_A) + a Y, Y tiré exactement par rejet
meixner_draw <- function(p) {
  function(t_inf, delta, t, t_sup) {
    z <- (delta - p[["m"]] * (t_sup - t_inf)) / p[["a"]]
    y <- rmeixner_bridge_y(p[["d"]] * (t - t_inf), p[["d"]] * (t_sup - t), z)
    c(p[["m"]] * (t - t_inf) + p[["a"]] * y, p[["m"]] * (t_sup - t) + p[["a"]] * (z - y))
  }
}

lh <- function(alpha, y) {
  n <- max(length(alpha), length(y))
  2 * Re(gsl::lngamma_complex(rep_len(alpha, n), rep_len(y, n)))
}

# Taux d'acceptation exact de l'enveloppe de Cauchy : la constante de normalisation de la loi cible est
# 2 pi B(2 alpha1, 2 alpha2) h_{alpha1 + alpha2}(z). Le choix de l'algorithme ne dépend que des
# paramètres : le tirage reste exact.
rmeixner_bridge_y <- function(alpha1, alpha2, z) {
  s <- min(alpha1, alpha2)
  log_rate <- log(2) + lbeta(2 * alpha1, 2 * alpha2) + lh(alpha1 + alpha2, z) -
    log(s) - 2 * lgamma(alpha1) - 2 * lgamma(alpha2)
  if (log_rate >= log(0.05)) rmeixner_bridge_cauchy(alpha1, alpha2, z) else rmeixner_bridge_cells(alpha1, alpha2, z)
}

rmeixner_bridge_cauchy <- function(alpha1, alpha2, z) {
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

# Enveloppe constante par cases de largeur 0.1 (raffinées près des pics en 0 et z) et queues de Cauchy :
# h_alpha(y) = |Gamma(alpha + i y)|^2 décroît en |y|, le taux d'acceptation reste bon pour tous les paramètres.
rmeixner_bridge_cells <- function(alpha1, alpha2, z) {
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
