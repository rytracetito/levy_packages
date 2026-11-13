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
#' (Devroye, 1986, section II.3): a value is proposed from a simpler density,
#' the envelope, and accepted or rejected so that the accepted values follow
#' exactly the distribution of \eqn{Y}. The envelope is built in three steps.
#'
#' 1. **A bound of the gamma function.** For \eqn{\alpha > 0},
#'    \deqn{|\Gamma(\alpha + i y)|^2 \le \Gamma(\alpha)^2 \frac{\alpha^2}{\alpha^2 + y^2}
#'    \Big(1 + \frac{y^2}{(\alpha + 1)^2}\Big)^{\alpha + 1/2}
#'    \exp\Big(-2 |y| \arctan \frac{|y|}{\alpha + 1}\Big),}
#'    a bound close to \eqn{|\Gamma(\alpha + i y)|^2} for all \eqn{y}, which
#'    uses only elementary functions.
#' 2. **A piecewise exponential envelope.** The product of the bounds for
#'    \eqn{h_{\alpha_1}(y)} and \eqn{h_{\alpha_2}(z - y)} is bounded above by an
#'    exponential function on each of a set of intervals, which are halved
#'    until this exponential function is close to the product. Beyond these
#'    intervals, the two tails are bounded by exponential functions.
#' 3. **Acceptance-rejection.** A piece of the envelope is chosen with a
#'    probability proportional to its area, a value is drawn from it by
#'    inversion of its distribution function, and this value is accepted with
#'    the probability given by the ratio of the density of \eqn{Y} to the
#'    envelope.
#'
#' **Guarantee.** For all the parameters, a proposed value is accepted with
#' probability at least 0.57, so that at most 1.75 values are proposed on
#' average for each draw. The cost of the envelope changes little with
#' \eqn{z}, and for \eqn{n} draws with the same endpoints it is built only
#' once.
#'
#' **Details.** The bound of step 1 follows from the product formula of the
#' gamma function (DLMF, equation 5.8.3) and the partial fraction expansion of
#' \eqn{\coth} (DLMF, equation 4.36.3). The function
#' \eqn{\log |\Gamma(\alpha + i y)|^2} is computed with
#' [gsl::lngamma_complex()], and only when the bound of step 1 and a similar
#' lower bound cannot decide the acceptance alone. When \eqn{\alpha_1} and
#' \eqn{\alpha_2} are both smaller than about \eqn{10^{-15} |z|}, the two peaks
#' of the density cannot be resolved with floating point numbers: the draws
#' stay exact, but the probability of acceptance can be lower.
#'
#' @inheritParams gamma_bridge
#' @param theta Vector \eqn{\Theta = (a, b, d, m)} of the parameters of the process,
#'   a numeric vector named `a`, `b`, `d` and `m` (scale, skewness, shape and
#'   location per unit of time), for instance the estimate returned by
#'   `meixner_mle()` of the package `levyfit`. The bridge uses only `a`, `d`
#'   and `m`.
#'
#' @inherit gamma_bridge return
#'
#' @references
#' Devroye, L. (1986). *Non-Uniform Random Variate Generation*. Springer.
#'
#' NIST Digital Library of Mathematical Functions (DLMF), release 1.2.8,
#' 2026-09-15, <https://dlmf.nist.gov/>: equation 5.8.3 (product formula
#' of the gamma function) and equation 4.36.3 (partial fraction expansion of
#' \eqn{\coth}).
#'
#' Schoutens, W. and Teugels, J. L. (1998). Lévy processes, polynomials and
#' martingales. *Communications in Statistics. Stochastic Models*, 14(1-2),
#' 335-349.
#'
#' @seealso
#' * [meixner_resampling()]: resampling of trajectories on a regular grid;
#' * [gamma_bridge()], [wiener_bridge()], [ig_bridge()], [vg_bridge()]: the bridges of the other processes.
#'
#' @examples
#' # 1. One draw of X_2, given X_1 = 3 and X_4 = 6 (a = 1, b = 0.3, d = 1.5, m = 0.2)
#' set.seed(1)
#' meixner_bridge(t_A = 1, x_A = 3, t_B = 4, x_B = 6, t = 2,
#'                theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2))
#'
#' # 2. Many draws: the histogram of the draws and, in red, the density of
#' #    the bridge, f_{t - t_A}(x - x_A) f_{t_B - t}(x_B - x) / f_{t_B - t_A}(x_B - x_A)
#' #    (see ?levyBridge), where f_s is the density of X_s
#' #    (see the section Model, computed on the log scale)
#' f <- function(x, s, a = 1, b = 0.3, d = 1.5, m = 0.2) {
#'   exp(2 * d * s * log(2 * cos(b / 2)) - log(2 * a * pi) - lgamma(2 * d * s) +
#'         b * (x - m * s) / a +
#'         2 * Re(gsl::lngamma_complex(rep(d * s, length(x)), (x - m * s) / a)))
#' }
#' dens <- function(x) f(x - 3, 1) * f(6 - x, 2) / f(3, 3)
#' y <- meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2), n = 50000)
#' hist(y, breaks = 50, freq = FALSE)
#' curve(dens, add = TRUE, col = "red")
#'
#' # The skewness b does not change the bridge
#' set.seed(2)
#' y1 <- meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2), n = 5)
#' set.seed(2)
#' y2 <- meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, b = -2, d = 1.5, m = 0.2), n = 5)
#' identical(y1, y2)
#'
#' # A bridge far from its typical value: x_B - x_A is more than 1000 times its
#' # standard deviation, and the draws stay fast and exact
#' system.time(y <- meixner_bridge(0, 0, 1, 2000, 0.5,
#'                                 theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2), n = 1000))
#' summary(y)
#'
#' # Errors: t outside (t_A, t_B), theta without all the parameters
#' try(meixner_bridge(1, 3, 4, 6, 5, theta = c(a = 1, b = 0.3, d = 1.5, m = 0.2)))
#' try(meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, d = 1.5, m = 0.2)))
#' try(meixner_bridge(1, 3, 4, 6, 2, theta = c(a = 1, b = 4, d = 1.5, m = 0.2)))
#'
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
  bridge_draws(t_A, x_A, t_B, x_B, t, n, meixner_draw(p), "Meixner", draw_n = meixner_draw(p))
}

#' Meixner process: resampling on a regular grid
#'
#' Resamples one or several trajectories of a Meixner process onto a regular
#' time grid: each grid point that is not an inspection time is drawn with
#' [meixner_bridge()] between the known levels around it. The grids and the
#' resampling algorithm are described in [levyBridge-package].
#'
#' @inheritParams gamma_resampling
#' @param theta Vector \eqn{\Theta = (a, b, d, m)} of the parameters of the process,
#'   a numeric vector named `a`, `b`, `d` and `m`, for instance the estimate
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
#' plot_bridge(r6, trajectory = 1:3)
#'
#' # 6. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- meixner_resampling(x, theta)[[1]]$X
#' set.seed(2); y2 <- meixner_resampling(x, theta)[[1]]$X
#' identical(y1, y2)
#'
#' # 7. With the parameters estimated by levyfit
#' set.seed(1)
#' plan_e <- DegradSim::inspection_times("decreasing", n = 20, n_traj = 10, tmax = 40)
#' xe <- DegradSim::sim_meixner(plan_e, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' theta_hat <- levyfit::meixner_mle(xe)
#' theta_hat
#' r7 <- meixner_resampling(xe, theta_hat)
#' attr(r7, "parameters")
#'
#' # 8. Four trajectories observed up to 8, 12, 45 and 58, each with its own
#' #    inspection times, resampled with "MultiBridge"
#' liste_t <- list(c(0, 3, 8), c(0, 2, 5, 9, 12),
#'                 c(0, 4, 7, 11, 15, 18, 22, 27, 33, 39, 45),
#'                 c(0, 2, 3, 5, 8, 9, 14, 15, 17, 20, 22, 25, 27, 30, 33, 46, 58))
#'
#' # Meixner process (a = 1, b = 0.3, d = 1.5, m = 0.2). With these 32
#' # increments, levyfit::meixner_mle() often finds no maximum inside the
#' # parameter space: the true parameters are used
#' set.seed(46)
#' x8 <- DegradSim::sim_meixner(liste_t, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' theta8 <- c(a = 1, b = 0.3, d = 1.5, m = 0.2)
#'
#' # Option 1: one grid for all the trajectories; each one stops at the last
#' # grid point before its last inspection
#' r8 <- meixner_resampling(x8, theta8, method = "MultiBridge")
#' attr(r8, "teq")
#' sapply(r8, function(d) max(d$t))
#' plot_bridge(r8, trajectory = 1:4)
#'
#' # Option 2: the grid ends exactly at 58; beyond its last inspection, each
#' # shorter trajectory is extended by increments of the process
#' r8b <- meixner_resampling(x8, theta8, method = "MultiBridge", option = 2)
#' attr(r8b, "teq")
#' sapply(r8b, function(d) max(d$t))
#' plot_bridge(r8b, trajectory = 1:4)
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
  function(t_inf, delta, t, t_sup, n = 1) {
    z <- (delta - p[["m"]] * (t_sup - t_inf)) / p[["a"]]
    y <- rmeixner_bridge_y(p[["d"]] * (t - t_inf), p[["d"]] * (t_sup - t), z, n)
    d <- cbind(p[["m"]] * (t - t_inf) + p[["a"]] * y, p[["m"]] * (t_sup - t) + p[["a"]] * (z - y))
    if (n == 1) d[1, ] else d
  }
}

lh <- function(alpha, y) {
  n <- max(length(alpha), length(y))
  2 * Re(gsl::lngamma_complex(rep_len(alpha, n), rep_len(y, n)))
}

# Y est tiré par acceptation-rejet avec l'enveloppe exponentielle par morceaux, construite une seule fois pour
# les n tirages. Si Y a la densité associée à (alpha1, alpha2, z), z - Y a celle associée à (alpha2, alpha1, z) :
# le pic le plus étroit est placé en 0, où les nombres flottants sont les plus fins.
rmeixner_bridge_y <- function(alpha1, alpha2, z, n = 1) {
  za <- abs(z)
  y <- if (alpha1 <= alpha2) rmeixner_bridge_envelope(alpha1, alpha2, za, n) else za - rmeixner_bridge_envelope(alpha2, alpha1, za, n)
  sign(z + (z == 0)) * y
}

# Majorant b_alpha de |Gamma(alpha + i y)|^2 / Gamma(alpha)^2, à un facteur au plus exp(1 / (6 (alpha + 1))) près
meixner_log_b <- function(alpha, y) {
  t <- abs(y)
  b <- alpha + 1
  2 * lgamma(alpha) + log(alpha^2 / (alpha^2 + t^2)) + (alpha + 0.5) * log1p(t^2 / b^2) - 2 * t * atan(t / b)
}

meixner_log_E <- function(alpha1, alpha2, z, y) meixner_log_b(alpha1, y) + meixner_log_b(alpha2, z - y)

# Test d'acceptation avec compression : log_E - correction <= log q <= log_E ; log q n'est calculé que si nécessaire
meixner_accept <- function(lu, alpha1, alpha2, z, y) {
  w <- z - y
  lE <- meixner_log_E(alpha1, alpha2, z, y)
  lmin <- lE - y^2 / (6 * (alpha1 + 1) * ((alpha1 + 1)^2 + y^2)) - w^2 / (6 * (alpha2 + 1) * ((alpha2 + 1)^2 + w^2))
  ok <- lu <= lmin
  j <- !ok & lu <= lE
  if (any(j)) ok[j] <- lu[j] <= lh(alpha1, y[j]) + lh(alpha2, w[j])
  ok
}

# Encadrement de t / (c^2 + t^2) pour t dans [tl, tu]
meixner_ratio_range <- function(cc, tl, tu) {
  fl <- tl / (cc^2 + tl^2)
  fu <- tu / (cc^2 + tu^2)
  list(lo = pmin(fl, fu), hi = ifelse(tl < cc & cc < tu, 1 / (2 * cc), pmax(fl, fu)))
}

# Encadrement de la dérivée de log b_alpha(t) pour t dans [tl, tu]
meixner_dlogb_range <- function(alpha, tl, tu) {
  b <- alpha + 1
  ra <- meixner_ratio_range(alpha, tl, tu)
  rb <- meixner_ratio_range(b, tl, tu)
  list(lo = -2 * ra$hi - rb$hi - 2 * atan(tu / b), hi = -2 * ra$lo - rb$lo - 2 * atan(tl / b))
}

# Encadrement de la dérivée de log E sur les cellules [lo, hi], qui ne contiennent ni 0 ni z en leur intérieur
meixner_dlogE_range <- function(alpha1, alpha2, z, lo, hi) {
  pos <- lo >= 0
  g1 <- meixner_dlogb_range(alpha1, ifelse(pos, lo, -hi), ifelse(pos, hi, -lo))
  wl <- z - hi
  wh <- z - lo
  pw <- wl >= 0
  g2 <- meixner_dlogb_range(alpha2, ifelse(pw, wl, -wh), ifelse(pw, wh, -wl))
  list(lo = ifelse(pos, g1$lo, -g1$hi) + ifelse(pw, -g2$hi, g2$lo),
       hi = ifelse(pos, g1$hi, -g1$lo) + ifelse(pw, -g2$lo, g2$hi))
}

# Enveloppe exponentielle par morceaux de la loi réduite, pour z >= 0 : cellules sur [-L, z + R] et deux queues
# exponentielles. Une cellule est coupée en deux tant que l'enveloppe y dépasse E d'un facteur supérieur à rho et
# que sa hauteur dépasse eps K / (z + L + R) ; la probabilité d'acceptation est au moins 1 / (rho e^(1/3) + 3 eps).
meixner_envelope <- function(alpha1, alpha2, z, rho = 1.25, eps = 1e-4) {
  ref <- log(2 * pi) + lbeta(2 * alpha1, 2 * alpha2) + lh(alpha1 + alpha2, z) + log(eps)
  F <- function(y) meixner_log_E(alpha1, alpha2, z, y)
  kappa <- function(alpha, t) 2 * atan(t / (alpha + 1))
  R <- 1
  while (F(z + R) - log(kappa(alpha1, z + R) + kappa(alpha2, R)) > ref) R <- 2 * R
  L <- 1
  while (F(-L) - log(kappa(alpha1, L) + kappa(alpha2, z + L)) > ref) L <- 2 * L
  geo <- function(s, M) if (M > s) s * 2^(0:floor(log2(M / s))) else numeric(0)
  c1 <- min(alpha1, 1) / 4
  c2 <- min(alpha2, 1) / 4
  pts <- c(-L, -geo(c1, L), 0, geo(c1, z / 2), z - geo(c2, z / 2), if (z > 0) z, z + geo(c2, R), z + R)
  pts <- sort(unique(pts[pts >= -L & pts <= z + R]))
  lo <- pts[-length(pts)]
  hi <- pts[-1]
  fl <- F(lo)
  fh <- F(hi)
  cel <- matrix(numeric(0), ncol = 4)
  repeat {
    w <- hi - lo
    sl <- (fh - fl) / w
    D <- meixner_dlogE_range(alpha1, alpha2, z, lo, hi)
    dp <- pmax(D$hi - sl, 0)
    dm <- pmax(sl - D$lo, 0)
    del <- ifelse(dp + dm > 0, dp * dm * w / (dp + dm), 0)
    mid <- (lo + hi) / 2
    cut <- 2 * del > log(rho) & pmax(fl, fh) + del > ref - log(z + L + R) & w > 1e-13 & mid > lo & mid < hi
    cel <- rbind(cel, cbind(lo[!cut], hi[!cut], fl[!cut] + del[!cut], sl[!cut]))
    if (!any(cut)) break
    lo <- lo[cut]
    hi <- hi[cut]
    fl <- fl[cut]
    fh <- fh[cut]
    mid <- (lo + hi) / 2
    fm <- F(mid)
    lo <- c(lo, mid)
    hi <- c(mid, hi)
    fl <- c(fl, fm)
    fh <- c(fm, fh)
  }
  w <- cel[, 2] - cel[, 1]
  sl <- cel[, 4]
  x <- sl * w
  lm <- cel[, 3] + log(w)
  ip <- x > 1e-12
  im <- x < -1e-12
  lm[ip] <- cel[ip, 3] + x[ip] + log(-expm1(-x[ip])) - log(sl[ip])
  lm[im] <- cel[im, 3] + log(-expm1(x[im])) - log(-sl[im])
  kR <- kappa(alpha1, z + R) + kappa(alpha2, R)
  kL <- kappa(alpha1, L) + kappa(alpha2, z + L)
  lm <- c(lm, F(z + R) - log(kR), F(-L) - log(kL))
  list(z = z, L = L, R = R, kL = kL, kR = kR, FL = F(-L), FR = F(z + R), cel = cel, prob = exp(lm - max(lm)))
}

# n propositions selon l'enveloppe, avec la valeur du log de l'enveloppe en chaque proposition
meixner_envelope_draw <- function(env, n) {
  nc <- nrow(env$cel)
  j <- sample.int(nc + 2, n, replace = TRUE, prob = env$prob)
  u <- stats::runif(n)
  y <- numeric(n)
  le <- numeric(n)
  ic <- j <= nc
  if (any(ic)) {
    cc <- env$cel[j[ic], , drop = FALSE]
    w <- cc[, 2] - cc[, 1]
    sl <- cc[, 4]
    x <- sl * w
    uu <- u[ic]
    t <- uu * w
    ip <- x > 1e-12
    im <- x < -1e-12
    t[ip] <- w[ip] + log(uu[ip] + (1 - uu[ip]) * exp(-x[ip])) / sl[ip]
    t[im] <- log1p(uu[im] * expm1(x[im])) / sl[im]
    y[ic] <- cc[, 1] + t
    le[ic] <- cc[, 3] + sl * t
  }
  iR <- j == nc + 1
  if (any(iR)) {
    e <- stats::rexp(sum(iR)) / env$kR
    y[iR] <- env$z + env$R + e
    le[iR] <- env$FR - env$kR * e
  }
  iL <- j == nc + 2
  if (any(iL)) {
    e <- stats::rexp(sum(iL)) / env$kL
    y[iL] <- -env$L - e
    le[iL] <- env$FL - env$kL * e
  }
  list(y = y, le = le)
}

# n tirages de la loi réduite pour z >= 0, avec une seule enveloppe
rmeixner_bridge_envelope <- function(alpha1, alpha2, z, n = 1) {
  env <- meixner_envelope(alpha1, alpha2, z)
  y <- numeric(0)
  while (length(y) < n) {
    k <- ceiling(1.5 * (n - length(y))) + 2
    prop <- meixner_envelope_draw(env, k)
    ok <- meixner_accept(log(stats::runif(k)) + prop$le, alpha1, alpha2, z, prop$y)
    y <- c(y, prop$y[ok])
  }
  y[seq_len(n)]
}
