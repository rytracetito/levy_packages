#' Inverse Gaussian bridge
#'
#' Draws the value at time \eqn{t} of an inverse Gaussian process known at two
#' times \eqn{t_A < t < t_B}, from the distribution of the inverse Gaussian
#' bridge. This function is used by [ig_resampling()] at each point of the
#' grid.
#'
#' @section Model:
#' The inverse Gaussian process with parameters \eqn{\mu > 0} and
#' \eqn{\eta > 0} has \eqn{X_0 = 0}, independent increments and
#' \eqn{X_t - X_s \sim IG(\mu (t - s), \eta (t - s)^2)} (mean
#' \eqn{\mu (t - s)}, shape \eqn{\eta (t - s)^2}). Its trajectories are
#' increasing. The vector of its parameters is \eqn{\Theta = (\mu, \eta)}.
#'
#' @section Distribution of the bridge:
#' Let \eqn{x_B > x_A}. Given \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B},
#' the ratio \eqn{W = (X_t - x_A) / (x_B - X_t)} has the density given by
#' Ribeiro and Webber (2003, eq. (27)), which is the mixture
#' \deqn{W \sim IG(m, \lambda) \text{ with probability } \frac{t_B - t}{t_B - t_A}, \qquad
#'   W \sim \frac{m^2}{IG(m, \lambda)} \text{ with probability } \frac{t - t_A}{t_B - t_A},}
#' where \eqn{m = (t - t_A) / (t_B - t)} and
#' \eqn{\lambda = \eta (t - t_A)^2 / (x_B - x_A)}; then
#' \eqn{X_t = x_A + (x_B - x_A)\, W / (1 + W)}. This distribution does not
#' depend on \eqn{\mu}. When \eqn{x_B = x_A} (two equal levels, which can
#' happen after rounding), the bridge is constant.
#'
#' The variable \eqn{W} is drawn from one Gaussian and one uniform variable,
#' by the method of Michael, Schucany and Haas (1976) as adapted by Ribeiro
#' and Webber (2003, sections 3.3 and 3.4).
#'
#' @inheritParams gamma_bridge
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(mu = , eta = )`, for instance the estimate returned by
#'   `ig_mle()` of the package `levyfit`. The bridge uses only `eta`.
#'
#' @inherit gamma_bridge return
#'
#' @references
#' Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
#' variates using transformations with multiple roots. *The American
#' Statistician*, 30(2), 88-90.
#' \doi{10.1080/00031305.1976.10479147}
#'
#' Ribeiro, C. and Webber, N. (2003). A Monte Carlo method for the normal
#' inverse Gaussian option valuation model using an inverse Gaussian bridge.
#' *Computing in Economics and Finance 2003*, 5.
#'
#' @seealso [ig_resampling()] to resample trajectories on a regular grid;
#'   [gamma_bridge()], [wiener_bridge()], [vg_bridge()], [meixner_bridge()].
#'
#' @examples
#' # 1. One draw of X_2, given X_1 = 3 and X_4 = 6 (mu = 1, eta = 3)
#' set.seed(1)
#' ig_bridge(t_A = 1, x_A = 3, t_B = 4, x_B = 6, t = 2, theta = c(mu = 1, eta = 3))
#'
#' # 2. Many draws: the values stay between 3 and 6
#' y <- ig_bridge(1, 3, 4, 6, 2, theta = c(mu = 1, eta = 3), n = 10000)
#' range(y)
#' hist(y, breaks = 50, freq = FALSE)
#'
#' # 3. Equal levels give a constant bridge; a decrease is impossible
#' ig_bridge(0, 1, 2, 1, 0.5, theta = c(mu = 1, eta = 3))
#' try(ig_bridge(0, 2, 2, 1, 0.5, theta = c(mu = 1, eta = 3)))
#'
#' # Errors: t outside (t_A, t_B), theta without all the parameters
#' try(ig_bridge(1, 3, 4, 6, 5, theta = c(mu = 1, eta = 3)))
#' try(ig_bridge(1, 3, 4, 6, 2, theta = c(eta = 3)))
#' @examplesIf requireNamespace("levyfit", quietly = TRUE)
#' # 4. Parameters estimated with levyfit on a simulated trajectory
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' x <- DegradSim::sim_ig(t, mu = 1, eta = 3)
#' est <- levyfit::ig_mle(x)
#' ig_bridge(1.2, x[2, 2], 2, x[3, 2], 1.5, theta = est, n = 3)
#' @export
ig_bridge <- function(t_A, x_A, t_B, x_B, t, theta, n = 1) {
  p <- check_theta(theta, c("mu", "eta"))
  bridge_draws(t_A, x_A, t_B, x_B, t, n, ig_draw(p[["eta"]]), "inverse Gaussian", increasing = TRUE)
}

#' Inverse Gaussian process: resampling on a regular grid
#'
#' Resamples one or several trajectories of an inverse Gaussian process onto
#' a regular time grid: each grid point that is not an inspection time is
#' drawn with [ig_bridge()] between the known levels around it. The grids and
#' the resampling algorithm are described in [levyBridge-package].
#'
#' @inheritParams gamma_resampling
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(mu = , eta = )`, for instance the estimate returned by
#'   `ig_mle()` of the package `levyfit`. The bridge uses only `eta`; `mu` is
#'   used by the option 2 of the method `"MultiBridge"`, to extend the
#'   trajectories.
#'
#' @inherit gamma_resampling return
#'
#' @seealso [ig_bridge()] for the distribution of the bridge,
#'   [levyBridge-package] for the grids and the algorithm, [plot_bridge()].
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (mu = 1, eta = 3)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- DegradSim::sim_ig(t, mu = 1, eta = 3)
#' theta <- c(mu = 1, eta = 3)
#' r <- ig_resampling(x, theta)
#' r
#'
#' # The resampled trajectory, on the regular grid 0, 2.25, 4.5, ..., 22.5
#' r[[1]]
#' attr(r, "teq")
#' plot_bridge(r)
#'
#' # 2. One trajectory given as a data frame (times, then levels): the grid is
#' #    0, 2, 4, 6, 8, 10, and the inspection times 2 and 6 keep their levels
#' x2 <- data.frame(t = c(0, 2, 3, 6, 7, 10), X = c(0, 1.4, 2.0, 4.1, 4.6, 6.8))
#' ig_resampling(x2, theta)[[1]]
#'
#' # Equal observed levels give a constant bridge; a decrease is impossible
#' ig_resampling(data.frame(t = c(0, 1, 2.5, 3), X = c(0, 1, 1, 2)), theta)[[1]]
#' try(ig_resampling(data.frame(t = c(0, 1, 2, 3), X = c(0, 2, 1.5, 3)), theta))
#'
#' # 3. Method "ExtBridge": the times have 1 decimal place, so the grid
#' #    0, 0.1, ..., 22.5 contains all of them
#' r3 <- ig_resampling(x, theta, method = "ExtBridge")
#' attr(r3, "teq")
#' all(t %in% r3[[1]]$t)
#' plot_bridge(r3)
#'
#' # digits = 0 truncates the times to integers, with a warning
#' attr(ig_resampling(x, theta, method = "ExtBridge", digits = 0), "teq")
#'
#' # 4. Five trajectories with the same inspection times (a matrix)
#' xm <- DegradSim::sim_ig(t, mu = 1, eta = 3, n_traj = 5)
#' head(xm)
#' r4 <- ig_resampling(xm, theta)
#' attr(r4, "teq")
#' plot_bridge(r4, trajectory = 1:5)
#'
#' # 5. Three trajectories, each with its own inspection times (a list)
#' plan <- DegradSim::inspection_times("decreasing", n = 10, n_traj = 3,
#'                                     tmax = c(20, 30, 40), digits = 1)
#' xs <- DegradSim::sim_ig(plan, mu = 1, eta = 3)
#'
#' # One grid and one step per trajectory
#' attr(ig_resampling(xs, theta), "teq")
#' attr(ig_resampling(xs, theta, method = "ExtBridge"), "teq")
#'
#' # The same step for all the trajectories
#' attr(ig_resampling(xs, theta, method = "ExtBridge", common_step = TRUE), "teq")
#' r5 <- ig_resampling(xs, theta, method = "MultiBridge")
#' attr(r5, "teq")
#' sapply(r5, nrow)
#' plot_bridge(r5, trajectory = 1:3)
#'
#' # Option 2: all the trajectories up to 40; beyond its last inspection, each
#' #    one is extended by increments of the process (this uses mu)
#' r6 <- ig_resampling(xs, theta, method = "MultiBridge", option = 2)
#' sapply(r6, nrow)
#'
#' # 6. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- ig_resampling(x, theta)[[1]]$X
#' set.seed(2); y2 <- ig_resampling(x, theta)[[1]]$X
#' identical(y1, y2)
#' @examplesIf requireNamespace("levyfit", quietly = TRUE)
#' # 7. With the parameters estimated by levyfit
#' theta_hat <- levyfit::ig_mle(xs)
#' theta_hat
#' r7 <- ig_resampling(xs, theta_hat)
#' attr(r7, "parameters")
#' @export
ig_resampling <- function(x, theta, method = "Bridge", digits = NULL, common_step = FALSE, option = 1) {
  p <- check_theta(theta, c("mu", "eta"))
  # IG(mu s, eta s^2) par la méthode de Michael, Schucany et Haas ; petite racine = m^2 / grande racine
  sim <- function(s) {
    m <- p[["mu"]] * s
    lam <- p[["eta"]] * s^2
    q <- stats::rnorm(1)^2
    w2 <- m + m^2 * q / (2 * lam) + m / (2 * lam) * sqrt(4 * m * lam * q + m^2 * q^2)
    w1 <- m^2 / w2
    if (stats::runif(1) <= m / (m + w1)) w1 else w2
  }
  bridge_paths(x, method, ig_draw(p[["eta"]]), "inverse Gaussian", p, increasing = TRUE, digits = digits,
               common_step = common_step, option = option, sim = sim)
}

# Pont inverse gaussien (Ribeiro et Webber, 2003) : accroissements à gauche et à droite de t
ig_draw <- function(eta) {
  function(t_inf, delta, t, t_sup) {
    if (delta == 0) return(c(0, 0))
    tau1 <- t - t_inf
    tau2 <- t_sup - t
    m <- tau1 / tau2
    lam <- eta * tau1^2 / delta
    q <- stats::rnorm(1)^2
    w2 <- m + m^2 * q / (2 * lam) + m / (2 * lam) * sqrt(4 * m * lam * q + m^2 * q^2)
    w1 <- m^2 / w2
    w <- if (stats::runif(1) <= tau1 * (1 + w1) / ((tau1 + tau2) * (m + w1))) w1 else w2
    delta * c(w / (1 + w), 1 / (1 + w))
  }
}
