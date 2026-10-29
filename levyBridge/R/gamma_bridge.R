#' Gamma bridge
#'
#' Draws the value at time \eqn{t} of a gamma process known at two times
#' \eqn{t_A < t < t_B}, from the distribution of the gamma bridge. This
#' function is used by [gamma_resampling()] at each point of the grid.
#'
#' @section Model:
#' The gamma process with shape \eqn{a > 0} and rate \eqn{b > 0} has
#' \eqn{X_0 = 0}, independent increments and
#' \eqn{X_t - X_s \sim \Gamma(a (t - s), b)} for \eqn{s < t} (shape
#' \eqn{a (t - s)}, rate \eqn{b}). Its trajectories are increasing. The
#' vector of its parameters is \eqn{\Theta = (a, b)}.
#'
#' @section Distribution of the bridge:
#' Given \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B},
#' \deqn{X_t = x_A + (x_B - x_A)\, Y, \qquad Y \sim \mathrm{Beta}\big(a (t - t_A),\ a (t_B - t)\big)}
#' (Grall-Maës, 2012, eq. (3)). This distribution does not depend on the rate
#' \eqn{b}. When \eqn{x_B = x_A} (two equal levels, which can happen after
#' rounding), the bridge is constant.
#'
#' @param t_A,x_A Time and level of the known point on the left.
#' @param t_B,x_B Time and level of the known point on the right,
#'   \eqn{t_B > t_A}.
#' @param t Time of the draw, \eqn{t_A < t < t_B}.
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(a = , b = )`, for instance the estimate returned by
#'   `gamma_mle()` of the package `levyfit`. The bridge uses only `a`.
#' @param n Number of independent draws, a positive integer (default 1).
#'
#' @return A numeric vector of `n` independent draws of \eqn{X_t} given
#'   \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B}.
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' @seealso [gamma_resampling()] to resample trajectories on a regular grid;
#'   [wiener_bridge()], [ig_bridge()], [vg_bridge()], [meixner_bridge()].
#'
#' @examples
#' # 1. One draw of X_2, given X_1 = 3 and X_4 = 6 (a = 2, b = 3)
#' set.seed(1)
#' gamma_bridge(t_A = 1, x_A = 3, t_B = 4, x_B = 6, t = 2, theta = c(a = 2, b = 3))
#'
#' # 2. Many draws: X_2 = 3 + 3 Y with Y ~ Beta(2, 4)
#' y <- gamma_bridge(1, 3, 4, 6, 2, theta = c(a = 2, b = 3), n = 10000)
#' hist(y, breaks = 50, freq = FALSE)
#' curve(dbeta((x - 3) / 3, 2, 4) / 3, add = TRUE, col = "red")
#'
#' # The rate b does not change the bridge
#' set.seed(2); y1 <- gamma_bridge(1, 3, 4, 6, 2, theta = c(a = 2, b = 3), n = 5)
#' set.seed(2); y2 <- gamma_bridge(1, 3, 4, 6, 2, theta = c(a = 2, b = 100), n = 5)
#' identical(y1, y2)
#'
#' # 3. Equal levels give a constant bridge; a decrease is impossible
#' gamma_bridge(0, 1, 2, 1, 0.5, theta = c(a = 2, b = 3))
#' try(gamma_bridge(0, 2, 2, 1, 0.5, theta = c(a = 2, b = 3)))
#' @examplesIf requireNamespace("levyfit", quietly = TRUE)
#' # 4. Parameters estimated with levyfit on a simulated trajectory
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3)
#' est <- levyfit::gamma_mle(x)
#' est
#' # Three draws at t = 1.5, between the inspections 1.2 and 2
#' gamma_bridge(1.2, x[2, 2], 2, x[3, 2], 1.5, theta = est, n = 3)
#' @export
gamma_bridge <- function(t_A, x_A, t_B, x_B, t, theta, n = 1) {
  p <- check_theta(theta, c("a", "b"))
  bridge_draws(t_A, x_A, t_B, x_B, t, n, gamma_draw(p[["a"]]), "gamma", increasing = TRUE)
}

#' Gamma process: resampling on a regular grid
#'
#' Resamples one or several trajectories of a gamma process onto a regular
#' time grid: each grid point that is not an inspection time is drawn with
#' [gamma_bridge()] between the known levels around it. The grids and the
#' resampling algorithm are described in [levyBridge-package].
#'
#' @param x The observed trajectories: a data frame whose first column
#'   contains the inspection times (starting at 0, strictly increasing) and
#'   second column the degradation levels; a list of such data frames, one
#'   per trajectory; or a matrix whose first column contains the times and
#'   the other columns the levels, one column per trajectory. The results of
#'   the simulation functions of `DegradSim` are accepted. See the section
#'   Data of [levyBridge-package].
#' @param theta Vector \eqn{\Theta} of the parameters of the process, a named
#'   numeric vector `c(a = , b = )`, for instance the estimate returned by
#'   `gamma_mle()` of the package `levyfit`. The bridge uses only `a`; `b` is
#'   used by the option 2 of the method `"MultiBridge"`, to extend the
#'   trajectories.
#' @param method Grid of the resampling: `"Bridge"` (default), `"ExtBridge"`
#'   or `"MultiBridge"`, see the section Resampling methods of
#'   [levyBridge-package].
#' @param digits Used only with `method = "ExtBridge"`: number of decimal
#'   places of the inspection times, which sets the step of the grid. `NULL`
#'   (default) takes the number of decimal places of the times.
#' @param common_step Used only with `method = "ExtBridge"`: `FALSE`
#'   (default) gives each trajectory its own step, `TRUE` gives the same step
#'   to all the trajectories.
#' @param option Used only with `method = "MultiBridge"`: `1` (default) to
#'   resample each trajectory up to its own last inspection, `2` to resample
#'   all the trajectories up to the last inspection time over all of them,
#'   extending the shorter ones by increments of the process.
#'
#' @return An object of class `"levybridge"`: a list with one data frame per
#'   trajectory, with the columns `t` (regular grid) and `X` (resampled
#'   levels), and the attributes `teq` (step of each grid), `dx` (increments
#'   between consecutive grid points), `observed`, `model`, `parameters` and
#'   `method`. See the section Results of [levyBridge-package].
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' @seealso [gamma_bridge()] for the distribution of the bridge,
#'   [levyBridge-package] for the grids and the algorithm, [plot_bridge()].
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (a = 2, b = 3)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3)
#' r <- gamma_resampling(x, theta = c(a = 2, b = 3))
#' r
#'
#' # The resampled trajectory, on the regular grid 0, 2.25, 4.5, ..., 22.5
#' r[[1]]
#' attr(r, "teq")
#' plot_bridge(r)
#'
#' # 2. A grid point equal to an inspection time keeps the observed level:
#' #    the grid is 0, 2, 4, 6, 8, 10 and 2 and 6 are inspection times
#' x2 <- data.frame(t = c(0, 2, 3, 6, 7, 10), X = c(0, 1.4, 2.0, 4.1, 4.6, 6.8))
#' gamma_resampling(x2, theta = c(a = 2, b = 3))[[1]]
#'
#' # 3. Method "ExtBridge": the times have 1 decimal place, so the grid
#' #    0, 0.1, ..., 22.5 contains all of them
#' r3 <- gamma_resampling(x, theta = c(a = 2, b = 3), method = "ExtBridge")
#' attr(r3, "teq")
#' all(t %in% r3[[1]]$t)
#' plot_bridge(r3)
#'
#' # digits = 0 truncates the times to integers, with a warning
#' attr(gamma_resampling(x, theta = c(a = 2, b = 3), method = "ExtBridge", digits = 0), "teq")
#'
#' # 4. Three trajectories, each with its own inspection times
#' plan <- DegradSim::inspection_times("decreasing", n = 10, n_traj = 3,
#'                                     tmax = c(20, 30, 40), digits = 1)
#' xs <- DegradSim::sim_gamma(plan, a = 2, b = 3)
#'
#' # One grid per trajectory
#' attr(gamma_resampling(xs, theta = c(a = 2, b = 3)), "teq")
#' attr(gamma_resampling(xs, theta = c(a = 2, b = 3), method = "ExtBridge"), "teq")
#'
#' # The same step for all the trajectories
#' attr(gamma_resampling(xs, theta = c(a = 2, b = 3), method = "ExtBridge", common_step = TRUE), "teq")
#' r5 <- gamma_resampling(xs, theta = c(a = 2, b = 3), method = "MultiBridge")
#' attr(r5, "teq")
#' sapply(r5, nrow)
#' plot_bridge(r5, trajectory = 1:3)
#'
#' # Option 2: all the trajectories up to 40, extended with the rate b
#' r6 <- gamma_resampling(xs, theta = c(a = 2, b = 3), method = "MultiBridge", option = 2)
#' sapply(r6, nrow)
#'
#' # 5. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- gamma_resampling(x, theta = c(a = 2, b = 3))[[1]]$X
#' set.seed(2); y2 <- gamma_resampling(x, theta = c(a = 2, b = 3))[[1]]$X
#' identical(y1, y2)
#' @examplesIf requireNamespace("levyfit", quietly = TRUE)
#' # 6. With the parameters estimated by levyfit
#' r_est <- gamma_resampling(xs, theta = levyfit::gamma_mle(xs))
#' attr(r_est, "parameters")
#' @export
gamma_resampling <- function(x, theta, method = "Bridge", digits = NULL, common_step = FALSE, option = 1) {
  p <- check_theta(theta, c("a", "b"))
  sim <- function(s) exp(log_rgamma(p[["a"]] * s) - log(p[["b"]]))
  bridge_paths(x, method, gamma_draw(p[["a"]]), "gamma", p, increasing = TRUE, digits = digits,
               common_step = common_step, option = option, sim = sim)
}

# Pont gamma : accroissements à gauche et à droite de t, delta * (Y, 1 - Y) avec Y ~ Beta
gamma_draw <- function(a) {
  function(t_inf, delta, t, t_sup) delta * rbeta_pair(a * (t - t_inf), a * (t_sup - t))
}
