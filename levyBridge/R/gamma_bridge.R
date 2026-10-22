#' Gamma bridge: resampling on a regular grid
#'
#' Resamples one or several trajectories of a homogeneous gamma process,
#' observed at irregular inspection times, onto a regular time grid, by
#' exact simulation of the gamma bridge between the known values.
#'
#' @section Model:
#' The homogeneous gamma process with shape \eqn{a > 0} and rate \eqn{b > 0}
#' has \eqn{X_0 = 0}, independent increments and
#' \eqn{X_t - X_s \sim \Gamma(a (t - s), b)} for \eqn{s < t} (shape
#' \eqn{a (t - s)}, rate \eqn{b}). Its trajectories are increasing.
#'
#' @section Law of the bridge:
#' Let \eqn{t_A < t < t_B}, \eqn{\tau_1 = t - t_A}, \eqn{\tau_2 = t_B - t} and
#' \eqn{\Delta = x_B - x_A}. Given \eqn{X_{t_A} = x_A} and
#' \eqn{X_{t_B} = x_B},
#' \deqn{X_t = x_A + \Delta\, Y, \qquad Y \sim \mathrm{Beta}(a \tau_1, a \tau_2).}
#' The law of the bridge does not depend on the rate \eqn{b}: only the shape
#' \eqn{a} is needed. When \eqn{\Delta = 0} (two equal observed values, which
#' can happen after rounding), the bridge is constant.
#'
#' @section Resampling algorithm (method "Bridge"):
#' This is the interpolation algorithm of Grall-Maës (2012), used for all the
#' processes of the package. Let \eqn{0 = t_0 < t_1 < \dots < t_N} be the
#' inspection times of a trajectory and \eqn{x_0, \dots, x_N} the observed
#' levels. The regular grid is \eqn{t'_i = i\, \Delta t_{eq}},
#' \eqn{i = 0, \dots, N}, with \eqn{\Delta t_{eq} = t_N / N}: it has as many
#' points as observations and the same end points.
#' 1. \eqn{x'_0 = x_0} and \eqn{x'_N = x_N}.
#' 2. For \eqn{i = 1, \dots, N - 1}, in increasing order:
#'    * if \eqn{t'_i} is an inspection time, \eqn{x'_i} is the observed value;
#'    * otherwise, \eqn{x'_i} is drawn from the bridge between the nearest
#'      known point on the left (an inspection or the grid point
#'      \eqn{t'_{i-1}} already drawn, whichever is closer) and the nearest
#'      inspection on the right.
#'
#' By the Markov property of the bridge, the vector
#' \eqn{(x'_1, \dots, x'_{N-1})} follows exactly the conditional law of
#' \eqn{(X_{t'_1}, \dots, X_{t'_{N-1}})} given the observations. Hence, if the
#' data come from the model and the parameters are the true ones, the
#' reconstructed increments \eqn{x'_i - x'_{i-1}} (see [increments()]) are
#' independent with the law of \eqn{X_{\Delta t_{eq}}}: they can be tested
#' with an i.i.d. goodness-of-fit test. A grid point that equals an inspection
#' time up to rounding errors (relative difference below \eqn{10^{-9}}) is
#' set to that inspection time.
#'
#' **Precision.** Each draw splits the increment \eqn{\Delta} of the current
#' interval into a left part (up to the grid point) and a right part (up to
#' the next inspection), both computed directly: here \eqn{\Delta\, Y} and
#' \eqn{\Delta (1 - Y)}, with \eqn{Y = G_1 / (G_1 + G_2)},
#' \eqn{G_1 \sim \Gamma(a \tau_1)}, \eqn{G_2 \sim \Gamma(a \tau_2)} drawn on
#' the log scale. The reconstructed increments ([increments()]) are these
#' parts, not differences of levels: a very small increment keeps all its
#' significant digits instead of being rounded to 0, even for very small
#' shapes \eqn{a \tau}.
#'
#' @section Method "ExtBridge":
#' Extension of the method of Grall-Maës whose regular grid contains all the
#' inspection times (Tito and Paroissin). The inspection times are multiples of
#' \eqn{10^{-d}}, where \eqn{d} is the number of decimal places (argument
#' `digits`), and the time step is the largest one such that every inspection
#' time is a grid point:
#' \deqn{h = 10^{-d}\, \mathrm{gcd}\big(10^d t_1, \dots, 10^d t_N\big), \qquad
#'   t'_k = k\, h, \quad k = 0, \dots, t_N / h.}
#' The grid is then covered as in the method "Bridge": every inspection time
#' keeps its observed value, and the other grid points are drawn from the
#' bridge between the last known value on the left and the next inspection.
#' The reconstructed increments are independent with the law of \eqn{X_h},
#' and the pseudo-sample contains all the observations. The grid has
#' \eqn{t_N / h \ge N} steps, with equality if and only if the inspection
#' times are regular; then the method returns the observations themselves.
#'
#' With `digits = NULL` (default), \eqn{d} is the smallest number of decimal
#' places, between 0 and 6, such that all the inspection times are multiples
#' of \eqn{10^{-d}}. If there is none (for instance times drawn with
#' `runif()`), the times are rounded to 2 decimal places, with a warning; give
#' `digits` to choose another rounding. When `digits` is given and some times
#' have more decimal places, they are rounded to `digits` decimal places, with
#' a warning; two times equal after rounding give an error. A warning is also
#' given when the grid has more than \eqn{10^5} steps, since the resampling
#' then takes time. The grid and the time step are computed for each
#' trajectory; with `common_step = TRUE`, \eqn{d} and \eqn{h} are computed
#' from the inspection times of all the trajectories, so that all the grids
#' have the same step.
#'
#' @section Method "MultiBridge":
#' Version of the method "Bridge" for several trajectories, with **one
#' regular grid common to all the trajectories** (Tito and Paroissin). With
#' the method "Bridge", each trajectory has its own time step
#' \eqn{t_N / N}: the reconstructed increments of two trajectories then have
#' different laws and cannot be tested together. With the method
#' "MultiBridge", all the reconstructed increments have the law of \eqn{X_h}
#' for the same \eqn{h}: they form a single i.i.d. sample under the model,
#' and one goodness-of-fit test can be applied to all the trajectories
#' together.
#'
#' **Time step.** For \eqn{n} trajectories, the trajectory \eqn{i} being
#' observed on \eqn{[0, t_{i, max}]} with \eqn{m_i} increments,
#' \deqn{h = \frac{t_{1, max} + \dots + t_{n, max}}{m_1 + \dots + m_n},}
#' the total observed duration divided by the total number of increments.
#' The number of reconstructed increments is then about the number of
#' observed increments, as with the method "Bridge"; for one trajectory,
#' \eqn{h = t_N / N} is the step of the method "Bridge". The grid is
#' \eqn{t'_k = k h}; as with the method "Bridge", it does not necessarily
#' contain the inspection times (a grid point equal to an inspection time
#' keeps the observed value).
#'
#' **Option 1** (`option = 1`, default): each trajectory is resampled up to
#' its own last inspection, on the grid points \eqn{k h \le t_{i, max}}. Each
#' grid point is drawn from the bridge between the last known value on the
#' left and the next inspection, exactly as with the method "Bridge". The end
#' of the trajectory, between its last grid point and \eqn{t_{i, max}}, is not
#' used. A long trajectory gives many reconstructed increments, a short one
#' few (none when \eqn{t_{i, max} < h}, with a warning). Under the model, with
#' the true parameters, all the reconstructed increments are independent with
#' the law of \eqn{X_h}: the pooled sample is exactly i.i.d. This is the
#' recommended option.
#'
#' **Option 2** (`option = 2`): all the trajectories are resampled up to the
#' last inspection time over all the trajectories, \eqn{\max_i t_{i, max}}, so
#' that they all have the same number of grid points. Up to
#' \eqn{t_{i, max}}, as in option 1; beyond, there is no observation left and
#' the trajectory is extended by increments of the process drawn with the
#' given parameters (not a bridge): the first one starts from the last
#' observation \eqn{x_{i, N}}. This needs all the parameters of the process
#' (for instance \eqn{b} for the gamma process, see the arguments). Under the
#' model, with the true parameters, the pooled sample is also i.i.d. with the
#' law of \eqn{X_h}; but the extended values do not come from the data: they
#' follow the null model by construction and dilute the information of the
#' observations, which lowers the power of a test, all the more as the
#' trajectories have different durations. Use it only for trajectories of
#' similar durations.
#'
#' **When to use "MultiBridge".** For several trajectories observed at
#' different times, when one test on all the trajectories together is
#' wanted. With one trajectory, both options give the method "Bridge". The
#' method "ExtBridge" with `common_step = TRUE` also gives a common step, with
#' a grid that contains all the inspection times, but it can have many more
#' grid points.
#'
#' @param x Degradation levels of one trajectory; or a list of data frames,
#'   one per trajectory, whose first column contains the inspection times and
#'   second column the degradation levels (for instance the result of a
#'   `levysim` simulation function).
#' @param t Inspection times of the trajectory, starting at 0 and strictly
#'   increasing. Not used when `x` is a list.
#' @param a Shape parameter of the gamma process, \eqn{a > 0}, for instance
#'   estimated with the function `gamma_mle()` of the package `levyfit`.
#' @param method `"Bridge"` (default) for the algorithm of Grall-Maës (2012),
#'   with as many grid points as observations and one grid per trajectory;
#'   `"ExtBridge"` for its extension, whose grid contains all the inspection
#'   times (see the section Method "ExtBridge"); `"MultiBridge"` for several
#'   trajectories on one common grid (see the section Method "MultiBridge").
#' @param digits Used only with `method = "ExtBridge"`: number of decimal
#'   places of the inspection times, which sets the time step of the grid.
#'   `NULL` (default) finds it automatically.
#' @param common_step Used only with `method = "ExtBridge"` and several
#'   trajectories: `FALSE` (default) gives each trajectory its own time step
#'   \eqn{h}; `TRUE` gives all the trajectories the same time step, the
#'   largest one such that every inspection time of every trajectory is a grid
#'   point. With a common step, the reconstructed increments of all the
#'   trajectories have the same law under the model and can be tested
#'   together.
#' @param option Used only with `method = "MultiBridge"`: `1` (default) to
#'   resample each trajectory up to its own last inspection, `2` to resample
#'   all the trajectories up to the last inspection time over all the
#'   trajectories, extending the shorter ones by increments of the process (see
#'   the section Method "MultiBridge").
#' @param b Rate parameter of the gamma process, \eqn{b > 0}; used only with
#'   `method = "MultiBridge"` and `option = 2`, to extend the trajectories.
#'
#' @return An object of class `"levybridge"`: a list with one data frame per
#'   trajectory, with columns `t` (regular grid) and `X` (resampled levels).
#'   Its attributes are `model`, `parameters`, `method`, `teq` (time step of
#'   the grid of each trajectory) and `observed` (the observed data, one data
#'   frame per trajectory; with the method "ExtBridge", the times are the
#'   rounded ones when a rounding was needed). Use [increments()] for the reconstructed
#'   increments and [plot_bridge()] to compare observed and resampled values.
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' Avramidis, A. N., L'Ecuyer, P. and Tremblay, P.-A. (2003). Efficient
#' simulation of gamma and variance-gamma processes. *Proceedings of the 2003
#' Winter Simulation Conference*, 319-326.
#'
#' @seealso
#' * [increments()]: the reconstructed increments;
#' * [plot_bridge()]: observed and resampled values;
#' * [wiener_bridge()]: bridge of the Wiener process;
#' * [ig_bridge()]: bridge of the inverse Gaussian process;
#' * [vg_bridge()]: bridge of the variance-gamma process;
#' * [meixner_bridge()]: bridge of the Meixner process.
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (a = 2, b = 3)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- c(0, cumsum(rgamma(10, shape = 2 * diff(t), rate = 3)))
#' set.seed(1)
#' r <- gamma_bridge(x, t, a = 2)
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
#' x2 <- c(0, 1.4, 2.0, 4.1, 4.6, 6.8)
#' gamma_bridge(x2, t2, a = 2)[[1]]
#'
#' # 5. One trajectory given as a data frame (times, then levels)
#' gamma_bridge(data.frame(t = t2, X = x2), a = 2)
#'
#' # 6. Several trajectories, each with its own inspection times:
#' #    one regular grid and one time step per trajectory
#' liste <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.8, 2.1, 4.3)),
#'               data.frame(t = c(0, 2, 4, 5, 9), X = c(0, 1.2, 2.9, 3.1, 6.2)))
#' r2 <- gamma_bridge(liste, a = 2)
#' r2
#' attr(r2, "teq")
#' increments(r2)
#' plot_bridge(r2, trajectory = 1:2)
#'
#' # 7. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- gamma_bridge(x, t, a = 2)[[1]]$X
#' set.seed(2); y2 <- gamma_bridge(x, t, a = 2)[[1]]$X
#' identical(y1, y2)
#'
#' # 8. Equal observed values (after rounding) give a constant bridge,
#' #    a decrease is impossible for this increasing process
#' gamma_bridge(c(0, 1, 1, 2), c(0, 1, 2.5, 3), a = 2)[[1]]
#' try(gamma_bridge(c(0, 2, 1.5, 3), c(0, 1, 2, 3), a = 2))
#'
#' # 9. Method "ExtBridge": the times have 1 decimal place, the gcd of the
#' #    gaps is 0.1, and the grid 0, 0.1, ..., 22.5 contains all of them
#' r3 <- gamma_bridge(x, t, a = 2, method = "ExtBridge")
#' r3
#' attr(r3, "teq")
#' nrow(r3[[1]])
#' all(t %in% r3[[1]]$t)
#' plot_bridge(r3)
#'
#' # Times with 2 decimal places give a grid of step 0.01, here 0.05
#' t4 <- c(0, 0.35, 1.2, 1.85, 2.5)
#' x4 <- c(0, 0.3, 1.1, 1.4, 2.2)
#' attr(gamma_bridge(x4, t4, a = 2, method = "ExtBridge"), "teq")
#'
#' # digits = 1 rounds the times to 1 decimal place, with a warning
#' attr(gamma_bridge(x4, t4, a = 2, method = "ExtBridge", digits = 1), "teq")
#'
#' # Regular times: the grid is the inspection times, the data are unchanged
#' gamma_bridge(c(0, 1.2, 2.9, 3.5), c(0, 2, 4, 6), a = 2, method = "ExtBridge")[[1]]
#'
#' # Several trajectories with a common step: gcd of all the gaps, here 0.25
#' liste2 <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.8, 2.1, 4.3)),
#'                data.frame(t = c(0, 0.5, 1.25, 2), X = c(0, 0.4, 1.1, 1.6)))
#' attr(gamma_bridge(liste2, a = 2, method = "ExtBridge"), "teq")
#' attr(gamma_bridge(liste2, a = 2, method = "ExtBridge", common_step = TRUE), "teq")
#'
#' # Very small shapes: the reconstructed increments are tiny but never 0
#' min(increments(gamma_bridge(x, t, a = 0.05, method = "ExtBridge"))$dx)
#'
#' # Method "MultiBridge": four trajectories with their own times, one common grid
#' liste_t <- list(c(0, 3, 8), c(0, 2, 5, 9, 12), c(0, 4, 7, 11, 15, 18, 22, 27, 33, 39, 45),
#'                 c(0, 2, 3, 5, 8, 9, 14, 15, 17, 20, 22, 25, 27, 30, 33, 46, 58))
#' liste4 <- lapply(liste_t, function(tk) data.frame(t = tk, X = c(0, cumsum(rgamma(length(tk) - 1, 2 * diff(tk), 3)))))
#'
#' # Option 1: each trajectory up to its last inspection;
#' # step h = (8 + 12 + 45 + 58) / (2 + 4 + 10 + 16) = 3.84 for all
#' r5 <- gamma_bridge(liste4, a = 2, method = "MultiBridge")
#' r5
#' attr(r5, "teq")
#' sapply(r5, nrow)
#' plot_bridge(r5, trajectory = 1:4)
#'
#' # The reconstructed increments of the four trajectories have the same law Gamma(a h, b):
#' # one sample for one test
#' d5 <- increments(r5)
#' ks.test(d5$dx, "pgamma", shape = 2 * attr(r5, "teq")[1], rate = 3)
#'
#' # Option 2: all the trajectories up to 58, the shorter ones extended by gamma
#' # increments, which needs the rate b
#' r6 <- gamma_bridge(liste4, a = 2, method = "MultiBridge", option = 2, b = 3)
#' sapply(r6, nrow)
#' try(gamma_bridge(liste4, a = 2, method = "MultiBridge", option = 2))
#'
#' # With one trajectory, "MultiBridge" is the method "Bridge"
#' attr(gamma_bridge(x, t, a = 2, method = "MultiBridge"), "teq")
#'
#' # 10. Error: the first time must be 0
#' try(gamma_bridge(x, t + 1, a = 2))
#' @examplesIf requireNamespace("levysim", quietly = TRUE) && requireNamespace("levyfit", quietly = TRUE)
#' # 11. Ten trajectories simulated with levysim, inspected more and more often
#' plan <- levysim::inspection_times("increasing", n = 20, n_traj = 10, tmax = 40)
#' sims <- levysim::sim_gamma(plan, a = 2, b = 3)
#'
#' # Parameters estimated with levyfit
#' a_hat <- levyfit::gamma_mle(sims)[["a"]]
#' r <- gamma_bridge(sims, a = a_hat)
#' r
#' plot_bridge(r, trajectory = 1:4)
#'
#' # Reconstructed increments of all the trajectories
#' d <- increments(r)
#' head(d)
#' # Under the model, the increments follow Gamma(a * teq, b)
#' ks.test(d$dx, "pgamma", shape = 2 * d$teq[1], rate = 3)
#' @export
gamma_bridge <- function(x, t = NULL, a, method = "Bridge", digits = NULL, common_step = FALSE, option = 1,
                         b = NULL) {
  check_par(a = a)
  draw <- function(t_inf, delta, t, t_sup) delta * rbeta_pair(a * (t - t_inf), a * (t_sup - t))
  sim <- NULL
  if (!is.null(b)) {
    check_par(b = b)
    sim <- function(s) exp(log_rgamma(a * s) - log(b))
  }
  bridge_paths(x, t, method, draw, "gamma", c(a = unname(a)), increasing = TRUE, digits = digits,
               common_step = common_step, option = option, sim = sim)
}
