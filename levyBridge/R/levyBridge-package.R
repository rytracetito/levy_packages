#' levyBridge: stochastic bridges for Levy degradation processes
#'
#' `levyBridge` simulates exactly the bridges of Levy degradation processes,
#' and uses them to resample degradation trajectories, observed at regular or
#' irregular inspection times, onto a regular time grid.
#'
#' The package has two kinds of functions, for each process:
#' * the **bridge** functions `gamma_bridge()`, `wiener_bridge()`,
#'   `ig_bridge()`, `vg_bridge()` and `meixner_bridge()` draw the value
#'   \eqn{X_t} of the process at one time \eqn{t}, given its levels at two
#'   times \eqn{t_A < t < t_B};
#' * the **resampling** functions `gamma_resampling()`, `wiener_resampling()`,
#'   `ig_resampling()`, `vg_resampling()` and `meixner_resampling()` resample
#'   whole trajectories onto a regular grid, by calling the bridge function at
#'   each grid point.
#'
#' All these functions take the vector \eqn{\Theta} of all the parameters of
#' the process, as a named vector `theta`, for instance `c(a = 2, b = 3)` for
#' the gamma process.
#'
#' @section Bridge of a Levy process:
#' Let \eqn{X} be a Levy process such that \eqn{X_s} has the density
#' \eqn{f_s}, and let \eqn{t_A < t < t_B}. Given \eqn{X_{t_A} = x_A} and
#' \eqn{X_{t_B} = x_B}, the value \eqn{X_t} has the density
#' \deqn{f_{t | t_A, t_B}(x) = \frac{f_{t - t_A}(x - x_A)\, f_{t_B - t}(x_B - x)}{f_{t_B - t_A}(x_B - x_A)}}
#' (Hoyle, Hughston and Macrina, 2011). The bridge functions draw exactly
#' from this distribution:
#'
#' | Process | Bridge | Resampling | \eqn{\Theta} (`theta`) | Used by the bridge | Simulation of the bridge |
#' |---|---|---|---|---|---|
#' | gamma | [gamma_bridge()] | [gamma_resampling()] | `a`, `b` | `a` | beta variable |
#' | Wiener | [wiener_bridge()] | [wiener_resampling()] | `mu`, `sigma` | `sigma` | Gaussian variable |
#' | inverse Gaussian | [ig_bridge()] | [ig_resampling()] | `mu`, `eta` | `eta` | IG or scaled inverse IG |
#' | variance-gamma | [vg_bridge()] | [vg_resampling()] | `theta`, `sigma`, `nu` | all | GIG, beta, Gaussian |
#' | Meixner | [meixner_bridge()] | [meixner_resampling()] | `a`, `b`, `d`, `m` | `a`, `d`, `m` | acceptance-rejection |
#'
#' The distribution of the bridge does not always depend on all the
#' parameters: the parameters that it does not use are ignored by the bridge
#' functions, and used by the resampling functions only to extend the
#' trajectories (option 2 of the method `"MultiBridge"`). The help page of
#' each bridge function gives the distribution of its bridge and the way it
#' is drawn.
#'
#' @section Data:
#' The argument `x` of the resampling functions accepts three formats:
#' * **one trajectory**: a data frame whose first column contains the
#'   inspection times and second column the degradation levels;
#' * **several trajectories, each with its own inspection times**: a list of
#'   such data frames, one per trajectory, for instance the result of a
#'   simulation function of `DegradSim` when the trajectories have different
#'   inspection times;
#' * **several trajectories with the same inspection times**: a matrix whose
#'   first column contains the times and the other columns the levels, one
#'   column per trajectory, for instance the result of a simulation function
#'   of `DegradSim` when all the trajectories have the same inspection times.
#'
#' The columns are read by position, not by name. The inspection times of a
#' trajectory start at 0 and are strictly increasing.
#'
#' @section Resampling algorithm:
#' Let \eqn{0 = t_0 < t_1 < \dots < t_N} be the inspection times of a
#' trajectory and \eqn{x_0, \dots, x_N} its observed levels. The method
#' (section Resampling methods) gives a regular grid
#' \eqn{0 = t'_0 < t'_1 < \dots}, and the resampled levels \eqn{x'_k} are
#' computed in increasing order of \eqn{k}:
#' * if \eqn{t'_k} is an inspection time, \eqn{x'_k} is the observed level;
#' * otherwise, \eqn{x'_k} is drawn by the bridge function between the last
#'   known point on the left (the inspection time or the grid point
#'   \eqn{t'_{k-1}}, whichever is closer to \eqn{t'_k}) and the next
#'   inspection time on the right.
#'
#' This is the interpolation algorithm of Grall-Maës (2012), proposed for the
#' gamma process and used here for all the processes. Since a Levy process is
#' a Markov process, the value at \eqn{t'_k} given the observations and the
#' levels already drawn depends only on the two known points around
#' \eqn{t'_k}: the resampled levels follow the conditional distribution of
#' the process at the grid points given the observations.
#'
#' A grid point that differs from an inspection time only by a rounding
#' error of the computer (relative difference below \eqn{10^{-9}}) is set to
#' that inspection time.
#'
#' @section Resampling methods:
#' Every resampling function has an argument `method`, which chooses the grid.
#'
#' **`"Bridge"`** (default). The grid of Grall-Maës (2012):
#' \deqn{t'_k = k\, \frac{t_N}{N}, \qquad k = 0, \dots, N,}
#' with as many points as observations and the same end points. Each
#' trajectory has its own grid, of step \eqn{t_N / N}.
#'
#' **`"ExtBridge"`**. A grid that contains all the inspection times. The
#' times are taken with \eqn{d} decimal places (argument `digits`) and
#' \deqn{h = 10^{-d}\, \mathrm{gcd}\big(\lfloor 10^d t_1 \rfloor, \dots, \lfloor 10^d t_N \rfloor\big),
#'   \qquad t'_k = k\, h, \quad k = 0, \dots, t_N / h,}
#' where \eqn{\lfloor u \rfloor} is the largest integer not greater than
#' \eqn{u}. The step \eqn{h} is the largest one such that every inspection
#' time, truncated to \eqn{d} decimal places, is a grid point: all the observations are kept. The grid
#' has \eqn{t_N / h \ge N} steps, with equality if and only if the
#' inspection times are regular; the grid is then the inspection times and
#' the data are returned unchanged.
#'
#' With `digits = NULL` (default), \eqn{d} is the number of decimal places
#' of the inspection times, that is the smallest integer \eqn{d} such that
#' all the times are multiples of \eqn{10^{-d}}, up to rounding errors of the
#' computer. This package looks for \eqn{d} between 0 and 6, a limit chosen to
#' avoid very fine grids; if the times have more decimal places, an error
#' asks for `digits`. When `digits` is given and some times have more decimal
#' places, they are truncated to `digits` decimal places, with a warning;
#' two times equal after truncation give an error. A warning is also given when the
#' grid has more than \eqn{10^5} steps, since the resampling then takes time.
#'
#' Each trajectory has its own step \eqn{h}. With `common_step = TRUE`,
#' \eqn{d} and \eqn{h} are computed from the inspection times of all the
#' trajectories, which then have the same step.
#'
#' **`"MultiBridge"`**. One step common to all the trajectories: the
#' trajectory \eqn{i}, observed on \eqn{[0, t_{i, max}]} with \eqn{m_i}
#' increments, is resampled on \eqn{t'_k = k h} with
#' \deqn{h = \frac{t_{1, max} + \dots + t_{n, max}}{m_1 + \dots + m_n},}
#' the total observed duration divided by the total number of increments.
#' Grall-Maës (2012, eq. (7)) proposes this step for two trajectories. The
#' number of reconstructed increments is then about the number of observed
#' ones; for one trajectory, \eqn{h = t_N / N} and the method is `"Bridge"`.
#' The grid does not necessarily contain the inspection times.
#' * `option = 1` (default): the trajectory \eqn{i} is resampled on the grid
#'   points \eqn{k h \le t_{i, max}}, that is up to its own last inspection.
#'   The end of the trajectory, between its last grid point and
#'   \eqn{t_{i, max}}, is not used. A trajectory shorter than \eqn{h} gives
#'   no increment, with a warning.
#' * `option = 2`: all the trajectories are resampled up to the last
#'   inspection time over all of them, \eqn{\max_i t_{i, max}}. Beyond its
#'   last inspection, a trajectory is extended by increments of the process
#'   drawn with the given parameters, starting from its last observed level.
#'   This needs all the parameters of the process (for instance the rate `b`
#'   of the gamma process), and the extended values do not come from the
#'   data.
#'
#' **Choice of the method.**
#' * `"Bridge"`: one grid per trajectory, with as many points as
#'   observations; the steps of two trajectories are different in general.
#' * `"ExtBridge"`: all the observations are kept on the grid, which can have
#'   many more points than observations.
#' * `"MultiBridge"`, or `"ExtBridge"` with `common_step = TRUE`: the same
#'   step for all the trajectories.
#'
#' @section Results:
#' A resampling function returns an object of class `"levybridge"`, a list
#' with one data frame per trajectory:
#' * `res[[k]]` is the \eqn{k}-th trajectory on its regular grid, with the
#'   columns `t` (grid) and `X` (resampled levels);
#' * `attr(res, "teq")` gives the step of the grid of each trajectory;
#' * `attr(res, "dx")` gives, for each trajectory, the increments between
#'   consecutive grid points. They are computed directly by the bridges, not
#'   as differences of the resampled levels: a very small increment keeps its
#'   significant digits instead of being rounded to 0;
#' * `attr(res, "observed")` gives the observed data, one data frame per
#'   trajectory (with the method `"ExtBridge"`, the truncated times when a
#'   truncation was needed);
#' * `attr(res, "model")`, `attr(res, "parameters")` and
#'   `attr(res, "method")` recall the call.
#'
#' Printing `res` shows the model, the parameters, the steps and the first
#' resampled levels; [plot_bridge()] compares the observed and resampled
#' levels. The bridge values are random: use `set.seed()` before a
#' resampling to make it reproducible.
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' Hoyle, E., Hughston, L. P. and Macrina, A. (2011). Lévy random bridges and
#' the modelling of financial information. *Stochastic Processes and their
#' Applications*, 121(4), 856-884.
#'
#' @examples
#' # 1. The bridge: draws of X_2 given X_1 = 3 and X_4 = 6 (gamma, a = 2, b = 3)
#' set.seed(1)
#' gamma_bridge(t_A = 1, x_A = 3, t_B = 4, x_B = 6, t = 2, theta = c(a = 2, b = 3), n = 5)
#'
#' # 2. Resampling of one gamma trajectory observed at 11 irregular times
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3)
#'
#' # on the regular grid 0, 2.25, 4.5, ..., 22.5
#' r <- gamma_resampling(x, theta = c(a = 2, b = 3))
#' r
#' r[[1]]
#' attr(r, "teq")
#' attr(r, "dx")
#' plot_bridge(r)
#'
#' # 3. Three trajectories, each with its own inspection times (a list)
#' plan <- DegradSim::inspection_times("decreasing", n = 10, n_traj = 3,
#'                                     tmax = c(20, 30, 40), digits = 1)
#' xs <- DegradSim::sim_gamma(plan, a = 2, b = 3)
#'
#' # One grid per trajectory, then the two methods with a common step
#' theta <- c(a = 2, b = 3)
#' attr(gamma_resampling(xs, theta), "teq")
#' attr(gamma_resampling(xs, theta, method = "ExtBridge", common_step = TRUE), "teq")
#' attr(gamma_resampling(xs, theta, method = "MultiBridge"), "teq")
#' plot_bridge(gamma_resampling(xs, theta, method = "MultiBridge"), trajectory = 1:3)
#' @importFrom DegradSim inspection_times
#' @keywords package
"_PACKAGE"
