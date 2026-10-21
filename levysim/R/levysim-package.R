#' levysim: simulation of degradation processes
#'
#' `levysim` simulates exactly one or several degradation trajectories of six
#' stochastic processes, homogeneous (Lévy processes) or not, at inspection
#' times that can be irregular and differ from one trajectory to another. The
#' simulated trajectories can be plotted and extracted as a data frame.
#'
#' @section Functions of the package:
#' Simulation (click on a function name for its help page):
#'
#' | Function | Process | Parameters | Functions of time |
#' |---|---|---|---|
#' | [sim_gamma()] | gamma | shape `a`, rate `b` | `A` |
#' | [sim_wiener()] | Wiener with drift | `mu`, `sigma` | `Lambda`, `tau` |
#' | [sim_ig()] | inverse Gaussian | `mu`, `eta` | `Lambda` |
#' | [sim_vg()] | variance-gamma | `theta`, `sigma`, `nu` | `Lambda` |
#' | [sim_meixner()] | Meixner | `a`, `b`, `d`, `m` | `Lambda` |
#' | [sim_igam()] | inverse gamma | `theta`, `alpha`, `beta` | `eta_inv` |
#'
#' The functions of time are the identity by default: the processes are then
#' homogeneous, that is Lévy processes (except the inverse gamma process, which
#' is never a Lévy process).
#'
#' Simulation methods:
#' * gamma: [stats::rgamma()]; Wiener: [stats::rnorm()];
#' * inverse Gaussian: method of Michael, Schucany and Haas;
#' * variance-gamma: subordination, or difference of two gamma processes
#'   (argument `method`);
#' * Meixner: numerical inversion of the distribution function (`Runuran`);
#' * inverse gamma: exact passage algorithm.
#'
#' Results and inspection schedules:
#' * [plot_traj()]: plot of the simulated trajectories (`ggplot2`);
#' * [plot_plan()]: plot of the inspection times, one row per trajectory
#'   (`ggplot2`);
#' * [traj_df()]: simulated data as a data frame with columns `id`, `t`, `X`,
#'   and optionally the time steps `dt` and increments `dX`;
#' * [inspection_times()]: inspection times with a regular, increasing,
#'   decreasing or alternating number of inspections over time.
#'
#' @section Getting started:
#' A simulation always follows the same three steps.
#' 1. **Choose the inspection times** `t`. Three forms are accepted by every
#'    `sim_` function:
#'    * a vector, for instance `seq(0, 10, by = 0.5)`: all the trajectories
#'      share these times, and their number is given by `n_traj`;
#'    * a list of vectors, one per trajectory, for instance
#'      `list(c(0, 1, 3), c(0, 2, 4, 6))`;
#'    * the result of [inspection_times()], for instance with the pattern
#'      `"increasing"` for trajectories inspected more and more often.
#'
#'    The times of each trajectory start at 0 and are strictly increasing.
#' 2. **Simulate** with the function of the chosen process and its
#'    parameters, for instance [sim_gamma()] for a gamma process.
#'    For a non-homogeneous process, give the function of time, for instance
#'    `A = function(t) t^1.5`.
#' 3. **Store and use the result** `x`. Its form depends on the times:
#'    * same times for all the trajectories: a matrix with the column `t`
#'      and one column per trajectory, `X1`, ..., `Xn`; the \eqn{k}-th
#'      trajectory is `x[, "Xk"]`;
#'    * different times: a list with one data frame (`t`, `X`) per
#'      trajectory; the \eqn{k}-th trajectory is `x[[k]]`.
#'
#'    In both cases, `x` prints the model, the parameters and the first
#'    simulated values, `traj_df(x)` gives all the data in one data frame
#'    (`id`, `t`, `X`, and `dt`, `dX` with `increments = TRUE`) and
#'    `plot_traj(x)` plots the trajectories.
#'
#' Use `set.seed()` before a simulation to make it reproducible.
#'
#' @section Exact simulation:
#' For the first five processes, which have independent increments, the joint
#' distribution at the inspection times is determined by the distributions of
#' the increments: the simulation draws independent increments with the
#' prescribed distributions and cumulates them, without any time
#' discretization. The inverse gamma process, whose increments are dependent,
#' uses a specific exact algorithm.
#'
#' @examples
#' # 1. Inspection times: 10 trajectories, inspected more and more often
#' plan <- inspection_times("increasing", n = 20, n_traj = 10, tmax = 100)
#' plot_plan(plan)
#'
#' # 2. Simulation of a homogeneous gamma process
#' set.seed(1)
#' x <- sim_gamma(plan, a = 2, b = 3)
#'
#' # 3. Use of the result: the same times for all the trajectories give a
#' #    matrix (column t, then X1, ..., X10)
#' x
#' x[, "X1"]
#' head(traj_df(x, increments = TRUE))
#' plot_traj(x)
#' @keywords package
"_PACKAGE"

utils::globalVariables(c("t", "X", "id"))
