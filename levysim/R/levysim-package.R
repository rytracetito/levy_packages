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
#' | Function | Process | Parameters | Functions of time | Method |
#' |---|---|---|---|---|
#' | [sim_gamma()] | gamma | shape `a`, rate `b` | `A` | `rgamma` |
#' | [sim_wiener()] | Wiener with drift | `mu`, `sigma` | `Lambda`, `tau` | `rnorm` |
#' | [sim_ig()] | inverse Gaussian | `mu`, `eta` | `Lambda` | Michael, Schucany and Haas |
#' | [sim_vg()] | variance-gamma | `theta`, `sigma`, `nu` | `Lambda` | subordination, or difference of two gamma processes (`method`) |
#' | [sim_meixner()] | Meixner | `a`, `b`, `d`, `m` | `Lambda` | numerical inversion (`Runuran`) |
#' | [sim_igam()] | inverse gamma | `theta`, `alpha`, `beta` | `eta_inv` | exact passage algorithm |
#'
#' The functions of time are the identity by default: the processes are then
#' homogeneous, that is Lévy processes (except the inverse gamma process, which
#' is never a Lévy process).
#'
#' Results and inspection schedules:
#'
#' | Function | Use |
#' |---|---|
#' | [plot_traj()] | plot of the simulated trajectories (`ggplot2`) |
#' | [traj_df()] | simulated data as a data frame with columns `id`, `t`, `X` |
#' | [inspection_times()] | inspection times with a regular, increasing, decreasing or alternating number of inspections over time |
#'
#' @section Getting started:
#' A simulation always follows the same three steps.
#' 1. **Choose the inspection times** `t`. Three forms are accepted by every
#'    `sim_` function:
#'    * a vector, for instance `seq(0, 10, by = 0.5)`: all the trajectories
#'      share these times, and their number is given by `n_traj`;
#'    * a list of vectors, one per trajectory, for instance
#'      `list(c(0, 1, 3), c(0, 2, 4, 6))`;
#'    * the result of [inspection_times()], for instance
#'      `inspection_times("increasing", n = 20, n_traj = 5, tmax = 10)` for
#'      five trajectories inspected more and more often.
#'
#'    The times of each trajectory start at 0 and are strictly increasing.
#' 2. **Simulate** with the function of the chosen process and its
#'    parameters, for instance `x <- sim_gamma(t, a = 2, b = 3, n_traj = 10)`.
#'    For a non-homogeneous process, give the function of time, for instance
#'    `A = function(t) t^1.5`.
#' 3. **Use the result** `x`:
#'    * `x` prints the model, the parameters and the first simulated values;
#'    * `plot_traj(x)` plots the trajectories;
#'    * `traj_df(x)` gives all the data in one data frame (`id`, `t`, `X`);
#'    * `x[[k]]` is the \eqn{k}-th trajectory, a data frame with columns `t`
#'      and `X`.
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
#'
#' # 2. Simulation of a homogeneous gamma process
#' set.seed(1)
#' x <- sim_gamma(plan, a = 2, b = 3)
#'
#' # 3. Use of the result
#' x
#' plot_traj(x)
#' head(traj_df(x))
#' x[[1]]
#' @keywords package
"_PACKAGE"

utils::globalVariables(c("t", "X", "id"))
