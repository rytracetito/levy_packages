#' DegradSim: simulation of Lévy and non-homogeneous Lévy degradation processes
#'
#' `DegradSim` simulates exactly, at given inspection times, one or several
#' degradation trajectories of two families of processes:
#' 1. **Lévy processes** (homogeneous case): gamma, Wiener, inverse Gaussian,
#'    variance-gamma and Meixner processes;
#' 2. **non-homogeneous Lévy processes** obtained by a time change
#'    \eqn{X_t = Y_{A(t)}} of a Lévy process \eqn{Y}: gamma, Wiener and
#'    inverse Gaussian processes.
#'
#' The inspection times can be irregular and differ from one trajectory to
#' another; [inspection_times()] builds them with a chosen pattern. The
#' simulated trajectories can be printed, plotted and extracted as a data
#' frame.
#'
#' @section Lévy processes:
#' **Definition** (Sato, 2001, p. 3). A real stochastic process
#' \eqn{(X_t)_{t \ge 0}} is a Lévy process if
#' 1. it has independent increments: for all \eqn{n \ge 1} and
#'    \eqn{0 \le t_0 < t_1 < \dots < t_n}, the random variables
#'    \eqn{X_{t_0}, X_{t_1} - X_{t_0}, \dots, X_{t_n} - X_{t_{n-1}}} are
#'    independent;
#' 2. it starts at the origin: \eqn{X_0 = 0} almost surely;
#' 3. it is time homogeneous: the distribution of
#'    \eqn{(X_{s+t} - X_s)_{t \ge 0}} does not depend on \eqn{s};
#' 4. it is stochastically continuous: for all \eqn{\varepsilon > 0},
#'    \eqn{P(|X_{s+t} - X_s| > \varepsilon) \to 0} as \eqn{t \to 0};
#' 5. its trajectories \eqn{t \mapsto X_t} are right-continuous with left
#'    limits (càdlàg), almost surely.
#'
#' **Infinitely divisible distribution** (Sato, 2001, p. 4). A distribution
#' \eqn{\mu} on \eqn{\mathbb{R}} is infinitely divisible if, for every integer
#' \eqn{n \ge 1}, there is a distribution \eqn{\mu_n} such that
#' \eqn{\mu = \mu_n^{n*}}, the \eqn{n}-fold convolution of \eqn{\mu_n}: a
#' random variable with distribution \eqn{\mu} has the distribution of the sum
#' of \eqn{n} independent random variables with distribution \eqn{\mu_n}. For a
#' Lévy process, the distribution of \eqn{X_t} is infinitely divisible for
#' every \eqn{t}; conversely, every infinitely divisible distribution is the
#' distribution of \eqn{X_1} for a Lévy process, unique in law (Sato, 2001,
#' Theorems 1.1 and 1.2). The gamma, normal, inverse Gaussian, variance-gamma
#' and Meixner distributions are infinitely divisible.
#'
#' The simulation functions, with the distribution of an increment over a time
#' \eqn{\Delta t = t - s}:
#'
#' | Function | Process | Parameters | \eqn{X_t - X_s} |
#' |---|---|---|---|
#' | [sim_gamma()] | gamma | `a`, `b` | \eqn{\Gamma(a \Delta t, b)}, rate \eqn{b} |
#' | [sim_wiener()] | Wiener with drift | `mu`, `sigma` | \eqn{\mathcal{N}(\mu \Delta t, \sigma^2 \Delta t)} |
#' | [sim_ig()] | inverse Gaussian | `mu`, `eta` | \eqn{IG(\mu \Delta t, \eta \Delta t^2)} |
#' | [sim_vg()] | variance-gamma | `theta`, `sigma`, `nu` | \eqn{VG(\theta, \sigma, \nu)} over \eqn{\Delta t} |
#' | [sim_meixner()] | Meixner | `a`, `b`, `d`, `m` | \eqn{Meixner(a, b, d \Delta t, m \Delta t)} |
#'
#' @section Non-homogeneous Lévy processes:
#' **Time change.** Let \eqn{A : [0, +\infty[ \to [0, +\infty[} be a
#' non-decreasing function, right-continuous with left limits, with
#' \eqn{A(0) = 0} and \eqn{A(t) \to +\infty} as \eqn{t \to +\infty}, and let
#' \eqn{Y} be a homogeneous Lévy process. The process
#' \deqn{X_t = Y_{A(t)}, \qquad t \ge 0,}
#' is a non-homogeneous Lévy process (Kahle et al., 2016, p. xix-xx; Cha and
#' Mercier, 2019, section 2): \eqn{X_0 = 0}, its increments are independent
#' and \eqn{X_t - X_s} has the distribution of \eqn{Y_{A(t) - A(s)}}. Its
#' increments are not stationary when \eqn{A} is not linear. When \eqn{A} is
#' continuous, \eqn{X} is also stochastically continuous: it is an additive
#' process in the sense of Sato (2001, p. 3), that is a process satisfying
#' the conditions 1, 2, 4 and 5 above.
#'
#' **The time change in the package.** \eqn{A} is given as the argument `A`
#' of [sim_gamma()], [sim_wiener()] and [sim_ig()], a vectorized R function
#' of `t` (`A(t)` returns one value per time), for instance
#' `function(t) 2 * t^1.5`. Each of these functions checks that
#' \eqn{A(0) = 0} and that \eqn{A} is finite and non-decreasing at the
#' inspection times, and raises an error otherwise. The conditions
#' "right-continuous with left limits" and \eqn{A(t) \to +\infty} cannot be
#' checked on a finite set of times: they are left to the user. If \eqn{A} is
#' constant between two inspection times, the increment of \eqn{X} is 0.
#'
#' **Identifiability.** One parameter of \eqn{Y} (\eqn{a} for the gamma
#' process, \eqn{\mu} for the Wiener and inverse Gaussian processes) is set to
#' 1 and contained in \eqn{A}. With \eqn{\Delta A = A(t) - A(s)}:
#'
#' | Function | Parameters | \eqn{X_t - X_s} |
#' |---|---|---|
#' | [sim_gamma()] | `b`, `A` | \eqn{\Gamma(\Delta A, b)} |
#' | [sim_wiener()] | `sigma`, `A` | \eqn{\mathcal{N}(\Delta A, \sigma^2 \Delta A)} |
#' | [sim_ig()] | `eta`, `A` | \eqn{IG(\Delta A, \eta \Delta A^2)} |
#'
#' The time-scaled Wiener process is \eqn{X_t = A(t) + \sigma B_{A(t)}}, with
#' \eqn{B} a standard Brownian motion. The homogeneous case is
#' \eqn{A(t) = a t} for the gamma process, and \eqn{A(t) = \mu t} for the
#' Wiener and inverse Gaussian processes, with \eqn{\sigma / \sqrt{\mu}}
#' instead of \eqn{\sigma} and \eqn{\eta / \mu^2} instead of \eqn{\eta}.
#'
#' @section Getting started:
#' A simulation always follows the same three steps.
#' 1. **Choose the inspection times** `t`. Three forms are accepted by every
#'    `sim_` function:
#'    * a numeric vector, for instance `seq(0, 10, by = 0.5)`: all the
#'      trajectories share these times, and their number is given by
#'      `n_traj`;
#'    * a list of numeric vectors, one per trajectory, for instance
#'      `list(c(0, 1, 3), c(0, 2, 4, 6))`;
#'    * the result of [inspection_times()], for instance with the pattern
#'      `"increasing"` for trajectories inspected more and more often.
#'
#'    The times of each trajectory must be finite, start at 0 and be
#'    strictly increasing, with at least two times (0 and one inspection).
#'    These conditions are checked by every `sim_` function, also for the
#'    result of [inspection_times()], and an error explains which condition
#'    fails.
#' 2. **Simulate** with the function of the chosen process and its
#'    parameters, for instance `sim_gamma(t, a = 2, b = 3)` for a homogeneous
#'    gamma process. For a non-homogeneous process, give the time change `A`
#'    (see the section Non-homogeneous Lévy processes) instead of the fixed
#'    parameter, for instance `sim_gamma(t, b = 3, A = function(t) 2 * t^1.5)`.
#' 3. **Store and use the result** `x`, an object of class `DegradSim`. Its
#'    form depends on the times:
#'    * same times for all the trajectories: a matrix with the column `t`
#'      and one column per trajectory, `X1`, ..., `Xn`; the levels of the
#'      \eqn{k}-th trajectory are `x[, k + 1]` or `x[, "Xk"]`;
#'    * different times: a list with one data frame (`t`, `X`) per
#'      trajectory; the \eqn{k}-th trajectory is `x[[k]]`, with times
#'      `x[[k]]$t` (or `x[[k]][, 1]`) and levels `x[[k]]$X` (or
#'      `x[[k]][, 2]`).
#'
#'    In both cases, `x` prints the model, the parameters, the time change if
#'    any and the first simulated values; `traj_df(x)` gives all the data in
#'    one data frame (`id`, `t`, `X`, and `dt`, `dX` with
#'    `increments = TRUE`) and `plot_traj(x)` plots the trajectories.
#'
#' Use `set.seed()` before a simulation to make it reproducible.
#'
#' Other functions:
#' * [inspection_times()]: inspection times following a chosen pattern
#'   (regular, increasing, decreasing, or alternating between frequent and
#'   rare inspections);
#' * [plot_plan()]: plot of the inspection times, one row per trajectory
#'   (`ggplot2`);
#' * [plot_traj()]: plot of the simulated trajectories (`ggplot2`);
#' * [traj_df()]: simulated data as a data frame.
#'
#' @section Exact simulation at the inspection times (random walk approximation):
#' All the simulated processes have independent increments with known
#' distributions. The values at the inspection times
#' \eqn{0 = t_0 < t_1 < \dots < t_n} are therefore simulated by the random walk
#' approximation (RWA) of Kahle et al. (2016, algorithms 1.7 and 2.32):
#' independent increments \eqn{\Delta X_i} with the distribution of
#' \eqn{X_{t_i} - X_{t_{i-1}}} are drawn, then \eqn{X_{t_0} = 0} and
#' \eqn{X_{t_i} = \Delta X_1 + \dots + \Delta X_i}. The vector
#' \eqn{(X_{t_1}, \dots, X_{t_n})} has exactly the distribution of the process
#' at the inspection times; the name "approximation" refers only to the path
#' between two inspection times, which is not simulated. Algorithm 2.32 is
#' written for the gamma process, but the same scheme applies to every process
#' with independent increments. The way each increment is drawn is given in
#' the help page of each simulation function.
#'
#' @references
#' **Lévy processes and time change**
#'
#' Sato, K. (2001). Basic results on Lévy processes. In O. E.
#' Barndorff-Nielsen, T. Mikosch and S. I. Resnick (eds.), *Lévy Processes:
#' Theory and Applications*, 3-37. Birkhäuser, Boston.
#' \doi{10.1007/978-1-4612-0197-7}
#'
#' Kahle, W., Mercier, S. and Paroissin, C. (2016). *Degradation Processes in
#' Reliability*. ISTE and Wiley, London. ISBN 978-1-84821-888-8.
#'
#' Cha, J. H. and Mercier, S. (2019). Transformed Lévy processes as
#' state-dependent wear models. *Advances in Applied Probability*, 51(2),
#' 468-486. \doi{10.1017/apr.2019.21}
#'
#' **Processes**
#'
#' Gamma process: van Noortwijk, J. M. (2009). A survey of the application of
#' gamma processes in maintenance. *Reliability Engineering & System Safety*,
#' 94(1), 2-21. \doi{10.1016/j.ress.2007.03.019}
#'
#' Wiener process: Whitmore, G. A. (1995). Estimating degradation by a Wiener
#' diffusion process subject to measurement error. *Lifetime Data Analysis*,
#' 1, 307-319.
#'
#' Inverse Gaussian process: Wang, X. and Xu, D. (2010). An inverse Gaussian
#' process model for degradation data. *Technometrics*, 52(2), 188-197.
#' \doi{10.1198/TECH.2009.08197}
#'
#' Variance-gamma process: Madan, D. B., Carr, P. P. and Chang, E. C. (1998).
#' The variance gamma process and option pricing. *European Finance Review*,
#' 2(1), 79-105. \doi{10.1023/A:1009703431535}
#'
#' Meixner process: Schoutens, W. and Teugels, J. L. (1998). Lévy processes,
#' polynomials and martingales. *Communications in Statistics. Stochastic
#' Models*, 14(1-2), 335-349.
#'
#' **Simulation methods**
#'
#' Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
#' variates using transformations with multiple roots. *The American
#' Statistician*, 30(2), 88-90. \doi{10.1080/00031305.1976.10479147}
#'
#' Hörmann, W. and Leydold, J. (2003). Continuous random variate generation by
#' fast numerical inversion. *ACM Transactions on Modeling and Computer
#' Simulation*, 13(4), 347-362. \doi{10.1145/945511.945517}
#'
#' @examples
#' # 1. Inspection times: 10 trajectories, inspected more and more often
#' plan <- inspection_times("increasing", n = 20, n_traj = 10, tmax = 100)
#' plot_plan(plan)
#'
#' # 2. Simulation of a homogeneous gamma process (a Lévy process)
#' set.seed(1)
#' x <- sim_gamma(plan, a = 0.2, b = 1)
#'
#' # 3. Use of the result: the same times for all the trajectories give a
#' #    matrix (column t, then X1, ..., X10)
#' x
#' x[, "X1"]
#' head(traj_df(x, increments = TRUE))
#' plot_traj(x)
#'
#' # 4. Non-homogeneous gamma process with the time change A(t) = 0.02 t^1.5
#' y <- sim_gamma(plan, b = 1, A = function(t) 0.02 * t^1.5)
#' y
#' plot_traj(y, "Non-homogeneous gamma process")
#'
#' # 5. Different inspection times: a list of data frames (t, X)
#' z <- sim_wiener(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), mu = 0.5, sigma = 1.5)
#' z[[2]]
#' z[[2]]$X
#' @keywords package
"_PACKAGE"

utils::globalVariables(c("t", "X", "id"))
