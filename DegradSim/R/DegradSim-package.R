#' DegradSim: simulation of Levy and time-scaled Levy degradation processes
#'
#' `DegradSim` simulates exactly, at given inspection times, one or several
#' degradation trajectories of two families of processes:
#' 1. **Levy processes**: gamma, Wiener, inverse Gaussian, variance-gamma,
#'    Meixner, normal inverse Gaussian and compound Poisson processes;
#' 2. **time-scaled Levy processes** \eqn{X_t = Y_{A(t)}}, where \eqn{Y} is a
#'    Levy process and \eqn{A} a time-scaling function: time-scaled gamma,
#'    Wiener and inverse Gaussian processes.
#'
#' The inspection times can be irregular and differ from one trajectory to
#' another; [inspection_times()] builds them with a chosen pattern. The
#' simulated trajectories can be printed, plotted and extracted as a data
#' frame.
#'
#' @section Levy processes:
#' **Definition** (Sato, 2001, p. 3). A real stochastic process
#' \eqn{(X_t)_{t \ge 0}} is a Levy process if
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
#' Levy process, the distribution of \eqn{X_t} is infinitely divisible for
#' every \eqn{t}; conversely, every infinitely divisible distribution is the
#' distribution of \eqn{X_1} for a Levy process, unique in distribution (Sato, 2001,
#' Theorems 1.1 and 1.2).
#'
#' The simulation functions, with the distribution of an increment over a time
#' \eqn{\Delta t = t - s}:
#'
#' | Function | Process | Parameters | \eqn{X_t - X_s} |
#' |---|---|---|---|
#' | [sim_gamma()] | gamma | `a`, `b` | \eqn{\Gamma(a \Delta t, b)} |
#' | [sim_wiener()] | Wiener | `mu`, `sigma` | \eqn{\mathcal{N}(\mu \Delta t, \sigma^2 \Delta t)} |
#' | [sim_ig()] | inverse Gaussian | `mu`, `eta` | \eqn{IG(\mu \Delta t, \eta \Delta t^2)} |
#' | [sim_vg()] | variance-gamma | `theta`, `sigma`, `nu` | \eqn{\theta G_{\Delta t} + \sigma B_{G_{\Delta t}}} |
#' | [sim_meixner()] | Meixner | `a`, `b`, `d`, `m` | \eqn{Meixner(a, b, d \Delta t, m \Delta t)} |
#' | [sim_nig()] | NIG | `alpha`, `beta`, `delta`, `mu` | \eqn{NIG(\alpha, \beta, \delta \Delta t, \mu \Delta t)} |
#' | [sim_cpoisson()] | compound Poisson | `lambda`, `rjump` | sum of \eqn{Poisson(\lambda \Delta t)} jumps |
#'
#' For the variance-gamma process, \eqn{B} is a standard Brownian motion and
#' \eqn{G} an independent gamma process with
#' \eqn{G_{\Delta t} \sim \Gamma(\Delta t / \nu, 1 / \nu)} (see [sim_vg()]).
#'
#' @section Other Levy processes:
#' There are other Levy processes. `DegradSim` does not simulate them, except
#' the Poisson process, given by [sim_cpoisson()] with jumps equal to 1:
#' * the **Poisson process**, increments \eqn{Poisson(\lambda \Delta t)}
#'   (Hoyle, 2010, section 2.11);
#' * the **Cauchy process**, increments \eqn{Cauchy(0, \gamma \Delta t)}, of
#'   infinite mean (Hoyle, 2010, section 2.9);
#' * the **\eqn{\alpha}-stable processes**, of infinite variance for
#'   \eqn{\alpha < 2}; for \eqn{\alpha < 1}, the stable subordinator is
#'   increasing (Samorodnitsky and Taqqu, 1994; exact simulation by Chambers,
#'   Mallows and Stuck, 1976);
#' * the **generalized hyperbolic (GH) processes**, whose distribution at time
#'   1 is generalized hyperbolic (Barndorff-Nielsen, 1977); the NIG and the
#'   variance-gamma processes are special or limiting cases;
#' * the **tempered stable processes**, stable processes whose large jumps are
#'   tempered by an exponential factor (Koponen, 1995; Rosiński, 2007); the
#'   gamma and inverse Gaussian processes are tempered stable subordinators;
#' * the **Levy jump-diffusion processes**, a Brownian motion with drift plus
#'   a compound Poisson process, with normal jumps (Merton, 1976) or double
#'   exponential jumps (Kou, 2002);
#' * the **CGMY process** (Carr, Geman, Madan and Yor, 2002), a tempered
#'   stable process that contains the variance-gamma process as a special
#'   case.
#'
#' @section Time-scaled Levy processes:
#' **Time-scaling function.** Let \eqn{A : [0, +\infty[ \to [0, +\infty[} be
#' a non-decreasing function, right-continuous with left limits, with
#' \eqn{A(0) = 0} and \eqn{A(t) \to +\infty} as \eqn{t \to +\infty} (a
#' time-scaling function, Cha and Mercier, 2019, Example 1), and let \eqn{Y}
#' be a Levy process. The process
#' \deqn{X_t = Y_{A(t)}, \qquad t \ge 0,}
#' is a time-scaled Levy process (Kahle et al., 2016, p. xix-xx; Cha and
#' Mercier, 2019, section 2): \eqn{X_0 = 0}, its increments are independent
#' and \eqn{X_t - X_s} has the distribution of \eqn{Y_{A(t) - A(s)}}. Its
#' increments are not stationary when \eqn{A} is not linear: \eqn{X} is not a
#' Levy process, but a non-homogeneous Levy process. When \eqn{A} is
#' continuous, \eqn{X} is also stochastically continuous: it is an additive
#' process in the sense of Sato (2001, p. 3), that is a process satisfying
#' the conditions 1, 2, 4 and 5 above.
#'
#' **The time-scaling function in the package.** The function \eqn{A} is
#' given as the argument `A` of the functions [sim_gamma()], [sim_wiener()]
#' and [sim_ig()]. It is a vectorized R function of `t`: `A(t)` returns one
#' value per time, for instance with `function(t) 2 * t^1.5`. Each of these
#' functions checks that
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
#' \eqn{B} a standard Brownian motion.
#'
#' **Linear time-scaling function.** With \eqn{A(t) = a t}, the time-scaled
#' gamma process is the gamma process with shape \eqn{a}. With
#' \eqn{A(t) = \mu t}, the time-scaled Wiener process
#' \eqn{\mu t + \sigma B_{\mu t}} has the distribution of
#' \eqn{\mu t + \sigma \sqrt{\mu}\, B_t}, since \eqn{(B_{\mu t})_{t \ge 0}}
#' has the distribution of \eqn{(\sqrt{\mu}\, B_t)_{t \ge 0}}: it is the
#' Wiener process with drift \eqn{\mu} and diffusion coefficient
#' \eqn{\sigma \sqrt{\mu}}. In the same way, with \eqn{A(t) = \mu t}, the
#' increments of the time-scaled inverse Gaussian process are
#' \eqn{IG(\mu \Delta t, \eta \mu^2 \Delta t^2)}: it is the inverse Gaussian
#' process with parameters \eqn{\mu} and \eqn{\eta \mu^2}.
#'
#' @section Getting started:
#' A simulation always follows the same three steps.
#' 1. **Choose the inspection times** `t`. Three forms are accepted by every
#'    `sim_` function:
#'    * a numeric vector, for instance the component `time` of the result of
#'      [inspection_times()]: all the trajectories share these times, and
#'      their number is given by `n_traj`;
#'    * a list of numeric vectors, one per trajectory, for instance
#'      `list(c(0, 1, 3), c(0, 2, 4, 6))`;
#'    * the data frame returned by [inspection_times()], which gives one plan
#'      per trajectory, for instance with the pattern `"decreasing"` for
#'      trajectories inspected more and more often (`n_traj` is then ignored).
#'
#'    The times of each trajectory must be finite, start at 0 and be
#'    strictly increasing, with at least two times (0 and one inspection).
#'    These conditions are checked by every `sim_` function, also for the
#'    result of [inspection_times()], and an error explains which condition
#'    fails.
#' 2. **Simulate** with the function of the chosen process and its
#'    parameters; for a gamma process, for instance,
#'    `sim_gamma(t, a = 2, b = 3)`. For a time-scaled Levy process, give the time-scaling function
#'    `A` (see the section Time-scaled Levy processes) instead of the fixed
#'    parameter, for instance `sim_gamma(t, b = 3, A = function(t) 2 * t^1.5)`.
#' 3. **Store and use the result** `x`, an object of class `DegradSim`. Its
#'    form depends on the times:
#'    * **same times for all the trajectories**: a matrix with the column `t`
#'      (inspection times, starting at 0) and one column per trajectory, `X1`,
#'      ..., `Xn` (levels, equal to 0 at `t = 0`). The levels of the
#'      \eqn{k}-th trajectory are the column `x[, "Xk"]`, which is also the
#'      column number \eqn{k + 1}, `x[, k + 1]`;
#'    * **different times**: a list with one data frame per trajectory. The
#'      data frame `x[[k]]` of the \eqn{k}-th trajectory has two columns: `t`
#'      (its inspection times) and `X` (its levels). A column can be taken by
#'      its name or by its position: `x[[k]]$t` is the same as `x[[k]][, 1]`
#'      (first column), and `x[[k]]$X` is the same as `x[[k]][, 2]` (second
#'      column).
#'
#'    In both cases, `x` prints the model, the parameters, the time-scaling
#'    function if any and the first simulated values; `traj_df(x)` gives all
#'    the data in one data frame (`id`, `t`, `X`, and `dt`, `dX` with
#'    `increments = TRUE`) and `plot_traj(x)` plots the trajectories. See the
#'    examples below.
#'
#' Use `set.seed()` before a simulation to make it reproducible.
#'
#' Other functions:
#' * [inspection_times()]: inspection times following a chosen pattern
#'   (regular, with increasing or decreasing gaps, or alternating between frequent and
#'   rare inspections);
#' * [plot_plan()]: plot of the inspection times, one row per trajectory
#'   (`ggplot2`);
#' * [plot_traj()]: plot of the simulated trajectories (`ggplot2`);
#' * [traj_df()]: simulated data as a data frame.
#'
#' @section Simulation methods:
#' The argument `method` of the `sim_` functions chooses among the three
#' simulation methods of Kahle et al. (2016) for the Wiener process (section
#' 1.3) and the gamma process (section 2.5). Only `"RWA"` is currently
#' available.
#'
#' **Random walk approximation** (`method = "RWA"`, default). All the
#' simulated processes have independent increments with known distributions.
#' The values at the inspection times \eqn{0 = t_0 < t_1 < \dots < t_n} are
#' therefore simulated by the random walk approximation (Kahle et al., 2016,
#' section 1.3.1, algorithm 1.7, p. 6, for the Brownian motion, and section
#' 2.5.1, algorithm 2.32, p. 81, for the gamma process):
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
#' **Bridge sampling** (`method = "Bridge"`; Kahle et al., 2016, section
#' 1.3.2, algorithm 1.9, p. 8, for the Brownian motion, and section 2.5.2,
#' algorithm 2.34, p. 85, for the gamma process). The value at
#' the end \eqn{T} is drawn first, then the value at the middle of
#' \eqn{[0, T]} given the two end values, by the bridge of the process, and
#' so on, each interval being divided in half at each step: the values at the
#' times \eqn{i T / 2^n} have exactly the distribution of the process.
#'
#' **Truncated series representation** (`method = "Series"`, argument
#' `type`). The path is approximated on the whole interval
#' \eqn{[0, T]} by a truncated series: the Karhunen-Loève expansion of the
#' Brownian motion (Kahle et al., 2016, section 1.3.3, p. 9, equation \[1.1\]),
#' and, for the gamma process, compound Poisson process approximations
#' obtained by truncating the series representations of Bondesson, of
#' rejection and of thinning (section 2.5.3, algorithm 2.35, p. 87, based on
#' theorem 2.31).
#'
#' @references
#' **Levy processes and time-scaled Levy processes**
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
#' Normal inverse Gaussian process: Barndorff-Nielsen, O. E. (1997). Normal
#' inverse Gaussian distributions and stochastic volatility modelling.
#' *Scandinavian Journal of Statistics*, 24(1), 1-13.
#' \doi{10.1111/1467-9469.00045}
#'
#' Compound Poisson process: Kahle et al. (2016), chapter 3, above.
#'
#' Hoyle, A. E. V. (2010). *Information-Based Models for Finance and
#' Insurance*. PhD thesis, Imperial College London. arXiv:1010.0829.
#'
#' **Other Levy processes**
#'
#' Barndorff-Nielsen, O. E. (1977). Exponentially decreasing distributions for
#' the logarithm of particle size. *Proceedings of the Royal Society of London.
#' Series A*, 353, 401-419.
#'
#' Carr, P., Geman, H., Madan, D. B. and Yor, M. (2002). The fine structure of
#' asset returns: an empirical investigation. *The Journal of Business*,
#' 75(2), 305-332. \doi{10.1086/338705}
#'
#' Chambers, J. M., Mallows, C. L. and Stuck, B. W. (1976). A method for
#' simulating stable random variables. *Journal of the American Statistical
#' Association*, 71(354), 340-344.
#'
#' Koponen, I. (1995). Analytic approach to the problem of convergence of
#' truncated Lévy flights towards the Gaussian stochastic process. *Physical
#' Review E*, 52, 1197-1199.
#'
#' Kou, S. G. (2002). A jump-diffusion model for option pricing. *Management
#' Science*, 48(8), 1086-1101. \doi{10.1287/mnsc.48.8.1086.166}
#'
#' Merton, R. C. (1976). Option pricing when underlying stock returns are
#' discontinuous. *Journal of Financial Economics*, 3(1-2), 125-144.
#' \doi{10.1016/0304-405X(76)90022-2}
#'
#' Rosiński, J. (2007). Tempering stable processes. *Stochastic Processes and
#' their Applications*, 117(6), 677-707.
#'
#' Samorodnitsky, G. and Taqqu, M. S. (1994). *Stable Non-Gaussian Random
#' Processes*. Chapman and Hall, New York.
#'
#' **Simulation methods**
#'
#' Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
#' variates using transformations with multiple roots. *The American
#' Statistician*, 30(2), 88-90.\cr
#' \doi{10.1080/00031305.1976.10479147}
#'
#' Hörmann, W. and Leydold, J. (2003). Continuous random variate generation by
#' fast numerical inversion. *ACM Transactions on Modeling and Computer
#' Simulation*, 13(4), 347-362. \doi{10.1145/945511.945517}
#'
#' @examples
#' # 1. Inspection times: 10 trajectories, inspected more and more often
#' plan <- inspection_times("decreasing", n = 20, n_traj = 10, tmax = 100)
#' plot_plan(plan)
#'
#' # 2. Simulation of a gamma process (a Levy process)
#' set.seed(1)
#' x <- sim_gamma(plan, a = 0.2, b = 1)
#'
#' # 3. Same times for all the trajectories: a matrix (column t, then X1, ..., X10)
#' x
#' x[, "t"]      # inspection times
#' x[, "X1"]     # levels of trajectory 1
#' x[, 2]        # the same column, by its position
#' head(traj_df(x, increments = TRUE))
#' plot_traj(x)
#'
#' # 4. Time-scaled gamma process with A(t) = 0.02 t^1.5
#' y <- sim_gamma(plan, b = 1, A = function(t) 0.02 * t^1.5)
#' y
#' plot_traj(y, "Time-scaled gamma process")
#'
#' # 5. Different inspection times: a list of data frames (t, X)
#' z <- sim_gamma(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8)), a = 2, b = 3)
#' z[[2]]        # trajectory 2: data frame with columns t and X
#' z[[2]]$t      # its inspection times, the same as z[[2]][, 1]
#' z[[2]]$X      # its levels, the same as z[[2]][, 2]
#' z[[2]][, 2]
#' @keywords package
"_PACKAGE"

utils::globalVariables(c("t", "X", "id"))
