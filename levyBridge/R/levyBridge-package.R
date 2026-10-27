#' levyBridge: stochastic bridges for Levy degradation processes
#'
#' `levyBridge` resamples degradation trajectories observed at irregular
#' inspection times onto a regular time grid, by exact simulation of the
#' bridge of the process between the known values. Under the model, the
#' reconstructed increments are independent and identically distributed: they
#' are the input of the goodness-of-fit tests for Levy degradation processes.
#'
#' @section Bridge of a Levy process:
#' Let \eqn{X} be a Levy process such that \eqn{X_s} has the density
#' \eqn{f_s}, and \eqn{t_A < t < t_B}. Given \eqn{X_{t_A} = x_A} and
#' \eqn{X_{t_B} = x_B}, \eqn{X_t} has the density
#' \deqn{f_{t | t_A, t_B}(x) = \frac{f_{t - t_A}(x - x_A)\, f_{t_B - t}(x_B - x)}{f_{t_B - t_A}(x_B - x_A)}}
#' (Hoyle, Hughston and Macrina, 2011). Each `_bridge` function draws exactly
#' from this law, for one process:
#'
#' | Function | Process | Parameters | Simulation of the bridge |
#' |---|---|---|---|
#' | [gamma_bridge()] | gamma | shape `a` | beta variable |
#' | [wiener_bridge()] | Wiener | `sigma` | Gaussian variable |
#' | [ig_bridge()] | inverse Gaussian | `eta` | IG or scaled inverse IG |
#' | [vg_bridge()] | variance-gamma | `theta`, `sigma`, `nu` | GIG, beta, Gaussian |
#' | [meixner_bridge()] | Meixner | `a`, `d`, `m` | rejection, Cauchy envelope |
#'
#' The law of a bridge never depends on all the parameters of the process: the
#' gamma bridge does not depend on the rate `b`, the Brownian bridge does not
#' depend on the drift `mu`, the inverse Gaussian bridge does not depend on
#' `mu`, and the Meixner bridge does not depend on `b`. Click on a function
#' name for the law of its bridge and the details of the simulation.
#'
#' @section Resampling methods:
#' Every `_bridge` function has an argument `method`:
#' * `"Bridge"` (default): interpolation algorithm of Grall-Maës (2012); see
#'   also the help page of [gamma_bridge()]. A trajectory observed at
#'   \eqn{0 = t_0 < t_1 < \dots < t_N} is resampled on the regular grid
#'   \eqn{t'_k = k\, t_N / N}, \eqn{k = 0, \dots, N}: as many points as
#'   observations, with the time step \eqn{\Delta t_{eq} = t_N / N}. The grid
#'   is covered from left to right; each grid point is drawn from the bridge
#'   between the last known value before it (an observation or the previous
#'   grid point, already drawn) and the next observation. A grid point equal
#'   to an inspection time keeps the observed value;
#' * `"ExtBridge"`: extension of this method whose grid contains all the
#'   inspection times. With \eqn{d} the number of decimal places of the times
#'   (argument `digits`, found automatically by default), the time step is
#'   \eqn{h = 10^{-d}\, \mathrm{gcd}(10^d t_1, \dots, 10^d t_N)}, the largest
#'   one such that every inspection time is a grid point. The grid is covered
#'   as for `"Bridge"`, so the pseudo-sample contains all the observations.
#'   The grid has \eqn{t_N / h \ge N} steps, with equality for regular times.
#'   With several trajectories, `common_step = TRUE` gives them the same step
#'   \eqn{h} (gcd of all the gaps of all the trajectories);
#' * `"MultiBridge"`: the method `"Bridge"` for several trajectories, on one
#'   regular grid common to all of them, of step
#'   \eqn{h = (t_{1, max} + \dots + t_{n, max}) / (m_1 + \dots + m_n)} (total
#'   observed duration divided by the total number of increments; \eqn{t_N / N}
#'   for one trajectory). With `option = 1` (default) each trajectory is
#'   resampled up to its own last inspection; with `option = 2` all the
#'   trajectories go up to the last inspection time over all of them, the
#'   shorter ones being extended by increments of the process (this needs all
#'   the parameters, and these extended values do not come from the data).
#'
#' The three methods use the same bridges and the same algorithm: only the
#' grid differs.
#'
#' **Which method?**
#' * one trajectory: `"Bridge"`, or `"ExtBridge"` to keep all the observations
#'   on the grid;
#' * several trajectories, one test per trajectory: `"Bridge"` or
#'   `"ExtBridge"`;
#' * several trajectories, one test on all of them together: the
#'   reconstructed increments must have the same law, hence the same time step:
#'   `"MultiBridge"` (grid of Grall-Maës type, about as many points as
#'   observations) or `"ExtBridge"` with `common_step = TRUE` (grid containing
#'   all the inspection times, often many more points).
#'
#' @section Data:
#' All the `_bridge` functions accept the same formats:
#' * **one trajectory**: a vector of levels `x` and a vector of inspection
#'   times `t`, starting at 0, or a data frame whose first column contains the
#'   times and second column the levels;
#' * **several trajectories with the same inspection times**: the matrix
#'   returned by the `DegradSim` simulation functions, with columns `t`, `X1`,
#'   `X2`, ...;
#' * **several trajectories, each with its own inspection times**: a list of
#'   data frames, one per trajectory (first column the times, second column
#'   the levels), for instance the list returned by `DegradSim` with one plan
#'   per trajectory.
#'
#' With the methods "Bridge" and "ExtBridge", each trajectory gets its own
#' regular grid and its own time step. With the method "MultiBridge", or
#' "ExtBridge" with `common_step = TRUE`, all the trajectories get the same
#' time step: their reconstructed increments then have the same law and can be
#' tested together.
#'
#' @section Results:
#' A `_bridge` function returns an object of class `"levybridge"`:
#' * `res` prints the model, the parameters, the time steps and the first
#'   resampled values;
#' * `res[[k]]` is the \eqn{k}-th trajectory on its regular grid, a data
#'   frame with columns `t` and `X`;
#' * `attr(res, "teq")` gives the time step of the grid of each trajectory;
#' * [increments()] gives the reconstructed increments, one row per increment
#'   with the columns `id` (trajectory), `teq` and `dx`: they are the values
#'   to be tested;
#' * [plot_bridge()] compares the observed and resampled values (`ggplot2`).
#'
#' The bridge values are random: use `set.seed()` before a resampling to make
#' it reproducible.
#'
#' @section Goodness-of-fit test, step by step:
#' The package covers steps 2 to 4 of the methodology of goodness-of-fit tests
#' for a parametric family of Levy processes:
#' 1. **Estimate** the parameters from the observed trajectories, for instance
#'    with the package `levyfit` (`gamma_mle()`, `meixner_mle()`, ...);
#' 2. **Law of the bridge**: it is given in the help page of each `_bridge`
#'    function;
#' 3. **Simulation of the bridge**: done by the `_bridge` function, with the
#'    estimated parameters;
#' 4. **Pseudo-sample on a regular grid**: the result of the `_bridge`
#'    function, with the method "Bridge", "ExtBridge" or "MultiBridge";
#' 5. **Choose a test** for an i.i.d. sample of the law of
#'    \eqn{X_{\Delta t_{eq}}}, for instance Kolmogorov-Smirnov, Cramér-von
#'    Mises or Anderson-Darling;
#' 6. **Apply the test** to the reconstructed increments `increments(res)$dx`.
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
#' Tito, A. R. and Paroissin, C. (in preparation).
#'
#' @examples
#' # 1. One gamma trajectory observed at irregular times
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- c(0, cumsum(rgamma(10, shape = 2 * diff(t), rate = 3)))
#'
#' # Resampling on the regular grid 0, 2.25, 4.5, ..., 22.5
#' r <- gamma_bridge(x, t, a = 2)
#' r
#' r[[1]]
#' attr(r, "teq")
#' increments(r)
#' plot_bridge(r)
#'
#' # 2. The same trajectory with the other processes (parameters of the bridge only)
#' wiener_bridge(x, t, sigma = 1.5)[[1]]
#' ig_bridge(x, t, eta = 3)[[1]]
#' vg_bridge(x, t, theta = 0.5, sigma = 1.5, nu = 2)[[1]]
#' meixner_bridge(x, t, a = 1, d = 1.5, m = 0.2)[[1]]
#'
#' # 3. Several trajectories, each with its own inspection times (a list)
#' liste <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.4, 2.3, 2.9)),
#'               data.frame(t = c(0, 2, 4, 5, 9), X = c(0, 1.6, 2.8, 3.1, 4.2)))
#' r2 <- gamma_bridge(liste, a = 2)
#' attr(r2, "teq")
#' increments(r2)
#' @examplesIf requireNamespace("DegradSim", quietly = TRUE) && requireNamespace("levyfit", quietly = TRUE)
#' # 4. Complete use: simulation (DegradSim), estimation (levyfit), bridge
#' plan <- DegradSim::inspection_times("decreasing", n = 20, n_traj = 10, tmax = 40)
#' sims <- DegradSim::sim_gamma(plan, a = 2, b = 3)
#'
#' # Step 1: estimation of the parameters
#' est <- levyfit::gamma_mle(sims)
#' est
#'
#' # Steps 3 and 4: resampling with the estimated shape
#' r3 <- gamma_bridge(sims, a = est[["a"]])
#' plot_bridge(r3, trajectory = 1:4)
#'
#' # Step 6: the reconstructed increments, to be tested against the law
#' # Gamma(a teq, b) with the estimated parameters
#' d <- increments(r3)
#' head(d)
#' @keywords package
"_PACKAGE"
