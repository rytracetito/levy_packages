#' levyBridge: stochastic bridges for Lévy degradation processes
#'
#' `levyBridge` resamples degradation trajectories observed at irregular
#' inspection times onto a regular time grid, by exact simulation of the
#' bridge of the process between the known values. Under the model, the
#' reconstructed increments are independent and identically distributed: they
#' are the input of the goodness-of-fit tests for Lévy degradation processes.
#'
#' @section Functions of the package:
#' Resampling (click on a function name for its help page):
#'
#' | Function | Process | Parameters needed | Law of the bridge |
#' |---|---|---|---|
#' | [gamma_bridge()] | gamma | shape `a` | beta |
#' | [wiener_bridge()] | Wiener with drift | `sigma` | Gaussian |
#' | [ig_bridge()] | inverse Gaussian | `eta` | mixture of an IG variable and its inverse |
#' | [vg_bridge()] | variance-gamma | `theta`, `sigma`, `nu` | GIG, beta and Gaussian variables |
#' | [meixner_bridge()] | Meixner | `a`, `d`, `m` | acceptance-rejection |
#'
#' The law of a bridge never depends on all the parameters of the process: for
#' instance the gamma bridge does not depend on the rate `b`, and the Brownian
#' bridge does not depend on the drift `mu`.
#'
#' Results:
#'
#' | Function | Use |
#' |---|---|
#' | [increments()] | reconstructed increments, one row per increment, with the time step `teq` |
#' | [plot_bridge()] | observed values and resampled values on the regular grid (`ggplot2`) |
#'
#' @section Methods:
#' Every `_bridge` function has an argument `method`:
#' * `"Bridge"` (default): interpolation algorithm of Grall-Maës (2012),
#'   described in [gamma_bridge()]. The regular grid has as many points as
#'   observations, with time step \eqn{\Delta t_{eq} = t_N / N};
#' * `"ExtBridge"`: extension of this method that integrates the
#'   observations. Not available yet.
#'
#' @section Getting started:
#' 1. **Prepare the data**: either one trajectory, given by its levels `x`
#'    and its inspection times `t` (starting at 0), or several trajectories,
#'    given as a list of data frames whose first column contains the times and
#'    second column the levels (for instance the result of a `levysim`
#'    simulation function).
#' 2. **Resample** with the bridge of the chosen model and its parameters, for
#'    instance `res <- gamma_bridge(x, t, a = 2)`. The parameters can be
#'    estimated with the package `levyfit`.
#' 3. **Use the result** `res`:
#'    * `res` prints the model, the parameters, the time step and the first
#'      resampled values;
#'    * `res[[k]]` is the \eqn{k}-th trajectory on the regular grid, a data
#'      frame with columns `t` and `X`;
#'    * `increments(res)` gives the reconstructed increments, to be tested;
#'    * `plot_bridge(res)` compares the observed and resampled values.
#'
#' Use `set.seed()` before a resampling to make it reproducible: the bridge
#' values are random.
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
#' # 1. One gamma trajectory observed at irregular times
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- c(0, cumsum(rgamma(10, shape = 2 * diff(t), rate = 3)))
#'
#' # 2. Resampling on the regular grid 0, 2.25, 4.5, ..., 22.5
#' r <- gamma_bridge(x, t, a = 2)
#'
#' # 3. Use of the result
#' r
#' r[[1]]
#' increments(r)
#' plot_bridge(r)
#' @keywords package
"_PACKAGE"
