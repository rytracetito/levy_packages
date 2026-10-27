#' levygof: goodness-of-fit tests for Levy degradation processes
#'
#' `levygof` tests whether degradation trajectories come from a parametric
#' Levy process: gamma (`"HGP"`), Wiener (`"HWP"`), inverse Gaussian
#' (`"InvGaus"`), variance-gamma (`"VG"`) or Meixner (`"Meixner"`), with
#' unknown parameters. The main function is [gof.levy()].
#'
#' @section Methodology:
#' A Levy process has independent and stationary increments: the increments
#' over a same time step are i.i.d. The test proceeds in six steps.
#'
#' 1. **Estimate** the parameters from the observed trajectories by maximum
#'    likelihood (package `levyfit`).
#' 2. **Law of the bridge** of the process between two observations (package
#'    `levyBridge`).
#' 3. **Simulate the bridge** with the estimated parameters.
#' 4. **Pseudo-sample on a regular grid**: each trajectory is resampled on a
#'    regular grid, with the method `"Bridge"` (Grall-Maës, 2012, one grid per
#'    trajectory), `"ExtBridge"` (grid containing all the inspection times,
#'    step common to all the trajectories) or `"MultiBridge"` (one grid of
#'    Grall-Maës type common to all the trajectories). With several
#'    trajectories, `"ExtBridge"` and `"MultiBridge"` give one i.i.d. sample
#'    for all the trajectories together. The grid, the step and the drawing of
#'    the bridge of each method are detailed in the section Resampling methods
#'    of [gof.levy()].
#' 5. **Choose a test** for an i.i.d. sample of the law of the increments,
#'    see [tests_disponibles()].
#' 6. **Apply the test** to the reconstructed increments.
#'
#' With regular inspections the observed increments are already i.i.d.:
#' steps 1 to 4 are skipped. [gof.levy()] does all the steps;
#' [test_adequation()] does steps 5 and 6 for a sample of increments.
#'
#' @section Data:
#' [gof.levy()] accepts the formats of `levyfit` and `levyBridge`: a data
#' frame (times, levels) for one trajectory, a list of such data frames for
#' several trajectories with their own times, or a matrix (times, levels of
#' each trajectory), as returned by the simulation functions of `DegradSim`.
#'
#' @section Results:
#' * `res` prints the model, the type of inspections, the time step, the
#'   estimated parameters and one row per test;
#' * `res$results` is the data frame of the tests (statistic, p-value,
#'   critical value, decision);
#' * `res$bridge` contains the resampled trajectories and `res$increments`
#'   the tested increments.
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' Tito, A. R. and Paroissin, C. (in preparation).
#'
#' @examples
#' set.seed(1)
#' # A gamma trajectory observed at irregular times
#' t <- c(0, cumsum(round(runif(30, 0.5, 3.5), 1)))
#' x <- c(0, cumsum(rgamma(30, shape = 2 * diff(t), rate = 3)))
#'
#' # Is it a gamma process? Is it a Wiener process?
#' gof.levy(data.frame(t, x), "HGP", test = c("KS", "AD"), B = 100)
#' gof.levy(data.frame(t, x), "HWP")
#'
#' # The tests of each model
#' tests_disponibles()
#' @keywords internal
#' @importFrom Runuran ur pinvd.new udmeixner
"_PACKAGE"
