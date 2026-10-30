#' levygof: goodness-of-fit tests for Levy degradation processes
#'
#' `levygof` tests whether degradation trajectories come from a parametric
#' Levy process with unknown parameters:
#' * gamma (`"HGP"`);
#' * Wiener (`"HWP"`);
#' * inverse Gaussian (`"InvGaus"`);
#' * variance-gamma (`"VG"`);
#' * Meixner (`"Meixner"`).
#'
#' The main function is [gof.levy()].
#'
#' @section Functions of the package:
#' Click on a function name for its help page:
#'
#' | Function | Input | Role |
#' |---|---|---|
#' | [gof.levy()] | one or several trajectories | estimation, bridges if needed, tests |
#' | [test_adequation()] | i.i.d. sample of increments | tests only (steps 5 and 6 below) |
#' | [tests_disponibles()] | name of a model, or nothing | list of the tests of each model |
#'
#' Inspection times following a chosen pattern, and simulated trajectories
#' to try the tests, are given by the package `DegradSim`
#' (`DegradSim::inspection_times()` and its simulation functions).
#'
#' @section Methodology:
#' A Levy process has independent and stationary increments: the increments
#' over a same time step are i.i.d. The test proceeds in six steps.
#'
#' 1. **Estimate** the parameters from the observed trajectories by maximum
#'    likelihood (package `levyfit`).
#' 2. **Distribution of the bridge** of the process between two observations (package
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
#' 5. **Choose a test** for an i.i.d. sample of the distribution of the
#'    increments (list of the tests: [tests_disponibles()]).
#' 6. **Apply the test** to the reconstructed increments.
#'
#' With regular inspections the observed increments are already i.i.d.:
#' steps 1 to 4 are skipped. The function [gof.levy()] does all the steps,
#' and the function [test_adequation()] does steps 5 and 6 for a sample of
#' increments.
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
#' @examples
#' # 1. A gamma trajectory (a = 2, b = 3) with irregular inspections,
#' #    more and more frequent
#' set.seed(1)
#' t <- DegradSim::inspection_times("decreasing", n = 30, tmax = 60, digits = 1)$time
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3)
#'
#' # 2. Is it a gamma process? Is it a Wiener process?
#' gof.levy(x, "HGP", test = c("KS", "AD"), B = 100)
#' gof.levy(x, "HWP")
#'
#' # 3. The tests of each model
#' tests_disponibles()
#' @keywords package
#' @importFrom Runuran ur pinvd.new udmeixner
#' @importFrom DegradSim inspection_times
"_PACKAGE"
