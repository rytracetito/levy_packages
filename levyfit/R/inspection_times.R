#' Inspection times with a chosen pattern
#'
#' Generates inspection times on \eqn{[0, t_{max}]} for one or several
#' trajectories, with a number of inspections per unit of time that is
#' constant, increasing, decreasing or alternating.
#'
#' @details
#' The times are the values of an increasing function \eqn{f} on a regular grid
#' \eqn{x_i = i\, x_{max}/n}, \eqn{i = 0, \dots, n}, with \eqn{x_{max} = 10},
#' rescaled to \eqn{[0, t_{max}]}:
#' \deqn{t_i = t_{max} \frac{f(x_i) - f(0)}{f(x_{max}) - f(0)}.}
#' The gap between two inspections is about \eqn{f'(x)\, x_{max}/n}: the
#' inspections are frequent where \eqn{f} is flat and rare where \eqn{f} is
#' steep. The patterns are
#' * `"regular"`: \eqn{f(x) = x}, constant gaps;
#' * `"increasing"`: \eqn{f(x) = \ln(x + 1)}, more and more inspections over time;
#' * `"decreasing"`: \eqn{f(x) = -\ln(x_{max} + 1 - x)}, fewer and fewer
#'   inspections, mirror image of `"increasing"`;
#' * `"inc-dec"`: \eqn{f(x) = x + 0.8\, \frac{x_{max}}{2\pi} \sin(2\pi x / x_{max})},
#'   the number of inspections increases then decreases;
#' * `"dec-inc"`: same with a minus sign, the number of inspections decreases
#'   then increases;
#' * `"inc-dec-inc"`: \eqn{f(x) = x + 0.8\, \frac{x_{max}}{3\pi} \sin(3\pi x / x_{max})},
#'   the number of inspections increases, decreases, then increases;
#' * `"dec-inc-dec"`: same with a minus sign.
#'
#' Any other strictly increasing function on \eqn{[0, 10]} can be given as
#' `pattern`.
#'
#' @param pattern Name of the pattern (see Details) or a strictly increasing
#'   function on \eqn{[0, 10]}.
#' @param n Number of inspections after time 0, for all trajectories or one
#'   per trajectory.
#' @param n_traj Number of trajectories.
#' @param tmax Time of the last inspection, for all trajectories or one per
#'   trajectory.
#' @param digits Number of decimal places of the times, `NULL` for no rounding.
#'   An error is raised if rounding makes two times equal.
#'
#' @return A data frame with columns `id` (trajectory) and `time`, time 0
#'   included.
#'
#' @examples
#' # Regular inspections: 10 inspections after time 0, up to 100
#' inspection_times("regular", n = 10, tmax = 100)
#'
#' # More and more frequent inspections, times rounded to integers
#' inspection_times("increasing", n = 10, digits = 0)
#'
#' # Fewer and fewer inspections
#' inspection_times("decreasing", n = 10, digits = 0)
#'
#' # Frequent in the middle of the period, rare at the start and at the end
#' inspection_times("inc-dec", n = 20, digits = 1)
#'
#' # Frequent at the start and at the end, rare in the middle
#' inspection_times("dec-inc", n = 20, digits = 1)
#'
#' # Three trajectories with 10, 15 and 20 inspections, up to 50, 100 and 150
#' plan <- inspection_times("increasing", n = c(10, 15, 20), n_traj = 3,
#'                          tmax = c(50, 100, 150), digits = 1)
#' plan
#'
#' # Times of the second trajectory only
#' plan$time[plan$id == 2]
#'
#' # Any strictly increasing function on [0, 10]: with sqrt, the inspections
#' # are rare at the start and more and more frequent
#' inspection_times(sqrt, n = 10, digits = 0)
#' @export
inspection_times <- function(pattern = "regular", n = 10,
                             n_traj = 1, tmax = 100, digits = NULL) {

  xmax <- 10

  if (is.function(pattern)) {
    f <- pattern
  } else {
    f <- switch(pattern,
      "regular"     = function(x) x,
      "increasing"  = function(x) log(x + 1),
      "decreasing"  = function(x) -log(xmax + 1 - x),
      "inc-dec"     = function(x) x + 0.8 * xmax / (2 * pi) * sin(2 * pi * x / xmax),
      "dec-inc"     = function(x) x - 0.8 * xmax / (2 * pi) * sin(2 * pi * x / xmax),
      "inc-dec-inc" = function(x) x + 0.8 * xmax / (3 * pi) * sin(3 * pi * x / xmax),
      "dec-inc-dec" = function(x) x - 0.8 * xmax / (3 * pi) * sin(3 * pi * x / xmax),
      stop("unknown pattern: ", pattern))
  }

  n <- rep_len(n, n_traj)
  tmax <- rep_len(tmax, n_traj)

  resultat <- NULL
  for (j in 1:n_traj) {
    x.pts <- seq(from = 0, to = xmax, length.out = n[j] + 1)
    tps <- f(x.pts)
    tps <- tmax[j] * (tps - tps[1]) / (tps[n[j] + 1] - tps[1])
    if (!is.null(digits)) tps <- round(tps, digits)
    if (any(is.na(tps)) || any(diff(tps) <= 0)) stop("two times are equal or decreasing: f is not increasing or the rounding is too strong")
    resultat <- rbind(resultat, data.frame(id = j, time = tps))
  }
  resultat
}
