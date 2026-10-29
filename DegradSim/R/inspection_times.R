#' Inspection times following a chosen pattern
#'
#' Generates inspection times on \eqn{[0, t_{max}]} for one or several
#' trajectories, with gaps between two successive inspections that are
#' constant, increasing, decreasing, or that vary between the beginning, the
#' middle and the end of the period. The result can be given directly to the
#' simulation functions, such as [sim_gamma()], and plotted with [plot_plan()].
#'
#' @details
#' **Construction.** The times are the values of a strictly increasing
#' function \eqn{f} on the regular grid \eqn{x_i = i / n},
#' \eqn{i = 0, \dots, n}, of \eqn{[0, 1]}, rescaled to \eqn{[0, t_{max}]}:
#' \deqn{t_i = t_{max} \frac{f(x_i) - f(0)}{f(1) - f(0)}, \qquad i = 0, \dots, n.}
#' So \eqn{t_0 = 0}, \eqn{t_n = t_{max}}, and the gap between two successive
#' inspections is about proportional to the slope \eqn{f'(x_i)}: the
#' inspections are frequent where \eqn{f} is flat and rare where \eqn{f} is
#' steep.
#'
#' **Patterns.**
#' * `"regular"`: \eqn{f(x) = x}, constant gaps \eqn{t_{max} / n};
#' * `"decreasing"`: \eqn{f(x) = \ln(1 + 10 x)}, decreasing gaps: inspections
#'   rare at the beginning and more and more frequent; the first gap is about
#'   \eqn{f'(0) / f'(1) = 11} times the last one;
#' * `"increasing"`: \eqn{f(x) = -\ln(11 - 10 x)}, increasing gaps:
#'   inspections frequent at the beginning and rarer and rarer, mirror image
#'   of `"decreasing"`;
#' * `"dec-inc"`: \eqn{f(x) = x + \frac{a}{2\pi} \sin(2\pi x)}, gaps that
#'   decrease then increase: inspections rare at the beginning, frequent in
#'   the middle, rare at the end;
#' * `"inc-dec"`: \eqn{f(x) = x - \frac{a}{2\pi} \sin(2\pi x)}, gaps that
#'   increase then decrease: inspections frequent at the beginning, rare in
#'   the middle, frequent at the end.
#'
#' The names describe the gaps between successive inspections:
#' `"decreasing"` means gaps that decrease, so more and more frequent
#' inspections; `"dec-inc"` means gaps that decrease then increase.
#'
#' **Amplitude `a` of `"dec-inc"` and `"inc-dec"`.** Here
#' \eqn{f'(x) = 1 \pm a \cos(2\pi x)}: the gaps vary between
#' \eqn{1 - a} and \eqn{1 + a} times the mean gap \eqn{t_{max} / n}, so that
#' the largest gap is \eqn{(1 + a) / (1 - a)} times the smallest one. With the
#' default \eqn{a = 0.5}, this ratio is 3 and the gaps change progressively;
#' with \eqn{a = 0.8}, it is 9 and the frequent inspections gather in a short
#' time; \eqn{a = 0} gives regular inspections.
#'
#' **Rounding.** With `digits`, the times are rounded to `digits` decimal
#' places. If the rounding makes two successive times equal, the repeated
#' times are removed, so that the trajectory has fewer than `n` inspections
#' after time 0, and a warning gives the number of inspections that are kept.
#'
#' **Other patterns.** Any strictly increasing R function on \eqn{[0, 1]} can
#' be given as `pattern`, for instance `sqrt`; an error is raised if it is not
#' strictly increasing on the grid.
#'
#' @param pattern Pattern of the inspection times (see Details):
#'   * `"regular"` (default), `"decreasing"` or `"increasing"`;
#'   * `"dec-inc"` or `"inc-dec"`;
#'   * or a strictly increasing R function on \eqn{[0, 1]}.
#' @param n Number of inspections after time 0, a positive integer. A single
#'   value gives the same number of inspections to all the trajectories; a
#'   vector of length `n_traj` gives one number per trajectory.
#' @param n_traj Number of trajectories, a positive integer.
#' @param tmax Time of the last inspection, \eqn{t_{max} > 0}. A single value
#'   gives the same last time to all the trajectories; a vector of length
#'   `n_traj` gives one value per trajectory.
#' @param digits Number of decimal places of the times, a non-negative
#'   integer, or `NULL` (default) for no rounding. Times made equal by the
#'   rounding are removed, with a warning (see Details).
#' @param a Amplitude of the patterns `"dec-inc"` and `"inc-dec"`, in
#'   \eqn{[0, 1[}: it sets the contrast between frequent and rare inspections,
#'   from 0 (regular inspections) to 1 excluded; default 0.5. Not used by the
#'   other patterns.
#'
#' @return A data frame with one row per inspection and two columns: `id`
#'   (number of the trajectory, 1 to `n_traj`) and `time` (inspection time).
#'   Each trajectory starts at time 0 and ends at its `tmax`, with strictly
#'   increasing times. The times of the trajectory \eqn{k} are
#'   `plan$time[plan$id == k]`.
#'
#' @seealso [plot_plan()] to plot the inspection times, [sim_gamma()] and the
#'   other simulation functions, which accept the result as inspection times.
#'
#' @examples
#' # 1. Regular inspections: 10 inspections after time 0, up to 100
#' inspection_times("regular", n = 10, tmax = 100)
#'
#' # 2. More and more frequent inspections, times rounded to integers
#' inspection_times("decreasing", n = 10, digits = 0)
#'
#' # 3. Fewer and fewer inspections
#' inspection_times("increasing", n = 10, digits = 0)
#'
#' # 4. Rare at the beginning, frequent in the middle, rare at the end
#' inspection_times("dec-inc", n = 20, digits = 1)
#'
#' # 5. Frequent at the beginning and at the end, rare in the middle
#' inspection_times("inc-dec", n = 20, digits = 1)
#'
#' # 6. Stronger contrast between frequent and rare inspections
#' inspection_times("inc-dec", n = 20, digits = 1, a = 0.8)
#'
#' # Gaps between inspections: from about 2.5 to 7.5 with the default amplitude
#' range(diff(inspection_times("inc-dec", n = 20)$time))
#'
#' # 7. Three trajectories with 10, 15 and 20 inspections, up to 50, 100 and 150
#' plan <- inspection_times("decreasing", n = c(10, 15, 20), n_traj = 3,
#'                          tmax = c(50, 100, 150), digits = 1)
#' plan
#'
#' # Times of the second trajectory only
#' plan$time[plan$id == 2]
#'
#' # 8. Same number of inspections and same tmax: the same times for the
#' #    three trajectories
#' inspection_times("decreasing", n = 5, n_traj = 3, tmax = 10)
#'
#' # 9. Rounding too strong: the repeated times are removed with a warning
#' inspection_times("increasing", n = 30, tmax = 20, digits = 0)
#'
#' # 10. Any strictly increasing function on [0, 1]: with sqrt, the
#' #     inspections are rare at the beginning and more and more frequent
#' inspection_times(sqrt, n = 10, digits = 0)
#'
#' # 11. Mixture of the increasing and decreasing patterns: frequent
#' #     inspections around a point that moves to the right when w increases
#' #     (0 < w < 1)
#' w <- 0.8
#' f <- function(x) w * log(1 + 10 * x) - (1 - w) * log(11 - 10 * x)
#' inspection_times(f, n = 20, digits = 1)
#'
#' # 12. Plot the inspection times, one row per trajectory, then simulate
#' plan <- inspection_times("inc-dec", n = c(10, 20, 30), n_traj = 3, tmax = 100)
#' plot_plan(plan)
#' x <- sim_gamma(plan, a = 0.2, b = 1)
#' plot_traj(x)
#' @export
inspection_times <- function(pattern = "regular", n = 10,
                             n_traj = 1, tmax = 100, digits = NULL, a = 0.5) {
  is_count <- function(v) is.numeric(v) && all(is.finite(v)) && all(v >= 1) && all(v == round(v))
  if (length(n_traj) != 1 || !is_count(n_traj)) stop("n_traj must be a positive integer")
  if (!length(n) %in% c(1, n_traj) || !is_count(n)) {
    stop("n must be a positive integer, or a vector of n_traj positive integers")
  }
  if (!length(tmax) %in% c(1, n_traj) || !is.numeric(tmax) || any(!is.finite(tmax)) || any(tmax <= 0)) {
    stop("tmax must be a positive number, or a vector of n_traj positive numbers")
  }
  if (!is.null(digits) && (length(digits) != 1 || !is.numeric(digits) || digits < 0 || digits != round(digits))) {
    stop("digits must be NULL or a non-negative integer")
  }
  if (!is.numeric(a) || length(a) != 1 || is.na(a) || a < 0 || a >= 1) {
    stop("a must be a number in [0, 1)")
  }

  if (is.function(pattern)) {
    f <- pattern
  } else {
    f <- switch(pattern,
      "regular"    = function(x) x,
      "decreasing" = function(x) log(1 + 10 * x),
      "increasing" = function(x) -log(11 - 10 * x),
      "dec-inc"    = function(x) x + a / (2 * pi) * sin(2 * pi * x),
      "inc-dec"    = function(x) x - a / (2 * pi) * sin(2 * pi * x),
      stop("unknown pattern: ", pattern,
           " (use \"regular\", \"increasing\", \"decreasing\", \"inc-dec\", \"dec-inc\" or a function)"))
  }

  n <- rep_len(n, n_traj)
  tmax <- rep_len(tmax, n_traj)

  resultat <- NULL
  for (j in 1:n_traj) {
    fx <- f(seq(from = 0, to = 1, length.out = n[j] + 1))
    if (length(fx) != n[j] + 1 || any(!is.finite(fx)) || any(diff(fx) <= 0)) {
      stop("the pattern function must be vectorized and strictly increasing on [0, 1]")
    }
    tps <- tmax[j] * (fx - fx[1]) / (fx[n[j] + 1] - fx[1])
    if (!is.null(digits)) {
      tps <- unique(round(tps, digits))
      if (length(tps) < 2) stop("the rounding is too strong: all the times of trajectory ", j, " are rounded to 0")
      if (length(tps) < n[j] + 1) {
        warning("trajectory ", j, ": the rounding to ", digits, " decimal place(s) made some times equal; ",
                length(tps) - 1, " inspections after time 0 are kept instead of ", n[j], call. = FALSE)
      }
    }
    resultat <- rbind(resultat, data.frame(id = j, time = tps))
  }
  resultat
}
