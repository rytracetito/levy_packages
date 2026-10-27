#' Simulation of the compound Poisson process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' a compound Poisson process, a Levy process: a shock model where
#' each shock adds a random amount of degradation.
#'
#' @section Model:
#' Let \eqn{N} be a Poisson process with intensity \eqn{\lambda > 0} (the
#' times of the shocks) and \eqn{J_1, J_2, \dots} independent random variables
#' with the same distribution, independent of \eqn{N} (the sizes of the
#' shocks). The compound Poisson process is
#' \deqn{X_t = \sum_{k = 1}^{N_t} J_k, \qquad t \ge 0,}
#' with \eqn{X_t = 0} when \eqn{N_t = 0}. It is a Levy process: \eqn{X_0 = 0},
#' and for \eqn{0 \le s < t}, \eqn{X_t - X_s} is the sum of
#' \eqn{N_t - N_s \sim Poisson(\lambda (t - s))} independent jumps. Hence
#' \eqn{E[X_t] = \lambda t\, E[J]} and \eqn{Var(X_t) = \lambda t\, E[J^2]}.
#' The trajectories are piecewise constant, with a finite number of jumps on
#' every interval: an increment is 0 with probability at least
#' \eqn{e^{-\lambda (t - s)}} (no shock). The process is non-decreasing when
#' the jumps are non-negative, which is the usual case for a degradation
#' (cumulative shock model, Kahle et al., 2016, chapter 3). With jumps equal
#' to 1, it is the Poisson process \eqn{N} itself.
#'
#' @section Exact simulation:
#' Random walk approximation, described in [DegradSim-package]: with
#' \eqn{\Delta t_i = t_i - t_{i-1}}, the increments are drawn independently,
#' then cumulated:
#' 1. draw the numbers of shocks \eqn{N_i \sim Poisson(\lambda \Delta t_i)}
#'    with [stats::rpois()];
#' 2. draw \eqn{N_1 + \dots + N_n} jumps with `rjump`;
#' 3. \eqn{\Delta X_i} is the sum of the \eqn{N_i} jumps of the interval
#'    \eqn{]t_{i-1}, t_i]} (0 if \eqn{N_i = 0}).
#'
#' @inheritParams sim_gamma
#' @param lambda Intensity of the shocks, \eqn{\lambda > 0}: mean number of
#'   shocks per unit of time.
#' @param rjump Distribution of the jumps: an R function of `n` that returns
#'   `n` independent jumps, for instance `function(n) rgamma(n, 2, 3)`. The
#'   default `function(n) rexp(n, 1)` gives exponential jumps of mean 1, and
#'   `function(n) rep(1, n)` gives the Poisson process.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Kahle, W., Mercier, S. and Paroissin, C. (2016). *Degradation Processes in
#' Reliability*. ISTE and Wiley, London. ISBN 978-1-84821-888-8.
#'
#' Hoyle, A. E. V. (2010). *Information-Based Models for Finance and
#' Insurance*. PhD thesis, Imperial College London. arXiv:1010.0829.
#'
#' Sato, K. (2001). Basic results on Lévy processes. In O. E.
#' Barndorff-Nielsen, T. Mikosch and S. I. Resnick (eds.), *Lévy Processes:
#' Theory and Applications*, 3-37. Birkhäuser, Boston.
#' \doi{10.1007/978-1-4612-0197-7}
#'
#' @seealso [sim_gamma()], [plot_traj()], [traj_df()], [inspection_times()],
#'   [plot_plan()]
#'
#' @examples
#' # 1. One trajectory with 200 regular inspections up to 10 (step 0.05), exponential
#' #    jumps of mean 1 (default)
#' t <- inspection_times("regular", n = 200, tmax = 10)$time
#' set.seed(1)
#' x <- sim_cpoisson(t, lambda = 2)
#' x
#' plot_traj(x)
#'
#' # 2. Ten trajectories with gamma jumps of mean 2/3
#' x <- sim_cpoisson(t, lambda = 2, rjump = function(n) rgamma(n, 2, 3), n_traj = 10)
#' plot_traj(x)
#'
#' # Same times for all the trajectories: a matrix with the column t,
#' # then one column per trajectory (X1, X2, ...)
#' head(x[, 1:4])
#' x[, "X2"]
#'
#' # 3. Rare and large shocks, then frequent and small shocks: same mean
#' x1 <- sim_cpoisson(t, lambda = 0.2, rjump = function(n) rexp(n, 0.1), n_traj = 10)
#' x2 <- sim_cpoisson(t, lambda = 20, rjump = function(n) rexp(n, 10), n_traj = 10)
#' plot_traj(x1, "Rare and large shocks")
#' plot_traj(x2, "Frequent and small shocks")
#'
#' # 4. The Poisson process: jumps equal to 1
#' x <- sim_cpoisson(t, lambda = 2, rjump = function(n) rep(1, n), n_traj = 10)
#' plot_traj(x, "Poisson process")
#'
#' # 5. Irregular inspections, more and more frequent over time
#' plan <- inspection_times("decreasing", n = 50, n_traj = 5, tmax = 10)
#' x <- sim_cpoisson(plan, lambda = 2)
#' plot_traj(x)
#'
#' # 6. Different inspection times for each trajectory: a list of data frames
#' x <- sim_cpoisson(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), lambda = 2)
#' x[[2]]
#' x[[2]]$X
#'
#' # 7. All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#'
#' # 8. Increments over the constant step h = 1, with exponential jumps of
#' #    rate 1: an increment is 0 with probability exp(-lambda h) (no shock);
#' #    otherwise its density is sum_k P(N = k) dgamma(y, k, 1) / (1 - exp(-lambda h))
#' t8 <- inspection_times("regular", n = 100, tmax = 100)$time
#' x <- sim_cpoisson(t8, lambda = 2, n_traj = 10)
#' dx <- traj_df(x, increments = TRUE)$dX
#' dx <- dx[!is.na(dx)]
#' c(zero_observed = mean(dx == 0), zero_theory = exp(-2))
#' dpos <- function(y) {
#'   vapply(y, function(v) sum(dpois(1:30, 2) * dgamma(v, 1:30, 1)), 0) / (1 - exp(-2))
#' }
#' ggplot2::ggplot(data.frame(dx = dx[dx > 0]), ggplot2::aes(dx)) +
#'   ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
#'                           bins = 30, boundary = 0, fill = "grey80", colour = "grey40") +
#'   ggplot2::stat_function(fun = dpos, colour = "red") +
#'   ggplot2::labs(x = "dX > 0", y = "Density", title = "Non-zero increments, compound Poisson")
#' @export
sim_cpoisson <- function(t, lambda, rjump = function(n) stats::rexp(n, 1), n_traj = 1) {
  if (!is.numeric(lambda) || length(lambda) != 1 || !(lambda > 0) || !is.finite(lambda)) {
    stop("lambda must be a positive number")
  }
  if (!is.function(rjump)) stop("rjump must be an R function of n, for instance function(n) rexp(n, 1)")
  one_path <- function(t) {
    N <- stats::rpois(length(t) - 1, lambda * diff(t))
    dX <- numeric(length(N))
    if (sum(N) > 0) {
      J <- rjump(sum(N))
      if (!is.numeric(J) || length(J) != sum(N) || any(!is.finite(J))) {
        stop("rjump(n) must return n finite numbers")
      }
      pos <- N > 0
      dX[pos] <- rowsum(J, rep(seq_along(N), N))[, 1]
    }
    c(0, cumsum(dX))
  }
  x <- simulate_paths(t, n_traj, one_path, "compound Poisson", c(lambda = lambda))
  attr(x, "jumps") <- paste(deparse(rjump), collapse = " ")
  x
}
