#' Simulation of the inverse Gaussian process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' a homogeneous inverse Gaussian (IG) process (a Lévy process) or of a
#' non-homogeneous IG process, obtained by a time change \eqn{A}.
#'
#' @section Model:
#' \eqn{IG(m, \lambda)} denotes the inverse Gaussian distribution with mean
#' \eqn{m > 0} and shape \eqn{\lambda > 0}, of density
#' \deqn{\sqrt{\frac{\lambda}{2 \pi x^3}} \exp\Big(-\frac{\lambda (x - m)^2}{2 m^2 x}\Big), \qquad x > 0,}
#' and variance \eqn{m^3 / \lambda}.
#'
#' **Homogeneous IG process** (`A = NULL`, default). Let \eqn{\mu > 0} and
#' \eqn{\eta > 0}. The process satisfies \eqn{X(0) = 0}, has independent
#' increments, and for \eqn{0 \le s < t},
#' \deqn{X(t) - X(s) \sim IG\big(\mu (t - s),\ \eta (t - s)^2\big),}
#' with the parametrization of Wang and Xu (2010). Then \eqn{E[X(t)] = \mu t}
#' and \eqn{Var(X(t)) = \mu^3 t / \eta}. It is a Lévy process, increasing and
#' purely made of jumps.
#'
#' **Non-homogeneous IG process** (`A` given). Let \eqn{A} be a time change: a
#' non-decreasing function, right-continuous with left limits, with
#' \eqn{A(0) = 0} and \eqn{A(t) \to +\infty} as \eqn{t \to +\infty}. The
#' process is \eqn{X(t) = Y(A(t))}, where \eqn{Y} is the homogeneous IG
#' process with \eqn{\mu = 1}, so that
#' \deqn{X(t) - X(s) \sim IG\big(A(t) - A(s),\ \eta\, (A(t) - A(s))^2\big),}
#' \eqn{E[X(t)] = A(t)} and \eqn{Var(X(t)) = A(t) / \eta}. This is the
#' non-homogeneous IG process with mean function \eqn{A} of Cha and Mercier
#' (2019, Example 3), whose parameter \eqn{b} is \eqn{\eta}.
#'
#' **Why \eqn{\mu = 1} in the non-homogeneous case.** With a mean
#' \eqn{\mu \Delta A}, replacing \eqn{(A, \mu, \eta)} by
#' \eqn{(c A, \mu / c, \eta / c^2)}, \eqn{c > 0}, gives the same distribution
#' \eqn{IG(\mu \Delta A, \eta \Delta A^2)} of the increments: \eqn{\mu},
#' \eqn{\eta} and \eqn{A} could not be identified separately. The mean is therefore contained in \eqn{A}. With
#' \eqn{A(t) = \mu t}, the non-homogeneous process has increments
#' \eqn{IG(\mu \Delta t, \eta \mu^2 \Delta t^2)}: the homogeneous process with
#' parameters \eqn{(\mu, \eta)} is obtained with \eqn{A(t) = \mu t} and
#' \eqn{\eta / \mu^2}.
#'
#' @section Exact simulation at the inspection times (random walk approximation):
#' The increments over disjoint intervals are independent with known
#' distributions, so the values at the inspection times
#' \eqn{0 = t_0 < t_1 < \dots < t_n} are simulated exactly by the random walk
#' approximation (RWA) of Kahle et al. (2016, section 2.5.1, written there for the gamma process): independent
#' increments are drawn, then cumulated.
#' 1. \eqn{m_i = \mu \Delta t_i} and \eqn{\lambda_i = \eta \Delta t_i^2}
#'    (homogeneous case), or \eqn{m_i = \Delta A_i} and
#'    \eqn{\lambda_i = \eta \Delta A_i^2} (non-homogeneous case), with
#'    \eqn{\Delta t_i = t_i - t_{i-1}} and \eqn{\Delta A_i = A(t_i) - A(t_{i-1})}.
#' 2. Draw independently \eqn{\Delta X_i \sim IG(m_i, \lambda_i)} by the method
#'    of Michael, Schucany and Haas (1976), based on
#'    \eqn{\lambda (X - m)^2 / (m^2 X) \sim \chi^2_1}: draw
#'    \eqn{Z \sim N(0, 1)}, set \eqn{\nu = Z^2},
#'    \eqn{x_1 = m + m^2 \nu / (2\lambda) - m \sqrt{4 m \lambda \nu + m^2 \nu^2} / (2\lambda)},
#'    and return \eqn{x_1} with probability \eqn{m / (m + x_1)}, otherwise
#'    \eqn{m^2 / x_1}.
#' 3. \eqn{X(t_0) = 0} and \eqn{X(t_i) = \sum_{j \le i} \Delta X_j}.
#'
#' The values at the inspection times have exactly the distribution of the
#' process; only the path between two inspection times is not simulated. If
#' \eqn{A} is constant between two inspection times, the increment is 0.
#'
#' @inheritSection sim_gamma Conditions on the inspection times and on A
#'
#' @inheritParams sim_gamma
#' @param mu Mean degradation per unit of time of the homogeneous process,
#'   \eqn{\mu > 0}. Not given when `A` is given (the mean of \eqn{X(t)} is then
#'   \eqn{A(t)}).
#' @param eta Shape parameter, \eqn{\eta > 0}: the larger \eqn{\eta}, the less
#'   variable the increments. The variance of \eqn{X(t)} is
#'   \eqn{\mu^3 t / \eta} (homogeneous case) or \eqn{A(t) / \eta}
#'   (non-homogeneous case).
#' @param A `NULL` (default) for the homogeneous process, or the time change
#'   for the non-homogeneous process: a vectorized R function of `t`,
#'   non-decreasing with \eqn{A(0) = 0}, for instance
#'   `function(t) 0.5 * t^1.5`.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Cha, J. H. and Mercier, S. (2019). Transformed Lévy processes as
#' state-dependent wear models. *Advances in Applied Probability*, 51(2),
#' 468-486. \doi{10.1017/apr.2019.21}
#'
#' Kahle, W., Mercier, S. and Paroissin, C. (2016). *Degradation Processes in
#' Reliability*. ISTE and Wiley, London. ISBN 978-1-84821-888-8.
#'
#' Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
#' variates using transformations with multiple roots. *The American
#' Statistician*, 30(2), 88-90. \doi{10.1080/00031305.1976.10479147}
#'
#' Wang, X. and Xu, D. (2010). An inverse Gaussian process model for
#' degradation data. *Technometrics*, 52(2), 188-197.
#' \doi{10.1198/TECH.2009.08197}
#'
#' @seealso [sim_gamma()] and [sim_wiener()] for the other processes with a
#'   time change, [plot_traj()], [traj_df()], [inspection_times()],
#'   [plot_plan()]
#'
#' @examples
#' # 1. Homogeneous process: one trajectory inspected every 0.5 time unit
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_ig(t, mu = 0.5, eta = 2)
#' x
#' plot_traj(x)
#'
#' # 2. Ten trajectories with the same inspection times
#' x <- sim_ig(t, mu = 0.5, eta = 2, n_traj = 10)
#' plot_traj(x)
#'
#' # Same times for all the trajectories: a matrix with the column t,
#' # then one column per trajectory (X1, X2, ...)
#' head(x[, 1:4])
#' x[, "X2"]
#'
#' # 3. Irregular inspections, less and less frequent over time
#' plan <- inspection_times("decreasing", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_ig(plan, mu = 0.5, eta = 2)
#' plot_traj(x)
#'
#' # 4. Different inspection times for each trajectory: a list of data frames
#' x <- sim_ig(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), mu = 0.5, eta = 2)
#' x[[2]]
#' x[[2]]$X
#'
#' # 5. Non-homogeneous IG process with mean function A(t) = 0.5 t^1.5
#' #    (Cha and Mercier, 2019, Example 3): X(t) - X(s) ~ IG(dA, eta dA^2)
#' x <- sim_ig(t, eta = 2, A = function(t) 0.5 * t^1.5, n_traj = 10)
#' x
#' plot_traj(x, "Non-homogeneous inverse Gaussian process")
#'
#' # 6. The homogeneous process (mu, eta) is A(t) = mu t with eta / mu^2
#' set.seed(2)
#' y1 <- sim_ig(t, mu = 0.5, eta = 2)
#' set.seed(2)
#' y2 <- sim_ig(t, eta = 2 / 0.5^2, A = function(t) 0.5 * t)
#' all.equal(unclass(y1)[, 2], unclass(y2)[, 2])
#'
#' # 7. All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#' @export
sim_ig <- function(t, mu, eta, A = NULL, n_traj = 1) {
  if (missing(eta) || !is.numeric(eta) || length(eta) != 1 || !(eta > 0)) stop("eta must be a positive number")
  if (is.null(A)) {
    if (missing(mu)) stop("mu is missing: give mu (homogeneous process) or A (non-homogeneous process)")
    if (!is.numeric(mu) || length(mu) != 1 || !(mu > 0)) stop("mu must be a positive number")
    clock <- function(t) diff(t)
    drift <- mu
    model <- "inverse Gaussian"
    parameters <- c(mu = mu, eta = eta)
  } else {
    if (!missing(mu)) stop("give mu or A, not both: with A, the mean of X(t) is A(t) (mu = 1)")
    clock <- function(t) time_change(A, t)
    drift <- 1
    model <- "non-homogeneous inverse Gaussian"
    parameters <- c(eta = eta)
  }
  one_path <- function(t) {
    dA <- clock(t)
    dX <- numeric(length(dA))
    pos <- dA > 0
    dX[pos] <- rig(sum(pos), drift * dA[pos], eta * dA[pos]^2)
    c(0, cumsum(dX))
  }
  simulate_paths(t, n_traj, one_path, model, parameters, A)
}

rig <- function(n, m, lambda) {
  nu <- stats::rnorm(n)^2
  x1 <- m + m^2 * nu / (2 * lambda) - m / (2 * lambda) * sqrt(4 * m * lambda * nu + m^2 * nu^2)
  ifelse(stats::runif(n) <= m / (m + x1), x1, m^2 / x1)
}
