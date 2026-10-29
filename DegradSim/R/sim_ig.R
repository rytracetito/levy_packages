#' Simulation of the inverse Gaussian process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' an inverse Gaussian (IG) process (a Levy process) or of a time-scaled IG
#' process, obtained with a time-scaling function \eqn{A}.
#'
#' @section Model:
#' \eqn{IG(m, \lambda)} denotes the inverse Gaussian distribution with mean
#' \eqn{m > 0} and shape \eqn{\lambda > 0}, of density
#' \deqn{\sqrt{\frac{\lambda}{2 \pi x^3}} \exp\Big(-\frac{\lambda (x - m)^2}{2 m^2 x}\Big), \qquad x > 0,}
#' and variance \eqn{m^3 / \lambda}.
#'
#' **IG process** (`A = NULL`, default). Let \eqn{\mu > 0} and
#' \eqn{\eta > 0}. The process satisfies \eqn{X_0 = 0}, has independent
#' increments, and for \eqn{0 \le s < t},
#' \deqn{X_t - X_s \sim IG\big(\mu (t - s),\ \eta (t - s)^2\big),}
#' with the parametrization of Wang and Xu (2010). Then \eqn{E[X_t] = \mu t}
#' and \eqn{Var(X_t) = \mu^3 t / \eta}. It is a Levy process, increasing and
#' purely made of jumps.
#'
#' **Time-scaled IG process** (`A` given). With the time-scaling function
#' \eqn{A} defined in [DegradSim-package], \eqn{X_t = Y_{A(t)}}, where \eqn{Y}
#' is the IG process with \eqn{\mu = 1} (identifiability, see
#' [DegradSim-package]):
#' \deqn{X_t - X_s \sim IG\big(A(t) - A(s),\ \eta\, (A(t) - A(s))^2\big),}
#' \eqn{E[X_t] = A(t)} and \eqn{Var(X_t) = A(t) / \eta}. This is the
#' time-scaled IG process with mean function \eqn{A} of Cha and Mercier
#' (2019, Example 3), whose parameter \eqn{b} is \eqn{\eta}.
#'
#' @section Exact simulation:
#' Random walk approximation, described in [DegradSim-package]: the
#' increments \eqn{\Delta X_i \sim IG(m_i, \lambda_i)}, with
#' \eqn{m_i = \mu \Delta t_i}, \eqn{\lambda_i = \eta \Delta t_i^2}
#' (IG process) or \eqn{m_i = \Delta A_i},
#' \eqn{\lambda_i = \eta \Delta A_i^2} (time-scaled IG process), are drawn
#' independently, then cumulated. Each increment \eqn{IG(m, \lambda)} is drawn
#' by the method of Michael, Schucany and Haas (1976), based on
#' \eqn{\lambda (X - m)^2 / (m^2 X) \sim \chi^2_1}: draw
#' \eqn{Z \sim \mathcal{N}(0, 1)}, set \eqn{\nu = Z^2},
#' \eqn{x_1 = m + m^2 \nu / (2\lambda) - m \sqrt{4 m \lambda \nu + m^2 \nu^2} / (2\lambda)},
#' and return \eqn{x_1} with probability \eqn{m / (m + x_1)}, otherwise
#' \eqn{m^2 / x_1}.
#'
#' @inheritParams sim_gamma
#' @param mu Mean degradation per unit of time of the IG process,
#'   \eqn{\mu > 0}. Not given when `A` is given (the mean of \eqn{X_t} is then
#'   \eqn{A(t)}).
#' @param eta Shape parameter, \eqn{\eta > 0}: the larger \eqn{\eta}, the less
#'   variable the increments. The variance of \eqn{X_t} is
#'   \eqn{\mu^3 t / \eta} (IG process) or \eqn{A(t) / \eta} (time-scaled
#'   IG process).
#'
#' @inherit sim_gamma return
#'
#' @references
#' Cha, J. H. and Mercier, S. (2019). Transformed Lévy processes as
#' state-dependent wear models. *Advances in Applied Probability*, 51(2),
#' 468-486. \doi{10.1017/apr.2019.21}
#'
#' Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
#' variates using transformations with multiple roots. *The American
#' Statistician*, 30(2), 88-90.\cr
#' \doi{10.1080/00031305.1976.10479147}
#'
#' Wang, X. and Xu, D. (2010). An inverse Gaussian process model for
#' degradation data. *Technometrics*, 52(2), 188-197.
#' \doi{10.1198/TECH.2009.08197}
#'
#' @seealso
#' * [sim_gamma()], [sim_wiener()]: the other time-scaled processes;
#' * [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()].
#'
#' @examples
#' # 1. One trajectory with 200 regular inspections up to 10 (step 0.05)
#' t <- inspection_times("regular", n = 200, tmax = 10)$time
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
#' plan <- inspection_times("increasing", n = 50, n_traj = 5, tmax = 10)
#' x <- sim_ig(plan, mu = 0.5, eta = 2)
#' plot_traj(x)
#'
#' # 4. Different inspection times for each trajectory: a list of data frames
#' x <- sim_ig(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), mu = 0.5, eta = 2)
#' x[[2]]
#' x[[2]]$X
#'
#' # 5. Time-scaled IG process with mean function A(t) = 0.5 t^1.5
#' #    (Cha and Mercier, 2019, Example 3): X_t - X_s ~ IG(dA, eta dA^2)
#' x <- sim_ig(t, eta = 2, A = function(t) 0.5 * t^1.5, n_traj = 10)
#' x
#' plot_traj(x, "Time-scaled inverse Gaussian process")
#'
#' # 6. With A(t) = mu t, the time-scaled IG process with eta is the IG
#' #    process with parameters mu and eta mu^2
#' set.seed(2)
#' y1 <- sim_ig(t, eta = 2, A = function(t) 0.5 * t)
#' set.seed(2)
#' y2 <- sim_ig(t, mu = 0.5, eta = 2 * 0.5^2)
#' all.equal(unclass(y1)[, 2], unclass(y2)[, 2])
#'
#' # 7. All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#'
#' # 8. Histogram of the increments, with the density of the increments in red.
#' #    Constant step h = 0.5: the increments are i.i.d.
#' #    with distribution IG(m, lambda), m = mu h and lambda = eta h^2
#' t8 <- inspection_times("regular", n = 50, tmax = 25)$time
#' x <- sim_ig(t8, mu = 0.5, eta = 2, n_traj = 20)
#' dx <- traj_df(x, increments = TRUE)$dX
#' dx <- dx[!is.na(dx)]
#' m <- 0.5 * 0.5
#' lambda <- 2 * 0.5^2
#' dig <- function(y) sqrt(lambda / (2 * pi * y^3)) * exp(-lambda * (y - m)^2 / (2 * m^2 * y))
#' ggplot2::ggplot(data.frame(dx), ggplot2::aes(dx)) +
#'   ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
#'                           bins = 30, boundary = 0, fill = "grey80", colour = "grey40") +
#'   ggplot2::stat_function(fun = dig, colour = "red") +
#'   ggplot2::labs(x = "dX", y = "Density",
#'                 title = "Increments of the inverse Gaussian process")
#' @export
sim_ig <- function(t, mu, eta, A = NULL, method = "RWA", type = NULL, n_traj = 1) {
  check_sim_method(method, type)
  if (missing(eta) || !is.numeric(eta) || length(eta) != 1 || !(eta > 0)) stop("eta must be a positive number")
  if (is.null(A)) {
    if (missing(mu)) stop("mu is missing: give mu (IG process) or A (time-scaled IG process)")
    if (!is.numeric(mu) || length(mu) != 1 || !(mu > 0)) stop("mu must be a positive number")
    clock <- function(t) diff(t)
    drift <- mu
    model <- "inverse Gaussian"
    parameters <- c(mu = mu, eta = eta)
  } else {
    if (!missing(mu)) stop("give mu or A, not both: with A, the mean of X_t is A(t) (mu = 1)")
    clock <- function(t) time_change(A, t)
    drift <- 1
    model <- "time-scaled inverse Gaussian"
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
