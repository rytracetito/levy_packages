#' Maximum likelihood estimation for the gamma process
#'
#' Computes the maximum likelihood estimator (MLE) of the shape \eqn{a} and
#' the rate \eqn{b} of the homogeneous gamma process from one or several
#' trajectories observed at possibly irregular inspection times. The MLE always exists
#' (with probability one), is unique, and is the root of an equation in one
#' variable.
#'
#' @section Model and data:
#' The gamma process is increasing, with independent increments
#' \eqn{X_t - X_s \sim \Gamma(a (t - s), b)} (shape \eqn{a(t - s)}, rate
#' \eqn{b}), so that \eqn{E[X_t] = a t / b} and \eqn{Var(X_t) = a t / b^2}.
#' The trajectories are independent and have the same law; each one is
#' observed at its own inspection times. Let
#' \eqn{(\Delta t_j, \Delta x_j)}, \eqn{j = 1, \dots, m}, be the durations and
#' degradation increments (\eqn{\Delta x_j > 0}) between two consecutive
#' inspections of a same trajectory, gathered over all the trajectories,
#' \eqn{T = \sum_j \Delta t_j}, \eqn{S = \sum_j \Delta x_j},
#' \eqn{U = \sum_j \Delta t_j \ln \Delta x_j}, the speeds
#' \eqn{r_j = \Delta x_j / \Delta t_j} and the weights
#' \eqn{w_j = \Delta t_j / T}.
#'
#' @section Likelihood:
#' The increments being independent, the log-likelihood is
#' \deqn{\ell(a, b) = \sum_j \big[(a \Delta t_j - 1) \ln \Delta x_j - \ln \Gamma(a \Delta t_j)
#'   + a \Delta t_j \ln b - b \Delta x_j\big].}
#' It is a two-parameter exponential family with sufficient statistic
#' \eqn{(U, S)}.
#'
#' @section Maximum likelihood estimator:
#' 1. The score in \eqn{b} is \eqn{\partial \ell / \partial b = a T / b - S}:
#'    for a fixed \eqn{a}, \eqn{\ell} is maximal at \eqn{\hat b(a) = a T / S}.
#' 2. The profile log-likelihood \eqn{\ell_p(a) = \ell(a, \hat b(a))} has
#'    derivative
#'    \deqn{g(a) = U - \sum_j \Delta t_j\, \psi(a \Delta t_j) + T \ln \frac{a T}{S},}
#'    where \eqn{\psi} is the digamma function.
#' 3. Since \eqn{\psi'(x) > 1/x} for \eqn{x > 0}, \eqn{g} is strictly
#'    decreasing. Moreover \eqn{g(a) \to +\infty} as \eqn{a \to 0} and, by
#'    Jensen's inequality,
#'    \eqn{g(a) \to T \big(\sum_j w_j \ln r_j - \ln \sum_j w_j r_j\big) < 0}
#'    as \eqn{a \to \infty}, unless all the speeds are equal.
#' 4. Hence \eqn{g} has a unique root \eqn{\hat a}, and
#'    \eqn{\hat b = \hat a T / S}. If all the speeds are equal (in particular
#'    if \eqn{m = 1}), the MLE does not exist (`NA` with a warning); this has
#'    probability zero when \eqn{m \ge 2}.
#'
#' The root is computed by [stats::uniroot()] in \eqn{\ln a}, starting from the
#' method of moments: with \eqn{\hat r = S / T} and
#' \eqn{\hat\sigma^2 = \sum_j (\Delta x_j - \hat r \Delta t_j)^2 / (T - \sum_j \Delta t_j^2 / T)},
#' \eqn{a_{MM} = \hat r^2 / \hat\sigma^2}. Since \eqn{g} is monotone, the
#' search converges whatever the starting value.
#'
#' @section Bias:
#' No unbiased version of \eqn{\hat a} and \eqn{\hat b} is available. With few
#' increments (about ten) both are biased upwards, by about 40 percent; the
#' bias vanishes as the number of increments grows. Their ratio
#' \eqn{\hat a / \hat b = S / T} is unbiased.
#'
#' @inheritSection wiener_mle Confidence intervals and estimation failures
#'
#' @inheritParams wiener_mle
#' @param x Degradation levels of one trajectory, with strictly positive increments;
#'   or a list of data frames, one per trajectory, whose first column contains
#'   the inspection times and second column the degradation levels.
#'   The increments can also be given directly with [increments_data()].
#' @param t Inspection times of the trajectory, strictly increasing. Not used
#'   when `x` is a list.
#'
#' @return A named vector with `a` \eqn{= \hat a} and `b` \eqn{= \hat b}.
#'
#'   With `ci = TRUE`, a data frame with one row per parameter and the columns
#'   `parameter`, `estimate`, `lower` and `upper` (bounds of the confidence
#'   interval). When the estimator does not exist, the estimates (and bounds)
#'   are `NA`, with a warning.
#'
#' @references
#' van Noortwijk, J. M. (2009). A survey of the application of gamma processes
#' in maintenance. *Reliability Engineering & System Safety*, 94(1), 2-21.
#' \doi{10.1016/j.ress.2007.03.019}
#'
#' @seealso [wiener_mle()], [ig_mle()]
#'
#' @examplesIf requireNamespace("DegradSim", quietly = TRUE)
#' # 1. One trajectory inspected every time unit up to 50 (a = 2, b = 3)
#' t <- seq(0, 50, by = 1)
#' set.seed(1)
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3)
#' gamma_mle(x)
#'
#' # The same trajectory given as two vectors: degradation levels, then times
#' gamma_mle(x[, "X1"], t)
#'
#' # 2. One trajectory with irregular inspections, more and more frequent
#' t <- DegradSim::inspection_times("increasing", n = 50, tmax = 50)$time
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3)
#' gamma_mle(x)
#'
#' # 3. Ten trajectories with the same inspection times
#' t <- seq(0, 50, by = 5)
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3, n_traj = 10)
#' gamma_mle(x)
#'
#' # 4. Ten trajectories, each with its own inspection times
#' plan <- DegradSim::inspection_times("decreasing", n = 5:14, n_traj = 10, tmax = 50)
#' x <- DegradSim::sim_gamma(plan, a = 2, b = 3)
#' gamma_mle(x)
#'
#' # 5. Estimates with 95 % confidence intervals
#' gamma_mle(x, ci = TRUE)
#'
#' # 6. Your own data: a list of data frames, one per trajectory,
#' #    with the times in the first column and the levels in the second
#' liste <- list(data.frame(t = c(0, 10, 20, 30), X = c(0, 6.1, 13.4, 19.8)),
#'               data.frame(t = c(0, 15, 30), X = c(0, 9.7, 20.5)))
#' gamma_mle(liste)
#'
#' # The gamma process only increases: a decrease stops the estimation
#' try(gamma_mle(c(0, 2, 1.5, 3), c(0, 1, 2, 3)))
#' @export
gamma_mle <- function(x, t = NULL, ci = FALSE, level = 0.95) {
  estimate("gamma", x, t, ci, level)
}

gamma_est <- function(x, t = NULL) {
  d <- increments(x, t)
  gamma_fit(d$dt, d$dx)
}

gamma_fit <- function(dt, dx) {
  if (any(dx <= 0)) stop("increments must be strictly positive")
  r <- dx / dt
  if (all(abs(r - r[1]) <= 1e-12 * r[1])) no_mle("the MLE does not exist: all speeds dx/dt are equal")

  Tt <- sum(dt)
  S <- sum(dx)
  U <- sum(dt * log(dx))

  r_hat <- S / Tt
  s2 <- sum((dx - r_hat * dt)^2) / (Tt - sum(dt^2) / Tt)
  a_mm <- r_hat^2 / s2

  g <- function(u) U - sum(dt * digamma(exp(u) * dt)) + Tt * log(exp(u) * Tt / S)
  a <- exp(uniroot(g, c(log(a_mm) - 1, log(a_mm) + 1), extendInt = "downX", tol = 1e-12)$root)

  c(a = a, b = a * Tt / S)
}
