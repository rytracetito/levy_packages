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
#'    if \eqn{m = 1}), the MLE does not exist and the function stops; this has
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
#' @param x Degradation levels of one trajectory, with strictly positive increments;
#'   or a list of data frames, one per trajectory, whose first column contains
#'   the inspection times and second column the degradation levels.
#' @param t Inspection times of the trajectory, strictly increasing. Not used
#'   when `x` is a list.
#'
#' @return A named vector with `a` \eqn{= \hat a} and `b` \eqn{= \hat b}.
#'
#' @references
#' van Noortwijk, J. M. (2009). A survey of the application of gamma processes
#' in maintenance. *Reliability Engineering & System Safety*, 94(1), 2-21.
#' \doi{10.1016/j.ress.2007.03.019}
#'
#' @seealso [wiener_mle()], [ig_mle()], [inspection_times()]
#'
#' @examples
#' t <- inspection_times("decreasing", n = 30, tmax = 50)$time
#' set.seed(1)
#' x <- c(0, cumsum(rgamma(30, shape = 2 * diff(t), rate = 3)))
#' gamma_mle(x, t)
#' @export
gamma_mle <- function(x, t = NULL) {
  d <- increments(x, t)
  dt <- d$dt
  dx <- d$dx
  if (any(dx <= 0)) stop("increments must be strictly positive")
  r <- dx / dt
  if (all(abs(r - r[1]) <= 1e-12 * r[1])) stop("the MLE does not exist: all speeds dx/dt are equal")

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
