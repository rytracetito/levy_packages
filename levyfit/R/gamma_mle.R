#' Maximum likelihood estimation for the gamma process
#'
#' Computes the maximum likelihood estimator (MLE) of the shape \eqn{a} and
#' the rate \eqn{b} of the gamma process from one or several trajectories
#' observed at regular or irregular inspection times. The MLE always exists
#' (with probability one), is unique, and is the root of an equation in one
#' variable.
#'
#' @section Model and data:
#' The gamma process is increasing, with independent increments
#' \eqn{X_t - X_s \sim \Gamma(a (t - s), b)}, with shape \eqn{a(t - s)} and rate
#' \eqn{b} (Kahle, Mercier and Paroissin, 2016, definition 2.24 with
#' \eqn{A(t) = a t}), so that \eqn{E[X_t] = a t / b} and
#' \eqn{Var(X_t) = a t / b^2} (Kahle et al., 2016, section 2.7). The density of
#' an increment over a duration \eqn{\Delta t} is (Kahle et al., 2016,
#' definition 2.1)
#' \deqn{f_{\Delta t}(x; a, b) = \frac{b^{a \Delta t}}{\Gamma(a \Delta t)}\, x^{a \Delta t - 1} e^{-b x}, \qquad x > 0.}
#' The data are the \eqn{n} increments \eqn{(\Delta t_j, \Delta x_j)},
#' \eqn{j = 1, \dots, n}, with \eqn{\Delta x_j > 0}, pooled over the
#' trajectories, with \eqn{T = \sum_{j=1}^{n} \Delta t_j} and
#' \eqn{S = \sum_{j=1}^{n} \Delta x_j} (see the section Notation of
#' [levyfit-package]). Let also
#' \deqn{U = \sum_{j=1}^{n} \Delta t_j \ln \Delta x_j, \qquad
#'   r_j = \frac{\Delta x_j}{\Delta t_j} \ \text{(speeds)}, \qquad
#'   w_j = \frac{\Delta t_j}{T} \ \text{(weights)}.}
#'
#' @section Likelihood:
#' The increments being independent, the log-likelihood is (Kahle et al.,
#' 2016, equation \[2.42\] with \eqn{\Delta A_{i,j}(\theta) = a s_{i,j}})
#' \deqn{\ell(a, b) = \sum_{j=1}^{n} \big[(a \Delta t_j - 1) \ln \Delta x_j - \ln \Gamma(a \Delta t_j)
#'   + a \Delta t_j \ln b - b \Delta x_j\big].}
#' It is a two-parameter exponential family with sufficient statistic
#' \eqn{(U, S)}, whatever the durations \eqn{\Delta t_j}, equal or not. With
#' equal durations \eqn{\Delta t}, the increments are i.i.d.
#' \eqn{\Gamma(a \Delta t, b)} and \eqn{(U, S)} is, up to the factor
#' \eqn{\Delta t}, the usual sufficient statistic
#' \eqn{(\sum_{j=1}^{n} \ln \Delta x_j, \sum_{j=1}^{n} \Delta x_j)} of an
#' i.i.d. gamma sample.
#'
#' @section Maximum likelihood estimator:
#' The score in \eqn{b} is \eqn{\partial \ell / \partial b = a T / b - S}, so
#' that, for a fixed \eqn{a}, \eqn{\ell} is maximal at
#' \eqn{\hat b(a) = a T / S}. The profile log-likelihood
#' \eqn{\ell_p(a) = \ell(a, \hat b(a))} has derivative
#' \deqn{g(a) = U - \sum_{j=1}^{n} \Delta t_j\, \psi(a \Delta t_j) + T \ln \frac{a T}{S},}
#' where \eqn{\psi} is the digamma function. The function
#' \eqn{x \mapsto \ln x - \psi(x)} is strictly completely monotonic on
#' \eqn{(0, \infty)} (Alzer, 1997, Theorem 1 with \eqn{\alpha = 0}), hence
#' strictly decreasing, so that \eqn{\psi'(x) > 1/x} for \eqn{x > 0}; therefore
#' \eqn{g'(a) = -\sum_{j=1}^{n} \Delta t_j^2 \psi'(a \Delta t_j) + T / a < 0}
#' and \eqn{g} is continuous and strictly decreasing on \eqn{(0, \infty)}.
#' Moreover, by the inequality
#' \eqn{1/(2x) < \ln x - \psi(x) < 1/x}, \eqn{x > 0} (Alzer, 1997,
#' equation (2.2)), \eqn{g(a) > g_0 + n / (2a)} with \eqn{g_0} not depending
#' on \eqn{a}, so that \eqn{g(a) \to +\infty} as \eqn{a \to 0}, and
#' \eqn{\psi(x) - \ln x \to 0} as \eqn{x \to \infty}, so that, by Jensen's
#' inequality,
#' \deqn{\lim_{a \to \infty} g(a) = T \Big(\sum_{j=1}^{n} w_j \ln r_j - \ln \sum_{j=1}^{n} w_j r_j\Big) < 0,}
#' unless all the speeds \eqn{r_j} are equal. By the intermediate value
#' theorem, \eqn{g} has a root \eqn{\hat a}, unique since \eqn{g} is strictly
#' decreasing, and \eqn{\hat b = \hat a T / S}. These are the equations
#' \[2.43\] and \[2.44\] of Kahle, Mercier and Paroissin (2016, section 2.7.1.1).
#' With equal durations \eqn{\Delta t}, the equation \eqn{g(a) = 0} becomes
#' equation \[2.45a\] of Kahle et al. (2016), whose common step \eqn{\delta}
#' is \eqn{\Delta t} here:
#' \deqn{\psi(a \Delta t) - \ln(a \Delta t) = \frac{1}{n} \sum_{j=1}^{n} \ln \Delta x_j
#'   - \ln\Big(\frac{1}{n} \sum_{j=1}^{n} \Delta x_j\Big),}
#' the logarithm of the ratio of the geometric and arithmetic means of the
#' increments. If all the speeds are equal (in particular if \eqn{n = 1}), the
#' MLE does not exist and the function returns `NA` with a warning; this has
#' probability zero when \eqn{n \ge 2}.
#'
#' The root is computed by [stats::uniroot()] in \eqn{\ln a}. Kahle, Mercier
#' and Paroissin (2016, remark 2.50) recommend a moment estimator as starting
#' point. The mean and the variance of the gamma process are \eqn{a t / b} and
#' \eqn{a t / b^2} (Kahle et al., 2016, section 2.7), that is
#' \eqn{\mu = a / b} and \eqn{\nu = a / b^2} per unit of time, so that
#' \eqn{a = \mu^2 / \nu}. The search starts from
#' \eqn{a_{MM} = \hat\mu^2 / \hat\nu}, where
#' \deqn{\hat\mu = \frac{\sum_{j=1}^{n} \Delta x_j}{\sum_{j=1}^{n} \Delta t_j} = \frac{S}{T}, \qquad
#'   \hat\nu = \frac{1}{n} \sum_{j=1}^{n} \frac{(\Delta x_j - \hat\mu \Delta t_j)^2}{\Delta t_j}}
#' are the estimators of the mean and of the variance per unit of time of
#' Kahle et al. (2016, section 1.6.1, equation \[1.19\]), given there for the
#' Wiener process. Since the root of \eqn{g} is unique, it does not depend on
#' this starting point, which only makes the computation faster.
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
#' Kahle, W., Mercier, S. and Paroissin, C. (2016). *Degradation Processes in
#' Reliability*. ISTE and Wiley, London. Definitions 2.1 and 2.24, section
#' 1.6.1 (equation \[1.19\]) and section 2.7 (equations \[2.42\] to \[2.45\]
#' and remark 2.50).
#'
#' Alzer, H. (1997). On some inequalities for the gamma and psi functions.
#' *Mathematics of Computation*, 66(217), 373-389.
#' \doi{10.1090/S0025-5718-97-00807-7}. Equation (2.2) and Theorem 1.
#'
#' van Noortwijk, J. M. (2009). A survey of the application of gamma processes
#' in maintenance. *Reliability Engineering & System Safety*, 94(1), 2-21.
#' \doi{10.1016/j.ress.2007.03.019}
#'
#' @seealso [wiener_mle()], [ig_mle()]
#'
#' @examplesIf requireNamespace("DegradSim", quietly = TRUE)
#' # 1. One trajectory inspected every time unit up to 50 (a = 2, b = 3)
#' t <- DegradSim::inspection_times("regular", n = 50, tmax = 50)$time
#' set.seed(1)
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3)
#' gamma_mle(x)
#'
#' # The same trajectory given as two vectors: degradation levels, then times
#' gamma_mle(x[, "X1"], t)
#'
#' # 2. One trajectory with irregular inspections, more and more frequent
#' t <- DegradSim::inspection_times("decreasing", n = 50, tmax = 50)$time
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3)
#' gamma_mle(x)
#'
#' # 3. Ten trajectories with the same inspection times
#' t <- DegradSim::inspection_times("regular", n = 10, tmax = 50)$time
#' x <- DegradSim::sim_gamma(t, a = 2, b = 3, n_traj = 10)
#' gamma_mle(x)
#'
#' # 4. Ten trajectories, each with its own inspection times
#' plan <- DegradSim::inspection_times("increasing", n = 5:14, n_traj = 10, tmax = 50)
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

  mu_hat <- S / Tt
  nu_hat <- mean((dx - mu_hat * dt)^2 / dt)
  a_mm <- mu_hat^2 / nu_hat

  g <- function(u) U - sum(dt * digamma(exp(u) * dt)) + Tt * log(exp(u) * Tt / S)
  a <- exp(uniroot(g, c(log(a_mm) - 1, log(a_mm) + 1), extendInt = "downX", tol = 1e-12)$root)

  c(a = a, b = a * Tt / S)
}
