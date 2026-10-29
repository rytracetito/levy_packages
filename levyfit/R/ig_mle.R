#' Maximum likelihood estimation for the inverse Gaussian process
#'
#' Computes the maximum likelihood estimator (MLE) of the parameters \eqn{\mu}
#' and \eqn{\eta} of the inverse Gaussian (IG) process from one or several
#' trajectories observed at regular or irregular inspection times. The MLE is
#' obtained in closed form,
#' without numerical optimization.
#'
#' @section Model and data:
#' \eqn{IG(m, \lambda)} denotes the inverse Gaussian distribution with mean
#' \eqn{m} and shape \eqn{\lambda}. The IG process is increasing, with
#' independent increments \eqn{X_t - X_s \sim IG(\mu (t - s), \eta (t - s)^2)},
#' so that \eqn{E[X_t] = \mu t} and \eqn{Var(X_t) = \mu^3 t / \eta}. The
#' density of an increment over a duration \eqn{\Delta t} is
#' \deqn{f_{\Delta t}(x; \mu, \eta) = \Delta t \sqrt{\frac{\eta}{2 \pi x^3}}
#'   \exp\Big(-\frac{\eta (x - \mu \Delta t)^2}{2 \mu^2 x}\Big), \qquad x > 0.}
#' This is the inverse Gaussian process of Wang and Xu (2010, section 1 and
#' equation (1)) with mean function \eqn{\Lambda(t) = \mu t}: they write the
#' increments \eqn{IG(\Delta\Lambda, \eta_{WX} \Delta\Lambda^2)}, with
#' \eqn{\Delta\Lambda = \Lambda(t) - \Lambda(s)}, so that the parameter
#' \eqn{\eta} used here is their \eqn{\eta_{WX}} multiplied by \eqn{\mu^2}:
#' \eqn{\eta = \eta_{WX} \mu^2}.
#' The data are the \eqn{n} increments \eqn{(\Delta t_j, \Delta x_j)},
#' \eqn{j = 1, \dots, n}, with \eqn{\Delta x_j > 0}, pooled over the
#' trajectories, with \eqn{T = \sum_{j=1}^{n} \Delta t_j} and
#' \eqn{S = \sum_{j=1}^{n} \Delta x_j} (see the section Notation of
#' [levyfit-package]).
#'
#' @section Likelihood:
#' The increments being independent (Wang and Xu, 2010, section 2,
#' equation (3)), with \eqn{\delta = 1 / \mu}, the log-likelihood is
#' \deqn{\ell(\mu, \eta) = c_0 + \frac{n}{2} \ln \eta
#'   - \frac{\eta}{2} \sum_{j=1}^{n} \frac{(\delta \Delta x_j - \Delta t_j)^2}{\Delta x_j},
#'   \qquad c_0 = \sum_{j=1}^{n} \Big[\ln \Delta t_j - \frac{1}{2} \ln(2\pi) - \frac{3}{2} \ln \Delta x_j\Big].}
#'
#' @section Maximum likelihood estimator:
#' Expanding the square gives
#' \deqn{\sum_{j=1}^{n} \frac{(\delta \Delta x_j - \Delta t_j)^2}{\Delta x_j} = \frac{(\delta S - T)^2}{S} + \hat Q,
#'   \qquad \hat Q = \sum_{j=1}^{n} \frac{\Delta t_j^2}{\Delta x_j} - \frac{T^2}{S} \ge 0,}
#' with \eqn{\hat Q = 0} if and only if the speeds
#' \eqn{\Delta x_j / \Delta t_j} are all equal. For every \eqn{\eta},
#' \eqn{\ell} is therefore maximal at \eqn{\delta = T / S}, that is
#' \eqn{\hat\mu = S / T}. It remains to maximize
#' \eqn{\frac{n}{2} \ln \eta - \frac{\eta}{2} \hat Q}, which gives
#' \eqn{\hat\eta = n / \hat Q}. The MLE exists and is unique if and only if
#' \eqn{\hat Q > 0}, which has probability one as soon as \eqn{n \ge 2}.
#' Otherwise the likelihood is unbounded and the function returns `NA` with a
#' warning.
#'
#' These are the estimators of Tweedie (1957, section 3, equations (13) and
#' (14)), who considers independent observations \eqn{x_j \sim IG(\mu, \lambda_0 w_j)}
#' with known weights \eqn{w_j}. A change of scale multiplies the mean and the
#' shape of an inverse Gaussian variable by the same factor (Tweedie, 1957,
#' section 1, p. 364), so that the speeds \eqn{x_j = \Delta x_j / \Delta t_j}
#' follow \eqn{IG(\mu, \eta \Delta t_j)}: this is his setting with
#' \eqn{w_j = \Delta t_j} and \eqn{\lambda_0 = \eta}. His estimators are then
#' \eqn{\hat\mu = \sum_{j=1}^{n} w_j x_j / \sum_{j=1}^{n} w_j = S / T} and
#' \eqn{1 / \hat\lambda_0 = \frac{1}{n} \sum_{j=1}^{n} w_j (1 / x_j - 1 / \hat\mu) = \hat Q / n}.
#'
#' @section Bias:
#' With the same correspondence, Tweedie (1957, equation (15) and section 4,
#' equations (21) and (22), p. 367-369) shows that, for \eqn{n \ge 2},
#' \deqn{\hat\mu \sim IG(\mu, \eta T), \qquad \eta \hat Q \sim \chi^2_{n-1},}
#' and that \eqn{\hat\mu} and \eqn{\hat Q} are independent. Hence
#' \eqn{\hat\mu} is unbiased and, since \eqn{\hat\eta = n \eta / V} with
#' \eqn{V = \eta \hat Q \sim \chi^2_{n-1}} and
#' \eqn{E[1/V] = 1/(k - 2)} for \eqn{V \sim \chi^2_k},
#' \eqn{E[\hat\eta] = \frac{n}{n - 3} \eta} for \eqn{n > 3}: \eqn{\hat\eta} is
#' biased upwards. The function also returns its unbiased version
#' \deqn{\tilde\eta = \frac{n - 3}{\hat Q} = \frac{n - 3}{n} \hat\eta,}
#' which requires at least four increments (`NA` otherwise). With equal
#' durations \eqn{\Delta t}, the speeds are i.i.d. \eqn{IG(\mu, \eta \Delta t)}
#' and \eqn{\tilde\eta \Delta t} is the uniformly minimum variance unbiased
#' estimator \eqn{(n - 3) / V} of their shape parameter given by Seshadri
#' (1999, section 2.1, Table 2.1, \eqn{n \ge 4}).
#'
#' @inheritParams wiener_mle
#' @param x Degradation levels of one trajectory, with strictly positive increments;
#'   or a list of data frames, one per trajectory, whose first column contains
#'   the inspection times and second column the degradation levels.
#'   The increments can also be given directly with [increments_data()].
#' @param t Inspection times of the trajectory, strictly increasing. Not used
#'   when `x` is a list.
#'
#' @return A named vector with the maximum likelihood estimates
#'   `mu` \eqn{= \hat\mu} and `eta` \eqn{= \hat\eta}, and the unbiased
#'   estimate `eta_unbiased` \eqn{= \tilde\eta} of \eqn{\eta}.
#'
#'   With `ci = TRUE`, a data frame with one row per parameter and the columns
#'   `parameter`, `estimate`, `lower` and `upper` (bounds of the confidence
#'   interval). When the estimator does not exist, the estimates (and bounds)
#'   are `NA`, with a warning.
#'
#' @references
#' Tweedie, M. C. K. (1957). Statistical properties of inverse Gaussian
#' distributions. I. *The Annals of Mathematical Statistics*, 28(2), 362-377.
#' \doi{10.1214/aoms/1177706964}. Sections 1, 3 and 4, equations (13) to (15),
#' (21) and (22).
#'
#' Seshadri, V. (1999). *The Inverse Gaussian Distribution: Statistical Theory
#' and Applications*. Lecture Notes in Statistics, 137. Springer, New York.
#' \doi{10.1007/978-1-4612-1456-4}. Section 2.1, Table 2.1.
#'
#' Wang, X. and Xu, D. (2010). An inverse Gaussian process model for
#' degradation data. *Technometrics*, 52(2), 188-197.
#' \doi{10.1198/TECH.2009.08197}. Sections 1 and 2, equations (1) and (3).
#'
#' Ye, Z.-S. and Chen, N. (2014). The inverse Gaussian process as a
#' degradation model. *Technometrics*, 56(3), 302-311.
#' \doi{10.1080/00401706.2013.830074}
#'
#' @seealso [wiener_mle()], [gamma_mle()]
#'
#' @examplesIf requireNamespace("DegradSim", quietly = TRUE)
#' # 1. One trajectory inspected every time unit up to 50 (mu = 0.5, eta = 2)
#' t <- DegradSim::inspection_times("regular", n = 50, tmax = 50)$time
#' set.seed(1)
#' x <- DegradSim::sim_ig(t, mu = 0.5, eta = 2)
#' ig_mle(x)
#'
#' # The same trajectory given as two vectors: degradation levels, then times
#' ig_mle(x[, "X1"], t)
#'
#' # 2. One trajectory with irregular inspections, more and more frequent
#' t <- DegradSim::inspection_times("decreasing", n = 50, tmax = 50)$time
#' x <- DegradSim::sim_ig(t, mu = 0.5, eta = 2)
#' ig_mle(x)
#'
#' # 3. Ten trajectories with the same inspection times
#' t <- DegradSim::inspection_times("regular", n = 10, tmax = 50)$time
#' x <- DegradSim::sim_ig(t, mu = 0.5, eta = 2, n_traj = 10)
#' ig_mle(x)
#'
#' # 4. Ten trajectories, each with its own inspection times
#' plan <- DegradSim::inspection_times("increasing", n = 5:14, n_traj = 10, tmax = 50)
#' x <- DegradSim::sim_ig(plan, mu = 0.5, eta = 2)
#' ig_mle(x)
#'
#' # 5. Estimates with 95 % confidence intervals
#' ig_mle(x, ci = TRUE)
#'
#' # 6. Your own data: a list of data frames, one per trajectory,
#' #    with the times in the first column and the levels in the second
#' liste <- list(data.frame(t = c(0, 10, 20, 30), X = c(0, 5.2, 9.8, 15.3)),
#'               data.frame(t = c(0, 15, 30), X = c(0, 7.4, 14.6)))
#' ig_mle(liste)
#'
#' # The inverse Gaussian process only increases: a decrease stops the estimation
#' try(ig_mle(c(0, 2, 1.5, 3), c(0, 1, 2, 3)))
#' @export
ig_mle <- function(x, t = NULL, ci = FALSE, level = 0.95) {
  estimate("ig", x, t, ci, level)
}

ig_est <- function(x, t = NULL) {
  d <- increments(x, t)
  dt <- d$dt
  dx <- d$dx
  if (any(dx <= 0)) stop("increments must be strictly positive")

  m <- length(dx)
  Tt <- sum(dt)
  S <- sum(dx)
  Q <- sum(dt^2 / dx) - Tt^2 / S
  if (m < 2 || Q <= 1e-12 * sum(dt^2 / dx)) no_mle("the MLE does not exist: all speeds dx/dt are equal")

  c(mu = S / Tt, eta = m / Q, eta_unbiased = if (m > 3) (m - 3) / Q else NA)
}
