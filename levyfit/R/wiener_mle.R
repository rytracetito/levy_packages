#' Maximum likelihood estimation for the Wiener process
#'
#' Computes the maximum likelihood estimator (MLE) of the drift \eqn{\mu} and
#' the variance \eqn{\sigma^2} of the Wiener process
#' \eqn{X_t = \mu t + \sigma B_t}, where \eqn{B} is a standard Brownian motion,
#' from one or several trajectories observed at possibly irregular inspection times. The MLE
#' is obtained in closed form by maximizing the likelihood directly, without
#' numerical optimization.
#'
#' @section Model and data:
#' The increments are independent and
#' \eqn{X_t - X_s \sim N(\mu (t - s), \sigma^2 (t - s))} for \eqn{s < t}; they
#' can be negative. The trajectories are independent and have the same law;
#' each one is observed at its own inspection times. Let
#' \eqn{(\Delta t_j, \Delta x_j)}, \eqn{j = 1, \dots, m}, be the durations and
#' degradation increments between two consecutive inspections of a same
#' trajectory, gathered over all the trajectories, \eqn{T = \sum_j \Delta t_j},
#' \eqn{S = \sum_j \Delta x_j} and \eqn{r_j = \Delta x_j / \Delta t_j} (speeds).
#' The likelihood is a product over all these independent increments, so the
#' formulas below hold for one or several trajectories.
#'
#' @section Likelihood:
#' The increments being independent and Gaussian, the log-likelihood is
#' \deqn{\ell(\mu, \sigma^2) = -\frac{m}{2} \ln(2\pi) - \frac{1}{2} \sum_j \ln \Delta t_j
#'   - \frac{m}{2} \ln \sigma^2 - \frac{Q(\mu)}{2 \sigma^2},
#'   \qquad Q(\mu) = \sum_j \frac{(\Delta x_j - \mu \Delta t_j)^2}{\Delta t_j}.}
#'
#' @section Maximum likelihood estimator:
#' 1. Completing the square gives \eqn{Q(\mu) = T (\mu - \hat\mu)^2 + \hat Q},
#'    with \eqn{\hat\mu = S / T} and
#'    \eqn{\hat Q = Q(\hat\mu) = \sum_j \Delta t_j (r_j - \hat\mu)^2}.
#' 2. For every \eqn{\sigma^2}, \eqn{\ell} is therefore maximal in \eqn{\mu} at
#'    \eqn{\hat\mu = S / T}, the speeds averaged with weights \eqn{\Delta t_j / T}.
#' 3. The remaining function \eqn{-\frac{m}{2} \ln \sigma^2 - \hat Q / (2 \sigma^2)}
#'    is maximal at \eqn{\hat\sigma^2 = \hat Q / m}.
#' 4. The MLE exists and is unique if and only if \eqn{\hat Q > 0}, that is if
#'    the speeds \eqn{r_j} are not all equal, which has probability one as soon
#'    as \eqn{m \ge 2}. Otherwise (in particular if \eqn{m = 1}) the likelihood
#'    is unbounded and the function stops.
#'
#' With one trajectory, \eqn{\hat\mu} only depends on the first and last
#' observations: the intermediate inspections only bring information on
#' \eqn{\sigma^2}.
#'
#' @section Bias:
#' Writing the standardized increments as a Gaussian vector and applying an
#' orthogonal rotation shows that, for \eqn{m \ge 2},
#' \deqn{\hat\mu \sim N(\mu, \sigma^2 / T), \qquad \hat Q / \sigma^2 \sim \chi^2_{m-1},}
#' and that \eqn{\hat\mu} and \eqn{\hat Q} are independent. Hence
#' \eqn{\hat\mu} is unbiased, while
#' \eqn{E[\hat\sigma^2] = \frac{m - 1}{m} \sigma^2}: \eqn{\hat\sigma^2} is
#' biased downwards. The function also returns its unbiased version
#' \deqn{s^2 = \frac{\hat Q}{m - 1} = \frac{m}{m - 1} \hat\sigma^2.}
#'
#' @param x Degradation levels of one trajectory; or a list of data frames,
#'   one per trajectory, whose first column contains the inspection times and
#'   second column the degradation levels.
#' @param t Inspection times of the trajectory, strictly increasing. Not used
#'   when `x` is a list.
#'
#' @return A named vector with the maximum likelihood estimates
#'   `mu` \eqn{= \hat\mu} and `sigma2` \eqn{= \hat\sigma^2}, and the unbiased
#'   estimate `s2` \eqn{= s^2} of \eqn{\sigma^2}.
#'
#' @seealso [gamma_mle()], [ig_mle()], [vg_mle()], [meixner_mle()],
#'   [igam_mle()], [inspection_times()]
#'
#' @examplesIf requireNamespace("levysim", quietly = TRUE)
#' # 1. One trajectory inspected every time unit up to 50 (mu = 0.5, sigma2 = 2.25)
#' t <- seq(0, 50, by = 1)
#' set.seed(1)
#' x <- levysim::sim_wiener(t, mu = 0.5, sigma = 1.5)
#' wiener_mle(x)
#'
#' # The same trajectory given as two vectors: degradation levels, then times
#' wiener_mle(x[[1]]$X, t)
#'
#' # 2. One trajectory with irregular inspections, more and more frequent
#' t <- inspection_times("increasing", n = 50, tmax = 50)$time
#' x <- levysim::sim_wiener(t, mu = 0.5, sigma = 1.5)
#' wiener_mle(x)
#'
#' # 3. Ten trajectories with the same inspection times
#' t <- seq(0, 50, by = 5)
#' x <- levysim::sim_wiener(t, mu = 0.5, sigma = 1.5, n_traj = 10)
#' wiener_mle(x)
#'
#' # 4. Ten trajectories, each with its own inspection times
#' plan <- inspection_times("decreasing", n = 5:14, n_traj = 10, tmax = 50)
#' x <- levysim::sim_wiener(plan, mu = 0.5, sigma = 1.5)
#' wiener_mle(x)
#'
#' # 5. Estimates with 95 % confidence intervals
#' confint_mle(x, model = "wiener")
#'
#' # 6. Your own data: a list of data frames, one per trajectory,
#' #    with the times in the first column and the levels in the second
#' liste <- list(data.frame(t = c(0, 10, 20, 30), X = c(0, 4.2, 11.5, 13.9)),
#'               data.frame(t = c(0, 15, 30), X = c(0, 9.1, 12.8)))
#' wiener_mle(liste)
#'
#' # sigma2 is the MLE and s2 its unbiased version: with few increments,
#' # here 3, they differ clearly
#' wiener_mle(c(0, 1.2, 0.7, 2.5), c(0, 1, 2, 3))
#' @export
wiener_mle <- function(x, t = NULL) {
  d <- increments(x, t)
  dt <- d$dt
  dx <- d$dx

  m <- length(dx)
  mu <- sum(dx) / sum(dt)
  Q <- sum((dx - mu * dt)^2 / dt)
  if (m < 2 || Q <= 1e-12 * sum(dx^2 / dt)) stop("the MLE does not exist: all speeds dx/dt are equal")

  c(mu = mu, sigma2 = Q / m, s2 = Q / (m - 1))
}
