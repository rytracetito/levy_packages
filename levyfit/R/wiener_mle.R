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
#' @examples
#' t <- inspection_times("increasing", n = 20, tmax = 50)$time
#' set.seed(1)
#' x <- c(0, cumsum(0.5 * diff(t) + 1.5 * sqrt(diff(t)) * rnorm(20)))
#' wiener_mle(x, t)
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
