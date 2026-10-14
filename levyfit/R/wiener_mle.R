#' Maximum likelihood estimation for the Wiener process
#'
#' Computes the maximum likelihood estimator (MLE) of the drift \eqn{\mu} and
#' the diffusion coefficient \eqn{\sigma} of the Wiener process
#' \eqn{X_t = \mu t + \sigma B_t}, where \eqn{B} is a standard Brownian motion,
#' from one trajectory observed at possibly irregular inspection times. The
#' estimator is explicit and its exact distribution is known, which gives exact
#' confidence intervals.
#'
#' @section Model and data:
#' The increments are independent and
#' \eqn{X_t - X_s \sim N(\mu (t - s), \sigma^2 (t - s))} for \eqn{s < t}; they
#' can be negative. From the observations \eqn{(t_j, x_j)}, \eqn{j = 0, \dots, m},
#' let \eqn{\Delta t_j = t_j - t_{j-1}}, \eqn{\Delta x_j = x_j - x_{j-1}},
#' \eqn{T = \sum_j \Delta t_j}, \eqn{S = \sum_j \Delta x_j} and
#' \eqn{r_j = \Delta x_j / \Delta t_j} (speeds).
#'
#' @section Likelihood:
#' The increments being independent and Gaussian, the log-likelihood is
#' \deqn{\ell(\mu, \sigma^2) = -\frac{m}{2} \ln(2\pi) - \frac{1}{2} \sum_j \ln \Delta t_j
#'   - \frac{m}{2} \ln \sigma^2 - \frac{Q(\mu)}{2 \sigma^2},
#'   \qquad Q(\mu) = \sum_j \frac{(\Delta x_j - \mu \Delta t_j)^2}{\Delta t_j}.}
#'
#' @section Derivation of the estimator:
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
#' @section Exact distributions and confidence intervals:
#' Writing the standardized increments as a Gaussian vector and applying an
#' orthogonal rotation shows that, for \eqn{m \ge 2},
#' \deqn{\hat\mu \sim N(\mu, \sigma^2 / T), \qquad \hat Q / \sigma^2 \sim \chi^2_{m-1},}
#' and that \eqn{\hat\mu} and \eqn{\hat\sigma^2} are independent. Consequently
#' \eqn{s^2 = \hat Q / (m - 1)} is unbiased,
#' \eqn{(\hat\mu - \mu) / \sqrt{s^2 / T}} follows a Student distribution with
#' \eqn{m - 1} degrees of freedom, and the intervals
#' \deqn{\hat\mu \pm t_{m-1, 1-\alpha/2} \sqrt{s^2 / T}, \qquad
#'   \Big[\frac{\hat Q}{\chi^2_{m-1, 1-\alpha/2}}, \frac{\hat Q}{\chi^2_{m-1, \alpha/2}}\Big]}
#' have exact level \eqn{1 - \alpha} for \eqn{\mu} and \eqn{\sigma^2}; the
#' interval for \eqn{\sigma} is obtained by taking square roots.
#'
#' @section Standard errors:
#' The Fisher information is
#' \eqn{I(\mu, \sigma^2) = \mathrm{diag}(T / \sigma^2, m / (2 \sigma^4))}, hence
#' \eqn{se(\hat\mu) = \hat\sigma / \sqrt{T}},
#' \eqn{se(\hat\sigma^2) = \hat\sigma^2 \sqrt{2 / m}} and, by the delta method,
#' \eqn{se(\hat\sigma) = \hat\sigma / \sqrt{2m}}. `vcov` is
#' \eqn{I^{-1}} at the estimate. With few increments, prefer the exact
#' intervals to Wald intervals.
#'
#' @param x Degradation levels.
#' @param t Inspection times, strictly increasing.
#' @param level Confidence level \eqn{1 - \alpha} of the exact intervals.
#'
#' @return A list with the estimates `mu`, `sigma`, `sigma2`, the unbiased
#'   estimate `s2` of \eqn{\sigma^2}, the standard errors `se`, the covariance
#'   matrix `vcov` of \eqn{(\hat\mu, \hat\sigma^2)}, the maximum log-likelihood
#'   `loglik`, the exact confidence intervals `ci_mu`, `ci_sigma2`, `ci_sigma`
#'   and the number of increments `n_incr`.
#'
#' @seealso [gamma_mle()], [ig_mle()], [inspection_times()]
#'
#' @examples
#' t <- inspection_times("increasing", n = 20, tmax = 50)$time
#' set.seed(1)
#' x <- c(0, cumsum(0.5 * diff(t) + 1.5 * sqrt(diff(t)) * rnorm(20)))
#' fit <- wiener_mle(x, t)
#' c(mu = fit$mu, sigma = fit$sigma)
#' fit$ci_mu
#' fit$ci_sigma
#' @export
wiener_mle <- function(x, t, level = 0.95) {
  dt <- diff(t)
  dx <- diff(x)
  if (any(dt <= 0)) stop("inspection times must be strictly increasing")

  m <- length(dx)
  Tt <- sum(dt)
  S <- sum(dx)
  mu <- S / Tt
  Q <- sum((dx - mu * dt)^2 / dt)
  if (m < 2 || Q <= 1e-12 * sum(dx^2 / dt)) stop("the MLE does not exist: all speeds dx/dt are equal")

  sigma2 <- Q / m
  s2 <- Q / (m - 1)
  loglik <- sum(dnorm(dx, mean = mu * dt, sd = sqrt(sigma2 * dt), log = TRUE))

  se <- c(mu = sqrt(sigma2 / Tt), sigma2 = sigma2 * sqrt(2 / m), sigma = sqrt(sigma2 / (2 * m)))
  vcov <- diag(c(sigma2 / Tt, 2 * sigma2^2 / m))
  dimnames(vcov) <- list(c("mu", "sigma2"), c("mu", "sigma2"))

  alpha <- 1 - level
  ci_mu <- mu + c(-1, 1) * qt(1 - alpha / 2, m - 1) * sqrt(s2 / Tt)
  ci_sigma2 <- Q / qchisq(c(1 - alpha / 2, alpha / 2), m - 1)

  list(mu = mu, sigma = sqrt(sigma2), sigma2 = sigma2, s2 = s2, se = se, vcov = vcov, loglik = loglik,
       ci_mu = ci_mu, ci_sigma2 = ci_sigma2, ci_sigma = sqrt(ci_sigma2), n_incr = m)
}
