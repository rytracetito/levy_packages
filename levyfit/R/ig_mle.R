#' Estimation of the parameters of the inverse Gaussian process
#'
#' Estimates the parameters \eqn{\mu} and \eqn{\eta} of the homogeneous
#' inverse Gaussian (IG) process from one trajectory observed at possibly
#' irregular inspection times. Both estimators are explicit and unbiased.
#'
#' @section Model and data:
#' \eqn{IG(m, \lambda)} denotes the inverse Gaussian distribution with mean
#' \eqn{m} and shape \eqn{\lambda}, with density
#' \deqn{\sqrt{\frac{\lambda}{2 \pi x^3}} \exp\Big(-\frac{\lambda (x - m)^2}{2 m^2 x}\Big), \qquad x > 0.}
#' The IG process is increasing, with independent increments
#' \eqn{X_t - X_s \sim IG(\mu (t - s), \eta (t - s)^2)}, so that
#' \eqn{E[X_t] = \mu t} and \eqn{Var(X_t) = \mu^3 t / \eta}. From the
#' observations \eqn{(t_j, x_j)}, \eqn{j = 0, \dots, m}, let
#' \eqn{\Delta t_j = t_j - t_{j-1}}, \eqn{\Delta x_j = x_j - x_{j-1} > 0},
#' \eqn{T = \sum_j \Delta t_j} and \eqn{S = \sum_j \Delta x_j}.
#'
#' @section Likelihood:
#' With \eqn{\delta = 1 / \mu}, the log-likelihood is
#' \deqn{\ell(\mu, \eta) = c_0 + \frac{m}{2} \ln \eta
#'   - \frac{\eta}{2} \sum_j \frac{(\delta \Delta x_j - \Delta t_j)^2}{\Delta x_j},
#'   \qquad c_0 = \sum_j \Big[\ln \Delta t_j - \frac{1}{2} \ln(2\pi) - \frac{3}{2} \ln \Delta x_j\Big].}
#'
#' @section Maximum likelihood estimator:
#' 1. Expanding the square gives
#'    \deqn{\sum_j \frac{(\delta \Delta x_j - \Delta t_j)^2}{\Delta x_j} = \frac{(\delta S - T)^2}{S} + \hat Q,
#'      \qquad \hat Q = \sum_j \frac{\Delta t_j^2}{\Delta x_j} - \frac{T^2}{S} \ge 0,}
#'    with \eqn{\hat Q = 0} if and only if the speeds
#'    \eqn{\Delta x_j / \Delta t_j} are all equal.
#' 2. For every \eqn{\eta}, \eqn{\ell} is therefore maximal at
#'    \eqn{\delta = T / S}, that is \eqn{\hat\mu = S / T}.
#' 3. The remaining function \eqn{\frac{m}{2} \ln \eta - \frac{\eta}{2} \hat Q}
#'    is maximal at \eqn{\hat\eta = m / \hat Q}.
#' 4. The MLE exists and is unique if and only if \eqn{\hat Q > 0}, which has
#'    probability one as soon as \eqn{m \ge 2}. Otherwise the likelihood is
#'    unbounded and the function stops.
#'
#' @section Unbiased estimators:
#' For \eqn{m \ge 2},
#' \deqn{\hat\mu \sim IG(\mu, \eta T), \qquad \eta \hat Q \sim \chi^2_{m-1},}
#' and \eqn{\hat Q} is independent of \eqn{S}: the conditional density of the
#' increments given \eqn{S} does not depend on \eqn{\mu}, and its exponent only
#' involves \eqn{\hat Q} (this extends a result of Tweedie (1957) to irregular
#' inspection times). Hence \eqn{\hat\mu} is unbiased, while, since
#' \eqn{E[1/V] = 1/(k - 2)} for \eqn{V \sim \chi^2_k},
#' \eqn{E[\hat\eta] = \frac{m}{m - 3} \eta} for \eqn{m > 3}: the maximum
#' likelihood estimator of \eqn{\eta} is biased upwards. The function returns
#' its unbiased version
#' \deqn{\tilde\eta = \frac{m - 3}{\hat Q},}
#' which requires at least four increments.
#'
#' @param x Degradation levels, with strictly positive increments.
#' @param t Inspection times, strictly increasing.
#'
#' @return A named vector with `mu` \eqn{= \hat\mu} and
#'   `eta` \eqn{= \tilde\eta}.
#'
#' @references
#' Tweedie, M. C. K. (1957). Statistical properties of inverse Gaussian
#' distributions. I. *The Annals of Mathematical Statistics*, 28(2), 362-377.
#' \doi{10.1214/aoms/1177706964}
#'
#' Wang, X. and Xu, D. (2010). An inverse Gaussian process model for
#' degradation data. *Technometrics*, 52(2), 188-197.
#' \doi{10.1198/TECH.2009.08197}
#'
#' Ye, Z.-S. and Chen, N. (2014). The inverse Gaussian process as a
#' degradation model. *Technometrics*, 56(3), 302-311.
#' \doi{10.1080/00401706.2013.830074}
#'
#' @seealso [wiener_mle()], [gamma_mle()], [inspection_times()]
#'
#' @examples
#' t <- inspection_times("inc-dec", n = 30, tmax = 50)$time
#' rig <- function(n, mean, shape) {
#'   nu <- rnorm(n)^2
#'   y <- mean + mean^2 * nu / (2 * shape) -
#'     mean / (2 * shape) * sqrt(4 * mean * shape * nu + mean^2 * nu^2)
#'   ifelse(runif(n) <= mean / (mean + y), y, mean^2 / y)
#' }
#' set.seed(1)
#' x <- c(0, cumsum(rig(30, 0.5 * diff(t), 2 * diff(t)^2)))
#' ig_mle(x, t)
#' @export
ig_mle <- function(x, t) {
  dt <- diff(t)
  dx <- diff(x)
  if (any(dt <= 0)) stop("inspection times must be strictly increasing")
  if (any(dx <= 0)) stop("increments must be strictly positive")

  m <- length(dx)
  if (m < 4) stop("at least four increments are needed")
  Tt <- sum(dt)
  S <- sum(dx)
  Q <- sum(dt^2 / dx) - Tt^2 / S
  if (Q <= 1e-12 * sum(dt^2 / dx)) stop("the MLE does not exist: all speeds dx/dt are equal")

  c(mu = S / Tt, eta = (m - 3) / Q)
}
