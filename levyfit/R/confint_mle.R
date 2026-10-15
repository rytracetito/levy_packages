#' Confidence intervals for the maximum likelihood estimators
#'
#' Computes the maximum likelihood estimates of a degradation model with the
#' corresponding `_mle` function, and asymptotic confidence intervals for its
#' parameters, from one or several trajectories.
#'
#' @section Method:
#' Let \eqn{\hat\psi} be the MLE and \eqn{\ell} the log-likelihood. The
#' observed Fisher information is the matrix
#' \eqn{J = -\partial^2 \ell / \partial \psi \partial \psi^\top} evaluated at
#' \eqn{\hat\psi}, computed by numerical differentiation ([stats::optimHess()]).
#' The standard errors are the square roots of the diagonal of \eqn{J^{-1}}.
#' With \eqn{z} the quantile of order \eqn{(1 + level) / 2} of the standard
#' normal distribution, the interval of a parameter \eqn{\psi_k} is
#' * \eqn{\hat\psi_k \pm z\, se_k} for a real parameter
#'   (`mu` of the Wiener process, `theta` of the variance-gamma process,
#'   `b` and `m` of the Meixner process);
#' * \eqn{\hat\psi_k \exp(\pm z\, se_k / \hat\psi_k)} for a positive
#'   parameter: this is the interval \eqn{\ln \hat\psi_k \pm z\, se_k / \hat\psi_k}
#'   for \eqn{\ln \psi_k} transformed back, which always stays positive.
#'
#' These intervals are asymptotic: their coverage is close to `level` when the
#' number of increments (or of trajectories, for the inverse gamma process) is
#' large. The unbiased estimates `s2` (Wiener) and `eta_unbiased` (inverse
#' Gaussian) are given without interval: use the interval of `sigma2` or `eta`.
#' When \eqn{J} is not positive definite (estimate close to a boundary of the
#' parameter space), the bounds are `NA`.
#'
#' @inheritParams wiener_mle
#' @param model Model to estimate: `"gamma"`, `"wiener"`, `"ig"` (inverse
#'   Gaussian), `"vg"` (variance-gamma), `"meixner"` or `"igam"` (inverse
#'   gamma).
#' @param level Confidence level, between 0 and 1 (default 0.95).
#'
#' @return A data frame with one row per parameter and four columns:
#'   * `parameter`: name of the parameter, as returned by the `_mle` function;
#'   * `estimate`: maximum likelihood estimate (or unbiased estimate for `s2`
#'     and `eta_unbiased`);
#'   * `lower`, `upper`: bounds of the confidence interval.
#'
#' @seealso [gamma_mle()], [wiener_mle()], [ig_mle()], [vg_mle()],
#'   [meixner_mle()], [igam_mle()]
#'
#' @examples
#' # One trajectory of a gamma process (a = 2, b = 3)
#' t <- inspection_times("regular", n = 50, tmax = 100)$time
#' set.seed(1)
#' x <- c(0, cumsum(rgamma(50, shape = 2 * diff(t), rate = 3)))
#' confint_mle(x, t, model = "gamma")
#'
#' # 99 % intervals
#' confint_mle(x, t, model = "gamma", level = 0.99)
#'
#' # Ten trajectories of a Wiener process (mu = 0.5, sigma2 = 2.25)
#' t <- inspection_times("increasing", n = 20, tmax = 100)$time
#' liste <- lapply(1:10, function(i)
#'   data.frame(t = t, X = c(0, cumsum(0.5 * diff(t) + 1.5 * sqrt(diff(t)) * rnorm(20)))))
#' confint_mle(liste, model = "wiener")
#' @export
confint_mle <- function(x, t = NULL, model, level = 0.95) {
  model <- match.arg(model, c("gamma", "wiener", "ig", "vg", "meixner", "igam"))
  if (level <= 0 || level >= 1) stop("level must be between 0 and 1")

  fit <- switch(model, gamma = gamma_mle(x, t), wiener = wiener_mle(x, t), ig = ig_mle(x, t),
                vg = vg_mle(x, t), meixner = meixner_mle(x, t), igam = igam_mle(x, t))
  par <- switch(model, wiener = fit[c("mu", "sigma2")], ig = fit[c("mu", "eta")], fit)
  real <- names(par) %in% c("mu", "theta", "b", "m") & model != "ig"

  ll <- loglik_fun(x, t, model)
  H <- optimHess(par, function(p) -ll(p), control = list(parscale = pmax(abs(par), 1e-8)))
  V <- tryCatch(solve(H), error = function(e) matrix(NA, length(par), length(par)))
  se <- sqrt(ifelse(diag(V) > 0, diag(V), NA))

  z <- qnorm((1 + level) / 2)
  lower <- ifelse(real, par - z * se, par * exp(-z * se / par))
  upper <- ifelse(real, par + z * se, par * exp(z * se / par))
  out <- data.frame(parameter = names(fit), estimate = unname(fit), lower = NA_real_, upper = NA_real_)
  out$lower[match(names(par), out$parameter)] <- lower
  out$upper[match(names(par), out$parameter)] <- upper
  out
}

loglik_fun <- function(x, t, model) {
  if (model == "igam") {
    traj <- as_trajectories(x, t)
    W <- lapply(traj, function(d) d$x)
    tt <- lapply(traj, function(d) d$t)
    return(function(p) igam_loglik(W, tt, p[1], p[2], p[3]))
  }
  d <- increments(x, t)
  dt <- d$dt
  dx <- d$dx
  switch(model,
    gamma = function(p) sum(dgamma(dx, shape = p[1] * dt, rate = p[2], log = TRUE)),
    wiener = function(p) sum(dnorm(dx, p[1] * dt, sqrt(p[2] * dt), log = TRUE)),
    ig = function(p) {
      mu <- p[1] * dt
      la <- p[2] * dt^2
      sum(0.5 * log(la / (2 * pi * dx^3)) - la * (dx - mu)^2 / (2 * mu^2 * dx))
    },
    vg = function(p) vg_loglik(dx, dt, p[1], p[2], p[3]),
    meixner = function(p) meixner_loglik(dx, dt, p[1], p[2], p[3], p[4]))
}
