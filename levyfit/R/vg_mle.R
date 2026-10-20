#' Maximum likelihood estimation for the variance-gamma process
#'
#' Computes the maximum likelihood estimator (MLE) of the parameters
#' \eqn{(\theta, \sigma, \nu)} of the variance-gamma (VG) process from one or several
#' trajectories observed at possibly irregular inspection times. There is no
#' closed form: the MLE is computed by an EM algorithm, whose steps are explicit
#' or reduce to an equation with a unique root, followed by a quasi-Newton
#' refinement.
#'
#' @section Model and data:
#' Let \eqn{\gamma} be a gamma process such that \eqn{\gamma_s} has mean
#' \eqn{s} and variance \eqn{\nu s} (shape \eqn{s / \nu}, scale \eqn{\nu}), and
#' \eqn{W} an independent standard Brownian motion. The VG process is
#' \eqn{X_s = \theta \gamma_s + \sigma W_{\gamma_s}}, a Lévy process that can
#' decrease. With \eqn{\lambda_s = s / \nu},
#' \eqn{c = \sqrt{\theta^2 + 2 \sigma^2 / \nu}} and \eqn{\kappa = c / \sigma^2},
#' \eqn{X_s} has the density
#' \deqn{f_s(x) = \frac{2 e^{\theta x / \sigma^2}}{\sqrt{2\pi}\, \sigma\, \nu^{\lambda_s} \Gamma(\lambda_s)}
#'   \Big(\frac{|x|}{c}\Big)^{\lambda_s - 1/2} K_{|\lambda_s - 1/2|}(\kappa |x|), \qquad x \ne 0,}
#' where \eqn{K_\alpha} is the modified Bessel function of the second kind. The
#' trajectories are independent and have the same law; each one is observed at
#' its own inspection times. Let \eqn{(\Delta t_j, \Delta x_j)},
#' \eqn{j = 1, \dots, m}, be the durations and degradation increments
#' (\eqn{\Delta x_j \ne 0}) between two consecutive inspections of a same
#' trajectory, gathered over all the trajectories.
#'
#' @section Likelihood:
#' The increments being independent,
#' \eqn{\ell(\theta, \sigma, \nu) = \sum_j \ln f_{\Delta t_j}(\Delta x_j)}.
#' The derivatives of \eqn{K_\alpha} with respect to its order have no closed
#' form, so the score equations cannot be solved explicitly.
#'
#' @section Boundary behaviour:
#' As \eqn{\nu \to 0}, \eqn{\ell(\theta, \sigma, \nu)} tends to the
#' log-likelihood of the Wiener process with drift \eqn{\theta} and diffusion
#' \eqn{\sigma}; as \eqn{\nu \to \infty}, \eqn{\ell \to -\infty}. The supremum
#' of the likelihood can therefore be reached at the Brownian limit
#' \eqn{\nu = 0}: the MLE does not always exist (about one case in four with a
#' single trajectory of ten increments). The function then returns `NA` with a
#' warning and the Wiener model, see [wiener_mle()], should be used.
#'
#' As \eqn{\sigma \to 0}, \eqn{X_s \to \theta \gamma_s}: if \eqn{\theta > 0},
#' the VG process tends to the gamma process with shape \eqn{a = 1 / \nu} and
#' rate \eqn{b = 1 / (\theta \nu)}, and when all the increments are positive
#' \eqn{\ell} tends to the gamma log-likelihood. Hence
#' \eqn{\sup \ell \ge \ell_\Gamma(\hat a, \hat b)}, the maximum of the gamma
#' log-likelihood. When all the increments are positive and the VG estimate
#' does not exceed this value, the supremum is reached at the gamma limit
#' \eqn{\sigma = 0}: the function returns `NA` with a warning and the gamma
#' model, see [gamma_mle()], should be used. This typically happens with monotone
#' degradation data.
#'
#' @section EM algorithm:
#' Let \eqn{G_j} be the increment of the hidden gamma time on the \eqn{j}-th
#' interval. Given \eqn{\Delta x_j}, \eqn{G_j} follows a generalized inverse
#' Gaussian distribution \eqn{GIG(p_j, \chi_j, \psi)} with
#' \eqn{p_j = \Delta t_j / \nu - 1/2}, \eqn{\chi_j = \Delta x_j^2 / \sigma^2} and
#' \eqn{\psi = \theta^2 / \sigma^2 + 2 / \nu}, whose moments
#' \eqn{e^+_j = E[G_j]}, \eqn{e^-_j = E[1/G_j]} and \eqn{l_j = E[\ln G_j]} are
#' ratios of Bessel functions.
#' 1. E step: compute \eqn{e^+_j}, \eqn{e^-_j}, \eqn{l_j} at the current
#'    parameters.
#' 2. M step for \eqn{(\theta, \sigma^2)}, explicit and unique: with
#'    \eqn{S = \sum_j \Delta x_j} and \eqn{E^+ = \sum_j e^+_j},
#'    \deqn{\theta = \frac{S}{E^+}, \qquad
#'      \sigma^2 = \frac{1}{m} \Big(\sum_j \Delta x_j^2 e^-_j - \frac{S^2}{E^+}\Big).}
#' 3. M step for \eqn{\nu}: with \eqn{a = 1 / \nu}, the function
#'    \deqn{q(a) = \sum_j \big[\Delta t_j l_j - e^+_j + \Delta t_j \ln a + \Delta t_j - \Delta t_j \psi(a \Delta t_j)\big]}
#'    (\eqn{\psi} the digamma function) is strictly decreasing, from
#'    \eqn{+\infty} to a negative limit: its unique root gives \eqn{\nu}.
#' 4. Repeat until the relative change of \eqn{\ell} is below \eqn{10^{-10}}.
#'
#' Each iteration increases the likelihood (Dempster, Laird and Rubin, 1977).
#' The EM limit is then refined by L-BFGS-B on
#' \eqn{(\theta, \ln \sigma, \ln \nu)}; if this refinement fails or does not
#' increase the likelihood (for instance with data on a very small scale), the
#' EM limit is returned. Starting values: the Wiener MLE for
#' \eqn{(\theta, \sigma)} and \eqn{\nu = T / m}. The algorithm stops at the
#' Brownian limit when \eqn{\nu < \max_j \Delta t_j / 300}.
#'
#' @inheritSection wiener_mle Confidence intervals and estimation failures
#'
#' @inheritParams wiener_mle
#' @param x Degradation levels of one trajectory, with non-zero increments;
#'   or a list of data frames, one per trajectory, whose first column contains
#'   the inspection times and second column the degradation levels.
#' @param t Inspection times of the trajectory, strictly increasing. Not used
#'   when `x` is a list.
#'
#' @return A named vector with the maximum likelihood estimates `theta`,
#'   `sigma` and `nu`.
#'
#'   With `ci = TRUE`, a data frame with one row per parameter and the columns
#'   `parameter`, `estimate`, `lower` and `upper` (bounds of the confidence
#'   interval). When the estimator does not exist, the estimates (and bounds)
#'   are `NA`, with a warning.
#'
#' @references
#' Madan, D. B., Carr, P. P. and Chang, E. C. (1998). The variance gamma
#' process and option pricing. *European Finance Review*, 2(1), 79-105.
#' \doi{10.1023/A:1009703431535}
#'
#' Belhaj Salem, M., Fouladirad, M. and Deloux, E. (2022). Variance Gamma
#' process as degradation model for prognosis and imperfect maintenance of
#' centrifugal pumps. *Reliability Engineering & System Safety*, 223, 108417.
#' \doi{10.1016/j.ress.2022.108417}
#'
#' Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977). Maximum likelihood
#' from incomplete data via the EM algorithm. *Journal of the Royal Statistical
#' Society, Series B*, 39(1), 1-22. \doi{10.1111/j.2517-6161.1977.tb01600.x}
#'
#' @seealso [wiener_mle()], [meixner_mle()], [inspection_times()]
#'
#' @examplesIf requireNamespace("levysim", quietly = TRUE)
#' # 1. One trajectory inspected every 2 time units up to 400
#' #    (theta = 0.5, sigma = 1.5, nu = 2)
#' t <- seq(0, 400, by = 2)
#' set.seed(1)
#' x <- levysim::sim_vg(t, theta = 0.5, sigma = 1.5, nu = 2)
#' vg_mle(x)
#'
#' # The same trajectory given as two vectors: degradation levels, then times
#' vg_mle(x[, "X1"], t)
#'
#' # 2. Twenty trajectories with the same inspection times
#' t <- seq(0, 40, by = 2)
#' x <- levysim::sim_vg(t, theta = 0.5, sigma = 1.5, nu = 2, n_traj = 20)
#' vg_mle(x)
#'
#' # 3. Twenty trajectories, each with its own inspection times
#' plan <- inspection_times("dec-inc", n = 15:34, n_traj = 20, tmax = 40)
#' x <- levysim::sim_vg(plan, theta = 0.5, sigma = 1.5, nu = 2)
#' vg_mle(x)
#'
#' # 4. Estimates with 95 % confidence intervals
#' vg_mle(x, ci = TRUE)
#'
#' # 5. Increasing data, here from a gamma process: the likelihood is maximal
#' #    at the gamma limit (sigma = 0): warning suggesting gamma_mle(), NA estimates
#' y <- levysim::sim_gamma(seq(0, 100, by = 1), a = 2, b = 3)
#' vg_mle(y)
#' gamma_mle(y)
#' @export
vg_mle <- function(x, t = NULL, ci = FALSE, level = 0.95) {
  estimate("vg", x, t, ci, level)
}

vg_est <- function(x, t = NULL) {
  d <- increments(x, t)
  dt <- d$dt
  dx <- d$dx
  if (any(dx == 0)) stop("zero increment: log|dx| is not defined")
  m <- length(dx)
  if (m < 3) stop("at least three increments are needed to estimate three parameters")

  brownian <- "nu tends to 0: the likelihood is maximal at the Brownian limit, use wiener_mle()"
  theta <- sum(dx) / sum(dt)
  sigma <- sqrt(sum((dx - theta * dt)^2 / dt) / m)
  nu <- mean(dt)
  nu_min <- max(dt) / 300
  ll <- vg_loglik(dx, dt, theta, sigma, nu)

  for (it in 1:5000) {
    g <- gig_moments(dt / nu - 0.5, dx^2 / sigma^2, theta^2 / sigma^2 + 2 / nu)
    if (!all(is.finite(unlist(g))) || nu < nu_min) no_mle(brownian)
    theta <- sum(dx) / sum(g$EG)
    sigma <- sqrt(mean(dx^2 * g$EinvG - 2 * theta * dx + theta^2 * g$EG))
    q <- function(la) sum(dt * g$ElogG - g$EG + dt * la + dt - dt * digamma(exp(la) * dt))
    if (q(log(1 / nu_min)) > 0) no_mle(brownian)
    nu <- exp(-uniroot(q, interval = c(-log(1e6 * max(dt)), log(1 / nu_min)), tol = 1e-12)$root)
    ll_new <- vg_loglik(dx, dt, theta, sigma, nu)
    if (!is.finite(ll_new)) no_mle(brownian)
    if (ll_new - ll < -1e-8 * abs(ll)) warning("the log-likelihood decreased: numerical problem")
    conv <- abs(ll_new - ll) <= 1e-10 * (abs(ll) + 1e-10)
    ll <- ll_new
    if (conv) break
  }

  nll <- function(p) {
    v <- -vg_loglik(dx, dt, p[1], exp(p[2]), exp(p[3]))
    if (is.finite(v)) v else 1e10
  }
  s0 <- sd(dx / sqrt(dt))
  lower <- c(-Inf, log(s0 / 100), log(max(dt) / 300))
  upper <- c(Inf, log(100 * s0), log(1000 * max(dt)))
  p_em <- pmin(pmax(c(theta, log(sigma), log(nu)), lower), upper)
  opt <- tryCatch(optim(p_em, nll, method = "L-BFGS-B", lower = lower, upper = upper,
                        control = list(factr = 1, maxit = 1000)),
                  error = function(e) NULL)
  p <- if (!is.null(opt) && opt$value <= nll(p_em)) opt$par else p_em

  if (all(dx > 0)) {
    g <- gamma_fit(dt, dx)
    ll_gamma <- sum(dgamma(dx, shape = g[["a"]] * dt, rate = g[["b"]], log = TRUE))
    if (-nll(p) <= ll_gamma + 1e-8 * abs(ll_gamma)) {
      no_mle("sigma tends to 0: the likelihood is maximal at the gamma limit, use gamma_mle()")
    }
  }
  c(theta = p[1], sigma = exp(p[2]), nu = exp(p[3]))
}

vg_loglik <- function(dx, dt, theta, sigma, nu) {
  lam <- dt / nu
  cc <- sqrt(theta^2 + 2 * sigma^2 / nu)
  ts <- theta * sign(dx)
  d <- ifelse(ts > 0, (2 * sigma^2 / nu) / (cc + ts), cc - ts)
  z <- cc * abs(dx) / sigma^2
  sum(log(2) - 0.5 * log(2 * pi) - log(sigma) - lam * log(nu) - lgamma(lam) +
        (lam - 0.5) * (log(abs(dx)) - log(cc)) + log_besselKx(z, lam - 0.5) - abs(dx) * d / sigma^2)
}

log_besselKx <- function(z, p) {
  n <- max(length(z), length(p))
  z <- rep_len(z, n)
  p <- rep_len(abs(p), n)
  out <- rep(NA_real_, n)
  small <- p <= 1e4
  out[small] <- log(besselK(z[small], p[small], expon.scaled = TRUE))
  debye <- function(z, p) {
    s <- sqrt(1 + (z / p)^2)
    0.5 * log(pi / (2 * p)) - p * log((z / p) / (1 + s)) - 0.5 * log(s) - p^2 / (z + p * s)
  }
  for (i in which(!is.finite(out) & is.finite(z) & z > 0 & is.finite(p))) {
    out[i] <- if (p[i] > 1e4) debye(z[i], p[i]) else {
      g <- function(u) -2 * z[i] * sinh(u / 2)^2 + p[i] * u + log1p(exp(-2 * p[i] * u)) - log(2)
      opt <- optimize(g, c(0, asinh(p[i] / z[i]) + 1), maximum = TRUE, tol = 1e-12)
      w <- 12 / sqrt(z[i] * cosh(opt$maximum))
      I <- integrate(function(u) exp(g(u) - opt$objective), max(0, opt$maximum - w), opt$maximum + w,
                     rel.tol = 1e-12, subdivisions = 1000, stop.on.error = FALSE)$value
      opt$objective + log(I)
    }
  }
  out
}

gig_moments <- function(p, chi, psi, h = 1e-5) {
  om <- sqrt(chi * psi)
  eta <- sqrt(chi / psi)
  lk <- log_besselKx(om, p)
  list(EG = eta * exp(log_besselKx(om, p + 1) - lk),
       EinvG = exp(log_besselKx(om, p - 1) - lk) / eta,
       ElogG = log(eta) + (log_besselKx(om, p + h) - log_besselKx(om, p - h)) / (2 * h))
}
