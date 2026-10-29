#' Maximum likelihood estimation for the variance-gamma process
#'
#' Computes the maximum likelihood estimator (MLE) of the parameters
#' \eqn{(\theta, \sigma, \nu)} of the variance-gamma (VG) process from one or
#' several trajectories observed at regular or irregular inspection times.
#' There is no closed form: the MLE is computed by an EM algorithm, whose
#' steps are explicit or reduce to an equation with a unique root, followed by
#' a quasi-Newton refinement.
#'
#' @section Model and data:
#' Let \eqn{\gamma} be a gamma process such that \eqn{\gamma_s} has mean
#' \eqn{s} and variance \eqn{\nu s} (shape \eqn{s / \nu}, rate \eqn{1 / \nu}),
#' and \eqn{W} an independent standard Brownian motion. The VG process is
#' \deqn{X_s = \theta \gamma_s + \sigma W_{\gamma_s}:}
#' a Brownian motion with drift \eqn{\theta} and diffusion \eqn{\sigma},
#' observed at the random time \eqn{\gamma_s} instead of the time \eqn{s}
#' (Seneta, 2004, equations (1) and (2); Bee, Dickson and Santi, 2018,
#' section 2, where the increments of the random time over one unit of time
#' follow a gamma distribution with shape and rate both equal to their
#' \eqn{\nu}, that is \eqn{1 / \nu} here). It is a Levy process that can
#' decrease.
#'
#' Given the random time, an increment over a duration \eqn{\Delta t} is
#' \eqn{\theta G + \sigma \sqrt{G} Z}, with \eqn{G \sim \Gamma(\Delta t / \nu, 1 / \nu)}
#' and \eqn{Z \sim \mathcal{N}(0, 1)} independent. This is the representation
#' (2.18) of Fischer, Gaunt and Sarantsev (2023) with \eqn{S = 2 G / \nu},
#' so that the increment follows their distribution
#' \eqn{VG(r, \theta_F, \sigma_F, 0)} with
#' \eqn{r = 2 \Delta t / \nu}, \eqn{\theta_F = \theta \nu / 2} and
#' \eqn{\sigma_F = \sigma \sqrt{\nu / 2}}. Their density (1.1) then gives, with
#' \eqn{\lambda = \Delta t / \nu}, \eqn{c = \sqrt{\theta^2 + 2 \sigma^2 / \nu}}
#' and \eqn{\kappa = c / \sigma^2},
#' \deqn{f_{\Delta t}(x; \theta, \sigma, \nu) = \frac{2 e^{\theta x / \sigma^2}}{\sqrt{2\pi}\, \sigma\, \nu^{\lambda} \Gamma(\lambda)}
#'   \Big(\frac{|x|}{c}\Big)^{\lambda - 1/2} K_{|\lambda - 1/2|}(\kappa |x|), \qquad x \ne 0,}
#' where \eqn{K_p} is the modified Bessel function of the second kind of order
#' \eqn{p}, defined by
#' \eqn{K_p(x) = \int_0^\infty e^{-x \cosh t} \cosh(p t)\, dt} (Fischer et al.,
#' 2023, appendix A), so that \eqn{K_{-p} = K_p}. For
#' \eqn{\Delta t = 1}, this is the density (1) of Bee et al. (2018) with their
#' location \eqn{\mu = 0}, their \eqn{\gamma} equal to \eqn{\theta} and their
#' \eqn{\nu} equal to \eqn{1 / \nu}.
#'
#' The data are the \eqn{n} increments \eqn{(\Delta t_j, \Delta x_j)},
#' \eqn{j = 1, \dots, n}, with \eqn{\Delta x_j \ne 0}, pooled over the
#' trajectories, with \eqn{T = \sum_{j=1}^{n} \Delta t_j} and
#' \eqn{S = \sum_{j=1}^{n} \Delta x_j} (see the section Notation of
#' [levyfit-package]).
#'
#' @section Likelihood:
#' The increments being independent, the log-likelihood is
#' \deqn{\ell(\theta, \sigma, \nu) = \sum_{j=1}^{n} \ln f_{\Delta t_j}(\Delta x_j; \theta, \sigma, \nu),}
#' which is equation (4) of Bee et al. (2018) for \eqn{\Delta t_j = 1}. As in
#' their section 3.2, the likelihood is maximized by an EM algorithm that
#' treats the random times as missing data.
#'
#' @section When the estimator does not exist:
#' The VG model contains two other models as limits, and the likelihood can be
#' maximal at one of them; the function then returns `NA` with a warning that
#' names the model to use instead.
#' * \eqn{\nu \to 0}: with the parameters above, the cumulant generating
#'   function of an increment (Fischer et al., 2023, equation (2.13)) is
#'   \eqn{-(\Delta t / \nu) \ln(1 - \theta \nu u - \sigma^2 \nu u^2 / 2)}, which
#'   tends to \eqn{\Delta t (\theta u + \sigma^2 u^2 / 2)}, the cumulant
#'   generating function of \eqn{\mathcal{N}(\theta \Delta t, \sigma^2 \Delta t)}:
#'   the VG process tends to the Wiener process with drift \eqn{\theta} and
#'   diffusion \eqn{\sigma} (convergence in distribution, by Levy's continuity
#'   theorem as in Fischer et al., 2023, section 2.5.1). Use [wiener_mle()].
#' * \eqn{\sigma \to 0}: when \eqn{\theta > 0}, an increment tends in
#'   distribution to the gamma distribution with shape \eqn{\Delta t / \nu} and
#'   rate \eqn{1 / (\theta \nu)} (Fischer et al., 2023, section 2.5.1, item 2,
#'   with their \eqn{r} equal to \eqn{\Delta t / \nu} and their \eqn{\lambda}
#'   equal to \eqn{1 / (\theta \nu)}): the VG process tends to the gamma
#'   process with shape \eqn{1 / \nu} and rate \eqn{1 / (\theta \nu)}. When all
#'   the increments are positive, the function compares the VG likelihood with
#'   the maximum of the gamma likelihood; if the VG likelihood is not larger,
#'   the function considers that the maximum is reached at this limit. Use
#'   [gamma_mle()].
#'
#' @section EM algorithm:
#' The algorithm is the one of Bee et al. (2018, section 3.2), who follow
#' Nitithumbundit and Chan (2015), written here for increments over unequal
#' durations and without location parameter.
#'
#' **Idea.** Let \eqn{G_j} be the increment of the random time \eqn{\gamma}
#' on the \eqn{j}-th interval, which is not observed. If the \eqn{G_j} were
#' known, the estimation would be easy, since
#' \deqn{\Delta x_j \mid G_j \sim \mathcal{N}(\theta G_j, \sigma^2 G_j), \qquad
#'   G_j \sim \Gamma(\Delta t_j / \nu,\ 1 / \nu),}
#' and the log-likelihood of the complete data \eqn{(\Delta x_j, G_j)},
#' \eqn{j = 1, \dots, n}, only involves \eqn{G_j}, \eqn{1 / G_j} and
#' \eqn{\ln G_j} (Bee et al., 2018, equation (5)). The EM algorithm replaces
#' these three quantities by their expectations given the observed
#' increments, computed at the current parameters (E step), and then
#' maximizes the complete log-likelihood (M step). The likelihood increases at
#' each iteration (Bee et al., 2018, section 3.2).
#'
#' **Distribution of the hidden times.** Given \eqn{\Delta x_j}, \eqn{G_j}
#' follows a generalized inverse Gaussian distribution \eqn{GIG(p_j, \chi_j, \omega)}
#' (Bee et al., 2018, section 3.2, after Embrechts, 1983),
#' with density proportional to \eqn{g^{p_j - 1} \exp(-(\chi_j / g + \omega g) / 2)},
#' \eqn{g > 0}, where
#' \deqn{p_j = \frac{\Delta t_j}{\nu} - \frac{1}{2}, \qquad
#'   \chi_j = \frac{\Delta x_j^2}{\sigma^2}, \qquad
#'   \omega = \frac{\theta^2}{\sigma^2} + \frac{2}{\nu}.}
#' With \eqn{z_j = \sqrt{\chi_j \omega}} and \eqn{\rho_j = \sqrt{\chi_j / \omega}},
#' its moments are ratios of Bessel functions (Bee et al., 2018, equations (6)
#' to (8)):
#' \deqn{e^+_j = E[G_j] = \rho_j \frac{K_{p_j + 1}(z_j)}{K_{p_j}(z_j)}, \qquad
#'   e^-_j = E[1 / G_j] = \frac{1}{\rho_j} \frac{K_{p_j - 1}(z_j)}{K_{p_j}(z_j)},}
#' \deqn{l_j = E[\ln G_j] = \ln \rho_j + \frac{\partial}{\partial p} \ln K_p(z_j) \Big|_{p = p_j},}
#' the last derivative being computed numerically, as in Bee et al. (2018).
#'
#' **Steps.**
#' 1. E step: compute \eqn{e^+_j}, \eqn{e^-_j} and \eqn{l_j},
#'    \eqn{j = 1, \dots, n}, at the current parameters.
#' 2. M step for \eqn{(\theta, \sigma^2)}, explicit (Bee et al., 2018,
#'    equations (10) and (11) with their location \eqn{\mu = 0}): with
#'    \eqn{E^+ = \sum_{j=1}^{n} e^+_j},
#'    \deqn{\theta = \frac{S}{E^+}, \qquad
#'      \sigma^2 = \frac{1}{n} \Big(\sum_{j=1}^{n} \Delta x_j^2 e^-_j - \frac{S^2}{E^+}\Big).}
#' 3. M step for \eqn{\nu}: with \eqn{a = 1 / \nu} and \eqn{\psi} the digamma
#'    function, the derivative of the complete log-likelihood in \eqn{a} is
#'    \deqn{q(a) = \sum_{j=1}^{n} \big[\Delta t_j l_j - e^+_j + \Delta t_j \ln a + \Delta t_j - \Delta t_j\, \psi(a \Delta t_j)\big],}
#'    which is the derivative given in the CM-step 2a of Bee et al. (2018)
#'    when \eqn{\Delta t_j = 1}. Since \eqn{\psi'(x) > 1/x} (Alzer, 1997,
#'    Theorem 1), \eqn{q} is strictly decreasing; by the inequality (2.2) of
#'    Alzer (1997) it tends to \eqn{+\infty} as \eqn{a \to 0} and, by Jensen's
#'    inequality (\eqn{l_j < \ln e^+_j}) and \eqn{\ln u \le u - 1}, to a
#'    negative limit as
#'    \eqn{a \to \infty}: its unique root gives \eqn{\nu = 1 / a}.
#' 4. Repeat steps 1 to 3 until the increase of \eqn{\ell} becomes small.
#'
#' **Start and end.** The algorithm starts from the Wiener MLE for
#' \eqn{(\theta, \sigma)} and from \eqn{\nu = T / n}. The maximization is then
#' finished by L-BFGS-B on \eqn{(\theta, \ln \sigma, \ln \nu)}, started at the
#' last EM iterate. If this refinement fails or does not increase the
#' likelihood, the EM iterate is returned.
#'
#' **Numerical detail.** The Bessel functions \eqn{K_p(z)} are computed on the
#' log scale; when `besselK()` overflows (large order \eqn{p}, that is small
#' \eqn{\nu}, and small \eqn{z}), \eqn{\ln K_p(z)} is computed exactly by the
#' recurrence on the order \eqn{K_{p+1}(z) = K_{p-1}(z) + (2p/z) K_p(z)}
#' (Olver et al., 2010, equation 10.29.1, with \eqn{\mathcal{Z}_\nu = e^{\nu \pi i} K_\nu}
#' as defined in section 10.25(ii)).
#'
#' @inheritParams wiener_mle
#' @param x Degradation levels of one trajectory, with non-zero increments;
#'   or a list of data frames, one per trajectory, whose first column contains
#'   the inspection times and second column the degradation levels.
#'   The increments can also be given directly with [increments_data()].
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
#' Alzer, H. (1997). On some inequalities for the gamma and psi functions.
#' *Mathematics of Computation*, 66(217), 373-389.
#' \doi{10.1090/S0025-5718-97-00807-7}. Equation (2.2) and Theorem 1.
#'
#' Bee, M., Dickson, M. M. and Santi, F. (2018). Likelihood-based risk
#' estimation for variance-gamma models. *Statistical Methods & Applications*,
#' 27(1), 69-89. \doi{10.1007/s10260-017-0393-z}. Sections 2 and 3.2,
#' equations (1) and (4) to (11).
#'
#' Belhaj Salem, M., Fouladirad, M. and Deloux, E. (2022). Variance Gamma
#' process as degradation model for prognosis and imperfect maintenance of
#' centrifugal pumps. *Reliability Engineering & System Safety*, 223, 108417.
#' \doi{10.1016/j.ress.2022.108417}
#'
#' Fischer, A., Gaunt, R. E. and Sarantsev, A. (2023). The variance-gamma
#' distribution: a review. arXiv:2303.05615. Equations (1.1), (2.13) and
#' (2.18), section 2.5.1 and appendix A.
#'
#' Nitithumbundit, T. and Chan, J. S. K. (2015). An ECM algorithm for skewed
#' multivariate variance gamma distribution in normal mean-variance
#' representation. arXiv:1504.01239.
#'
#' Seneta, E. (2004). Fitting the variance-gamma model to financial data.
#' *Journal of Applied Probability*, 41A, 177-187.
#' \doi{10.1017/S0021900200112288}. Equations (1) and (2).
#'
#' Olver, F. W. J., Lozier, D. W., Boisvert, R. F. and Clark, C. W. (eds.)
#' (2010). *NIST Handbook of Mathematical Functions*. Cambridge University
#' Press. Online version: *NIST Digital Library of Mathematical Functions*,
#' \url{https://dlmf.nist.gov/10.29}, equation 10.29.1.
#'
#' @seealso [wiener_mle()], [meixner_mle()]
#'
#' @examplesIf requireNamespace("DegradSim", quietly = TRUE)
#' # 1. One trajectory inspected every 2 time units up to 400
#' #    (theta = 0.5, sigma = 1.5, nu = 2)
#' t <- DegradSim::inspection_times("regular", n = 200, tmax = 400)$time
#' set.seed(1)
#' x <- DegradSim::sim_vg(t, theta = 0.5, sigma = 1.5, nu = 2)
#' vg_mle(x)
#'
#' # The same trajectory given as two vectors: degradation levels, then times
#' vg_mle(x[, "X1"], t)
#'
#' # 2. Twenty trajectories with the same inspection times
#' t <- DegradSim::inspection_times("regular", n = 20, tmax = 40)$time
#' x <- DegradSim::sim_vg(t, theta = 0.5, sigma = 1.5, nu = 2, n_traj = 20)
#' vg_mle(x)
#'
#' # 3. Twenty trajectories, each with its own inspection times
#' plan <- DegradSim::inspection_times("inc-dec", n = 15:34, n_traj = 20, tmax = 40)
#' x <- DegradSim::sim_vg(plan, theta = 0.5, sigma = 1.5, nu = 2)
#' vg_mle(x)
#'
#' # 4. Estimates with 95 % confidence intervals
#' vg_mle(x, ci = TRUE)
#'
#' # 5. Increasing data, here from a gamma process: the likelihood is maximal
#' #    at the gamma limit (sigma = 0): warning suggesting gamma_mle(), NA estimates
#' t100 <- DegradSim::inspection_times("regular", n = 100, tmax = 100)$time
#' y <- DegradSim::sim_gamma(t100, a = 2, b = 3)
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
    # l'EM s'arrête dès qu'il ralentit : L-BFGS-B termine la maximisation, bien plus vite près de nu = 0
    slow <- it >= 10 && ll_new - ll <= 1e-6 * (abs(ll) + 1)
    ll <- ll_new
    if (slow) break
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
  if (p[3] <= lower[3] + 1e-6) no_mle(brownian)

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

# log(K_p(z) exp(z)) par récurrence ascendante sur l'ordre, en logarithmes :
# K_{v+1}(z) = K_{v-1}(z) + (2 v / z) K_v(z), à partir de l'ordre fractionnaire v0 = p - floor(p)
log_besselKx_rec <- function(z, p) {
  n0 <- floor(p)
  v0 <- p - n0
  l <- log(besselK(z, v0, expon.scaled = TRUE))
  r <- besselK(z, v0 + 1, expon.scaled = TRUE) / besselK(z, v0, expon.scaled = TRUE)
  for (k in seq_len(max(n0))) {
    i <- n0 >= k
    l[i] <- l[i] + log(r[i])
    r <- 1 / r + 2 * (v0 + k) / z
  }
  l
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
  # débordement de besselK (grand ordre, petit z) : récurrence exacte sur l'ordre
  bad <- which(!is.finite(out) & is.finite(z) & z > 0 & p <= 1e4)
  if (length(bad) > 0) out[bad] <- log_besselKx_rec(z[bad], p[bad])
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
