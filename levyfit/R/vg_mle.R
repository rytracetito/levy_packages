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
#' (Madan and Seneta, 1990, section III.A, equation (7), for \eqn{\theta = 0};
#' Seneta, 2004, equations (1) and (2); Bee, Dickson and Santi, 2018,
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
#' 2023, appendix A), so that \eqn{K_{-p} = K_p}. Some authors call
#' \eqn{K_p} the modified Bessel function of the third kind; Seneta (2004,
#' section 3) points out this ambiguity of terminology. For
#' \eqn{\Delta t = 1}, this is the density (15) of Seneta (2004) with his
#' location \eqn{c = 0}, and the density (1) of Bee et al. (2018) with their
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
#' which is equation (4) of Bee et al. (2018) for \eqn{\Delta t_j = 1}.
#'
#' **Why an EM algorithm.** The density contains the modified Bessel function
#' \eqn{K_p}, whose order depends on \eqn{\nu}, and the likelihood cannot be
#' maximized in closed form. Direct maximization is likely to be
#' computationally expensive (Madan and Seneta, 1990, section II.F), the
#' log-likelihood can have many local optima (Belhaj Salem, Fouladirad and
#' Deloux, 2022, section 2.2), and Seneta (2004, section 3) reports, for the
#' model with a location parameter, a likelihood surface that appeared flat
#' near the moment estimates. Standard optimization algorithms
#' applied directly to this likelihood, such as the Nelder-Mead algorithm of
#' [stats::optim()], are very sensitive to the starting values and give
#' highly unstable estimates (Bee et al., 2018, section 1, after Cervellera
#' and Tucci). On the
#' other hand, the maximization would be easy if the random times were known:
#' this is why the MLE is computed by an EM algorithm, which treats the
#' random times as missing data (Bee et al., 2018, sections 1 and 3.2).
#'
#' @section When the estimator does not exist:
#' The VG model has two other models as limits, which are not in its
#' parameter space (\eqn{\nu > 0}, \eqn{\sigma > 0}). The likelihood can
#' increase towards one of these limits, without a maximum inside the
#' parameter space: the MLE then does not exist, and the function returns
#' `NA` with a warning that names the model to use instead.
#' * \eqn{\nu \to 0}: with the parameters above, the characteristic function
#'   of an increment (Fischer et al., 2023, equation (2.12)) is
#'   \eqn{\varphi(u) = (1 - i \theta \nu u + \sigma^2 \nu u^2 / 2)^{-\Delta t / \nu}},
#'   which tends to \eqn{\exp(\Delta t (i \theta u - \sigma^2 u^2 / 2))}, the
#'   characteristic function of \eqn{\mathcal{N}(\theta \Delta t, \sigma^2 \Delta t)}.
#'   By Levy's continuity theorem (van der Vaart, 1998, Theorem 2.13), each
#'   increment converges in distribution to this normal distribution: the VG
#'   process tends to the Wiener process with drift \eqn{\theta} and
#'   diffusion \eqn{\sigma}. For \eqn{\theta = 0}, Seneta (2004, section 3)
#'   gives the same limit. The function concludes that the likelihood
#'   increases towards this limit when \eqn{\nu} falls below
#'   \eqn{\max_j \Delta t_j / 300} during the EM algorithm, or when L-BFGS-B
#'   stops on this bound. The value 300 is a choice of this implementation,
#'   without a source. Its meaning: for \eqn{\theta = 0}, the kurtosis of an
#'   increment over a duration \eqn{\Delta t} is \eqn{3 (1 + \nu / \Delta t)}
#'   (Madan and Seneta, 1990, section III.C), so that at the threshold the
#'   excess kurtosis of the increment over the longest duration is
#'   \eqn{3 / 300 = 0.01}: this increment is then almost normal. Use
#'   [wiener_mle()].
#' * \eqn{\sigma \to 0}: when \eqn{\theta > 0}, an increment tends in
#'   distribution to the gamma distribution with shape \eqn{\Delta t / \nu} and
#'   rate \eqn{1 / (\theta \nu)} (Fischer et al., 2023, section 2.5.1, item 2,
#'   with their \eqn{r} equal to \eqn{\Delta t / \nu} and their \eqn{\lambda}
#'   equal to \eqn{1 / (\theta \nu)}): the VG process tends to the gamma
#'   process with shape \eqn{1 / \nu} and rate \eqn{1 / (\theta \nu)}. The
#'   density converges too: applying the limiting form (A.65) of \eqn{K_p} to
#'   the density (1.1) of Fischer et al. (2023) shows that, as
#'   \eqn{\sigma \to 0}, the VG density tends to this gamma density at every
#'   \eqn{x > 0}, and to 0 at every \eqn{x < 0}. The VG likelihood therefore
#'   tends to the gamma likelihood. When all the increments are positive, the
#'   function compares the VG likelihood at the estimate with the maximum of
#'   the gamma likelihood; if the VG likelihood is not larger, the function
#'   concludes that the likelihood increases towards this limit. Use
#'   [gamma_mle()].
#'
#' When a computation of the EM algorithm gives a value that is not finite,
#' the function returns `NA` with a warning that reports this numerical
#' failure, without naming a limit.
#'
#' @section EM algorithm:
#' **The EM algorithm in general** (Bee et al., 2018, section 3.2). The data
#' are completed by unobserved variables (the missing data), chosen so that
#' the log-likelihood of the complete data is simple. Starting from initial
#' values of the parameters, the algorithm repeats two steps:
#' * E step (expectation): compute the expected value of the complete
#'   log-likelihood given the observed data, at the current value of the
#'   parameters;
#' * M step (maximization): maximize this expected value with respect to the
#'   parameters, which gives their new value.
#'
#' The likelihood of the observed data increases at each iteration (Bee et
#' al., 2018, section 3.2).
#'
#' **Application to the VG process.** The algorithm used here is the one of
#' Bee et al. (2018, section 3.2), who follow Nitithumbundit and Chan (2015),
#' written for increments over unequal durations and without location
#' parameter. The correspondence with the general algorithm is the following.
#' * Observed data: the increments \eqn{\Delta x_j}, \eqn{j = 1, \dots, n}.
#' * Missing data: the increments \eqn{G_j} of the random time \eqn{\gamma} on
#'   the \eqn{j}-th interval.
#' * Complete log-likelihood: since
#'   \deqn{\Delta x_j \mid G_j \sim \mathcal{N}(\theta G_j, \sigma^2 G_j), \qquad
#'     G_j \sim \Gamma(\Delta t_j / \nu,\ 1 / \nu),}
#'   it is the sum of normal and gamma log-densities, and it only involves
#'   \eqn{G_j}, \eqn{1 / G_j} and \eqn{\ln G_j} (Bee et al., 2018, equation
#'   (5)).
#' * E step: compute the expectations \eqn{e^+_j}, \eqn{e^-_j} and \eqn{l_j}
#'   of \eqn{G_j}, \eqn{1 / G_j} and \eqn{\ln G_j} given \eqn{\Delta x_j}
#'   (see below).
#' * M step: update \eqn{(\theta, \sigma^2)}, then \eqn{\nu} (see Steps
#'   below).
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
#' **Start and end.** As recommended by Seneta (2004, section 3), the
#' algorithm starts from moment estimates.
#' * \eqn{\theta} and \eqn{\sigma}: the Wiener MLE (see [wiener_mle()]),
#'   \eqn{\theta_0 = S / T} and
#'   \eqn{\sigma_0^2 = \frac{1}{n} \sum_{j=1}^{n} (\Delta x_j - \theta_0 \Delta t_j)^2 / \Delta t_j},
#'   which estimate the mean and the variance of the process per unit of
#'   time. Seneta (2004, equations (23) and (24)) uses the same estimates for
#'   \eqn{\theta} and \eqn{\sigma^2} when the terms in \eqn{\theta^2},
#'   \eqn{\theta^3}, \eqn{\theta^4} are ignored (his location \eqn{c} being 0
#'   here).
#' * \eqn{\nu}: the kurtosis of \eqn{X_1} is \eqn{3 (1 + \nu)} when these terms
#'   are ignored (Seneta, 2004, equation (26); exactly when \eqn{\theta = 0},
#'   Madan and Seneta, 1990, sections II.B and II.E), so that
#'   \eqn{\nu_0 = \hat\gamma_2 / 3}, where \eqn{\hat\gamma_2 = k_4 / k_2^2} is
#'   the excess kurtosis of \eqn{X_1} estimated with the residuals
#'   \eqn{e_j = \Delta x_j - \theta_0 \Delta t_j} by
#'   \eqn{k_2 = \frac{1}{T} \sum_{j=1}^{n} e_j^2} and
#'   \eqn{k_4 = \frac{1}{T} \big(\sum_{j=1}^{n} e_j^4 - 3 k_2^2 \sum_{j=1}^{n} \Delta t_j^2\big)}.
#'   These estimators of the cumulants per unit of time, which take the
#'   unequal durations into account, are an extension made in this
#'   implementation: they use the fact that the cumulants of an increment are
#'   proportional to its duration \eqn{\Delta t}, the logarithm of its
#'   characteristic function (see When the estimator does not exist) being
#'   proportional to \eqn{\Delta t}. When \eqn{\nu_0} is not larger than the threshold
#'   \eqn{\max_j \Delta t_j / 300} (see When the estimator does not exist),
#'   in particular when \eqn{\hat\gamma_2 \le 0}, the moment estimate is not
#'   usable and the algorithm starts from \eqn{\nu_0 = T / n}, a choice of
#'   this implementation, without a source.
#'
#' The maximization is then
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
#' \doi{10.1016/j.ress.2022.108417}. Section 2.2.
#'
#' Fischer, A., Gaunt, R. E. and Sarantsev, A. (2023). The variance-gamma
#' distribution: a review. arXiv:2303.05615. Equations (1.1), (2.12) and
#' (2.18), section 2.5.1 and appendix A (limiting form (A.65)).
#'
#' Madan, D. B. and Seneta, E. (1990). The variance gamma (V.G.) model for
#' share market returns. *Journal of Business*, 63(4), 511-524.
#' \url{https://www.jstor.org/stable/2353303}. Sections II.B, II.E, II.F, III.A and III.C.
#'
#' Nitithumbundit, T. and Chan, J. S. K. (2015). An ECM algorithm for skewed
#' multivariate variance gamma distribution in normal mean-variance
#' representation. arXiv:1504.01239.
#'
#' van der Vaart, A. W. (1998). *Asymptotic Statistics*. Cambridge University
#' Press. \doi{10.1017/CBO9780511802256}. Theorem 2.13.
#'
#' Seneta, E. (2004). Fitting the variance-gamma model to financial data.
#' *Journal of Applied Probability*, 41A, 177-187.
#' \doi{10.1017/S0021900200112288}. Equations (1), (2), (15) and (23) to
#' (26), section 3.
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
#' # 5. Increasing data, here from a gamma process: the likelihood increases
#' #    towards the gamma limit (sigma = 0): warning suggesting gamma_mle(), NA estimates
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

  brownian <- "nu tends to 0: the likelihood increases towards the Brownian limit, use wiener_mle()"
  failure <- "numerical failure of the EM algorithm (a value is not finite)"
  theta <- sum(dx) / sum(dt)
  sigma <- sqrt(sum((dx - theta * dt)^2 / dt) / m)
  nu_min <- max(dt) / 300
  e <- dx - theta * dt
  k2 <- sum(e^2) / sum(dt)
  k4 <- (sum(e^4) - 3 * k2^2 * sum(dt^2)) / sum(dt)
  nu <- k4 / k2^2 / 3
  if (!is.finite(nu) || nu <= nu_min) nu <- mean(dt)
  ll <- vg_loglik(dx, dt, theta, sigma, nu)

  for (it in 1:5000) {
    g <- gig_moments(dt / nu - 0.5, dx^2 / sigma^2, theta^2 / sigma^2 + 2 / nu)
    if (nu < nu_min) no_mle(brownian)
    if (!all(is.finite(unlist(g)))) no_mle(failure)
    theta <- sum(dx) / sum(g$EG)
    sigma <- sqrt(mean(dx^2 * g$EinvG - 2 * theta * dx + theta^2 * g$EG))
    q <- function(la) sum(dt * g$ElogG - g$EG + dt * la + dt - dt * digamma(exp(la) * dt))
    if (q(log(1 / nu_min)) > 0) no_mle(brownian)
    nu <- exp(-uniroot(q, interval = c(-log(1e6 * max(dt)), log(1 / nu_min)), tol = 1e-12)$root)
    ll_new <- vg_loglik(dx, dt, theta, sigma, nu)
    if (!is.finite(ll_new)) no_mle(failure)
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
      no_mle("sigma tends to 0: the likelihood increases towards the gamma limit, use gamma_mle()")
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
