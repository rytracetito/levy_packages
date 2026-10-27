#' levyfit: maximum likelihood estimation for degradation processes
#'
#' `levyfit` estimates the parameters of stochastic degradation models from one
#' or several degradation trajectories observed at inspection times, which can
#' be irregular and differ from one trajectory to another. Each estimator is the
#' maximum likelihood estimator (MLE); the help page of each function details
#' the main mathematical steps of its derivation.
#'
#' @section Functions of the package:
#' Estimation of the parameters (click on a function name for its help page):
#'
#' | Function | Process | Parameters | Computation |
#' |---|---|---|---|
#' | [wiener_mle()] | Wiener with drift | `mu`, `sigma2`, `s2` | closed form |
#' | [gamma_mle()] | gamma | shape `a`, rate `b` | root of one equation |
#' | [ig_mle()] | inverse Gaussian | `mu`, `eta`, `eta_unbiased` | closed form |
#' | [vg_mle()] | variance-gamma | `theta`, `sigma`, `nu` | EM, then L-BFGS-B |
#' | [meixner_mle()] | Meixner | `a`, `b`, `d`, `m` | profile likelihood |
#'
#' `s2` and `eta_unbiased` are the unbiased versions of `sigma2` and `eta`.
#' Every estimation function returns the point estimates and, with
#' `ci = TRUE`, asymptotic confidence intervals (see [wiener_mle()]). When the
#' estimator does not exist for the data, it gives a warning and returns `NA`
#' estimates instead of stopping.
#'
#' Inspection times following a chosen pattern, and simulated trajectories
#' to try the estimators, are given by the package `DegradSim`
#' (`DegradSim::inspection_times()` and its simulation functions).
#'
#' @section Data:
#' Every estimation function accepts the data in three forms:
#' * one trajectory: `f(x, t)`, with `x` the degradation levels and `t` the
#'   inspection times;
#' * several trajectories with the same inspection times: `f(x)`, with `x` the
#'   matrix returned by the `DegradSim` simulation functions (columns `t`, `X1`,
#'   `X2`, ...);
#' * several trajectories, each with its own inspection times: `f(liste)`,
#'   with `liste` a list of data frames, one per trajectory, whose first column
#'   contains the inspection times and second column the degradation levels.
#'
#' The trajectories are assumed independent and with the same law. The
#' increments are computed within each trajectory and then pooled.
#'
#' @section Computation of the estimators:
#' With \eqn{(\Delta t_j, \Delta x_j)}, \eqn{j = 1, \dots, m}, the durations and
#' increments between two consecutive inspections, pooled over the
#' trajectories, the log-likelihood of a Lévy process is
#' \eqn{\ell(\psi) = \sum_j \ln f_{\Delta t_j}(\Delta x_j; \psi)}.
#' * **Closed form**: Wiener and inverse Gaussian. The MLE is explicit, with
#'   \eqn{\hat\mu = S / T}, \eqn{S} the total degradation and \eqn{T} the total
#'   observed time.
#' * **Root of one equation**: gamma. \eqn{\hat b = \hat a T / S} is explicit
#'   and \eqn{\hat a} is the unique root of a decreasing function, found by
#'   [stats::uniroot()] from the method of moments estimator; the result does
#'   not depend on this starting point.
#' * **Numerical maximization**: variance-gamma and Meixner.
#'   * Variance-gamma: EM algorithm, the gamma clock being the hidden
#'     variable (each step is explicit or a unique root), started from the
#'     Wiener MLE, then L-BFGS-B ([stats::optim()]) with a numerical gradient.
#'   * Meixner: L-BFGS-B on the likelihood profiled in \eqn{b} (\eqn{\hat b} is
#'     explicit for given \eqn{a}, \eqn{d}, \eqn{m}), with an explicit gradient
#'     and several starting points (`n_start`) built from the empirical mean
#'     and variance; the best interior maximum is kept.
#'
#' When the supremum of the likelihood is reached at the boundary of the
#' parameter space (for instance the Brownian limit of the variance-gamma
#' process), the MLE does not exist: the function gives a warning and returns
#' `NA`. The help page of each function details its derivation.
#'
#' @section Confidence intervals:
#' With `ci = TRUE`, the intervals are Wald intervals. The observed Fisher
#' information \eqn{J = -\nabla^2 \ell(\hat\psi)} is computed numerically by
#' [stats::optimHess()] for every model, and the standard error of
#' \eqn{\hat\psi_k} is \eqn{\widehat{se}_k = \sqrt{(J^{-1})_{kk}}}. With
#' \eqn{z} the quantile of order \eqn{(1 + level) / 2} of the standard normal
#' distribution:
#' * real parameter (`mu` of the Wiener process, `theta` of the variance-gamma
#'   process, `b` and `m` of the Meixner process):
#'   \eqn{\hat\psi_k \pm z\, \widehat{se}_k};
#' * positive parameter (all the others): the interval
#'   \eqn{\ln \hat\psi_k \pm z\, \widehat{se}_k / \hat\psi_k} for
#'   \eqn{\ln \psi_k} (delta method), transformed back:
#'   \eqn{[\hat\psi_k e^{-z \widehat{se}_k / \hat\psi_k},\ \hat\psi_k e^{z \widehat{se}_k / \hat\psi_k}]},
#'   which always stays positive.
#'
#' For the Wiener, inverse Gaussian and gamma processes, the numerical
#' information agrees with its closed form; for instance, for the Wiener
#' process, \eqn{\widehat{se}(\hat\mu) = \hat\sigma / \sqrt{T}} and
#' \eqn{\widehat{se}(\hat\sigma^2) = \hat\sigma^2 \sqrt{2 / m}}.
#'
#' @section Parametrizations:
#' * Wiener: \eqn{X_t - X_s \sim \mathcal{N}(\mu (t - s), \sigma^2 (t - s))}.
#' * Gamma: \eqn{X_t - X_s \sim \Gamma(a (t - s), b)}, shape \eqn{a(t - s)} and
#'   **rate** \eqn{b}, so that \eqn{E[X_t] = a t / b}. The scale is \eqn{1 / b}.
#' * Inverse Gaussian: \eqn{X_t - X_s \sim IG(\mu (t - s), \eta (t - s)^2)}.
#' * Variance-gamma: \eqn{X_s = \theta \gamma_s + \sigma W_{\gamma_s}}, with
#'   \eqn{\gamma_s} of mean \eqn{s} and variance \eqn{\nu s}.
#' * Meixner: \eqn{X_s \sim Meixner(a, b, d s, m s)}.
#'
#' @references
#' **Gamma process**: van Noortwijk, J. M. (2009). A survey of the application
#' of gamma processes in maintenance. *Reliability Engineering & System
#' Safety*, 94(1), 2-21. \doi{10.1016/j.ress.2007.03.019}
#'
#' **Wiener process**: Whitmore, G. A. (1995). Estimating degradation by a
#' Wiener diffusion process subject to measurement error. *Lifetime Data
#' Analysis*, 1, 307-319.
#'
#' **Inverse Gaussian process**: Wang, X. and Xu, D. (2010). An inverse
#' Gaussian process model for degradation data. *Technometrics*, 52(2),
#' 188-197. \doi{10.1198/TECH.2009.08197}
#'
#' **Variance-gamma process**: Madan, D. B., Carr, P. P. and Chang, E. C.
#' (1998). The variance gamma process and option pricing. *European Finance
#' Review*, 2(1), 79-105. \doi{10.1023/A:1009703431535}
#'
#' **Meixner process**: Schoutens, W. and Teugels, J. L. (1998). Lévy
#' processes, polynomials and martingales. *Communications in Statistics.
#' Stochastic Models*, 14(1-2), 335-349.
#'
#' The references of each method are given in the help page of each function.
#'
#' @examplesIf requireNamespace("DegradSim", quietly = TRUE)
#' # 1. Inspection times: 10 trajectories inspected more and more often
#' plan <- DegradSim::inspection_times("increasing", n = 20, n_traj = 10, tmax = 100)
#'
#' # 2. Simulated data, here a gamma process (a = 2, b = 3)
#' set.seed(1)
#' x <- DegradSim::sim_gamma(plan, a = 2, b = 3)
#'
#' # 3. Estimation of several models on the same data
#' gamma_mle(x)
#' wiener_mle(x)
#' ig_mle(x)
#'
#' # 4. Estimates with 95 % confidence intervals
#' gamma_mle(x, ci = TRUE)
#' @keywords package
"_PACKAGE"

#' @importFrom stats dgamma dnorm integrate optim optimHess optimize qnorm runif sd uniroot
NULL
