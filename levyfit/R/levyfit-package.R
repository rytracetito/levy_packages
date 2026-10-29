#' levyfit: maximum likelihood estimation for degradation processes
#'
#' `levyfit` estimates the parameters of stochastic degradation models from one
#' or several degradation trajectories observed at regular or irregular
#' inspection times. All the parameters are estimated by maximum likelihood
#' (maximum likelihood estimator, MLE). The help page of each function gives
#' the model, the likelihood and how the estimator is computed.
#'
#' @section Functions of the package:
#' Estimation of the parameters:
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
#' `ci = TRUE`, asymptotic confidence intervals. [increments_data()] gives the
#' data by their increments
#' instead of their levels.
#'
#' Inspection times following a chosen pattern, and simulated trajectories
#' to try the estimators, are given by the package `DegradSim`
#' (`DegradSim::inspection_times()` and its simulation functions).
#'
#' @section Data:
#' Every estimation function accepts the data in four forms (shown here with
#' [gamma_mle()]; the call is the same for all the estimation functions):
#' * one trajectory: `gamma_mle(x, t)`, with `x` the degradation levels and
#'   `t` the inspection times;
#' * several trajectories with the same inspection times: `gamma_mle(x)`,
#'   with `x` the matrix returned by the `DegradSim` simulation functions
#'   (columns `t`, `X1`, `X2`, ...);
#' * several trajectories, each with its own inspection times:
#'   `gamma_mle(liste)`, with `liste` a list of data frames, one per
#'   trajectory, whose first column contains the inspection times and second
#'   column the degradation levels;
#' * the increments themselves: `gamma_mle(increments_data(dx, dt))`, see
#'   [increments_data()].
#'
#' The trajectories are assumed independent and with the same distribution.
#' The increments are computed within each trajectory and then pooled.
#'
#' @section Notation:
#' The same notation is used in all the help pages of the package.
#' * \eqn{(\Delta t_j, \Delta x_j)}, \eqn{j = 1, \dots, n}: the durations and
#'   the degradation increments between two consecutive inspections of a same
#'   trajectory, pooled over all the trajectories, so that \eqn{n} is the
#'   total number of increments. The durations \eqn{\Delta t_j} can be equal
#'   (regular inspections) or not: all the formulas hold in both cases.
#' * \eqn{T = \sum_{j=1}^{n} \Delta t_j}, the total observed time, and
#'   \eqn{S = \sum_{j=1}^{n} \Delta x_j}, the total degradation.
#' * \eqn{\theta}: the vector of the parameters of the Levy process, for
#'   instance \eqn{\theta = (a, b)} for the gamma process.
#' * \eqn{f_{\Delta t}(x; \theta)}: the density of an increment of the
#'   process over a duration \eqn{\Delta t}.
#' * \eqn{\psi} and \eqn{\psi'}: the digamma and trigamma functions, the first
#'   and second derivatives of \eqn{\ln \Gamma}.
#'
#' @section Computation of the estimators:
#' The increments of a Levy process being independent, the log-likelihood is
#' \deqn{\ell(\theta) = \sum_{j=1}^{n} \ln f_{\Delta t_j}(\Delta x_j; \theta).}
#' * **Closed form**: Wiener and inverse Gaussian. The MLE is explicit, with
#'   \eqn{\hat\mu = S / T}.
#' * **Root of one equation**: gamma. The estimator \eqn{\hat b = \hat a T / S} is explicit
#'   and \eqn{\hat a} is the unique root of a decreasing function, found by
#'   [stats::uniroot()] from the method of moments estimator; the result does
#'   not depend on this starting point.
#' * **Numerical maximization**: variance-gamma and Meixner.
#'   * Variance-gamma: the process is a Brownian motion observed at a random
#'     time given by a gamma process. The EM algorithm treats this random time
#'     as a hidden variable (each step is explicit or a unique root); it starts
#'     from the Wiener MLE and is followed by L-BFGS-B ([stats::optim()]).
#'   * Meixner: L-BFGS-B ([stats::optim()]) on the likelihood profiled in
#'     \eqn{b} (\eqn{\hat b} is explicit for given \eqn{a}, \eqn{d}, \eqn{m}),
#'     with the exact gradient, from several starting points: the moment
#'     estimates when they exist, and points built from the empirical mean
#'     and variance (`n_start`); the solution of largest likelihood is kept.
#'
#' @section Confidence intervals:
#' With `ci = TRUE`, the intervals are asymptotic intervals based on the
#' asymptotic normality of the MLE. For i.i.d. observations, that is equal
#' durations, the MLE is asymptotically normal with covariance matrix the
#' inverse of the Fisher information (van der Vaart, 1998, Theorem 5.39);
#' with unequal durations the increments are independent but not identically
#' distributed, and the same normal approximation is used. Let
#' \eqn{\hat\theta} be the MLE. The Fisher information is estimated by the
#' observed information
#' \eqn{J = -\partial^2 \ell / \partial \theta \partial \theta^\top} taken at
#' \eqn{\hat\theta} (Grigoletto and Provasi, 2009, section 4.2), computed
#' numerically by [stats::optimHess()], except for the Meixner process whose
#' Hessian is explicit (see [meixner_mle()]). The standard error of the
#' \eqn{k}-th parameter is \eqn{se_k = \sqrt{(J^{-1})_{kk}}}; replacing the
#' unknown standard deviation by a consistent estimate keeps the asymptotic
#' level of the interval (van der Vaart, 1998, section 3.2). Let
#' \eqn{z = \Phi^{-1}((1 + level) / 2)}, where \eqn{\Phi} is the distribution
#' function of the standard normal distribution: \eqn{(1 + level) / 2} is the
#' order of the quantile for a two-sided interval, for instance
#' \eqn{z = 1.96} for `level = 0.95`. The interval is
#' * \eqn{\hat\theta_k \pm z\, se_k} for a real parameter (`mu` of the
#'   Wiener process, `theta` of the variance-gamma process, `b` and `m` of the
#'   Meixner process);
#' * \eqn{[\hat\theta_k e^{-z\, se_k / \hat\theta_k},\ \hat\theta_k e^{z\, se_k / \hat\theta_k}]}
#'   for a positive parameter (all the others): this is the interval
#'   \eqn{\ln \hat\theta_k \pm z\, se_k / \hat\theta_k} for
#'   \eqn{\ln \theta_k}, whose standard error \eqn{se_k / \hat\theta_k} is
#'   given by the delta method (van der Vaart, 1998, Theorem 3.1); since
#'   \eqn{\ln} is monotone, this interval is transformed back into an interval
#'   for \eqn{\theta_k} (van der Vaart, 1998, section 3.2), which always stays
#'   positive.
#'
#' The unbiased estimates `s2` (Wiener) and `eta_unbiased` (inverse Gaussian)
#' are given without interval. When \eqn{J} is not positive definite, the
#' bounds are `NA`.
#'
#' For the Wiener process, the information matrix of Kahle, Mercier and
#' Paroissin (2016, equation \[1.17\]), written for the \eqn{n} increments,
#' gives \eqn{se(\hat\mu) = \hat\sigma / \sqrt{T}} and
#' \eqn{se(\hat\sigma^2) = \hat\sigma^2 \sqrt{2 / n}}.
#'
#' @section Parametrizations:
#' * Wiener: \eqn{X_t - X_s \sim \mathcal{N}(\mu (t - s), \sigma^2 (t - s))}
#'   (Kahle et al., 2016, definition 1.11).
#' * Gamma: \eqn{X_t - X_s \sim \Gamma(a (t - s), b)}, with shape
#'   \eqn{a(t - s)} and rate \eqn{b} (scale \eqn{1 / b}), so that
#'   \eqn{E[X_t] = a t / b} (Kahle et al., 2016, definition 2.24 and
#'   section 2.7).
#' * Inverse Gaussian: \eqn{X_t - X_s \sim IG(\mu (t - s), \eta (t - s)^2)},
#'   the process of Wang and Xu (2010) with \eqn{\Lambda(t) = \mu t} (see
#'   [ig_mle()] for the correspondence of the parameters).
#' * Variance-gamma: \eqn{X_s = \theta \gamma_s + \sigma W_{\gamma_s}}, with
#'   \eqn{\gamma_s} of mean \eqn{s} and variance \eqn{\nu s} (Seneta, 2004;
#'   Bee, Dickson and Santi, 2018).
#' * Meixner: \eqn{X_s \sim Meixner(a, b, d s, m s)}, the Levy process
#'   associated with the Meixner distribution of Grigoletto and Provasi
#'   (2009) (see [meixner_mle()]).
#'
#' @references
#' **Degradation processes**: Kahle, W., Mercier, S. and Paroissin, C. (2016).
#' *Degradation Processes in Reliability*. ISTE and Wiley, London. Chapter 1
#' (Wiener process) and chapter 2 (gamma process, in particular section 2.7 on
#' the maximum likelihood and moment estimators).
#'
#' **Gamma process**: van Noortwijk, J. M. (2009). A survey of the application
#' of gamma processes in maintenance. *Reliability Engineering & System
#' Safety*, 94(1), 2-21. \doi{10.1016/j.ress.2007.03.019}. Section 3.7.1
#' (maximum likelihood estimation).
#'
#' **Wiener process**: Whitmore, G. A. (1995). Estimating degradation by a
#' Wiener diffusion process subject to measurement error. *Lifetime Data
#' Analysis*, 1, 307-319.
#'
#' **Inverse Gaussian process**: Wang, X. and Xu, D. (2010). An inverse
#' Gaussian process model for degradation data. *Technometrics*, 52(2),
#' 188-197. \doi{10.1198/TECH.2009.08197}
#'
#' **Variance-gamma process**: Seneta, E. (2004). Fitting the variance-gamma
#' model to financial data. *Journal of Applied Probability*, 41A, 177-187.
#' \doi{10.1017/S0021900200112288}. Bee, M., Dickson, M. M. and Santi, F.
#' (2018). Likelihood-based risk estimation for variance-gamma models.
#' *Statistical Methods & Applications*, 27(1), 69-89.
#' \doi{10.1007/s10260-017-0393-z}
#'
#' **Meixner process**: Grigoletto, M. and Provasi, C. (2009). Simulation and
#' estimation of the Meixner distribution. *Communications in Statistics -
#' Simulation and Computation*, 38(1), 58-77. \doi{10.1080/03610910802395679}
#'
#' **Confidence intervals**: van der Vaart, A. W. (1998). *Asymptotic
#' Statistics*. Cambridge University Press. \doi{10.1017/CBO9780511802256}.
#' Theorem 3.1, section 3.2 and Theorem 5.39.
#'
#' The references of each method are given in the help page of each function.
#'
#' @examplesIf requireNamespace("DegradSim", quietly = TRUE)
#' # 1. Inspection times: 10 trajectories inspected more and more often
#' plan <- DegradSim::inspection_times("decreasing", n = 20, n_traj = 10, tmax = 100)
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
