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
#' | Function | Model | Parameters | Method |
#' |---|---|---|---|
#' | [wiener_mle()] | Wiener process with drift | `mu`, `sigma2` (and unbiased `s2`) | closed form |
#' | [gamma_mle()] | homogeneous gamma process | shape `a`, rate `b` | root of one equation |
#' | [ig_mle()] | inverse Gaussian process | `mu`, `eta` (and unbiased `eta_unbiased`) | closed form |
#' | [vg_mle()] | variance-gamma process | `theta`, `sigma`, `nu` | EM algorithm then quasi-Newton |
#' | [meixner_mle()] | Meixner process | `a`, `b`, `d`, `m` | profile likelihood, several starts |
#' | [igam_mle()] | inverse gamma process | `theta`, `alpha`, `beta` | exact likelihood by forward filtering |
#'
#' Confidence intervals:
#'
#' | Function | Use |
#' |---|---|
#' | [confint_mle()] | estimates and asymptotic confidence intervals of the parameters of any of the six models |
#'
#' Inspection schedules:
#'
#' | Function | Use |
#' |---|---|
#' | [inspection_times()] | inspection times on \eqn{[0, t_{max}]} for one or several trajectories, with a regular, increasing, decreasing or alternating number of inspections over time |
#'
#' @section Data:
#' Every estimation function accepts the data in two forms:
#' * one trajectory: `f(x, t)`, with `x` the degradation levels and `t` the
#'   inspection times;
#' * several trajectories: `f(liste)`, with `liste` a list of data frames, one
#'   per trajectory, whose first column contains the inspection times and
#'   second column the degradation levels.
#'
#' The trajectories are assumed independent and with the same law. For the
#' Lévy processes, the increments are computed within each trajectory and then
#' pooled; for the inverse gamma process, which is not a Lévy process, the
#' log-likelihoods of the trajectories are added.
#'
#' @section Parametrizations:
#' * Wiener: \eqn{X_t - X_s \sim N(\mu (t - s), \sigma^2 (t - s))}.
#' * Gamma: \eqn{X_t - X_s \sim \Gamma(a (t - s), b)}, shape \eqn{a(t - s)} and
#'   **rate** \eqn{b}, so that \eqn{E[X_t] = a t / b}. The scale is \eqn{1 / b}.
#' * Inverse Gaussian: \eqn{X_t - X_s \sim IG(\mu (t - s), \eta (t - s)^2)}.
#' * Variance-gamma: \eqn{X_s = \theta \gamma_s + \sigma W_{\gamma_s}}, with
#'   \eqn{\gamma_s} of mean \eqn{s} and variance \eqn{\nu s}.
#' * Meixner: \eqn{X_s \sim Meixner(a, b, d s, m s)}.
#' * Inverse gamma: \eqn{W(t) = \inf\{w : T(w) > t\}}, with \eqn{T} a gamma
#'   process of shape \eqn{(w / \alpha)^\beta} and scale \eqn{\theta}.
#'
#' @examplesIf requireNamespace("levysim", quietly = TRUE)
#' # 1. Inspection times: 10 trajectories inspected more and more often
#' plan <- inspection_times("increasing", n = 20, n_traj = 10, tmax = 100)
#'
#' # 2. Simulated data, here a gamma process (a = 2, b = 3)
#' set.seed(1)
#' x <- levysim::sim_gamma(plan, a = 2, b = 3)
#'
#' # 3. Estimation of several models on the same data
#' gamma_mle(x)
#' wiener_mle(x)
#' ig_mle(x)
#'
#' # 4. Estimates with 95 % confidence intervals
#' confint_mle(x, model = "gamma")
#' @keywords package
"_PACKAGE"

#' @importFrom stats coef dgamma dnorm integrate lm optim optimHess optimize pgamma qnorm runif sd uniroot
NULL
