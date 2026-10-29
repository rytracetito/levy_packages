#' Maximum likelihood estimation for the Meixner process
#'
#' Computes the maximum likelihood estimator of the four parameters
#' \eqn{(a, b, d, m)} of the Meixner process from one or several trajectories
#' observed at regular or irregular inspection times. The likelihood is
#' maximized by a Newton-type algorithm started at the moment estimates.
#'
#' @section Source and notation:
#' The results on the Meixner distribution used below come from Grigoletto
#' and Provasi (2009), who treat an i.i.d. sample of this
#' distribution. The distribution written here \eqn{Meixner(a, b, d, m)} is
#' their \eqn{MXN(m, a, b, d)}: \eqn{m} is a location parameter, \eqn{a > 0} a
#' scale parameter, and \eqn{b \in (-\pi, \pi)} and \eqn{d > 0} are shape
#' parameters (Grigoletto and Provasi, 2009, section 2); the same notation is used by Schoutens
#' (2002), whose Meixner process is used here. The extension of the results
#' of Grigoletto and Provasi (2009) to increments over unequal durations, and the algorithm, are
#' choices of this implementation; they are marked as such.
#'
#' @section Model and data:
#' The density of \eqn{Meixner(a, b, d, m)} is (Grigoletto and Provasi, 2009, section 2)
#' \deqn{f(x) = \frac{(2 \cos(b/2))^{2 d}}{2 a \pi\, \Gamma(2 d)}
#'   \exp\Big(\frac{b (x - m)}{a}\Big) \Big|\Gamma\Big(d + i\, \frac{x - m}{a}\Big)\Big|^2,
#'   \qquad x \in (-\infty, \infty),}
#' where \eqn{i} is the imaginary unit. Its characteristic function is
#' \eqn{e^{i m u} \big(\cos(b/2) / \cosh((a u - i b)/2)\big)^{2 d}}
#' (Schoutens, 2002, section 2), so that the distribution is infinitely
#' divisible. Its mean and variance are \eqn{m + a d \tan(b/2)} and
#' \eqn{a^2 d / (1 + \cos b)} (Grigoletto and Provasi, 2009, equation (1)).
#'
#' The Meixner process \eqn{X} is the Levy process associated with this
#' distribution: \eqn{X_0 = 0}, its increments are independent and
#' stationary, and \eqn{X_s \sim Meixner(a, b, d s, m s)} (Schoutens, 2002,
#' section 2, where the increment over \eqn{[s, s + t]} has the characteristic
#' function of \eqn{Meixner(a, b, d, m)} raised to the power \eqn{t}): an
#' increment over a duration
#' \eqn{\Delta t} has the density \eqn{f_{\Delta t}} obtained by replacing
#' \eqn{d} by \eqn{d \Delta t} and \eqn{m} by \eqn{m \Delta t} in \eqn{f}. The
#' mean and the variance of the process per unit of time are
#' \eqn{\kappa_1 = m + a d \tan(b/2)} and \eqn{\kappa_2 = a^2 d / (1 + \cos b)}.
#'
#' The data are the \eqn{n} increments \eqn{(\Delta t_j, \Delta x_j)},
#' \eqn{j = 1, \dots, n}, pooled over the trajectories, with
#' \eqn{T = \sum_{j=1}^{n} \Delta t_j} and \eqn{S = \sum_{j=1}^{n} \Delta x_j}
#' (see the section Notation of [levyfit-package]). To shorten the formulas, let
#' \deqn{s_j = d \Delta t_j, \qquad y_j = \frac{\Delta x_j - m \Delta t_j}{a}, \qquad j = 1, \dots, n;}
#' \eqn{y_j} is the standardized value \eqn{z_j = (x_j - m)/a} of Grigoletto and Provasi (2009)
#' (section 4.2) for an increment over \eqn{\Delta t_j}.
#'
#' @section Likelihood and score:
#' The increments being independent, the log-likelihood is
#' \deqn{\ell(a, b, d, m) = \sum_{j=1}^{n} \Big[2 s_j \ln(2 \cos(b/2)) - \ln(2 a \pi) - \ln \Gamma(2 s_j)
#'   + b y_j + \ln |\Gamma(s_j + i y_j)|^2\Big].}
#' With equal durations \eqn{\Delta t_j = 1}, it is \eqn{n} times the average
#' log-likelihood \eqn{\ell_n(\theta)} of Grigoletto and Provasi (2009, equation (6)); with unequal
#' durations, it is its extension to the process.
#'
#' The first derivatives (the score) are explicit (Grigoletto and Provasi, 2009, appendix A). Grigoletto and Provasi (2009)
#' write them with \eqn{\Psi_j^+ = \psi(d + i z_j) + \psi(d - i z_j)} and
#' \eqn{\Psi_j^- = \psi(d + i z_j) - \psi(d - i z_j)}, \eqn{\psi} being the
#' digamma function. Here, with \eqn{\psi_j = \psi(s_j + i y_j)} and its real
#' and imaginary parts \eqn{\mathrm{Re}\,\psi_j} and \eqn{\mathrm{Im}\,\psi_j},
#' these are \eqn{\Psi_j^+ = 2\, \mathrm{Re}\,\psi_j} and
#' \eqn{i\, \Psi_j^- = -2\, \mathrm{Im}\,\psi_j}, and the score is
#' \deqn{\frac{\partial \ell}{\partial a} = \frac{1}{a} \sum_{j=1}^{n} \big(-1 - b y_j + 2 y_j\, \mathrm{Im}\,\psi_j\big),}
#' \deqn{\frac{\partial \ell}{\partial b} = \sum_{j=1}^{n} y_j - \tan(b/2) \sum_{j=1}^{n} s_j,}
#' \deqn{\frac{\partial \ell}{\partial d} = \sum_{j=1}^{n} \Delta t_j \big(2 \ln(2 \cos(b/2)) - 2 \psi(2 s_j) + 2\, \mathrm{Re}\,\psi_j\big),}
#' \deqn{\frac{\partial \ell}{\partial m} = \frac{1}{a} \sum_{j=1}^{n} \Delta t_j \big(2\, \mathrm{Im}\,\psi_j - b\big).}
#'
#' @section The estimate of b:
#' Setting \eqn{\partial \ell / \partial b = 0} shows that the estimate of
#' \eqn{b} is a function of the estimates of the other parameters:
#' \deqn{\hat b = 2 \arctan \frac{S - \hat m T}{\hat a \hat d T},}
#' which is the lemma of Grigoletto and Provasi (2009, section 4.2),
#' \eqn{\hat b = 2 \tan^{-1}(\bar z / \hat d)}, written for increments over
#' unequal durations. It makes the expected total increment
#' \eqn{\kappa_1 T} equal to \eqn{S}.
#'
#' @section Method of moments:
#' With unequal durations (extension of Grigoletto and Provasi, 2009, section 4.1), the mean and the
#' higher cumulants of the process per unit of time are estimated by
#' \eqn{k_1 = S / T} and, with the residuals
#' \eqn{e_j = \Delta x_j - k_1 \Delta t_j},
#' \deqn{k_2 = \frac{1}{T} \sum_{j=1}^{n} e_j^2, \qquad
#'   k_3 = \frac{1}{T} \sum_{j=1}^{n} e_j^3, \qquad
#'   k_4 = \frac{1}{T} \Big(\sum_{j=1}^{n} e_j^4 - 3 k_2^2 \sum_{j=1}^{n} \Delta t_j^2\Big).}
#' Let \eqn{\gamma_1 = k_3 / k_2^{3/2}} and \eqn{\gamma_2 = k_4 / k_2^2} be
#' the skewness and the excess kurtosis of \eqn{X_1}. Since the skewness of
#' \eqn{Meixner(a, b, d, m)} is \eqn{\sin b / \sqrt{d (1 + \cos b)}} and its
#' kurtosis \eqn{3 + (2 - \cos b) / d} (Grigoletto and Provasi, 2009, equation (1)), that is
#' \eqn{\gamma_1^2 = (1 - \cos b) / d} and \eqn{\gamma_2 = (2 - \cos b) / d},
#' the moment estimates are (Grigoletto and Provasi, 2009, section 4.1, with the excess kurtosis)
#' \deqn{d = \frac{1}{\gamma_2 - \gamma_1^2}, \quad
#'   b = \mathrm{sign}(\gamma_1) \arccos(2 - d \gamma_2), \quad
#'   a = \sqrt{k_2 (1 + \cos b) / d}, \quad m = k_1 - a d \tan(b/2).}
#' They exist if and only if \eqn{d > 0} and \eqn{|2 - d \gamma_2| < 1}, that
#' is \eqn{\gamma_2 > \frac{3}{2} \gamma_1^2} and \eqn{\gamma_1 \ne 0}. With
#' equal durations \eqn{\Delta t}, they are the estimates of Grigoletto and Provasi (2009, section 4.1)
#' computed on the increments, with the variance divided by \eqn{n}, and with
#' \eqn{d} and \eqn{m} divided by \eqn{\Delta t}.
#'
#' @section Algorithm:
#' The likelihood is maximized as recommended by Grigoletto and Provasi
#' (2009, section 4.2): by a Newton-type algorithm that uses the first and
#' second derivatives of \eqn{\ell} (see Likelihood and score, and Hessian),
#' started at the moment estimates.
#' 1. Compute the moment estimates (see Method of moments). When they do not
#'    exist, there is no starting point: the function returns `NA` with a
#'    warning.
#' 2. From this starting point, maximize \eqn{\ell} with [stats::nlminb()],
#'    which uses the exact gradient and the exact Hessian. To stay in the
#'    parameter space \eqn{a > 0}, \eqn{d > 0}, \eqn{b \in (-\pi, \pi)}, the
#'    maximization is done in \eqn{(\ln a, b, \ln d, m)}, with \eqn{b} kept in
#'    \eqn{(-\pi, \pi)}.
#' 3. The solution is a local maximum when \eqn{\ell} has a zero gradient and a
#'    negative definite Hessian there. If [stats::nlminb()] does not converge,
#'    or if the Hessian at the solution is not negative definite, the function
#'    returns `NA` with a warning.
#'
#' Grigoletto and Provasi (2009, section 4.2) find that the convergence of the estimator to its
#' asymptotic normal distribution is very slow, with a speed inversely related
#' to \eqn{d}, and that inference based on asymptotic normality must be carefully
#' evaluated even for medium samples.
#'
#' @section Hessian:
#' The second derivatives of \eqn{\ell} are explicit (Grigoletto and Provasi, 2009, appendix A, where
#' they are written with \eqn{\Upsilon_j^\pm = \psi'(d + i z_j) \pm \psi'(d - i z_j)},
#' \eqn{\psi'} being the trigamma function). Here, with
#' \eqn{\psi'_j = \psi'(s_j + i y_j)},
#' \deqn{G_j = b - 2\, \mathrm{Im}\,\psi_j, \qquad R_j = \mathrm{Re}\,\psi'_j, \qquad I_j = \mathrm{Im}\,\psi'_j,
#'   \qquad j = 1, \dots, n,}
#' so that \eqn{\Upsilon_j^+ = 2 R_j} and \eqn{\Upsilon_j^- = 2 i I_j}, they are
#' \deqn{\frac{\partial^2 \ell}{\partial a^2} = \frac{1}{a^2} \sum_{j=1}^{n} \big(1 - 2 R_j y_j^2 + 2 G_j y_j\big), \qquad
#'   \frac{\partial^2 \ell}{\partial a \partial b} = -\frac{1}{a} \sum_{j=1}^{n} y_j, \qquad
#'   \frac{\partial^2 \ell}{\partial a \partial d} = \frac{2}{a} \sum_{j=1}^{n} I_j y_j \Delta t_j,}
#' \deqn{\frac{\partial^2 \ell}{\partial a \partial m} = \frac{1}{a^2} \sum_{j=1}^{n} \Delta t_j \big(G_j - 2 R_j y_j\big), \qquad
#'   \frac{\partial^2 \ell}{\partial b^2} = -\frac{1}{1 + \cos b} \sum_{j=1}^{n} s_j,}
#' \deqn{\frac{\partial^2 \ell}{\partial b \partial d} = -T \tan(b/2), \qquad
#'   \frac{\partial^2 \ell}{\partial b \partial m} = -\frac{T}{a}, \qquad
#'   \frac{\partial^2 \ell}{\partial d^2} = \sum_{j=1}^{n} \Delta t_j^2 \big(2 R_j - 4 \psi'(2 s_j)\big),}
#' \deqn{\frac{\partial^2 \ell}{\partial d \partial m} = \frac{2}{a} \sum_{j=1}^{n} I_j \Delta t_j^2, \qquad
#'   \frac{\partial^2 \ell}{\partial m^2} = -\frac{2}{a^2} \sum_{j=1}^{n} R_j \Delta t_j^2.}
#' With `ci = TRUE`, the observed Fisher information is minus this Hessian at
#' the estimate, as proposed by Grigoletto and Provasi (2009, section 4.2), instead of the numerical
#' Hessian used by the other models.
#'
#' @inheritParams wiener_mle
#' @param x Degradation levels of one trajectory; or a list of data frames,
#'   one per trajectory, whose first column contains the inspection times and
#'   second column the degradation levels.
#'   The increments can also be given directly with [increments_data()].
#' @param t Inspection times of the trajectory, strictly increasing. Not used
#'   when `x` is a list.
#'
#' @return A named vector with the estimates `a`, `b`, `d` and `m`.
#'
#'   With `ci = TRUE`, a data frame with one row per parameter and the columns
#'   `parameter`, `estimate`, `lower` and `upper` (bounds of the confidence
#'   interval). When the estimator does not exist, the estimates (and bounds)
#'   are `NA`, with a warning.
#'
#' @references
#' Grigoletto, M. and Provasi, C. (2009). Simulation and estimation of the
#' Meixner distribution. *Communications in Statistics - Simulation and
#' Computation*, 38(1), 58-77. \doi{10.1080/03610910802395679}
#'
#' Schoutens, W. (2002). The Meixner process: theory and applications in
#' finance. EURANDOM report 2002-004. Section 2.
#'
#' @seealso [vg_mle()]
#'
#' @examplesIf requireNamespace("DegradSim", quietly = TRUE)
#' # 1. One trajectory with 5000 inspections up to 500 (step 0.1)
#' #    (a = 1, b = 0.3, d = 1.5, m = 0.2)
#' t <- DegradSim::inspection_times("regular", n = 5000, tmax = 500)$time
#' set.seed(1)
#' x <- DegradSim::sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' meixner_mle(x)
#'
#' # The same trajectory given as two vectors: degradation levels, then times
#' meixner_mle(x[, "X1"], t)
#'
#' # 2. The same duration with a step 1 (500 increments); comparison of the
#' #    mean and the variance of an increment over one time unit, true and
#' #    fitted
#' t1 <- DegradSim::inspection_times("regular", n = 500, tmax = 500)$time
#' set.seed(1)
#' x1 <- DegradSim::sim_meixner(t1, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' p <- meixner_mle(x1)
#' p
#' c(true = 0.2 + 1.5 * tan(0.15), fitted = p[["m"]] + p[["a"]] * p[["d"]] * tan(p[["b"]] / 2))
#' c(true = 1.5 / (2 * cos(0.15)^2), fitted = p[["a"]]^2 * p[["d"]] / (2 * cos(p[["b"]] / 2)^2))
#'
#' # 3. Twenty trajectories with the same inspection times (step 0.1)
#' t <- DegradSim::inspection_times("regular", n = 200, tmax = 20)$time
#' x <- DegradSim::sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2, n_traj = 20)
#' meixner_mle(x)
#'
#' # 4. Estimates with 95 % confidence intervals
#' meixner_mle(x, ci = TRUE)
#'
#' # 5. Ten increments over long intervals: when no acceptable maximum is
#' #    found, the function gives a warning and NA estimates
#' t10 <- DegradSim::inspection_times("regular", n = 10, tmax = 100)$time
#' y <- DegradSim::sim_meixner(t10, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' meixner_mle(y)
#' @export
meixner_mle <- function(x, t = NULL, ci = FALSE, level = 0.95) {
  estimate("meixner", x, t, ci, level)
}

meixner_est <- function(x, t = NULL) {
  d <- increments(x, t)
  dt <- d$dt
  dx <- d$dx
  if (length(dx) < 4) stop("at least four increments are needed to estimate four parameters")

  mom <- meixner_moments(dx, dt)
  if (is.null(mom)) no_mle("the moment estimates, used as starting point, do not exist")

  theta <- function(p) c(a = exp(p[1]), b = p[2], d = exp(p[3]), m = p[4])
  nll <- function(p) {
    q <- theta(p)
    v <- -meixner_loglik(dx, dt, q[["a"]], q[["b"]], q[["d"]], q[["m"]])
    if (is.finite(v)) v else Inf
  }
  ngr <- function(p) {
    q <- theta(p)
    s <- meixner_score(dx, dt, q[["a"]], q[["b"]], q[["d"]], q[["m"]])
    -c(q[["a"]] * s[["a"]], s[["b"]], q[["d"]] * s[["d"]], s[["m"]])
  }
  nhe <- function(p) {
    q <- theta(p)
    s <- meixner_score(dx, dt, q[["a"]], q[["b"]], q[["d"]], q[["m"]])
    j <- c(q[["a"]], 1, q[["d"]], 1)
    h <- meixner_hessian(dx, dt, q[["a"]], q[["b"]], q[["d"]], q[["m"]]) * outer(j, j)
    h[1, 1] <- h[1, 1] + q[["a"]] * s[["a"]]
    h[3, 3] <- h[3, 3] + q[["d"]] * s[["d"]]
    -h
  }

  start <- c(log(mom[["a"]]), mom[["b"]], log(mom[["d"]]), mom[["m"]])
  fit <- tryCatch(nlminb(start, nll, ngr, nhe, lower = c(-Inf, -pi, -Inf, -Inf), upper = c(Inf, pi, Inf, Inf),
                         control = list(eval.max = 1000, iter.max = 500)),
                  error = function(e) NULL)
  if (is.null(fit) || fit$convergence != 0) no_mle("the Newton-type algorithm did not converge")

  q <- theta(fit$par)
  h <- meixner_hessian(dx, dt, q[["a"]], q[["b"]], q[["d"]], q[["m"]])
  if (any(!is.finite(h)) || max(eigen(h, symmetric = TRUE, only.values = TRUE)$values) >= 0) {
    no_mle("the solution is not a local maximum (Hessian not negative definite)")
  }
  q
}

meixner_loglik <- function(dx, dt, a, b, d, m) {
  y <- (dx - m * dt) / a
  sum(2 * d * dt * log(2 * cos(b / 2)) - log(2 * a * pi) - lgamma(2 * d * dt) + b * y + log_abs_gamma2(d * dt, y))
}

meixner_score <- function(dx, dt, a, b, d, m) {
  y <- (dx - m * dt) / a
  ps <- complex_digamma(complex(real = d * dt, imaginary = y))
  c(a = sum(-1 / a - b * y / a + 2 * Im(ps) * y / a),
    b = sum(-d * dt * tan(b / 2) + y),
    d = sum(dt * (2 * log(2 * cos(b / 2)) - 2 * digamma(2 * d * dt) + 2 * Re(ps))),
    m = sum(dt / a * (2 * Im(ps) - b)))
}

log_abs_gamma2 <- function(x, y) {
  n <- max(length(x), length(y))
  2 * Re(gsl::lngamma_complex(rep_len(x, n), rep_len(y, n)))
}

complex_digamma <- function(z) {
  n <- pmax(0, ceiling(10 - Re(z)))
  s <- complex(length(z))
  for (k in seq_len(max(n, 0)) - 1) s <- s + ifelse(k < n, 1 / (z + k), 0)
  w <- z + n
  w2 <- 1 / w^2
  log(w) - 1 / (2 * w) - w2 * (1 / 12 - w2 * (1 / 120 - w2 * (1 / 252 - w2 * (1 / 240 - w2 / 132)))) - s
}

meixner_moments <- function(dx, dt) {
  Tt <- sum(dt)
  k1 <- sum(dx) / Tt
  r <- dx - k1 * dt
  k2 <- sum(r^2) / Tt
  k3 <- sum(r^3) / Tt
  k4 <- (sum(r^4) - 3 * k2^2 * sum(dt^2)) / Tt
  g1 <- k3 / k2^1.5
  g2 <- k4 / k2^2
  d <- 1 / (g2 - g1^2)
  cb <- 2 - d * g2
  if (!is.finite(d) || d <= 0 || !is.finite(cb) || abs(cb) >= 1) return(NULL)
  b <- sign(g1) * acos(cb)
  a <- sqrt(k2 * (1 + cb) / d)
  c(a = a, b = b, d = d, m = k1 - a * d * tan(b / 2))
}

complex_trigamma <- function(z) {
  n <- pmax(0, ceiling(10 - Re(z)))
  s <- complex(length(z))
  for (k in seq_len(max(n, 0)) - 1) s <- s + ifelse(k < n, 1 / (z + k)^2, 0)
  w <- z + n
  w2 <- 1 / w^2
  s + 1 / w + w2 / 2 + w2 / w * (1 / 6 - w2 * (1 / 30 - w2 * (1 / 42 - w2 * (1 / 30 - w2 * 5 / 66))))
}

meixner_hessian <- function(dx, dt, a, b, d, m) {
  y <- (dx - m * dt) / a
  z <- complex(real = d * dt, imaginary = y)
  ps <- complex_digamma(z)
  p1 <- complex_trigamma(z)
  Gy <- -2 * Im(ps)
  Gyy <- -2 * Re(p1)
  Guy <- -2 * Im(p1)
  Guu <- 2 * Re(p1)
  ya <- -y / a
  ym <- -dt / a
  h <- matrix(0, 4, 4, dimnames = list(c("a", "b", "d", "m"), c("a", "b", "d", "m")))
  h["a", "a"] <- sum(1 / a^2 + Gyy * ya^2 + (b + Gy) * 2 * y / a^2)
  h["a", "b"] <- sum(ya)
  h["a", "d"] <- sum(Guy * dt * ya)
  h["a", "m"] <- sum(Gyy * ya * ym + (b + Gy) * dt / a^2)
  h["b", "b"] <- -sum(d * dt) / (2 * cos(b / 2)^2)
  h["b", "d"] <- -sum(dt) * tan(b / 2)
  h["b", "m"] <- sum(ym)
  h["d", "d"] <- sum(dt^2 * (Guu - 4 * trigamma(2 * d * dt)))
  h["d", "m"] <- sum(Guy * dt * ym)
  h["m", "m"] <- sum(Gyy * ym^2)
  h[lower.tri(h)] <- t(h)[lower.tri(h)]
  h
}
