#' Maximum likelihood estimation for the Meixner process
#'
#' Computes the maximum likelihood estimator of the four parameters
#' \eqn{(a, b, d, m)} of the Meixner process from one or several trajectories observed at
#' possibly irregular inspection times. The likelihood is unbounded, so the
#' estimator is defined as the interior local maximum of largest likelihood.
#'
#' @section Model and data:
#' For \eqn{a > 0}, \eqn{-\pi < b < \pi}, \eqn{d > 0} and real \eqn{m}, the
#' \eqn{Meixner(a, b, d, m)} distribution has density
#' \deqn{f(x) = \frac{(2 \cos(b/2))^{2d}}{2 a \pi\, \Gamma(2d)}
#'   \exp\Big(\frac{b (x - m)}{a}\Big) \Big|\Gamma\Big(d + i \frac{x - m}{a}\Big)\Big|^2, \qquad x \in R.}
#' The Meixner process is the Lévy process such that
#' \eqn{X_s \sim Meixner(a, b, d s, m s)}. Its mean and variance per unit of
#' time are \eqn{\kappa_1 = m + a d \tan(b/2)} and
#' \eqn{\kappa_2 = a^2 d / (2 \cos^2(b/2))}. The trajectories are independent
#' and have the same law; each one is observed at its own inspection times. Let
#' \eqn{(\Delta t_j, \Delta x_j)}, \eqn{j = 1, \dots, N}, be the durations and
#' degradation increments between two consecutive inspections of a same
#' trajectory, gathered over all the trajectories,
#' \eqn{T = \sum_j \Delta t_j}, \eqn{S = \sum_j \Delta x_j},
#' \eqn{s_j = d \Delta t_j} and \eqn{y_j = (\Delta x_j - m \Delta t_j) / a}.
#'
#' @section Likelihood and score:
#' The increments being independent,
#' \deqn{\ell(a, b, d, m) = \sum_j \big[2 s_j \ln(2 \cos(b/2)) - \ln(2 a \pi) - \ln \Gamma(2 s_j)
#'   + b y_j + \ln |\Gamma(s_j + i y_j)|^2\big].}
#' With \eqn{\psi_j = \psi(s_j + i y_j)} (digamma function of a complex
#' argument), the score is explicit:
#' \eqn{\partial_a \ell = \sum_j (-1 - b y_j + 2 y_j \Im \psi_j) / a},
#' \eqn{\partial_b \ell = \sum_j (y_j - s_j \tan(b/2))},
#' \eqn{\partial_d \ell = \sum_j \Delta t_j (2 \ln(2 \cos(b/2)) - 2 \psi(2 s_j) + 2 \Re \psi_j)}
#' and \eqn{\partial_m \ell = \sum_j \Delta t_j (2 \Im \psi_j - b) / a}.
#'
#' @section Profile in b:
#' For fixed \eqn{(a, d, m)}, \eqn{\ell} is strictly concave in \eqn{b}, with
#' unique maximizer
#' \deqn{\hat b(a, d, m) = 2 \arctan \frac{S - m T}{a d T},}
#' which matches the expected total increment \eqn{\kappa_1 T} with \eqn{S}.
#' The maximization is therefore done in dimension 3, on the profile
#' log-likelihood \eqn{\ell_p(a, d, m) = \ell(a, \hat b(a, d, m), d, m)}.
#'
#' @section The likelihood is unbounded:
#' As \eqn{a \to 0} with \eqn{b = \pi - \lambda a}, the Meixner distribution
#' tends to a shifted gamma distribution. Placing the shift \eqn{m} on the
#' smallest observed speed \eqn{\Delta x_j / \Delta t_j} and taking
#' \eqn{2 d \Delta t_j < 1} makes \eqn{\ell \to +\infty}: the global MLE never
#' exists. The estimator is the interior local maximum (root of the score with
#' negative definite Hessian) of largest likelihood. Solutions with
#' \eqn{\pi - |b| < 0.1} (shifted gamma limit), with
#' \eqn{d \max_j \Delta t_j > 1000} (Gaussian limit) or on a bound are
#' discarded as degenerate.
#'
#' @section Algorithm:
#' 1. Compute the mean \eqn{c_1 = S / T} and the variance
#'    \eqn{c_2 = \sum_j (\Delta x_j - c_1 \Delta t_j)^2 / (T - \sum_j \Delta t_j^2 / T)}
#'    per unit of time.
#' 2. Draw `n_start` starting points: \eqn{b_0} in \eqn{(-1.5, 1.5)} and
#'    \eqn{d_0} around \eqn{1 / \overline{\Delta t}} (the first one is
#'    \eqn{b_0 = 0}, \eqn{d_0 = 1 / \overline{\Delta t}}), then
#'    \eqn{a_0 = \sqrt{2 c_2 \cos^2(b_0/2) / d_0}} and
#'    \eqn{m_0 = c_1 - a_0 d_0 \tan(b_0/2)}, so that the starting point matches
#'    \eqn{c_1} and \eqn{c_2}.
#' 3. From each start, maximize \eqn{\ell_p} by L-BFGS-B on
#'    \eqn{(\ln a, \ln d, m)} with the exact gradient and wide bounds relative
#'    to the scale of the data.
#' 4. Discard the degenerate solutions and keep the one of largest likelihood;
#'    a warning is given if a degenerate solution had a larger likelihood.
#'
#' The parameters \eqn{a} and \eqn{d} are strongly correlated and \eqn{d} is
#' poorly identified: hundreds of increments are needed for a stable
#' estimation, and with few increments there is often no interior maximum, in
#' which case the function stops.
#'
#' @param x Degradation levels of one trajectory; or a list of data frames,
#'   one per trajectory, whose first column contains the inspection times and
#'   second column the degradation levels.
#' @param t Inspection times of the trajectory, strictly increasing. Not used
#'   when `x` is a list.
#' @param n_start Number of starting points.
#'
#' @return A named vector with the estimates `a`, `b`, `d` and `m`.
#'
#' @references
#' Schoutens, W. and Teugels, J. L. (1998). Lévy processes, polynomials and
#' martingales. *Communications in Statistics. Stochastic Models*, 14(1-2),
#' 335-349.
#'
#' Schoutens, W. (2002). The Meixner process: theory and applications in
#' finance. EURANDOM report 2002-004.
#'
#' @seealso [vg_mle()], [inspection_times()]
#'
#' @examplesIf requireNamespace("levysim", quietly = TRUE)
#' # The Meixner model needs many small increments: over long intervals its
#' # increments are almost Gaussian and its parameters are hard to separate.
#'
#' # 1. One trajectory inspected every time unit up to 500
#' #    (a = 1, b = 0.3, d = 1.5, m = 0.2)
#' t <- seq(0, 500, by = 1)
#' set.seed(1)
#' x <- levysim::sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' meixner_mle(x)
#'
#' # The same trajectory given as two vectors: degradation levels, then times
#' meixner_mle(x[, "X1"], t)
#'
#' # 2. Twenty trajectories with the same inspection times
#' t <- seq(0, 20, by = 1)
#' x <- levysim::sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2, n_traj = 20)
#' meixner_mle(x)
#'
#' # 3. Estimates with 95 % confidence intervals
#' confint_mle(x, model = "meixner")
#'
#' # 4. Few increments over long intervals: the maximum is often at the
#' #    boundary (Gaussian or shifted gamma limit) and the function stops
#' y <- levysim::sim_meixner(seq(0, 100, by = 10), a = 1, b = 0.3, d = 1.5, m = 0.2)
#' try(meixner_mle(y))
#' @export
meixner_mle <- function(x, t = NULL, n_start = 5) {
  d <- increments(x, t)
  dt <- d$dt
  dx <- d$dx
  if (length(dx) < 4) stop("at least four increments are needed to estimate four parameters")
  Tt <- sum(dt)

  b_hat <- function(a, d, m) 2 * atan(sum((dx - m * dt) / a) / (d * Tt))
  nlp <- function(p) {
    a <- exp(p[1]); d <- exp(p[2])
    -meixner_loglik(dx, dt, a, b_hat(a, d, p[3]), d, p[3])
  }
  grp <- function(p) {
    a <- exp(p[1]); d <- exp(p[2])
    g <- meixner_score(dx, dt, a, b_hat(a, d, p[3]), d, p[3])
    -c(a * g[["a"]], d * g[["d"]], g[["m"]])
  }

  c1 <- sum(dx) / Tt
  s0 <- sqrt(sum((dx - c1 * dt)^2) / (Tt - sum(dt^2) / Tt))
  lower <- c(log(1e-3 * s0), log(1e-3 / max(dt)), c1 - 100 * s0)
  upper <- c(log(1e3 * s0), log(1e3 / min(dt)), c1 + 100 * s0)
  fits <- lapply(meixner_starts(dx, dt, n_start), function(st)
    tryCatch(optim(st, nlp, grp, method = "L-BFGS-B", lower = lower, upper = upper,
                   control = list(factr = 1e4, maxit = 2000)),
             error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) stop("the optimization failed from all starting points")

  deg <- vapply(fits, function(f) {
    a <- exp(f$par[1]); d <- exp(f$par[2])
    pi - abs(b_hat(a, d, f$par[3])) < 0.1 || d * max(dt) > 1e3 ||
      any(abs(f$par - lower) < 1e-6) || any(abs(f$par - upper) < 1e-6)
  }, TRUE)
  if (all(deg)) stop("no interior maximum: all solutions are at the boundary (shifted gamma or Gaussian limit)")
  ll <- -vapply(fits, function(f) f$value, 0)
  best <- fits[!deg][[which.max(ll[!deg])]]
  if (any(deg) && max(ll[deg]) > -best$value) {
    warning("a boundary solution has a larger likelihood than the interior maximum returned")
  }
  if (best$convergence != 0) warning("optim did not converge (code ", best$convergence, ")")

  a <- exp(best$par[1]); d <- exp(best$par[2]); m <- best$par[3]
  c(a = a, b = b_hat(a, d, m), d = d, m = m)
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

meixner_starts <- function(dx, dt, n_start) {
  Tt <- sum(dt)
  c1 <- sum(dx) / Tt
  c2 <- sum((dx - c1 * dt)^2) / (Tt - sum(dt^2) / Tt)
  b0 <- c(0, runif(n_start - 1, -1.5, 1.5))
  d0 <- c(1, exp(runif(n_start - 1, log(0.2), log(5)))) / mean(dt)
  lapply(seq_len(n_start), function(k) {
    a0 <- sqrt(2 * c2 * cos(b0[k] / 2)^2 / d0[k])
    c(log(a0), log(d0[k]), c1 - a0 * d0[k] * tan(b0[k] / 2))
  })
}
