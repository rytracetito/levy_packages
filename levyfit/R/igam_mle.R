#' Maximum likelihood estimation for the inverse gamma process
#'
#' Computes the maximum likelihood estimator of the parameters
#' \eqn{(\theta, \alpha, \beta)} of the inverse gamma process from one or
#' several trajectories observed at inspection times. The exact likelihood, which accounts
#' for equal observations, is computed by forward filtering over the hidden
#' positions of the underlying gamma process, then maximized numerically.
#'
#' @section Model and data:
#' Let \eqn{T(w)} be a gamma process in the level \eqn{w}, with increments
#' \eqn{T(w_2) - T(w_1) \sim Gamma(\eta(w_2) - \eta(w_1))} with scale
#' \eqn{\theta}, and shape function \eqn{\eta(w) = (w / \alpha)^\beta}. The
#' inverse gamma process is \eqn{W(t) = \inf\{w : T(w) > t\}}. It is increasing
#' but it is neither a Lévy process nor a Markov process: the likelihood is not
#' a product of densities of increments, and two consecutive observations can
#' be equal with positive probability. Writing \eqn{T(w) = \theta G(\eta(w))},
#' with \eqn{G} a standard gamma process, gives
#' \eqn{W(t) = \eta^{-1}(L(t / \theta))} with \eqn{L(z) = \inf\{u : G(u) > z\}}.
#' Each trajectory is observed at its own times
#' \eqn{0 = t_0 < t_1 < \dots < t_N}, with
#' \eqn{0 = w_0 < w_1 \le \dots \le w_N}. The trajectories are independent and
#' have the same law.
#'
#' @section Reduced data:
#' Let \eqn{z_j = t_j / \theta} and \eqn{u_j = \eta(w_j)}. Equal observations
#' are grouped into plateaus \eqn{k = 1, \dots, K}, with distinct levels
#' \eqn{\bar w_k}, first and last reduced times \eqn{\underline z_k} and
#' \eqn{\overline z_k}, and \eqn{s_k = \eta(\bar w_k) - \eta(\bar w_{k-1})}.
#' The parameters are identifiable: \eqn{E[T(w)] = \theta (w / \alpha)^\beta}
#' and \eqn{Var(T(w)) = \theta E[T(w)]}.
#'
#' @section Exact likelihood:
#' Let \eqn{g_s} be the density of \eqn{Gamma(s, 1)}, \eqn{E_1} the
#' exponential integral, and
#' \deqn{\kappa(s; r, d) = \int_0^r g_s(y) \frac{e^{-(d - y)}}{d - y} dy, \qquad
#'   \Pi(s; r, R) = \int_0^r g_s(y) E_1(R - y) dy.}
#' \eqn{\kappa} is the joint density of the passage level above \eqn{r} and of
#' the position after the jump for a standard gamma process, and
#' \eqn{\Pi(s; r, R)} the density of the passage level above \eqn{r} with a
#' jump that also exceeds \eqn{R}. The density of the distinct observed levels
#' is
#' \deqn{p(\bar w_1, \dots, \bar w_K) = \Big[\prod_k \eta'(\bar w_k)\Big] \Lambda,}
#' \deqn{\Lambda = \int \prod_{k=1}^{K-1} \kappa(s_k; \underline z_k - x_{k-1}, x_k - x_{k-1})\,
#'   \Pi(s_K; \underline z_K - x_{K-1}, \overline z_K - x_{K-1})\, dx_1 \cdots dx_{K-1},}
#' with \eqn{x_0 = 0} and \eqn{\overline z_k < x_k < \underline z_{k+1}}: the
#' \eqn{x_k} are the hidden positions of \eqn{G} after each jump. The
#' probability of equal observations is included through the plateaus. Since
#' the process is not a Lévy process, the increments of different trajectories
#' cannot be pooled: the log-likelihood of several trajectories is the sum of
#' the logarithms of their densities.
#'
#' @section Computation:
#' \eqn{\Lambda} is computed by forward filtering: each integral is computed by
#' the trapezoidal rule with 51 nodes, after the change of variable
#' \eqn{x = \overline z_k + (\underline z_{k+1} - \overline z_k) \phi(\omega)},
#' \eqn{\phi(\omega) = \omega^3 (10 - 15 \omega + 6 \omega^2)}, whose
#' derivative vanishes at both ends and smooths the two integrable
#' singularities of the integrand. The total mass of the density is then 1 up
#' to about \eqn{4 \times 10^{-5}}.
#'
#' @section Estimation:
#' 1. Starting point: regression of \eqn{\ln t_j} on \eqn{\ln w_j} (slope
#'    \eqn{\beta_0}, intercept \eqn{c_0}), \eqn{\theta_0} = mean of
#'    \eqn{(t - \hat t)^2 / \hat t}, \eqn{\alpha_0 = \exp((\ln \theta_0 - c_0) / \beta_0)}.
#' 2. L-BFGS-B on \eqn{(\ln \theta, \ln \alpha, \ln \beta)}, with bounds
#'    \eqn{\pm 4} around the start for \eqn{\theta} and \eqn{\alpha} and
#'    \eqn{\beta \in [0.05, 10]}.
#'
#' Neither the existence nor the uniqueness of the MLE is proven. The
#' likelihood of Guida and Pulcini (2013), a product of conditional densities,
#' is \eqn{-\infty} as soon as two observations are equal and is biased
#' otherwise; it is not used here.
#'
#' @param x Degradation levels of one trajectory, starting at 0 at time 0,
#'   non-decreasing, with \eqn{x_1 > 0}; or a list of data frames, one per
#'   trajectory, whose first column contains the inspection times and second
#'   column the degradation levels.
#' @param t Inspection times of the trajectory, strictly increasing, starting
#'   at 0. Not used when `x` is a list.
#'
#' @return A named vector with the estimates `theta`, `alpha` and `beta`.
#'
#' @references
#' Guida, M. and Pulcini, G. (2013). The inverse Gamma process: a family of
#' continuous stochastic models for describing state-dependent deterioration
#' phenomena. *Reliability Engineering & System Safety*, 120, 72-79.
#' \doi{10.1016/j.ress.2013.03.013}
#'
#' @seealso [gamma_mle()], [inspection_times()]
#'
#' @examplesIf requireNamespace("levysim", quietly = TRUE)
#' \donttest{
#' # 1. One trajectory inspected every 2 time units up to 20
#' #    (theta = 1, alpha = 1, beta = 1.5)
#' t <- seq(0, 20, by = 2)
#' set.seed(1)
#' x <- levysim::sim_igam(t, theta = 1, alpha = 1, beta = 1.5)
#' igam_mle(x)
#'
#' # The same trajectory given as two vectors: degradation levels, then times
#' igam_mle(x[[1]]$X, t)
#'
#' # 2. Ten trajectories with the same inspection times
#' x <- levysim::sim_igam(t, theta = 1, alpha = 1, beta = 1.5, n_traj = 10)
#' igam_mle(x)
#'
#' # 3. Estimates with 95 % confidence intervals
#' confint_mle(x, model = "igam")
#' }
#' @export
igam_mle <- function(x, t = NULL) {
  traj <- as_trajectories(x, t)
  for (d in traj) {
    if (d$t[1] != 0 || d$x[1] != 0) stop("each trajectory must start at t = 0 with x = 0")
    if (any(diff(d$x) < 0) || d$x[2] <= 0) stop("levels must be non-decreasing, with x[2] > 0")
  }

  W <- lapply(traj, function(d) d$x)
  tt <- lapply(traj, function(d) d$t)
  nll <- function(p) {
    v <- -igam_loglik(W, tt, exp(p[1]), exp(p[2]), exp(p[3]))
    if (is.finite(v)) v else 1e10
  }
  st <- igam_start(W, tt)
  lower <- c(st[1:2] - 4, log(0.05))
  upper <- c(st[1:2] + 4, log(10))
  st[3] <- min(max(st[3], lower[3] + 0.1), upper[3] - 0.1)
  opt <- optim(st, nll, method = "L-BFGS-B", lower = lower, upper = upper, control = list(factr = 1e7))
  if (opt$convergence != 0) warning("optim did not converge (code ", opt$convergence, ")")
  if (any(abs(opt$par - lower) < 1e-6 | abs(opt$par - upper) < 1e-6)) warning("the estimate is on a bound")

  c(theta = exp(opt$par[1]), alpha = exp(opt$par[2]), beta = exp(opt$par[3]))
}

igam_loglik <- function(W, t, theta, alpha, beta, M = 50) {
  sum(vapply(seq_along(W), function(i) {
    u <- (W[[i]] / alpha)^beta
    if (max(u) > 1e3 * (max(t[[i]]) / theta + 10)) return(-Inf)
    wd <- unique(W[[i]][-1])
    igam_loglik_red(t[[i]] / theta, u, M) + sum(log(beta / alpha) + (beta - 1) * log(wd / alpha))
  }, 0))
}

igam_loglik_red <- function(z, a, M = 50) {
  bl <- make_blocks(z, a)
  K <- nrow(bl) - 1
  s <- diff(bl$a)
  if (K == 1) return(log_pass_beyond(s[1], bl$z_first[2], bl$z_last[2]))
  om <- seq(0, 1, length.out = M + 1)
  phi <- om^3 * (10 - 15 * om + 6 * om^2)
  lw <- log(30 * om^2 * (1 - om)^2) + log(c(0.5, rep(1, M - 1), 0.5) / M)
  lc <- 0
  for (k in 1:(K - 1)) {
    lo <- bl$z_last[k + 1]
    D <- bl$z_first[k + 2] - lo
    x <- lo + D * phi
    la <- if (k == 1) log_kappa(s[1], bl$z_first[2], x) else {
      L <- outer(lp, rep(1, M + 1)) +
        matrix(log_kappa(s[k], rep(bl$z_first[k + 1] - xp, M + 1), rep(x, each = M + 1) - rep(xp, M + 1)), M + 1)
      mx <- apply(L, 2, max)
      mx[!is.finite(mx)] <- 0
      mx + log(colSums(exp(L - rep(mx, each = M + 1))))
    }
    v <- la + log(D) + lw
    v[!is.finite(v)] <- -Inf
    lc <- lc + max(v)
    lp <- v - max(v)
    xp <- x
  }
  lf <- lp + log_pass_beyond(s[K], bl$z_first[K + 1] - xp, bl$z_last[K + 1] - xp)
  mf <- max(lf[is.finite(lf)])
  lc + mf + log(sum(exp(lf - mf)))
}

igam_start <- function(W, t) {
  lw <- unlist(lapply(W, function(w) log(w[-1])))
  lt <- unlist(lapply(t, function(s) log(s[-1])))
  cf <- coef(lm(lt ~ lw))
  b0 <- cf[[2]]
  f0 <- exp(cf[[1]] + b0 * lw)
  th0 <- mean((exp(lt) - f0)^2 / f0)
  c(log(th0), (log(th0) - cf[[1]]) / b0, log(b0))
}

make_blocks <- function(z, a) {
  k <- cumsum(c(TRUE, diff(a) > 0))
  data.frame(a = tapply(a, k, `[`, 1), z_first = tapply(z, k, min), z_last = tapply(z, k, max))
}

gl_nodes <- function(K) {
  i <- seq_len(K - 1)
  bt <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, K, K)
  J[cbind(i, i + 1)] <- bt
  J[cbind(i + 1, i)] <- bt
  e <- eigen(J, symmetric = TRUE)
  list(x = (e$values + 1) / 2, w = e$vectors[1, ]^2)
}

GL20 <- gl_nodes(20)

logBs0 <- function(x, s) {
  n <- max(length(x), length(s))
  x <- rep_len(x, n)
  s <- rep_len(s, n)
  out <- rep(-Inf, n)
  ser <- function(xv, sv) {
    acc <- 1 / sv
    p <- rep(1, length(xv))
    for (k in seq_len(ceiling(-38 / log(max(xv))) + 2)) {
      p <- p * xv
      acc <- acc + p / (sv + k)
    }
    sv * log(xv) + log(acc)
  }
  lo <- x > 0 & x <= 0.6
  if (any(lo)) out[lo] <- ser(x[lo], s[lo])
  hi <- x > 0.6
  if (any(hi)) {
    xh <- x[hi]
    sh <- s[hi]
    v <- outer(1 - xh, GL20$x) + xh
    f <- (v^(sh - 1) - 1) / (1 - v)
    bad <- !is.finite(f)
    f[bad] <- (1 - sh)[row(f)[bad]]
    r <- -log1p(-xh) - digamma(sh) + digamma(1) - as.vector(f %*% GL20$w) * (1 - xh)
    cc <- r < 1e-6 * (-log1p(-xh) + abs(digamma(sh)) + 1)
    r <- log(pmax(r, 0))
    if (any(cc)) r[cc] <- ser(xh[cc], sh[cc])
    out[hi] <- r
  }
  out
}

E1 <- function(x) {
  out <- numeric(length(x))
  sm <- x <= 2
  if (any(sm)) {
    k <- 1:40
    out[sm] <- -0.5772156649015329 - log(x[sm]) - vapply(x[sm], function(t) sum((-t)^k / (k * factorial(k))), 0)
  }
  if (any(!sm)) {
    out[!sm] <- vapply(x[!sm], function(t) {
      b <- t + 1; cc <- 1e300; d <- 1 / b; h <- d
      for (i in 1:300) {
        an <- -i * i; b <- b + 2; d <- 1 / (an * d + b); cc <- b + an / cc
        del <- cc * d; h <- h * del
        if (abs(del - 1) < 1e-15) break
      }
      h * exp(-t)
    }, 0)
  }
  out
}

log_dz <- function(u, z, h = 1e-3) {
  n <- max(length(u), length(z))
  u <- rep_len(u, n)
  z <- rep_len(z, n)
  b <- log(u)
  lq <- function(b, z) pgamma(z, shape = exp(b), lower.tail = FALSE, log.p = TRUE)
  lp <- function(b, z) pgamma(z, shape = exp(b), lower.tail = TRUE, log.p = TRUE)
  der <- function(f, b, z) (-f(b + 2 * h, z) + 8 * f(b + h, z) - 8 * f(b - h, z) + f(b - 2 * h, z)) / (12 * h)
  vq <- lq(b, z)
  up <- vq < log(0.5)
  out <- numeric(n)
  if (any(up)) out[up] <- vq[up] + log(der(lq, b[up], z[up])) - b[up]
  if (any(!up)) out[!up] <- lp(b[!up], z[!up]) + log(-der(lp, b[!up], z[!up])) - b[!up]
  out
}

log_kappa <- function(s, r, d) {
  n <- max(length(s), length(r), length(d))
  s <- rep_len(s, n); r <- rep_len(r, n); d <- rep_len(d, n)
  out <- rep(-Inf, n)
  ok <- r > 0 & d > r
  if (any(ok)) out[ok] <- -d[ok] + (s[ok] - 1) * log(d[ok]) + logBs0(r[ok] / d[ok], s[ok]) - lgamma(s[ok])
  out
}

log_pass_beyond <- function(s, r, R) {
  n <- max(length(r), length(R))
  r <- rep_len(r, n); R <- rep_len(R, n)
  out <- rep(-Inf, n)
  eqr <- r > 0 & R - r <= 1e-12 * r
  if (any(eqr)) out[eqr] <- log_dz(s, r[eqr])
  gt <- r > 0 & !eqr
  if (any(gt)) {
    om <- GL20$x
    v <- 1 - (1 - om)^2
    out[gt] <- vapply(which(gt), function(i) {
      y <- r[i] * v^(1 / s)
      f <- exp(-y) * E1(R[i] - y) * 2 * (1 - om)
      s * log(r[i]) - log(s) - lgamma(s) + log(sum(f * GL20$w))
    }, 0)
  }
  out
}
