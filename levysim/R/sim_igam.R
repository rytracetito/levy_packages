#' Simulation of the inverse gamma process
#'
#' Simulates exactly one or several trajectories of the inverse gamma process
#' at given inspection times.
#'
#' @section Model:
#' Let \eqn{\eta} be continuous, strictly increasing, with \eqn{\eta(0) = 0},
#' and \eqn{\theta > 0}. Let \eqn{T(w)} be a gamma process indexed by the
#' degradation level \eqn{w}, with
#' \eqn{T(w_2) - T(w_1) \sim Ga(\eta(w_2) - \eta(w_1))} with scale
#' \eqn{\theta}: \eqn{T(w)} is the time needed to reach the level \eqn{w}. The
#' inverse gamma process is
#' \deqn{W(t) = \inf\{w \ge 0 : T(w) > t\}.}
#' Writing \eqn{T(w) = \theta G(\eta(w))} with \eqn{G} a standard gamma
#' process gives \eqn{W(t) = \eta^{-1}(L(t / \theta))}, with
#' \eqn{L(z) = \inf\{u : G(u) > z\}}, and the marginal distribution
#' \eqn{P(W(t) \le w) = Q(\eta(w), t / \theta)}, where \eqn{Q} is the
#' regularized upper incomplete gamma function. \eqn{W} is continuous,
#' increasing and constant on long plateaus: two inspections can give the same
#' value with positive probability. It is neither a Lévy process nor a Markov
#' process. By default \eqn{\eta(w) = (w / \alpha)^\beta}.
#'
#' @section Exact simulation:
#' The increments are not independent. In reduced scale \eqn{z_j = t_j / \theta},
#' the algorithm follows the position of \eqn{G} just after the last jump that
#' crossed a threshold; after such a jump, \eqn{G} starts afresh (strong Markov
#' property). For a standard gamma process starting at 0, the passage level
#' \eqn{s = L(r)} above \eqn{r} satisfies \eqn{P(L(r) \le s) = Q(s, r)}; given
#' \eqn{s}, the position before the jump \eqn{y \in (0, r)} has a density
#' proportional to \eqn{g_s(y) E_1(r - y)} (\eqn{g_s} the density of
#' \eqn{Ga(s, 1)}, \eqn{E_1} the exponential integral); given \eqn{y}, the jump
#' \eqn{J > r - y} satisfies \eqn{P(J > j) = E_1(j) / E_1(r - y)}.
#' 1. Start with level \eqn{a = 0} and position \eqn{pos = 0}.
#' 2. For each inspection \eqn{j}: if \eqn{pos > z_j}, the previous jump
#'    covers \eqn{z_j} and \eqn{a_j = a_{j-1}} (equal observations);
#'    otherwise, with \eqn{r = z_j - pos}, draw \eqn{s} by inversion of
#'    \eqn{Q(s, r)}, then \eqn{y}, then \eqn{J}, and set
#'    \eqn{a_j = a_{j-1} + s}, \eqn{pos = pos + y + J}.
#' 3. Return \eqn{W(t_j) = \eta^{-1}(a_j)}.
#'
#' The algorithm is exact up to numerical errors (quadrature on 401 nodes,
#' bisection, root finding).
#'
#' @inheritParams sim_gamma
#' @param theta Scale of the gamma process in level, \eqn{\theta > 0}: the
#'   larger \eqn{\theta}, the slower the degradation.
#' @param alpha,beta Parameters of the default shape function
#'   \eqn{\eta(w) = (w / \alpha)^\beta}, \eqn{\alpha > 0}, \eqn{\beta > 0}
#'   (default 1 for both).
#' @param eta_inv Optional: an R function giving the inverse
#'   \eqn{\eta^{-1}(u)} of another shape function, for instance
#'   `function(u) log(1 + u)` for \eqn{\eta(w) = e^w - 1}. It then replaces
#'   \eqn{(w / \alpha)^\beta}, and `alpha` and `beta` are ignored.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Guida, M. and Pulcini, G. (2013). The inverse Gamma process: a family of
#' continuous stochastic models for describing state-dependent deterioration
#' phenomena. *Reliability Engineering & System Safety*, 120, 72-79.
#' \doi{10.1016/j.ress.2013.03.013}
#'
#' @seealso [plot_traj()], [traj_df()], [inspection_times()]
#'
#' @examples
#' # One trajectory inspected every time unit up to 10
#' t <- seq(0, 10, by = 1)
#' set.seed(1)
#' x <- sim_igam(t, theta = 1, alpha = 1, beta = 1.5)
#' x
#' plot_traj(x)
#'
#' # Five trajectories with the same inspection times: note the plateaus
#' x <- sim_igam(t, theta = 1, alpha = 1, beta = 1.5, n_traj = 5)
#' plot_traj(x)
#'
#' # Irregular inspections, more and more frequent over time
#' plan <- inspection_times("increasing", n = 10, n_traj = 3, tmax = 10)
#' x <- sim_igam(plan, theta = 1, alpha = 1, beta = 1.5)
#' plot_traj(x)
#'
#' # Different inspection times for each trajectory
#' x <- sim_igam(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), theta = 1, alpha = 1, beta = 1.5)
#' x[[2]]
#'
#' # Another shape function: eta(w) = exp(w) - 1, given by its inverse
#' x <- sim_igam(t, theta = 1, eta_inv = function(u) log(1 + u), n_traj = 5)
#' plot_traj(x, "Inverse gamma process with eta(w) = exp(w) - 1")
#'
#' # All the simulated data in a single data frame
#' head(traj_df(x))
#' @export
sim_igam <- function(t, theta, alpha = 1, beta = 1, eta_inv = NULL, n_traj = 1) {
  if (theta <= 0 || alpha <= 0 || beta <= 0) stop("theta, alpha and beta must be positive")
  if (is.null(eta_inv)) eta_inv <- function(u) alpha * u^(1 / beta)
  one_path <- function(t) {
    z <- t / theta
    a <- numeric(length(z))
    pos <- 0
    for (j in seq_along(z)[-1]) {
      if (pos > z[j]) {
        a[j] <- a[j - 1]
        next
      }
      p <- r_pass(z[j] - pos)
      a[j] <- a[j - 1] + p[["s"]]
      pos <- pos + p[["x"]]
    }
    eta_inv(a)
  }
  simulate_paths(t, n_traj, one_path, "inverse gamma", c(theta = theta, alpha = alpha, beta = beta))
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

rtab <- function(xg, q) {
  q[!is.finite(q)] <- 0
  M <- length(xg) - 1
  h <- diff(xg)
  m <- h * (q[-1] + q[-(M + 1)]) / 2
  j <- findInterval(stats::runif(1) * sum(m), c(0, cumsum(m)), rightmost.closed = TRUE)
  v <- stats::runif(1)
  q0 <- q[j]; q1 <- q[j + 1]
  xg[j] + h[j] * (if (q0 + q1 == 0) v else v * (q0 + q1) / (q0 + sqrt(q0^2 + v * (q1^2 - q0^2))))
}

rgrid <- function(lw, M) {
  lw[!is.finite(lw)] <- -Inf
  rtab(seq(0, 1, length.out = M + 1), exp(lw - max(lw)))
}

r_under <- function(s, r, d = NA, R = NA, M = 400) {
  om <- seq(0, 1, length.out = M + 1)
  v <- 1 - (1 - om)^2
  u <- r * v^(1 / s)
  lw <- if (is.na(d)) -u + log(E1(pmax(R - u, 1e-300))) else -log(d - u)
  vv <- 1 - (1 - rgrid(lw + log(2 * (1 - om)), M))^2
  r * vv^(1 / s)
}

rincr_eta <- function(n, z) {
  u <- stats::runif(n)
  lo <- rep(-40, n)
  hi <- rep(log(z + 60 * sqrt(z + 1) + 60), n)
  for (it in 1:80) {
    mid <- (lo + hi) / 2
    up <- stats::pgamma(z, shape = exp(mid), lower.tail = FALSE) >= u
    hi[up] <- mid[up]
    lo[!up] <- mid[!up]
  }
  exp((lo + hi) / 2)
}

r_pass <- function(r) {
  s <- rincr_eta(1, r)
  y <- r_under(s, r, R = r)
  lt <- log(E1(r - y)) + log(stats::runif(1))
  j <- stats::uniroot(function(j) log(E1(j)) - lt, c(r - y, r - y + 50 + 10 * (r - y)), tol = 1e-12)$root
  c(s = s, x = y + j)
}
