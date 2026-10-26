#' Simulation of the Meixner process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' a homogeneous Meixner process, a Lévy process.
#'
#' @section Model:
#' The \eqn{Meixner(a, b, d, m)} distribution, with \eqn{a > 0},
#' \eqn{b \in ]-\pi, \pi[}, \eqn{d > 0} and \eqn{m} real, has density
#' \deqn{f(x) = \frac{(2 \cos(b/2))^{2d}}{2 a \pi\, \Gamma(2d)}
#'   \exp\Big(\frac{b (x - m)}{a}\Big) \Big|\Gamma\Big(d + i \frac{x - m}{a}\Big)\Big|^2,
#'   \qquad x \in \mathbb{R}.}
#' It is infinitely divisible: the sum of independent \eqn{Meixner(a, b, d_1, m_1)}
#' and \eqn{Meixner(a, b, d_2, m_2)} variables has the
#' \eqn{Meixner(a, b, d_1 + d_2, m_1 + m_2)} distribution.
#'
#' The Meixner process satisfies \eqn{X(0) = 0}, has independent increments,
#' and for \eqn{0 \le s < t},
#' \deqn{X(t) - X(s) \sim Meixner\big(a, b, d (t - s), m (t - s)\big).}
#' Hence \eqn{E[X(t)] = (m + a d \tan(b/2))\, t} and
#' \eqn{Var(X(t)) = a^2 d\, t / (2 \cos^2(b/2))}. It is a Lévy process whose
#' trajectories are not monotone, purely made of jumps, of infinite
#' variation: \eqn{a} is a scale parameter, \eqn{b} sets the asymmetry,
#' \eqn{d} the shape and \eqn{m} the location.
#'
#' @section Exact simulation at the inspection times (random walk approximation):
#' The increments over disjoint intervals are independent with known
#' distributions, so the values at the inspection times
#' \eqn{0 = t_0 < t_1 < \dots < t_n} are simulated by the random walk
#' approximation (RWA) of Kahle et al. (2016, section 2.5.1, written there for the gamma process):
#' 1. draw independently
#'    \eqn{\Delta X_i \sim Meixner(a, b, d \Delta t_i, m \Delta t_i)}, with
#'    \eqn{\Delta t_i = t_i - t_{i-1}};
#' 2. \eqn{X(t_0) = 0} and \eqn{X(t_i) = \sum_{j \le i} \Delta X_j}.
#'
#' There is no Meixner generator in `stats`: each increment is drawn by
#' numerical inversion of its distribution function (method PINV of the package
#' `Runuran`, built from the exact density, with an error in \eqn{u} of order
#' \eqn{10^{-10}}). One generator is built for each distinct time step
#' \eqn{\Delta t_i} and reused for all the trajectories. The simulation is
#' therefore exact up to this numerical error.
#'
#' @inheritSection sim_vg Conditions on the inspection times
#'
#' @inheritParams sim_gamma
#' @param a Scale parameter, \eqn{a > 0}.
#' @param b Asymmetry parameter, \eqn{-\pi < b < \pi}: \eqn{b = 0} gives
#'   symmetric increments, \eqn{b > 0} increments skewed to the right.
#' @param d Shape parameter per unit of time, \eqn{d > 0}.
#' @param m Location parameter per unit of time, real.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Kahle, W., Mercier, S. and Paroissin, C. (2016). *Degradation Processes in
#' Reliability*. ISTE and Wiley, London. ISBN 978-1-84821-888-8.
#'
#' Schoutens, W. (2003). *Lévy Processes in Finance: Pricing Financial
#' Derivatives*. Wiley, Chichester. ISBN 978-0-470-85156-2.
#'
#' Schoutens, W. and Teugels, J. L. (1998). Lévy processes, polynomials and
#' martingales. *Communications in Statistics. Stochastic Models*, 14(1-2),
#' 335-349.
#'
#' Hörmann, W. and Leydold, J. (2003). Continuous random variate generation by
#' fast numerical inversion. *ACM Transactions on Modeling and Computer
#' Simulation*, 13(4), 347-362. \doi{10.1145/945511.945517}
#'
#' @seealso [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()]
#'
#' @examples
#' # 1. One trajectory inspected every 0.5 time unit up to 10
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' x
#' plot_traj(x)
#'
#' # 2. Ten trajectories with the same inspection times
#' x <- sim_meixner(t, a = 1, b = 0.3, d = 1.5, m = 0.2, n_traj = 10)
#' plot_traj(x)
#'
#' # Same times for all the trajectories: a matrix with the column t,
#' # then one column per trajectory (X1, X2, ...)
#' head(x[, 1:4])
#' x[, "X2"]
#'
#' # 3. Irregular inspections, frequent in the middle of the period
#' plan <- inspection_times("inc-dec", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_meixner(plan, a = 1, b = 0.3, d = 1.5, m = 0.2)
#' plot_traj(x)
#'
#' # 4. Different inspection times for each trajectory: a list of data frames
#' x <- sim_meixner(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), a = 1, b = 0.3, d = 1.5, m = 0.2)
#' x[[2]]
#' x[[2]]$X
#'
#' # 5. Symmetric increments (b = 0) and no drift (m = 0)
#' x <- sim_meixner(t, a = 1, b = 0, d = 1.5, m = 0, n_traj = 10)
#' plot_traj(x, "Symmetric Meixner process")
#'
#' # 6. All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#' @export
sim_meixner <- function(t, a, b, d, m, n_traj = 1) {
  if (!is.numeric(a) || length(a) != 1 || !(a > 0)) stop("a must be a positive number")
  if (!is.numeric(b) || length(b) != 1 || !(abs(b) < pi)) stop("b must be in (-pi, pi)")
  if (!is.numeric(d) || length(d) != 1 || !(d > 0)) stop("d must be a positive number")
  if (!is.numeric(m) || length(m) != 1 || !is.finite(m)) stop("m must be a real number")
  generators <- new.env()
  draw <- function(s) {
    key <- format(s, digits = 15)
    if (is.null(generators[[key]])) {
      loi <- Runuran::udmeixner(alpha = a, beta = b, delta = d * s, mu = m * s)
      generators[[key]] <- Runuran::pinvd.new(loi)
    }
    Runuran::ur(generators[[key]], 1)
  }
  one_path <- function(t) c(0, cumsum(vapply(diff(t), draw, 0)))
  simulate_paths(t, n_traj, one_path, "Meixner", c(a = a, b = b, d = d, m = m))
}
