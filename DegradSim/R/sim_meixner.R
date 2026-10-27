#' Simulation of the Meixner process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' a Meixner process, a Lévy process.
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
#' The Meixner process satisfies \eqn{X_0 = 0}, has independent increments,
#' and for \eqn{0 \le s < t},
#' \deqn{X_t - X_s \sim Meixner\big(a, b, d (t - s), m (t - s)\big).}
#' Hence \eqn{E[X_t] = (m + a d \tan(b/2))\, t} and
#' \eqn{Var(X_t) = a^2 d\, t / (2 \cos^2(b/2))}. It is a Lévy process whose
#' trajectories are not monotone, purely made of jumps, of infinite
#' variation: \eqn{a} is a scale parameter, \eqn{b} sets the asymmetry,
#' \eqn{d} the shape and \eqn{m} the location.
#'
#' @section Exact simulation:
#' Random walk approximation, described in [DegradSim-package]: with
#' \eqn{\Delta t_i = t_i - t_{i-1}}, the increments
#' \deqn{\Delta X_i \sim Meixner(a, b, d \Delta t_i, m \Delta t_i)}
#' are drawn independently, then cumulated. There is no Meixner generator in `stats`: each increment is
#' drawn by numerical inversion of its distribution function (method PINV of
#' the package `Runuran`, built from the exact density, with an error in
#' \eqn{u} of order \eqn{10^{-10}}). One generator is built for each distinct
#' time step \eqn{\Delta t_i} and reused for all the trajectories. The
#' simulation is therefore exact up to this numerical error.
#'
#' @inheritParams sim_gamma
#' @param a Scale parameter, \eqn{a > 0}.
#' @param b Asymmetry parameter, \eqn{b \in ]-\pi, \pi[}: \eqn{b = 0} gives
#'   symmetric increments, \eqn{b > 0} increments skewed to the right.
#' @param d Shape parameter per unit of time, \eqn{d > 0}.
#' @param m Location parameter per unit of time, real.
#'
#' @inherit sim_gamma return
#'
#' @references
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
#' # 1. One trajectory with 200 regular inspections up to 10 (step 0.05)
#' t <- inspection_times("regular", n = 200, tmax = 10)$time
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
#' plan <- inspection_times("inc-dec", n = 50, n_traj = 5, tmax = 10)
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
#'
#' # 7. Histogram of the increments, with the density of the increments in red.
#' #    Constant step h = 1: the increments are i.i.d.
#' #    with distribution Meixner(a, b, d h, m h), density from Runuran::ud()
#' t7 <- inspection_times("regular", n = 50, tmax = 50)$time
#' x <- sim_meixner(t7, a = 1, b = 0.3, d = 1.5, m = 0.2, n_traj = 20)
#' dx <- traj_df(x, increments = TRUE)$dX
#' dx <- dx[!is.na(dx)]
#' law <- Runuran::udmeixner(alpha = 1, beta = 0.3, delta = 1.5, mu = 0.2)
#' ggplot2::ggplot(data.frame(dx), ggplot2::aes(dx)) +
#'   ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
#'                           bins = 30, boundary = 0, fill = "grey80", colour = "grey40") +
#'   ggplot2::stat_function(fun = function(y) Runuran::ud(law, y), colour = "red") +
#'   ggplot2::labs(x = "dX", y = "Density", title = "Increments of the Meixner process")
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
