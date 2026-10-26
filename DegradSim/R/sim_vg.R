#' Simulation of the variance-gamma process
#'
#' Simulates exactly, at given inspection times, one or several trajectories of
#' a homogeneous variance-gamma (VG) process, a Lévy process.
#'
#' @section Model:
#' Let \eqn{\theta} be real, \eqn{\sigma > 0} and \eqn{\nu > 0}. Let \eqn{G}
#' be a gamma process with
#' \eqn{G(t) - G(s) \sim \Gamma((t - s) / \nu, \text{rate } 1 / \nu)}, so that
#' \eqn{E[G(t)] = t} and \eqn{Var(G(t)) = \nu t}, and \eqn{B} an independent
#' standard Brownian motion. The VG process is the Brownian motion with drift
#' evaluated at the random clock \eqn{G}:
#' \deqn{X(t) = \theta G(t) + \sigma B(G(t)), \qquad t \ge 0.}
#' It is a Lévy process: \eqn{X(0) = 0}, its increments are independent and
#' stationary, \eqn{E[X(t)] = \theta t} and
#' \eqn{Var(X(t)) = (\sigma^2 + \theta^2 \nu)\, t}. The trajectories are not
#' monotone, purely made of jumps, of finite variation. As \eqn{\nu \to 0}, the
#' VG process tends to the Wiener process \eqn{\theta t + \sigma B(t)}.
#'
#' @section Exact simulation:
#' Random walk approximation, described in [DegradSim-package]: each
#' increment, over a time \eqn{\Delta t_i = t_i - t_{i-1}}, is drawn by one of
#' two exact methods (argument `method`), which give the same distribution,
#' then the increments are cumulated.
#'
#' **Method 1, subordination** (`method = "subordination"`, default).
#' Conditionally on the clock increment \eqn{\Delta G},
#' \eqn{X(t) - X(s) \sim N(\theta \Delta G, \sigma^2 \Delta G)}:
#' 1. draw independently \eqn{\Delta G_i \sim \Gamma(\Delta t_i / \nu)} with
#'    scale \eqn{\nu};
#' 2. \eqn{\Delta X_i = \theta \Delta G_i + \sigma \sqrt{\Delta G_i}\, Z_i},
#'    with \eqn{Z_i} independent \eqn{N(0, 1)}.
#'
#' This method follows the definition of the process and also gives the
#' random clock \eqn{G}.
#'
#' **Method 2, difference of two gamma processes** (`method = "difference"`).
#' Let
#' \deqn{\mu_\pm = \frac{1}{2} \sqrt{\theta^2 + \frac{2 \sigma^2}{\nu}} \pm \frac{\theta}{2} > 0,}
#' and \eqn{G_+}, \eqn{G_-} two independent gamma processes with shapes
#' \eqn{t / \nu} and scales \eqn{\nu \mu_\pm}. Then \eqn{X} has the same
#' distribution as \eqn{G_+ - G_-}: since \eqn{\mu_+ - \mu_- = \theta} and
#' \eqn{\nu^2 \mu_+ \mu_- = \sigma^2 \nu / 2}, the characteristic function of
#' \eqn{\Delta G_+ - \Delta G_-} over a time \eqn{\Delta t} is
#' \deqn{(1 - i u \nu \mu_+)^{-\Delta t / \nu} (1 + i u \nu \mu_-)^{-\Delta t / \nu}
#'   = \big(1 - i u \theta \nu + \sigma^2 \nu u^2 / 2\big)^{-\Delta t / \nu},}
#' which is the characteristic function of the VG increment (Madan et al.,
#' 1998).
#' 1. Draw independently \eqn{\Delta G^+_i \sim \Gamma(\Delta t_i / \nu)} with
#'    scale \eqn{\nu \mu_+} and \eqn{\Delta G^-_i \sim \Gamma(\Delta t_i / \nu)}
#'    with scale \eqn{\nu \mu_-};
#' 2. \eqn{\Delta X_i = \Delta G^+_i - \Delta G^-_i}.
#'
#' This method only uses gamma draws.
#'
#' @inheritParams sim_gamma
#' @param theta Drift of the subordinated Brownian motion, real: mean
#'   degradation per unit of time.
#' @param sigma Diffusion coefficient of the subordinated Brownian motion,
#'   \eqn{\sigma > 0} (standard deviation, not variance).
#' @param nu Variance rate of the gamma clock, \eqn{\nu > 0}: the larger
#'   \eqn{\nu}, the larger the jumps; as \eqn{\nu \to 0} the process becomes a
#'   Wiener process.
#' @param method Simulation algorithm: `"subordination"` (default) or
#'   `"difference"` (difference of two gamma processes). Both are exact and
#'   give the same distribution; see the section on the exact simulation.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Madan, D. B., Carr, P. P. and Chang, E. C. (1998). The variance gamma
#' process and option pricing. *European Finance Review*, 2(1), 79-105.
#' \doi{10.1023/A:1009703431535}
#'
#' @seealso [plot_traj()], [traj_df()], [inspection_times()], [plot_plan()]
#'
#' @examples
#' # 1. One trajectory inspected every 0.5 time unit up to 10
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_vg(t, theta = 0.5, sigma = 1.2, nu = 0.8)
#' x
#' plot_traj(x)
#'
#' # 2. Ten trajectories, simulated by subordination (default)
#' x <- sim_vg(t, theta = 0.5, sigma = 1.2, nu = 0.8, n_traj = 10)
#' plot_traj(x)
#'
#' # Same times for all the trajectories: a matrix with the column t,
#' # then one column per trajectory (X1, X2, ...)
#' head(x[, 1:4])
#' x[, "X2"]
#'
#' # 3. The same model simulated as a difference of two gamma processes
#' y <- sim_vg(t, theta = 0.5, sigma = 1.2, nu = 0.8, method = "difference", n_traj = 10)
#' plot_traj(y, "Variance-gamma process, difference of two gamma processes")
#'
#' # 4. Irregular inspections, more and more frequent over time
#' plan <- inspection_times("increasing", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_vg(plan, theta = 0.5, sigma = 1.2, nu = 0.8)
#' plot_traj(x)
#'
#' # 5. Different inspection times for each trajectory: a list of data frames
#' x <- sim_vg(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), theta = 0.5, sigma = 1.2, nu = 0.8)
#' x[[2]]
#' x[[2]]$X
#'
#' # 6. Small nu: the trajectories look like those of a Wiener process
#' x <- sim_vg(t, theta = 0.5, sigma = 1.2, nu = 0.01, n_traj = 10)
#' plot_traj(x, "Variance-gamma process with nu = 0.01")
#'
#' # 7. All the simulated data in a single data frame, with the increments
#' head(traj_df(x))
#' head(traj_df(x, increments = TRUE))
#' @export
sim_vg <- function(t, theta, sigma, nu, method = "subordination", n_traj = 1) {
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta)) stop("theta must be a real number")
  if (!is.numeric(sigma) || length(sigma) != 1 || !(sigma > 0)) stop("sigma must be a positive number")
  if (!is.numeric(nu) || length(nu) != 1 || !(nu > 0)) stop("nu must be a positive number")
  if (!identical(method, "subordination") && !identical(method, "difference")) {
    stop("method must be \"subordination\" or \"difference\"")
  }
  r <- sqrt(theta^2 + 2 * sigma^2 / nu) / 2
  one_path <- function(t) {
    dt <- diff(t)
    n <- length(dt)
    if (method == "subordination") {
      dG <- stats::rgamma(n, shape = dt / nu, scale = nu)
      dX <- theta * dG + sigma * sqrt(dG) * stats::rnorm(n)
    } else {
      dX <- stats::rgamma(n, shape = dt / nu, scale = nu * (r + theta / 2)) -
        stats::rgamma(n, shape = dt / nu, scale = nu * (r - theta / 2))
    }
    c(0, cumsum(dX))
  }
  simulate_paths(t, n_traj, one_path, "variance-gamma", c(theta = theta, sigma = sigma, nu = nu))
}
