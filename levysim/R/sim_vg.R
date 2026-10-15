#' Simulation of the variance-gamma process
#'
#' Simulates exactly one or several trajectories of a variance-gamma (VG)
#' process, homogeneous or not, at given inspection times, by subordination.
#'
#' @section Model:
#' Let \eqn{\Lambda} be continuous, strictly increasing, with
#' \eqn{\Lambda(0) = 0}, \eqn{\theta} real, \eqn{\sigma > 0} and \eqn{\nu > 0}.
#' Let \eqn{G} be a gamma process with
#' \eqn{G(t) - G(s) \sim Ga(\Delta\Lambda / \nu, \text{rate } 1 / \nu)}, so that
#' \eqn{E[G(t)] = \Lambda(t)} and \eqn{Var(G(t)) = \nu \Lambda(t)}, and \eqn{B}
#' an independent standard Brownian motion. The VG process is the Brownian
#' motion with drift evaluated at the random clock \eqn{G}:
#' \deqn{X(t) = \theta G(t) + \sigma B(G(t)).}
#' It has independent increments, \eqn{E[X(t)] = \theta \Lambda(t)} and
#' \eqn{Var(X(t)) = (\sigma^2 + \theta^2 \nu) \Lambda(t)}. The trajectories are
#' not monotone, purely made of jumps, of finite variation. With
#' \eqn{\Lambda(t) = t} (default), it is the homogeneous VG process, a Lévy
#' process. As
#' \eqn{\nu \to 0}, the VG process tends to the Wiener process
#' \eqn{\theta \Lambda(t) + \sigma B(\Lambda(t))}.
#'
#' @section Exact simulation 1: subordination (`method = "subordination"`):
#' Conditionally on the clock increment \eqn{\Delta G},
#' \eqn{X(t) - X(s) \sim N(\theta \Delta G, \sigma^2 \Delta G)}:
#' 1. Draw independently \eqn{\Delta G_i \sim Ga(\Delta\Lambda_i / \nu)} with
#'    scale \eqn{\nu}.
#' 2. \eqn{\Delta X_i = \theta \Delta G_i + \sigma \sqrt{\Delta G_i}\, Z_i},
#'    with \eqn{Z_i} independent \eqn{N(0, 1)}.
#' 3. Cumulative sum.
#'
#' This method also gives the random clock \eqn{G}, which is useful to
#' interpret the process.
#'
#' @section Exact simulation 2: difference of two gamma processes (`method = "difference"`):
#' Let
#' \deqn{\mu_\pm = \frac{1}{2} \sqrt{\theta^2 + \frac{2 \sigma^2}{\nu}} \pm \frac{\theta}{2} > 0,}
#' and \eqn{\Gamma_+}, \eqn{\Gamma_-} two independent gamma processes with
#' shape function \eqn{\Lambda(t) / \nu} and rates \eqn{1 / (\nu \mu_\pm)}.
#' Then \eqn{X} has the same distribution as \eqn{\Gamma_+ - \Gamma_-}: the
#' characteristic function of \eqn{\Delta\Gamma_+ - \Delta\Gamma_-} is
#' \deqn{(1 - i u \nu \mu_+)^{-\Delta\Lambda / \nu} (1 + i u \nu \mu_-)^{-\Delta\Lambda / \nu}
#'   = (1 - i u \nu (\mu_+ - \mu_-) + u^2 \nu^2 \mu_+ \mu_-)^{-\Delta\Lambda / \nu},}
#' with \eqn{\mu_+ - \mu_- = \theta} and \eqn{\nu^2 \mu_+ \mu_- = \sigma^2 \nu / 2},
#' which is the characteristic function
#' \eqn{(1 - i u \theta \nu + \sigma^2 \nu u^2 / 2)^{-\Delta\Lambda / \nu}} of
#' the VG increment.
#' 1. Draw independently \eqn{\Delta\Gamma^+_i \sim Ga(\Delta\Lambda_i / \nu)}
#'    with scale \eqn{\nu \mu_+} and
#'    \eqn{\Delta\Gamma^-_i \sim Ga(\Delta\Lambda_i / \nu)} with scale
#'    \eqn{\nu \mu_-}.
#' 2. \eqn{\Delta X_i = \Delta\Gamma^+_i - \Delta\Gamma^-_i}.
#' 3. Cumulative sum.
#'
#' This method only uses gamma draws. Both methods are exact and give the same
#' distribution.
#'
#' @inheritParams sim_gamma
#' @param theta Drift of the subordinated Brownian motion, real: mean
#'   degradation per unit of \eqn{\Lambda}.
#' @param sigma Diffusion coefficient of the subordinated Brownian motion,
#'   \eqn{\sigma > 0} (standard deviation, not variance).
#' @param nu Variance rate of the gamma clock, \eqn{\nu > 0}: the larger
#'   \eqn{\nu}, the larger the jumps; as \eqn{\nu \to 0} the process becomes a
#'   Wiener process.
#' @param Lambda Function of time: an R function of `t`, continuous and
#'   strictly increasing with \eqn{\Lambda(0) = 0}, for instance
#'   `function(t) t^1.5`. The default `function(t) t` gives the homogeneous
#'   process.
#' @param method Simulation algorithm: `"subordination"` (default) or
#'   `"difference"` (difference of two gamma processes). Both are exact and
#'   give the same distribution; see the sections above.
#'
#' @inherit sim_gamma return
#'
#' @references
#' Madan, D. B., Carr, P. P. and Chang, E. C. (1998). The variance gamma
#' process and option pricing. *European Finance Review*, 2(1), 79-105.
#' \doi{10.1023/A:1009703431535}
#'
#' @seealso [plot_traj()], [traj_df()], [inspection_times()]
#'
#' @examples
#' # One trajectory inspected every 0.5 time unit up to 10
#' t <- seq(0, 10, by = 0.5)
#' set.seed(1)
#' x <- sim_vg(t, theta = 0.5, sigma = 1.2, nu = 0.8)
#' x
#' plot_traj(x)
#'
#' # Ten trajectories, simulated by subordination (default)
#' x <- sim_vg(t, theta = 0.5, sigma = 1.2, nu = 0.8, n_traj = 10)
#' plot_traj(x)
#'
#' # The same model simulated as a difference of two gamma processes
#' y <- sim_vg(t, theta = 0.5, sigma = 1.2, nu = 0.8, method = "difference", n_traj = 10)
#' plot_traj(y, "Variance-gamma process, difference of two gamma processes")
#'
#' # Irregular inspections, more and more frequent over time
#' plan <- inspection_times("increasing", n = 20, n_traj = 5, tmax = 10)
#' x <- sim_vg(plan, theta = 0.5, sigma = 1.2, nu = 0.8)
#' plot_traj(x)
#'
#' # Different inspection times for each trajectory
#' x <- sim_vg(list(c(0, 1, 3, 6), c(0, 2, 4, 6, 8, 10)), theta = 0.5, sigma = 1.2, nu = 0.8)
#' x[[2]]
#'
#' # Non-homogeneous process with Lambda(t) = t^1.5
#' x <- sim_vg(t, theta = 0.5, sigma = 1.2, nu = 0.8, Lambda = function(t) t^1.5, n_traj = 10)
#' plot_traj(x, "Non-homogeneous variance-gamma process")
#'
#' # All the simulated data in a single data frame
#' head(traj_df(x))
#' @export
sim_vg <- function(t, theta, sigma, nu, Lambda = function(t) t, method = "subordination", n_traj = 1) {
  if (sigma <= 0 || nu <= 0) stop("sigma and nu must be positive")
  if (!method %in% c("subordination", "difference")) stop("method must be \"subordination\" or \"difference\"")
  r <- sqrt(theta^2 + 2 * sigma^2 / nu) / 2
  one_path <- function(t) {
    dL <- diff(Lambda(t))
    if (any(dL <= 0)) stop("Lambda must be strictly increasing")
    n <- length(dL)
    if (method == "subordination") {
      dG <- stats::rgamma(n, shape = dL / nu, scale = nu)
      dX <- theta * dG + sigma * sqrt(dG) * stats::rnorm(n)
    } else {
      dX <- stats::rgamma(n, shape = dL / nu, scale = nu * (r + theta / 2)) -
        stats::rgamma(n, shape = dL / nu, scale = nu * (r - theta / 2))
    }
    c(0, cumsum(dX))
  }
  simulate_paths(t, n_traj, one_path, "variance-gamma", c(theta = theta, sigma = sigma, nu = nu))
}
