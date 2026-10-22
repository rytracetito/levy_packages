#' Inverse Gaussian bridge: resampling on a regular grid
#'
#' Resamples one or several trajectories of a homogeneous inverse Gaussian
#' process, observed at irregular inspection times, onto a regular time grid,
#' by exact simulation of the inverse Gaussian bridge between the known
#' values.
#'
#' @section Model:
#' The homogeneous inverse Gaussian process with parameters \eqn{\mu > 0} and
#' \eqn{\eta > 0} has \eqn{X_0 = 0}, independent increments and
#' \eqn{X_t - X_s \sim IG(\mu (t - s), \eta (t - s)^2)} (mean
#' \eqn{\mu (t - s)}, shape \eqn{\eta (t - s)^2}). Its trajectories are
#' increasing.
#'
#' @section Law of the bridge:
#' Let \eqn{t_A < t < t_B}, \eqn{\tau_1 = t - t_A}, \eqn{\tau_2 = t_B - t},
#' \eqn{T = \tau_1 + \tau_2}, \eqn{\Delta = x_B - x_A > 0},
#' \eqn{m = \tau_1 / \tau_2} and \eqn{\lambda = \eta \tau_1^2 / \Delta}. Given
#' \eqn{X_{t_A} = x_A} and \eqn{X_{t_B} = x_B}, the ratio
#' \eqn{W = (X_t - x_A) / (x_B - X_t)} ("realized / remaining") is the
#' mixture
#' \deqn{W \sim IG(m, \lambda) \text{ with probability } \tau_2 / T, \qquad
#'   W \sim m^2 / IG(m, \lambda) \text{ with probability } \tau_1 / T,}
#' and \eqn{X_t = x_A + \Delta\, W / (1 + W)}. The law of the bridge does not
#' depend on \eqn{\mu}: only \eqn{\eta} is needed.
#'
#' The simulation combines this mixture with the method of Michael, Schucany
#' and Haas (1976), with one Gaussian and one uniform variable: with
#' \eqn{q = Z^2}, \eqn{Z \sim N(0, 1)}, the two roots
#' \deqn{w_2 = m + \frac{m^2 q}{2 \lambda} + \frac{m}{2 \lambda} \sqrt{4 m \lambda q + m^2 q^2}, \qquad w_1 = m^2 / w_2,}
#' and \eqn{W = w_1} with probability
#' \eqn{\tau_1 (1 + w_1) / (T (m + w_1))}, \eqn{W = w_2} otherwise. The small
#' root is computed as \eqn{m^2 / w_2} to avoid a cancellation. When
#' \eqn{\Delta = 0} (two equal observed values, which can happen after
#' rounding), the bridge is constant.
#'
#' @inheritSection gamma_bridge Resampling algorithm (method "Bridge")
#' @inheritSection gamma_bridge Method "ExtBridge"
#' @inheritSection gamma_bridge Method "MultiBridge"
#'
#' @inheritParams gamma_bridge
#' @param mu Mean degradation per unit of time, \eqn{\mu > 0}; used only with
#'   `method = "MultiBridge"` and `option = 2`, to extend the trajectories.
#' @param eta Shape parameter \eqn{\eta > 0}, for instance the estimate
#'   `ig_mle(x)[["eta"]]` of the package `levyfit`.
#'
#' @inherit gamma_bridge return
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
#' variates using transformations with multiple roots. *The American
#' Statistician*, 30(2), 88-90.\cr
#' \doi{10.1080/00031305.1976.10479147}
#'
#' Ribeiro, C. and Webber, N. (2003). A Monte Carlo method for the normal
#' inverse Gaussian option valuation model using an inverse Gaussian bridge.
#' *Computing in Economics and Finance 2003*, 5.
#'
#' @seealso
#' * [increments()]: the reconstructed increments;
#' * [plot_bridge()]: observed and resampled values;
#' * [gamma_bridge()]: bridge of the gamma process;
#' * [wiener_bridge()]: bridge of the Wiener process;
#' * [vg_bridge()]: bridge of the variance-gamma process;
#' * [meixner_bridge()]: bridge of the Meixner process.
#'
#' @examples
#' # 1. One trajectory observed at 11 irregular times (mu = 1, eta = 3)
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' x <- c(0, 1.31, 2.05, 3.62, 5.71, 7.15, 10.62, 12.83, 15.02, 18.94, 22.31)
#' set.seed(1)
#' r <- ig_bridge(x, t, eta = 3)
#' r
#'
#' # The resampled trajectory, on the regular grid 0, 2.25, 4.5, ..., 22.5
#' r[[1]]
#'
#' # Time step of the grid: t_N / N = 22.5 / 10
#' attr(r, "teq")
#'
#' # 2. Reconstructed increments, the input of the goodness-of-fit tests
#' increments(r)
#'
#' # 3. Observed and resampled values
#' plot_bridge(r)
#'
#' # 4. Grid points equal to inspection times keep the observed values:
#' #    here the grid is 0, 2, 4, 6, 8, 10 and 2 and 6 are inspection times
#' t2 <- c(0, 2, 3, 6, 7, 10)
#' x2 <- c(0, 1.4, 2.0, 4.1, 4.6, 6.8)
#' ig_bridge(x2, t2, eta = 3)[[1]]
#'
#' # 5. One trajectory given as a data frame (times, then levels)
#' ig_bridge(data.frame(t = t2, X = x2), eta = 3)
#'
#' # 6. Several trajectories, each with its own inspection times:
#' #    one regular grid and one time step per trajectory
#' liste <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.8, 2.1, 4.3)),
#'               data.frame(t = c(0, 2, 4, 5, 9), X = c(0, 1.2, 2.9, 3.1, 6.2)))
#' r2 <- ig_bridge(liste, eta = 3)
#' r2
#' attr(r2, "teq")
#' increments(r2)
#' plot_bridge(r2, trajectory = 1:2)
#'
#' # 7. The bridge values are random: set.seed() makes them reproducible
#' set.seed(2); y1 <- ig_bridge(x, t, eta = 3)[[1]]$X
#' set.seed(2); y2 <- ig_bridge(x, t, eta = 3)[[1]]$X
#' identical(y1, y2)
#'
#' # 8. Equal observed values (after rounding) give a constant bridge,
#' #    a decrease is impossible for this increasing process
#' ig_bridge(c(0, 1, 1, 2), c(0, 1, 2.5, 3), eta = 3)[[1]]
#' try(ig_bridge(c(0, 2, 1.5, 3), c(0, 1, 2, 3), eta = 3))
#'
#' # 9. Method "ExtBridge": the times have 1 decimal place, the grid
#' #    0, 0.1, ..., 22.5 contains all of them
#' r3 <- ig_bridge(x, t, eta = 3, method = "ExtBridge")
#' attr(r3, "teq")
#' all(t %in% r3[[1]]$t)
#' plot_bridge(r3)
#'
#' # Error: the first time must be 0
#' try(ig_bridge(x, t + 1, eta = 3))
#' @examplesIf requireNamespace("levysim", quietly = TRUE) && requireNamespace("levyfit", quietly = TRUE)
#' # Method "MultiBridge": several trajectories with their own times on one common grid,
#' # so that all the reconstructed increments have the same law
#' liste4 <- list(data.frame(t = c(0, 3, 8), X = c(0, 1.5, 4.2)),
#'                data.frame(t = c(0, 2, 5, 9, 12), X = c(0, 0.9, 2.8, 4.1, 6.5)),
#'                data.frame(t = c(0, 4, 7, 11, 15), X = c(0, 2.2, 3.1, 5.6, 7.4)))
#' r4 <- ig_bridge(liste4, eta = 3, method = "MultiBridge")
#' attr(r4, "teq")
#' sapply(r4, nrow)
#'
#' # 10. Ten trajectories simulated with levysim, inspected more and more often
#' plan <- levysim::inspection_times("increasing", n = 20, n_traj = 10, tmax = 40)
#' sims <- levysim::sim_ig(plan, mu = 1, eta = 3)
#'
#' # Parameters estimated with levyfit
#' eta_hat <- levyfit::ig_mle(sims)[["eta"]]
#' r <- ig_bridge(sims, eta = eta_hat)
#' r
#' plot_bridge(r, trajectory = 1:4)
#'
#' # Reconstructed increments of all the trajectories
#' d <- increments(r)
#' head(d)
#' # Under the model, the increments follow IG(mu * teq, eta * teq^2)
#' # (mean mu * teq, variance mu^3 * teq / eta)
#' c(mean(d$dx), 1 * d$teq[1], var(d$dx), d$teq[1] / 3)
#' @export
ig_bridge <- function(x, t = NULL, eta, method = "Bridge", digits = NULL, common_step = FALSE, option = 1,
                      mu = NULL) {
  check_par(eta = eta)
  sim <- NULL
  if (!is.null(mu)) {
    check_par(mu = mu)
    # IG(mu s, eta s^2) par la méthode de Michael, Schucany et Haas ; petite racine = m^2 / grande racine
    sim <- function(s) {
      m <- mu * s
      lam <- eta * s^2
      q <- stats::rnorm(1)^2
      w2 <- m + m^2 * q / (2 * lam) + m / (2 * lam) * sqrt(4 * m * lam * q + m^2 * q^2)
      w1 <- m^2 / w2
      if (stats::runif(1) <= m / (m + w1)) w1 else w2
    }
  }
  draw <- function(t_inf, delta, t, t_sup) {
    if (delta == 0) return(c(0, 0))
    tau1 <- t - t_inf
    tau2 <- t_sup - t
    m <- tau1 / tau2
    lam <- eta * tau1^2 / delta
    q <- stats::rnorm(1)^2
    w2 <- m + m^2 * q / (2 * lam) + m / (2 * lam) * sqrt(4 * m * lam * q + m^2 * q^2)
    w1 <- m^2 / w2
    w <- if (stats::runif(1) <= tau1 * (1 + w1) / ((tau1 + tau2) * (m + w1))) w1 else w2
    delta * c(w / (1 + w), 1 / (1 + w))
  }
  bridge_paths(x, t, method, draw, "inverse Gaussian", c(eta = unname(eta)), increasing = TRUE, digits = digits,
               common_step = common_step, option = option, sim = sim)
}
