as_trajectories <- function(x, t) {
  if (is.matrix(x)) x <- lapply(seq_len(ncol(x) - 1) + 1, function(k) data.frame(t = x[, 1], X = x[, k]))
  if (is.data.frame(x)) x <- list(x)
  if (is.list(x)) {
    traj <- lapply(x, function(d) list(t = d[[1]], x = d[[2]]))
  } else {
    traj <- list(list(t = t, x = x))
  }
  for (k in seq_along(traj)) {
    if (length(traj[[k]]$t) != length(traj[[k]]$x)) stop("trajectory ", k, ": t and x must have the same length")
    if (length(traj[[k]]$t) < 2) stop("trajectory ", k, ": at least two observations are needed")
    if (any(diff(traj[[k]]$t) <= 0)) stop("inspection times must be strictly increasing")
  }
  traj
}

increments <- function(x, t) {
  if (inherits(x, "levy_increments")) return(list(dt = x$dt, dx = x$dx))
  traj <- as_trajectories(x, t)
  list(dt = unlist(lapply(traj, function(d) diff(d$t))),
       dx = unlist(lapply(traj, function(d) diff(d$x))))
}

#' Data given by their increments
#'
#' Builds the data of the estimation functions ([gamma_mle()],
#' [wiener_mle()], [ig_mle()], [vg_mle()], [meixner_mle()]) from increments
#' instead of degradation levels: the increments \eqn{\Delta x_j} over the
#' durations \eqn{\Delta t_j}, \eqn{j = 1, \dots, n}, of one or several
#' trajectories. The likelihood of these processes depends on the data only
#' through the increments. This is the fourth form of the data accepted by the
#' estimation functions, after the three forms based on the degradation levels
#' (see the section Data of [levyfit-package]).
#'
#' Giving the increments directly keeps all their precision. In R,
#' `.Machine$double.eps` (about \eqn{2.2 \times 10^{-16}}) is the smallest
#' positive number \eqn{x} such that \eqn{1 + x \ne 1} (see [base::.Machine]), so
#' that an increment smaller than about \eqn{10^{-16}} times the level can be
#' rounded to 0 when computed as a difference of levels; this happens for
#' instance with the increments reconstructed by the bridges of the package
#' `levyBridge` on a fine grid.
#'
#' @param dx Numeric vector of increments.
#' @param dt Durations of the increments: one positive number (all the
#'   increments over the same duration, 1 by default) or a vector of the same
#'   length as `dx`.
#'
#' @return An object of class `"levy_increments"`, to be given as argument `x`
#'   of an estimation function (without `t`).
#'
#' @examples
#' set.seed(1)
#' # Gamma increments over a time step 0.1, some of them very small
#' dx <- rgamma(200, shape = 2 * 0.1, rate = 3)
#' min(dx)
#' gamma_mle(increments_data(dx, dt = 0.1))
#'
#' # Same estimate as from the levels, when no increment is lost
#' gamma_mle(c(0, cumsum(dx)), seq(0, 20, by = 0.1))
#'
#' # Increments over different durations
#' dt <- runif(50, 0.5, 2)
#' wiener_mle(increments_data(rnorm(50, 0.5 * dt, 1.5 * sqrt(dt)), dt))
#' @export
increments_data <- function(dx, dt = 1) {
  if (!is.numeric(dx) || length(dx) == 0 || !all(is.finite(dx))) stop("dx must be a vector of finite numbers")
  if (!is.numeric(dt) || !(length(dt) %in% c(1, length(dx))) || !all(is.finite(dt)) || any(dt <= 0)) {
    stop("dt must be one positive number or a vector of positive numbers of the same length as dx")
  }
  structure(list(dt = rep_len(dt, length(dx)), dx = dx), class = "levy_increments")
}
