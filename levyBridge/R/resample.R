as_trajectories <- function(x, t) {
  if (is.data.frame(x)) x <- list(x)
  if (is.list(x)) {
    traj <- lapply(x, function(d) list(t = d[[1]], x = d[[2]]))
  } else {
    if (is.null(t)) stop("give the inspection times t, or a list of data frames")
    traj <- list(list(t = t, x = x))
  }
  for (k in seq_along(traj)) {
    d <- traj[[k]]
    if (!is.numeric(d$t) || !is.numeric(d$x)) stop("trajectory ", k, ": times and levels must be numeric")
    if (length(d$t) != length(d$x)) stop("trajectory ", k, ": t and x must have the same length")
    if (length(d$t) < 2) stop("trajectory ", k, ": at least two observations are needed")
    if (anyNA(d$t) || anyNA(d$x) || any(!is.finite(d$t)) || any(!is.finite(d$x))) stop("trajectory ", k, ": missing or infinite values")
    if (d$t[1] != 0) stop("trajectory ", k, ": the first inspection time must be 0")
    if (any(diff(d$t) <= 0)) stop("trajectory ", k, ": inspection times must be strictly increasing")
  }
  traj
}

check_par <- function(..., positive = TRUE) {
  p <- list(...)
  for (nm in names(p)) {
    v <- p[[nm]]
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v)) stop(nm, " must be a single finite number")
    if (positive && v <= 0) stop(nm, " must be positive")
  }
}

check_method <- function(method) {
  method <- match.arg(method, c("Bridge", "ExtBridge"))
  if (method == "ExtBridge") stop("method \"ExtBridge\" is not available yet: use method = \"Bridge\"")
  method
}

regular_grid <- function(t) {
  n <- length(t) - 1
  g <- seq(0, t[n + 1], length.out = n + 1)
  tol <- 1e-9 * t[n + 1]
  for (i in seq_along(g)) {
    j <- which(abs(t - g[i]) <= tol)
    if (length(j) > 0) g[i] <- t[j[1]]
  }
  g
}

resample <- function(t, x, draw) {
  g <- regular_grid(t)
  n <- length(g) - 1
  y <- numeric(n + 1)
  y[1] <- x[1]
  y[n + 1] <- x[n + 1]
  for (i in seq_len(n - 1) + 1) {
    k <- findInterval(g[i], t)
    if (t[k] == g[i]) {
      y[i] <- x[k]
    } else {
      if (t[k] > g[i - 1]) {
        t_inf <- t[k]
        x_inf <- x[k]
      } else {
        t_inf <- g[i - 1]
        x_inf <- y[i - 1]
      }
      y[i] <- draw(t_inf, x_inf, g[i], t[k + 1], x[k + 1])
    }
  }
  data.frame(t = g, X = y)
}

bridge_paths <- function(x, t, method, draw, model, parameters, increasing = FALSE) {
  method <- check_method(method)
  traj <- as_trajectories(x, t)
  if (increasing) {
    for (k in seq_along(traj)) {
      if (any(diff(traj[[k]]$x) < 0)) stop("trajectory ", k, ": the ", model, " process is increasing, the levels must not decrease")
    }
  }
  paths <- lapply(traj, function(d) resample(d$t, d$x, draw))
  teq <- vapply(traj, function(d) d$t[length(d$t)] / (length(d$t) - 1), 0)
  observed <- lapply(traj, function(d) data.frame(t = d$t, X = d$x))
  structure(paths, class = "levybridge", model = model, parameters = parameters, method = method,
            teq = teq, observed = observed)
}

#' Reconstructed increments
#'
#' Extracts the increments of the trajectories resampled on the regular grid
#' by a `_bridge` function. Under the model, with the true parameters, these
#' increments are independent and have the law of an increment of the process
#' over a time step `teq` (Grall-Maës, 2012): they are the input of the
#' goodness-of-fit tests.
#'
#' @param x Result of a bridge function: [gamma_bridge()], [wiener_bridge()],
#'   [ig_bridge()], [vg_bridge()] or [meixner_bridge()].
#'
#' @return A data frame with one row per increment and three columns:
#'   * `id`: number of the trajectory (1, 2, ...);
#'   * `teq`: time step of the regular grid of this trajectory, \eqn{t_N / N};
#'   * `dx`: reconstructed increment \eqn{x'_i - x'_{i-1}}.
#'
#' @seealso [gamma_bridge()], [wiener_bridge()], [ig_bridge()], [vg_bridge()],
#'   [meixner_bridge()]
#'
#' @examples
#' # One gamma trajectory observed at irregular times
#' t <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)
#' set.seed(1)
#' x <- c(0, cumsum(rgamma(10, shape = 2 * diff(t), rate = 3)))
#' r <- gamma_bridge(x, t, a = 2)
#' increments(r)
#'
#' # Two trajectories with different times: one time step per trajectory
#' liste <- list(data.frame(t = c(0, 1, 3, 6), X = c(0, 0.8, 2.1, 4.3)),
#'               data.frame(t = c(0, 2, 4, 5, 9), X = c(0, 1.2, 2.9, 3.1, 6.2)))
#' increments(gamma_bridge(liste, a = 2))
#' @export
increments <- function(x) {
  if (!inherits(x, "levybridge")) stop("x must be the result of a _bridge function")
  teq <- attr(x, "teq")
  do.call(rbind, lapply(seq_along(x), function(i)
    data.frame(id = i, teq = teq[i], dx = diff(x[[i]]$X))))
}

#' @export
print.levybridge <- function(x, ...) {
  par <- attr(x, "parameters")
  cat(length(x), "trajectory(ies) resampled on a regular grid by the", attr(x, "model"), "bridge",
      paste0("(method \"", attr(x, "method"), "\")\n"))
  cat("parameters:", paste(names(par), signif(par, 4), sep = " = ", collapse = ", "), "\n")
  cat("time step teq:", paste(signif(unique(attr(x, "teq")), 4), collapse = ", "), "\n\n")
  d <- do.call(rbind, lapply(seq_along(x), function(i) data.frame(id = i, x[[i]])))
  print(utils::head(d))
  invisible(x)
}
