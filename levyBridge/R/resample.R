as_trajectories <- function(x, t) {
  if (is.matrix(x)) x <- lapply(seq_len(ncol(x) - 1) + 1, function(k) data.frame(t = x[, 1], X = x[, k]))
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
  match.arg(method, c("Bridge", "ExtBridge"))
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

gcd <- function(a, b) {
  while (b > 0) {
    r <- a %% b
    a <- b
    b <- r
  }
  a
}

places <- function(d) paste(d, if (d == 1) "decimal place" else "decimal places")

on_grid <- function(t, d) all(abs(t * 10^d - round(t * 10^d)) <= 1e-6)

ext_digits <- function(t, digits) {
  if (is.null(digits)) {
    d <- which(vapply(0:6, function(k) on_grid(t, k), TRUE))
    if (length(d) > 0) return(d[1] - 1)
    warning("the inspection times are not multiples of 10^-6: they are rounded to 2 decimal places ",
            "(choose the number of decimal places with digits)", call. = FALSE)
    return(2)
  }
  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) || digits < 0 || digits != round(digits)) {
    stop("digits must be a non-negative integer")
  }
  if (!on_grid(t, digits)) warning("the inspection times are rounded to ", places(digits), call. = FALSE)
  digits
}

# Grilles ExtBridge de plusieurs trajectoires ; common = TRUE : même pas h pour toutes
ext_grids <- function(ts, digits = NULL, common = FALSE) {
  if (!common) return(lapply(ts, function(t) ext_grids(list(t), digits, TRUE)[[1]]))
  digits <- ext_digits(unlist(ts), digits)
  ti <- lapply(ts, function(t) round(t * 10^digits))
  if (any(vapply(ti, function(v) any(diff(v) <= 0), TRUE))) {
    stop("two inspection times are equal after rounding to ", places(digits), ": increase digits")
  }
  h <- Reduce(gcd, unlist(lapply(ti, diff)))
  lapply(ti, function(v) {
    n <- v[length(v)] / h
    if (n > 1e5) warning("the grid has ", n, " time steps: the resampling can be slow", call. = FALSE)
    list(t = v / 10^digits, grid = (0:n) * h / 10^digits, h = h / 10^digits)
  })
}

ext_grid <- function(t, digits = NULL) ext_grids(list(t), digits, TRUE)[[1]]

# Rééchantillonnage sur la grille g. draw(t_inf, delta, t, t_sup) tire le pont entre t_inf et t_sup
# (accroissement total delta) au temps t et renvoie c(gauche, droite) : les accroissements sur
# [t_inf, t] et [t, t_sup]. Les accroissements de la grille sont calculés directement, sans
# différence de niveaux : ils gardent toute leur précision relative, même très petits.
resample <- function(t, x, draw, g = regular_grid(t)) {
  n <- length(g) - 1
  N <- length(t)
  y <- numeric(n + 1)
  dx <- numeric(n)
  y[1] <- x[1]
  k <- 1
  t_inf <- t[1]
  x_inf <- x[1]
  rest <- x[2] - x[1]
  for (i in seq_len(n)) {
    carry <- 0
    # inspections strictement entre deux points de grille (méthode Bridge)
    while (t[k + 1] < g[i + 1]) {
      carry <- carry + rest
      k <- k + 1
      t_inf <- t[k]
      x_inf <- x[k]
      rest <- x[k + 1] - x[k]
    }
    if (t[k + 1] == g[i + 1]) {
      dx[i] <- carry + rest
      y[i + 1] <- x[k + 1]
      k <- k + 1
      if (k < N) {
        t_inf <- t[k]
        x_inf <- x[k]
        rest <- x[k + 1] - x[k]
      }
    } else {
      d <- draw(t_inf, rest, g[i + 1], t[k + 1])
      dx[i] <- carry + d[1]
      # niveau calculé depuis l'extrémité la plus proche : il reste entre les deux observations
      y[i + 1] <- if (abs(d[1]) <= abs(d[2])) x_inf + d[1] else x[k + 1] - d[2]
      t_inf <- g[i + 1]
      x_inf <- y[i + 1]
      rest <- d[2]
    }
  }
  list(path = data.frame(t = g, X = y), dx = dx)
}

# log d'une variable Gamma(s, 1), sans dépassement de capacité pour s petit : G_s = G_{s+1} U^{1/s}
log_rgamma <- function(s) log(stats::rgamma(1, s + 1)) + log(stats::runif(1)) / s

# Y ~ Beta(p, q) et 1 - Y, chacun avec toute sa précision relative
rbeta_pair <- function(p, q) {
  l <- log_rgamma(p) - log_rgamma(q)
  c(1 / (1 + exp(-l)), 1 / (1 + exp(l)))
}

bridge_paths <- function(x, t, method, draw, model, parameters, increasing = FALSE, digits = NULL,
                         common_step = FALSE) {
  method <- check_method(method)
  if (!is.logical(common_step) || length(common_step) != 1 || is.na(common_step)) stop("common_step must be TRUE or FALSE")
  traj <- as_trajectories(x, t)
  if (increasing) {
    for (k in seq_along(traj)) {
      if (any(diff(traj[[k]]$x) < 0)) stop("trajectory ", k, ": the ", model, " process is increasing, the levels must not decrease")
    }
  }
  if (method == "Bridge") {
    res <- lapply(traj, function(d) resample(d$t, d$x, draw))
    teq <- vapply(traj, function(d) d$t[length(d$t)] / (length(d$t) - 1), 0)
  } else {
    grids <- ext_grids(lapply(traj, `[[`, "t"), digits, common_step)
    for (k in seq_along(traj)) traj[[k]]$t <- grids[[k]]$t
    res <- lapply(seq_along(traj), function(k) resample(traj[[k]]$t, traj[[k]]$x, draw, grids[[k]]$grid))
    teq <- vapply(grids, function(g) g$h, 0)
  }
  observed <- lapply(traj, function(d) data.frame(t = d$t, X = d$x))
  structure(lapply(res, `[[`, "path"), class = "levybridge", model = model, parameters = parameters, method = method,
            teq = teq, observed = observed, dx = lapply(res, `[[`, "dx"))
}

#' Reconstructed increments
#'
#' Extracts the increments of the trajectories resampled on the regular grid
#' by a `_bridge` function. Under the model, with the true parameters, these
#' increments are independent and have the law of an increment of the process
#' over a time step `teq` (Grall-Maës, 2012): they are the input of the
#' goodness-of-fit tests.
#'
#' @param x Result of a bridge function of the package, for instance
#'   [gamma_bridge()].
#'
#' @return A data frame with one row per increment and three columns:
#'   * `id`: number of the trajectory (1, 2, ...);
#'   * `teq`: time step of the regular grid of this trajectory, \eqn{t_N / N}
#'     for the method "Bridge" and \eqn{h} for the method "ExtBridge";
#'   * `dx`: reconstructed increment \eqn{x'_i - x'_{i-1}}.
#'
#' The increments are computed directly by the bridges, not as differences of
#' the resampled levels: a very small increment keeps its relative precision
#' instead of being rounded to 0 (the levels have about 16 significant
#' digits, so an increment smaller than about \eqn{10^{-16}} times the level
#' would be lost by a difference).
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
  dx <- attr(x, "dx")
  do.call(rbind, lapply(seq_along(x), function(i)
    data.frame(id = i, teq = teq[i], dx = if (is.null(dx)) diff(x[[i]]$X) else dx[[i]])))
}

#' @export
print.levybridge <- function(x, ...) {
  par <- attr(x, "parameters")
  cat(length(x), "trajectory(ies) resampled on a regular grid by the", attr(x, "model"), "bridge",
      paste0("(method \"", attr(x, "method"), "\")\n"))
  cat("parameters:", paste(names(par), signif(par, 4), sep = " = ", collapse = ", "), "\n")
  teq <- signif(unique(attr(x, "teq")), 4)
  if (length(teq) <= 4) {
    cat("time step teq:", paste(teq, collapse = ", "), "\n\n")
  } else {
    cat("time step teq: from", min(teq), "to", max(teq), "(one per trajectory, see attr(x, \"teq\"))\n\n")
  }
  d <- do.call(rbind, lapply(seq_along(x), function(i) data.frame(id = i, x[[i]])))
  print(utils::head(d))
  invisible(x)
}
