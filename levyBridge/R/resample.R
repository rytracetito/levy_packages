as_trajectories <- function(x) {
  if (is.matrix(x)) {
    if (ncol(x) < 2) stop("the matrix x must have the times in column 1 and the levels in the next columns")
    x <- lapply(seq_len(ncol(x) - 1) + 1, function(k) data.frame(t = x[, 1], X = x[, k]))
  }
  if (is.data.frame(x)) x <- list(x)
  if (!is.list(x) || length(x) == 0 || !all(vapply(x, is.data.frame, TRUE))) {
    stop("x must be a data frame (times, levels), a list of such data frames, or a matrix (times, levels...)")
  }
  traj <- lapply(seq_along(x), function(k) {
    if (ncol(x[[k]]) < 2) stop("trajectory ", k, ": the data frame must have the times in column 1 and the levels in column 2")
    list(t = x[[k]][[1]], x = x[[k]][[2]])
  })
  for (k in seq_along(traj)) {
    d <- traj[[k]]
    if (!is.numeric(d$t) || !is.numeric(d$x)) stop("trajectory ", k, ": times and levels must be numeric")
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

# theta : vecteur nommé de tous les paramètres du processus ; real = paramètres réels (les autres sont > 0)
check_theta <- function(theta, required, real = character(0)) {
  if (!is.numeric(theta) || is.null(names(theta))) {
    stop("theta must be a named numeric vector with the parameters ", paste(required, collapse = ", "), call. = FALSE)
  }
  miss <- setdiff(required, names(theta))
  if (length(miss) > 0) stop("theta must contain the parameter(s) ", paste(miss, collapse = ", "), call. = FALSE)
  p <- vapply(required, function(nm) unname(theta[[nm]]), 0)
  for (nm in required) {
    if (!is.finite(p[[nm]])) stop(nm, " must be a finite number", call. = FALSE)
    if (!(nm %in% real) && p[[nm]] <= 0) stop(nm, " must be positive", call. = FALSE)
  }
  p
}

# n tirages de X_t sachant X_{t_A} = x_A et X_{t_B} = x_B ; le niveau est calculé depuis l'extrémité la plus proche
bridge_draws <- function(t_A, x_A, t_B, x_B, t, n, draw, model, increasing = FALSE) {
  check_par(t_A = t_A, x_A = x_A, t_B = t_B, x_B = x_B, t = t, positive = FALSE)
  if (!(t_A < t && t < t_B)) stop("the times must satisfy t_A < t < t_B", call. = FALSE)
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n < 1 || n != round(n)) {
    stop("n must be a positive integer", call. = FALSE)
  }
  if (increasing && x_B < x_A) stop("the ", model, " process is increasing: x_B must not be smaller than x_A", call. = FALSE)
  vapply(seq_len(n), function(i) {
    d <- draw(t_A, x_B - x_A, t, t_B)
    if (abs(d[1]) <= abs(d[2])) x_A + d[1] else x_B - d[2]
  }, 0)
}

check_method <- function(method) {
  match.arg(method, c("Bridge", "ExtBridge", "MultiBridge"))
}

# Point de grille égal à une inspection à l'arrondi près : on prend l'instant observé
snap_grid <- function(g, t) {
  tol <- 1e-9 * max(t, g)
  for (i in seq_along(g)) {
    j <- which(abs(t - g[i]) <= tol)
    if (length(j) > 0) g[i] <- t[j[1]]
  }
  g
}

regular_grid <- function(t) {
  n <- length(t) - 1
  snap_grid(seq(0, t[n + 1], length.out = n + 1), t)
}

# Grilles MultiBridge : pas commun h = (somme des durées) / (somme des nombres d'accroissements) ;
# option 1 : chaque trajectoire jusqu'à son dernier instant ; option 2 : toutes jusqu'au dernier instant global
multi_grids <- function(ts, option) {
  t_max <- vapply(ts, function(t) t[length(t)], 0)
  h <- sum(t_max) / sum(lengths(ts) - 1)
  n_all <- floor(max(t_max) / h + 1e-9)
  lapply(seq_along(ts), function(k) {
    n <- if (option == 1) floor(t_max[k] / h + 1e-9) else n_all
    if (n == 0) warning("trajectory ", k, " is shorter than the time step ", signif(h, 4),
                        ": it gives no reconstructed increment", call. = FALSE)
    list(grid = snap_grid(seq(0, n * h, length.out = n + 1), ts[[k]]), h = h)
  })
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
    stop("the inspection times have more than 6 decimal places: give the number of decimal places with digits",
         call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) || digits < 0 || digits != round(digits)) {
    stop("digits must be a non-negative integer")
  }
  if (!on_grid(t, digits)) warning("the inspection times are truncated to ", places(digits), call. = FALSE)
  digits
}

# Grilles ExtBridge de plusieurs trajectoires ; common = TRUE : même pas h pour toutes
ext_grids <- function(ts, digits = NULL, common = FALSE) {
  if (!common) return(lapply(ts, function(t) ext_grids(list(t), digits, TRUE)[[1]]))
  digits <- ext_digits(unlist(ts), digits)
  # partie entière inférieure, à une erreur de calcul près (0.3 * 10 = 2.9999999999999996)
  ti <- lapply(ts, function(t) floor(t * 10^digits + 1e-6))
  if (any(vapply(ti, function(v) any(diff(v) <= 0), TRUE))) {
    stop("two inspection times are equal after truncation to ", places(digits), ": increase digits")
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
# Au-delà de la dernière observation (MultiBridge, option 2), sim(s) tire un accroissement du
# processus sur une durée s.
resample <- function(t, x, draw, g = regular_grid(t), sim = NULL) {
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
    # inspections strictement entre deux points de grille (méthodes Bridge et MultiBridge)
    while (k < N && t[k + 1] < g[i + 1]) {
      carry <- carry + rest
      k <- k + 1
      t_inf <- t[k]
      x_inf <- x[k]
      rest <- if (k < N) x[k + 1] - x[k] else 0
    }
    if (k == N) {
      # au-delà de la dernière observation : accroissement du processus, sans pont
      s <- sim(g[i + 1] - t_inf)
      dx[i] <- carry + s
      y[i + 1] <- x_inf + s
      t_inf <- g[i + 1]
      x_inf <- y[i + 1]
    } else if (t[k + 1] == g[i + 1]) {
      dx[i] <- carry + rest
      y[i + 1] <- x[k + 1]
      k <- k + 1
      t_inf <- t[k]
      x_inf <- x[k]
      rest <- if (k < N) x[k + 1] - x[k] else 0
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

bridge_paths <- function(x, method, draw, model, parameters, increasing = FALSE, digits = NULL,
                         common_step = FALSE, option = 1, sim = NULL) {
  method <- check_method(method)
  if (!is.logical(common_step) || length(common_step) != 1 || is.na(common_step)) stop("common_step must be TRUE or FALSE")
  if (!identical(option, 1) && !identical(option, 2) && !identical(option, 1L) && !identical(option, 2L)) {
    stop("option must be 1 or 2", call. = FALSE)
  }
  traj <- as_trajectories(x)
  if (increasing) {
    for (k in seq_along(traj)) {
      if (any(diff(traj[[k]]$x) < 0)) stop("trajectory ", k, ": the ", model, " process is increasing, the levels must not decrease")
    }
  }
  if (method == "Bridge") {
    res <- lapply(traj, function(d) resample(d$t, d$x, draw))
    teq <- vapply(traj, function(d) d$t[length(d$t)] / (length(d$t) - 1), 0)
  } else if (method == "MultiBridge") {
    grids <- multi_grids(lapply(traj, `[[`, "t"), option)
    res <- lapply(seq_along(traj), function(k) resample(traj[[k]]$t, traj[[k]]$x, draw, grids[[k]]$grid, sim))
    teq <- vapply(grids, function(g) g$h, 0)
  } else {
    grids <- ext_grids(lapply(traj, `[[`, "t"), digits, common_step)
    for (k in seq_along(traj)) traj[[k]]$t <- grids[[k]]$t
    res <- lapply(seq_along(traj), function(k) resample(traj[[k]]$t, traj[[k]]$x, draw, grids[[k]]$grid))
    teq <- vapply(grids, function(g) g$h, 0)
  }
  observed <- lapply(traj, function(d) data.frame(t = d$t, X = d$x))
  structure(lapply(res, `[[`, "path"), class = "levybridge", model = model, parameters = parameters,
            method = if (method == "MultiBridge") paste0("MultiBridge, option ", option) else method,
            teq = teq, observed = observed, dx = lapply(res, `[[`, "dx"))
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
