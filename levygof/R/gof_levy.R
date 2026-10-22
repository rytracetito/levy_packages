#' Goodness-of-fit test of a Lévy degradation process
#'
#' Tests whether one or several degradation trajectories come from a given
#' Lévy process (gamma, Wiener, inverse Gaussian, variance-gamma or Meixner),
#' with unknown parameters. With regular inspections the observed increments
#' are tested; with irregular inspections the parameters are estimated, the
#' trajectories are resampled on a regular grid by stochastic bridges, and
#' the reconstructed increments are tested.
#'
#' @param data The observed trajectories, in one of the formats of the
#'   packages `levyfit` and `levyBridge`:
#'   * one trajectory: a data frame whose first column contains the
#'     inspection times and second column the degradation levels;
#'   * several trajectories, each with its own inspection times: a list of
#'     such data frames;
#'   * several trajectories with the same inspection times: a matrix whose
#'     first column contains the times and the other columns the levels, for
#'     instance the result of a simulation function of `levysim`.
#'
#'   The columns are read by position, not by name. The inspection times of
#'   a trajectory must be strictly increasing; they do not have to start at
#'   0 (the first inspection is taken as origin).
#' @param model `"HGP"` (gamma process), `"HWP"` (Wiener process),
#'   `"InvGaus"` (inverse Gaussian process), `"VG"` (variance-gamma process)
#'   or `"Meixner"` (Meixner process). Upper or lower case does not matter.
#' @param test Names of the tests to apply; `NULL` (default) applies all the
#'   tests of the model, see [tests_disponibles()].
#' @param method Resampling method of the bridges, used with irregular
#'   inspections only: `"Bridge"` (default, Grall-Maës, 2012: regular grid
#'   with as many points as observations) or `"ExtBridge"` (regular grid
#'   containing all the inspection times, with a time step common to all the
#'   trajectories), see [levyBridge::gamma_bridge()].
#' @param alpha Level of the tests, strictly between 0 and 1.
#' @param B Number of bootstrap replications of the tests that use a
#'   parametric bootstrap.
#'
#' @return An object of class `"gof.levy"`, printed by `print()`: a list
#'   with
#'   * `model`, `method` (`NA` when no bridge was used), `alpha`, `B`;
#'   * `regular`: `TRUE` when all the trajectories are inspected with a
#'     constant time step;
#'   * `out_of_support`: `TRUE` when an observed increment is negative or
#'     zero under an increasing model (gamma, inverse Gaussian);
#'   * `parameters`: maximum likelihood estimates of the parameters of the
#'     process (`NULL` when no bridge was used);
#'   * `steps`: time step of the increments of each trajectory;
#'   * `bridge`: the resampled trajectories (list of data frames with
#'     columns `t` and `X`, on the time scale of the data), or `NULL`;
#'   * `increments`: data frame of the tested increments, with columns `id`
#'     (trajectory), `step` and `dx`;
#'   * `results`: data frame of the tests, with a column `group` (`"all"`
#'     for the test on all the trajectories together, or the number of the
#'     trajectory) followed by the columns of [test_adequation()].
#'
#' @section Method:
#' The null hypothesis H0 is: the trajectories are independent trajectories
#' of the Lévy process `model`, with unknown parameters.
#'
#' 1. **Regular inspections** (each trajectory inspected with a constant
#'    time step \eqn{\Delta}): the increments are i.i.d. with the law of
#'    \eqn{X_\Delta}; they are tested directly.
#' 2. **Irregular inspections**: the parameters are estimated by maximum
#'    likelihood on all the trajectories (package `levyfit`), then each
#'    trajectory is resampled on a regular grid by exact simulation of the
#'    bridge of the process with these parameters (package `levyBridge`,
#'    method `method`). Under H0 the reconstructed increments are i.i.d. with
#'    the law of \eqn{X_{t_{eq}}}, where \eqn{t_{eq}} is the time step of the
#'    grid; they are tested.
#' 3. **Gamma and inverse Gaussian models**: these processes are
#'    increasing. A negative or zero observed increment is impossible under
#'    H0, which is then rejected for sure (p-value 0), without estimation
#'    nor bridge.
#'
#' The tests are tests for i.i.d. samples with estimated parameters, see
#' [test_adequation()] and [tests_disponibles()].
#'
#' @section Several trajectories:
#' Each trajectory is tested alone (group `1`, `2`, ...). When all the
#' trajectories have the same time step, their increments have the same law
#' under H0: they are also tested together (group `"all"`). This is always
#' the case with the method `"ExtBridge"`, whose grids have a step common to
#' all the trajectories (the largest step such that every inspection time of
#' every trajectory is a grid point). With the method `"Bridge"` the step
#' \eqn{t_N / N} of each trajectory is generally different, and the
#' increments of different trajectories are then not tested together.
#'
#' @section Randomness:
#' The bridges and the bootstrap are random: use [set.seed()] before
#' `gof.levy()` to reproduce a result.
#'
#' @references
#' Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma
#' process. *Proceedings of the Institution of Mechanical Engineers, Part O:
#' Journal of Risk and Reliability*, 226(6), 624-634.
#' \doi{10.1177/1748006X12462522}
#'
#' @seealso [test_adequation()] for an i.i.d. sample of increments,
#'   [tests_disponibles()] for the list of the tests,
#'   [levyfit::gamma_mle()] for the estimation,
#'   [levyBridge::gamma_bridge()] for the bridges.
#'
#' @examples
#' set.seed(1)
#'
#' # 1. One gamma trajectory (a = 2, b = 3) inspected every 2 time units:
#' #    the observed increments are tested directly
#' t <- seq(0, 60, by = 2)
#' x <- c(0, cumsum(rgamma(30, shape = 2 * 2, rate = 3)))
#' gof.levy(data.frame(t, x), "HGP", B = 100)
#'
#' # 2. Irregular inspections (a = 10, b = 5): estimation, gamma bridge, test
#' t2 <- c(0, cumsum(round(runif(30, 0.5, 3.5), 1)))
#' x2 <- c(0, cumsum(rgamma(30, shape = 10 * diff(t2), rate = 5)))
#' res <- gof.levy(data.frame(t2, x2), "HGP", test = c("KS", "AD"), B = 100)
#' res
#' res$parameters
#' res$results
#'
#' # 3. The same trajectory with the method "ExtBridge"
#' gof.levy(data.frame(t2, x2), "HGP", test = c("KS", "AD"), method = "ExtBridge", B = 100)
#'
#' # 3b. A small shape (a = 2): the reconstructed increments can be tiny but stay positive
#' x4 <- c(0, cumsum(rgamma(30, shape = 2 * diff(t2), rate = 3)))
#' gof.levy(data.frame(t2, x4), "HGP", test = "AD", method = "ExtBridge", B = 100)
#'
#' # 4. Wiener model on the same data, at level 10 %
#' gof.levy(data.frame(t2, x2), "hwp", alpha = 0.1)
#'
#' # 5. A decreasing step: the gamma model is rejected for sure
#' x3 <- x2
#' x3[10] <- x3[9] - 0.1
#' gof.levy(data.frame(t2, x3), "HGP", test = "KS")
#'
#' # 6. Several trajectories with their own inspection times (a list):
#' #    one test per trajectory
#' liste <- lapply(1:3, function(k) {
#'   tk <- c(0, cumsum(round(runif(20, 0.5, 3.5), 1)))
#'   data.frame(t = tk, X = c(0, cumsum(rnorm(20, 0.5 * diff(tk), 1.5 * sqrt(diff(tk))))))
#' })
#' gof.levy(liste, "HWP", test = c("SW", "AD"))
#'
#' # 6b. The same trajectories with "ExtBridge": common step, so also one test on all
#' gof.levy(liste, "HWP", test = c("SW", "AD"), method = "ExtBridge")
#'
#' # 7. Several trajectories with the same irregular inspection times (matrix):
#' #    the grids have the same step, so the increments are also tested together
#' m <- cbind(t = t2, sapply(1:4, function(k) c(0, cumsum(rgamma(30, 10 * diff(t2), 5)))))
#' gof.levy(m, "HGP", test = "AD", B = 100)
#'
#' \donttest{
#' # 8. Inverse Gaussian and variance-gamma models (slower)
#' ig <- c(0, cumsum(statmod::rinvgauss(30, mean = diff(t2), shape = 3 * diff(t2)^2)))
#' gof.levy(data.frame(t2, ig), "InvGaus", test = c("KS", "AD", "HK2"), B = 100)
#'
#' g <- rgamma(30, shape = diff(t2) / 0.2, scale = 0.2)
#' vg <- c(0, cumsum(0.5 * g + 1.5 * sqrt(g) * rnorm(30)))
#' gof.levy(data.frame(t2, vg), "VG", test = "AD", B = 50)
#' }
#' @export
gof.levy <- function(data, model, test = NULL, method = "Bridge", alpha = 0.05, B = 200) {
  model <- model_name(model)
  test <- test_names(test, model)
  check_alpha_B(alpha, B)
  method <- match.arg(method, c("Bridge", "ExtBridge"))
  traj <- as_trajectories(data)

  dx_obs <- lapply(traj, function(d) diff(d$x))
  regular <- all(vapply(traj, function(d) is_regular(d$t), TRUE))
  out_of_support <- model %in% c("HGP", "InvGaus") && any(unlist(dx_obs) <= 0)
  par <- NULL
  bridge <- NULL

  if (out_of_support) {
    pieces <- list(list(step = NA_real_, dx = unlist(dx_obs)))
  } else if (regular) {
    pieces <- lapply(seq_along(traj), function(k) list(step = diff(traj[[k]]$t)[1], dx = dx_obs[[k]]))
  } else {
    par <- estimate_parameters(traj, model)
    r <- bridge_trajectories(traj, model, par, method)
    teq <- attr(r, "teq")
    inc <- levyBridge::increments(r)
    pieces <- lapply(seq_along(r), function(k) list(step = teq[k], dx = inc$dx[inc$id == k]))
    bridge <- lapply(seq_along(r), function(k) data.frame(t = r[[k]]$t + traj[[k]]$t0, X = r[[k]]$X))
    check_zero_increments(unlist(lapply(pieces, `[[`, "dx")), model, method)
  }

  structure(list(model = model, method = if (regular || out_of_support) NA_character_ else method,
                 alpha = alpha, B = B, regular = regular, out_of_support = out_of_support,
                 parameters = par, steps = vapply(pieces, `[[`, 0, "step"), bridge = bridge,
                 increments = do.call(rbind, lapply(seq_along(pieces), function(k)
                   data.frame(id = k, step = pieces[[k]]$step, dx = pieces[[k]]$dx))),
                 results = run_tests(pieces, model, test, alpha, B, out_of_support)),
            class = "gof.levy")
}

as_trajectories <- function(data) {
  if (is.matrix(data)) {
    if (ncol(data) < 2) stop("the matrix data must have the times in column 1 and the levels in the next columns", call. = FALSE)
    data <- lapply(seq_len(ncol(data) - 1) + 1, function(k) data.frame(t = data[, 1], X = data[, k]))
  }
  if (is.data.frame(data)) data <- list(data)
  if (!is.list(data) || length(data) == 0 || !all(vapply(data, is.data.frame, TRUE))) {
    stop("data must be a data frame (times, levels), a list of such data frames, or a matrix (times, levels...)", call. = FALSE)
  }
  traj <- lapply(seq_along(data), function(k) {
    d <- data[[k]]
    if (ncol(d) < 2) stop("trajectory ", k, ": the data frame must have the times in column 1 and the levels in column 2", call. = FALSE)
    t <- d[[1]]
    x <- d[[2]]
    if (!is.numeric(t) || !is.numeric(x)) stop("trajectory ", k, ": times and levels must be numeric", call. = FALSE)
    if (!all(is.finite(t)) || !all(is.finite(x))) stop("trajectory ", k, ": missing or infinite values", call. = FALSE)
    if (length(t) < 2) stop("trajectory ", k, ": at least two inspections are needed", call. = FALSE)
    if (any(diff(t) <= 0)) stop("trajectory ", k, ": the inspection times must be strictly increasing", call. = FALSE)
    list(t = t - t[1], x = x, t0 = t[1])
  })
  if (sum(vapply(traj, function(d) length(d$t) - 1, 0)) < 4) stop("at least 4 increments are needed", call. = FALSE)
  traj
}

is_regular <- function(t) {
  dt <- diff(t)
  all(abs(dt - dt[1]) <= 1e-6 * dt[1])
}

as_data_frames <- function(traj) lapply(traj, function(d) data.frame(t = d$t, X = d$x))

estimate_parameters <- function(traj, model) {
  d <- as_data_frames(traj)
  p <- suppressWarnings(switch(model,
    HGP     = levyfit::gamma_mle(d),
    HWP     = levyfit::wiener_mle(d),
    InvGaus = levyfit::ig_mle(d),
    VG      = levyfit::vg_mle(d),
    Meixner = levyfit::meixner_mle(d)
  ))
  if (anyNA(p)) {
    stop("estimation impossible for the model ", model, ", the bridge cannot be built: ", attr(p, "message"), call. = FALSE)
  }
  switch(model,
    HGP     = p[c("a", "b")],
    HWP     = c(mu = p[["mu"]], sigma = sqrt(p[["sigma2"]])),
    InvGaus = p[c("mu", "eta")],
    VG      = p[c("theta", "sigma", "nu")],
    Meixner = p[c("a", "b", "d", "m")]
  )
}

# Avec ExtBridge, pas commun à toutes les trajectoires : tous les accroissements reconstruits ont la même loi
bridge_trajectories <- function(traj, model, par, method) {
  d <- as_data_frames(traj)
  cs <- method == "ExtBridge"
  switch(model,
    HGP     = levyBridge::gamma_bridge(d, a = par[["a"]], method = method, common_step = cs),
    HWP     = levyBridge::wiener_bridge(d, sigma = par[["sigma"]], method = method, common_step = cs),
    InvGaus = levyBridge::ig_bridge(d, eta = par[["eta"]], method = method, common_step = cs),
    VG      = levyBridge::vg_bridge(d, theta = par[["theta"]], sigma = par[["sigma"]], nu = par[["nu"]], method = method,
                                    common_step = cs),
    Meixner = levyBridge::meixner_bridge(d, a = par[["a"]], d = par[["d"]], m = par[["m"]], method = method,
                                         common_step = cs)
  )
}

# Accroissement reconstruit nul : impossible sous ces modèles (dépassement de capacité extrême)
check_zero_increments <- function(dx, model, method) {
  if (any(dx == 0) && model %in% c("HGP", "InvGaus", "VG")) {
    stop(sum(dx == 0), " reconstructed increment(s) equal to 0 (below the smallest positive number): ",
         "the test cannot be computed", call. = FALSE)
  }
}

run_tests <- function(pieces, model, test, alpha, B, out_of_support) {
  safe <- function(dx, group) {
    r <- tryCatch(test_adequation(dx, model, test, alpha, B), error = function(e) e)
    if (inherits(r, "error")) r <- do.call(rbind, lapply(test, function(nom) result_row(nom, NA_real_, NA_real_, NA_real_, NA, conditionMessage(r))))
    cbind(group = group, r, stringsAsFactors = FALSE)
  }
  if (out_of_support) return(safe(pieces[[1]]$dx, "all"))
  steps <- vapply(pieces, `[[`, 0, "step")
  common <- all(abs(steps - steps[1]) <= 1e-8 * steps[1])
  res <- lapply(seq_along(pieces), function(k) safe(pieces[[k]]$dx, as.character(k)))
  if (length(pieces) > 1 && common) res <- c(list(safe(unlist(lapply(pieces, `[[`, "dx")), "all")), res)
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

#' @export
print.gof.levy <- function(x, ...) {
  n <- length(x$steps)
  cat("\nGoodness-of-fit test of the model", x$model, "\n")
  cat("Trajectories:", n, "| Increments tested:", nrow(x$increments), "| Level:", x$alpha, "| B:", x$B, "\n")
  cat("Inspections:",
      if (x$out_of_support) "a negative or zero increment is observed under an increasing model: H0 rejected for sure"
      else if (x$regular) "regular (observed increments)"
      else paste0("irregular (bridge, method \"", x$method, "\")"), "\n")
  if (!x$out_of_support) {
    st <- signif(unique(x$steps), 4)
    cat("Time step:", if (length(st) <= 4) paste(st, collapse = ", ") else paste("from", min(st), "to", max(st)), "\n")
  }
  if (!is.null(x$parameters)) {
    cat("Estimated parameters of the process:", paste(names(x$parameters), signif(x$parameters, 4), sep = " = ", collapse = ", "), "\n")
  }
  cat("\n")
  r <- x$results
  for (v in c("statistic", "p.value", "critical.value")) r[[v]] <- signif(r[[v]], 4)
  if (all(is.na(r$critical.value))) r$critical.value <- NULL
  if (all(r$message == "")) r$message <- NULL
  if (n == 1) r$group <- NULL
  print(r, row.names = FALSE)
  invisible(x)
}
