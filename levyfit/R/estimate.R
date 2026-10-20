no_mle <- function(message) {
  stop(structure(class = c("no_mle", "error", "condition"), list(message = message, call = NULL)))
}

parameter_names <- list(gamma = c("a", "b"), wiener = c("mu", "sigma2", "s2"),
                        ig = c("mu", "eta", "eta_unbiased"), vg = c("theta", "sigma", "nu"),
                        meixner = c("a", "b", "d", "m"), igam = c("theta", "alpha", "beta"))

estimate <- function(model, x, t, ci, level, ...) {
  if (!is.logical(ci) || length(ci) != 1 || is.na(ci)) stop("ci must be TRUE or FALSE")
  if (!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1) {
    stop("level must be a number between 0 and 1")
  }
  fit <- tryCatch(switch(model, gamma = gamma_est(x, t), wiener = wiener_est(x, t), ig = ig_est(x, t),
                         vg = vg_est(x, t), meixner = meixner_est(x, t, ...), igam = igam_est(x, t)),
                  no_mle = function(e) e)
  if (inherits(fit, "no_mle")) {
    warning("estimation impossible: ", conditionMessage(fit), call. = FALSE)
    out <- stats::setNames(rep(NA_real_, length(parameter_names[[model]])), parameter_names[[model]])
    if (ci) out <- data.frame(parameter = names(out), estimate = NA_real_, lower = NA_real_, upper = NA_real_)
    attr(out, "message") <- conditionMessage(fit)
    return(out)
  }
  if (!ci) return(fit)
  confidence_intervals(model, x, t, fit, level)
}

confidence_intervals <- function(model, x, t, fit, level) {
  par <- switch(model, wiener = fit[c("mu", "sigma2")], ig = fit[c("mu", "eta")], fit)
  real <- names(par) %in% switch(model, wiener = "mu", vg = "theta", meixner = c("b", "m"), character(0))

  ll <- loglik_fun(x, t, model)
  H <- optimHess(par, function(p) -ll(p), control = list(parscale = pmax(abs(par), 1e-8)))
  V <- tryCatch(solve(H), error = function(e) matrix(NA, length(par), length(par)))
  se <- sqrt(ifelse(diag(V) > 0, diag(V), NA))

  z <- qnorm((1 + level) / 2)
  lower <- ifelse(real, par - z * se, par * exp(-z * se / par))
  upper <- ifelse(real, par + z * se, par * exp(z * se / par))
  out <- data.frame(parameter = names(fit), estimate = unname(fit), lower = NA_real_, upper = NA_real_)
  out$lower[match(names(par), out$parameter)] <- lower
  out$upper[match(names(par), out$parameter)] <- upper
  out
}

loglik_fun <- function(x, t, model) {
  if (model == "igam") {
    traj <- as_trajectories(x, t)
    W <- lapply(traj, function(d) d$x)
    tt <- lapply(traj, function(d) d$t)
    return(function(p) igam_loglik(W, tt, p[1], p[2], p[3]))
  }
  d <- increments(x, t)
  dt <- d$dt
  dx <- d$dx
  switch(model,
    gamma = function(p) sum(dgamma(dx, shape = p[1] * dt, rate = p[2], log = TRUE)),
    wiener = function(p) sum(dnorm(dx, p[1] * dt, sqrt(p[2] * dt), log = TRUE)),
    ig = function(p) {
      mu <- p[1] * dt
      la <- p[2] * dt^2
      sum(0.5 * log(la / (2 * pi * dx^3)) - la * (dx - mu)^2 / (2 * mu^2 * dx))
    },
    vg = function(p) vg_loglik(dx, dt, p[1], p[2], p[3]),
    meixner = function(p) meixner_loglik(dx, dt, p[1], p[2], p[3], p[4]))
}
