as_trajectories <- function(x, t) {
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
  traj <- as_trajectories(x, t)
  list(dt = unlist(lapply(traj, function(d) diff(d$t))),
       dx = unlist(lapply(traj, function(d) diff(d$x))))
}
