incr <- function(r) {
  do.call(rbind, lapply(seq_along(r), function(i) data.frame(id = i, teq = attr(r, "teq")[i], dx = attr(r, "dx")[[i]])))
}
