remotes::install_github("rytracetito/levy_packages", subdir = "levyfit")
library(levyfit)

set.seed(1)
t1 <- inspection_times("increasing", n = 10, tmax = 22.5)$time
t2 <- inspection_times("decreasing", n = 10, tmax = 22.5)$time
t3 <- inspection_times("dec-inc", n = 10, tmax = 2000)$time
t4 <- inspection_times("regular", n = 10, tmax = 100)$time
R <- 1000

plot(x = t1, y = rep(1, length(t1)), pch = 19, ylab = "", main = "increasing, n = 10")
plot(x = t2, y = rep(1, length(t2)), pch = 19, ylab = "", main = "decreasing, n = 10")
plot(x = t3, y = rep(1, length(t3)), pch = 19, ylab = "", main = "dec-inc, n = 10")
plot(x = t4, y = rep(1, length(t4)), pch = 19, ylab = "", main = "regular, n = 10")

##############
### WIENER ###
##############

mu <- 0.5
sigma <- 1.5
sim_wiener <- function(t) c(0, cumsum(mu * diff(t) + sigma * sqrt(diff(t)) * rnorm(length(t) - 1)))

x <- sim_wiener(t1)
plot(x = t1, y = x, type = "b", pch = 19, main = "Wiener")
wiener_mle(x, t1)

x <- sim_wiener(t3)
wiener_mle(x, t3)

est <- matrix(NA, R, 2)
for (r in 1:R) est[r, ] <- wiener_mle(sim_wiener(t4), t4)
colMeans(est)
c(mu, sigma^2)

#############
### GAMMA ###
#############

a <- 2
b <- 3
sim_gamma <- function(t) c(0, cumsum(rgamma(length(t) - 1, shape = a * diff(t), rate = b)))

x <- sim_gamma(t1)
plot(x = t1, y = x, type = "b", pch = 19, main = "Gamma")
gamma_mle(x, t1)

x <- sim_gamma(t3)
gamma_mle(x, t3)

est <- matrix(NA, R, 2)
for (r in 1:R) est[r, ] <- gamma_mle(sim_gamma(t4), t4)
colMeans(est)
c(a, b)

########################
### INVERSE GAUSSIEN ###
########################

mu <- 0.5
eta <- 2
rig <- function(n, mean, shape) {
  nu <- rnorm(n)^2
  y <- mean + mean^2 * nu / (2 * shape) - mean / (2 * shape) * sqrt(4 * mean * shape * nu + mean^2 * nu^2)
  ifelse(runif(n) <= mean / (mean + y), y, mean^2 / y)
}
sim_ig <- function(t) c(0, cumsum(rig(length(t) - 1, mu * diff(t), eta * diff(t)^2)))

x <- sim_ig(t1)
plot(x = t1, y = x, type = "b", pch = 19, main = "Inverse gaussien")
ig_mle(x, t1)

x <- sim_ig(t3)
ig_mle(x, t3)

est <- matrix(NA, R, 2)
for (r in 1:R) est[r, ] <- ig_mle(sim_ig(t4), t4)
colMeans(est)
c(mu, eta)

###############################
### DONNEES NON VALABLES ###
###############################

try(wiener_mle(c(0, 1), c(0, 1)))
try(gamma_mle(c(0, 1, 1, 2), c(0, 1, 2, 3)))
try(gamma_mle(c(0, 1, 2, 3), c(0, 1, 2, 3)))
try(ig_mle(c(0, 1, 2), c(0, 1, 1)))
