devtools::install("C:/Users/ADMIN/Desktop/Ma Thèse/Thèse_1/Recherche_1/levy_packages/levyfit", upgrade = "never")
library(levyfit)

set.seed(1)
t <- inspection_times("increasing", n = 10, tmax = 22.5)$time
t_long <- inspection_times("dec-inc", n = 1000, tmax = 2000)$time
t_mc <- inspection_times("regular", n = 100, tmax = 100)$time
R <- 1000

##############
### WIENER ###
##############

mu <- 0.5
sigma <- 1.5
sim_wiener <- function(t) c(0, cumsum(mu * diff(t) + sigma * sqrt(diff(t)) * rnorm(length(t) - 1)))

x <- sim_wiener(t)
plot(x = t, y = x, type = "b", pch = 19, main = "Wiener")
fit <- wiener_mle(x, t)
fit

x <- sim_wiener(t_long)
fit <- wiener_mle(x, t_long)
c(mu = fit$mu, sigma = fit$sigma)

est <- matrix(NA, R, 2)
couvert <- matrix(NA, R, 2)
for (r in 1:R) {
  fit <- wiener_mle(sim_wiener(t_mc), t_mc)
  est[r, ] <- c(fit$mu, fit$sigma)
  couvert[r, ] <- c(fit$ci_mu[1] <= mu & mu <= fit$ci_mu[2],
                    fit$ci_sigma[1] <= sigma & sigma <= fit$ci_sigma[2])
}
colMeans(est)
colMeans(couvert)

#############
### GAMMA ###
#############

a <- 2
b <- 3
sim_gamma <- function(t) c(0, cumsum(rgamma(length(t) - 1, shape = a * diff(t), rate = b)))

x <- sim_gamma(t)
plot(x = t, y = x, type = "b", pch = 19, main = "Gamma")
fit <- gamma_mle(x, t)
fit

x <- sim_gamma(t_long)
fit <- gamma_mle(x, t_long)
c(a = fit$a, b = fit$b)

est <- matrix(NA, R, 2)
couvert <- matrix(NA, R, 2)
for (r in 1:R) {
  fit <- gamma_mle(sim_gamma(t_mc), t_mc)
  est[r, ] <- c(fit$a, fit$b)
  couvert[r, ] <- abs(c(fit$a, fit$b) - c(a, b)) <= 1.96 * fit$se
}
colMeans(est)
colMeans(couvert)

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

x <- sim_ig(t)
plot(x = t, y = x, type = "b", pch = 19, main = "Inverse gaussien")
fit <- ig_mle(x, t)
fit

x <- sim_ig(t_long)
fit <- ig_mle(x, t_long)
c(mu = fit$mu, eta = fit$eta)

est <- matrix(NA, R, 2)
couvert <- matrix(NA, R, 2)
for (r in 1:R) {
  fit <- ig_mle(sim_ig(t_mc), t_mc)
  est[r, ] <- c(fit$mu, fit$eta)
  couvert[r, ] <- c(fit$ci_mu[1] <= mu & mu <= fit$ci_mu[2],
                    fit$ci_eta[1] <= eta & eta <= fit$ci_eta[2])
}
colMeans(est)
colMeans(couvert)

###############################
### DONNEES NON VALABLES ###
###############################

try(wiener_mle(c(0, 1), c(0, 1)))
try(gamma_mle(c(0, 1, 1, 2), c(0, 1, 2, 3)))
try(gamma_mle(c(0, 1, 2, 3), c(0, 1, 2, 3)))
try(ig_mle(c(0, 1, 2), c(0, 1, 1)))
