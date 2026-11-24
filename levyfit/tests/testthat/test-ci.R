test_that("ci = TRUE gives the closed-form standard errors of the Wiener MLE", {
  set.seed(1)
  t <- seq(0, 100, by = 1)
  x <- c(0, cumsum(0.5 * diff(t) + 1.5 * sqrt(diff(t)) * rnorm(100)))
  r <- wiener_mle(x, t, ci = TRUE)
  f <- wiener_mle(x, t)
  z <- qnorm(0.975)
  expect_equal(r$upper[1] - r$estimate[1], z * sqrt(f[["sigma2"]] / 100), tolerance = 1e-4)
  expect_equal(log(r$upper[2] / r$estimate[2]), z * sqrt(2 / 100), tolerance = 1e-4)
  expect_true(is.na(r$lower[3]))
})

test_that("ci = TRUE intervals contain the estimate and widen with the level", {
  set.seed(2)
  t <- seq(0, 100, by = 2)
  x <- c(0, cumsum(rgamma(50, shape = 2 * diff(t), rate = 3)))
  r95 <- gamma_mle(x, t, ci = TRUE)
  r99 <- gamma_mle(x, t, ci = TRUE, level = 0.99)
  expect_equal(r95$estimate, unname(gamma_mle(x, t)))
  expect_true(all(r95$lower < r95$estimate & r95$estimate < r95$upper & r95$lower > 0))
  expect_true(all(r99$lower < r95$lower & r95$upper < r99$upper))
  expect_error(gamma_mle(x, t, ci = TRUE, level = 1.2))
})

test_that("without ci the estimates are a named vector, with ci a data frame", {
  set.seed(3)
  t <- seq(0, 50, by = 1)
  x <- c(0, cumsum(rgamma(50, shape = 2 * diff(t), rate = 3)))
  expect_named(gamma_mle(x, t), c("a", "b"))
  r <- gamma_mle(x, t, ci = TRUE)
  expect_s3_class(r, "data.frame")
  expect_named(r, c("parameter", "estimate", "lower", "upper"))
  expect_equal(r$estimate, unname(gamma_mle(x, t)))
  expect_error(gamma_mle(x, t, ci = "oui"), "ci must be")
  expect_error(gamma_mle(x, t, ci = TRUE, level = 0), "level")
})

test_that("an estimator that does not exist gives a warning and NA, not an error", {
  expect_warning(r <- wiener_mle(c(0, 1), c(0, 1)), "does not exist")
  expect_true(all(is.na(r)))
  expect_named(r, c("mu", "sigma2", "s2"))
  expect_match(attr(r, "message"), "does not exist")
  expect_warning(r <- gamma_mle(c(0, 1, 2, 3), c(0, 1, 2, 3), ci = TRUE), "does not exist")
  expect_true(all(is.na(r$estimate)) && all(is.na(r$lower)))
  set.seed(4)
  t <- seq(0, 100, by = 1)
  y <- c(0, cumsum(rgamma(100, shape = 2 * diff(t), rate = 3)))
  expect_warning(r <- vg_mle(y, t), "gamma limit")
  expect_true(all(is.na(r)))
  expect_error(gamma_mle(c(0, 2, 1.5, 3), c(0, 1, 2, 3)), "strictly positive")
})

test_that("positive parameters have intervals on the log scale, real ones symmetric intervals", {
  set.seed(4)
  t <- seq(0, 100, by = 1)
  x <- c(0, cumsum(rgamma(100, shape = 1.5 * diff(t), rate = 3)))
  r <- gamma_mle(x, t, ci = TRUE)
  expect_equal(r$lower * r$upper, r$estimate^2, tolerance = 1e-8)
  w <- wiener_mle(x, t, ci = TRUE)
  expect_equal(w$upper[1] - w$estimate[1], w$estimate[1] - w$lower[1], tolerance = 1e-8)
})

test_that("ci = TRUE works when the parameters are much smaller than 1e-3", {
  # Ordre de grandeur des données GaAs laser : mu = 0.002, eta = 5e-5, pas de 250 h
  set.seed(5)
  t <- seq(0, 4000, by = 250)
  m <- 2.04e-3 * diff(t)
  la <- 5.46e-5 * diff(t)^2
  liste <- lapply(1:15, function(i) {
    nu <- rnorm(16)^2
    x1 <- m + m^2 * nu / (2 * la) - m / (2 * la) * sqrt(4 * m * la * nu + m^2 * nu^2)
    dx <- ifelse(runif(16) <= m / (m + x1), x1, m^2 / x1)
    data.frame(t = t, X = c(0, cumsum(dx)))
  })
  f <- ig_mle(liste)
  r <- expect_no_warning(ig_mle(liste, ci = TRUE))
  z <- qnorm(0.975)
  T <- 15 * 4000
  expect_equal(log(r$upper[1] / r$estimate[1]), z * sqrt(f[["mu"]] / (f[["eta"]] * T)), tolerance = 1e-3)
  expect_equal(log(r$upper[2] / r$estimate[2]), z * sqrt(2 / 240), tolerance = 1e-3)

  # Wiener à petite échelle : mêmes intervalles relatifs que les données multipliées par 1000
  set.seed(6)
  x <- c(0, cumsum(0.5 * diff(seq(0, 100)) + 1.5 * rnorm(100)))
  r_small <- wiener_mle(x / 1000, 0:100, ci = TRUE)
  r_big <- wiener_mle(x, 0:100, ci = TRUE)
  expect_equal(r_small$lower[1:2] * c(1000, 1e6), r_big$lower[1:2], tolerance = 1e-4)
  expect_equal(r_small$upper[1:2] * c(1000, 1e6), r_big$upper[1:2], tolerance = 1e-4)
})
