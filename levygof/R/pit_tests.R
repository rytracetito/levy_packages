# Tests KS, CvM, AD et Watson avec paramètres estimés (VG, Meixner) et test VGE pour la loi gamma :
# p-valeur par bootstrap paramétrique.

stat_ks <- function(p) {
  m <- length(p$F); j <- seq_len(m)
  max(j / m - p$F, p$F - (j - 1) / m)
}

stat_cvm <- function(p) {
  m <- length(p$F); j <- seq_len(m)
  sum((p$F - (2 * j - 1) / (2 * m))^2) + 1 / (12 * m)
}

stat_ad <- function(p) {
  m <- length(p$F); j <- seq_len(m)
  -m - mean((2 * j - 1) * (log(p$F) + log(rev(p$S))))
}

stat_watson <- function(p) {
  stat_cvm(p) - length(p$F) * (mean(p$F) - 0.5)^2
}

pit_stats <- list(KS = stat_ks, CvM = stat_cvm, AD = stat_ad, Watson = stat_watson)

# F(y_(j)) et 1 - F(y_(j)) par intégration de la densité entre points triés consécutifs ;
# l'intégrale est coupée en `split` (singularité possible de la densité VG en 0)
pit_numeric <- function(y, dens, split = NULL) {
  y <- sort(y)
  n <- length(y)
  I <- function(lo, hi) {
    if (!is.null(split) && lo < split && split < hi) return(I(lo, split) + I(split, hi))
    stats::integrate(dens, lo, hi, rel.tol = 1e-8, subdivisions = 500L, stop.on.error = FALSE)$value
  }
  mid <- if (n > 1) vapply(2:n, function(j) if (y[j] > y[j - 1]) I(y[j - 1], y[j]) else 0, 0) else numeric(0)
  F <- pmax(I(-Inf, y[1]) + c(0, cumsum(mid)), 1e-300)
  S <- pmax(I(y[n], Inf) + rev(c(0, cumsum(rev(mid)))), 1e-300)
  list(F = F / (F + S), S = S / (F + S))
}

dvg_log <- function(x, theta, sigma, nu) {
  lam <- 1 / nu
  cc <- sqrt(theta^2 + 2 * sigma^2 / nu)
  ax <- pmax(abs(x), 1e-300)
  z <- cc * ax / sigma^2
  log(2) + theta * x / sigma^2 - 0.5 * log(2 * pi) - log(sigma) - lam * log(nu) - lgamma(lam) +
    (lam - 0.5) * (log(ax) - log(cc)) + log(besselK(z, abs(lam - 0.5), expon.scaled = TRUE)) - z
}

dmeixner_log <- function(x, a, b, d, m) {
  y <- (x - m) / a
  2 * d * log(2 * cos(b / 2)) - log(2 * a * pi) - lgamma(2 * d) + b * y +
    2 * Re(gsl::lngamma_complex(rep_len(d, length(y)), y))
}

# Lois sur un pas 1 : la famille est stable par changement de pas, le pas des accroissements est sans importance
# Ajustement sur les accroissements eux-mêmes (increments_data) : aucune perte de précision
families <- list(
  VG = list(
    fit = function(y) {
      p <- suppressWarnings(levyfit::vg_mle(levyfit::increments_data(y)))
      if (anyNA(p)) stop(attr(p, "message"), call. = FALSE)
      p
    },
    pit = function(y, p) pit_numeric(y, function(x) exp(dvg_log(x, p[["theta"]], p[["sigma"]], p[["nu"]])), split = 0),
    sim = function(n, p) {
      g <- stats::rgamma(n, 1 / p[["nu"]], scale = p[["nu"]])
      p[["theta"]] * g + p[["sigma"]] * sqrt(g) * stats::rnorm(n)
    }),
  Meixner = list(
    fit = function(y) {
      p <- suppressWarnings(levyfit::meixner_mle(levyfit::increments_data(y)))
      if (anyNA(p)) stop(attr(p, "message"), call. = FALSE)
      p
    },
    pit = function(y, p) pit_numeric(y, function(x) exp(dmeixner_log(x, p[["a"]], p[["b"]], p[["d"]], p[["m"]]))),
    sim = function(n, p) {
      gen <- Runuran::pinvd.new(Runuran::udmeixner(alpha = p[["a"]], beta = p[["b"]], delta = p[["d"]], mu = p[["m"]]))
      Runuran::ur(gen, n)
    })
)

boot_pvalue <- function(T_obs, T_boot) {
  ok <- !is.na(T_boot)
  if (!any(ok)) stop("all the bootstrap replications failed", call. = FALSE)
  (1 + sum(T_boot[ok] >= T_obs)) / (1 + sum(ok))
}

# Test `test` (KS, CvM, AD ou Watson) de la famille `model` (VG ou Meixner)
gof_pit <- function(dx, model, test, B) {
  fam <- families[[model]]
  stat <- pit_stats[[test]]
  par <- fam$fit(dx)
  T_obs <- stat(fam$pit(dx, par))
  T_boot <- vapply(seq_len(B), function(i) {
    y <- fam$sim(length(dx), par)
    tryCatch(stat(fam$pit(y, fam$fit(y))), error = function(e) NA_real_)
  }, 0)
  list(statistic = T_obs, p.value = boot_pvalue(T_obs, T_boot))
}

# Test VGE de la loi gamma : statistique de goft::gamma_test, p-valeur par bootstrap paramétrique
gof_gamma_vge <- function(dx, B) {
  fit <- goft::gamma_fit(dx)
  shape <- fit["shape", 1]
  scale <- fit["scale", 1]
  stat <- function(y) abs(as.numeric(suppressWarnings(goft::gamma_test(y))$statistic))
  T_obs <- stat(dx)
  T_boot <- vapply(seq_len(B), function(i) tryCatch(stat(stats::rgamma(length(dx), shape = shape, scale = scale)),
                                                   error = function(e) NA_real_), 0)
  list(statistic = T_obs, p.value = boot_pvalue(T_obs, T_boot))
}

# Test omnibus de D'Agostino-Pearson (K^2)
gof_dagostino_pearson <- function(dx) {
  r <- fBasics::dagoTest(dx)
  list(statistic = unname(r@test$statistic[1]), p.value = unname(r@test$p.value[1]))
}
