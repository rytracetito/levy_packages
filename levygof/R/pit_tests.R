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

# log(K_p(z) exp(z)) par récurrence ascendante sur l'ordre, en logarithmes :
# K_{v+1}(z) = K_{v-1}(z) + (2 v / z) K_v(z), à partir de l'ordre fractionnaire v0 = p - floor(p)
log_besselKx_rec <- function(z, p) {
  n0 <- floor(p)
  v0 <- p - n0
  l <- log(besselK(z, v0, expon.scaled = TRUE))
  r <- besselK(z, v0 + 1, expon.scaled = TRUE) / besselK(z, v0, expon.scaled = TRUE)
  for (k in seq_len(max(n0))) {
    i <- n0 >= k
    l[i] <- l[i] + log(r[i])
    r <- 1 / r + 2 * (v0 + k) / z
  }
  l
}

# log(K_p(z) exp(z)) sans dépassement de capacité (même fonction que dans levyfit) : besselK si possible,
# sinon développement de Debye (grand ordre) ou intégrale de la représentation K_p(z) = int exp(-z cosh u) cosh(pu) du
log_besselKx <- function(z, p) {
  n <- max(length(z), length(p))
  z <- rep_len(z, n)
  p <- rep_len(abs(p), n)
  out <- rep(NA_real_, n)
  small <- p <= 1e4
  out[small] <- log(besselK(z[small], p[small], expon.scaled = TRUE))
  debye <- function(z, p) {
    s <- sqrt(1 + (z / p)^2)
    0.5 * log(pi / (2 * p)) - p * log((z / p) / (1 + s)) - 0.5 * log(s) - p^2 / (z + p * s)
  }
  # débordement de besselK (grand ordre, petit z) : récurrence exacte sur l'ordre
  bad <- which(!is.finite(out) & is.finite(z) & z > 0 & p <= 1e4)
  if (length(bad) > 0) out[bad] <- log_besselKx_rec(z[bad], p[bad])
  for (i in which(!is.finite(out) & is.finite(z) & z > 0 & is.finite(p))) {
    out[i] <- if (p[i] > 1e4) debye(z[i], p[i]) else {
      g <- function(u) -2 * z[i] * sinh(u / 2)^2 + p[i] * u + log1p(exp(-2 * p[i] * u)) - log(2)
      opt <- stats::optimize(g, c(0, asinh(p[i] / z[i]) + 1), maximum = TRUE, tol = 1e-12)
      w <- 12 / sqrt(z[i] * cosh(opt$maximum))
      I <- stats::integrate(function(u) exp(g(u) - opt$objective), max(0, opt$maximum - w), opt$maximum + w,
                            rel.tol = 1e-12, subdivisions = 1000, stop.on.error = FALSE)$value
      opt$objective + log(I)
    }
  }
  out
}

dvg_log <- function(x, theta, sigma, nu) {
  lam <- 1 / nu
  cc <- sqrt(theta^2 + 2 * sigma^2 / nu)
  ax <- pmax(abs(x), 1e-300)
  z <- cc * ax / sigma^2
  log(2) + theta * x / sigma^2 - 0.5 * log(2 * pi) - log(sigma) - lam * log(nu) - lgamma(lam) +
    (lam - 0.5) * (log(ax) - log(cc)) + log_besselKx(z, lam - 0.5) - z
}

# Densité dont les valeurs non finies (sous-dépassement extrême) sont mises à 0, pour l'intégration numérique
safe_density <- function(logdens) {
  function(x) {
    v <- exp(logdens(x))
    v[!is.finite(v)] <- 0
    v
  }
}

dmeixner_log <- function(x, a, b, d, m) {
  y <- (x - m) / a
  2 * d * log(2 * cos(b / 2)) - log(2 * a * pi) - lgamma(2 * d) + b * y +
    2 * Re(gsl::lngamma_complex(rep_len(d, length(y)), y))
}

# Ajustement VG. Quand la vraisemblance croît sans fin vers nu = 0, la loi ajustée est la limite
# brownienne : la loi normale du maximum de vraisemblance, notée nu = 0
fit_vg <- function(d) {
  p <- suppressWarnings(levyfit::vg_mle(d))
  if (!anyNA(p)) return(p)
  if (!grepl("Brownian limit", attr(p, "message"))) stop(attr(p, "message"), call. = FALSE)
  w <- levyfit::wiener_mle(d)
  c(theta = w[["mu"]], sigma = sqrt(w[["sigma2"]]), nu = 0)
}

vg_limit_note <- "VG fitted at its Brownian limit (nu tends to 0): normal distribution"

pit_normal <- function(y, m, s) {
  y <- sort(y)
  list(F = pmax(stats::pnorm(y, m, s), 1e-300), S = pmax(stats::pnorm(y, m, s, lower.tail = FALSE), 1e-300))
}

# Lois sur un pas 1 : la famille est stable par changement de pas, le pas des accroissements est sans importance
# Ajustement sur les accroissements eux-mêmes (increments_data) : aucune perte de précision
families <- list(
  VG = list(
    fit = function(y) fit_vg(levyfit::increments_data(y)),
    pit = function(y, p) {
      if (p[["nu"]] == 0) return(pit_normal(y, p[["theta"]], p[["sigma"]]))
      pit_numeric(y, safe_density(function(x) dvg_log(x, p[["theta"]], p[["sigma"]], p[["nu"]])), split = 0)
    },
    sim = function(n, p) {
      if (p[["nu"]] == 0) return(p[["theta"]] + p[["sigma"]] * stats::rnorm(n))
      g <- stats::rgamma(n, 1 / p[["nu"]], scale = p[["nu"]])
      p[["theta"]] * g + p[["sigma"]] * sqrt(g) * stats::rnorm(n)
    }),
  Meixner = list(
    fit = function(y) {
      p <- suppressWarnings(levyfit::meixner_mle(levyfit::increments_data(y)))
      if (anyNA(p)) stop(attr(p, "message"), call. = FALSE)
      p
    },
    pit = function(y, p) pit_numeric(y, safe_density(function(x) dmeixner_log(x, p[["a"]], p[["b"]], p[["d"]], p[["m"]]))),
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

# Résultat d'un test par bootstrap : statistique, p-valeur et valeur critique (quantile 1 - alpha)
boot_result <- function(T_obs, T_boot, alpha) {
  list(statistic = T_obs, p.value = boot_pvalue(T_obs, T_boot),
       cv = unname(stats::quantile(T_boot, 1 - alpha, na.rm = TRUE)))
}

# Bootstrap commun : les tests d'un même modèle utilisent les mêmes B échantillons simulés (une
# simulation et une estimation par échantillon), ce qui ne change pas la loi de la p-valeur de chaque test.
# T_boot : matrice avec une ligne par test et une colonne par échantillon
boot_stats <- function(stats, B, one) {
  matrix(vapply(seq_len(B), function(i) tryCatch(one(), error = function(e) rep(NA_real_, length(stats))),
                numeric(length(stats))), nrow = length(stats))
}

# T_obs calculée test par test : une statistique impossible à calculer ne fait échouer que son test
obs_stats <- function(stats, f) lapply(stats, function(s) tryCatch(f(s), error = function(e) e))

# Tests de gofgamma (HME1, BE, KS, CvM, AD, Watson) : chaque test refait ligne pour ligne
# gofgamma::test.HME1, test.BE, test.KS, test.CM, test.AD et test.WA (avec gofgamma::crit.values), avec son
# propre bootstrap : même statistique, même valeur critique et même décision. Les statistiques simulées
# sont gardées pour donner aussi la p-valeur.
gamma_stats <- list(HME1 = function(y, k) gofgamma::HME1(y, k, 1), BE = function(y, k) gofgamma::BE(y, k, 1),
                    KS = function(y, k) gofgamma::KS(y, k), CvM = function(y, k) gofgamma::CM(y, k),
                    AD = function(y, k) gofgamma::AD(y, k), Watson = function(y, k) gofgamma::WA(y, k))

gofgamma_test <- function(dx, test, B, alpha) {
  f <- gamma_stats[[test]]
  par <- gofgamma::gamma_est(dx)
  n <- length(dx)
  y <- rep(0, B)
  for (j in seq_len(B)) {
    BS <- stats::rgamma(n, par[1], 1)
    k <- gofgamma::gamma_est(BS)[1]
    y[j] <- f(BS / (mean(BS) / k), k)
  }
  Tn <- sort(y)
  i <- floor((1 - alpha) * B)
  cv <- if (i >= 1 && i < B) Tn[i] + (1 - alpha) * (Tn[i + 1] - Tn[i]) else NA_real_
  T_obs <- f(dx / par[2], par[1])
  list(statistic = T_obs, cv = cv, Decision = T_obs > cv, p.value = mean(y >= T_obs))
}

# Tests de la loi normale (Wiener) : statistiques des packages nortest, stats, tseries et fBasics,
# p-valeur par bootstrap paramétrique commun. Les statistiques ne changent pas par changement de position
# et d'échelle : un échantillon normal simulé avec la moyenne et l'écart-type estimés, puis réestimé, a la
# même loi de statistique qu'un échantillon N(0, 1). SW rejette pour les petites valeurs de W : on teste -W.
normal_stats <- list(
  AD = function(x) unname(nortest::ad.test(x)$statistic),
  CvM = function(x) unname(suppressWarnings(nortest::cvm.test(x))$statistic),
  SW = function(x) -unname(stats::shapiro.test(x)$statistic),
  Lilliefors = function(x) unname(nortest::lillie.test(x)$statistic),
  JB = function(x) unname(tseries::jarque.bera.test(x)$statistic),
  DAgostinoPearson = function(x) unname(fBasics::dagoTest(x)@test$statistic[1]))

boot_normal <- function(dx, tests, B, alpha) {
  st <- normal_stats[tests]
  m <- mean(dx)
  s <- stats::sd(dx)
  T_obs <- obs_stats(st, function(f) f(dx))
  T_boot <- boot_stats(st, B, function() {
    x <- stats::rnorm(length(dx), m, s)
    vapply(st, function(f) tryCatch(f(x), error = function(e) NA_real_), 0)
  })
  stats::setNames(lapply(seq_along(st), function(j) {
    if (inherits(T_obs[[j]], "error")) return(T_obs[[j]])
    r <- boot_result(T_obs[[j]], T_boot[j, ], alpha)
    if (tests[j] == "SW") {
      r$statistic <- -r$statistic
      r$cv <- -r$cv
    }
    r
  }), tests)
}

# Tests KS, CvM, AD et Watson de la famille `model` (VG ou Meixner)
boot_pit <- function(dx, model, tests, B, alpha) {
  fam <- families[[model]]
  st <- pit_stats[tests]
  par <- fam$fit(dx)
  p0 <- fam$pit(dx, par)
  T_obs <- obs_stats(st, function(f) f(p0))
  T_boot <- boot_stats(st, B, function() {
    y <- fam$sim(length(dx), par)
    pp <- fam$pit(y, fam$fit(y))
    vapply(st, function(f) f(pp), 0)
  })
  note <- if (model == "VG" && par[["nu"]] == 0) vg_limit_note else NULL
  stats::setNames(lapply(seq_along(st), function(j) {
    if (inherits(T_obs[[j]], "error")) return(T_obs[[j]])
    r <- boot_result(T_obs[[j]], T_boot[j, ], alpha)
    r$note <- note
    r
  }), tests)
}

# Tests de gofIG (ABEV1, ABEV2, HK1, HK2, KS, CvM, AD) : statistiques de gofIG et même bootstrap
# paramétrique que gofIG (échantillons IG(1, phi), phi = lambda / mu du maximum de vraisemblance ;
# p-valeur : proportion des statistiques simulées strictement plus grandes)
# Tests VGE_ratio et VGE_gamma : statistiques de goft::ig_test (|T| pour le test bilatéral du rapport de
# variances), avec le même bootstrap. Toutes ces statistiques ne changent pas par changement d'échelle.
ig_stats <- list(ABEV1 = function(x) gofIG::ABEV1(x, 10, "MLE"), ABEV2 = function(x) gofIG::ABEV2(x, 10, "MLE"),
                 HK1 = function(x) gofIG::HK1(x, 0), HK2 = function(x) gofIG::HK2(x),
                 KS = function(x) gofIG::KS(x), CvM = function(x) gofIG::CM(x), AD = function(x) gofIG::AD(x),
                 VGE_ratio = function(x) abs(unname(goft::ig_test(x, method = "ratio")$statistic)),
                 VGE_gamma = function(x) unname(goft::ig_test(x, method = "transf")[[2]]$statistic))

# Loi IG(m, lambda) par la méthode de Michael, Schucany et Haas (1976)
rig <- function(n, m, lambda) {
  nu <- stats::rnorm(n)^2
  x1 <- m + m^2 * nu / (2 * lambda) - m / (2 * lambda) * sqrt(4 * m * lambda * nu + m^2 * nu^2)
  ifelse(stats::runif(n) <= m / (m + x1), x1, m^2 / x1)
}

boot_ig <- function(dx, tests, B) {
  st <- ig_stats[tests]
  mu <- mean(dx)
  phi <- (1 / mean(1 / dx - 1 / mu)) / mu
  T_obs <- obs_stats(st, function(f) f(dx))
  T_boot <- boot_stats(st, B, function() {
    x <- rig(length(dx), 1, phi)
    vapply(st, function(f) f(x), 0)
  })
  stats::setNames(lapply(seq_along(st), function(j) {
    if (inherits(T_obs[[j]], "error")) return(T_obs[[j]])
    ok <- !is.na(T_boot[j, ])
    if (!any(ok)) return(simpleError("all the bootstrap replications failed"))
    list(statistic = T_obs[[j]], p.value = mean(T_boot[j, ok] > T_obs[[j]]))
  }), tests)
}

# Tous les tests du modèle qui partagent un bootstrap, calculés ensemble : liste nommée par test
shared_tests <- function(dx, model, test, alpha, B) {
  groupe <- intersect(test, switch(model, HWP = names(normal_stats), InvGaus = names(ig_stats),
                                   VG = names(pit_stats), Meixner = names(pit_stats), character(0)))
  if (length(groupe) == 0) return(list())
  r <- tryCatch(switch(model, HWP = boot_normal(dx, groupe, B, alpha), InvGaus = boot_ig(dx, groupe, B),
                       boot_pit(dx, model, groupe, B, alpha)), error = function(e) e)
  if (inherits(r, "error")) r <- stats::setNames(rep(list(r), length(groupe)), groupe)
  r
}

# Test VGE de la loi gamma : statistique de goft::gamma_test, p-valeur par bootstrap paramétrique
gof_gamma_vge <- function(dx, B, alpha) {
  fit <- goft::gamma_fit(dx)
  shape <- fit["shape", 1]
  scale <- fit["scale", 1]
  stat <- function(y) abs(as.numeric(suppressWarnings(goft::gamma_test(y))$statistic))
  T_obs <- stat(dx)
  T_boot <- vapply(seq_len(B), function(i) tryCatch(stat(stats::rgamma(length(dx), shape = shape, scale = scale)),
                                                   error = function(e) NA_real_), 0)
  boot_result(T_obs, T_boot, alpha)
}
