library(levyfit)
library(levysim)

# Vraies valeurs, plans d'inspection et nombre de répétitions
vraies <- list(gamma   = c(a = 1.5, b = 3),
               wiener  = c(mu = 0.5, sigma2 = 2.25),
               ig      = c(mu = 0.5, eta = 2),
               vg      = c(theta = 0.5, sigma = 1.5, nu = 2),
               meixner = c(a = 1, b = 0.3, d = 1.5, m = 0.2))
noms <- c(gamma = "Gamma", wiener = "Wiener", ig = "Inverse gaussien",
          vg = "Variance-gamma", meixner = "Meixner")

plans <- c(regulier = "regular", croissant = "increasing", decroissant = "decreasing",
           mixte_bouts = "dec-inc", mixte_milieu = "inc-dec")
libelles <- c(regulier = "Régulier", croissant = "Croissant", decroissant = "Décroissant",
              mixte_bouts = "Mixte (bouts)", mixte_milieu = "Mixte (milieu)")

R <- 200
R_meixner <- 100

# Simulation et estimation
simuler <- function(modele, plan, p) {
  switch(modele,
    gamma   = sim_gamma(plan, p[["a"]], p[["b"]]),
    wiener  = sim_wiener(plan, p[["mu"]], sqrt(p[["sigma2"]])),
    ig      = sim_ig(plan, p[["mu"]], p[["eta"]]),
    vg      = sim_vg(plan, p[["theta"]], p[["sigma"]], p[["nu"]]),
    meixner = sim_meixner(plan, p[["a"]], p[["b"]], p[["d"]], p[["m"]]))
}

estimer <- function(modele, x) {
  f <- switch(modele, gamma = gamma_mle, wiener = wiener_mle, ig = ig_mle,
              vg = vg_mle, meixner = meixner_mle)
  suppressWarnings(f(x, ci = TRUE))
}

plan_de <- function(modele, config, nom_plan) {
  mx <- modele == "meixner"
  if (config == "une") {
    inspection_times(plans[[nom_plan]], n = if (mx) 499 else 99,
                     tmax = if (mx) 500 else 100)$time
  } else {
    inspection_times(plans[[nom_plan]], n = if (mx) 19 else 9, n_traj = 20,
                     tmax = if (mx) 20 else 100)
  }
}

une_cellule <- function(modele, config, nom_plan) {
  plan <- plan_de(modele, config, nom_plan)
  rep <- if (modele == "meixner") R_meixner else R
  do.call(rbind, lapply(seq_len(rep), function(r) {
    e <- estimer(modele, simuler(modele, plan, vraies[[modele]]))
    data.frame(modele, config, plan = nom_plan, r, e)
  }))
}

# Calcul en parallèle des 50 situations, ou relecture du fichier de résultats
fichier <- "test.mle_resultats.rds"
if (file.exists(fichier)) {
  res <- readRDS(fichier)
} else {
  cellules <- expand.grid(modele = names(vraies), config = c("une", "plusieurs"),
                          plan = names(plans), stringsAsFactors = FALSE)
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  parallel::clusterSetRNGStream(cl, 2026)
  parallel::clusterEvalQ(cl, {library(levyfit); library(levysim)})
  parallel::clusterExport(cl, c("vraies", "plans", "R", "R_meixner", "simuler",
                                "estimer", "plan_de", "une_cellule"))
  res <- do.call(rbind, parallel::parLapply(cl, split(cellules, seq_len(nrow(cellules))),
                                            function(k) une_cellule(k$modele, k$config, k$plan)))
  parallel::stopCluster(cl)
  saveRDS(res, fichier)
}

# Résumé des répétitions
resumer <- function(res) {
  garder <- mapply(function(m, p) p %in% names(vraies[[m]]), res$modele, res$parameter)
  res <- res[garder, ]
  cles <- unique(res[c("modele", "config", "plan", "parameter")])
  do.call(rbind, lapply(seq_len(nrow(cles)), function(i) {
    k <- cles[i, ]
    d <- res[res$modele == k$modele & res$config == k$config & res$plan == k$plan &
               res$parameter == k$parameter, ]
    v <- vraies[[k$modele]][[k$parameter]]
    ok <- !is.na(d$estimate)
    e <- d$estimate[ok]
    couv <- if (all(is.na(d$lower[ok]))) NA else
      100 * mean(d$lower[ok] <= v & v <= d$upper[ok], na.rm = TRUE)
    data.frame(k, vraie = v, moyenne = mean(e), biais = mean(e) - v,
               rmse = sqrt(mean((e - v)^2)), rmse_rel = 100 * sqrt(mean((e - v)^2)) / abs(v),
               couverture = couv, echecs = 100 * mean(!ok))
  }))
}
synthese <- resumer(res)

# Tableau d'un modèle pour une configuration
tableau <- function(modele, config) {
  d <- synthese[synthese$modele == modele & synthese$config == config, ]
  d <- d[order(match(d$parameter, names(vraies[[modele]])), match(d$plan, names(plans))), ]
  f <- function(v, k = 3) ifelse(is.na(v), "", formatC(v, digits = k, format = "fg", flag = "#"))
  t <- data.frame(Paramètre = ifelse(duplicated(d$parameter), "",
                                     paste0(d$parameter, " (", d$vraie, ")")),
                  Plan = libelles[d$plan], Moyenne = f(d$moyenne), Biais = f(d$biais),
                  RMSE = f(d$rmse),
                  `RMSE rel. %` = formatC(d$rmse_rel, digits = 1, format = "f"),
                  `Couv. %` = ifelse(is.na(d$couverture), "",
                                     formatC(d$couverture, digits = 0, format = "f")),
                  `Échecs %` = formatC(d$echecs, digits = 0, format = "f"), check.names = FALSE)
  tab <- knitr::kable(t, format = "latex", booktabs = TRUE, row.names = FALSE,
                      align = "llrrrrrr", linesep = c(rep("", length(plans) - 1), "\\addlinespace"))
  knitr::asis_output(paste0("\\begin{small}\n", tab, "\n\\end{small}"))
}
