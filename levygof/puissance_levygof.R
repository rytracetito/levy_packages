# Taux de rejet de gof.levy() : chaque processus simulé est testé contre les cinq modèles.
# Le résultat est enregistré dans puissance_levygof.rds, lu par exemple_levygof.Rmd.

library(levygof)
library(DegradSim)

dossier <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  getwd()
}
fichier <- file.path(dossier, "puissance_levygof.rds")

# N : nombre de répétitions pour les modèles rapides (HGP, HWP, InvGaus) et pour VG et Meixner
N_rapide <- 100
N_lent <- 50

plan <- inspection_times("dec-inc", n = 60, tmax = 120, digits = 1)

simuler <- function(processus) {
  switch(processus,
         gamma   = sim_gamma(plan$time, a = 1, b = 1, n_traj = 5),
         Wiener  = sim_wiener(plan$time, mu = 1, sigma = 1, n_traj = 5),
         IG      = sim_ig(plan$time, mu = 1, eta = 1, n_traj = 5),
         VG      = sim_vg(plan$time, theta = 1, sigma = 1, nu = 2, n_traj = 5),
         Meixner = sim_meixner(plan$time, a = 1, b = 1, d = 0.5, m = 0.2, n_traj = 5))
}

processus <- c("gamma", "Wiener", "IG", "VG", "Meixner")
modeles <- c("HGP", "HWP", "InvGaus", "VG", "Meixner")
test <- c(HGP = "VGE", HWP = "SW", InvGaus = "VGE_ratio", VG = "AD", Meixner = "AD")
B <- c(HGP = 100, HWP = 200, InvGaus = 200, VG = 30, Meixner = 30)

rejet <- matrix(NA, length(processus), length(modeles), dimnames = list(processus, modeles))
non_estime <- rejet

set.seed(2026)
debut <- Sys.time()
for (m in modeles) {
  N <- if (m %in% c("VG", "Meixner")) N_lent else N_rapide
  for (p in processus) {
    decision <- character(N)
    for (k in 1:N) {
      r <- tryCatch(gof.levy(simuler(p), m, test = test[[m]], B = B[[m]])$results,
                    error = function(e) NULL)
      decision[k] <- if (is.null(r)) "non estimable" else r$decision[r$group == "all"]
    }
    rejet[p, m] <- mean(decision == "reject H0")
    non_estime[p, m] <- mean(decision %in% c("non estimable", "not computed"))
    cat(p, "testé comme", m, ": rejet", rejet[p, m], "| non estimable", non_estime[p, m],
        "|", round(difftime(Sys.time(), debut, units = "mins"), 1), "min\n")
  }
}

saveRDS(list(rejet = rejet, non_estime = non_estime, N_rapide = N_rapide, N_lent = N_lent,
             test = test, B = B), fichier)
