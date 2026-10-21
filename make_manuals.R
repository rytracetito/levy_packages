# Manuels de référence (format CRAN) des packages, rangés dans le dossier manuels/
# À exécuter dans RStudio, ligne par ligne ou en entier.

packages <- c("levysim", "levyfit", "levyBridge")

racine <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  getwd()
}
sortie <- file.path(racine, "manuels")
dir.create(sortie, showWarnings = FALSE)

# R CMD Rd2pdf échoue si le chemin contient un accent ("Ma Thèse") :
# chaque package est copié dans le dossier temporaire de R, sans accent.
manuel <- function(pkg) {
  roxygen2::roxygenise(file.path(racine, pkg))

  copie <- file.path(tempdir(), pkg)
  unlink(copie, recursive = TRUE)
  dir.create(copie)
  file.copy(file.path(racine, pkg, c("DESCRIPTION", "NAMESPACE")), copie)
  file.copy(file.path(racine, pkg, "man"), copie, recursive = TRUE)

  # Champs Author et Maintainer, ajoutés d'habitude par R CMD build
  d <- read.dcf(file.path(copie, "DESCRIPTION"))
  auteurs <- eval(parse(text = d[, "Authors@R"]))
  cre <- auteurs[vapply(auteurs, function(p) "cre" %in% p$role, TRUE)]
  d <- cbind(d,
             Author = paste(format(auteurs, include = c("given", "family", "role"),
                                   braces = list(role = c("[", "]"))), collapse = ", "),
             Maintainer = format(cre, include = c("given", "family", "email")))
  write.dcf(d, file.path(copie, "DESCRIPTION"))

  pdf_tmp <- file.path(tempdir(), paste0(pkg, ".pdf"))
  system2(file.path(R.home("bin"), "R"),
          c("CMD", "Rd2pdf", shQuote(copie), paste0("--output=", shQuote(pdf_tmp)), "--force", "--no-preview"))
  file.copy(pdf_tmp, file.path(sortie, paste0(pkg, ".pdf")), overwrite = TRUE)
}

for (pkg in packages) manuel(pkg)

list.files(sortie)
