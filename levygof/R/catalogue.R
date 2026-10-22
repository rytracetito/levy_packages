catalogue <- list(
  HGP     = c("HME1", "BE", "VGE", "KS", "CvM", "AD", "Watson"),
  HWP     = c("AD", "CvM", "SW", "Lilliefors", "JB", "DAgostinoPearson"),
  InvGaus = c("ABEV1", "ABEV2", "HK1", "HK2", "KS", "CvM", "AD", "VGE_ratio", "VGE_gamma"),
  VG      = c("KS", "CvM", "AD", "Watson"),
  Meixner = c("KS", "CvM", "AD", "Watson")
)

#' Available goodness-of-fit tests
#'
#' Lists the tests that [gof.levy()] and [test_adequation()] can apply to
#' each model.
#'
#' @param model `NULL` (default) for all the models, or one model: `"HGP"`
#'   (gamma process), `"HWP"` (Wiener process), `"InvGaus"` (inverse Gaussian
#'   process), `"VG"` (variance-gamma process) or `"Meixner"` (Meixner
#'   process). Upper or lower case does not matter.
#'
#' @return With `model = NULL`, a named list with one vector of test names
#'   per model; otherwise the vector of test names of this model.
#'
#' @section The tests:
#' All the tests are tests of a composite hypothesis: the increments
#' \eqn{\Delta x_1, \dots, \Delta x_n} are i.i.d. with a law of the family of
#' the model, whose parameters are unknown and estimated on the sample.
#'
#' **Gamma process (`"HGP"`)**, increments \eqn{\Gamma(a \Delta, b)}:
#' * `"HME1"`: Henze, Meintanis and Ebner (2012), Laplace transform
#'   (`gofgamma::test.HME1`, tuning parameter 1);
#' * `"BE"`: Betsch and Ebner (2019), characterization by a fixed point
#'   (`gofgamma::test.BE`, tuning parameter 1);
#' * `"VGE"`: Villasenor and Gonzalez-Estrada (2015), ratio of two
#'   variance estimators (`goft::gamma_test`), with a p-value by
#'   parametric bootstrap;
#' * `"KS"`, `"CvM"`, `"AD"`, `"Watson"`: Kolmogorov-Smirnov,
#'   Cramer-von Mises, Anderson-Darling and Watson (`gofgamma`).
#'
#' The tests of `gofgamma` give the decision by comparison with a bootstrap
#' critical value, without p-value.
#'
#' **Wiener process (`"HWP"`)**, increments
#' \eqn{N(\mu \Delta, \sigma^2 \Delta)}: tests of normality with estimated
#' mean and variance.
#' * `"AD"`, `"CvM"`, `"Lilliefors"`: Anderson-Darling, Cramer-von Mises and
#'   Lilliefors (Kolmogorov-Smirnov) tests (`nortest`);
#' * `"SW"`: Shapiro-Wilk (`stats::shapiro.test`);
#' * `"JB"`: Jarque-Bera (`tseries::jarque.bera.test`);
#' * `"DAgostinoPearson"`: D'Agostino-Pearson omnibus test \eqn{K^2},
#'   combining skewness and kurtosis (`fBasics::dagoTest`).
#'
#' **Inverse Gaussian process (`"InvGaus"`)**, increments
#' \eqn{IG(\mu \Delta, \eta \Delta^2)}:
#' * `"ABEV1"`, `"ABEV2"`: Allison, Betsch, Ebner and Visagie (2022)
#'   (`gofIG`, tuning parameter 10, maximum likelihood);
#' * `"HK1"`, `"HK2"`: Henze and Klar (2002) (`gofIG`);
#' * `"KS"`, `"CvM"`, `"AD"`: Kolmogorov-Smirnov, Cramer-von Mises and
#'   Anderson-Darling (`gofIG`);
#' * `"VGE_ratio"`, `"VGE_gamma"`: Villasenor, Gonzalez-Estrada and Ochoa
#'   (2019), ratio of variances and transformation to gamma variables
#'   (`goft::ig_test`).
#'
#' **Variance-gamma (`"VG"`) and Meixner (`"Meixner"`) processes**:
#' `"KS"`, `"CvM"`, `"AD"`, `"Watson"`. The parameters are estimated by
#' maximum likelihood ([levyfit::vg_mle()], [levyfit::meixner_mle()]), the
#' distribution function is computed by numerical integration of the density,
#' and the p-value is obtained by parametric bootstrap.
#'
#' @seealso [gof.levy()], [test_adequation()]
#'
#' @examples
#' # All the models and their tests
#' tests_disponibles()
#'
#' # The tests of one model (upper or lower case)
#' tests_disponibles("HGP")
#' tests_disponibles("meixner")
#' @export
tests_disponibles <- function(model = NULL) {
  if (is.null(model)) return(catalogue)
  catalogue[[model_name(model)]]
}

model_name <- function(model) {
  if (!is.character(model) || length(model) != 1 || is.na(model)) stop("model must be a character string", call. = FALSE)
  k <- match(tolower(model), tolower(names(catalogue)))
  if (is.na(k)) stop("unknown model '", model, "'; available models: ", paste(names(catalogue), collapse = ", "), call. = FALSE)
  names(catalogue)[k]
}

test_names <- function(test, model) {
  dispo <- catalogue[[model]]
  if (is.null(test)) return(dispo)
  if (!is.character(test) || length(test) == 0 || anyNA(test)) stop("test must be a vector of test names", call. = FALSE)
  k <- match(tolower(test), tolower(dispo))
  if (anyNA(k)) {
    stop("unknown test for the model ", model, ": ", paste(test[is.na(k)], collapse = ", "),
         "; available tests: ", paste(dispo, collapse = ", "), call. = FALSE)
  }
  unique(dispo[k])
}

check_alpha_B <- function(alpha, B) {
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number strictly between 0 and 1", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1 || is.na(B) || B < 1 || B != round(B)) {
    stop("B must be an integer greater than or equal to 1", call. = FALSE)
  }
}
