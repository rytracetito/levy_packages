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
#'   (statistic `gofgamma::HME1`, tuning parameter 1);
#' * `"BE"`: Betsch and Ebner (2019), characterization by a fixed point
#'   (statistic `gofgamma::BE`, tuning parameter 1);
#' * `"VGE"`: Villasenor and Gonzalez-Estrada (2015), ratio of two
#'   variance estimators (`goft::gamma_test`), with a p-value by
#'   parametric bootstrap;
#' * `"KS"`, `"CvM"`, `"AD"`, `"Watson"`: Kolmogorov-Smirnov,
#'   Cramer-von Mises, Anderson-Darling and Watson (`gofgamma`).
#'
#' For the tests HME1, BE, KS, CvM, AD and Watson, the statistic and the
#' parametric bootstrap are those of the package `gofgamma` (the shape is
#' estimated by `gofgamma::gamma_est`, and each bootstrap sample of the
#' gamma law with this shape is estimated again); all the bootstrap
#' statistics are kept, which gives both the p-value and the critical value.
#'
#' **Wiener process (`"HWP"`)**, increments
#' \eqn{\mathcal{N}(\mu \Delta, \sigma^2 \Delta)}: tests of normality with estimated
#' mean and variance.
#' * `"AD"`, `"CvM"`, `"Lilliefors"`: Anderson-Darling, Cramer-von Mises and
#'   Lilliefors tests (`nortest`). The Lilliefors test is the
#'   Kolmogorov-Smirnov test with estimated mean and variance: `"KS"` is also
#'   accepted for it, and the result is labelled `"Lilliefors"`;
#' * `"SW"`: Shapiro-Wilk (`stats::shapiro.test`);
#' * `"JB"`: Jarque-Bera (`tseries::jarque.bera.test`);
#' * `"DAgostinoPearson"`: D'Agostino-Pearson omnibus test \eqn{K^2},
#'   combining skewness and kurtosis (`fBasics::dagoTest`).
#'
#' **Inverse Gaussian process (`"InvGaus"`)**, increments
#' \eqn{IG(\mu \Delta, \eta \Delta^2)}:
#' * `"ABEV1"`, `"ABEV2"`: Allison, Betsch, Ebner and Visagie (2022)
#'   (`gofIG`, tuning parameter 10, maximum likelihood), characterization of
#'   the distribution function (see below);
#' * `"HK1"`, `"HK2"`: Henze and Klar (2002) (`gofIG`), empirical Laplace
#'   transform (see below);
#' * `"KS"`, `"CvM"`, `"AD"`: Kolmogorov-Smirnov, Cramer-von Mises and
#'   Anderson-Darling (`gofIG`);
#' * `"VGE_ratio"`: Villasenor and Gonzalez-Estrada (2015), ratio of two
#'   variance estimators, asymptotic p-value, suited to large samples
#'   (`goft::ig_test`, method `"ratio"`);
#' * `"VGE_gamma"`: Villasenor and Gonzalez-Estrada (2015), transformation of
#'   the data to gamma variables of shape 1/2 and Anderson-Darling test
#'   (`goft::ig_test`, method `"transf"`).
#'
#' **Tests ABEV1, ABEV2, HK1 and HK2.** Let \eqn{X_1, \dots, X_n} be the
#' sample, \eqn{IG(m, \lambda)} the inverse Gaussian law of mean \eqn{m} and
#' shape \eqn{\lambda}, \eqn{\hat m = \bar X} and
#' \eqn{\hat\lambda = 1 / \overline{(1/X - 1/\bar X)}} the maximum likelihood
#' estimators. The four tests use the scaled data
#' \eqn{Y_j = X_j / \hat m} and \eqn{\hat\varphi = \hat\lambda / \hat m}: under
#' H0, \eqn{Y_j} is close to \eqn{IG(1, \varphi)}, and the law of the
#' statistic only depends on \eqn{\varphi}. Each statistic is a weighted
#' \eqn{L^2} distance, H0 is rejected for large values, and the p-value is
#' obtained by parametric bootstrap: `B` samples of \eqn{IG(1, \hat\varphi)}
#' of size \eqn{n}, the statistic being computed again (with new estimates)
#' on each one. `gofIG` gives only the p-value, so the column
#' `critical.value` is `NA`.
#'
#' * `"HK2"` compares the empirical Laplace transform
#'   \eqn{L_n(t) = \frac{1}{n} \sum_j e^{-t Y_j}} with the fitted one
#'   \eqn{L_{\hat\varphi}(t) = \exp(\hat\varphi (1 - \sqrt{1 + 2t/\hat\varphi}))}:
#'   \deqn{HK2 = n \int_0^\infty (L_n(t) - L_{\hat\varphi}(t))^2\, dt.}
#' * `"HK1"` uses the differential equation
#'   \eqn{L(t) + \sqrt{1 + 2t/\varphi}\, L'(t) = 0}, which characterizes the
#'   Laplace transform of \eqn{IG(1, \varphi)}; the statistic measures how
#'   far \eqn{L_n} is from satisfying it:
#'   \deqn{HK1 \propto n \int_0^\infty \left(L_n(t) + \sqrt{1 + 2t/\hat\varphi}\, L_n'(t)\right)^2 e^{-a t}\, dt,}
#'   with \eqn{a = 0}.
#' * `"ABEV1"` and `"ABEV2"` use a characterization of the distribution
#'   function: with the score
#'   \eqn{\tau(y) = -f'(y)/f(y) = \frac{1}{2}(\varphi + 3/y - \varphi/y^2)} of
#'   the density \eqn{f} of \eqn{IG(1, \varphi)}, a positive variable
#'   \eqn{Y} has this law if and only if
#'   \eqn{F(t) = E[\tau(Y) \min(Y, t)]} for all \eqn{t > 0}. With
#'   \eqn{F_n} the empirical distribution function,
#'   \deqn{ABEV = n \int_0^\infty \left(\frac{1}{n} \sum_j \hat\tau(Y_j) \min(Y_j, t) - F_n(t)\right)^2 w(t)\, dt,}
#'   with the weight \eqn{w(t) = e^{-a t}} for ABEV1 and
#'   \eqn{w(t) = e^{-a t^2}} for ABEV2, and \eqn{a = 10}.
#'
#' The four statistics are double sums over the sample, computed with R
#' loops: their cost grows as \eqn{n^2}, times \eqn{B + 1} for the bootstrap.
#' For \eqn{n = 300} and `B = 200`, one test takes about 8 seconds for HK2,
#' 2 minutes for ABEV1 and HK1, and 7 minutes for ABEV2. With many
#' increments (for instance the method `"ExtBridge"` of [gof.levy()]), prefer
#' HK2 or the other tests.
#'
#' **Variance-gamma (`"VG"`) and Meixner (`"Meixner"`) processes**:
#' `"KS"`, `"CvM"`, `"AD"`, `"Watson"` (Kolmogorov-Smirnov, Cramer-von Mises,
#' Anderson-Darling and Watson statistics, Stephens, 1986). The parameters are
#' estimated by maximum likelihood ([levyfit::vg_mle()],
#' [levyfit::meixner_mle()]), the distribution function is computed by
#' numerical integration of the density, and the p-value is obtained by
#' parametric bootstrap (Stute, Gonzalez Manteiga and Presedo Quindimil,
#' 1993). With a small sample, the maximum likelihood estimator of the
#' Meixner law may not exist (maximum on the boundary, towards a shifted
#' gamma or a Gaussian law): the test is then not computed, with this message.
#'
#' @references
#' **Gamma distribution**
#'
#' Henze, N., Meintanis, S. G. and Ebner, B. (2012). Goodness-of-fit tests
#' for the Gamma distribution based on the empirical Laplace transform.
#' *Communications in Statistics - Theory and Methods*, 41(9), 1543-1556.
#' \doi{10.1080/03610926.2010.542851} (tests HME1, KS, CvM, AD and Watson of
#' `gofgamma`)
#'
#' Betsch, S. and Ebner, B. (2019). A new characterization of the Gamma
#' distribution and associated goodness of fit tests. *Metrika*, 82(7),
#' 779-806. \doi{10.1007/s00184-019-00708-7} (test BE)
#'
#' Villasenor, J. A. and Gonzalez-Estrada, E. (2015). A variance ratio test of
#' fit for Gamma distributions. *Statistics and Probability Letters*, 96,
#' 281-286. \doi{10.1016/j.spl.2014.10.001} (test VGE)
#'
#' **Normal distribution (Wiener process)**
#'
#' Shapiro, S. S. and Wilk, M. B. (1965). An analysis of variance test for
#' normality (complete samples). *Biometrika*, 52(3-4), 591-611.
#' \doi{10.1093/biomet/52.3-4.591} (test SW)
#'
#' Royston, P. (1995). Remark AS R94: A remark on Algorithm AS 181: The W test
#' for normality. *Applied Statistics*, 44, 547-551. \doi{10.2307/2986146}
#' (computation of the test SW in `stats::shapiro.test`)
#'
#' Lilliefors, H. W. (1967). On the Kolmogorov-Smirnov test for normality with
#' mean and variance unknown. *Journal of the American Statistical
#' Association*, 62(318), 399-402. \doi{10.1080/01621459.1967.10482916}
#' (test Lilliefors)
#'
#' Dallal, G. E. and Wilkinson, L. (1986). An analytic approximation to the
#' distribution of Lilliefors' test for normality. *The American
#' Statistician*, 40, 294-296. (p-value of the test Lilliefors)
#'
#' Stephens, M. A. (1986). Tests based on EDF statistics. In D'Agostino, R. B.
#' and Stephens, M. A. (eds.), *Goodness-of-Fit Techniques*, Marcel Dekker,
#' New York. (tests AD and CvM with estimated parameters, and the EDF
#' statistics of the VG and Meixner tests)
#'
#' Jarque, C. M. and Bera, A. K. (1987). A test for normality of observations
#' and regression residuals. *International Statistical Review*, 55(2),
#' 163-172. \doi{10.2307/1403192} (test JB)
#'
#' D'Agostino, R. B. and Pearson, E. S. (1973). Tests for departure from
#' normality. Empirical results for the distributions of \eqn{b_2} and
#' \eqn{\sqrt{b_1}}. *Biometrika*, 60(3), 613-622.
#' \doi{10.1093/biomet/60.3.613} (test DAgostinoPearson)
#'
#' **Inverse Gaussian distribution**
#'
#' Allison, J. S., Betsch, S., Ebner, B. and Visagie, I. J. H. (2022). On
#' testing the adequacy of the inverse Gaussian distribution. *Mathematics*,
#' 10(3), 350. \doi{10.3390/math10030350} (tests ABEV1,
#' ABEV2, KS, CvM and AD of `gofIG`)
#'
#' Henze, N. and Klar, B. (2002). Goodness-of-fit tests for the inverse
#' Gaussian distribution based on the empirical Laplace transform. *Annals of
#' the Institute of Statistical Mathematics*, 54(2), 425-444.
#' \doi{10.1023/A:1022442506681} (tests HK1 and HK2)
#'
#' Villasenor, J. A. and Gonzalez-Estrada, E. (2015). Tests of fit for
#' inverse Gaussian distributions. *Statistics and Probability Letters*, 105,
#' 189-194. \doi{10.1016/j.spl.2015.06.017} (tests VGE_ratio and VGE_gamma)
#'
#' **EDF statistics and parametric bootstrap (VG and Meixner)**
#'
#' Anderson, T. W. and Darling, D. A. (1954). A test of goodness of fit.
#' *Journal of the American Statistical Association*, 49, 765-769.
#'
#' Watson, G. S. (1961). Goodness-of-fit tests on a circle. *Biometrika*,
#' 48(1-2), 109-114.
#'
#' Stute, W., Gonzalez Manteiga, W. and Presedo Quindimil, M. (1993).
#' Bootstrap based goodness-of-fit-tests. *Metrika*, 40, 243-256.
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

# Autres noms acceptés : pour la loi normale, le test de Kolmogorov-Smirnov à moyenne et variance estimées
# est le test de Lilliefors
aliases <- list(HWP = c(KS = "Lilliefors"))

test_names <- function(test, model) {
  dispo <- catalogue[[model]]
  if (is.null(test)) return(dispo)
  if (!is.character(test) || length(test) == 0 || anyNA(test)) stop("test must be a vector of test names", call. = FALSE)
  al <- aliases[[model]]
  if (!is.null(al)) {
    i <- match(tolower(test), tolower(names(al)))
    test[!is.na(i)] <- al[i[!is.na(i)]]
  }
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
