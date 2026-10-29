#' Goodness-of-fit tests for an i.i.d. sample of increments
#'
#' Applies goodness-of-fit tests to a sample of increments of a Levy process
#' observed with a constant time step: the increments are then independent
#' and identically distributed. This is the last step of [gof.levy()], also
#' useful alone when the increments are already available (regular
#' inspections, or increments reconstructed by the package `levyBridge`).
#'
#' @param dx Numeric vector of at least 3 increments, all over the same time
#'   step.
#' @param model `"HGP"` (gamma process), `"HWP"` (Wiener process),
#'   `"InvGaus"` (inverse Gaussian process), `"VG"` (variance-gamma process)
#'   or `"Meixner"` (Meixner process). Upper or lower case does not matter.
#' @param test Names of the tests to apply; `NULL` (default) applies all the
#'   tests of the model, see [tests_disponibles()].
#' @param alpha Level of the tests, strictly between 0 and 1.
#' @param B Number of bootstrap replications, for the tests whose p-value or
#'   critical value is computed by parametric bootstrap.
#'
#' @return A data frame with one row per test and the columns:
#'   * `test`: name of the test;
#'   * `statistic`: value of the test statistic;
#'   * `p.value`: p-value; the test rejects H0 when `p.value <= alpha`;
#'   * `critical.value`: for the tests whose p-value is computed by
#'     parametric bootstrap, the quantile of order `1 - alpha` of the
#'     bootstrap statistics (H0 is rejected when the statistic exceeds it);
#'     `NA` when it is not available (tests with an exact or asymptotic
#'     p-value, and tests of `gofIG`, which give only the p-value);
#'   * `decision`: `"reject H0"`, `"do not reject H0"`, or `"not computed"`
#'     when the test failed;
#'   * `message`: reason of a failure or of a certain rejection, `""`
#'     otherwise.
#'
#'   The column `message` is not printed when it is empty for all the tests;
#'   it stays in the data frame.
#'
#' @details
#' The null hypothesis is: the increments are i.i.d. with a distribution of the
#' family of the model (parameters unknown). Under the gamma and inverse
#' Gaussian models the increments are positive with probability one: a
#' negative or zero increment is impossible under H0, and all the tests then
#' reject H0 with a p-value 0, without computation.
#'
#' A test that fails (for instance a sample too small for the estimation)
#' does not stop the others: its row has the decision `"not computed"` and
#' the error message.
#'
#' The tests are described in [tests_disponibles()]. The bootstrap tests are
#' random: use [set.seed()] to reproduce a result.
#'
#' @inherit tests_disponibles references
#'
#' @seealso [gof.levy()] for a trajectory or several trajectories,
#'   [tests_disponibles()] for the list of the tests.
#'
#' @examples
#' set.seed(1)
#'
#' # 1. Gamma increments (time step 2, a = 2, b = 3): all the tests of the model
#' dx <- rgamma(40, shape = 2 * 2, rate = 3)
#' test_adequation(dx, "HGP", B = 100)
#'
#' # 2. Some tests only, at level 10 %
#' test_adequation(dx, "HGP", test = c("KS", "AD"), alpha = 0.1, B = 100)
#'
#' # 3. Normal increments: tests of the Wiener model
#' test_adequation(rnorm(50, mean = 1, sd = 2), "HWP")
#'
#' # "KS" is accepted for the Wiener model: it is the Lilliefors test
#' test_adequation(rnorm(50, mean = 1, sd = 2), "HWP", test = c("KS", "AD"))
#'
#' # 4. The same normal increments are negative: the gamma model is rejected for sure
#' test_adequation(rnorm(50), "HGP", test = c("KS", "AD"))
#'
#' # 5. Inverse Gaussian increments
#' dx_ig <- 1 / rgamma(40, 3) # not inverse Gaussian
#' test_adequation(dx_ig, "invgaus", test = c("KS", "AD", "VGE_ratio"), B = 100)
#'
#' \donttest{
#' # 6. Variance-gamma model (slower: estimation and numerical integration
#' #    in each bootstrap replication)
#' g <- rgamma(100, shape = 1, scale = 1)
#' dx_vg <- 0.5 * g + 1.5 * sqrt(g) * rnorm(100)
#' test_adequation(dx_vg, "VG", test = "AD", B = 30)
#'
#' # 7. Meixner model (a = 1, b = 1, d = 0.5, m = 0.2)
#' set.seed(1)
#' gen <- Runuran::pinvd.new(Runuran::udmeixner(alpha = 1, beta = 1, delta = 0.5, mu = 0.2))
#' test_adequation(Runuran::ur(gen, 100), "Meixner", test = c("KS", "AD"), B = 30)
#' }
#' @export
test_adequation <- function(dx, model, test = NULL, alpha = 0.05, B = 200) {
  model <- model_name(model)
  test <- test_names(test, model)
  check_alpha_B(alpha, B)
  if (!is.numeric(dx) || length(dx) < 3 || !all(is.finite(dx))) {
    stop("dx must contain at least 3 finite numeric increments", call. = FALSE)
  }

  if (model %in% c("HGP", "InvGaus") && any(dx <= 0)) {
    msg <- paste0(sum(dx <= 0), " negative or zero increment(s): impossible under H0 (increasing process)")
    return(tests_table(do.call(rbind, lapply(test, function(nom) result_row(nom, NA_real_, 0, NA_real_, TRUE, msg)))))
  }
  tests_table(do.call(rbind, lapply(test, function(nom) one_test(dx, model, nom, alpha, B))))
}

tests_table <- function(d) {
  rownames(d) <- NULL
  structure(d, class = c("levygof_tests", "data.frame"))
}

#' @export
print.levygof_tests <- function(x, ...) {
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  if (all(d$message == "")) d$message <- NULL
  print(d, ...)
  invisible(x)
}

result_row <- function(test, stat, p, cv, reject, message) {
  data.frame(test = test, statistic = stat, p.value = p, critical.value = cv,
             decision = if (is.na(reject)) "not computed" else if (reject) "reject H0" else "do not reject H0",
             message = message, stringsAsFactors = FALSE)
}

call_test <- function(dx, model, test, alpha, B) {
  switch(paste(model, test),
    "HGP HME1"             = gof_gamma_boot(dx, gofgamma::HME1, B, alpha, tuning = 1),
    "HGP BE"               = gof_gamma_boot(dx, gofgamma::BE, B, alpha, tuning = 1),
    "HGP VGE"              = gof_gamma_vge(dx, B, alpha),
    "HGP KS"               = gof_gamma_boot(dx, gofgamma::KS, B, alpha),
    "HGP CvM"              = gof_gamma_boot(dx, gofgamma::CM, B, alpha),
    "HGP AD"               = gof_gamma_boot(dx, gofgamma::AD, B, alpha),
    "HGP Watson"           = gof_gamma_boot(dx, gofgamma::WA, B, alpha),
    "HWP AD"               = nortest::ad.test(dx),
    "HWP CvM"              = suppressWarnings(nortest::cvm.test(dx)),
    "HWP SW"               = stats::shapiro.test(dx),
    "HWP Lilliefors"       = nortest::lillie.test(dx),
    "HWP JB"               = tseries::jarque.bera.test(dx),
    "HWP DAgostinoPearson" = gof_dagostino_pearson(dx),
    "InvGaus ABEV1"        = gofIG::test.ABEV1(dx, a = 10, meth = "MLE", B = B),
    "InvGaus ABEV2"        = gofIG::test.ABEV2(dx, a = 10, meth = "MLE", B = B),
    "InvGaus HK1"          = gofIG::test.HK1(dx, a = 0, B = B),
    "InvGaus HK2"          = gofIG::test.HK2(dx, B = B),
    "InvGaus KS"           = gofIG::test.KS(dx, B = B),
    "InvGaus CvM"          = gofIG::test.CM(dx, B = B),
    "InvGaus AD"           = gofIG::test.AD(dx, B = B),
    "InvGaus VGE_ratio"    = goft::ig_test(dx, method = "ratio"),
    "InvGaus VGE_gamma"    = goft::ig_test(dx, method = "transf")[[2]],
    gof_pit(dx, model, test, B, alpha)
  )
}

one_test <- function(dx, model, test, alpha, B) {
  r <- tryCatch(call_test(dx, model, test, alpha, B), error = function(e) e)
  if (inherits(r, "error")) return(result_row(test, NA_real_, NA_real_, NA_real_, NA, conditionMessage(r)))
  stat <- as.numeric(if (!is.null(r$T.value)) r$T.value else r$statistic)[1]
  p <- if (is.null(r$p.value)) NA_real_ else as.numeric(r$p.value)[1]
  cv <- if (is.null(r$cv)) NA_real_ else as.numeric(r$cv)[1]
  reject <- if (!is.na(p)) p <= alpha else if (is.logical(r$Decision) && length(r$Decision) == 1) r$Decision else NA
  result_row(test, stat, p, cv, reject, "")
}
