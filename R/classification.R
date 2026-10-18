#' Livingston and Lewis' "Effective Test Length".
#'
#' @description  According to Livingston and Lewis (1995), "The effective test length corresponding to a test score is the number of discrete, dichotomously scored, locally independent, equally difficult items required to produce a total score of the same reliability."
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param min The lower-bound (minimum possible value) of the observed-score distribution. Default is 0 (assuming observed scores represent proportions).
#' @param max The upper-bound (maximum possible value) of the observed-score distribution. Default is 1 (assuming observed scores represent proportions).
#' @param reliability The reliability of the observed scores (proportion of observed-score distribution variance shared with true-score distribution).
#' @return An estimate of the effective length of a test, given the stability of the observations it produces.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the reliability of this test was estimated to 0.7. To estimate and
#' # retrieve the effective test length using ETL():
#' ETL(mean = mean(testdata), variance = var(testdata), min = 0, max = 100,
#' reliability = .7)
#' @export
ETL <- function(mean, variance, min = 0, max = 1, reliability) {
  ((mean - min) * (max - mean) - reliability * variance) / (variance * (1 - reliability))
}

#' Model Implied Reliability from Livingston and Lewis' "Effective Test Length".
#'
#' @description  Calculate model-implied reliability given mean, variance, the minimum and maximum possible scores, and the effective test length.
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param min The lower-bound (minimum possible value) of the observed-score distribution. Default is 0 (assuming observed scores represent proportions).
#' @param max The upper-bound (maximum possible value) of the observed-score distribution. Default is 1 (assuming observed scores represent proportions).
#' @param ETL The effective test length as defined by Livingston and Lewis (1995).
#' @return An estimate of the reliability of a test, given the effective test length, mean, variance, and minimum and maximum possible scores of the observed-score distribution..
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # From the data-generating script above, the effective test length is 100.
#' # To estimate and retrieve the model-implied reliability using R.ETL():
#' R.ETL(mean = mean(testdata), variance = var(testdata), min = 0, max = 100,
#' ETL = 100)
#' @export
R.ETL <- function(mean, variance, min = 0, max = 1, ETL) {
  (ETL * variance - (mean - min) * (max - mean)) / ((ETL - 1) * variance)
}

#' An Implementation of the Livingston and Lewis (1995) Approach to Estimate Classification Consistency and Accuracy based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Livingston and Lewis approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to the binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores, or a list specifying parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location-, and \code{alpha} and \code{beta} for the shape parameters of the Beta true-score distribution, and \code{etl} for the effective test length (see documentation for the \code{ETL} function).
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param min The minimum value possible to attain on the test. Default is 0.
#' @param max The maximum value possible to attain on the test. Default is 1 (assumes that the values contained in \code{x} represents proportions of maximum credit).
#' @param cut The cutoff value for classifying observations into pass or fail categories.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param truecut Optional specification of a "true" cutoff. Useful for producing ROC curves (see documentation for the \code{LL.ROC()} function).
#' @param output Character vector indicating which types of statistics (i.e, accuracy and/or consistency) are to be computed and included in the output. Permissible values are \code{"accuracy"} and \code{"consistency"}.
#' @param failsafe Logical value indicating whether to engage the automatic fail-safe defaulting to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates (l < 0, u > 1, non-positive shape parameters, or non-finite estimates). Default is \code{TRUE} (i.e., the function will engage failsafe if the four-parameter Beta-distribution fitting-procedure produced impermissible estimates).
#' @param l If \code{true.model = "2P"} or \code{failsafe = TRUE}, the lower-bound location parameter to be used in the two-parameter fitting procedure. Default is 0 (i.e., the lower-bound of the Standard Beta distribution).
#' @param u If \code{true.model = "2P"} or \code{failsafe = TRUE}, the upper-bound location parameter to be used in the two-parameter fitting procedure. Default is 1 (i.e., the upper-bound of the Standard Beta distribution).
#' @param modelfit Allows for controlling the chi-square test for model fit. The argument takes either a vector of two values, or \code{NULL}. If set to \code{NULL}, the model-fit test is not executed. If a vector of values is supplied, the first value is to represent the initial number of bins the distribution of scores is to be divided in to. This value is set to a default of 100. If this default results in too few bins to conduct the chi-square test, this value can be made larger. The second value represents the minimum expected number of observations that the bins should consist of. In accordance with standard recommendations for chi-square tests, the default value is set to 10. Degrees of freedom are the number of bins minus one, minus the number of moments matched when fitting the true-score distribution (two for a two-parameter solution, four otherwise), following Lord (1965, p. 244). Lord (1965) grouped score points only as far as necessary to obtain expected frequencies of at least 1.0 per class; \code{modelfit = c(nbins, 1)} applies that rule (the default minimum of 10 gives a coarser, less sensitive test). Note that the test is carried out on observed scores rescaled and rounded to the effective test length, as the Livingston and Lewis model defines the observed-score distribution on that scale.
#' @return A list containing the estimated parameters necessary for the approach (i.e., the effective test-length and the beta distribution parameters), a chi-square test of model-fit, the confusion matrix containing estimated proportions of true/false pass/fail categorizations for a test, diagnostic performance statistics, and / or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @note It should be noted that this implementation differs from the original articulation of Livingston and Lewis (1995) in some respects. First, the procedure includes a number of diagnostic performance (accuracy) indices which the original procedure enables but that were not included. Second, the way consistency is calculated differs substantially from the original articulation of the procedure, which made use of a split-half approach. Rather, this implementation uses the approach to estimating classification consistency outlined by Hanson (1991).
#' @note A shiny application providing a GUI for this method is available at https://hthaa.shinyapps.io/shinybeta/ .
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. To estimate and
#' # retrieve the estimated parameters, confusion matrix, consistency and
#' # accuracy statistics using LL.CA():
#' LL.CA(x = testdata, reliability = .7, cut = 50, min = 0, max = 100)
#'
#' # Suppose the true-score parameter estimation procedure arrived at
#' # impermissible parameter estimates (i.e., l < 0, u > 1, alpha < 0, or
#' # beta < 0). For example:
#' set.seed(9)
#' testdata <- rbinom(100, 25, rBeta.4P(100, 0.25, 1, 5, 3))
#' Beta.tp.fit(testdata, 0, 25, 25, failsafe = TRUE)
#'
#' # Suppose further that you have good grounds for assuming that the lower-
#' # bound parameter is equal to 0.25 (e.g., the test consists of multiple-
#' # choice questions with four response options, leading to a 25% probability
#' # of guessing the correct answer per question), and good reason to believe
#' # that the upper-bound parameter is equal to 1 (i.e., there is no reason to
#' # believe that there are no members of the population who will attain a
#' # perfect score across all possible test-forms.) To set these lower and
#' # upper bounds for the fitting procedure in the LL.CA() function, set
#' # the argument true.model = "2p", and specify the location parameters
#' # l = 0.25 and u = 1:
#' LL.CA(testdata, 0.6287713, 12, 0, 25, true.model = "2p", l = 0.25, u = 1)
#'
#' # Alternatively to supplying scores to which a true-score distribution is
#' # to be fit, a list with true-score distribution parameter values can be
#' # supplied manually along with the effective test length (see documentation
#' # for the ETL() function), foregoing the need for actual data. The list
#' # entries must be named. "l" is the lower-bound and "u" the upper-bound
#' # location parameters of the true-score distribution, "alpha" and "beta" for
#' # the shape parameters, and "etl" for the effective test-length..
#' trueparams <- list("l" = 0.25, "u" = 0.75, "alpha" = 5, "beta" = 3, "etl" = 50)
#' LL.CA(x = trueparams, cut = 50, min = 0, max = 100)
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @references Lord. Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Lewis, Don and Burke, C. J. (1949). The Use and Misuse of the Chi-Square Test. Psychological Bulletin, 46(6).
#' @export
LL.CA <- function(x = NULL, reliability, cut, min = 0, max = 1, true.model = "4P", truecut = NULL, output = c("accuracy", "consistency"), failsafe = TRUE, l = 0, u = 1, modelfit = c("nbins" = 100, "minbin" = 10)) {
  model <- .true_model(true.model)
  from_scores <- !is.list(x)
  if (from_scores) {
    .warn_outside_scale(x, min, max)
    N <- ETL(mean(x), stats::var(x), min = min, max = max, reliability = reliability)
    fit <- .ll_fit(x, min, max, N, model, failsafe && model == "4P", l, u)
    params <- fit$params
    .warn_out_of_bounds(params)
  } else {
    params <- x
    N <- params$etl
  }
  if (is.null(truecut)) {
    truecut <- cut
  }
  N <- round(N)
  cut <- (cut - min) / (max - min)
  truecut <- (truecut - min) / (max - min)
  params[["etl_rounded"]] <- N
  out <- list(parameters = params)

  if (from_scores && !is.null(modelfit)) {
    out[["modelfit"]] <- .ll_modelfit((x - min) / (max - min), N, params, modelfit,
                                      n_params = if (model == "2P" || fit$engaged) 2 else 4)
  }

  # Observed scores on the effective-test-length scale are classified as
  # passing when they are at or above round(N * cut).
  observed_cut <- round(N * cut)
  below <- function(p) stats::pbinom(observed_cut - 1, N, p)
  above <- function(p) stats::pbinom(observed_cut - 1, N, p, lower.tail = FALSE)
  integrate_tp <- function(f, lower = 0, upper = 1) {
    .integrate_support(function(p) dBeta.4P(p, params$l, params$u, params$alpha, params$beta) * f(p),
                       params$l, params$u, lower, upper)
  }

  output <- tolower(output)
  if (any(startsWith(output, "a"))) {
    p.tp <- integrate_tp(above, truecut, 1)
    p.fp <- integrate_tp(above, 0, truecut)
    p.ff <- integrate_tp(below, truecut, 1)
    p.tf <- integrate_tp(below, 0, truecut)
    camat <- confmat(p.tf, p.tp, p.ff, p.fp, "prop")
    out[["confusionmatrix"]] <- camat
    out[["classification.accuracy"]] <- caStats(camat[1, 1], camat[1, 2], camat[2, 1], camat[2, 2])
  }
  if (any(startsWith(output, "c"))) {
    p.ij <- integrate_tp(function(p) below(p) * above(p))
    out <- c(out, .binary_consistency(ii = integrate_tp(function(p) below(p)^2), ij = p.ij, ji = p.ij,
                                      jj = integrate_tp(function(p) above(p)^2)))
  }
  out
}

# Fit the true-score distribution for the Livingston and Lewis approach.
.ll_fit <- function(x, min, max, etl, true.model, failsafe, l, u) {
  m <- HB.tsm((x - min) / (max - min) * etl, 4, etl, 0)
  fit <- .fit_true_score(m, true.model, failsafe, l, u)
  fit$params[["etl"]] <- etl
  fit
}

# Chi-square model fit for the Livingston and Lewis approach. `xs` are the
# observed scores on the proportion scale.
.ll_modelfit <- function(xs, N, params, modelfit, n_params) {
  scores <- round(xs * N)
  breaks <- seq(0, N, N / modelfit[1])
  nbins <- length(breaks) - 1
  # Integer score s falls in bin j when breaks[j] <= s < breaks[j + 1]; the
  # last bin also holds everything at or above its lower break. The same rule
  # is applied to observed and expected counts.
  bin_of <- function(s) findInterval(s, breaks[seq_len(nbins)])
  pmf <- .betabinom_pmf(0:N, N, params$l, params$u, params$alpha, params$beta)
  expected <- vapply(seq_len(nbins), function(j) sum(pmf[bin_of(0:N) == j]), numeric(1))
  expected <- expected / sum(expected) * length(xs)
  observed <- tabulate(bin_of(scores), nbins)
  .chisq_modelfit(expected, observed, modelfit[2], n_params)
}

#' An Extension of the Livingston and Lewis (1995) Approach to Estimate Classification Consistency and Accuracy for Multiple Classifications based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Livingston and Lewis approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to the binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores, or a list specifying parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location-, and \code{alpha} and \code{beta} for the shape parameters of the Beta true-score distribution, and \code{etl} for the effective test length (see documentation for the \code{ETL} function).
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param min The minimum value possible to attain on the test. Default is 0.
#' @param max The maximum value possible to attain on the test. Default is 1 (assumes that the values contained in \code{x} represents proportions of maximum credit).
#' @param cut A vector of cut-off values for classifying observations into two or more categories.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param failsafe Logical value indicating whether to engage the automatic fail-safe defaulting to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates. Default is \code{TRUE} (i.e., the function will engage failsafe if the four-parameter Beta-distribution fitting-procedure produced impermissible estimates).
#' @param l If \code{true.model = "2P"} or \code{failsafe = TRUE}, the lower-bound location parameter to be used in the two-parameter fitting procedure. Default is 0 (i.e., the lower-bound of the Standard Beta distribution).
#' @param u If \code{true.model = "2P"} or \code{failsafe = TRUE}, the upper-bound location parameter to be used in the two-parameter fitting procedure. Default is 1 (i.e., the upper-bound of the Standard Beta distribution).
#' @param modelfit Allows for controlling the chi-square test for model fit. The argument takes either a vector of two values, or \code{NULL}. If set to \code{NULL}, the model-fit test is not executed. If a vector of values is supplied, the first value is to represent the initial number of bins the distribution of scores is to be divided in to. This value is set to a default of 100. If this default results in too few bins to conduct the chi-square test, this value can be made larger. The second value represents the minimum expected number of observations that the bins should consist of. In accordance with standard recommendations for chi-square tests, the default value is set to 10. Degrees of freedom are the number of bins minus one, minus the number of moments matched when fitting the true-score distribution (two for a two-parameter solution, four otherwise), following Lord (1965, p. 244). Lord (1965) grouped score points only as far as necessary to obtain expected frequencies of at least 1.0 per class; \code{modelfit = c(nbins, 1)} applies that rule (the default minimum of 10 gives a coarser, less sensitive test). Note that the test is carried out on observed scores rescaled and rounded to the effective test length, as the Livingston and Lewis model defines the observed-score distribution on that scale.
#' @return A list containing the estimated parameters necessary for the approach (i.e., the effective test-length and the beta distribution parameters), a chi-square test of model-fit, the confusion matrix containing estimated proportions of true/false positive/negative categorizations for a test, diagnostic performance statistics, and/or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @note It should be noted that this implementation differs from the original articulation of Livingston and Lewis (1995) in some respects. First, the procedure includes a number of diagnostic performance (accuracy) indices which the original procedure enables but that were not included. Second, the way consistency is calculated differs substantially from the original articulation of the procedure, which made use of a split-half approach. Rather, this implementation uses the approach to estimating classification consistency outlined by Hanson (1991).
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.1, 0.95, 5, 3)
#' for (i in 1:100) {
#'   if (i == 1) {
#'     rawdata <- matrix(nrow = 1000, ncol = 100)
#'   }
#'   rawdata[, i] <- rbinom(1000, 1, p.success)
#' }
#'
#' # Suppose the cutoff value for being placed in the lower category is a score
#' # below 50, second lowest 60, then 70, 80, and 90. Using the cba() function
#' # to estimate the reliability of this test, to use the LL.CA.MC() function
#' # or estimating diagnostic performance and consistency indices of
#' # classifications when using several cut-points:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)
#'
#' # The output from this function can get quite verbose when operating with
#' # several cut-points. In order to retrieve only model parameter estimates:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)$parameters
#'
#' # To retrieve only the model-fit estimate:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)$modelfit
#'
#' # To retrieve only the diagnostic performance estimates:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)$accuracy
#'
#' # To retrieve only the classification consistency indices:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)$consistency
#'
#' # Alternatively, the MC.out.tabular() function can be used to organize the
#' # category-specific indices in a tabular format:
#' MC.out.tabular(LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100))
#'
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @references Lord. Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Lewis, Don and Burke, C. J. (1949). The Use and Misuse of the Chi-Square Test. Psychological Bulletin, 46(6).
#' @export
LL.CA.MC <- function(x = NULL, reliability, cut, min = 0, max = 1, true.model = "4P", failsafe = TRUE, l = 0, u = 1, modelfit = c("nbins" = 100, "minbin" = 10)) {
  model <- .true_model(true.model)
  from_scores <- !is.list(x)
  if (from_scores) {
    .warn_outside_scale(x, min, max)
    N <- ETL(mean(x), stats::var(x), min = min, max = max, reliability = reliability)
    fit <- .ll_fit(x, min, max, N, model, failsafe && model == "4P", l, u)
    params <- fit$params
    .warn_out_of_bounds(params)
  } else {
    params <- x
    N <- params$etl
  }
  labels <- cut
  cut <- (cut - min) / (max - min)
  N <- round(N)
  params[["etl_rounded"]] <- N
  out <- list(parameters = params)

  if (from_scores && !is.null(modelfit)) {
    out[["modelfit"]] <- .ll_modelfit((x - min) / (max - min), N, params, modelfit,
                                      n_params = if (model == "2P" || fit$engaged) 2 else 4)
  }

  pcut <- c(0, cut, 1)
  ocut <- round(pcut * N)
  K <- length(cut) + 1
  # Probability that the observed score falls in category j given true score p.
  category_prob <- function(p, j) {
    upper <- if (j == K) 1 else stats::pbinom(ocut[j + 1] - 1, N, p)
    upper - stats::pbinom(ocut[j] - 1, N, p)
  }
  integrate_tp <- function(f, lower = 0, upper = 1) {
    .integrate_support(function(p) dBeta.4P(p, params$l, params$u, params$alpha, params$beta) * f(p),
                       params$l, params$u, lower, upper)
  }

  # Accuracy: rows are observed categories, columns are true categories.
  camat <- matrix(nrow = K, ncol = K,
                  dimnames = list(.category_labels(labels, "Observed"), .category_labels(labels, "True    ")))
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      camat[j, i] <- integrate_tp(function(p) category_prob(p, j), pcut[i], pcut[i + 1])
    }
  }
  out[["accuracy"]] <- .mc_accuracy(camat / sum(camat))

  # Consistency: rows are the second administration, columns the first.
  ccmat <- matrix(nrow = K, ncol = K,
                  dimnames = list(.category_labels(labels, "2nd adm."), .category_labels(labels, "1st adm.")))
  for (i in seq_len(K)) {
    for (j in i:K) {
      ccmat[j, i] <- ccmat[i, j] <- integrate_tp(function(p) category_prob(p, j) * category_prob(p, i))
    }
  }
  out[["consistency"]] <- .mc_consistency(ccmat / sum(ccmat))
  out
}

#' Confusion matrix
#'
#' @description Organizes supplied values of true and false positives and negatives into a confusion matrix.
#' @param tp The frequency or rate of true-positive classifications.
#' @param tn The frequency or rate of true-negative classifications.
#' @param fp The frequency or rate of false-positive classifications.
#' @param fn The frequency or rate of false-negative classifications.
#' @param output Whether the returned output reflects frequencies or proportions. Defaults to returning frequencies.
#' @return A confusion matrix organizing the input values of true and false positive and negatives.
#' @examples
#' # Generate some true and observed conditions.
#' set.seed(1234)
#' true.ability <- rbeta(50, 4, 4)
#' true.category <- ifelse(true.ability < 0.5, 0, 1)
#' observed.score <- rbinom(50, 50, true.ability)
#' observed.category <- ifelse(observed.score < 25, 0, 1)
#' # Calculate the frequencies of true and false positives and negatives based on the true and
#' # observed conditions.
#' TP <- sum(ifelse(observed.category == 0 & true.category == 0, 1, 0))
#' FP <- sum(ifelse(observed.category == 0 & true.category != 0, 1, 0))
#' TN <- sum(ifelse(observed.category == 1 & true.category == 1, 1, 0))
#' FN <- sum(ifelse(observed.category == 1 & true.category != 1, 1, 0))
#' # Organize the above values in a confusion matrix using the confmat function:
#' confmat(tp = TP, fp = FP, tn = TN, fn = FN)
#' @export
confmat <- function(tp, tn, fp, fn, output = "freq") {
  mat <- matrix(c(tp, fp, tn, fn), nrow = 2,
                dimnames = list(c("True", "False"), c("Positive", "Negative")))
  mat <- cbind(mat, "Total" = rowSums(mat))
  mat <- rbind(mat, "Total" = colSums(mat))
  if (output != "freq") {
    mat <- mat / sum(tp, tn, fp, fn)
  }
  mat
}

#' Classification Accuracy Statistics.
#'
#' @description Provides a set of statistics often used for conveying information regarding the certainty of classifications based on tests.
#' @param tp The frequency or rate of true-positive classifications.
#' @param tn The frequency or rate of true-negative classifications.
#' @param fp The frequency or rate of false-positive classifications.
#' @param fn The frequency or rate of false-negative classifications.
#' @return A list of diagnostic performance statistics based on true/false positive/negative statistics. Specifically, the sensitivity, specificity, positive predictive value (PPV), negative predictive value (NPV), Youden's J. (Youden.J), and Accuracy.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. First, compute the
#' # estimated confusion matrix using LL.CA():
#' cmat <- LL.CA(x = testdata, reliability = 0.7, cut = 50, min = 0,
#' max = 100)$confusionmatrix
#'
#' # To estimate and retrieve diagnostic performance statistics using caStats(),
#' # feed it the appropriate entries of the confusion matrix.
#' caStats(tp = cmat["True", "Positive"], tn = cmat["True", "Negative"],
#' fp = cmat["False", "Positive"], fn = cmat["False", "Negative"])
#' @export
caStats <- function(tp, tn, fp, fn) {
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  list("Sensitivity" = sensitivity,
       "Specificity" = specificity,
       "PPV" = tp / (tp + fp),
       "NPV" = tn / (tn + fn),
       "Youden.J" = sensitivity + specificity - 1,
       "Accuracy" = (tp + tn) / (tp + tn + fp + fn))
}

#' Classification Consistency Statistics.
#'
#' @description Provides a set of statistics often used for conveying information regarding the consistency of classifications based on tests.
#' @param ii The frequency or rate of consistent classifications into category "i". Frequencies are converted to proportions before the statistics are computed.
#' @param ij The frequency or rate of inconsistent classifications into categories "i" and "j".
#' @param ji The frequency or rate of inconsistent classifications into categories "j" and "i".
#' @param jj The frequency or rate of consistent classifications into category "j".
#' @return A list of classification consistency statistics. Specifically, the coefficient of consistent classification (p), the coefficient of consistent classification by chance (p_c), the proportion of positive classifications due to chance (p_c_pos), the proportion of negative classifications due to chance (p_c_neg), and Cohen's Kappa coefficient.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. First, compute the
#' # estimated consistency matrix using LL.CA():
#' cmat <- LL.CA(x = testdata, reliability = .7, cut = 50, min = 0,
#' max = 100)$consistencymatrix
#'
#' # To estimate and retrieve consistency statistics using ccStats(),
#' # feed it the appropriate entries of the consistency matrix.
#' ccStats(ii = cmat["i", "i"], ij = cmat["i", "j"],
#' ji = cmat["j", "i"], jj = cmat["j", "j"])
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @export
ccStats <- function(ii, ij, ji, jj) {
  # Work with proportions so that frequencies and rates give the same result.
  total <- ii + ij + ji + jj
  ii <- ii / total
  ij <- ij / total
  ji <- ji / total
  jj <- jj / total
  p <- ii + jj
  p_c_pos <- (ii + ij) * (ii + ji)
  p_c_neg <- (ij + jj) * (ji + jj)
  p_c <- p_c_pos + p_c_neg
  list("p" = p, "p_c" = p_c, "p_c_pos" = p_c_pos, "p_c_neg" = p_c_neg, "Kappa" = (p - p_c) / (1 - p_c))
}

#' ROC curves for the Livingston and Lewis approach.
#'
#' @description Generate a ROC curve plotting the false-positive rate against the true-positive rate at different cut-off values across the observed-score scale.
#' @param x A vector of observed results.
#' @param min The minimum possible value to attain on the observed-score scale.
#' @param max The maximum value possible to attain on the test. Default is 1 (assumes that the values contained in \code{x} represents proportions of maximum credit).
#' @param reliability The reliability coefficient of the test.
#' @param truecut The true point along the x-scale that marks the categorization-threshold.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The \code{"4P"} method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the \code{"2P"} method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param failsafe If true-model == "4P": Whether to engage a fail-safe reverting to a two-parameter true-score distribution solution should the four-parameter fitting procedure produce impermissible results. Default is TRUE (engage fail-safe in the event of impermissible estimates).
#' @param l If \code{true.model == "2P"} or \code{failsafe == TRUE}: The lower-bound location parameter of the two-parameter true-score distribution solution.
#' @param u If \code{true.model == "2P"} or \code{failsafe == TRUE}: The upper-bound location parameter of the two-parameter true-score distribution solution.
#' @param AUC Logical. Calculate and include the area under the curve? Default is \code{FALSE}.
#' @param maxJ Logical. Mark the point along the curve where Youden's J statistic is maximized? Default is \code{FALSE}.
#' @param maxAcc Logical. Mark the point along the curve where the Accuracy statistic is maximized? Default is \code{FALSE}.
#' @param locate Ask the function to locate the cut-point at which sensitivity or NPV is greater than or equal to some value, or specificity or PPV is lesser than or equal to some value. Take as input a character-vector of length 2, with the first argument being which index is to be found (e.g., "sensitivity"), and the second argument the value to locate (e.g., "0.75"). For example: c("sensitivity", "0.75").
#' @param raw.out Give raw coordinates as output rather than plot? Default is \code{FALSE}
#' @param grainsize Specify the number of cutoff-points for which the ROC curve is to be calculated. The greater this number the greater the accuracy. Default is 100 points.
#' @return A plot tracing the ROC curve for the test, or matrix of coordinates if raw.out is \code{TRUE}.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' hist(testdata / 100, xlim = c(0, 1), freq = FALSE)
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct.
#' # Suppose further that the reliability of the test-scores were estimated to
#' # 0.75. To produce a plot with an ROC curve using LL.ROC(), along with the
#' # AUC statistics and the points at which Youden's J. is maximized:
#' LL.ROC(x = testdata, reliability = 0.7, truecut = 50, min = 0, max = 100,
#' AUC = TRUE, maxJ = TRUE)
#' # Or to locate the point at which accuracy is maximized:
#' LL.ROC(x = testdata, reliability = 0.7, truecut = 50, min = 0, max = 100,
#' maxAcc = TRUE)
#'
#' # Using the example data above, the function can be instructed to locate an
#' # operational cut-point at which sensitivity or specificity is equal to or
#' # greater than some specified value by specifying the "locate" argument with
#' # c("statistic", value). For example, to locate the operational cut-point at
#' # which sensitivity is first equal to or greater than 0.9:
#' LL.ROC(testdata, reliability = 0.7, min = 0, max = 100, truecut = 50,
#' locate = c("sensitivity", 0.9))
#' # For Negative Predictive value, the point at which it is equal or greater:
#' LL.ROC(testdata, reliability = 0.7, min = 0, max = 100, truecut = 50,
#' locate = c("NPV", 0.9))
#' # And so on for other statistics such as Specificity and Positive Predictive
#' # Value.
#' @export
LL.ROC <- function(x = NULL, reliability, min = 0, max = 1, truecut, true.model = "4P", failsafe = TRUE, l = 0, u = 1, AUC = FALSE, maxJ = FALSE, maxAcc = FALSE, locate = NULL, raw.out = FALSE, grainsize = 100) {
  if (!is.list(x)) {
    x <- Beta.tp.fit(x, min, max, reliability = reliability, true.model = true.model, failsafe = failsafe, l = l, u = u)
  }
  rocmat <- .roc_matrix(seq(min, max, (max - min) / grainsize), function(cut) {
    LL.CA(x = x, min = min, max = max, cut = cut, truecut = truecut, true.model = true.model,
          output = "a", l = l, u = u, modelfit = NULL)$classification.accuracy
  })
  if (raw.out) {
    return(rocmat)
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  .roc_plot(rocmat, truecut, AUC = AUC, maxJ = maxJ, maxAcc = maxAcc, locate = locate)
}

#' Area Under the ROC Curve.
#'
#' @description Given a vector of false-positive rates and a vector of true-positive rates, calculate the area under the Receiver Operator Characteristic (ROC) curve.
#' @param FPR Vector of False-Positive Rates.
#' @param TPR Vector of True-Positive Rates.
#' @return A value representing the area under the ROC curve.
#' @note Script originally retrieved and modified from https://blog.revolutionanalytics.com/2016/11/calculating-auc.html.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. To calculate the
#' # necessary (x, y) coordinates to compute the area under the curve statistic
#' # one can use the LL.ROC() function with the argument
#' # raw.out = TRUE.
#' coords <- LL.ROC(x = testdata, reliability = .7, truecut = 50, min = 0,
#' max = 100, raw.out = TRUE)
#'
#' # To calculate and retrieve the Area Under the Curve (AUC) with the AUC()
#' # function, feed it the raw coordinates calculated above.
#' AUC(coords[, "FPR"], coords[, "TPR"])
#' @export
AUC <- function(FPR, TPR) {
  # Trapezoidal rule over the (FPR, TPR) coordinates in the order given.
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR) / 2
}

#' Calculate Cronbach's Alpha reliability-coefficient from supplied variables.
#'
#' @description Calculates Cronbach's Alpha reliability coefficient of the sum-score.
#' @param x A data-frame or matrix of numerical values where rows represent respondents, and columns represent items.
#' @note Missing values are treated by passing \code{na.rm = TRUE} to the \code{var} function call.
#' @note Be aware that this function does not issue a warning if there are negative correlations between variables in the supplied data-set.
#' @return Cronbach's Alpha for the sum-score of supplied variables.
#' @references Cronbach, L.J. (1951). Coefficient alpha and the internal structure of tests. Psychometrika 16, 297--334. doi: 10.1007/BF02310555
#' @examples
#' # Generate some fictional data. Say 100 students take a 50-item long test
#' # where all items are equally difficult.
#' set.seed(1234)
#' p.success <- rBeta.4P(100, 0.25, 0.75, 5, 3)
#' for (i in 1:50) {
#'   if (i == 1) {
#'     rawdata <- matrix(nrow = 100, ncol = 50)
#'   }
#'   rawdata[, i] <- rbinom(100, 1, p.success)
#' }
#' # To calculate Cronbach's Alpha for this test:
#' cba(rawdata)
#' @export
cba <- function(x) {
  k <- ncol(x)
  S <- stats::var(x, na.rm = TRUE)
  k / (k - 1) * (1 - sum(diag(S)) / sum(S))
}

#' Calculate McDonald's Omega reliability-coefficient from supplied variables.
#'
#' @description Calculates McDonalds's Omega reliability-coefficient of the sum-score from the Spearman one-factor model using the procedure outlined in McDonald (1999).
#' @param x A data-frame or matrix of numerical values where rows represent respondents, and columns represent items.
#' @param fit Logical. Default is \code{FALSE}. If \code{TRUE}, the output changes from a vector containing the Omega reliability-estimate to a list containing additional detailed information concerning the fitted factor model.
#' @note Missing values are treated by passing \code{na.rm = TRUE} to the \code{var} function call and \code{use = "pairwise.complete.obs"} to the \code{cov} function call.
#' @note The function terminates with an error if there are negative covariance-matrix entries.
#' @return If \code{fit = FALSE}, A vector of length 1 containing the estimated McDonalds's Omega reliability-coefficient for the sum-score of the supplied variables. If \code{fit = TRUE}, a list containing the Omega-coefficient reliability-estimate as the first entry, followed by the goodness-of-fit index (GFI), a two-row matrix containing the estimated factor-loadings and error-variances, and the observed and fitted covariance-matrices and the discrepancy matrix.
#' @references McDonald, R. P. (1999). Test Theory: A Unified Treatment. Routledge.
#' @examples
#' # Generate some fictional data.
#' set.seed(1234)
#' rawdata <- matrix(rnorm(500), ncol = 5)
#' common <- rnorm(100)
#' rawdata <- apply(rawdata, 2, function(x) {x + common})
#'
#' # To estimate McDonald's Omega from this data:
#' mdo(rawdata)
#'
#' # To retrieve additional information such as the GFI fit-index and model-
#' # parameter estimates:
#' mdo(rawdata, fit = TRUE)
#' @export
mdo <- function(x, fit = FALSE) {
  n <- ncol(x)
  if (n < 3) {
    stop("At least three variables are required to estimate McDonald's Omega.")
  }
  omat <- stats::cov(x, use = "pairwise.complete.obs")
  if (any(omat <= 0)) {
    stop("Item-covariance(s) less than or equal to 0 detected. Consider reverse-scoring or excluding variables.")
  }
  vars <- apply(x, 2, stats::var, na.rm = TRUE)
  # Spearman one-factor model: for item i and any two other items a and b,
  # loading_i^2 = cov(i, a) * cov(i, b) / cov(a, b). Each loading is the mean
  # of these triad estimates over all pairs of other items.
  floads <- vapply(seq_len(n), function(i) {
    others <- seq_len(n)[-i]
    pairs <- which(upper.tri(diag(n - 1)), arr.ind = TRUE)
    a <- others[pairs[, 1]]
    b <- others[pairs[, 2]]
    mean(sqrt(omat[i, a] * omat[i, b] / omat[cbind(a, b)]))
  }, numeric(1))
  evars <- vars - floads^2
  Omega <- sum(floads)^2 / (sum(evars) + sum(floads)^2)
  if (!fit) {
    return(Omega)
  }
  imat <- outer(floads, floads)
  diag(imat) <- floads^2 + evars
  dmat <- omat - imat
  # Goodness-of-fit index (unweighted least squares form).
  gfi <- 1 - sum(dmat^2) / sum(omat^2)
  floads <- matrix(c(floads, evars), nrow = 2, byrow = TRUE, dimnames = list(c("Loadings", "Errors")))
  if (!is.null(colnames(x))) {
    colnames(floads) <- rownames(imat) <- colnames(imat) <- colnames(x)
  }
  list("Omega" = Omega,
       "GFI" = gfi,
       "Parameters" = floads,
       "Matrices" = list("Observed" = omat, "Fitted" = imat, "Discrepancy" = dmat))
}

#' Estimate Beta true-score distribution based on observed-score raw-moments and the effective test length.
#'
#' @description Estimator for the Beta true-score distribution shape-parameters from the observed-score distribution and Livingston and Lewis' effective test length. Returns a list with entries representing the lower- and upper shape parameters (l and u), and the shape parameters (alpha and beta) of the four-parameters beta distribution, and the effective test length.
#' @param x Vector of observed-scores.
#' @param min The minimum possible score to attain on the test.
#' @param max The maximum possible score to attain on the test.
#' @param etl The value of Livingston and Lewis' effective test length. See \code{?ETL()}. Not necessary to specify if reliability is supplied to the \code{reliability} argument.
#' @param reliability Optional specification of the test-score reliability coefficient. Used to compute the effective test length when \code{etl} is not supplied.
#' @param true.model The type of Beta distribution which is to be fit to the moments of the true-score distribution. Options are \code{"4P"} and \code{"2P"}, where \code{"4P"} refers to the four-parameter (with the same mean, variance, skewness, and kurtosis), and \code{"2P"} the two-parameter solution where both location-parameters are specified (with the same mean and variance).
#' @param failsafe Logical. Whether to revert to a fail-safe two-parameter solution should the four-parameter solution contain invalid parameter estimates.
#' @param l If \code{failsafe = TRUE} or \code{true.model = "2P"}: The lower-bound of the Beta distribution. Default is 0 (i.e., the lower-bound of the Standard, two-parameter Beta distribution).
#' @param u If \code{failsafe = TRUE} or \code{true.model = "2P"}: The upper-bound of the Beta distribution. Default is 1 (i.e., the upper-bound of the Standard, two-parameter Beta distribution).
#' @param output Option to specify true-score distribution moments as output if the value of the output argument does not equal \code{"parameters"}.
#' @return A list with the parameter values of a four-parameter Beta distribution. "l" is the lower location-parameter, "u" the upper location-parameter, "alpha" the first shape-parameter, "beta" the second shape-parameter, and "etl" the effective test length.
#' @references Hanson, B. A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series. Retrieved from https://files.eric.ed.gov/fulltext/ED344945.pdf
#' @references Lord, F. M. (1965). A strong true-score theory, with applications. Psychometrika. 30(3). pp. 239--270. doi: 10.1007/BF02289490
#' @references Rogosa, D. &  Finkelman, M. (2004). How Accurate Are the STAR Scores for Individual Students? An Interpretive Guide. Retrieved from http://statweb.stanford.edu/~rag/accguide/guide04.pdf
#' @examples
#' # Generate some fictional data. Say 1000 individuals take a 100-item test
#' # where all items are equally difficult, and the true-score distribution
#' # is a four-parameter Beta distribution with location parameters l = 0.25,
#' # u = 0.75, alpha = 5, and beta = 3:
#' set.seed(12)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#'
#' # Since this test contains items which are all equally difficult, the true
#' # effective test length (etl) is the actual test length. I.e., etl = 100.
#' # To estimate the four-parameter Beta distribution parameters underlying
#' # the draws from the binomial distribution:
#' Beta.tp.fit(testdata, 0, 100, 100)
#'
#' # Imagine a case where the fitting procedure produces an impermissible
#' # estimate (e.g., l < 0 or u > 1).
#' set.seed(1234)
#' testdata <- rbinom(1000, 50, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' Beta.tp.fit(testdata, 0, 50, 50)
#'
#' # This example produced an l-value estimate less than 0. One way of
#' # dealing with such an occurrence is to revert to a two-parameter
#' # model, specifying the l and u parameters and estimating the
#' # alpha and beta parameters necessary to produce a Beta distribution
#' # with the same mean and variance as the estimated true-score distribution.
#'
#' # Suppose you have good theoretical reasons to fix the l parameter at a
#' # value of 0.25 (e.g., the test is composed of multiple-choice questions
#' # with four response-options, resulting in a 25% chance of guessing the
#' # correct answer). The l-parameter could be specified to this theoretically
#' # justified value, and the u-parameter could be specified to be equal to the
#' # estimate above (u = 0.7256552) as such:
#' Beta.tp.fit(testdata, 0, 50, 50, true.model = "2P", l = 0.25, u = 0.7256552)
#' @export
Beta.tp.fit <- function(x, min, max, etl = NULL, reliability = NULL, true.model = "4P", failsafe = FALSE, l = 0, u = 1, output = "parameters") {
  if (is.null(etl)) {
    etl <- ETL(mean(x), stats::var(x), min, max, reliability)
  }
  if (output == "parameters") {
    return(.ll_fit(x, min, max, etl, true.model, failsafe, l, u)$params)
  }
  m <- HB.tsm((x - min) / (max - min) * etl, 4, etl, 0)
  s2 <- m[2] - m[1]^2
  s3 <- m[3] - 3 * m[1] * m[2] + 2 * m[1]^3
  s4 <- m[4] - 4 * m[1] * m[3] + 6 * m[1]^2 * m[2] - 3 * m[1]^4
  list("Raw" = list(m[1], m[2], m[3], m[4]),
       "Central" = list(0, s2, s3, s4),
       "Standardized" = list(0, 1, s3 / sqrt(s2)^3, s4 / sqrt(s2)^4))
}

#' Estimate Beta True-Score Distribution Based on Observed-Score Raw-Moments and Lord's k.
#'
#' @description Estimator for the Beta true-score distribution shape-parameters from the observed-score distribution and Lord's k. Returns a list with entries representing the lower- and upper shape parameters (l and u), and the shape parameters (alpha and beta) of the four-parameters beta distribution, as well as Lord's k and the test length.
#' @param x Vector of observed-scores.
#' @param N The test length.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param true.model The type of Beta distribution which is to be fit to the moments of the true-score distribution. Options are \code{"4P"} and \code{"2P"}, where \code{"4P"} refers to the four-parameter (with the same mean, variance, skewness, and kurtosis), and \code{"2P"} the two-parameter solution where both location-parameters are specified (with the same mean and variance).
#' @param failsafe Logical. Whether to revert to a fail-safe two-parameter solution should the four-parameter solution contain invalid parameter estimates.
#' @param l If \code{failsafe = TRUE} or \code{true.model = "2P"}: The lower-bound of the Beta distribution. Default is 0 (i.e., the lower-bound of the Standard, two-parameter Beta distribution).
#' @param u If \code{failsafe = TRUE} or \code{true.model = "2P"}: The upper-bound of the Beta distribution. Default is 1 (i.e., the upper-bound of the Standard, two-parameter Beta distribution).
#' @return A list with the parameter values of a four-parameter Beta distribution. "l" is the lower location-parameter, "u" the upper location-parameter, "alpha" the first shape-parameter, and "beta" the second shape-parameter. Also includes Lord's k and the test length.
#' @references Hanson, B. A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series. Retrieved from https://files.eric.ed.gov/fulltext/ED344945.pdf
#' @references Lord, F. M. (1965). A strong true-score theory, with applications. Psychometrika. 30(3). pp. 239--270. doi: 10.1007/BF02289490
#' @examples
#' # Generate some fictional data. Say 1000 individuals take a 100-item test
#' # where all items are equally difficult, and the true-score distribution
#' # is a four-parameter Beta distribution with location parameters l = 0.25,
#' # u = 0.75, alpha = 5, and beta = 3, and the error distribution is Binomial
#' # with Lord's k = 0:
#' set.seed(12)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#'
#' # To estimate the four-parameter Beta distribution parameters from this
#' # sample of observations:
#' HB.beta.tp.fit(testdata, 100, 0)
#' @export
HB.beta.tp.fit <- function(x, N, k, true.model = "4P", failsafe = FALSE, l = 0, u = 1) {
  .hb_fit(x, N, k, true.model, failsafe, l, u)$params
}

# Fit the true-score distribution for the Hanson and Brennan approach.
.hb_fit <- function(x, N, k, true.model, failsafe, l, u) {
  fit <- .fit_true_score(HB.tsm(x, 4, N, k), true.model, failsafe, l, u)
  fit$params[["k"]] <- k
  fit$params[["N"]] <- N
  fit
}

#' Descending (falling) factorial.
#'
#' @description Calculate the descending (or falling) factorial of a value \code{x} of order \code{r}.
#' @param x A value for which the descending factorial is to be calculated.
#' @param r The power \code{x} is to be raised to.
#' @return The descending factorial of value \code{x} raised to the \code{r}'th power.
#' @param method The method by which the descending factorials are to be calculated. Default is \code{"product"} which uses direct arithmetic. Alternative is \code{"gamma"} which calculates the ascending factorial using the Gamma function. The alternative method might be faster but might fail because the Gamma function is not defined for negative integers (returning Inf).
#' @export
#' @examples
#' # To calculate the 4th descending factorial for a value (e.g., 3.14):
#' dfac(x = 3.14, r = 4)
#'
#' # To calculate the 5th descending factorial for values 3.14, 2.72, and 0.58:
#' dfac(x = c(3.14, 2.72, 0.58), r = 5)
dfac <- function(x, r, method = "product") {
  if (method != "product") {
    return(gamma(x + 1) / gamma(x - r + 1))
  }
  # Values below the order are set to zero (as the falling factorial of an
  # integer smaller than r is zero).
  x <- ifelse(x < r, 0, x)
  if (r <= 1) {
    return(x^r)
  }
  out <- x
  for (i in seq_len(r - 1)) {
    out <- out * (x - i)
  }
  out
}

#' Ascending (rising) factorial.
#'
#' @description Calculate the ascending (or rising) factorial of a value \code{x} of order \code{r}.
#' @param x A value for which the ascending factorial is to be calculated.
#' @param r The power \code{x} is to be raised to.
#' @param method The method by which the descending factorials are to be calculated. Default is \code{"product"} which uses direct arithmetic. Alternative is \code{"gamma"} which calculates the descending factorial using the Gamma function. The alternative method might be faster but might fail because the Gamma function is not defined for negative integers (returning Inf).
#' @return The ascending factorial of value \code{x} raised to the \code{r}'th power.
#' @export
#' @examples
#' # To calculate the 4th ascending factorial for a value (e.g., 3.14):
#' afac(x = 3.14, r = 4)
#'
#' # To calculate the 5th ascending factorial for values 3.14, 2.72, and 0.58:
#' afac(x = c(3.14, 2.72, 0.58), r = 5)
afac <- function(x, r, method = "product") {
  if (method != "product") {
    return(gamma(x + r) / gamma(x))
  }
  if (r <= 1) {
    return(x^r)
  }
  out <- x
  for (i in seq_len(r - 1)) {
    out <- out * (x + i)
  }
  out
}

#' Proportional true-score distribution raw moments from Livingston and Lewis' effective test-score and effective test-length.
#'
#' @description An implementation of Lords (1965, p. 265) equation 37 for estimating the raw moments of the true-score distribution, modified to function for the Livingston and Lewis approach.
#' @param x The effective test-score of test-takers.
#' @param r The moment-order that is to be calculated (where 1 is the mean, 2 is the raw variance, 3 is the raw skewness, etc.).
#' @param n The effective test-length.
#' @param method The method by which the descending factorials are to be calculated. Default is \code{"product"} which uses direct arithmetic. Alternative is "gamma" which calculates the descending factorial using the Gamma function. The alternative method might be faster but might fail because the Gamma function is not defined for negative integers (returning Inf).
#' @references Lord, F. M. (1965). A strong true-score theory, with applications. Psychometrika. 30(3). pp. 239--270. doi: 10.1007/BF02289490
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @export
#' @examples
#' # Examine the raw moments of the underlying Beta distribution that is to provide the basis for
#' # observed-scores:
#' betamoments(alpha = 5, beta = 3, l = 0.25, u = 0.75, types = "raw")
#'
#' # Generate observed-scores from true-scores by passing the true-scores as binomial probabilities
#' # for the rbinom function.
#' set.seed(1234)
#' obs.scores <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' # Examine the raw moments of the observed-score distribution.
#' observedmoments(obs.scores, type = "raw")
#'
#' # First four estimated raw moment of the proportional true-score distribution from the observed-
#' # score distribution. As all items are equally difficult, the effective test-length is equal to
#' # the actual test-length.
#' tsm(x = obs.scores, r = 1, n = 100)
#' tsm(x = obs.scores, r = 2, n = 100)
#' tsm(x = obs.scores, r = 3, n = 100)
#' tsm(x = obs.scores, r = 4, n = 100)
#' # Which is fairly close to the true raw moments of the proportional true-score distribution
#' # calculated above.
tsm <- function(x, r, n, method = "product") {
  if (method == "product" && r == 1) {
    return(mean(x) / n)
  }
  # Lord (1965), eq. 37 with k = 0: E[p^r] = E[x^(r)] / n^(r), where the
  # falling factorial n^(r) is written as n^(2) * (n - 2)^(r - 2).
  mean(dfac(x, r, method)) / dfac(n - 2, r - 2, method) / dfac(n, 2, method)
}

#' Proportional True-Score Distribution Raw Moments for the Hanson-Brennan Approach to Classification Accuracy and Consistency.
#'
#' @description An implementation of Lords (1965, p. 265) equation 37 for estimating the raw moments of the true-score distribution.
#' @param x Vector of values representing sum-scores.
#' @param r The number of raw moments to be calculated.
#' @param N The number of test items (i.e., test length).
#' @param k Lord's k (see documentation for the \code{Lords.k()} function.
#' @export
#' @examples
#' # Generate some data under the Beta Compound-Binomial distribution, where the
#' # Compound Binomial distribution has 100 trials and Lord's k = 2, and the
#' # Beta distribution has location parameters l = .15 and u = .85, and shape
#' # parameters alpha = 6 and beta = 4:
#' obs <- rBetacBinom(1000, 100, 2, .15, .85, 6, 4)
#'
#' # To estimate the first four raw moments of the underlying Beta distribution:
#' HB.tsm(x = obs, r = 4, N = 100, k = 2)
HB.tsm <- function(x, r, N, k) {
  m <- numeric(r)
  m[1] <- mean(x) / N
  for (i in seq_len(r)[-1]) {
    m[i] <- (mean(dfac(x, i)) / dfac(N - 2, i - 2) + k * dfac(i, 2) * m[i - 1]) /
      (dfac(N, 2) + k * dfac(i, 2))
  }
  m
}

#' Tabular organization of accuracy and consistency output from the \code{LL.CA.MC()} function.
#'
#' @description Function that takes the output from the \code{LL.CA.MC()} function and organizes it in a table with accuracy and consistency indices represented by columns and categories as rows.
#' @param x The list-output from the \code{LL.CA.MC()} function.
#' @export
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.1, 0.95, 5, 3)
#' for (i in 1:100) {
#'   if (i == 1) {
#'     rawdata <- matrix(nrow = 1000, ncol = 100)
#'   }
#'   rawdata[, i] <- rbinom(1000, 1, p.success)
#' }
#'
#' # Estimate accuracy and consistency where the lowest category are scores
#' # below 50, second lowest 60, then 70, 80, and 90. Using the cba() function
#' # to estimate the reliability of this test, to use the LL.CA.MC() function
#' # or estimating diagnostic performance and consistency indices of
#' # classifications when using several cut-points:
#' output <- LL.CA.MC(rowSums(rawdata), cba(rawdata), seq(50, 90, 10), 0, 100)
#'
#' # As this output can get quite verbose as the number of categories increase,
#' # the MC.out.tabular() function can be used to organize the output more
#' # concisely in a tabular format.
#' MC.out.tabular(output)
MC.out.tabular <- function(x) {
  categories <- x$accuracy$specific
  rows <- lapply(seq_along(categories), function(i) {
    cm <- categories[[i]]$confusionmatrix
    c("TP" = cm["True", "Positive"], "FP" = cm["False", "Positive"],
      "TN" = cm["True", "Negative"], "FN" = cm["False", "Negative"],
      unlist(categories[[i]]$statistics),
      unlist(x$consistency$specific[[i]]$statistics))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- paste0("Category.", seq_along(categories))
  out
}

#' Graphical presentation of model fit for the Beta-Binomial classification accuracy and consistency model.
#'
#' @description Tool for visually gauging the discrepancy between the observed and model-implied frequencies of observed-scores.
#' @param x The output object from the \code{LL.CA()}, \code{LL.MC.CA()}, \code{HB.CA()}, or \code{HB.CA.MC()} functions.
#' @param x.tickat The points along the x-axis that bins are to be labeled. Default is \code{NULL} (places a tick for each of the bins).
#' @param y.tickat The points along the y-axis where frequencies are to be labelled. Default is \code{NULL}.
#' @param y.lim The limits of the y-axis (frequencies). Useful for keeping the scale equal across several plots.
#' @param main.lab The main label (title) of the plot.
#' @param x.lab The label for the x-axis (the bins).
#' @param y.lab The label for the y-axis (the frequencies).
#' @param x.grid Control the vertical grid-lines of the plot. Takes \code{NULL}, \code{NA}, or a vector of values as input. If \code{NULL}, grid-lines are drawn automatically for each bin. If \code{NA}, no grid-lines are drawn. If a vector of values are supplied, lines are drawn at each value provided along the x-axis.
#' @param y.grid Control the horizontal grid-lines of the plot. Takes \code{NULL}, \code{NA}, or a vector of values as input. If \code{NULL}, grid-lines are drawn automatically for each frequency (i.e., increments of 1). If \code{NA}, no grid-lines are drawn. If a vector of values are supplied, lines are drawn at each value provided along the y-axis.
#' @export
#' @examples
#' # Generate some data. 1000 respondents taking 100 item test:
#' set.seed(060121)
#' p.success <- rBeta.4P(1000, 0.25, 0.75, 5, 3)
#' for (i in 1:100) {
#'   if (i == 1) {
#'    rawdata <- matrix(nrow = 1000, ncol = 100)
#'  }
#'  rawdata[, i] <- rbinom(1000, 1, p.success)
#' }
#'
#' # Analyse the accuracy and consistency of the test and store the object:
#' out <- LL.CA(x = rowSums(rawdata), reliability = cba(rawdata), cut = 50,
#' min = 0, max = 100, modelfit = c(nbins = 20, minbin = 1))
#'
#' # Feed the object to the mdlfit.gfx() function:
#' mdlfit.gfx(out)
#'
#' # Given the number of observations, the y-axis ticks are a bit crowded. We
#' # can make it look less crowded by changing the number of ticks, labels, and
#' # the grid-lines:
#' mdlfit.gfx(out, y.tickat = seq(0, 250, 25), y.lim = c(0, 250),
#' y.grid = seq(0, 250, 12.5))
mdlfit.gfx <- function(x, x.tickat = NULL, y.tickat = NULL, y.lim = NULL, main.lab = "Observed vs. Expected Frequencies",  x.lab = "Bins", y.lab = "Frequency", x.grid = NULL, y.grid = NULL) {
  tab <- x$modelfit$contingencytable
  bins <- seq_len(ncol(tab))
  if (is.null(y.lim)) y.lim <- c(0, max(tab) * 1.3)
  if (is.null(x.tickat)) x.tickat <- bins
  default_y_ticks <- seq(0, ceiling(y.lim[2]), ceiling(y.lim[2] / 10))
  if (is.null(y.tickat)) y.tickat <- default_y_ticks
  graphics::plot(NULL, xlim = range(bins), ylim = y.lim, ylab = y.lab, xlab = x.lab, main = main.lab, axes = FALSE)
  # NULL draws the default grid, NA draws none, and a vector draws lines at
  # the given positions.
  if (is.null(x.grid)) x.grid <- bins
  if (!is.na(x.grid[1])) graphics::abline(v = x.grid, col = "lightgrey", lty = 3)
  if (is.null(y.grid)) y.grid <- default_y_ticks
  if (!is.na(y.grid[1])) graphics::abline(h = y.grid, col = "lightgrey", lty = 3)
  graphics::lines(bins, tab[1, ], type = "o", col = "grey", lwd = 3, lty = 1, pch = 16, cex = 1.5)
  graphics::lines(bins, tab[2, ], type = "o", lwd = 3, pch = 1, cex = 1.5)
  graphics::box()
  graphics::axis(1, at = x.tickat, labels = FALSE)
  graphics::axis(2, at = y.tickat)
  graphics::legend("topright", legend = c("Observed", "Expected"), col = c("black", "darkgrey"),
                   lty = c(1, 1), lwd = c(3, 3), pch = c(1, 16), bty = "n", cex = 1.5)
  p_label <- if (x$modelfit$pvalue < 0.001) "p < 0.001" else paste("p = ", round(x$modelfit$pvalue, 3))
  graphics::legend("topleft", bty = "n", cex = 1.5,
                   legend = c(paste("chi-square = ", round(x$modelfit$chisquared, 2)),
                              paste("df = ", x$modelfit$df), p_label))
}



#' Function for estimating "Lord's k" for Lord's two-term approximation to the compound binomial distribution.
#'
#' @description Calculates Lord's k.
#' @param x A vector of observed-scores.
#' @param N The test length.
#' @param reliability The test-score reliability coefficient.
#' @return A value representing Lord's k
#' @export
#' @examples
#' # Generate some fictional data. Say 100 students take a 50-item long test
#' # where all items are equally difficult (i.e., where the true Lord's k = 0).
#' set.seed(1234)
#' p.success <- rBeta.4P(100, 0.25, 0.75, 5, 3)
#' for (i in 1:50) {
#'   if (i == 1) {
#'     rawdata <- matrix(nrow = 100, ncol = 50)
#'   }
#'   rawdata[, i] <- rbinom(100, 1, p.success)
#' }
#'
#' # Estimate the reliability of these scores with Cronbach's Alpha:
#' reliability <- cba(rawdata)
#'
#' # Estimate Lord's k using Lords.k():
#' Lords.k(rowSums(rawdata), 50, reliability)
Lords.k <- function(x, N, reliability) {
  mu <- mean(x)
  sigma2 <- stats::var(x)
  sigma2_true <- sigma2 * reliability  # observed minus error variance
  num <- N * ((N - 1) * sigma2_true - N * sigma2 + mu * (N - mu))
  den <- 2 * (mu * (N - mu) - sigma2_true)
  num / den
}


#' An Implementation of the Hanson and Brennan Approach to Estimate Classification Consistency and Accuracy based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Hanson and Brennan approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to a two-term approximation of the compound binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores, or a list specifying parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location-, and \code{alpha} and \code{beta} for the shape parameters of the Beta true-score distribution, and \code{k} for the "Lord's k" parameter (see documentation for the \code{Lords.k} function).
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param testlength The total number of test items (or maximum possible score). Must be an integer.
#' @param cut The cutoff value for classifying observations into above/below categories. Observed sum scores at or above \code{cut} are classified as above the cut.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param truecut Optional specification of a "true" cutoff. Useful for producing ROC curves (see documentation for the \code{HB.ROC()} function).
#' @param output Character vector indicating which types of statistics (i.e, accuracy and/or consistency) are to be computed and included in the output. Permissible values are \code{"accuracy"} and \code{"consistency"}.
#' @param failsafe Logical value indicating whether to engage the automatic fail-safe defaulting to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates. Default is \code{TRUE} (i.e., the function will engage failsafe if the four-parameter Beta-distribution fitting-procedure produced impermissible estimates).
#' @param l If \code{true.model = "2P"} or \code{failsafe = TRUE}, the lower-bound location parameter to be used in the two-parameter fitting procedure. Default is 0 (i.e., the lower-bound of the Standard Beta distribution).
#' @param u If \code{true.model = "2P"} or \code{failsafe = TRUE}, the upper-bound location parameter to be used in the two-parameter fitting procedure. Default is 1 (i.e., the upper-bound of the Standard Beta distribution).
#' @param modelfit Allows for controlling the chi-square test for model fit by setting the minimum bin-size for expected observations. Can alternatively be set to \code{NULL} to forego model-fit testing (speeding up the function). In accordance with standard recommendations for chi-square tests the default input to this argument is 10. Expected and observed frequencies are compared at each possible sum score, merging adjacent scores only where the expected frequency falls below this minimum. Degrees of freedom are the number of bins minus one, minus the number of moments matched when fitting the true-score distribution (two for a two-parameter solution, four otherwise); no degree of freedom is deducted for Lord's k. With \code{modelfit = 1} this reproduces the procedure of Lord (1965, pp. 244, 249), which yields \eqn{n - 4} degrees of freedom for an ungrouped four-parameter fit.
#' @note This implementation of the Hanson-Brennan approach is much slower than the implementation of the Livingston and Lewis approach, as there is no native implementation of Lord's two-term approximation to the Compound-Binomial distribution in R. This implementation uses a "brute-force" method of computing the cumulative probabilities from the compound-Binomial distribution, which will by necessity be more resource intensive.
#' @return A list containing the estimated parameters necessary for the approach (i.e., the effective test-length and the beta distribution parameters), a chi-square test of model-fit, the confusion matrix containing estimated proportions of true/false pass/fail categorizations for a test, diagnostic performance statistics, and / or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 50.
#' # Generate some fictional data. Say, 1000 individuals take a 20-item test.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
#'  for (i in 1:20) {
#'    if (i == 1) {
#'      rawdata <- matrix(nrow = 1000, ncol = 20)
#'      }
#'    rawdata[, i] <- rbinom(1000, 1, p.success)
#'  }
#'
#' # Suppose the cutoff value for attaining a pass is 10 items correct, and
#' # that the reliability of this test was estimated using the Cronbach's Alpha
#' # estimator. To estimate and retrieve the estimated parameters, confusion and
#' # consistency matrices, and accuracy and consistency indices using HB.CA():
#' HB.CA(x = rowSums(rawdata), reliability = cba(rawdata), cut = 10,
#' testlength = 20)
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @references Lord. Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Lewis, Don and Burke, C. J. (1949). The Use and Misuse of the Chi-Square Test. Psychological Bulletin, 46(6).
#' @export
HB.CA <- function(x = NULL, reliability, cut, testlength, true.model = "4P", truecut = NULL, output = c("accuracy", "consistency"), failsafe = TRUE, l = 0, u = 1, modelfit = 10) {
  model <- .true_model(true.model)
  from_scores <- !is.list(x)
  if (from_scores) {
    fit <- .hb_fit(x, testlength, Lords.k(x, testlength, reliability), model, failsafe && model == "4P", l, u)
    params <- fit$params
    .warn_out_of_bounds(params)
  } else {
    params <- x
  }
  if (is.null(truecut)) {
    truecut <- cut
  }
  N <- params$N
  k <- params$k
  truecut <- truecut / N
  out <- list(parameters = params)

  if (from_scores && !is.null(modelfit)) {
    out[["modelfit"]] <- .hb_modelfit(x, params, modelfit,
                                      n_params = if (model == "2P" || fit$engaged) 2 else 4)
  }

  # Observed scores at or above `cut` are classified as above the cut.
  below <- function(p) .pcbinom_tail(cut, N, k, p, lower.tail = TRUE)
  above <- function(p) .pcbinom_tail(cut, N, k, p, lower.tail = FALSE)
  integrate_tp <- function(f, lower = 0, upper = 1) {
    .integrate_support(function(p) dBeta.4P(p, params$l, params$u, params$alpha, params$beta) * f(p),
                       params$l, params$u, lower, upper)
  }

  output <- tolower(output)
  if (any(startsWith(output, "a"))) {
    p.tp <- integrate_tp(below, 0, truecut)
    p.fp <- integrate_tp(below, truecut, 1)
    p.tn <- integrate_tp(above, truecut, 1)
    p.fn <- integrate_tp(above, 0, truecut)
    camat <- confmat(p.tp, p.tn, p.fp, p.fn, "freq")
    out[["confusionmatrix"]] <- camat
    out[["classification.accuracy"]] <- caStats(camat[1, 1], camat[1, 2], camat[2, 1], camat[2, 2])
  }
  if (any(startsWith(output, "c"))) {
    p.ij <- integrate_tp(function(p) below(p) * above(p))
    out <- c(out, .binary_consistency(ii = integrate_tp(function(p) below(p)^2), ij = p.ij, ji = p.ij,
                                      jj = integrate_tp(function(p) above(p)^2)))
  }
  out
}

# Chi-square model fit for the Hanson and Brennan approach: one bin per
# possible sum score before merging.
.hb_modelfit <- function(x, params, minbin, n_params) {
  N <- params$N
  expected <- vapply(0:N, function(i) {
    .integrate_support(function(p) dBeta.4P(p, params$l, params$u, params$alpha, params$beta) * dcBinom(i, N, params$k, p),
                       params$l, params$u, 0, 1)
  }, numeric(1)) * length(x)
  observed <- vapply(0:N, function(i) sum(x == i), numeric(1))
  .chisq_modelfit(expected, observed, minbin, n_params, rownames = c("Expected", "Observed"))
}
#' ROC curves for the Hanson and Brennan approach.
#'
#' @description Generate a ROC curve plotting the false-positive rate against the true-positive rate at different cut-off values across the observed-score scale.
#' @param x A vector of observed results (sum scores) or a list of parameter values (see documentation for the \code{HB.beta.tp.fit() function}.
#' @param reliability The reliability coefficient of the test.
#' @param testlength The total number of test items (or maximum possible score). Must be an integer.
#' @param truecut The point along the x-scale that marks true category membership.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The \code{"4P"} method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the \code{"2P"} method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param failsafe If true-model == "4P": Whether to engage a fail-safe reverting to a two-parameter true-score distribution solution should the four-parameter fitting procedure produce impermissible results. Default is TRUE (engage fail-safe in the event of impermissible estimates).
#' @param l If \code{true.model == "2P"} or \code{failsafe == TRUE}: The lower-bound location parameter of the two-parameter true-score distribution solution.
#' @param u If \code{true.model == "2P"} or \code{failsafe == TRUE}: The upper-bound location parameter of the two-parameter true-score distribution solution.
#' @param AUC Logical. Calculate and include the area under the curve? Default is \code{FALSE}.
#' @param maxJ Logical. Mark the point along the curve where Youden's J statistic is maximized? Default is \code{FALSE}.
#' @param maxAcc Logical. Mark the point along the curve where the Accuracy statistic is maximized? Default is \code{FALSE}.
#' @param locate Ask the function to locate the cut-point at which sensitivity or NPV is greater than or equal to some value, or specificity or PPV is lesser than or equal to some value. Take as input a character-vector of length 2, with the first argument being which index is to be found (e.g., "sensitivity"), and the second argument the value to locate (e.g., "0.75"). For example: c("sensitivity", "0.75").
#' @param raw.out Give raw coordinates as output rather than plot? Default is \code{FALSE}
#' @param grainsize Specify the number of cutoff-points for which the ROC curve is to be calculated. The greater this number the greater the accuracy. Default is set to the stated test length (N).
#' @note This implementation of the Hanson-Brennan approach is much slower than the implementation of the Livingston and Lewis approach, as there is no native implementation of Lord's two-term approximation to the Compound-Binomial distribution in R. This implementation uses a "brute-force" method of computing the cumulative probabilities from the compound-Binomial distribution, which will by necessity be more resource intensive.
#' @return A plot tracing the ROC curve for the test, or matrix of coordinates if raw.out is \code{TRUE}.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 50.
#' # Generate some fictional data. Say, 1000 individuals take a 20-item test.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
#'  for (i in 1:20) {
#'    if (i == 1) {
#'      rawdata <- matrix(nrow = 1000, ncol = 20)
#'      }
#'    rawdata[, i] <- rbinom(1000, 1, p.success)
#'  }
#'
#' # Suppose the cutoff value for attaining a pass is 10 items correct, and
#' # that the reliability of this test was estimated using the Cronbach's Alpha
#' # estimator. To draw the ROC-graph and locate the points at which Youden's J
#' # and Accuracy are maximized:
#' HB.ROC(rowSums(rawdata), cba(rawdata), 20, 10, maxAcc = TRUE, maxJ = TRUE)
#'
#' # For further examples regarding how to use the locate argument to locate
#' # points at which various criteria are satisfied, see documentation for the
#' # LL.ROC() function.
#' @export
HB.ROC <- function(x = NULL, reliability, testlength, truecut, true.model = "4P", failsafe = TRUE, l = 0, u = 1, AUC = FALSE, maxJ = FALSE, maxAcc = FALSE, locate = NULL, raw.out = FALSE, grainsize = testlength) {
  if (!is.list(x)) {
    x <- HB.beta.tp.fit(x, testlength, Lords.k(x, testlength, reliability), true.model = true.model,
                        failsafe = failsafe, l = l, u = u)
  }
  rocmat <- .roc_matrix(seq(0, testlength, testlength / grainsize), function(cut) {
    HB.CA(x = x, cut = cut, truecut = truecut, true.model = true.model,
          output = "a", l = l, u = u, modelfit = NULL)$classification.accuracy
  })
  if (raw.out) {
    return(rocmat)
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  .roc_plot(rocmat, truecut, AUC = AUC, maxJ = maxJ, maxAcc = maxAcc, locate = locate)
}



#' An Extension of the Hanson and Brennan Approach to Estimate Classification Consistency and Accuracy for Multiple Classifications based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Hanson and Brennan approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to the binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores, or a list specifying parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location-, and \code{alpha} and \code{beta} for the shape parameters of the Beta true-score distribution, and \code{k} for the "Lord's k" parameter (see documentation for the \code{Lords.k} function).
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param testlength The total number of test items (or maximum possible score). Must be an integer.
#' @param cut A vector of cut-off values for classifying observations into two or more categories.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param failsafe Logical value indicating whether to engage the automatic fail-safe defaulting to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates. Default is \code{TRUE} (i.e., the function will engage failsafe if the four-parameter Beta-distribution fitting-procedure produced impermissible estimates).
#' @param l If \code{true.model = "2P"} or \code{failsafe = TRUE}, the lower-bound location parameter to be used in the two-parameter fitting procedure. Default is 0 (i.e., the lower-bound of the Standard Beta distribution).
#' @param u If \code{true.model = "2P"} or \code{failsafe = TRUE}, the upper-bound location parameter to be used in the two-parameter fitting procedure. Default is 1 (i.e., the upper-bound of the Standard Beta distribution).
#' @param modelfit Allows for controlling the chi-square test for model fit by setting the minimum bin-size for expected observations. Can alternatively be set to \code{NULL} to forego model-fit testing (speeding up the function). In accordance with standard recommendations for chi-square tests the default input to this argument is 10. Expected and observed frequencies are compared at each possible sum score, merging adjacent scores only where the expected frequency falls below this minimum. Degrees of freedom are the number of bins minus one, minus the number of moments matched when fitting the true-score distribution (two for a two-parameter solution, four otherwise); no degree of freedom is deducted for Lord's k. With \code{modelfit = 1} this reproduces the procedure of Lord (1965, pp. 244, 249), which yields \eqn{n - 4} degrees of freedom for an ungrouped four-parameter fit.
#' @return A list containing the estimated parameters necessary for the approach (i.e., Lord's k, test-length, and the true-score Beta distribution parameters), a chi-square test of model-fit, the confusion matrix containing estimated proportions of true/false positive/negative categorizations for a test, diagnostic performance statistics, and/or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @note This implementation of the Hanson-Brennan approach is much slower than the implementation of the Livingston and Lewis approach, as there is no native implementation of Lord's two-term approximation to the Compound-Binomial distribution in R. This implementation uses a "brute-force" method of computing the cumulative probabilities from the compound-Binomial distribution, which will by necessity be more resource intensive.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a 20-item test.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
#'  for (i in 1:20) {
#'    if (i == 1) {
#'      rawdata <- matrix(nrow = 1000, ncol = 20)
#'      }
#'    rawdata[, i] <- rbinom(1000, 1, p.success)
#'  }
#'
#' # Suppose the cutoff value for attaining a pass is 10 items correct, and
#' # that the reliability of this test was estimated using the Cronbach's Alpha
#' # estimator. To estimate and retrieve the estimated parameters, confusion and
#' # consistency matrices, and accuracy and consistency indices using HB.CA():
#' (output <- HB.CA.MC(x = rowSums(rawdata), reliability = cba(rawdata),
#' cut = c(8, 12), testlength = 20))
#'
#' # The output for this function can get quite verbose as more categories are
#' # included. The output from the function can be fed to the MC.out.tabular()
#' # function in order to organize the output in a tabular format.
#' MC.out.tabular(output)
#'
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @references Lord. Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Lewis, Don and Burke, C. J. (1949). The Use and Misuse of the Chi-Square Test. Psychological Bulletin, 46(6).
#' @export
HB.CA.MC <- function(x = NULL, reliability, cut, testlength, true.model = "4P", failsafe = TRUE, l = 0, u = 1, modelfit = 10) {
  model <- .true_model(true.model)
  from_scores <- !is.list(x)
  if (from_scores) {
    fit <- .hb_fit(x, testlength, Lords.k(x, testlength, reliability), model, failsafe && model == "4P", l, u)
    params <- fit$params
    .warn_out_of_bounds(params)
  } else {
    params <- x
  }
  N <- params$N
  k <- params$k
  out <- list(parameters = params)

  if (from_scores && !is.null(modelfit)) {
    out[["modelfit"]] <- .hb_modelfit(x, params, modelfit,
                                      n_params = if (model == "2P" || fit$engaged) 2 else 4)
  }

  truecut <- c(0, cut / N, 1)
  bounds <- c(0, cut)
  K <- length(cut) + 1
  # Probability that the observed score falls in category j given true score p.
  category_prob <- function(p, j) {
    upper <- if (j == K) 1 else .pcbinom_tail(bounds[j + 1], N, k, p)
    upper - .pcbinom_tail(bounds[j], N, k, p)
  }
  integrate_tp <- function(f, lower = 0, upper = 1) {
    .integrate_support(function(p) dBeta.4P(p, params$l, params$u, params$alpha, params$beta) * f(p),
                       params$l, params$u, lower, upper)
  }

  # Accuracy: rows are observed categories, columns are true categories.
  camat <- matrix(nrow = K, ncol = K,
                  dimnames = list(.category_labels(cut, "Observed"), .category_labels(cut, "True    ")))
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      camat[j, i] <- integrate_tp(function(p) category_prob(p, j), truecut[i], truecut[i + 1])
    }
  }
  out[["accuracy"]] <- .mc_accuracy(camat / sum(camat))

  # Consistency: rows are the second administration, columns the first.
  ccmat <- matrix(nrow = K, ncol = K,
                  dimnames = list(.category_labels(cut, "2nd adm."), .category_labels(cut, "1st adm.")))
  for (i in seq_len(K)) {
    for (j in i:K) {
      ccmat[j, i] <- ccmat[i, j] <- integrate_tp(function(p) category_prob(p, j) * category_prob(p, i))
    }
  }
  out[["consistency"]] <- .mc_consistency(ccmat / sum(ccmat))
  out
}
