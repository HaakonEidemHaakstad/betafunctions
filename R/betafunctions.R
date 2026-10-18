#' Compute Moments of Two-to-Four Parameter Beta Probability Density Distributions.
#'
#' @description Computes Raw, Central, or Standardized moment properties of defined Beta probability density distributions.
#' @param alpha The alpha shape parameter.
#' @param beta The beta shape parameter.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param types A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To compute the first four raw, central, and standardized
#' # moments of this distribution using betamoments():
#' betamoments(alpha = 5, beta = 3, l = 0.25, u = 0.75,
#' types = c("raw", "central", "standardized"), orders = 4)
#' @references Hanson, B. A (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series.
#' @return A list of moment types, each a list of moment orders.
#' @export
betamoments <- function(alpha, beta, l = 0, u = 1, types = c("raw", "central", "standardized"), orders = 4) {
  # Closed-form moments. With Y ~ Beta(alpha, beta) on [0, 1] and
  # X = l + (u - l) * Y, E[Y^k] = prod_{j = 0}^{k - 1} (alpha + j) / (alpha + beta + j),
  # and the moments of X follow from binomial expansions.
  i <- seq_len(orders)
  j <- seq_len(orders) - 1
  ey <- c(1, cumprod((alpha + j) / (alpha + beta + j)))  # E[Y^0], ..., E[Y^orders]
  width <- u - l
  central_y <- function(r) {
    k <- 0:r
    sum(choose(r, k) * ey[k + 1] * (-ey[2])^(r - k))
  }
  .moment_list(
    types,
    raw = function() vapply(i, function(r) {
      k <- 0:r
      sum(choose(r, k) * l^(r - k) * width^k * ey[k + 1])
    }, numeric(1)),
    central = function() width^i * vapply(i, central_y, numeric(1)),
    standardized = function() vapply(i, central_y, numeric(1)) / sqrt(central_y(2))^i
  )
}

#' Compute Mode of Two- and Four-Parameter Beta Probability Density distribution.
#'
#' @description Computes the mode of a Beta distribution with specified shape- and location parameters.
#' @param alpha The alpha shape parameter of the Probability Density Distribution.
#' @param beta The beta shape parameter of the Probability Density Distribution.
#' @param l The first (lower) location parameter of a four-parameter distribution. Default set to \code{0}.
#' @param u The second (upper) location parameter of a four-parameter distribution. Default set to \code{1}.
#' @examples
#' # To calculate the mode of a two-parameter (standard) Beta distribution with
#' # shape parameters alpha = 5 and beta = 3:
#' betamode(alpha = 5, beta = 3)
#'
#' # To calculate the mode of a four-parameter Beta distribution with shape
#' # parameters alpha = 5 and beta = 3, and location parameters l = 25 and
#' # u = 150:
#' betamode(alpha = 5, beta = 3, l = 25, u = 150)
#' @export
betamode <- function(alpha, beta, l = 0, u = 1) {
  if (any(alpha <= 1 | beta <= 1)) {
    warning("The mode is only uniquely defined (as computed here) for alpha > 1 and beta > 1.")
  }
  (alpha - 1) / (alpha + beta - 2) * (u - l) + l
}

#' Compute Median of Two- and Four-Parameter Beta Probability Density distribution.
#'
#' @description Computes the median of a Beta distribution with specified shape- and location parameters.
#' @note Versions prior to 1.9.1 returned the closed-form approximation \code{(alpha - 1/3) / (alpha + beta - 2/3)}, which is only accurate for \code{alpha, beta >= 1}. The exact median is now returned.
#' @param alpha The alpha shape parameter.
#' @param beta The beta shape parameter.
#' @param l The first (lower) location parameter. Default set to \code{0}.
#' @param u The second (upper) location parameter. Default set to \code{1}.
#' @examples
#' # To calculate the median of a two-parameter (standard) Beta distribution with
#' # shape parameters alpha = 5 and beta = 3:
#' betamedian(alpha = 5, beta = 3)
#'
#' # To calculate the median of a four-parameter Beta distribution with shape
#' # parameters alpha = 5 and beta = 3, and location parameters l = 25 and
#' # u = 150:
#' betamedian(alpha = 5, beta = 3, l = 25, u = 150)
#' @export
betamedian <- function(alpha, beta, l = 0, u = 1) {
  stats::qbeta(0.5, alpha, beta) * (u - l) + l
}

#' Compute Moments of Binomial Probability Mass Functions.
#'
#' @description Computes Raw, Central, or Standardized moment properties of defined Binomial probability mass functions.
#' @param n Number of Binomial trials
#' @param p Probability of success per trial.
#' @param types A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @examples
#' # Assume some variable follows a Binomial distribution with number of trials
#' # equal to 100 and a probability of success on each trial of 0.75. To compute
#' # the first four raw, central, and standardized moments of this distribution
#' # using binomialmoments():
#' binomialmoments(n = 100, p = 0.75, types = c("raw", "central",
#' "standardized"), orders = 4)
#'
#' # To only compute the (e.g.) standardized moments:
#' binomialmoments(n = 100, p = 0.75, types = "standardized")
#'
#' # To compute moments beyond the fourth order (e.g., the sixth):
#' binomialmoments(n = 100, p = 0.75, orders = 6)
#' @return A list of moment types, each a list of moment orders.
#' @export
binomialmoments <- function(n, p, types = c("raw", "central", "standardized"), orders = 4) {
  .pmf_moments(0:n, stats::dbinom(0:n, n, p), types, orders)
}

#' Compute Moments of Beta-Binomial Probability Mass Functions.
#'
#' @description Computes Raw, Central, or Standardized moment properties of defined Beta-Binomial probability mass functions.
#' @param N Number of trials.
#' @param l The first (lower) location-parameter of the Beta distribution.
#' @param u The second (upper) location-parameter of the Beta distribution.
#' @param alpha The alpha (first) shape-parameter of the Beta distribution.
#' @param beta The beta (second) shape-parameter of the Beta-distribution.
#' @param types A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @examples
#' # Assume 100 observations of a discrete variable with probabilities of
#' # positive outcomes adhering to a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = .95, and shape parameters a = 5 and
#' # b = 3. To compute the first four raw, central, and standardized moments of
#' # this distrubution using betabinomialmoments():
#' betabinomialmoments(N = 100, l = .25, u = .95, alpha = 5, beta = 3,
#' types = c("raw", "central", "standardized"), orders = 4)
#' @references Hanson, B. A (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series.
#' @return A list of moment types, each a list of moment orders.
#' @export
betabinomialmoments <- function(N, l, u, alpha, beta, types = c("raw", "central", "standardized"), orders = 4) {
  if (l < 0 || u > 1) {
    warning("Beta-distribution location parameter(s) out of bounds (l < 0 or u > 1).")
  }
  .pmf_moments(0:N, .betabinom_pmf(0:N, N, l, u, alpha, beta), types, orders)
}

#' Compute Moments of Observed Value Distribution.
#'
#' @description Computes Raw, Central, or Standardized moment properties of a vector of observed scores.
#' @param x A vector of values, the distribution of which moments are to be calculated.
#' @param type A character vector determining which moment-types are to be calculated. Permissible values are \code{"raw"}, \code{"central"}, and \code{"standardized"}.
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @param correct Logical. Whether to include bias correction in estimation of orders. Default is \code{TRUE}. For central moments, \code{TRUE} divides the sums of deviations by \code{n - 1} rather than \code{n}. For standardized moments, \code{TRUE} scales deviations by the sample standard deviation (\code{sd()}, denominator \code{n - 1}) and \code{FALSE} by the population standard deviation (denominator \code{n}). (Prior to version 1.9.1, \code{correct} had no effect on the standardized moments.)
#' @return A list of moment types, each a list of moment orders.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # To compute the first four raw, central, and standardized moments for this
#' # distribution of observed scores using observedmoments():
#' observedmoments(x = testdata, type = c("raw", "central", "standardized"),
#' orders = 4, correct = TRUE)
#' @export
observedmoments <- function(x, type = c("raw", "central", "standardized"),  orders = 4, correct = TRUE) {
  x <- x[!is.na(x)]
  n <- length(x)
  i <- seq_len(orders)
  dev <- x - mean(x)
  central_sums <- vapply(i, function(r) sum(dev^r), numeric(1))
  # With correct = TRUE the deviations are scaled by the (n - 1) standard
  # deviation; with correct = FALSE by the (n) population standard deviation.
  scale <- if (correct) stats::sd(x) else sqrt(sum(dev^2) / n)
  .moment_list(
    type,
    raw = function() vapply(i, function(r) sum(x^r) / n, numeric(1)),
    central = function() central_sums / (if (correct) n - 1 else n),
    standardized = function() central_sums / n / scale^i
  )
}


#' Alpha Shape-Parameter Given Location-Parameters, Mean, and Variance a Four-Parameter Beta Probability Density Distribution.
#'
#' @description Calculates the Beta value required to produce a Beta probability density distribution with defined moments and parameters. Be advised that not all combinations of moments and parameters can be satisfied (e.g., specifying mean, variance, skewness and kurtosis uniquely determines both location-parameters, meaning that the value of the lower-location parameter will take on which ever value it must, and cannot be specified).
#' @param mean The mean (first raw moment) of the target Standard Beta probability density distribution.
#' @param variance The variance (second central moment) of the target Standard Beta probability density distribution.
#' @param l The lower-bound location parameter of the Beta distribution. Default is 0 (as it is for the Standard Beta distribution).
#' @param u The upper-bound location parameter of the Beta distribution. Default is 1 (as it is for the Standard Beta distribution).
#' @param sd Optional alternative to specifying \code{var}. The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Alpha shape-parameter in order to produce a  Beta probability density distribution with the target mean and variance, given specified lower- and upper bounds of the Beta distribution.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0, rescaled to proportion
#' # of maximum.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3)) / 100
#' hist(testdata, xlim = c(0, 1))
#'
#' # To find the alpha shape-parameter of a Standard (two-parameter) Beta
#' # distribution with the same mean and variance as the observed-score
#' # distribution using AMS():
#' AMS(mean(testdata), var(testdata))
#' @export
AMS <- function(mean, variance, l = 0, u = 1, sd = NULL) {
  if (!is.null(sd)) {
    variance <- sd^2
  }
  # Method of moments on the standardized [0, 1] scale:
  # alpha = m * (m * (1 - m) / s2 - 1), with m and s2 the standardized mean
  # and variance.
  m <- (mean - l) / (u - l)
  s2 <- variance / (u - l)^2
  m * (m * (1 - m) / s2 - 1)
}

#' Beta Shape-Parameter Given Location-Parameters, Mean, and Variance of a Four-Parameter Beta Probability Density Distribution.
#'
#' @description Calculates the Beta value required to produce a Beta probability density distribution with defined moments and parameters. Be advised that not all combinations of moments and parameters can be satisfied (e.g., specifying mean, variance, skewness and kurtosis uniquely determines both location-parameters, meaning that the value of the lower-location parameter will take on which ever value it must, and cannot be specified).
#' @param mean The mean (first raw moment) of the target Standard Beta probability density distribution.
#' @param variance The variance (second central moment) of the target Standard Beta probability density distribution.
#' @param l The lower-bound location parameter of the Beta distribution. Default is 0 (as it is for the Standard Beta distribution).
#' @param u The upper-bound location parameter of the Beta distribution. Default is 1 (as it is for the Standard Beta distribution).
#' @param sd Optional alternative to specifying \code{var}. The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance, given specified lower- and upper bounds of the Beta distribution.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0, rescaled to proportion
#' # of maximum.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3)) / 100
#' hist(testdata, xlim = c(0, 1))
#'
#' # To find the beta shape-parameter of a Standard (two-parameter) Beta
#' # distribution with the same mean and variance as the observed-score
#' # distribution using BMS():
#' BMS(mean(testdata), var(testdata))
#'
#' # To find the beta shape-parameter of a four-parameter Beta
#' # distribution with specified lower- and upper-bounds of l = 0.25 and
#' # u = 0.75 using BMS:
#' BMS(mean(testdata), var(testdata), 0.25, 0.75)
#' @export
BMS <- function(mean, variance, l = 0, u = 1, sd = NULL) {
  if (!is.null(sd)) {
    variance <- sd^2
  }
  # Method of moments on the standardized [0, 1] scale:
  # beta = (1 - m) * (m * (1 - m) / s2 - 1).
  m <- (mean - l) / (u - l)
  s2 <- variance / (u - l)^2
  (1 - m) * (m * (1 - m) / s2 - 1)
}

#' Lower Location Parameter Given Shape Parameters, Mean, Variance, and Upper Location Parameter of a Four-Parameter Beta Probability Density Distribution.
#'
#' @description Calculates the lower-bound value required to produce a Beta probability density distribution with defined moments and parameters. Be advised that not all combinations of moments and parameters can be satisfied (e.g., specifying mean, variance, skewness and kurtosis uniquely determines both location-parameters, meaning that the value of the lower-location parameter will take on which ever value it must, and cannot be specified).
#' @param alpha The alpha (first) shape-parameter of the target Beta probability density distribution.
#' @param beta The beta (second) shape-parameter of the target Beta probability density distribution.
#' @param mean The mean (first raw moment) of the target Standard Beta probability density distribution.
#' @param variance The variance (second central moment) of the target Standard Beta probability density distribution.
#' @param u The upper-bound of the Beta distribution. Default is NULL (i.e., does not take a specified u-parameter into account). If supplied, \code{l} is determined by \code{alpha}, \code{beta}, \code{mean}, and \code{u} alone, and \code{variance} is not used.
#' @param sd Optional alternative to specifying \code{var}. The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta lower location-parameter (\code{l}) in order to produce a Beta probability density distribution with the target moments and parameters.
#' @examples
#' # Generate some fictional data.
#' set.seed(1234)
#' testdata <- rBeta.4P(100000, 0.25, 0.75, 5, 3)
#' hist(testdata, xlim = c(0, 1), freq = FALSE)
#'
#' # Suppose you know three of the four necessary parameters to fit a four-
#' # parameter Beta distribution (i. e., u = 0.75, alpha = 5, beta = 3) to this
#' # data. To find the value for the necessary l parameter, estimate the mean
#' # and variance of the distribution:
#' M <- mean(testdata)
#' S2 <- var(testdata)
#'
#' # To find the l parameter necessary to produce a four-parameter Beta
#' # distribution with the target mean, variance, and u, alpha, and beta
#' # parameters using the LABMSU() function:
#' (l <- LABMSU(alpha = 5, beta = 3, mean = M, variance = S2, u = 0.75))
#' curve(dBeta.4P(x, l, .75, 5, 3), add = TRUE, lwd = 2)
#' @export
LABMSU <- function(alpha = NULL, beta = NULL, u = NULL, mean = NULL, variance = NULL,  sd = NULL) {
  if (!is.null(sd)) {
    variance <- sd^2
  }
  if (is.null(alpha) || is.null(beta) || is.null(mean)) {
    warning("Not enough information.")
    return(NULL)
  }
  if (!is.null(u)) {
    # The mean alone determines l once u and the shape parameters are known:
    # mean = l + (u - l) * alpha / (alpha + beta).
    return((alpha * mean - alpha * u + beta * mean) / beta)
  }
  if (!is.null(variance)) {
    return(mean - alpha * sqrt(variance * (alpha + beta + 1)) / sqrt(alpha * beta))
  }
  warning("Not enough information.")
  NULL
}

#' Upper Location Parameter Given Shape Parameters, Mean, Variance, and Lower Location Parameter of a Four-Parameter Beta Probability Density Distribution.
#'
#' @description Calculates the upper-bound value required to produce a Beta probability density distribution with defined moments and parameters. Be advised that not all combinations of moments and parameters can be satisfied (e.g., specifying mean, variance, skewness and kurtosis uniquely determines both location-parameters, meaning that the value of the upper-location parameter will take on which ever value it must, and cannot be specified).
#' @param alpha The alpha shape-parameter of the target Beta probability density distribution.
#' @param beta The beta shape-parameter of the target Beta probability density distribution.
#' @param mean The mean (first raw moment) of the target Standard Beta probability density distribution.
#' @param variance The variance (second central moment) of the target Standard Beta probability density distribution.
#' @param l The lower-bound of the Beta distribution. Default is NULL (i.e., does not take a specified l-parameter into account). If supplied, \code{u} is determined by \code{alpha}, \code{beta}, \code{mean}, and \code{l} alone, and \code{variance} is not used.
#' @param sd Optional alternative to specifying \code{var}. The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta upper location-parameter (\code{u}) in order to produce a Beta probability density distribution with the target moments and parameters.
#' @examples
#' # Generate some fictional data.
#' set.seed(1234)
#' testdata <- rBeta.4P(100000, 0.25, 0.75, 5, 3)
#' hist(testdata, xlim = c(0, 1), freq = FALSE)
#'
#' # Suppose you know three of the four necessary parameters to fit a four-
#' # parameter Beta distribution (i. e., l = 0.25, alpha = 5, beta = 3) to this
#' # data. To find the value for the necessary u parameter, estimate the mean
#' # and variance of the distribution:
#' M <- mean(testdata)
#' S2 <- var(testdata)
#'
#' # To find the u parameter necessary to produce a four-parameter Beta
#' # distribution with the target mean, variance, and l, alpha, and beta
#' # parameters using the UABMSL() function:
#' (u <- UABMSL(alpha = 5, beta = 3, mean = M, variance = S2, l = 0.25))
#' curve(dBeta.4P(x, 0.25, u, 5, 3), add = TRUE, lwd = 2)
#' @export
UABMSL <- function(alpha = NULL, beta = NULL, mean = NULL, variance = NULL, l = NULL, sd = NULL) {
  if (!is.null(sd)) {
    variance <- sd^2
  }
  if (is.null(alpha) || is.null(beta) || is.null(mean)) {
    warning("Not enough information.")
    return(NULL)
  }
  if (!is.null(l)) {
    # The mean alone determines u once l and the shape parameters are known:
    # mean = l + (u - l) * alpha / (alpha + beta).
    return(beta * (mean - l) / alpha + mean)
  }
  if (!is.null(variance)) {
    return(mean + beta * sqrt(variance * (alpha + beta + 1)) / sqrt(alpha * beta))
  }
  warning("Not enough information.")
  NULL
}

#' Probability of Some Specific Observation under the Beta Probability Density Distribution with Specific Location Parameters, Mean, and Variance.
#'
#' @description Calculates the probability of some specific observation falling under a specified interval  ([0, x] or [x, 1]) under the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param q A specific point on the x-axis of the Standard Beta probability density distribution with a defined mean and variance.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param variance The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @param lower.tail Whether the density that should be considered is between the lower-end (i.e., [0 -> x]) or the higher-end of the distribution (i.e., [x -> 1]).
#' @param l The lower-bound location parameter. Default set to 0 (the standard Beta distribution).
#' @param u The upper-bound location parameter. Default set to 1 (the standard Beta distribution).
#' @return A value representing the probability of a random draw from the Standard Beta probability density distribution with a defined mean and variance being from one of two defined intervals (i.e., [0 -> x] or [x -> 1]).
#' @examples
#' # To compute the proportion of the density under the lower-end tail of a
#' # point along the Standard (two-parameter) Probability Density Distribution
#' # (e.g., 0.5) with mean of 0.6 and variance of 0.04:
#' pBetaMS(q = 0.5, mean = 0.6, variance = 0.04)
#'
#' # To compute the proportion of the density under the lower-end tail of a
#' # point along the Four-Parameter Beta Probability Density Distribution
#' # (e.g., 50) with mean of 60 and variance of 400, and lower-bound of 0 and
#' # upper-bound of 100:
#' pBetaMS(q = 50, mean = 60, variance = 400, l = 0, u = 100)
#' @export
pBetaMS <- function(q, mean, variance = NULL, sd = NULL, lower.tail = TRUE, l = 0, u = 1) {
  variance <- .resolve_variance(variance, sd)
  pBeta.4P(q, l, u, AMS(mean, variance, l, u), BMS(mean, variance, l, u), lower.tail = lower.tail)
}

#' Density Under a Specific Point of the Beta Probability Density Distribution with Specific Location Parameters, Mean, and Variance.
#'
#' @description Calculates the density under specific points of the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param x A specific point on the x-axis of the Standard Beta Probability Density Distribution.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param variance The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @param l The lower-bound location parameter. Default set to 0 (the standard Beta distribution).
#' @param u The upper-bound location parameter. Default set to 1 (the standard Beta distribution).
#' @return The probability density at \code{x} under the Beta distribution with the specified mean, variance, and location parameters.
#' @examples
#' # To compute the density at a specific point (e.g., 0.5) along the Standard
#' # (two-parameter) Probability Density Distribution with mean of 0.6 and variance of 0.04:
#' dBetaMS(x = 0.5, mean = 0.6, variance = 0.04)
#'
#' # To compute the density at a specific point (e.g., 50) along the four-
#' # parameter Beta distribution with a mean of 60, variance of 400, and lower-
#' # bound of 0 and upper-bound of 100:
#' dBetaMS(x = 50, mean = 60, variance = 400, l = 0, u = 100)
#' @export
dBetaMS <- function(x, mean, variance = NULL, sd = NULL, l = 0, u = 1) {
  variance <- .resolve_variance(variance, sd)
  dBeta.4P(x, l, u, AMS(mean, variance, l, u), BMS(mean, variance, l, u))
}

#' Quantile Containing Specific Proportion of the Distribution, Given a Specific Probability of the Beta Probability Density Distribution with Specific Mean and Variance.
#'
#' @description Calculates the quantile corresponding to a specific probability of some observation falling within the [0, x] (\code{lt = TRUE}) or [x, 1] (\code{lt = FALSE}) interval under the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param p A value of probability marking the point of the Y-axis to correspond to the X-axis.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param variance The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @param lower.tail Logical. Specifies which end of the tail for which to calculate quantile. Default is \code{TRUE} (meaning, find q for lower tail.)
#' @param l The lower-bound location parameter. Default set to 0 (the standard Beta distribution).
#' @param u The upper-bound location parameter. Default set to 1 (the standard Beta distribution).
#' @return A numeric value representing the quantile for which the specified proportion of observations fall within.
#' @examples
#' # To compute the quantile at a specific point (e.g., 0.5) along the Standard
#' # (two-parameter) Probability Density Distribution with mean of 0.6 and variance of 0.04:
#' qBetaMS(p = 0.5, mean = 0.6, variance = 0.04)
#'
#' # To compute the quantile at a specific points(e.g., 0.5) along the four-
#' # parameter Beta distribution with a mean of 60, variance of 400, and lower-
#' # bound of 0 and upper-bound of 100:
#' qBetaMS(p = 0.5, mean = 60, variance = 400, l = 0, u = 100)
#' @export
qBetaMS <- function(p, mean, variance = NULL, sd = NULL, lower.tail = TRUE, l = 0, u = 1) {
  variance <- .resolve_variance(variance, sd)
  qBeta.4P(p, l, u, AMS(mean, variance, l, u), BMS(mean, variance, l, u), lower.tail = lower.tail)
}

#' Random Draw from the Beta Probability Density Distribution With Specific Mean and Variance.
#'
#' @description Draws random samples of observations from the Standard Beta probability density distribution with defined mean and variance.
#' @param n Number of observations to be drawn from under the Standard Beta Probability Density Distribution.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param variance The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard probability density distribution.
#' @param l The lower-bound location parameter. Default set to 0 (the standard Beta distribution).
#' @param u The upper-bound location parameter. Default set to 1 (the standard Beta distribution).
#' @return A vector of length \code{n}, each value representing a random draw from the Standard Beta probability density distribution with defined mean and variance.
#' @examples
#' # To draw a random sample of 100 values from a Standard Beta distribution
#' # with a mean of 0.6 and variance = 0.04:
#' rBetaMS(n = 100, mean = 0.6, variance = 0.04)
#'
#' # To draw a random sample of 100 values from a Four-Parameter Beta
#' # distribution with a mean of 60 and variance = 400, and a lower-bound of
#' # 0 and an upper-bound of 100:
#' rBetaMS(n = 100, mean = 60, variance = 400, l = 0, u = 100)
#' @export
rBetaMS <- function(n, mean, variance = NULL, sd = NULL, l = 0, u = 1) {
  variance <- .resolve_variance(variance, sd)
  rBeta.4P(n, l, u, AMS(mean, variance, l, u), BMS(mean, variance, l, u))
}

#' Coordinate Generation for Marking an Area Under the Curve for the Beta Probability Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Beta probability density functions.
#' @param from The point of the \code{x}-axis from where to start producing \code{y}-density values.
#' @param to The point of the x-axis to where y-density values are to be produced.
#' @param by The resolution (or spacing) at which to produce y-density values.
#' @param alpha The alpha (first) shape-parameter value for the Standard Beta probability density distribution.
#' @param beta The beta (second) shape-parameter for the Standard Beta probability density distribution.
#' @param l The lower-bound location parameter of the Beta distribution.
#' @param u The upper-bound location parameter of the Beta distribution.
#' @return A two-column matrix with density-values of y to plot against corresponding location values of x.
#' @examples
#' # To box in an area under a four-parameter Beta distribution with location
#' # parameters l = .25 and u = .75, and shape parameters alpha = 5 and
#' # rbeta = 3, from 0.4 to 0.6:
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 7))
#' coords <- Beta.gfx.poly.pdf(from = 0.4, to = 0.6, by = 0.001, alpha = 5,
#' beta = 3, l = 0.25, u = 0.75)
#' polygon(coords)
#' @export
Beta.gfx.poly.pdf <- function(from, to, by, alpha, beta, l = 0, u = 1) {
  inner <- seq(from, to, by)
  x <- c(from, inner, to)
  y <- c(0, dBeta.4P(inner, l, u, alpha, beta), 0)
  cbind(x, y)
}

#' Coordinate Generation for Marking an Area Under the Curve for the Beta Quantile Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Beta probability quantile functions.
#' @param from The point of the \code{x}-axis from where to start producing \code{y}-quantile values.
#' @param to The point of the \code{x}-axis to where \code{y}-quantile values are to be produced.
#' @param by The resolution (or spacing) at which to produce \code{y}-density values.
#' @param alpha The alpha shape-parameter value for the Standard Beta probability distribution.
#' @param beta The beta shape-parameter for the Standard Beta probability distribution.
#' @param l The lower-bound location parameter of the Beta distribution.
#' @param u The upper-bound location parameter of the Beta distribution.
#' @return A two-column matrix with quantile-values of \code{y} to plot against corresponding location values of \code{x}.
#' @examples
#' # To box in an area under a four-parameter Beta quantile distribution with
#' # location parameters l = .25 and u = 75, and shape parameters alpha = 5 and
#' # beta = 3, from .4 to .6:
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 1))
#' coords <- Beta.gfx.poly.qdf(from = 0.4, to = 0.6, by = 0.001, alpha = 5,
#' beta = 3, l = 0.25, u = 0.75)
#' polygon(coords)
#' @export
Beta.gfx.poly.qdf <- function(from, to, by, alpha, beta, l = 0, u = 1) {
  inner <- seq(from, to, by)
  x <- c(from, inner, to)
  y <- c(l, qBeta.4P(inner, l, u, alpha, beta), l)
  cbind(x, y)
}

#' Coordinate Generation for Marking an Area Under the Curve for the Beta Cumulative Probability Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Beta cumulative probability functions.
#' @param from The point of the \code{x}-axis from where to start producing \code{y}-density values.
#' @param to The point of the \code{x}-axis to where \code{y}-density values are to be produced.
#' @param by The resolution (or spacing) at which to produce \code{y}-density values.
#' @param alpha The alpha shape-parameter value for the Standard Beta cumulative probability distribution.
#' @param beta The beta shape-parameter for the Standard Beta cumulative probability distribution.
#' @param l The lower-bound location parameter of the Beta distribution.
#' @param u The upper-bound location parameter of the Beta distribution.
#' @return A two-column matrix with cumulative probability-values of y to plot against corresponding location values of \code{x}.
#' @examples
#' # To box in an area under a four-parameter Beta cumulative distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters
#' # alpha = 5 and beta = 3, from 0.4 to 0.6:
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 1))
#' coords <- Beta.gfx.poly.cdf(from = 0.4, to = 0.6, by = 0.001, alpha = 5,
#' beta = 3, l = 0.25, u = 0.75)
#' polygon(coords)
#' @export
Beta.gfx.poly.cdf <- function(from, to, by, alpha, beta, l = 0, u = 1) {
  inner <- seq(from, to, by)
  x <- c(from, inner, to)
  y <- c(0, pBeta.4P(inner, l, u, alpha, beta), 0)
  cbind(x, y)
}

#' Most Likely True Alpha Value Given Observed Outcome.
#'
#' @description Given a fitted Standard (two-parameter) Beta Distribution, return the alpha shape-parameter value where the observed mean becomes the mode.
#' @param alpha Observed alpha-parameter value for fitted Standard Beta Probability Density Distribution.
#' @param beta Observed beta-parameter value for fitted Standard Beta Probability Density Distribution.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @return The Alpha shape-parameter value for the Standard Beta probability density distribution where the observed mean is the expected mode.
#' @examples
#' # Assuming a prior Standard (two-parameter) Beta distribution is fit, which
#' # yield an alpha parameter of 10 and a beta parameter of 8, calculate the
#' # true-alpha parameter most likely to have produced the observations:
#' MLA(a = 10, b = 8)
#' @export
MLA <- function(alpha, beta, x = NULL, n = NULL) {
  if (is.null(x) || is.null(n)) {
    # Use the mean of the fitted distribution as the observed proportion.
    n <- alpha + beta
    x <- alpha / n
  }
  x * (n - 2) + 1
}

#' Most Likely True Beta Value Given Observed Outcome.
#'
#' @description Assuming a prior standard (two-parameter) Beta Distribution, return the beta shape-parameter value where the observed mean becomes the mode.
#' @param alpha Observed alpha-parameter value for fitted Standard Beta Probability Density Distribution.
#' @param beta Observed beta-parameter value for fitted Standard Beta Probability Density Distribution.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @examples
#' # Assuming a prior Standard (two-parameter) Beta distribution is fit, which
#' # yield an alpha parameter of 10 and a beta parameter of 8, calculate the
#' # true-beta parameter most likely to have produced the observations:
#' MLB(a = 10, b = 8)
#' @return The Beta shape-parameter value for the Standard Beta probability density distribution where the observed mean is the expected mode.
#' @export
MLB <- function(alpha, beta, x = NULL, n = NULL) {
  if (is.null(x) || is.null(n)) {
    n <- alpha + beta
    x <- alpha / n
  }
  n - (x * (n - 2) + 1)
}

#' Most Likely Mean of the Standard Beta Probability Density Distribution, Given that the Observation is Considered the Most Likely Observation of the Standard Beta Probability Density Distribution (i.e., the mode).
#'
#' @description Assuming a prior Standard (two-parameter) Beta Distribution, returns the expected mean of the distribution under the assumption that the observed value is the most likely value of the distribution.
#' @param alpha Observed alpha value for fitted Standard Beta Probability Density Distribution.
#' @param beta Observed beta value for fitted Standard Beta Probability Density Distribution.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @examples
#' # Assuming a prior Standard (two-parameter) Beta distribution is fit, which
#' # yield an alpha parameter of 10 and a beta parameter of 8, calculate the
#' # true-mean most likely to have produced the observations:
#' MLM(a = 10, b = 8)
#' @return The expected mean of the Standard Beta probability density distribution, for which the observed mean is the most likely value.
#' @export
MLM <- function(alpha, beta, x = NULL, n = NULL) {
  if (is.null(x) || is.null(n)) {
    n <- alpha + beta
    x <- alpha / n
  }
  (x * (n - 2) + 1) / n
}

#' Probability Density under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Gives the density at desired values of \code{x} under the Four-Parameter Beta Probability Density Distribution.
#' @param x Value of \code{x}.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @return The value for the probability density at specified values of \code{x}.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To compute the probability density at a specific point of
#' # the distribution (e.g., 0.5) using dBeta.4P():
#' dBeta.4P(x = 0.5, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
dBeta.4P <- function(x, l, u, alpha, beta) {
  stats::dbeta((x - l) / (u - l), alpha, beta) / (u - l)
}

#' Random Number Generation under the Beta-Binomial Probability Mass Distribution.
#'
#' @param n Number of draws.
#' @param N Number of trials.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The alpha (first) shape parameter.
#' @param beta The beta (second) shape parameter.
#' @return A vector with length \code{n} of random values drawn from the Beta-Binomial Distribution.
#' @examples
#' # To draw a sample of 50 values from a Beta-Binomial distribution with
#' # number of trials = 100, and with success-probabilities drawn from a
#' # Four-Parameter Beta distribution with location parameters l = 0.25 and
#' # u = 0.95, and shape-parameters alpha = 5 and beta = 3:
#' rBetaBinom(n = 50, N = 100, l = 0.25, u = 0.95, alpha = 5, beta = 3)
#' @export
rBetaBinom <- function(n, N, l, u, alpha, beta) {
  sample(0:N, size = n, prob = .betabinom_pmf(0:N, N, l, u, alpha, beta), replace = TRUE)
}

#' Probability Mass under the Beta-Binomial Probability-Mass Distribution.
#'
#' @description Gives the density at \code{x} under the Beta-Binomial PMF.
#' @param x Value of \code{x} (a specific number of successes).
#' @param N The total number of trials.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @return The value for the probability mass at \code{x} given the specified Beta-Binomial distribution.
#' @examples
#' # Assume some variable follows a Beta-Binomial distribution with 100 number
#' # of trials, and with probabilities of successful trials drawn from a four-
#' # parameter Beta distribution with location parameters l = 0.25 and u = 0.75
#' # and shape parameters alpha = 5 and beta = 3. To compute the probability
#' # density at a specific point of the distribution (e.g., 50):
#' dBetaBinom(x = 50, N = 100, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
dBetaBinom <- function(x, N, l, u, alpha, beta) {
  .betabinom_pmf(x, N, l, u, alpha, beta)
}

#' Cumulative Probability Function under the Beta-Binomial Probability Distribution.
#'
#' @description Function for calculating the proportion of observations up to a specifiable quantile under the Beta-Binomial Probability Distribution.
#' @param q The quantile or a vector of quantiles for which the proportion is to be calculated.
#' @param N The total number of trials.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lower.tail Logical. If \code{TRUE} (default), probabilities are P[X < q], otherwise P[X >= q]. Note that this differs from base-R \code{pbinom()}, which gives P[X <= q].
#' @return A vector of proportions of observations falling under specified quantiles under the Beta-Binomial distribution.
#' @examples
#' # Assume some variable follows a Beta-Binomial distribution with number of
#' # trials = 50, and probabilities of successful trials are drawn from a four-
#' # parameter Beta distribution with location parameters l = 0.25 and u =
#' # 0.75, and shape parameters alpha = 5 and beta = 3. To compute the
#' # cumulative probability at a specific point of the distribution (e.g., 25):
#' pBetaBinom(q = 25, N = 50, l = .25, u = .75, alpha = 5, beta = 3)
#' @export
pBetaBinom <- function(q, N, l, u, alpha, beta, lower.tail = TRUE) {
  .discrete_tail(q, .betabinom_pmf(0:N, N, l, u, alpha, beta), lower.tail = lower.tail)
}

#' Random Number Generation under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Function for generating random numbers from a specified Four-Parameter Beta Distribution.
#' @param n Number of draws.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The alpha (first) shape parameter.
#' @param beta The beta (second) shape parameter.
#' @return A vector with length \code{n} of random values drawn from the Four-Parameter Beta Distribution.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To draw a random value from this distribution using
#' # rBeta.4P():
#' rBeta.4P(n = 1, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
rBeta.4P <- function(n, l, u, alpha, beta) {
  stats::rbeta(n, alpha, beta) * (u - l) + l
}

#' Cumulative Probability Function under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Function for calculating the proportion of observations up to a specifiable quantile under the Four-Parameter Beta Distribution.
#' @param q The quantile or a vector of quantiles for which the proportion is to be calculated.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lower.tail Whether the proportion to be calculated is to be under the lower or upper tail. Default is \code{TRUE} (lower tail).
#' @return A vector of proportions of observations falling under specified quantiles under the four-parameter Beta distribution.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To compute the cumulative probability at a specific point of
#' # the distribution (e.g., 0.5)
#' # using pBeta.4P():
#' pBeta.4P(q = 0.5, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
pBeta.4P <- function(q, l, u, alpha, beta, lower.tail = TRUE) {
  stats::pbeta((q - l) / (u - l), alpha, beta, lower.tail = lower.tail)
}

#' Quantile Given Probability Under the Four-Parameter Beta Distribution.
#'
#' @description Function for calculating the quantile (i.e., value of \code{x}) for a given proportion (i.e., the value of \code{y}) under the Four-Parameter Beta Distribution.
#' @param p A vector (or single value) of proportions or probabilities for which the corresponding value of \code{x} (i.e., the quantiles) are to be calculated.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lower.tail Logical. If \code{TRUE} (default), returns \code{x} such that P[X <= x] = \code{p}; otherwise \code{x} such that P[X > x] = \code{p}.
#' @note Prior to version 1.9.1, \code{lower.tail = FALSE} returned \code{1 - qbeta(p, alpha, beta)} on the [l, u] scale, which is only correct when \code{alpha == beta}.
#' @return A vector of quantiles for specified probabilities or proportions of observations under the four-parameter Beta distribution.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To compute the quantile at a specific point of the
#' # distribution (e.g., 0.5) using qBeta.4P():
#' qBeta.4P(p = 0.5, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
qBeta.4P <- function(p, l, u, alpha, beta, lower.tail = TRUE) {
  stats::qbeta(p, alpha, beta, lower.tail = lower.tail) * (u - l) + l
}

#' Method of Moment Estimates of Shape- and Location Parameters of the Four-Parameter Beta Distribution.
#'
#' @description An implementation of the method of moments estimation of four-parameter Beta distribution parameters presented by Hanson (1991). Given a vector of values, calculates the shape- and location parameters required to produce a four-parameter Beta distribution with the same mean, variance, skewness and kurtosis (i.e., the first four moments) as the observed-score distribution.
#' @param scores A vector of values to which the four-parameter Beta distribution is to be fitted.
#' @param mean If scores are not supplied: specification of the mean for the target four-parameter Beta distribution.
#' @param variance If scores are not supplied: specification of the variance for the target four-parameter Beta distribution.
#' @param skewness If scores are not supplied: specification of the skewness for the target four-parameter Beta distribution.
#' @param kurtosis If scores are not supplied: specification of the kurtosis for the target four-parameter Beta distribution.
#' @return A list of parameter-values required to produce a four-parameter Beta distribution with the same first four moments as the observed distribution.
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100), freq = FALSE)
#'
#' # To fit and retrieve the parameters for a four-parameter Beta distribution
#' # to the observed-score distribution using Beta.4p.fit():
#' (params.4p <- Beta.4p.fit(testdata))
#' curve(dBeta.4P(x, params.4p$l, params.4p$u, params.4p$alpha, params.4p$beta), add = TRUE)
#' @export
Beta.4p.fit <- function(scores, mean = NULL, variance = NULL, skewness = NULL, kurtosis = NULL) {
  if (!is.null(mean) && !is.null(variance) && !is.null(skewness) && !is.null(kurtosis)) {
    m1 <- mean
    s2 <- variance
    g3 <- skewness
    g4 <- kurtosis
  } else {
    if (missing(scores) || is.null(scores)) {
      stop("Supply either `scores`, or all of `mean`, `variance`, `skewness`, and `kurtosis`.")
    }
    m <- observedmoments(scores)
    m1 <- m$raw[[1]]
    s2 <- m$central[[2]]
    g3 <- m$standardized[[3]]
    g4 <- m$standardized[[4]]
  }
  # Hanson (1991), equations for the four-parameter Beta method of moments.
  r <- 6 * (g4 - g3^2 - 1) / (6 + 3 * g3^2 - 2 * g4)
  root <- sqrt(1 - 24 * (r + 1) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))
  larger <- r / 2 * (1 + root)
  smaller <- r / 2 * (1 - root)
  # Negative skew means alpha > beta.
  if (g3 < 0) {
    a <- larger
    b <- smaller
  } else {
    a <- smaller
    b <- larger
  }
  spread <- sqrt(s2 * (a + b + 1)) / sqrt(a * b)
  list("alpha" = a, "beta" = b, "l" = m1 - a * spread, "u" = m1 + b * spread)
}

#' Method of Moment Estimates of Shape-Parameters of the Two-Parameter (Standard) Beta Distribution.
#'
#' @description An implementation of the method of moments estimation of two-parameter Beta distribution parameters. Given a vector of values, calculates the shape parameters required to produce a two-parameter Beta distribution with the same mean and variance (i.e., the first two moments) as the observed-score distribution.
#' @param scores A vector of values to which the two-parameter Beta distribution is to be fitted. The values ought to fall within the [0, 1] interval.
#' @param mean The mean of the target Beta distribution. Alternative to feeding the function raw scores.
#' @param variance The variance of the target Beta distribution. Alternative to feeding the function raw scores.
#' @param l Optional specification of a lower-bound parameter of the Beta distribution. Default is 0 (i.e., the lower-bound of the Standard two-parameter Beta distribution).
#' @param u Optional specification of an upper-bound parameter of the Beta distribution. Default is 1 (i.e., the lower-bound of the Standard two-parameter Beta distribution).
#' @return A list of parameter-values required to produce a Standard two-parameter Beta distribution with the same first two moments as the observed distribution.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3)) / 100
#' hist(testdata, xlim = c(0, 1), freq = FALSE)
#'
#' # To fit and retrieve the parameters for a two-parameter Beta distribution
#' # to the observed-score distribution using Beta.2p.fit():
#' (params.2p <- Beta.2p.fit(testdata))
#' curve(dbeta(x, params.2p$alpha, params.2p$beta), add = TRUE)
#' @export
Beta.2p.fit <- function(scores = NULL, mean = NULL, variance = NULL, l = 0, u = 1) {
  if (!is.null(scores)) {
    if (max(scores) > u || min(scores) < l) {
      warning("Input values outside the range of the specified location parameters of the Beta distribution (i.e., there are values falling outside the [",
              l, ", ", u, "] interval).")
    }
    mean <- base::mean(scores)
    variance <- stats::var(scores)
  }
  list("alpha" = AMS(mean, variance, l, u), "beta" = BMS(mean, variance, l, u), "l" = l, "u" = u)
}

#' An implementation of the Beta-density Compound Cumulative Binomial Distribution.
#'
#' @description The Beta Compound Binomial distribution: The product of the four-parameter Beta probability density function and the binomial cumulative probability mass function. Used in the Livingston and Lewis approach to classification accuracy and consistency, the output can be interpreted as the population density of passing scores produced at "x" (a value of true-score).
#' @param x x-axis input for which \code{p} (proportion or probability) is to be computed.
#' @param l The lower-bound of the four-parameter Beta distribution.
#' @param u The upper-bound of the four-parameter Beta distribution.
#' @param alpha The alpha shape-parameter of the Beta distribution.
#' @param beta The beta shape-parameter of the Beta distribution.
#' @param n The number of trials for the Binomial distribution.
#' @param c The "true-cut" (proportion) of the Binomial distribution.
#' @param lower.tail Logical. Whether to compute the lower or upper tail of the Binomial distribution. Default is \code{FALSE} (i.e., upper tail).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @note The Binomial distribution cut-point is up-to but not including, unlike the standard behaviour of base-R pbinom() function.
#' @examples
#' # Given a four-parameter Beta distribution with parameters l = 0.25, u = 0.75,
#' # alpha = 5, and beta = 3, and a Binomial error distribution with number of
#' # trials (n) = 10 and a cutoff-point (c) at 50% correct (i.e., proportion correct
#' # of 0.5), the population density of passing scores produced at true-score
#' # (x) = 0 can be calculated as:
#' dBeta.pBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5)
#'
#' # Conversely, the density of failing scores produced at x can be calculated
#' # by passing the additional argument "lower.tail = TRUE" to the function.
#' # That is:
#' dBeta.pBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5,
#' lower.tail = TRUE)
#'
#' #By integration, the population proportion of (e.g.) passing scores in some
#' #region of the true-score distribution (e.g. between 0.25 and 0.5) can be
#' #calculated as:
#' integrate(function(x) { dBeta.pBinom(x, 0.25, .75, 5, 3, 10, 0.5) },
#' lower = 0.25, upper = 0.5)
#' @export
dBeta.pBinom <- function(x, l, u, alpha, beta, n, c, lower.tail = FALSE) {
  dBeta.4P(x, l, u, alpha, beta) * stats::pbinom(round(n * c) - 1, round(n), x, lower.tail = lower.tail)
}

#' An implementation of a Beta-density Compound Cumulative Gamma-Binomial Distribution.
#'
#' @description The Beta Compound Binomial distribution: The product of the four-parameter Beta probability density function and the binomial cumulative probability mass function. Used in the Livingston and Lewis approach to classification accuracy and consistency, the output can be interpreted as the population density of passing scores produced at "x" (a value of true-score).
#' @param x x-axis input for which \code{p} (proportion or probability) is to be computed.
#' @param l The lower-bound of the four-parameter Beta distribution.
#' @param u The upper-bound of the four-parameter Beta distribution.
#' @param alpha The alpha shape-parameter of the four-parameter Beta distribution.
#' @param beta The beta shape-parameter of the four-parameter Beta distribution.
#' @param n The number of "trials" for the Gamma-Binomial distribution.
#' @param c The "true-cut" (proportion) on the Gamma-Binomial distribution. Need not be an integer (unlike Binomial distribution).
#' @param lower.tail Logical. Whether to compute the lower or upper tail of the Binomial distribution. Default is \code{FALSE} (i.e., upper tail).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @export
#' @examples
#' # Given a four-parameter Beta distribution with parameters l = 0.25, u = 0.75,
#' # alpha = 5, and beta = 3, and a Binomial error distribution with number of
#' # trials (n) = 10 and a cutoff-point (c) at 50% correct (i.e., proportion correct
#' # of 0.5), the population density of passing scores produced at true-score
#' # (x) = 0 can be calculated as:
#' dBeta.pGammaBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5)
#'
#' # Conversely, the density of failing scores produced at x can be calculated
#' # by passing the additional argument "lower.tail = TRUE" to the function.
#' # That is:
#' dBeta.pGammaBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10.1, c = 0.5,
#' lower.tail = TRUE)
#'
#' #By integration, the population proportion of (e.g.) passing scores in some
#' #region of the true-score distribution (e.g. between 0.25 and 0.5) can be
#' #calculated as:
#' integrate(function(x) { dBeta.pGammaBinom(x, 0.25, 0.75, 5, 3, 10, 0.5) },
#' lower = 0.25, upper = 0.5)
dBeta.pGammaBinom <- function(x, l, u, alpha, beta, n, c, lower.tail = FALSE) {
  density <- dBeta.4P(x, l, u, alpha, beta)
  out <- density
  # pGammaBinom() integrates numerically, so only evaluate it where the
  # density is non-zero.
  inside <- !is.na(density) & density != 0
  out[inside] <- density[inside] * pGammaBinom(n * c, n, x[inside], lower.tail = lower.tail)
  out
}

#' An implementation of a Beta-density Compound Cumulative-Beta Distribution.
#'
#' @description The Beta Compound Beta distribution: The product of the four-parameter Beta probability density function and the Beta cumulative probability function. Used in the Livingston and Lewis approach to classification accuracy and consistency, the output can be interpreted as the population density of passing scores produced at "x" (a value of true-score).
#' @param x x-axis input for which \code{p} (proportion or probability) is to be computed.
#' @param l The lower-bound of the four-parameter Beta distribution.
#' @param u The upper-bound of the four-parameter Beta distribution.
#' @param alpha The alpha shape-parameter of the Beta density distribution.
#' @param beta The beta shape-parameter of the Beta density distribution.
#' @param n The number of trials for the Beta cumulative probability distribution.
#' @param c The "true-cut" (proportion) of on the Beta cumulative probability distribution.
#' @param lower.tail Logical. Whether to compute the lower or upper tail of the Beta cumulative probability distribution. Default is \code{FALSE} (i.e., upper tail).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @examples
#' # Given a four-parameter Beta distribution with parameters l = 0.25, u = 0.75,
#' # alpha = 5, and beta = 3, and a Beta error distribution with number of
#' # trials (n) = 10 and a cutoff-point (c) at 50% correct (i.e., proportion correct
#' # of 0.5), the population density of passing scores produced at true-score
#' # (x) = 0.5 can be calculated as:
#' dBeta.pBeta(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5)
#'
#' # Conversely, the density of failing scores produced at x can be calculated
#' # by passing the additional argument "lower.tail = TRUE" to the function.
#' # That is:
#' dBeta.pBeta(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5,
#' lower.tail = TRUE)
#'
#' # By integration, the population proportion of (e.g.) passing scores in some
#' # region of the true-score distribution (e.g. between 0.25 and 0.5) can be
#' # calculated as:
#' integrate(function(x) { dBeta.pBeta(x, 0.25, 0.75, 5, 3, 10, 0.5) },
#' lower = 0.25, upper = 0.5)
#' @export
dBeta.pBeta <- function(x, l, u, alpha, beta, n, c, lower.tail = FALSE) {
  dBeta.4P(x, l, u, alpha, beta) * stats::pbeta(c, x * n, (1 - x) * n, lower.tail = lower.tail)
}

#' Gamma-extended Binomial coefficient (choose function).
#'
#' @description Extends the Binomial coefficient for positive non-integers (including 0) by employing the Gamma rather than the factorial function.
#' @param n In Binomial terms, the number of Binomial "trials". Need not be an integer.
#' @param k In Binomial terms, the number of successful "trials". Need not be an integer.
#' @note Not defined for negative integers.
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @examples
#' # Compare choose function with gchoose function for integers:
#' gchoose(c(8, 9, 10), c(3, 4, 5)) == choose(c(8, 9, 10), c(3, 4, 5))
#'
#' # The gchoose function also works for non-integers:
#' gchoose(10.5, 7.5)
#' @export
gchoose <- function(n, k) {
  len <- max(length(n), length(k))
  a <- rep_len(n + 1, len)
  b <- rep_len(k + 1, len)
  c <- rep_len(n - k + 1, len)
  out <- gamma(a) / (gamma(b) * gamma(c))
  # gamma() overflows for arguments above ~171, giving Inf / Inf = NaN. Where
  # that happens and all arguments are positive, compute on the log scale.
  overflow <- !is.na(a) & !is.na(b) & !is.na(c) & pmax(a, b, c) > 170 & a > 0 & b > 0 & c > 0
  out[overflow] <- exp(lgamma(a[overflow]) - lgamma(b[overflow]) - lgamma(c[overflow]))
  out
}

#' Cumulative probability density function under the Gamma-extended Binomial distribution.
#'
#' @description Extends the cumulative Binomial probability mass function to positive non-integers, effectively turning the mass-function into a density-function.
#' @param q Vector of quantiles.
#' @param size Number of "trials" (zero or more). Need not be integer.
#' @param prob Probability of "success" on each "trial". Need not be integer. \code{q} and \code{prob} are recycled to a common length.
#' @param lower.tail Logical. If TRUE (default), probabilities are P[X<x], otherwise, P[X >= x]. Note that this differs from base-R \code{binom()} functions.
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @examples
#' # Assume some variable follows a Gamma-Binomial  distribution with
#' # "number of trials" = 10.5 and probability of "success" for each "trial"
#' # = 0.75, to compute the cumulative probability to attain a "number of
#' # success" below a specific point (e.g., less than 7.5 "successes":
#' pGammaBinom(q = 7.5, size = 10.5, prob = 0.75)
#'
#' # Conversely, to attain a value at or above 7.5:
#' pGammaBinom(q = 7.5, size = 10.5, prob = 0.75, lower.tail = FALSE)
#' @export
pGammaBinom <- function(q, size, prob, lower.tail = TRUE) {
  n <- max(length(q), length(prob))
  q <- rep_len(q, n)
  prob <- rep_len(prob, n)
  vapply(seq_len(n), function(i) .gammabinom_cdf(size, prob[i])(q[i], lower.tail), numeric(1))
}

# Returns a function computing the cumulative probability of the
# Gamma-extended Binomial distribution, with the normalizing integral computed
# once.
.gammabinom_cdf <- function(size, prob) {
  kernel <- function(y) dGammaBinom(y, size, prob)
  total <- stats::integrate(kernel, lower = 0, upper = size)$value
  function(q, lower.tail = TRUE) {
    q <- min(max(q, 0), size)
    if (lower.tail) {
      if (q == 0) return(0)
      stats::integrate(kernel, lower = 0, upper = q)$value / total
    } else {
      if (q == size) return(0)
      stats::integrate(kernel, lower = q, upper = size)$value / total
    }
  }
}

#' Probability density function under the Gamma-extended Binomial distribution.
#'
#' @param x Vector of quantiles.
#' @param size Number of "trials" (zero or more). Need not be integer.
#' @param prob Probability of "success" on each "trial". Need not be integer.
#' @param nc Whether to include a normalizing constant making sure that the sum of the distribution's density is 1.
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @examples
#' # Assume some variable follows a Gamma-Binomial distribution with
#' # "number of trials" = 10.5 and probability of "success" for each "trial"
#' # = 0.75, to compute the probability density to attain a "number of success"
#' # at a specific point (e.g., 7.5 "successes"):
#' dGammaBinom(x = 7.5, size = 10.5, prob = 0.75)
#'
#' # Including a normalizing constant (then diverges from binomial dist.):
#' dGammaBinom(x = 7.5, size = 10.5, prob = 0.75, nc = TRUE)
#' dGammaBinom(x = 7, size = 10, prob = 0.75) == dbinom(7, 10, 0.75)
#' dGammaBinom(x = 7, size = 10, prob = 0.75, nc = TRUE) == dbinom(7, 10, 0.75)
#' @export
dGammaBinom <- function(x, size, prob, nc = FALSE) {
  kernel <- function(y) {
    len <- max(length(y), length(prob))
    y <- rep_len(y, len)
    pr <- rep_len(prob, len)
    out <- numeric(len)
    out[is.na(y) | is.na(pr)] <- NA
    inside <- !is.na(out) & y >= 0 & y <= size
    out[inside] <- gchoose(size, y[inside]) * pr[inside]^y[inside] * (1 - pr[inside])^(size - y[inside])
    # For large `size` the coefficient overflows while the powers underflow;
    # compute those on the log scale instead.
    if (size > 170) {
      interior <- inside & pr > 0 & pr < 1
      yi <- y[interior]
      out[interior] <- exp(lgamma(size + 1) - lgamma(yi + 1) - lgamma(size - yi + 1) +
                             yi * log(pr[interior]) + (size - yi) * log1p(-pr[interior]))
    }
    out
  }
  if (nc) {
    kernel(x) / stats::integrate(kernel, lower = 0, upper = size)$value
  } else {
    kernel(x)
  }
}

#' Random number generation under the Gamma-extended Binomial distribution.
#'
#' @param n Number of observations.
#' @param size Number of "trials" (zero or more). Need not be integer.
#' @param prob Probability of "success" on each "trial". Need not be integer.
#' @param precision The convergence tolerance used by \code{qGammaBinom()} when computing the quantile. Default is 1e-4. Tuning this value will impact the time it takes for the search algorithm to arrive at an estimate.
#' @note Calls \code{qGammaBinom()}, which makes the random draw slower than what one might be used to (since \code{qGammaBinom()} calls \code{pGammaBinom()} and employs a search-algorithm to find the appropriate value down to a specifiable level of precision).
#' @examples
#' # Assume some variable follows a Gamma-Binomial distribution with
#' # "number of trials" = 10.5 and probability of "success" for each "trial"
#' # = 0.75 To draw a random value from this distribution:
#' rGammaBinom(n = 1, size = 10, prob = 0.75)
#' @export
rGammaBinom <- function(n, size, prob, precision = 1e-4) {
  qGammaBinom(stats::runif(n, 0, 1), size = size, prob = prob, precision = precision)
}

#' Quantile function for the Gamma-extended Binomial distribution.
#'
#' @param p Vector of probabilities.
#' @param size Number of "trials" (zero or more, including positive non-integers).
#' @param prob Probability of success on each "trial".
#' @param lower.tail Logical. If TRUE (default), probabilities are P[X < x], otherwise P[X > x].
#' @param precision The convergence tolerance of the root-finding algorithm (passed to \code{stats::uniroot()} as \code{tol}). Default is 1e-7. Tuning this value will impact the time it takes for the search algorithm to arrive at an estimate.
#' @note This function uses root-finding (\code{stats::uniroot()}) on the cumulative distribution function to find the number of successes corresponding to the specified quantile(s).
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @examples
#' # For a Gamma-extended Binomial distribution with number of trials = 10 and
#' # probability of success per trial of 0.75, calculate the number of success-
#' # ful trials at or below the 25% quantile:
#' qGammaBinom(p = 0.25, size = 10, prob = 0.75)
#'
#' # Conversely, for a Gamma-extended Binomial distribution with number of
#' # trials = 10 and probability of success per trial of 0.75, calculate the
#' # number of successful trials at or above the 25% quantile:
#' qGammaBinom(p = 0.25, size = 10, prob = 0.75, lower.tail = FALSE)
#' @export
qGammaBinom <- function(p, size, prob, lower.tail = TRUE, precision = 1e-7) {
  cdf <- .gammabinom_cdf(size, prob)
  vapply(p, function(target) {
    if (is.na(target) || target < 0 || target > 1) {
      warning("Probabilities must be within [0, 1]; returning NaN.")
      return(NaN)
    }
    # The quantile at probability 0 or 1 is an end-point of the support; which
    # end depends on the tail.
    if (target == 0) return(if (lower.tail) 0 else size)
    if (target == 1) return(if (lower.tail) size else 0)
    stats::uniroot(function(x) cdf(x, lower.tail) - target,
                   lower = 0, upper = size, tol = precision)$root
  }, numeric(1))
}


#' Probability Mass function for Lord's Two-Term Approximation to the Compound Binomial Distribution.
#'
#' @description Gives the density at \code{x} under Lord's two-term approximation to the compound Binomial PMF.
#' @param x Value of \code{x} (a specific number of successes).
#' @param N The total number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param p Probability of success for each trial.
#' @export
#' @examples
#' # Assume some variable follows a compound Binomial distribution with 100
#' # trials, a 50% probability of success on each trial, and Lord's k = 1. To
#' # compute the probability density at a specific point of the distribution
#' # (e.g., 50):
#' dcBinom(x = 50, N = 100, k = 1, p = .5)
dcBinom <- function(x, N, k, p) {
  stats::dbinom(x, N, p) -
    k * p * (1 - p) * (stats::dbinom(x, N - 2, p) - 2 * stats::dbinom(x - 1, N - 2, p) + stats::dbinom(x - 2, N - 2, p))
}

#' Cumulative Probability Mass function for Lord's Two-Term Approximation to the Compound Binomial Distribution.
#'
#' @description Function for calculating the proportion of observations up to a specifiable quantile under Lord's two-term approximation to the compound Binomial distribution.
#' @param q The quantile or vector of quantiles for which the proportion is to be calculated.
#' @param N Total number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param p Probability of success for each trial.
#' @param lower.tail Logical. If TRUE (default), probabilities are P[X<x], otherwise, P[X >= x]. Note that this differs from base-R \code{binom()} functions.
#' @return A vector of cumulative probabilities.
#' @export
#' @examples
#' # Assume some variable follows a compound Binomial distribution with 100
#' # trials, a 50% probability of success on each trial, and Lord's k = 1. To
#' # compute the cumulative probability at a specific point of the distribution
#' # (e.g., 50):
#' pcBinom(q = 50, N = 100, k = 1, p = .5)
pcBinom <- function(q, N, k, p, lower.tail = TRUE) {
  .pcbinom_tail(q, N, k, p, lower.tail = lower.tail)
}

#' Random Number Generation under Lord's Two-Term Approximation to the Compound Binomial Distribution.
#'
#' @description Random Number Generation under Lord's Two-Term Approximation to the Compound Binomial Distribution.
#' @param n Number of draws.
#' @param N Number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param p Probability of success for each trial.
#' @note For larger values of \code{k}, the distribution can yield negative probabilities. This function handles such occurrences by adding the absolute value of the minimum probability to all observations if there are any negative probabilities and then normalize the distribution so that the total density is equal to 1. If \code{p} is a vector it is recycled to length \code{n}, and one value is drawn per element.
#' @export
#' @examples
#' # To draw a sample of 50 values from a Compound-Binomial distribution with
#' # number of trials = 100, a 50% probability of success for each trial, and
#' # Lord's k = 1:
#' set.seed(1234)
#' rcBinom(n = 50, N = 100, k = 1, p = .5)
#'
#' # To draw values where the probabilities vary for each draw:
#' rcBinom(n = 50, N = 100, k = 1, p = runif(50))
rcBinom <- function(n, N, k, p) {
  if (length(p) == 1) {
    return(sample(0:N, size = n, prob = .fix_weights(dcBinom(0:N, N, k, p)), replace = TRUE))
  }
  p <- rep_len(p, n)
  negative <- FALSE
  draws <- vapply(p, function(pi) {
    weights <- dcBinom(0:N, N, k, pi)
    if (any(weights < 0)) {
      negative <<- TRUE
      weights <- weights + abs(min(weights))
    }
    sample(0:N, size = 1, prob = weights / sum(weights), replace = TRUE)
  }, integer(1))
  if (negative) {
    warning("Density function returned negative probabilities.")
  }
  draws
}

#' Probability Mass function for Lord's Beta Compound Binomial Distribution.
#'
#' @description Gives the density at \code{x} under the Beta Compound-Binomial distribution, where the Compound-Binomial distribution is Lord's two-term approximation.
#' @param x Value of \code{x} (a specific number of successes).
#' @param N Number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param l The lower-bound location parameter of the four-parameter Beta distribution.
#' @param u The upper-bound location parameter of the four-parameter Beta distribution.
#' @param alpha The first shape-parameter of the four-parameter Beta distribution.
#' @param beta The second shape-parameter of the four-parameter Beta distribution.
#' @return The probability mass at \code{x}.
#' @export
#' @examples
#' # Assume some variable follows a Beta compound Binomial distribution with 100
#' # trials, Lord's k = 1, and probabilities of successful trials drawn from a
#' # four-parameter Beta distribution with location-parameters l = .15 and u =
#' # .85, and shape parameters alpha = 6 and beta = 4. To compute the
#' # probability density at a specific point of the distribution (e.g., 50):
#' dBetacBinom(x = 50, N = 100, k = 1, l = .15, u = .85, alpha = 6, beta = 4)
dBetacBinom <- function(x, N, k, l, u, alpha, beta) {
  vapply(x, function(xi) {
    .integrate_support(function(p) dcBinom(xi, N, k, p) * dBeta.4P(p, l, u, alpha, beta), l, u, 0, 1)
  }, numeric(1))
}

#' Random Number Generation under Lord's Beta Compound-Binomial Distribution.
#'
#' @description Random number generation under Lord's Beta Compound-Binomial distribution, where the Compound-Binomial distribution is Lord's two-term approximation.
#' @param x Number of draws.
#' @param N Number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param l The lower-bound location parameter of the four-parameter Beta distribution.
#' @param u The upper-bound location parameter of the four-parameter Beta distribution.
#' @param alpha The first shape-parameter of the four-parameter Beta distribution.
#' @param beta The second shape-parameter of the four-parameter Beta distribution.
#' @note For larger values of \code{k}, the distribution can yield negative probabilities. As in \code{rcBinom()}, these are handled by shifting all probabilities up by the absolute value of the minimum and re-normalizing, with a warning.
#' @export
#' @examples
#' # To draw a sample of 50 values from a Beta Compound-Binomial distribution
#' # with number of trials = 100, Lord's k = 1, and probabilities of successful
#' # trials drawn from a four-parameter Beta distribution with location-
#' # parameters l = .15 and u = .85, and shape parameters alpha = 6 and
#' # beta = 4:
#' rBetacBinom(x = 50, N = 100, k = 1, l = .15, u = .85, alpha = 6, beta = 4)
rBetacBinom <- function(x, N, k, l, u, alpha, beta) {
  sample(0:N, size = x, prob = .fix_weights(dBetacBinom(0:N, N, k, l, u, alpha, beta)), replace = TRUE)
}
