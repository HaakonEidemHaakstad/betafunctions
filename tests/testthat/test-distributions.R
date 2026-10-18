test_that("four-parameter Beta d/p/q/r agree with base R on the rescaled scale", {
  x <- c(-1, 0.25, 0.3, 0.5, 0.75, 0.8)
  expect_equal(dBeta.4P(x, 0.25, 0.75, 5, 3), dbeta((x - 0.25) / 0.5, 5, 3) / 0.5)
  expect_equal(pBeta.4P(0.5, 0.25, 0.75, 5, 3), pbeta(0.5, 5, 3))
  expect_equal(pBeta.4P(0.5, 0.25, 0.75, 5, 3, lower.tail = FALSE), pbeta(0.5, 5, 3, lower.tail = FALSE))
  expect_equal(qBeta.4P(0.3, 0.25, 0.75, 5, 3), qbeta(0.3, 5, 3) * 0.5 + 0.25)
  set.seed(1)
  r <- rBeta.4P(1000, 0.25, 0.75, 5, 3)
  expect_true(all(r >= 0.25 & r <= 0.75))
})

test_that("qBeta.4P upper tail inverts pBeta.4P upper tail (asymmetric shapes)", {
  q <- qBeta.4P(0.1, 0.25, 0.75, 5, 3, lower.tail = FALSE)
  expect_equal(pBeta.4P(q, 0.25, 0.75, 5, 3, lower.tail = FALSE), 0.1)
})

test_that("mean/variance parameterization recovers the moments", {
  a <- AMS(60, 400, 0, 100)
  b <- BMS(60, 400, 0, 100)
  m <- betamoments(a, b, 0, 100, types = c("raw", "central"), orders = 2)
  expect_equal(m$raw[[1]], 60)
  expect_equal(m$central[[2]], 400)
  expect_equal(AMS(60, sd = 20, l = 0, u = 100), a)
  expect_equal(qBetaMS(pBetaMS(55, 60, 400, l = 0, u = 100), 60, 400, l = 0, u = 100), 55)
  expect_warning(dBetaMS(0.5, 0.6, variance = 0.04, sd = 0.1), "Nonequivalent")
  expect_error(dBetaMS(0.5, 0.6), "variance")
})

test_that("location parameters are recovered from shape parameters and moments", {
  m <- betamoments(5, 3, 0.25, 0.75, types = c("raw", "central"), orders = 2)
  mu <- m$raw[[1]]
  s2 <- m$central[[2]]
  expect_equal(LABMSU(5, 3, u = 0.75, mean = mu), 0.25)
  expect_equal(LABMSU(5, 3, mean = mu, variance = s2), 0.25)
  expect_equal(UABMSL(5, 3, mean = mu, l = 0.25), 0.75)
  expect_equal(UABMSL(5, 3, mean = mu, variance = s2), 0.75)
  expect_warning(UABMSL(5, 3), "Not enough information")
  expect_warning(LABMSU(5, 3), "Not enough information")
})

test_that("moment functions return lists of moment orders", {
  m <- betamoments(5, 3, 0.25, 0.75)
  expect_named(m, c("raw", "central", "standardized"))
  expect_type(m$raw, "list")
  expect_length(m$raw, 4)
  expect_equal(m$standardized[[2]], 1)

  b <- binomialmoments(100, 0.75)
  expect_equal(b$raw[[1]], 75)
  expect_equal(b$central[[2]], 100 * 0.75 * 0.25)

  bb <- betabinomialmoments(30, 0, 1, 2, 3)
  expect_equal(bb$raw[[1]], 30 * 2 / 5)

  x <- c(1, 2, 2, 3, 7, NA)
  om <- observedmoments(x)
  expect_equal(om$raw[[1]], mean(x, na.rm = TRUE))
  expect_equal(om$central[[2]], var(x, na.rm = TRUE))
  # `correct = FALSE` uses the population standard deviation.
  y <- x[!is.na(x)]
  expect_equal(observedmoments(x, correct = FALSE)$standardized[[2]], 1)
  expect_equal(observedmoments(x, correct = TRUE)$standardized[[2]], 1)
  # Unbiased third central moment and k-statistic skewness/kurtosis.
  n <- length(y)
  d <- y - mean(y)
  expect_equal(om$central[[3]], n * sum(d^3) / ((n - 1) * (n - 2)))
  k2 <- var(y)
  k4 <- (n * (n + 1) * sum(d^4) - 3 * (n - 1) * sum(d^2)^2) / ((n - 1) * (n - 2) * (n - 3))
  expect_equal(om$standardized[[3]], om$central[[3]] / k2^1.5)
  expect_equal(om$standardized[[4]], k4 / k2^2 + 3)
  expect_equal(observedmoments(x, correct = FALSE)$central[[3]], mean(d^3))
  expect_length(observedmoments(x, orders = 6)$central, 6)
})

test_that("betamedian is exact and betamode warns outside its domain", {
  expect_equal(betamedian(2, 8), qbeta(0.5, 2, 8))
  expect_equal(betamedian(5, 3, 25, 150), qbeta(0.5, 5, 3) * 125 + 25)
  expect_equal(betamode(5, 3), 4 / 6)
  expect_warning(betamode(0.5, 3), "alpha > 1")
})

test_that("Beta-Binomial pmf and cdf are accurate and use the P[X < q] convention", {
  # Closed form for the standard Beta-Binomial.
  N <- 10
  x <- 0:N
  closed <- choose(N, x) * beta(x + 2, N - x + 3) / beta(2, 3)
  expect_equal(dBetaBinom(x, N, 0, 1, 2, 3), closed)

  # Four-parameter case: a proper distribution, compared to a tight integration.
  pmf <- dBetaBinom(0:50, 50, 0.25, 0.75, 5, 3)
  expect_true(all(pmf >= 0))
  expect_equal(sum(pmf), 1)
  ref <- integrate(function(p) dbinom(20, 50, p) * dBeta.4P(p, 0.25, 0.75, 5, 3), 0.25, 0.75, rel.tol = 1e-12)$value
  expect_equal(pmf[21], ref)

  expect_equal(pBetaBinom(25, 50, 0.25, 0.75, 5, 3), sum(pmf[1:25]))
  expect_equal(pBetaBinom(25, 50, 0.25, 0.75, 5, 3, lower.tail = FALSE), sum(pmf[26:51]))
  expect_equal(pBetaBinom(25.5, 50, 0.25, 0.75, 5, 3), sum(pmf[1:26]))
  expect_equal(pBetaBinom(c(-1, 0, 51), 50, 0.25, 0.75, 5, 3), c(0, 0, 1))
})

test_that("Lord's two-term compound Binomial functions are consistent", {
  pmf <- dcBinom(0:10, 10, 1.2, 0.4)
  expect_equal(sum(pmf), 1)
  expect_equal(pcBinom(c(0, 3, 10, 11), 10, 1.2, 0.4), c(0, sum(pmf[1:3]), sum(pmf[1:10]), 1))
  expect_equal(pcBinom(3.5, 10, 1.2, 0.4), sum(pmf[1:4]))
  expect_equal(pcBinom(3, 10, 1.2, 0.4, lower.tail = FALSE), sum(pmf[4:11]))
  # k = 0 reduces to the Binomial.
  expect_equal(dcBinom(0:10, 10, 0, 0.4), dbinom(0:10, 10, 0.4))
})

test_that("rcBinom and rBetacBinom draw without spurious warnings", {
  set.seed(1)
  expect_silent(rcBinom(10, 20, 0.5, 0.5))
  expect_silent(rBetacBinom(10, 20, 0.5, 0.15, 0.85, 6, 4))
  expect_length(rcBinom(2, 10, 0.1, c(0.2, 0.4, 0.6)), 2)
  expect_length(rcBinom(5, 10, 0.1, c(0.2, 0.4)), 5)
  expect_warning(w <- rcBinom(200, 10, 0.5, 0.2), "set to zero")
  expect_true(all(dcBinom(w, 10, 0.5, 0.2) > 0))
})

test_that("Gamma-extended Binomial functions still work but are deprecated", {
  expect_equal(gchoose(c(8, 9, 10), c(3, 4, 5)), choose(c(8, 9, 10), c(3, 4, 5)))
  expect_equal(gchoose(200, 100), choose(200, 100))
  expect_warning(d <- dGammaBinom(0:10, 10, 0.75), "deprecated")
  expect_equal(d, dbinom(0:10, 10, 0.75))
  expect_warning(pGammaBinom(7, 10, 0.75), "deprecated")
  expect_warning(qGammaBinom(0.5, 10, 0.75), "deprecated")
  expect_warning(rGammaBinom(1, 10, 0.75), "deprecated")
  expect_warning(dBeta.pGammaBinom(0.5, 0.25, 0.75, 5, 3, 10, 0.5), "deprecated")
  expect_warning(dBeta.pBeta(0.5, 0.25, 0.75, 5, 3, 10, 0.5), "dBeta.pBinom")
  expect_warning(MLA(10, 8), "deprecated")

  # Behavior of the (internal) implementations.
  expect_equal(.dgammabinom(100, 200, 0.5), dbinom(100, 200, 0.5))
  expect_equal(.dgammabinom(c(-1, 11), 10, 0.75), c(0, 0))
  expect_length(.pgammabinom(7.5, 10.5, c(0.5, 0.75)), 2)
  expect_equal(.pgammabinom(7.5, 10.5, 0.75) + .pgammabinom(7.5, 10.5, 0.75, lower.tail = FALSE), 1)
  q <- .qgammabinom(0.25, 10, 0.75)
  expect_equal(.pgammabinom(q, 10, 0.75), 0.25, tolerance = 1e-5)
  qu <- .qgammabinom(0.25, 10, 0.75, lower.tail = FALSE)
  expect_equal(.pgammabinom(qu, 10, 0.75, lower.tail = FALSE), 0.25, tolerance = 1e-5)
  expect_equal(.qgammabinom(c(0, 1), 10, 0.75), c(0, 10))
  expect_equal(.qgammabinom(c(0, 1), 10, 0.75, lower.tail = FALSE), c(10, 0))
})

test_that("method-of-moments fits", {
  p4 <- Beta.4p.fit(mean = 0.5, variance = 0.01, skewness = -0.3, kurtosis = 2.6)
  m <- betamoments(p4$alpha, p4$beta, p4$l, p4$u)
  expect_equal(m$raw[[1]], 0.5)
  expect_equal(m$central[[2]], 0.01)
  expect_equal(m$standardized[[3]], -0.3)
  expect_equal(m$standardized[[4]], 2.6)
  expect_error(Beta.4p.fit(mean = 0.5), "Supply either")

  p2 <- Beta.2p.fit(mean = 0.5, variance = 0.01, l = 0.2, u = 0.9)
  expect_equal(p2$alpha, AMS(0.5, 0.01, 0.2, 0.9))
})

test_that("plotting coordinate helpers close the polygon", {
  coords <- Beta.gfx.poly.pdf(0.4, 0.6, 0.1, 5, 3, 0.25, 0.75)
  expect_equal(colnames(coords), c("x", "y"))
  expect_equal(coords[, "y"], c(0, dBeta.4P(seq(0.4, 0.6, 0.1), 0.25, 0.75, 5, 3), 0))
  expect_equal(unname(Beta.gfx.poly.qdf(0.4, 0.6, 0.1, 5, 3, 0.25, 0.75)[1, "y"]), 0.25)
})
