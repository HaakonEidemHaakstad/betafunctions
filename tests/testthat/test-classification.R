make_scores <- function(n = 1000, items = 100, seed = 1234) {
  set.seed(seed)
  rbinom(n, items, rBeta.4P(n, 0.25, 0.75, 5, 3))
}

test_that("ETL and R.ETL are inverses", {
  etl <- ETL(50, 100, 0, 100, 0.7)
  expect_equal(R.ETL(50, 100, 0, 100, etl), 0.7)
})

test_that("confmat, caStats and ccStats", {
  cm <- confmat(1, 2, 3, 4)
  expect_equal(dimnames(cm), list(c("True", "False", "Total"), c("Positive", "Negative", "Total")))
  expect_equal(cm["Total", "Total"], 10)
  expect_equal(sum(confmat(1, 2, 3, 4, "prop")[1:2, 1:2]), 1)
  s <- caStats(tp = 40, tn = 30, fp = 10, fn = 20)
  expect_equal(s$Sensitivity, 40 / 60)
  expect_equal(s$Accuracy, 0.7)
  # Frequencies and proportions give the same consistency statistics.
  expect_equal(ccStats(40, 10, 10, 40), ccStats(0.4, 0.1, 0.1, 0.4))
  expect_equal(ccStats(0.4, 0.1, 0.1, 0.4)$Kappa, 0.6)
})

test_that("LL.CA output structure and internal consistency", {
  x <- make_scores()
  out <- LL.CA(x, reliability = 0.7, cut = 50, min = 0, max = 100)
  expect_named(out, c("parameters", "modelfit", "confusionmatrix", "classification.accuracy",
                      "consistencymatrix", "classification.consistency"))
  expect_equal(sum(out$confusionmatrix[1:2, 1:2]), 1, tolerance = 1e-6)
  expect_equal(sum(out$consistencymatrix), 1)
  expect_equal(out$consistencymatrix["i", "j"], out$consistencymatrix["j", "i"])

  # Binary and multiple-category versions agree for a single cut-point,
  # including cut-points that do not fall on an integer of the effective test
  # length.
  for (cut in c(50, 37)) {
    b <- LL.CA(x, 0.7, cut, 0, 100, modelfit = NULL)
    m <- LL.CA.MC(x, 0.7, cut, 0, 100, modelfit = NULL)
    expect_equal(unname(m$consistency$overall$consistencymatrix), unname(b$consistencymatrix), tolerance = 1e-6)
    expect_equal(m$consistency$overall$statistics$Kappa, b$classification.consistency$Kappa, tolerance = 1e-6)
    expect_equal(m$accuracy$overall$accuracy, b$classification.accuracy$Accuracy, tolerance = 1e-6)
  }
})

test_that("LL.CA model fit bins observed and expected scores identically", {
  x <- make_scores()
  out <- LL.CA(x, reliability = 0.7, cut = 50, min = 0, max = 100, modelfit = c(nbins = 100, minbin = 10))
  # The effective test length is not a multiple of the number of bins, so
  # bins have fractional widths.
  expect_false(out$parameters$etl_rounded %% 100 == 0)
  tab <- out$modelfit$contingencytable
  expect_equal(sum(tab[2, ]), length(x))
  expect_equal(sum(tab[1, ]), length(x))
  expect_true(all(tab[1, ] >= 10))
  expect_equal(out$modelfit$df, ncol(tab) - 4)
  expect_equal(out$modelfit$pvalue, pchisq(out$modelfit$chisquared, out$modelfit$df, lower.tail = FALSE))

  two <- LL.CA(x, 0.7, 50, 0, 100, true.model = "2P")
  expect_equal(two$modelfit$df, ncol(two$modelfit$contingencytable) - 2)
  expect_equal(two$modelfit$pvalue,
               pchisq(two$modelfit$chisquared, two$modelfit$df, lower.tail = FALSE))
})

test_that("LL.CA accepts a parameter list", {
  tp <- list(l = 0.25, u = 0.75, alpha = 5, beta = 3, etl = 50)
  out <- LL.CA(tp, cut = 50, min = 0, max = 100)
  expect_null(out$modelfit)
  expect_equal(out$parameters$etl_rounded, 50)
})

test_that("LL.CA.MC category labels use the original score scale", {
  x <- make_scores(items = 100) + 100
  out <- LL.CA.MC(x, 0.7, c(150, 160), min = 100, max = 200, modelfit = NULL)
  expect_equal(rownames(out$accuracy$overall$confusionmatrix),
               c("Observed < 150", " >= 150 & < 160", " >= 160"))
  expect_equal(sum(out$accuracy$overall$confusionmatrix), 1)
  tab <- MC.out.tabular(out)
  expect_equal(dim(tab), c(3, 13))
  expect_equal(colnames(tab)[1:5], c("TP", "FP", "TN", "FN", "Sensitivity"))
})

test_that("Hanson-Brennan functions", {
  set.seed(1234)
  p <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
  raw <- sapply(1:20, function(i) rbinom(1000, 1, p))
  x <- rowSums(raw)
  rel <- cba(raw)
  out <- HB.CA(x, rel, cut = 10, testlength = 20)
  expect_named(out, c("parameters", "modelfit", "confusionmatrix", "classification.accuracy",
                      "consistencymatrix", "classification.consistency"))
  expect_equal(rownames(out$modelfit$contingencytable), c("Expected", "Observed"))
  expect_equal(sum(out$modelfit$contingencytable[2, ]), 1000)

  mc <- HB.CA.MC(x, rel, cut = 10, testlength = 20, modelfit = NULL)
  expect_equal(unname(mc$consistency$overall$consistencymatrix), unname(out$consistencymatrix), tolerance = 1e-6)

  # A cut at zero places nobody below the cut.
  zero <- HB.CA(out$parameters, cut = 0, truecut = 10, output = "a")
  expect_equal(zero$classification.accuracy$Sensitivity, 0)

  roc <- HB.ROC(x, rel, 20, 10, raw.out = TRUE)
  expect_equal(colnames(roc), c("FPR", "TPR", "Youden.J", "Cutoff", "Accuracy", "PPV", "NPV"))
  expect_equal(nrow(roc), 21)
})

test_that("ROC helpers", {
  x <- make_scores()
  roc <- LL.ROC(x, 0.7, 0, 100, 50, raw.out = TRUE, grainsize = 20)
  expect_equal(nrow(roc), 21)
  expect_true(AUC(roc[, "FPR"], roc[, "TPR"]) > 0.5)
  expect_equal(AUC(c(0, 1), c(0, 1)), 0.5)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(LL.ROC(x, 0.7, 0, 100, 50, AUC = TRUE, maxJ = TRUE, maxAcc = TRUE,
                       locate = c("specificity", 0.8), grainsize = 20))
  out <- LL.CA(x, 0.7, 50, 0, 100, modelfit = c(nbins = 20, minbin = 1))
  expect_silent(mdlfit.gfx(out))
})

test_that("reliability coefficients", {
  set.seed(1234)
  x <- matrix(rnorm(500), ncol = 5) + rnorm(100)
  colnames(x) <- letters[1:5]
  w <- mdo(x)
  expect_true(w > 0 && w < 1)
  f <- mdo(x, fit = TRUE)
  expect_equal(f$Omega, w)
  expect_equal(colnames(f$Parameters), letters[1:5])
  # The GFI does not depend on the scale of the data.
  expect_equal(mdo(x * 10, fit = TRUE)$GFI, f$GFI)
  expect_error(mdo(x[, 1:2]), "three")
  expect_true(cba(x) > 0)
})

test_that("factorials and true-score moments", {
  expect_equal(dfac(5, 3), 60)
  expect_equal(afac(5, 3), 210)
  expect_equal(dfac(c(3.14, 7), 3, "gamma"), gamma(c(4.14, 8)) / gamma(c(1.14, 5)))
  x <- make_scores()
  expect_equal(tsm(x, 1, 100), mean(x) / 100)
  expect_equal(HB.tsm(x, 4, 100, 0)[2], tsm(x, 2, 100))
  fit <- Beta.tp.fit(x, 0, 100, 100)
  expect_named(fit, c("alpha", "beta", "l", "u", "etl"))
  hb <- HB.beta.tp.fit(x, 100, 0)
  expect_equal(hb[c("alpha", "beta", "l", "u")], fit[c("alpha", "beta", "l", "u")])
  expect_error(Beta.tp.fit(x, 0, 100, 100, true.model = "3P"), "true.model")
})
