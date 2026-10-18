# Examples distributed with MULT-CLASS v3.0 (Lee, 2008): item data and the
# printed program output (five decimals).
# The printed output has five decimals.
near <- function(actual, expected) expect_lt(max(abs(actual - expected)), 5.01e-6)

mc_dich <- function() {
  matrix(c(
    1, 1, 0, 0, 1, 1, 1, 1, 0, 1,
    1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
    1, 1, 0, 1, 1, 1, 1, 0, 0, 1,
    0, 1, 0, 0, 1, 1, 0, 0, 0, 1,
    1, 1, 0, 0, 1, 1, 1, 0, 1, 1,
    1, 1, 0, 0, 1, 1, 1, 0, 0, 1,
    1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
    0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
    1, 1, 0, 0, 1, 1, 0, 0, 0, 1,
    1, 1, 1, 0, 0, 1, 1, 0, 0, 0
  ), nrow = 10, byrow = TRUE)
}
mc_mix <- function() {
  cbind(
  matrix(c(
    2, 2, 2, 2, 2, 4, 2, 2,
    2, 6, 4, 6, 4, 4, 4, 2,
    4, 4, 4, 6, 6, 2, 2, 2,
    8, 6, 8, 8, 4, 8, 2, 4,
    4, 8, 4, 8, 4, 2, 2, 2,
    4, 2, 2, 8, 2, 4, 2, 4,
    4, 2, 2, 4, 4, 2, 2, 2,
    2, 6, 6, 6, 4, 4, 2, 4,
    6, 6, 8, 8, 8, 8, 8, 8,
    2, 2, 2, 6, 2, 4, 2, 2
  ), nrow = 10, byrow = TRUE),
  matrix(c(
    0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
    1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1,
    1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1,
    1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
    1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0,
    1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1,
    1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
    1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1
  ), nrow = 10, byrow = TRUE))
}

test_that("MN.CA.MC reproduces the MULT-CLASS dichotomous example (no correction)", {
  out <- MN.CA.MC(mc_dich(), 5, correct = FALSE, conditional = TRUE)
  near(out$parameters$observed.variance, 3.64)
  near(out$parameters$error.variance, 2.33333)
  near(out$parameters$reliability, 0.35897)
  near(out$consistency$overall$statistics$p, 0.75198)
  near(out$consistency$overall$statistics$p_c, 0.55809)
  near(out$consistency$overall$statistics$Kappa, 0.43877)
  near(out$accuracy$overall$accuracy, 0.81711)
  near(out$accuracy$overall$overclassification, 0.07666)
  near(out$accuracy$overall$underclassification, 0.10623)
  expect_equal(round(out$conditional$phi, 5),
               c(0.90979, 0.53543, 0.90979, 0.53543, 0.90979, 0.72279, 0.99971, 0.93656, 0.53028, 0.53028))
  expect_equal(round(out$conditional$gamma, 5),
               c(0.95265, 0.63310, 0.95265, 0.63310, 0.95265, 0.83376, 0.99985, 0.96721, 0.62305, 0.62305))
})

test_that("MN.CA.MC reproduces the MULT-CLASS mixed-format example (with correction)", {
  out <- MN.CA.MC(mc_mix(), c(30, 49, 65, 90), sets = list(1:8, 9:28), weights = c(1, 2),
                  scorepoints = list(c(2, 4, 6, 8), c(0, 1)), correct = TRUE, conditional = TRUE)
  near(out$parameters$observed.variance, 354.56)
  near(out$parameters$error.variance, 39.43459)
  near(out$parameters$reliability, 0.88878)
  near(out$consistency$overall$statistics$p, 0.62527)
  near(out$consistency$overall$statistics$p_c, 0.29565)
  near(out$consistency$overall$statistics$Kappa, 0.46799)
  near(out$accuracy$overall$accuracy, 0.72649)
  near(out$accuracy$overall$overclassification, 0.15871)
  near(out$accuracy$overall$underclassification, 0.11480)
  expect_equal(out$conditional$total, c(28, 46, 52, 82, 46, 56, 40, 52, 96, 54))
  expect_equal(round(out$conditional$overclassification, 5),
               c(0.52511, 0.33653, 0.02188, 0.11891, 0.36633, 0.09854, 0.07025, 0.02113, 0, 0.02843))
  expect_equal(round(out$conditional$underclassification, 5),
               c(0, 0.00102, 0.30997, 0.01927, 0.00709, 0.15574, 0.01172, 0.30883, 0.15406, 0.18032))
  # Same structure as LL.CA.MC(), so MC.out.tabular() works.
  expect_equal(dim(MC.out.tabular(out)), c(5, 13))
  expect_equal(rownames(out$accuracy$overall$confusionmatrix)[1], "Observed < 30")
})

test_that("MN.CA agrees with MN.CA.MC for a single cut-point", {
  x <- mc_dich()
  for (correct in c(TRUE, FALSE)) {
    b <- MN.CA(x, 5, correct = correct)
    m <- MN.CA.MC(x, 5, correct = correct)
    expect_named(b, c("parameters", "confusionmatrix", "classification.accuracy",
                      "consistencymatrix", "classification.consistency"))
    expect_equal(unname(b$consistencymatrix), unname(m$consistency$overall$consistencymatrix))
    expect_equal(b$classification.consistency$Kappa, m$consistency$overall$statistics$Kappa)
    expect_equal(b$classification.accuracy$Accuracy, m$accuracy$overall$accuracy)
    expect_equal(sum(b$confusionmatrix[1:2, 1:2]), 1)
  }
  expect_null(MN.CA(x, 5)$conditional)
  expect_s3_class(MN.CA(x, 5, conditional = TRUE)$conditional, "data.frame")
})

test_that("multinomial model inputs are validated", {
  x <- mc_dich()
  expect_error(MN.CA(x, c(4, 6)), "MN.CA.MC")
  expect_error(MN.CA(x, 5, sets = list(1:5)), "exactly one item set")
  expect_error(MN.CA(x, 5, sets = list(1, 2:10)), "at least two items")
  expect_error(MN.CA(x, 5, scorepoints = list(c(0, 2))), "declared score points")
  expect_error(MN.CA.MC(x, c(6, 4)), "increasing")
  x[1, 1] <- NA
  expect_error(MN.CA(x, 5), "Missing")
  colnames(x) <- paste0("i", 1:10)
  y <- mc_dich()
  colnames(y) <- paste0("i", 1:10)
  expect_equal(MN.CA(y, 5, sets = list(paste0("i", 1:5), paste0("i", 6:10)))$parameters$sets, 2)
})

test_that("MN.ROC", {
  set.seed(1234)
  p <- rBeta.4P(500, 0.15, 0.85, 6, 4)
  x <- sapply(1:20, function(i) rbinom(500, 1, p))
  roc <- MN.ROC(x, truecut = 10, raw.out = TRUE, grainsize = 20)
  expect_equal(colnames(roc), c("FPR", "TPR", "Youden.J", "Cutoff", "Accuracy", "PPV", "NPV"))
  expect_equal(nrow(roc), 21)
  expect_equal(unname(roc[11, "Accuracy"]), MN.CA(x, 10)$classification.accuracy$Accuracy)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(MN.ROC(x, truecut = 10, grainsize = 20, AUC = TRUE, maxJ = TRUE))
})
