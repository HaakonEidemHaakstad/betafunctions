# betafunctions

Functions for two- and four-parameter Beta distributions and for estimating the accuracy and consistency of classifications based on test scores.

## Installation

```r
install.packages("betafunctions")

# Development version:
# install.packages("remotes")
remotes::install_github("HaakonEidemHaakstad/betafunctions")
```

## Classification accuracy and consistency

Three single-administration approaches are implemented. Each has a function for one cut-point (`.CA`), one for several cut-points (`.CA.MC`), and an ROC function (`.ROC`):

| Approach | Functions | Input | Model |
|---|---|---|---|
| Livingston and Lewis (1995) | `LL.CA()`, `LL.CA.MC()`, `LL.ROC()` | Total scores and a reliability coefficient | Four- (or two-) parameter Beta true scores, Binomial errors on the effective-test-length scale |
| Hanson and Brennan (Hanson, 1991) | `HB.CA()`, `HB.CA.MC()`, `HB.ROC()` | Sum scores, test length and a reliability coefficient | Four- (or two-) parameter Beta true scores, Lord's two-term compound-Binomial errors |
| Multinomial / compound multinomial (Lee, 2007; Lee, Brennan & Wan, 2009) | `MN.CA()`, `MN.CA.MC()`, `MN.ROC()` | Item scores, optionally in weighted item sets | Examinee-level multinomial errors; no true-score distribution |

```r
library(betafunctions)

# 1000 examinees, 20 dichotomous items.
set.seed(1234)
p <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
items <- sapply(1:20, function(i) rbinom(1000, 1, p))
scores <- rowSums(items)

# Livingston and Lewis, with coefficient alpha as the reliability estimate.
LL.CA(scores, reliability = cba(items), cut = 10, min = 0, max = 20)

# Hanson and Brennan, three categories, with results conditional on each score.
HB.CA.MC(scores, reliability = cba(items), cut = c(8, 12), testlength = 20,
         conditional = TRUE)

# Compound multinomial (MULT-CLASS), with Brennan and Lee's bias correction.
MN.CA(items, cut = 10)

# Category-specific indices in one table.
MC.out.tabular(MN.CA.MC(items, cut = c(8, 12)))
```

The `.CA` functions return the confusion matrix and accuracy statistics (sensitivity, specificity, predictive values, Youden's J, accuracy) and the consistency matrix with agreement, chance agreement and Kappa. The LL and HB functions also report a chi-square test of model fit. `mdlfit.gfx()` plots observed against expected frequencies.

## Distributions and estimation

- Four-parameter Beta distribution: `dBeta.4P()`, `pBeta.4P()`, `qBeta.4P()`, `rBeta.4P()`, and the mean-variance parameterization `dBetaMS()`, `pBetaMS()`, `qBetaMS()`, `rBetaMS()`.
- Beta-Binomial distribution: `dBetaBinom()`, `pBetaBinom()`, `rBetaBinom()`.
- Lord's two-term compound Binomial: `dcBinom()`, `pcBinom()`, `rcBinom()`, `dBetacBinom()`, `rBetacBinom()`.
- Moments: `betamoments()`, `binomialmoments()`, `betabinomialmoments()`, `observedmoments()`.
- Fitting: `Beta.4p.fit()`, `Beta.2p.fit()`, and true-score fits `Beta.tp.fit()` (LL) and `HB.beta.tp.fit()` (HB).
- Reliability: Cronbach's alpha (`cba()`) and McDonald's omega (`mdo()`, maximum likelihood by default); effective test length (`ETL()`) and Lord's k (`Lords.k()`).

Functions deprecated in version 2.0.0 are listed in `?"betafunctions-deprecated"`.

## References

- Hanson, B. A. (1991). *Method of moments estimates for the four-parameter beta compound binomial model and the calculation of classification consistency indexes*. ACT Research Report Series.
- Lee, W.-C. (2007). Multinomial and compound multinomial error models for tests with complex item scoring. *Applied Psychological Measurement, 31*(4), 255–274.
- Lee, W.-C., Brennan, R. L., & Wan, L. (2009). Classification consistency and accuracy for complex assessments under the compound multinomial model. *Applied Psychological Measurement, 33*(5), 374–390.
- Livingston, S. A., & Lewis, C. (1995). Estimating the consistency and accuracy of classifications based on test scores. *Journal of Educational Measurement, 32*(2), 179–197.
- Lord, F. M. (1965). A strong true-score theory, with applications. *Psychometrika, 30*(3), 239–270.

## License

CC0. Bug reports and suggestions: <https://github.com/HaakonEidemHaakstad/betafunctions/issues>.
