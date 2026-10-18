## Submission

This is a major update (1.9.0 -> 2.0.0): an audit and refactor of all functions, bug fixes, new functions for the compound multinomial approach to classification accuracy and consistency, and deprecation of functions without a theoretical basis (they still work, with a deprecation warning). See NEWS.md.

## Test environments

* Local: Ubuntu 24.04, R 4.3.3 (`R CMD check --as-cran`).
* TODO before submitting: win-builder (R-devel and R-release), R-hub / mac-builder.

## R CMD check results

0 errors | 0 warnings | 0 notes

(Local check run without internet access, so the incoming URL and system-clock checks were skipped.)

## Reverse dependencies

There is one reverse dependency, IRTest, which imports only `rBeta.4P()`. That function's arguments, implementation and output are unchanged in this release.
