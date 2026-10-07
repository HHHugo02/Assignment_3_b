# Assignment 3B: EM for mixtures of t-distributions

Presentation: `presentation.html` (open in a browser; press `S` for speaker notes,
`F` for full screen). Source: `presentation.qmd`.

The presentation is built up step by step:

1. A simple implementation (density scale, `dt()`)
2. Is it correct? (theory, `numDeriv`, simulation, `optim()`, `teigen`, `mclust`)
3. Where is it slow? (`Rprof`, line profiling, call counts)

The first version of the whole project (log scale, Rcpp, Fisher information,
robustness study, ...) is kept in `old/`.

## Files

| File | Content |
|------|---------|
| `R/tmix-common.R` | shared basics: t-density, mixture density, simulation, starting value |
| `R/em.R` | the EM driver `em()` (shared) |
| `R/tmix-simple.R` | version 1: log-likelihood, E-step, `Q()`, `grad_Q()` as plain functions |
| `R/gem.R` | version 1: the GEM step (`create_gem_step(x, nu)`) and `gem()` |
| `R/tmix-classes.R` | version 2: S3 classes `tmix` (parent), `tmix_density` and `tmix_log` |
| `R/gem-classes.R` | version 2: the GEM step for any `tmix` model and `fit()` |
| `R/gauss-mix.R` | two-component Gaussian mixture fitted by EM (`fit_gauss()`), for the outlier comparison |
| `R/profiling-tools.R` | `count_calls()`, `profile_expr()`, `profile_lines()` |
| `R/profiling-demo.R` | interactive `profvis` demo for the exam |
| `tests/test-simple.R` | tests for version 1 |
| `tests/test-classes.R` | tests for version 2 (each check runs for both classes) |

## Run

From this folder in R. Version 2 (S3 classes):

```r
source("R/load.R")
testthat::test_file("tests/test-classes.R")
x <- rtmix(400, c(0.6, 0, 3, 1, 0.7), c(5, 10))
fit(tmix_log(x, c(5, 10)))
fit(tmix_density(x, c(5, 10)))
```

Version 1 (plain functions, used in the first part of the presentation):

```r
source("R/load-v1.R")
gem(x, c(5, 10))
```

Rebuild the slides:

```
quarto render presentation.qmd
```

Required packages: ggplot2, patchwork, testthat, numDeriv, bench, profvis,
teigen, mclust.
