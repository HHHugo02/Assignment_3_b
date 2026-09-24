# Assignment 3B: EM for mixtures of t-distributions

Presentation: `presentation.html` (open in a browser; press `S` for speaker notes,
`F` for full screen). Source: `presentation.qmd`.

## Files

| File | Content |
|------|---------|
| `R/tmix-model.R` | t log-density, mixture density, simulation, log-likelihood, E-step, `Q()`, `grad_Q()`, `Q_aug()` |
| `R/em-algorithms.R` | generic driver `em()`, convergence criteria, GEM step, augmented EM (readable + optimized), exact EM, Gaussian mixture EM |
| `R/direct-methods.R` | gradient ascent and BFGS on the marginal log-likelihood |
| `R/tracer.R` | small tracer (records parameters and time, like `CSwR::tracer()`) |
| `R/fisher.R` | Fisher information three ways + spectral radius |
| `R/tmix-class.R` | S3 class `tmix` with `print`, `plot`, `coef`, `logLik`, `vcov`, `predict` |
| `R/v1-naive.R` | the first versions (kept to show the development) |
| `R/profiling-demo.R` | interactive `profvis` demo for the exam |
| `src/em_aug_step.cpp` | Rcpp version of the augmented EM step |
| `tests/test-tmix.R` | testthat tests |

## Run

From this folder in R:

```r
source("R/load.R")
testthat::test_file("tests/test-tmix.R")
fit <- tmix(x = rtmix(400, c(0.6, 0, 3, 1, 0.7), c(5, 10)), nu = c(5, 10))
fit; plot(fit)
```

Rebuild the slides (takes about 1-2 minutes):

```
quarto render presentation.qmd
```

Required packages: ggplot2, patchwork, testthat, numDeriv, bench, profvis, Rcpp
(+ Rtools for Rcpp).
