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
| `R/tmix-simple.R` | t-density, mixture density, simulation, log-likelihood, E-step, `Q()`, `grad_Q()` |
| `R/gem.R` | the GEM step (`create_gem_step()`), the EM driver `em()` and `gem()` |
| `R/profiling-tools.R` | `count_calls()`, `profile_expr()`, `profile_lines()` |
| `R/profiling-demo.R` | interactive `profvis` demo for the exam |
| `tests/test-simple.R` | testthat tests |

## Run

From this folder in R:

```r
source("R/load.R")
testthat::test_file("tests/test-simple.R")
x <- rtmix(400, c(0.6, 0, 3, 1, 0.7), c(5, 10))
gem(x, c(5, 10))
```

Rebuild the slides:

```
quarto render presentation.qmd
```

Required packages: ggplot2, patchwork, testthat, numDeriv, bench, profvis,
teigen, mclust.
