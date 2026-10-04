library(testthat)

# Tests for version 1 (plain functions)
if (!exists("gem", mode = "function")) {
  root <- if (file.exists("R/load-v1.R")) "." else ".."
  owd <- setwd(root)
  source("R/load-v1.R")
  setwd(owd)
}

set.seed(42)
nu <- c(5, 10)
par_true <- c(0.6, 0, 3, 1, 0.7)
x <- rtmix(400, par_true, nu)
par0 <- default_start(x)
par_rand <- c(0.4, -0.5, 2.5, 1.3, 0.9)

ptmix <- function(q, par, nu) {
  par[1] * pt((q - par[2]) / par[4], nu[1]) +
    (1 - par[1]) * pt((q - par[3]) / par[5], nu[2])
}

# Building blocks ---------------------------------------------------------------

test_that("t density matches the formula from the assignment", {
  f_formula <- function(x, mu, sigma, nu) {
    gamma((nu + 1) / 2) / (sqrt(pi * nu * sigma^2) * gamma(nu / 2)) *
      (1 + (x - mu)^2 / (nu * sigma^2))^(-(nu + 1) / 2)
  }
  z <- seq(-10, 10, length.out = 101)
  expect_equal(dt_ls(z, 1, 2, 3), f_formula(z, 1, 2, 3))
})

test_that("mixture density integrates to 1", {
  expect_equal(integrate(dtmix, -Inf, Inf, par = par_true, nu = nu)$value, 1,
    tolerance = 1e-6
  )
})

test_that("simulated data follow the mixture distribution (KS test)", {
  set.seed(3)
  xs <- rtmix(20000, par_true, nu)
  expect_gt(ks.test(xs, ptmix, par = par_true, nu = nu)$p.value, 0.01)
})

test_that("E-step probabilities are calibrated", {
  set.seed(3)
  d <- sim_tmix(20000, par_true, nu)
  gam <- e_step(par_true, d$x, nu)
  bin <- cut(gam, seq(0, 1, 0.1), include.lowest = TRUE)
  mean_gam <- tapply(gam, bin, mean)
  frac_z1 <- tapply(d$z == 1, bin, mean)
  se <- sqrt(mean_gam * (1 - mean_gam) / tapply(gam, bin, length))
  expect_true(all(abs(frac_z1 - mean_gam) < 4 * se + 0.01))
})

# Derivatives -------------------------------------------------------------------

test_that("gradient of Q matches numerical differentiation", {
  gam <- e_step(par0, x, nu)
  num <- numDeriv::grad(function(p) Q(p, gam, x, nu), par_rand)
  expect_equal(grad_Q(par_rand, gam, x, nu), num, tolerance = 1e-6)
})

test_that("Fisher's identity: grad_Q(par | par) is the score", {
  num <- numDeriv::grad(function(p) loglik(p, x, nu), par_rand)
  expect_equal(grad_Q(par_rand, e_step(par_rand, x, nu), x, nu), num, tolerance = 1e-6)
})

# The algorithm -----------------------------------------------------------------

test_that("each GEM step increases Q and the log-likelihood", {
  step <- create_gem_step(x, nu)
  par <- par0
  for (i in 1:30) {
    gam <- e_step(par, x, nu)
    par_new <- step(par)
    expect_gte(Q(par_new, gam, x, nu), Q(par, gam, x, nu))
    expect_gte(loglik(par_new, x, nu), loglik(par, x, nu) - 1e-10)
    par <- par_new
  }
})

test_that("GEM and optim() find the same maximum", {
  par_gem <- gem(x, nu, eps = 1e-14)
  nll <- function(eta) {
    -loglik(c(plogis(eta[1]), eta[2:3], exp(eta[4:5])), x, nu)
  }
  eta <- optim(c(qlogis(par0[1]), par0[2:3], log(par0[4:5])), nll,
    method = "BFGS", control = list(reltol = 1e-14, maxit = 1000)
  )$par
  par_optim <- c(plogis(eta[1]), eta[2:3], exp(eta[4:5]))
  expect_equal(as.vector(par_gem), par_optim, tolerance = 1e-4)
  expect_lt(max(abs(grad_Q(par_gem, e_step(par_gem, x, nu), x, nu))), 1e-3)
})

# Against R packages ------------------------------------------------------------

test_that("GEM agrees with teigen when nu1 = nu2", {
  skip_if_not_installed("teigen")
  set.seed(1)
  x5 <- rtmix(1000, par_true, c(5, 5))
  par_gem <- gem(x5, c(5, 5), eps = 1e-14)
  tg <- teigen::teigen(x5,
    Gs = 2, models = "univUU", dfstart = 5, dfupdate = FALSE,
    scale = FALSE, verbose = FALSE, eps = c(1e-12, 1e-12)
  )
  o <- order(tg$parameters$mean) # teigen may order components differently
  par_tg <- c(
    tg$parameters$pig[o][1], tg$parameters$mean[o],
    sqrt(tg$parameters$sigma[1, 1, o])
  )
  expect_equal(as.vector(par_gem), par_tg, tolerance = 1e-4)
  expect_equal(loglik(par_gem, x5, c(5, 5)), tg$logl, tolerance = 1e-8)
})

test_that("GEM agrees with mclust in the Gaussian limit (nu very large)", {
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(library(mclust))
  set.seed(2)
  xg <- c(rnorm(600, 0, 1), rnorm(400, 3, 0.7))
  par_gem <- gem(xg, c(1e6, 1e6), eps = 1e-16, maxit = 5000)
  mc <- Mclust(xg,
    G = 2, modelNames = "V", verbose = FALSE,
    control = emControl(tol = c(1e-12, sqrt(.Machine$double.eps)))
  )
  par_mc <- unname(c(
    mc$parameters$pro[1], mc$parameters$mean,
    sqrt(mc$parameters$variance$sigmasq)
  ))
  expect_equal(as.vector(par_gem), par_mc, tolerance = 1e-4)
})

test_that("estimates are close to the truth for a large sample", {
  set.seed(4)
  x_big <- rtmix(20000, par_true, nu)
  expect_equal(as.vector(gem(x_big, nu)), par_true, tolerance = 0.05)
})
