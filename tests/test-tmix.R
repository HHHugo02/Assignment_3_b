library(testthat)

if (!exists("tmix", mode = "function")) {
  root <- if (file.exists("R/load.R")) "." else ".."
  owd <- setwd(root)
  source("R/load.R")
  setwd(owd)
}

set.seed(2024)
nu <- c(5, 10)
par_true <- c(0.6, 0, 3, 1, 0.7)
x <- rtmix(400, par_true, nu)
par0 <- default_start(x)
par_rand <- c(0.4, -0.5, 2.5, 1.3, 0.9)

test_that("t log-density matches dt() and the mixture density integrates to 1", {
  z <- seq(-10, 10, length.out = 101)
  expect_equal(ldt_ls(z, 1, 2, 3), dt((z - 1) / 2, 3, log = TRUE) - log(2))
  expect_equal(integrate(dtmix, -Inf, Inf, par = par_true, nu = nu)$value, 1,
    tolerance = 1e-6
  )
})

test_that("gradient of Q matches numerical differentiation", {
  num <- numDeriv::grad(function(p) Q(p, par0, x, nu), par_rand)
  expect_equal(grad_Q(par_rand, par0, x, nu), num, tolerance = 1e-6)
})

test_that("Fisher's identity: grad_Q(par | par) is the score", {
  num <- numDeriv::grad(function(p) loglik(p, x, nu), par_rand)
  expect_equal(grad_loglik(par_rand, x, nu), num, tolerance = 1e-6)
})

test_that("the three implementations of the augmented EM step agree", {
  expected <- create_em_aug_step_v1(x, nu)(par_rand)
  expect_equal(create_em_aug_step(x, nu)(par_rand), expected)
  expect_equal(create_em_aug_step_fast(x, nu)(par_rand), expected)
  expect_equal(em_aug_step_cpp(x, par_rand, nu), expected)
})

test_that("GEM steps increase Q, and all EM variants increase the log-likelihood", {
  gem_step <- create_gem_step(x, nu)
  par <- par0
  for (i in 1:20) {
    par_new <- gem_step(par)
    expect_gte(Q(par_new, par, x, nu), Q(par, par, x, nu))
    par <- par_new
  }
  for (step in list(gem_step, create_em_aug_step(x, nu), create_em_exact_step(x, nu))) {
    tr <- trace_run(function(cb) em(par0, step, maxit = 30, cb = cb), x, nu, "")
    expect_true(all(diff(tr$loglik) >= -1e-10))
  }
})

test_that("all methods converge to the same stationary point", {
  fits <- sapply(tmix_methods, function(m) tmix(x, nu, par0, method = m, eps = 1e-14)$par)
  expect_equal(unname(fits[, "gem"]), unname(fits[, "bfgs"]), tolerance = 1e-5)
  for (m in tmix_methods) {
    expect_equal(fits[, m], fits[, "em_aug"], tolerance = 1e-5)
  }
  par_hat <- fits[, "em_aug"]
  expect_lt(max(abs(grad_loglik(par_hat, x, nu))), 1e-3)
  # ... which is a fixed point of the EM map
  expect_equal(create_em_aug_step(x, nu)(par_hat), unname(par_hat), tolerance = 1e-6)
})

test_that("estimates are close to the truth for a large sample", {
  x_big <- rtmix(20000, par_true, nu)
  fit <- tmix(x_big, nu, default_start(x_big), method = "em_aug_cpp")
  expect_equal(unname(coef(fit)), par_true, tolerance = 0.05)
})

test_that("the three ways of computing the Fisher information agree", {
  par_hat <- tmix(x, nu, par0, eps = 1e-16)$par
  i1 <- fisher_score(par_hat, x, nu)
  expect_equal(fisher_louis(par_hat, x, nu)$iX, i1, tolerance = 1e-5)
  expect_equal(fisher_em_map(par_hat, x, nu)$iX, i1, tolerance = 1e-5)
})
