library(testthat)

# Tests for version 2 (S3 classes): every check is run for both implementations
if (!exists("tmix_log", mode = "function")) {
  root <- if (file.exists("R/load.R")) "." else ".."
  owd <- setwd(root)
  source("R/load.R")
  setwd(owd)
}

set.seed(42)
nu <- c(5, 10)
par_true <- c(0.6, 0, 3, 1, 0.7)
x <- rtmix(400, par_true, nu)
par0 <- default_start(x)
par_rand <- c(0.4, -0.5, 2.5, 1.3, 0.9)

models <- list(density = tmix_density(x, nu), log = tmix_log(x, nu))

test_that("both implementations are subclasses of tmix", {
  expect_s3_class(models$density, c("tmix_density", "tmix"), exact = TRUE)
  expect_s3_class(models$log, c("tmix_log", "tmix"), exact = TRUE)
})

# The same correctness checks for each implementation -----------------------------

for (name in names(models)) {
  model <- models[[name]]

  test_that(paste0(name, ": gradient of Q matches numerical differentiation"), {
    gam <- e_step(model, par0)
    num <- numDeriv::grad(function(p) Q(model, p, gam), par_rand)
    expect_equal(grad_Q(model, par_rand, gam), num, tolerance = 1e-6)
  })

  test_that(paste0(name, ": Fisher's identity"), {
    num <- numDeriv::grad(function(p) loglik(model, p), par_rand)
    expect_equal(grad_Q(model, par_rand, e_step(model, par_rand)), num, tolerance = 1e-6)
  })

  test_that(paste0(name, ": each GEM step increases Q and the log-likelihood"), {
    step <- create_gem_step(model)
    par <- par0
    for (i in 1:30) {
      gam <- e_step(model, par)
      par_new <- step(par)
      expect_gte(Q(model, par_new, gam), Q(model, par, gam))
      expect_gte(loglik(model, par_new), loglik(model, par) - 1e-10)
      par <- par_new
    }
  })

  test_that(paste0(name, ": GEM and optim() find the same maximum"), {
    par_gem <- fit(model, eps = 1e-14)
    nll <- function(eta) -loglik(model, c(plogis(eta[1]), eta[2:3], exp(eta[4:5])))
    eta <- optim(c(qlogis(par0[1]), par0[2:3], log(par0[4:5])), nll,
      method = "BFGS", control = list(reltol = 1e-14, maxit = 1000)
    )$par
    expect_equal(as.vector(par_gem), c(plogis(eta[1]), eta[2:3], exp(eta[4:5])), tolerance = 1e-4)
  })

  test_that(paste0(name, ": GEM agrees with teigen when nu1 = nu2"), {
    skip_if_not_installed("teigen")
    set.seed(1)
    x5 <- rtmix(1000, par_true, c(5, 5))
    model5 <- if (name == "density") tmix_density(x5, c(5, 5)) else tmix_log(x5, c(5, 5))
    par_gem <- fit(model5, eps = 1e-14)
    tg <- teigen::teigen(x5,
      Gs = 2, models = "univUU", dfstart = 5, dfupdate = FALSE,
      scale = FALSE, verbose = FALSE, eps = c(1e-12, 1e-12)
    )
    o <- order(tg$parameters$mean)
    par_tg <- c(tg$parameters$pig[o][1], tg$parameters$mean[o], sqrt(tg$parameters$sigma[1, 1, o]))
    expect_equal(as.vector(par_gem), par_tg, tolerance = 1e-4)
    expect_equal(loglik(model5, par_gem), tg$logl, tolerance = 1e-8)
  })
}

# The two implementations against each other ------------------------------------

test_that("the two implementations compute the same building blocks", {
  d <- models$density
  l <- models$log
  gam <- e_step(d, par0)
  expect_equal(loglik(l, par_rand), loglik(d, par_rand))
  expect_equal(e_step(l, par_rand), e_step(d, par_rand))
  expect_equal(Q(l, par_rand, gam), Q(d, par_rand, gam))
  expect_equal(grad_Q(l, par_rand, gam), grad_Q(d, par_rand, gam))
})

test_that("the two implementations give the same fit", {
  expect_equal(fit(models$log), fit(models$density), tolerance = 1e-8)
})

test_that("the log scale survives an extreme outlier where the density scale fails", {
  x_out <- c(x, -1e60)
  expect_true(is.nan(tail(e_step(tmix_density(x_out, nu), par_true), 1)))
  expect_equal(tail(e_step(tmix_log(x_out, nu), par_true), 1), 1)
})

# The Gaussian mixture used in the outlier comparison ----------------------------

test_that("Gaussian EM: the log-likelihood never decreases and the fit is a stationary point", {
  step <- create_em_gauss_step(x)
  par <- par0
  ll <- loglik_gauss(par, x)
  for (i in 1:50) {
    par <- step(par)
    ll <- c(ll, loglik_gauss(par, x))
  }
  expect_true(all(diff(ll) >= -1e-8))
  expect_equal(numDeriv::grad(function(p) loglik_gauss(p, x), fit_gauss(x, eps = 1e-14)), rep(0, 5), tolerance = 1e-3)
})
