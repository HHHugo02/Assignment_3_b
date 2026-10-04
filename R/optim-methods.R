# Alternatives to the GEM algorithm: R's general-purpose optimizers -------------
#
# They maximize the observed log-likelihood directly, for any "tmix" model.
# To avoid constraints the optimizers work in the unconstrained parametrization
#   eta = (logit p, mu1, mu2, log sigma1, log sigma2).

eta_to_par <- function(eta) c(plogis(eta[1]), eta[2:3], exp(eta[4:5]))
par_to_eta <- function(par) c(qlogis(par[1]), par[2:3], log(par[4:5]))

# Negative log-likelihood and its gradient as functions of eta.
# The gradient of the log-likelihood comes from Fisher's identity,
# grad loglik(par) = grad_Q(par | par), and the chain rule
# d par / d eta = (p (1 - p), 1, 1, sigma1, sigma2).
neg_loglik_eta <- function(model) {
  force(model)
  function(eta) -loglik(model, eta_to_par(eta))
}

neg_score_eta <- function(model) {
  force(model)
  function(eta) {
    par <- eta_to_par(eta)
    -grad_Q(model, par, e_step(model, par)) *
      c(par[1] * (1 - par[1]), 1, 1, par[4], par[5])
  }
}

# Fit with optim(); Nelder-Mead does not use the gradient
fit_optim <- function(model, method = "BFGS", par0 = default_start(model$x),
                      reltol = 1e-10, maxit = 5000) {
  gr <- if (method == "Nelder-Mead") NULL else neg_score_eta(model)
  control <- if (method == "L-BFGS-B") {
    list(factr = reltol / .Machine$double.eps, maxit = maxit)
  } else {
    list(reltol = reltol, maxit = maxit)
  }
  opt <- optim(par_to_eta(par0), neg_loglik_eta(model), gr, method = method, control = control)
  par <- setNames(eta_to_par(opt$par), par_names)
  attr(par, "iter") <- unname(opt$counts[1]) # number of function evaluations
  par
}

# Fit with nlminb() (a quasi-Newton method from the PORT library)
fit_nlminb <- function(model, par0 = default_start(model$x), rel.tol = 1e-10) {
  opt <- nlminb(par_to_eta(par0), neg_loglik_eta(model), neg_score_eta(model),
    control = list(rel.tol = rel.tol, eval.max = 5000, iter.max = 5000)
  )
  par <- setNames(eta_to_par(opt$par), par_names)
  attr(par, "iter") <- unname(opt$evaluations[1])
  par
}
