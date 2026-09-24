# Alternatives: optimize the marginal log-likelihood directly ------------------

# Gradient ascent on loglik with backtracking, written as a step function so it
# can be run by the same driver `em()` (and the same convergence criterion).
# The gradient is computed with Fisher's identity: grad_Q(par | par).
create_gd_step <- function(x, nu, gamma0 = 1, d = 0.5, c = 0.1) {
  force(x)
  force(nu)
  gamma_last <- gamma0

  function(par) {
    val <- loglik(par, x, nu)
    g <- grad_loglik(par, x, nu)
    gamma <- gamma_last / d

    repeat {
      par_new <- par + gamma * g
      if (loglik(par_new, x, nu) >= val + c * gamma * sum(g^2)) {
        break
      }
      gamma <- d * gamma
      if (gamma < 1e-14) {
        return(par)
      }
    }
    gamma_last <<- gamma
    par_new
  }
}

# Quasi-Newton (BFGS via optim) in the unconstrained parametrization
#   eta = (logit p, mu1, mu2, log sigma1, log sigma2).
eta_to_par <- function(eta) c(plogis(eta[1]), eta[2:3], exp(eta[4:5]))
par_to_eta <- function(par) c(qlogis(par[1]), par[2:3], log(par[4:5]))

fit_bfgs <- function(par, x, nu, reltol = 1e-12, maxit = 1000, cb = NULL) {
  fn <- function(eta) {
    par <- eta_to_par(eta)
    if (!is.null(cb)) {
      cb() # the tracer picks up `par` from this frame
    }
    -loglik(par, x, nu)
  }
  gr <- function(eta) {
    par <- eta_to_par(eta)
    # chain rule: d par / d eta = (p(1 - p), 1, 1, sigma1, sigma2)
    -grad_loglik(par, x, nu) * c(par[1] * (1 - par[1]), 1, 1, par[4], par[5])
  }
  opt <- optim(par_to_eta(par), fn, gr,
    method = "BFGS",
    control = list(reltol = reltol, maxit = maxit)
  )
  par <- setNames(eta_to_par(opt$par), par_names)
  attr(par, "iter") <- opt$counts[["gradient"]]
  par
}
