# Version 1: the generalized EM algorithm with plain functions -----------------

# One GEM step: E-step, exact M-step for p, and one gradient step on Q for
# (mu1, mu2, sigma1, sigma2). The step size is found by backtracking: start at
# step0 and halve until Q has increased enough (Armijo condition), so that
# Q(par_new | par) >= Q(par | par) is guaranteed.
create_gem_step <- function(x, nu, step0 = 1, d = 0.5, c = 0.1) {
  force(x)
  force(nu)

  function(par) {
    gam <- e_step(par, x, nu) # E-step
    par[1] <- mean(gam) # M-step for p (closed form)

    g <- grad_Q(par, gam, x, nu)
    g[1] <- 0 # p is already optimal
    Q_old <- Q(par, gam, x, nu)
    step <- step0

    repeat {
      par_new <- par + step * g
      if (par_new[4] > 0 && par_new[5] > 0 &&
        Q(par_new, gam, x, nu) >= Q_old + c * step * sum(g^2)) {
        break
      }
      step <- d * step
      if (step < 1e-14) {
        return(par) # no ascent possible: stationary point
      }
    }
    par_new
  }
}

# Convenience wrapper: fit the mixture with the GEM algorithm
gem <- function(x, nu, par0 = default_start(x), ...) {
  em(par0, create_gem_step(x, nu), ...)
}
