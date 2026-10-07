# Version 2: the GEM algorithm written against the generics --------------------
#
# Same algorithm as version 1 (gem.R), but it only talks to the model through
# e_step(), Q() and grad_Q(), so it works for every subclass of "tmix".
#
# The backtracking starts at step0 = 2^-5 instead of 1: the accepted steps are
# far below 1 (they scale like 1 / N), so the first halvings were wasted Q calls.

create_gem_step <- function(model, step0 = 2^-5, d = 0.5, c = 0.1) {
  force(model)

  function(par) {
    gam <- e_step(model, par) # E-step
    par[1] <- mean(gam) # M-step for p (closed form)

    g <- grad_Q(model, par, gam)
    g[1] <- 0 # p is already optimal
    Q_old <- Q(model, par, gam)
    step <- step0

    repeat {
      par_new <- par + step * g
      if (par_new[4] > 0 && par_new[5] > 0 &&
        Q(model, par_new, gam) >= Q_old + c * step * sum(g^2)) {
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

fit <- function(model, ...) UseMethod("fit")

# Shared method: fit any "tmix" model with the GEM algorithm
fit.tmix <- function(model, par0 = default_start(model$x), eps = 1e-10, maxit = 1000) {
  em(par0, create_gem_step(model), eps, maxit)
}
