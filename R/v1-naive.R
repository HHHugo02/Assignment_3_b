# Version 1: my first attempts (kept to show the development) -----------------

# GEM with a FIXED step size for the gradient step on Q (no backtracking).
# Problem: no guarantee that Q increases, so this is not a GEM algorithm.
create_gem_step_v1 <- function(x, nu, gamma = 1e-3) {
  force(x)
  force(nu)
  function(par) {
    gam <- e_step(par, x, nu)
    par[1] <- mean(gam)
    g <- grad_Q(par, par, x, nu, gam)
    g[1] <- 0
    par + gamma * g
  }
}

# First version of the augmented EM step: densities from dt(), everything
# recomputed in every iteration, and a data frame to hold intermediate results.
create_em_aug_step_v1 <- function(x, nu) {
  function(par) {
    d <- data.frame(x = x)
    d$f1 <- par[1] * dt((x - par[2]) / par[4], df = nu[1]) / par[4]
    d$f2 <- (1 - par[1]) * dt((x - par[3]) / par[5], df = nu[2]) / par[5]
    d$gam <- d$f1 / (d$f1 + d$f2)
    d$w1 <- (nu[1] + 1) / (nu[1] + ((x - par[2]) / par[4])^2)
    d$w2 <- (nu[2] + 1) / (nu[2] + ((x - par[3]) / par[5])^2)
    mu1 <- sum(d$gam * d$w1 * d$x) / sum(d$gam * d$w1)
    mu2 <- sum((1 - d$gam) * d$w2 * d$x) / sum((1 - d$gam) * d$w2)
    s1 <- sqrt(sum(d$gam * d$w1 * (d$x - mu1)^2) / sum(d$gam))
    s2 <- sqrt(sum((1 - d$gam) * d$w2 * (d$x - mu2)^2) / sum(1 - d$gam))
    c(mean(d$gam), mu1, mu2, s1, s2)
  }
}
