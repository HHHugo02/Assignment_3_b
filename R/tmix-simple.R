# Version 1: the simple implementation as plain functions ----------------------
#
# Everything is computed directly on the density scale with R's dt()
# (dt_ls() and dtmix() are in tmix-common.R).

# Observed-data log-likelihood
loglik <- function(par, x, nu) {
  sum(log(dtmix(x, par, nu)))
}

# E-step: gamma_i = P(Z_i = 1 | X_i = x_i) by Bayes' formula
e_step <- function(par, x, nu) {
  a <- par[1] * dt_ls(x, par[2], par[4], nu[1])
  b <- (1 - par[1]) * dt_ls(x, par[3], par[5], nu[2])
  a / (a + b)
}

# Q(par | par0), where gam = e_step(par0, x, nu) carries everything from par0
Q <- function(par, gam, x, nu) {
  sum(
    gam * log(par[1] * dt_ls(x, par[2], par[4], nu[1])) +
      (1 - gam) * log((1 - par[1]) * dt_ls(x, par[3], par[5], nu[2]))
  )
}

# Gradient of Q w.r.t. par. With r = x - mu and w = (nu + 1) / (nu + r^2 / sigma^2):
#   dQ/dmu    = sum gam * w * r / sigma^2
#   dQ/dsigma = sum gam * (w * r^2 / sigma^2 - 1) / sigma
grad_Q <- function(par, gam, x, nu) {
  p <- par[1]
  s1 <- par[4]
  s2 <- par[5]
  r1 <- x - par[2]
  r2 <- x - par[3]
  w1 <- (nu[1] + 1) / (nu[1] + r1^2 / s1^2)
  w2 <- (nu[2] + 1) / (nu[2] + r2^2 / s2^2)
  c(
    sum(gam) / p - sum(1 - gam) / (1 - p),
    sum(gam * w1 * r1) / s1^2,
    sum((1 - gam) * w2 * r2) / s2^2,
    sum(gam * (w1 * r1^2 / s1^2 - 1)) / s1,
    sum((1 - gam) * (w2 * r2^2 / s2^2 - 1)) / s2
  )
}
