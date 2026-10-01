# Two-component mixture of t-distributions: the simple implementation ---------
#
# Parameter vector: par = c(p, mu1, mu2, sigma1, sigma2)
# Shape parameters: nu  = c(nu1, nu2) (fixed, not estimated)
#
# Everything is computed directly on the density scale with R's dt().

par_names <- c("p", "mu1", "mu2", "sigma1", "sigma2")

# Density of the location-scale t-distribution
dt_ls <- function(x, mu, sigma, nu) {
  dt((x - mu) / sigma, df = nu) / sigma
}

# Mixture density
dtmix <- function(x, par, nu) {
  par[1] * dt_ls(x, par[2], par[4], nu[1]) +
    (1 - par[1]) * dt_ls(x, par[3], par[5], nu[2])
}

# Simulate N observations together with their component labels z
sim_tmix <- function(N, par, nu) {
  z <- ifelse(runif(N) < par[1], 1, 2)
  x <- ifelse(
    z == 1,
    par[2] + par[4] * rt(N, nu[1]),
    par[3] + par[5] * rt(N, nu[2])
  )
  data.frame(x = x, z = z)
}

rtmix <- function(N, par, nu) sim_tmix(N, par, nu)$x

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

# Data-driven starting value: p = 1/2, quartiles for mu, sd(x) / 2 for sigma
default_start <- function(x) {
  q <- quantile(x, c(0.25, 0.75), names = FALSE)
  c(0.5, q[1], q[2], sd(x) / 2, sd(x) / 2)
}
