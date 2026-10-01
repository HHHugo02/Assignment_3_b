# Two-component mixture of t-distributions ------------------------------------
#
# Parameter vector: par = c(p, mu1, mu2, sigma1, sigma2)
# Shape parameters: nu  = c(nu1, nu2) (fixed, not estimated)

par_names <- c("p", "mu1", "mu2", "sigma1", "sigma2")

valid_par <- function(par) {
  par[1] > 0 && par[1] < 1 && par[4] > 0 && par[5] > 0
}

# Log-density of the location-scale t-distribution
ldt_ls <- function(x, mu, sigma, nu) {
  lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(pi * nu) - log(sigma) -
    (nu + 1) / 2 * log1p((x - mu)^2 / (nu * sigma^2))
}

dtmix <- function(x, par, nu) {
  par[1] * exp(ldt_ls(x, par[2], par[4], nu[1])) +
    (1 - par[1]) * exp(ldt_ls(x, par[3], par[5], nu[2]))
}

rtmix <- function(N, par, nu) {
  z <- rbinom(N, 1, par[1]) == 1
  x <- numeric(N)
  x[z] <- par[2] + par[4] * rt(sum(z), nu[1])
  x[!z] <- par[3] + par[5] * rt(sum(!z), nu[2])
  x
}

# Log of the two weighted component densities, an N x 2 matrix
log_comp <- function(par, x, nu) {
  cbind(
    log(par[1]) + ldt_ls(x, par[2], par[4], nu[1]),
    log1p(-par[1]) + ldt_ls(x, par[3], par[5], nu[2])
  )
}

# Observed-data (marginal) log-likelihood, computed with log-sum-exp
loglik <- function(par, x, nu) {
  if (!valid_par(par)) {
    return(-Inf)
  }
  lc <- log_comp(par, x, nu)
  m <- pmax(lc[, 1], lc[, 2])
  sum(m + log(exp(lc[, 1] - m) + exp(lc[, 2] - m)))
}

# E-step: responsibilities gamma_i = P(Z_i = 1 | X_i = x_i, par)
e_step <- function(par, x, nu) {
  lc <- log_comp(par, x, nu)
  plogis(lc[, 1] - lc[, 2])
}

# The Q-function with only the component labels Z as missing data:
#   Q(par | par0) = sum_i gamma_i (log p + log f1(x_i))
#                   + (1 - gamma_i) (log(1 - p) + log f2(x_i)),
# where gamma_i are computed at par0.
Q <- function(par, par0, x, nu, gam = e_step(par0, x, nu)) {
  if (!valid_par(par)) {
    return(-Inf)
  }
  sum(
    gam * (log(par[1]) + ldt_ls(x, par[2], par[4], nu[1])) +
      (1 - gam) * (log1p(-par[1]) + ldt_ls(x, par[3], par[5], nu[2]))
  )
}

# Gradient of Q w.r.t. the first argument, par = (p, mu1, mu2, sigma1, sigma2).
# With r = x - mu and w = (nu + 1) / (nu + r^2 / sigma^2):
#   dQ/dmu    = sum gamma * w * r / sigma^2
#   dQ/dsigma = sum gamma * (w * r^2 / sigma^2 - 1) / sigma
grad_Q <- function(par, par0, x, nu, gam = e_step(par0, x, nu)) {
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

# Fisher's identity: grad loglik(par) = grad_1 Q(par | par)
grad_loglik <- function(par, x, nu) grad_Q(par, par, x, nu)

# Q-function for the *augmented* complete data (X, Z, W), where
# W | Z = j ~ Gamma(nu_j / 2, rate = nu_j / 2) and X | W, Z = j ~ N(mu_j, sigma_j^2 / W).
# Terms that do not depend on par are dropped.
e_step_aug <- function(par, x, nu) {
  gam <- e_step(par, x, nu)
  w1 <- (nu[1] + 1) / (nu[1] + (x - par[2])^2 / par[4]^2)
  w2 <- (nu[2] + 1) / (nu[2] + (x - par[3])^2 / par[5]^2)
  list(gam = gam, w1 = w1, w2 = w2)
}

Q_aug <- function(par, par0, x, nu, e = e_step_aug(par0, x, nu)) {
  if (!valid_par(par)) {
    return(-Inf)
  }
  sum(
    e$gam * (log(par[1]) - log(par[4]) - e$w1 * (x - par[2])^2 / (2 * par[4]^2)) +
      (1 - e$gam) * (log1p(-par[1]) - log(par[5]) - e$w2 * (x - par[3])^2 / (2 * par[5]^2))
  )
}

# Sensible data-driven starting value
default_start <- function(x) {
  q <- quantile(x, c(0.25, 0.75), names = FALSE)
  c(0.5, q[1], q[2], sd(x) / 2, sd(x) / 2)
}
