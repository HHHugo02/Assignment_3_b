# Version 2: the building blocks as S3 classes ----------------------------------
#
# A model object holds the data x and the shape parameters nu. The building
# blocks loglik(), e_step(), Q() and grad_Q() are generics, so the GEM
# algorithm (gem-classes.R) works unchanged for every implementation:
#
#   "tmix"              parent class: shared methods (grad_Q, fit, print)
#     "tmix_density"    simple implementation on the density scale with dt()
#     "tmix_log"        log scale: log-densities written out, plogis(), log-sum-exp

# Constructors -------------------------------------------------------------------

new_tmix <- function(x, nu, subclass, ...) {
  stopifnot(is.numeric(x), length(nu) == 2, all(nu > 0))
  structure(list(x = x, nu = nu, ...), class = c(subclass, "tmix"))
}

tmix_density <- function(x, nu) new_tmix(x, nu, "tmix_density")

tmix_log <- function(x, nu) {
  # log of the normalizing constant Gamma((nu + 1) / 2) / (Gamma(nu / 2) sqrt(pi nu)),
  # computed once here instead of in every call
  log_const <- lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(pi * nu)
  new_tmix(x, nu, "tmix_log", log_const = log_const)
}

# Generics -----------------------------------------------------------------------

loglik <- function(model, par) UseMethod("loglik")
e_step <- function(model, par) UseMethod("e_step")
Q <- function(model, par, gam) UseMethod("Q")
grad_Q <- function(model, par, gam) UseMethod("grad_Q")

# Shared methods (parent class "tmix") -------------------------------------------

# The gradient is the same formula for both implementations.
# With r = x - mu and w = (nu + 1) / (nu + r^2 / sigma^2):
#   dQ/dmu    = sum gam * w * r / sigma^2
#   dQ/dsigma = sum gam * (w * r^2 / sigma^2 - 1) / sigma
grad_Q.tmix <- function(model, par, gam) {
  x <- model$x
  nu <- model$nu
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

print.tmix <- function(x, ...) {
  cat(
    "Two-component t-mixture model, implementation: ", class(x)[1], "\n",
    "N = ", length(x$x), ", nu = (", paste(x$nu, collapse = ", "), ")\n",
    sep = ""
  )
  invisible(x)
}

# Density scale: "tmix_density" ----------------------------------------------------

loglik.tmix_density <- function(model, par) {
  sum(log(dtmix(model$x, par, model$nu)))
}

e_step.tmix_density <- function(model, par) {
  x <- model$x
  nu <- model$nu
  a <- par[1] * dt_ls(x, par[2], par[4], nu[1])
  b <- (1 - par[1]) * dt_ls(x, par[3], par[5], nu[2])
  a / (a + b)
}

Q.tmix_density <- function(model, par, gam) {
  x <- model$x
  nu <- model$nu
  sum(
    gam * log(par[1] * dt_ls(x, par[2], par[4], nu[1])) +
      (1 - gam) * log((1 - par[1]) * dt_ls(x, par[3], par[5], nu[2]))
  )
}

# Log scale: "tmix_log" -------------------------------------------------------------

# Log-density of the location-scale t-distribution, written out:
#   log f(x) = log_const - log(sigma) - (nu + 1) / 2 * log(1 + (x - mu)^2 / (nu sigma^2))
ldt_ls <- function(x, mu, sigma, nu, log_const) {
  log_const - log(sigma) - (nu + 1) / 2 * log1p((x - mu)^2 / (nu * sigma^2))
}

# a = log(p f1(x_i)) and b = log((1 - p) f2(x_i))
log_comp <- function(model, par) {
  x <- model$x
  nu <- model$nu
  k <- model$log_const
  list(
    a = log(par[1]) + ldt_ls(x, par[2], par[4], nu[1], k[1]),
    b = log1p(-par[1]) + ldt_ls(x, par[3], par[5], nu[2], k[2])
  )
}

# log(e^a + e^b) with the log-sum-exp trick
loglik.tmix_log <- function(model, par) {
  lc <- log_comp(model, par)
  m <- pmax(lc$a, lc$b)
  sum(m + log(exp(lc$a - m) + exp(lc$b - m)))
}

# gamma = e^a / (e^a + e^b) = 1 / (1 + e^(b - a)) = plogis(a - b)
e_step.tmix_log <- function(model, par) {
  lc <- log_comp(model, par)
  plogis(lc$a - lc$b)
}

# Q only needs the log-densities: no dt() and no log() of a density
Q.tmix_log <- function(model, par, gam) {
  lc <- log_comp(model, par)
  sum(gam * lc$a + (1 - gam) * lc$b)
}
