# Two-component mixture of t-distributions: shared basics ----------------------
#
# Parameter vector: par = c(p, mu1, mu2, sigma1, sigma2)
# Shape parameters: nu  = c(nu1, nu2) (fixed, not estimated)
#
# Used by both the function version (tmix-simple.R) and the S3 version
# (tmix-classes.R).

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

# Data-driven starting value: p = 1/2, quartiles for mu, sd(x) / 2 for sigma
default_start <- function(x) {
  q <- quantile(x, c(0.25, 0.75), names = FALSE)
  c(0.5, q[1], q[2], sd(x) / 2, sd(x) / 2)
}
