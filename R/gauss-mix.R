# Two-component Gaussian mixture, for comparison with the t-mixture ------------
#
# Same parameter vector par = c(p, mu1, mu2, sigma1, sigma2). Here the M-step
# has a closed form (weighted means and standard deviations), so this is an
# ordinary EM algorithm, run with the same driver em().

loglik_gauss <- function(par, x) {
  sum(log(
    par[1] * dnorm(x, par[2], par[4]) + (1 - par[1]) * dnorm(x, par[3], par[5])
  ))
}

create_em_gauss_step <- function(x) {
  force(x)

  function(par) {
    # E-step on log scale: gam = a / (a + b) = plogis(log a - log b)
    a <- log(par[1]) + dnorm(x, par[2], par[4], log = TRUE)
    b <- log1p(-par[1]) + dnorm(x, par[3], par[5], log = TRUE)
    gam <- plogis(a - b)

    # M-step (closed form)
    N1 <- sum(gam)
    N2 <- length(x) - N1
    mu1 <- sum(gam * x) / N1
    mu2 <- sum((1 - gam) * x) / N2
    c(
      N1 / length(x), mu1, mu2,
      sqrt(sum(gam * (x - mu1)^2) / N1),
      sqrt(sum((1 - gam) * (x - mu2)^2) / N2)
    )
  }
}

# Convenience wrapper: fit the Gaussian mixture with the EM algorithm
fit_gauss <- function(x, par0 = default_start(x), eps = 1e-10, maxit = 1000) {
  em(par0, create_em_gauss_step(x), eps, maxit)
}
