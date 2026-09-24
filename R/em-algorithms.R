# Generic fixed-point driver and step-function factories ----------------------
#
# All algorithms below are written as a map  par -> new par  ("step function")
# and run by the same driver `em()`. This makes it easy to use *exactly* the
# same convergence criterion for every algorithm.

# Convergence criteria (function factories) -----------------------------------

# Criterion from the lecture slides: small relative change in the parameters
par_conv <- function(eps = 1e-10) {
  force(eps)
  function(par, par0) sum((par - par0)^2) <= eps * (sum(par^2) + eps)
}

# Criterion on the gradient of the marginal log-likelihood (per observation).
# Thanks to Fisher's identity it only costs one extra E-step.
grad_conv <- function(x, nu, tol = 1e-6) {
  force(x)
  force(nu)
  force(tol)
  N <- length(x)
  function(par, par0) sqrt(sum(grad_loglik(par, x, nu)^2)) / N <= tol
}

# The driver ------------------------------------------------------------------

em <- function(par, em_step, eps = 1e-10, maxit = 5000, cb = NULL,
               conv = par_conv(eps)) {
  for (i in seq_len(maxit)) {
    par0 <- par
    par <- em_step(par)

    if (!is.null(cb)) {
      cb()
    }

    if (conv(par, par0)) {
      break
    }
  }
  par <- setNames(as.vector(par), par_names)
  attr(par, "iter") <- i
  par
}

# 1. Generalized EM: exact p-update + ONE gradient step on Q ------------------
#
# The M-step for (mu, sigma) has no closed form. A GEM step only has to
# *increase* Q, so we take one gradient ascent step with backtracking
# (Armijo condition) on Q( . | par0).
#
#  warm = FALSE: backtracking starts from gamma0 in every iteration.
#  warm = TRUE : backtracking starts from the last accepted step size / d
#                (found after profiling; see the presentation).

create_gem_step <- function(x, nu, gamma0 = 1, d = 0.5, c = 0.1, warm = FALSE) {
  force(x)
  force(nu)
  gamma_last <- gamma0

  function(par) {
    gam <- e_step(par, x, nu) # E-step
    par[1] <- mean(gam) # exact M-step for p
    Q_val <- Q(par, par, x, nu, gam)
    g <- grad_Q(par, par, x, nu, gam)
    g[1] <- 0 # p is already optimal
    gamma <- if (warm) gamma_last / d else gamma0

    repeat {
      par_new <- par + gamma * g
      if (Q(par_new, par, x, nu, gam) >= Q_val + c * gamma * sum(g^2)) {
        break
      }
      gamma <- d * gamma
      if (gamma < 1e-14) {
        return(par) # no ascent possible: we are at a stationary point
      }
    }
    gamma_last <<- gamma
    par_new
  }
}

# 2. EM for the augmented complete data (X, Z, W): closed-form M-step --------
#
# E-step: gamma_i and w_ij = E(W | X = x_i, Z = j) = (nu_j + 1) / (nu_j + d_ij)
# M-step: weighted means and variances.

m_step_aug <- function(e, x) {
  a1 <- e$gam * e$w1
  a2 <- (1 - e$gam) * e$w2
  mu1 <- sum(a1 * x) / sum(a1)
  mu2 <- sum(a2 * x) / sum(a2)
  N1 <- sum(e$gam)
  N2 <- length(x) - N1
  c(
    N1 / length(x),
    mu1,
    mu2,
    sqrt(sum(a1 * (x - mu1)^2) / N1),
    sqrt(sum(a2 * (x - mu2)^2) / N2)
  )
}

create_em_aug_step <- function(x, nu) {
  force(x)
  force(nu)
  function(par) m_step_aug(e_step_aug(par, x, nu), x)
}

# Same algorithm after profiling (plogis, cbind and repeated lgamma/log calls
# were the hot spots): constants precomputed, the standardized squared
# residuals d are reused for both the E-step and the weights, and only one
# exp() per observation is needed for the responsibilities.
create_em_aug_step_fast <- function(x, nu) {
  force(x)
  N <- length(x)
  lk <- lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(pi * nu)
  h <- (nu + 1) / 2
  nu1 <- nu[1]
  nu2 <- nu[2]

  function(par) {
    r1 <- x - par[2]
    r2 <- x - par[3]
    d1 <- r1 * r1 / (par[4] * par[4])
    d2 <- r2 * r2 / (par[5] * par[5])
    # log(weighted f2) - log(weighted f1)
    const <- log1p(-par[1]) + lk[2] - log(par[5]) - log(par[1]) - lk[1] + log(par[4])
    gam <- 1 / (1 + exp(const - h[2] * log1p(d2 / nu2) + h[1] * log1p(d1 / nu1)))
    a1 <- gam * (nu1 + 1) / (nu1 + d1)
    a2 <- (1 - gam) * (nu2 + 1) / (nu2 + d2)
    mu1 <- sum(a1 * x) / sum(a1)
    mu2 <- sum(a2 * x) / sum(a2)
    N1 <- sum(gam)
    c(
      N1 / N, mu1, mu2,
      sqrt(sum(a1 * (x - mu1)^2) / N1),
      sqrt(sum(a2 * (x - mu2)^2) / (N - N1))
    )
  }
}

# 3. "Exact" EM with only Z missing: M-step solved numerically by BFGS -------

create_em_exact_step <- function(x, nu, reltol = 1e-12) {
  force(x)
  force(nu)

  # maximize sum_i g_i log f(x_i | mu, sigma, nu_j) over (mu, log sigma)
  weighted_t_mle <- function(mu, sigma, g, nu_j) {
    fn <- function(eta) -sum(g * ldt_ls(x, eta[1], exp(eta[2]), nu_j))
    gr <- function(eta) {
      s <- exp(eta[2])
      r <- x - eta[1]
      w <- (nu_j + 1) / (nu_j + r^2 / s^2)
      -c(sum(g * w * r) / s^2, sum(g * (w * r^2 / s^2 - 1)))
    }
    eta <- optim(c(mu, log(sigma)), fn, gr,
      method = "BFGS",
      control = list(reltol = reltol)
    )$par
    c(eta[1], exp(eta[2]))
  }

  function(par) {
    gam <- e_step(par, x, nu)
    c1 <- weighted_t_mle(par[2], par[4], gam, nu[1])
    c2 <- weighted_t_mle(par[3], par[5], 1 - gam, nu[2])
    c(mean(gam), c1[1], c2[1], c1[2], c2[2])
  }
}

# 4. EM for a two-component *Gaussian* mixture (for the robustness study) ----

loglik_gauss <- function(par, x) {
  sum(log(
    par[1] * dnorm(x, par[2], par[4]) + (1 - par[1]) * dnorm(x, par[3], par[5])
  ))
}

create_em_gauss_step <- function(x) {
  force(x)
  function(par) {
    a <- log(par[1]) + dnorm(x, par[2], par[4], log = TRUE)
    b <- log1p(-par[1]) + dnorm(x, par[3], par[5], log = TRUE)
    gam <- plogis(a - b)
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
