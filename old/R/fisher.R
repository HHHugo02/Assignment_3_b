# Observed Fisher information: three ways (lecture 10) ------------------------
#
# Note: numDeriv's functions have an argument called `x`, so the data are passed
# via closures instead of through `...`.

# 1. Gradient identity: differentiate the score grad_Q(par | par) numerically
fisher_score <- function(par, x, nu) {
  -numDeriv::jacobian(function(p) grad_loglik(p, x, nu), par)
}

# 2. Information identity (Louis): i_X = i_Y - i_{Y|X}, complete data (X, Z)
fisher_louis <- function(par, x, nu) {
  iY <- -numDeriv::jacobian(function(p) grad_Q(p, par, x, nu), par)
  iYX <- numDeriv::jacobian(function(p0) grad_Q(par, p0, x, nu), par)
  list(iY = iY, iYX = iYX, iX = iY - iYX)
}

# 3. EM map: i_X = (I - D Phi^T) i_Y, here with the augmented EM (X, Z, W)
fisher_em_map <- function(par, x, nu) {
  Phi <- create_em_aug_step(x, nu)
  DPhi <- numDeriv::jacobian(Phi, par)
  iY <- -numDeriv::hessian(function(p) Q_aug(p, par, x, nu), par)
  list(DPhi = DPhi, iY = iY, iX = (diag(length(par)) - t(DPhi)) %*% iY)
}

# Rate of linear convergence = spectral radius of D Phi = i_{Y|X} i_Y^{-1}
spectral_radius <- function(A) max(Mod(eigen(A, only.values = TRUE)$values))
