# EM driver from the lecture, with the same stopping rule:
# stop when ||par - par0||^2 <= eps (||par||^2 + eps)
#
# `em_step` is any function par -> new par, so the driver is shared by the
# function version and the S3 version.
em <- function(par, em_step, eps = 1e-10, maxit = 1000) {
  for (i in seq_len(maxit)) {
    par0 <- par
    par <- em_step(par)

    if (sum((par - par0)^2) <= eps * (sum(par^2) + eps)) {
      break
    }
  }
  par <- setNames(as.vector(par), par_names)
  attr(par, "iter") <- i
  par
}
