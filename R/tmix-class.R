# A small S3 interface on top of the algorithms --------------------------------

tmix_methods <- c("em_aug", "em_aug_fast", "em_aug_cpp", "gem", "em_exact", "gd", "bfgs")

tmix <- function(x, nu = c(5, 5), par0 = default_start(x),
                 method = "em_aug", eps = 1e-10, maxit = 5000, conv = NULL) {
  method <- match.arg(method, tmix_methods)
  if (is.null(conv)) {
    conv <- par_conv(eps)
  }
  time <- system.time({
    par <- if (method == "bfgs") {
      fit_bfgs(par0, x, nu)
    } else {
      step <- switch(method,
        em_aug      = create_em_aug_step(x, nu),
        em_aug_fast = create_em_aug_step_fast(x, nu),
        em_aug_cpp  = function(par) em_aug_step_cpp(x, par, nu),
        gem         = create_gem_step(x, nu, warm = TRUE),
        em_exact    = create_em_exact_step(x, nu),
        gd          = create_gd_step(x, nu)
      )
      em(par0, step, maxit = maxit, conv = conv)
    }
  })[["elapsed"]]

  structure(
    list(
      par = setNames(as.vector(par), par_names),
      iter = attr(par, "iter"),
      loglik = loglik(par, x, nu),
      method = method, nu = nu, x = x, time = time
    ),
    class = "tmix"
  )
}

coef.tmix <- function(object, ...) object$par

logLik.tmix <- function(object, ...) {
  structure(object$loglik, df = 5, nobs = length(object$x), class = "logLik")
}

vcov.tmix <- function(object, ...) {
  solve(fisher_score(object$par, object$x, object$nu))
}

predict.tmix <- function(object, newdata = object$x, ...) {
  e_step(object$par, newdata, object$nu) # P(Z = 1 | x)
}

print.tmix <- function(x, ...) {
  se <- sqrt(diag(vcov(x)))
  cat("Two-component t-mixture, nu = (", paste(x$nu, collapse = ", "), ")\n", sep = "")
  cat("Method:", x$method, "| iterations:", x$iter,
    "| loglik:", format(x$loglik, digits = 8), "\n\n")
  print(rbind(estimate = x$par, se = se), digits = 4)
  invisible(x)
}

plot.tmix <- function(x, ...) {
  grid <- seq(min(x$x), max(x$x), length.out = 400)
  p <- x$par
  df_fit <- data.frame(
    x = rep(grid, 3),
    density = c(
      dtmix(grid, p, x$nu),
      p[1] * exp(ldt_ls(grid, p[2], p[4], x$nu[1])),
      (1 - p[1]) * exp(ldt_ls(grid, p[3], p[5], x$nu[2]))
    ),
    curve = rep(c("mixture", "component 1", "component 2"), each = length(grid))
  )
  ggplot2::ggplot(data.frame(x = x$x), ggplot2::aes(x)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 40, fill = "grey85", colour = "white"
    ) +
    ggplot2::geom_line(
      data = df_fit,
      ggplot2::aes(x, density, colour = curve, linetype = curve),
      linewidth = 0.9
    ) +
    ggplot2::scale_linetype_manual(values = c(2, 2, 1)) +
    ggplot2::labs(colour = NULL, linetype = NULL, y = "density")
}
