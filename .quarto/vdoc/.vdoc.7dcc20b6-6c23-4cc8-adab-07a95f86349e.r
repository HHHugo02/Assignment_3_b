#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: setup
#| include: false

library(ggplot2)
library(testthat)
source("R/load-v1.R")
theme_set(theme_minimal(base_size = 15))
options(digits = 5, width = 110)

# print the source of a function (with comments) in a code chunk
show_fun <- function(name) {
  paste0(name, " <- ", paste(as.character(attr(get(name), "srcref")), collapse = "\n"))
}
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: sim-data
set.seed(42)
nu <- c(5, 10)
par_true <- c(p = 0.6, mu1 = 0, mu2 = 3,
              sigma1 = 1, sigma2 = 0.7)
x <- rtmix(400, par_true, nu)
par0 <- default_start(x)
round(par0, 3)
#
#
#
#
#
#
#| label: plot-data
#| echo: false
#| fig-width: 6.5
#| fig-height: 4.3
grid <- seq(min(x), max(x), length.out = 400)
ggplot(data.frame(x), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey85", colour = "white") +
  geom_line(data = data.frame(x = grid, d = dtmix(grid, par_true, nu)), aes(x, d), linewidth = 1) +
  labs(y = "density", title = "N = 400, true density")
```
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: fixed-steps
#| echo: false
gam <- e_step(par0, x, nu)
par1 <- par0
par1[1] <- mean(gam)                       # exact M-step for p
g <- grad_Q(par1, gam, x, nu)
g[1] <- 0
change_in_Q <- function(s) Q(par1 + s * g, gam, x, nu) - Q(par1, gam, x, nu)
#
#
#
#| label: plot-fixed-steps
#| echo: false
#| fig-height: 3.2
step_df <- data.frame(s = seq(0, 0.0085, length.out = 300))
step_df$change <- sapply(step_df$s, change_in_Q)
s_best <- optimize(change_in_Q, c(0, 0.0085), maximum = TRUE)$maximum
s_zero <- uniroot(change_in_Q, c(s_best, 0.0085))$root
ggplot(step_df, aes(s, change)) +
  annotate("rect", xmin = s_zero, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#d95f02", alpha = 0.12) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(linewidth = 1.1, colour = "#2c7fb8") +
  geom_vline(xintercept = s_best, linetype = 2, colour = "grey40") +
  annotate("point", x = s_best, y = change_in_Q(s_best), size = 3) +
  annotate("text", x = s_best, y = change_in_Q(s_best) + 0.9, label = "best step", size = 4.5) +
  annotate("text", x = 0.0004, y = -2.2, label = "too short:\nbarely moves", size = 4.5, hjust = 0) +
  annotate("text", x = s_zero + 0.0002, y = 2.2, label = "too long:\nQ decreases", size = 4.5, hjust = 0, colour = "#d95f02") +
  labs(x = "step size s", y = "change in Q")
#
#
#
#
#
#
#
#
#
#
#
#| label: run-gem
fit <- gem(x, nu)
fit
loglik(fit, x, nu)
```
#
#
#
#| label: plot-fit
#| echo: false
#| fig-width: 6.5
#| fig-height: 4.3
p <- fit
df_fit <- data.frame(
  x = rep(grid, 4),
  density = c(dtmix(grid, p, nu), p[1] * dt_ls(grid, p[2], p[4], nu[1]), (1 - p[1]) * dt_ls(grid, p[3], p[5], nu[2]),
              dtmix(grid, par_true, nu)),
  curve = rep(c("mixture", "component 1", "component 2", "true density"), each = length(grid))
)
ggplot(data.frame(x), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey85", colour = "white") +
  geom_line(data = df_fit, aes(x, density, colour = curve, linetype = curve), linewidth = 0.9) +
  scale_colour_manual(values = c("#F8766D", "#00BA38", "#619CFF", "black")) +
  scale_linetype_manual(values = c(2, 2, 1, 3)) +
  labs(colour = NULL, linetype = NULL, y = "density", title = "Fitted mixture and true density")
```
#
#
#
#
#
#
#
#
#
#
#
#
#| label: test-optim
#| echo: false
par_hat <- gem(x, nu, eps = 1e-14)
nll <- function(eta) -loglik(c(plogis(eta[1]), eta[2:3], exp(eta[4:5])), x, nu)
eta <- optim(c(qlogis(par0[1]), par0[2:3], log(par0[4:5])), nll,
             method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))$par
par_optim <- c(plogis(eta[1]), eta[2:3], exp(eta[4:5]))

optim_tab <- rbind(
  GEM = c(par_hat, loglik(par_hat, x, nu)),
  `optim() (BFGS)` = c(par_optim, loglik(par_optim, x, nu))
)
optim_fmt <- rbind(
  formatC(optim_tab, format = "f", digits = 4),
  Difference = formatC(abs(optim_tab[1, ] - optim_tab[2, ]), format = "e", digits = 1)
)
knitr::kable(
  optim_fmt,
  col.names = c("$p$", "$\\mu_1$", "$\\mu_2$", "$\\sigma_1$", "$\\sigma_2$", "$\\ell(\\hat\\theta)$"),
  align = "r"
)
#
#
#
#
#
#
#
#
#
#| label: prof-setup
#| include: false
set.seed(7)
x_big <- rtmix(5000, par_true, nu)
fit_big <- gem(x_big, nu)
n_iter <- attr(fit_big, "iter")
bm_fit <- bench::mark(gem(x_big, nu), iterations = 5)

prof_file <- profile_expr(for (i in 1:5) gem(x_big, nu))
prof <- summaryRprof(prof_file)
#
#
#
#
#
#| label: prof-calls
#| echo: false
calls <- count_calls(gem(x_big, nu), c("e_step", "Q", "grad_Q", "dt_ls"))

gam_big <- e_step(fit_big, x_big, nu)
cost <- bench::mark(
  e_step = e_step(fit_big, x_big, nu),
  Q      = Q(fit_big, gam_big, x_big, nu),
  grad_Q = grad_Q(fit_big, gam_big, x_big, nu),
  dt_ls  = dt_ls(x_big, 0, 1, 5),
  check = FALSE
)
per_call <- setNames(as.numeric(cost$median), as.character(cost$expression))[names(calls)]
total_ms <- 1000 * calls * per_call
calls_tab <- data.frame(
  fun   = paste0("`", names(calls), "()`"),
  calls = formatC(calls, format = "d"),
  per_iter = formatC(calls / n_iter, format = "f", digits = 1),
  ms_call  = formatC(1000 * per_call, format = "f", digits = 2),
  total    = formatC(total_ms, format = "f", digits = 0),
  share    = paste0(round(100 * total_ms / (1000 * as.numeric(bm_fit$median))), "%")
)
calls_tab <- calls_tab[c(2, 1, 3, 4), ]
calls_tab[1, ] <- paste0("**", calls_tab[1, ], "**")
calls_tab$fun[4] <- paste0(calls_tab$fun[4], " *(inside `Q()` and `e_step()`)*")
knitr::kable(
  calls_tab, row.names = FALSE, align = c("l", "r", "r", "r", "r", "r"),
  col.names = c("Function", "Calls per fit", "Calls per iteration", "Time per call (ms)",
                "Total time (ms)", "Share of the fit")
)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: start-points
#| echo: false
# Replay the GEM iterations and recover the accepted step size s in each one:
# par_new = par + s * g, so s = <par_new - par, g> / ||g||^2.
step_sizes <- function(x, nu, par = default_start(x), eps = 1e-10, maxit = 1000) {
  gem_step <- create_gem_step(x, nu)
  out <- list()
  for (i in seq_len(maxit)) {
    gam <- e_step(par, x, nu)
    par1 <- par
    par1[1] <- mean(gam)
    g <- grad_Q(par1, gam, x, nu)
    g[1] <- 0
    par_new <- gem_step(par)
    out[[i]] <- data.frame(iter = i, step = sum((par_new - par1) * g) / sum(g^2), grad = sqrt(sum(g^2)),
                           mu1 = par_new[2], mu2 = par_new[3])
    done <- sum((par_new - par)^2) <= eps * (sum(par_new^2) + eps)
    par <- par_new
    if (done) break
  }
  do.call(rbind, out)
}

# 50 starting points: the default one and 49 random ones
set.seed(1)
starts <- c(list(par0), replicate(
  49, c(runif(1, 0.1, 0.9), runif(2, min(x), max(x)), sd(x) * runif(2, 0.2, 2)), simplify = FALSE
))

# Where a fit ends: at the maximum found from the default start, with p at 0 or 1, or at the mirror image
end_levels <- c("global maximum", "labels swapped", "one component")
end_of <- function(est) {
  if (loglik(par_hat, x, nu) - loglik(est, x, nu) < 1e-3) end_levels[1]
  else if (est[1] < 0.01 || est[1] > 0.99) end_levels[3]
  else end_levels[2]
}

# Fit from every starting point with the backtracking started at step0
fit_starts <- function(step0) {
  do.call(rbind, lapply(starts, function(par) {
    n_Q <- count_calls(est <- em(par, create_gem_step(x, nu, step0 = step0)), "Q")
    data.frame(n_Q = n_Q[["Q"]], iter = attr(est, "iter"), end = end_of(est), loglik = loglik(est, x, nu))
  }))
}
runs <- lapply(list(`1` = 1, `2^-5` = 2^-5, `2^-8` = 2^-8), fit_starts)
n_end <- table(factor(runs[["1"]]$end, levels = end_levels))
ll_swapped <- runs[["1"]]$loglik[runs[["1"]]$end == "labels swapped"]
stopifnot(diff(range(ll_swapped)) < 1e-3) # the swapped fits are one and the same local maximum

paths <- do.call(rbind, lapply(seq_along(starts), function(k) {
  rbind(data.frame(start = k, iter = 0, step = NA, grad = NA, mu1 = starts[[k]][2], mu2 = starts[[k]][3]),
        cbind(start = k, step_sizes(x, nu, starts[[k]])))
}))
paths$end <- factor(runs[["1"]]$end[paths$start], levels = end_levels)
accepted <- paths[paths$iter > 0, ]
accepted$log2_step <- round(log2(accepted$step))
n_long <- sum(accepted$log2_step > -5)
end_cols <- setNames(c("#2c7fb8", "#d95f02", "#1b9e77"), end_levels)
#
#
#
#
#
#
#
#| label: plot-start-paths
#| echo: false
#| fig-width: 6
#| fig-height: 3.2
ggplot(paths, aes(mu1, mu2, group = start, colour = end)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 3) +
  geom_path(linewidth = 0.4, alpha = 0.7) +
  geom_point(data = paths[paths$iter == 0, ], shape = 21, fill = "white", size = 2) +
  geom_path(data = paths[paths$start == 1, ], colour = "black", linewidth = 1) +
  geom_point(data = paths[paths$start == 1 & paths$iter == 0, ], colour = "black", size = 3) +
  annotate("point", x = par_hat[2], y = par_hat[3], shape = 4, size = 4, stroke = 1.5) +
  annotate("segment", x = -2.3, y = 5.5, xend = par0[2] - 0.12, yend = par0[3] + 0.15, colour = "grey30") +
  annotate("text", x = -4.4, y = 5.7, label = "default start", size = 4.5, hjust = 0) +
  annotate("text", x = 5.8, y = 3.9, label = "mu[1] == mu[2]", parse = TRUE, size = 4.5, hjust = 1, colour = "grey40") +
  scale_x_continuous(breaks = seq(-4, 6, 2)) +
  scale_colour_manual(values = end_cols) +
  labs(x = expression(mu[1]), y = expression(mu[2]), colour = "ends in") +
  theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: plot-start-steps
#| echo: false
#| fig-width: 6
#| fig-height: 3.2
ggplot(accepted, aes(log2_step, fill = end)) +
  geom_bar(width = 0.8, colour = "white", linewidth = 0.4) +
  geom_vline(xintercept = -4.5, linetype = 2) +
  annotate("text", x = -4.3, y = Inf, vjust = 1.5, hjust = 0, size = 4.5, parse = TRUE, label = "s > 2^-5") +
  annotate("text", x = -4.3, y = Inf, vjust = 4, hjust = 0, size = 4.5,
           label = paste(n_long, "of", nrow(accepted))) +
  scale_x_continuous(breaks = -12:0, labels = function(k) parse(text = paste0("2^", k))) +
  scale_fill_manual(values = end_cols) +
  labs(x = "accepted step size s", y = "iterations", fill = "ends in") +
  theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        legend.key.size = unit(0.9, "lines"))
```
#
#
#
#
#
#| label: tab-start-steps
#| echo: false
start_tab <- data.frame(
  step0 = c("$1$", "$2^{-5}$", "$2^{-8}$"),
  n_Q   = formatC(sapply(runs, function(r) sum(r$n_Q)), format = "d"),
  iter  = formatC(sapply(runs, function(r) sum(r$iter)), format = "d")
)
start_tab[2, ] <- paste0("**", start_tab[2, ], "**")
knitr::kable(
  start_tab, row.names = FALSE, align = c("l", "r", "r"),
  col.names = c("First step", "Evaluations of $Q$", "Iterations")
)
```
#
#
#
#
#
#
#| label: load-s3
#| include: false
fit_v1 <- fit # the version 1 estimate, before fit() becomes a generic
# Version 2 of the code: the building blocks become S3 generics (replaces the plain functions)
source("R/load.R")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: s3-fit
model_d <- tmix_density(x, nu)   # density scale, dt()
model_l <- tmix_log(x, nu)       # log scale

rbind(density = fit(model_d),
      log     = fit(model_l))
#
#
#
#
#
#
#
#
#| label: s3-agree
par_test <- c(0.4, -0.5, 2.5, 1.3, 0.9)
gam_test <- e_step(model_d, par0)

test_that("the log class computes the same as the density class", {
  expect_equal(loglik(model_l, par_test), loglik(model_d, par_test))
  expect_equal(e_step(model_l, par_test), e_step(model_d, par_test))
  expect_equal(Q(model_l, par_test, gam_test), Q(model_d, par_test, gam_test))
  expect_equal(fit(model_l), fit(model_d), tolerance = 1e-8)
})
#
#
#
#
#
#
#
#| label: s3-bench
#| echo: false
model_density <- tmix_density(x_big, nu)
model_log <- tmix_log(x_big, nu)
same_fit <- function(a, b) isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = 1e-8))
bm_classes <- bench::mark(density = fit(model_density), log = fit(model_log),
                          check = same_fit, iterations = 10)

par_b <- fit(model_log)
gam_b <- e_step(model_log, par_b)
blocks <- bench::mark(
  Q_density = Q(model_density, par_b, gam_b),           Q_log = Q(model_log, par_b, gam_b),
  e_step_density = e_step(model_density, par_b),        e_step_log = e_step(model_log, par_b),
  loglik_density = loglik(model_density, par_b),        loglik_log = loglik(model_log, par_b),
  check = FALSE
)
med <- setNames(1e6 * as.numeric(blocks$median), as.character(blocks$expression))

ms_density <- c(1000 * as.numeric(bm_classes$median[1]), med[c("Q_density", "e_step_density", "loglik_density")] / 1000)
ms_log     <- c(1000 * as.numeric(bm_classes$median[2]), med[c("Q_log", "e_step_log", "loglik_log")] / 1000)
fmt_ms <- function(v) ifelse(v >= 10, formatC(v, format = "f", digits = 0), formatC(v, format = "f", digits = 2))
faster_tab <- data.frame(
  what    = c("Whole fit", "One call to `Q()`", "One call to `e_step()`", "One call to `loglik()`"),
  density = fmt_ms(ms_density),
  log     = fmt_ms(ms_log),
  speed   = paste0(formatC(ms_density / ms_log, format = "f", digits = 1), " ×")
)
faster_tab[1, ] <- paste0("**", faster_tab[1, ], "**")
knitr::kable(
  faster_tab, row.names = FALSE, align = c("l", "r", "r", "r"),
  col.names = c("", "`tmix_density` (ms)", "`tmix_log` (ms)", "Speed-up")
)
#
#
#
#
#
#
#
#
#
#| label: s3-profile
#| echo: false
prof_file_log <- profile_expr(for (i in 1:15) fit(model_log))
prof_log <- summaryRprof(prof_file_log)
calls_log <- count_calls(fit(model_log), c("e_step", "Q", "grad_Q"))

share <- function(p, f) { v <- p$by.total[dQuote(f, FALSE), "total.pct"]; if (is.na(v)) 0 else v }
now_tab <- data.frame(
  fun = c("`Q()`", "`e_step()`", "`grad_Q()`"),
  calls_before = formatC(calls[c("Q", "e_step", "grad_Q")], format = "d"),
  calls_now    = formatC(calls_log[c("Q", "e_step", "grad_Q")], format = "d"),
  share_before = paste0(round(c(share(prof, "Q"), share(prof, "e_step"), share(prof, "grad_Q"))), "%"),
  share_now    = paste0(round(c(share(prof_log, "Q.tmix_log"), share(prof_log, "e_step.tmix_log"),
                                share(prof_log, "grad_Q.tmix"))), "%")
)
now_tab[1, ] <- paste0("**", now_tab[1, ], "**")
knitr::kable(
  now_tab, row.names = FALSE, align = c("l", "r", "r", "r", "r"),
  col.names = c("Function", "Calls per fit, version 1", "Calls per fit, version 2",
                "Share of time, version 1", "Share of time, version 2")
)
#
#
#
#| label: plot-prof-compare
#| echo: false
#| fig-height: 3.2
prof_bars <- function(p, funs, labels, version) {
  tab <- p$by.total[dQuote(funs, FALSE), c("total.pct", "self.pct")]
  tab[is.na(tab)] <- 0
  data.frame(
    version = version,
    fun  = rep(labels, 2),
    type = rep(c("total", "self"), each = length(funs)),
    pct  = c(tab$total.pct, tab$self.pct)
  )
}
prof_cmp <- rbind(
  prof_bars(prof, c("Q", "dt_ls", "dt", "e_step", "grad_Q"),
            c("Q()", "dt_ls()", "dt()", "e_step()", "grad_Q()"), "Version 1: density scale"),
  prof_bars(prof_log, c("Q.tmix_log", "ldt_ls", "e_step.tmix_log", "grad_Q.tmix"),
            c("Q()", "ldt_ls()", "e_step()", "grad_Q()"), "Version 2: log scale")
)
prof_cmp$fun <- factor(prof_cmp$fun, levels = rev(c("Q()", "dt_ls()", "dt()", "ldt_ls()", "e_step()", "grad_Q()")))
prof_cmp$type <- factor(prof_cmp$type, levels = c("self", "total"))
ggplot(prof_cmp, aes(x = pct, y = fun, fill = type)) +
  geom_col(position = "dodge") +
  facet_wrap(~version, scales = "free_y") +
  scale_fill_manual(values = c(total = "grey70", self = "steelblue"), breaks = c("total", "self")) +
  labs(x = "share of the run time (%)", y = NULL, fill = NULL)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: methods-accuracy
#| echo: false
methods <- list(
  GEM           = function(m) fit(m),
  BFGS          = function(m) fit_optim(m, "BFGS"),
  `L-BFGS-B`    = function(m) fit_optim(m, "L-BFGS-B"),
  CG            = function(m) fit_optim(m, "CG"),
  `Nelder-Mead` = function(m) fit_optim(m, "Nelder-Mead"),
  nlminb        = function(m) fit_nlminb(m)
)
fits <- lapply(methods, function(f) f(model_log))
ll <- sapply(fits, function(p) loglik(model_log, p))
max_tab <- data.frame(
  method = names(fits),
  formatC(t(sapply(fits, as.vector)), format = "f", digits = 4),
  evals  = formatC(sapply(fits, attr, "iter"), format = "d"),
  gap    = ifelse(max(ll) - ll == 0, "0", formatC(max(ll) - ll, format = "e", digits = 1))
)
max_tab$method[1] <- "**GEM (ours)**"
knitr::kable(
  max_tab, row.names = FALSE, align = c("l", rep("r", 7)),
  col.names = c("Method", "$p$", "$\\mu_1$", "$\\mu_2$", "$\\sigma_1$", "$\\sigma_2$",
                "Evaluations", "Distance to best $\\ell$")
)
#
#
#
#
#
#
#
#| label: press-methods
#| echo: false
press <- bench::press(N = c(500, 2000, 8000, 32000), {
  set.seed(1)
  model <- tmix_log(rtmix(N, par_true, nu), nu)
  bench::mark(
    GEM           = fit(model),
    BFGS          = fit_optim(model, "BFGS"),
    `L-BFGS-B`    = fit_optim(model, "L-BFGS-B"),
    CG            = fit_optim(model, "CG"),
    `Nelder-Mead` = fit_optim(model, "Nelder-Mead"),
    nlminb        = fit_nlminb(model),
    check = FALSE, min_iterations = 3
  )
})
#
#
#
#| label: plot-press-methods
#| echo: false
#| fig-height: 3.6
press_df <- data.frame(method = as.character(press$expression), N = press$N,
                       ms = 1000 * as.numeric(press$median))
ggplot(press_df, aes(N, ms, colour = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_x_log10(breaks = unique(press_df$N)) + scale_y_log10() +
  labs(y = "median time (ms)", colour = NULL)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: gauss-one
#| echo: false
set.seed(1)
x_cont <- c(x, runif(10, 8, 20) * sample(c(-1, 1), 10, replace = TRUE))
one_tab <- rbind(
  `True parameters`              = par_true,
  `$t$ mixture, clean`           = fit(tmix_log(x, nu), par0 = par_true),
  `Gaussian mixture, clean`      = fit_gauss(x, par0 = par_true),
  `$t$ mixture, 10 outliers`     = fit(tmix_log(x_cont, nu), par0 = par_true),
  `Gaussian mixture, 10 outliers` = fit_gauss(x_cont, par0 = par_true)
)
knitr::kable(
  one_tab, digits = 2, align = "r",
  col.names = c("$p$", "$\\mu_1$", "$\\mu_2$", "$\\sigma_1$", "$\\sigma_2$")
)
#
#
#
#
#
#| label: plot-gauss-fit
#| echo: false
#| fig-height: 4.4
dgmix <- function(x, par) par[1] * dnorm(x, par[2], par[4]) + (1 - par[1]) * dnorm(x, par[3], par[5])
fit_t_cont <- one_tab[4, ]
fit_g_cont <- one_tab[5, ]
grid_cont <- seq(min(x_cont) - 1, max(x_cont) + 1, length.out = 1000)
df_cont <- data.frame(
  x = rep(grid_cont, 5),
  density = c(
    dgmix(grid_cont, fit_g_cont),
    fit_g_cont[1] * dnorm(grid_cont, fit_g_cont[2], fit_g_cont[4]),
    (1 - fit_g_cont[1]) * dnorm(grid_cont, fit_g_cont[3], fit_g_cont[5]),
    dtmix(grid_cont, fit_t_cont, nu),
    dtmix(grid_cont, par_true, nu)
  ),
  curve = factor(rep(c("Gaussian mixture", "Gaussian component 1", "Gaussian component 2",
                       "t mixture", "true density"), each = length(grid_cont)),
                 levels = c("Gaussian mixture", "Gaussian component 1", "Gaussian component 2",
                            "t mixture", "true density"))
)
ggplot(data.frame(x = x_cont), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "grey85", colour = "white") +
  geom_rug(data = data.frame(x = tail(x_cont, 10)), colour = "#d95f02", linewidth = 1,
           length = unit(0.06, "npc")) +
  geom_line(data = df_cont, aes(x, density, colour = curve, linetype = curve), linewidth = 0.9) +
  scale_colour_manual(values = c("#d95f02", "#d95f02", "#d95f02", "#2c7fb8", "black")) +
  scale_linetype_manual(values = c(1, 2, 3, 1, 2)) +
  labs(colour = NULL, linetype = NULL, y = "density",
       title = "N = 400 and 10 outliers (orange marks)")
#
#
#
#
#
#
#
#
#
#
#
#| label: gauss-outliers
#| echo: false
set.seed(3)
k_out <- c(0, 2, 5, 10, 20, 40)
out_res <- do.call(rbind, lapply(1:20, function(rep) {
  x0 <- rtmix(400, par_true, nu)
  out <- runif(max(k_out), 8, 20) * sample(c(-1, 1), max(k_out), replace = TRUE)
  do.call(rbind, lapply(k_out, function(k) {
    xc <- c(x0, out[seq_len(k)])
    rbind(
      data.frame(k, model = "t mixture", t(fit(tmix_log(xc, nu), par0 = par_true))),
      data.frame(k, model = "Gaussian mixture", t(fit_gauss(xc, par0 = par_true)))
    )
  }))
}))
out_med <- aggregate(cbind(mu1, mu2, sigma1, sigma2) ~ k + model, out_res, median)
out_long <- reshape(out_med, direction = "long", varying = par_names[-1], v.names = "estimate",
                    timevar = "parameter", times = par_names[-1])
med_at <- function(model, k, par) out_med[out_med$model == model & out_med$k == k, par]
#
#
#
#| label: plot-gauss-outliers
#| echo: false
#| fig-height: 3.4
ggplot(out_long, aes(k, estimate, colour = model)) +
  geom_hline(data = data.frame(parameter = par_names[-1], value = par_true[-1]),
             aes(yintercept = value), linetype = 2) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~parameter, nrow = 1, scales = "free_y") +
  labs(x = "number of outliers", y = "median estimate", colour = NULL) +
  theme(legend.position = "bottom")
#
#
#
#
#
#
#
#
#
#
#
#
