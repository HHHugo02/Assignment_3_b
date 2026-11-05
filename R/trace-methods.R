# Convergence traces, so all methods are compared on the same criterion ---------
#
# The methods stop on different criteria: the GEM on the change in par, optim()
# on the relative change in the function value. Their times to
# convergence are therefore not comparable. Instead every method records each
# point it visits together with the elapsed time. The log-likelihood of the
# points is computed afterwards, outside the timing, and a method is judged by
# when its best log-likelihood so far first comes within a given gap of the
# maximum. The methods are run with very tight tolerances, so they run past
# every gap we look at.

now <- function() as.numeric(Sys.time())

# GEM: one point per iteration
trace_gem <- function(model, par0 = default_start(model$x), eps = 1e-20, maxit = 20000) {
  gem_step <- create_gem_step(model)
  pars <- vector("list", maxit + 1)
  times <- numeric(maxit + 1)
  pars[[1]] <- par0
  par <- par0
  t0 <- now()
  for (i in seq_len(maxit)) {
    par0 <- par
    par <- gem_step(par)
    times[i + 1] <- now() - t0
    pars[[i + 1]] <- par
    if (sum((par - par0)^2) <= eps * (sum(par^2) + eps)) {
      break
    }
  }
  list(time = times[seq_len(i + 1)], par = do.call(rbind, pars[seq_len(i + 1)]))
}

# optim(): one point per evaluation of the objective. The
# objective is wrapped so it records eta and the time after each evaluation
trace_optimizer <- function(model, method, par0 = default_start(model$x), tol = 1e-16) {
  f <- neg_loglik_eta(model)
  gr <- neg_score_eta(model)
  rec_eta <- list()
  rec_time <- numeric()
  t0 <- now()
  fn <- function(eta) {
    value <- f(eta)
    rec_eta[[length(rec_eta) + 1]] <<- eta
    rec_time[length(rec_time) + 1] <<- now() - t0
    value
  }
  eta0 <- par_to_eta(par0)
  control <- if (method == "L-BFGS-B") list(factr = 1, maxit = 20000) else list(reltol = tol, maxit = 20000)
  optim(eta0, fn, if (method == "Nelder-Mead") NULL else gr, method = method, control = control)
  list(time = c(0, rec_time), par = rbind(par0, t(vapply(rec_eta, eta_to_par, numeric(5)))))
}

trace_method <- function(model, method, par0 = default_start(model$x)) {
  if (method == "GEM") trace_gem(model, par0) else trace_optimizer(model, method, par0)
}

# Best log-likelihood so far at each recorded time. The line searches of the
# quasi-Newton methods also evaluate trial points that are worse, hence cummax
best_loglik <- function(model, tr) {
  ll <- apply(tr$par, 1, function(p) loglik(model, p))
  ll[!is.finite(ll)] <- -Inf
  data.frame(time = tr$time, loglik = cummax(ll))
}

# First time the best log-likelihood is within `gap` of `ll_max` (NA if never)
time_to_gap <- function(best, ll_max, gap) {
  k <- which(ll_max - best$loglik <= gap)
  if (length(k) == 0) NA else best$time[k[1]]
}
