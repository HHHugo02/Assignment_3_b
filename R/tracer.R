# A minimal tracer in the spirit of CSwR::tracer() ----------------------------
#
# tr <- make_tracer("par"); em(par0, step, cb = tr$tracer); summary(tr)
#
# Each call stores the variables `fields` from the calling frame and the time
# spent *outside* the tracer since the previous call, so the tracing overhead
# is not counted in the run time.

make_tracer <- function(fields = "par") {
  store <- list()
  elapsed <- numeric()
  n <- 0
  t_last <- as.numeric(bench::hires_time())

  tracer <- function() {
    t_now <- as.numeric(bench::hires_time())
    n <<- n + 1
    elapsed[n] <<- t_now - t_last
    store[[n]] <<- unlist(mget(fields, envir = parent.frame()))
    t_last <<- as.numeric(bench::hires_time())
    invisible(NULL)
  }

  # restart the clock (call right before running the algorithm)
  start <- function() {
    store <<- list()
    elapsed <<- numeric()
    n <<- 0
    t_last <<- as.numeric(bench::hires_time())
    invisible(NULL)
  }

  structure(
    list(tracer = tracer, start = start, env = environment()),
    class = "tmix_tracer"
  )
}

summary.tmix_tracer <- function(object, ...) {
  env <- object$env
  vals <- do.call(rbind, env$store)
  data.frame(n = seq_len(env$n), time = cumsum(env$elapsed), vals)
}

# Trace an algorithm and add the log-likelihood (computed *after* the run)
trace_run <- function(run, x, nu, label) {
  tr <- make_tracer("par")
  tr$start()
  run(tr$tracer)
  tr_df <- summary(tr)
  pars <- as.matrix(tr_df[, -(1:2)])
  tr_df$loglik <- apply(pars, 1, loglik, x = x, nu = nu)
  tr_df$method <- label
  tr_df
}
