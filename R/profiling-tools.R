# Helpers for the profiling analysis --------------------------------------------

# Count how many times the functions `fun_names` are called while `expr` runs.
# Each function is temporarily replaced in the global environment by a wrapper
# that increments a counter and then calls the original function.
count_calls <- function(expr, fun_names) {
  counts <- setNames(integer(length(fun_names)), fun_names)
  originals <- mget(fun_names, envir = globalenv())
  on.exit(
    for (f in fun_names) assign(f, originals[[f]], envir = globalenv())
  )
  for (f in fun_names) {
    local({
      name <- f
      original <- originals[[f]]
      assign(name, function(...) {
        counts[name] <<- counts[name] + 1L
        original(...)
      }, envir = globalenv())
    })
  }
  force(expr)
  counts
}

# Run Rprof (with line profiling) on `expr` and return the file with the samples
profile_expr <- function(expr, interval = 0.005) {
  prof_file <- tempfile()
  Rprof(prof_file, interval = interval, line.profiling = TRUE)
  force(expr)
  Rprof(NULL)
  prof_file
}

# Line-level profile: one row per source line, with the function the line
# belongs to and the code on that line. Rprof reports lines as "file.R#n";
# the code is read from the R/ folder.
profile_lines <- function(prof_file, n = 8) {
  tab <- summaryRprof(prof_file, lines = "show")$by.total
  loc <- rownames(tab)
  info <- t(vapply(loc, function(l) {
    parts <- strsplit(l, "#", fixed = TRUE)[[1]]
    path <- file.path("R", parts[1])
    if (length(parts) != 2 || !file.exists(path)) {
      return(c("", ""))
    }
    src <- readLines(path)
    i <- as.integer(parts[2])
    starts <- grep("^[A-Za-z_.][A-Za-z0-9_.]* <- function", src[seq_len(i)])
    fun <- if (length(starts) > 0) sub(" <- function.*", "", src[max(starts)]) else ""
    c(fun, trimws(src[i]))
  }, character(2)))
  out <- data.frame(
    line = loc, in_function = info[, 1], code = info[, 2],
    total.pct = tab$total.pct, self.pct = tab$self.pct
  )
  out <- out[out$code != "" & !startsWith(out$line, "profiling-tools.R"), ]
  rownames(out) <- NULL
  head(out, n)
}
