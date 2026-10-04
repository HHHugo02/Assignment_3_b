# Source version 1 of the code: plain functions (run from the project folder).
# The S3 version is loaded with R/load.R instead.
for (f in c("tmix-common.R", "tmix-simple.R", "em.R", "gem.R", "profiling-tools.R")) {
  source(file.path("R", f), keep.source = TRUE)
}
