# Source the S3 version of the code (run from the project folder).
# Version 1 (plain functions) is loaded with R/load-v1.R instead.
for (f in c("tmix-common.R", "tmix-classes.R", "em.R", "gem-classes.R", "optim-methods.R",
            "profiling-tools.R")) {
  source(file.path("R", f), keep.source = TRUE)
}
