# Source all code (run from the project folder)
for (f in c("tmix-simple.R", "gem.R", "profiling-tools.R")) {
  source(file.path("R", f), keep.source = TRUE)
}
