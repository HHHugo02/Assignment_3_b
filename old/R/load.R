# Source all code for the assignment (run from the project folder)
for (f in c("tmix-model.R", "em-algorithms.R", "direct-methods.R", "tracer.R",
            "fisher.R", "tmix-class.R", "v1-naive.R")) {
  source(file.path("R", f), keep.source = TRUE)
}
Rcpp::sourceCpp(file.path("src", "em_aug_step.cpp"))
