# Interactive profiling demo (run in RStudio from the project folder) ----------
source("R/load-v1.R")
library(profvis)

set.seed(7)
nu <- c(5, 10)
x_big <- rtmix(5000, c(0.6, 0, 3, 1, 0.7), nu)

# The flame graph shows Q() -> dt_ls() -> dt() inside the line search of
# create_gem_step() as the dominating part of the run time.
profvis(for (i in 1:5) gem(x_big, nu))
