# Interactive profiling demo (run in RStudio from the project folder) ----------
source("R/load.R")
library(profvis)

set.seed(7)
nu <- c(5, 10)
x_big <- rtmix(5000, c(0.6, 0, 3, 1, 0.7), nu)
p0_big <- default_start(x_big)

# GEM with cold-started line search: most time is spent in Q() -> ldt_ls()
profvis(for (i in 1:10) em(p0_big, create_gem_step(x_big, nu, warm = FALSE)))

# GEM with warm-started line search
profvis(for (i in 1:10) em(p0_big, create_gem_step(x_big, nu, warm = TRUE)))

# Augmented EM: plogis(), ldt_ls() and cbind() are the hot spots
profvis(for (i in 1:10) em(p0_big, create_em_aug_step(x_big, nu), maxit = 100))
