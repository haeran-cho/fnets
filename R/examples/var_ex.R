library(fnets)

set.seed(123)
n <- 500
p <- 50
idio <- sim.var(n, p)
x <- idio$data

fv <- fnets.var(x,
  center = TRUE, method = "lasso", var.order = 1,
  tuning.args = list(tuning = "cv", n.folds = 1, path.length = 10, do.plot = TRUE),
  n.cores = 2
)
norm(fv$beta - t(idio$A), "F") / norm(t(idio$A), "F")
