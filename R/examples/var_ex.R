library(fnets)

set.seed(123)
n <- 500
p <- 50
idio <- sim.var(n, p)
x <- idio$data

fv <- fit.var(x, center = TRUE, method = 'lasso', var.order = 1)
norm(fv$beta - t(idio$A), 'F')/norm(t(idio$A), 'F')
