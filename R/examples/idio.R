require(fnets)
require(doParallel)
require(lpSolve)
set.seed(222)
n <- 200
p<- 100
chi <- sim.factor.M2(n,p)
xi <- sim.idio(n,p)
sample.data <- chi$data + xi$data
fit.var(sample.data, .1, idio.method = "lasso")
