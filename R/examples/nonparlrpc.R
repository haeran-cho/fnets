#nonpar.lrpc
require(doParallel)
require(lpSolve)
set.seed(222)
n <- 200
p<- 100
chi <- sim.factor.M2(n,p)
xi <- sim.idio(n,p)
sample.data <- chi$data + xi$data
model <- fnets(sample.data, q=2, idio.method = "lasso")
nonpar.lrpc(model, sample.data, 1, n.cores = 1)
