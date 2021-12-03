#nonpar.pcn

require(doParallel)
require(lpSolve)
model <- fnets(sample.data, q=2, idio.method = "lasso")
nonpar.pcn(model, sample.data, 1, n.cores = 1)
