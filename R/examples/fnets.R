set.seed(222)
n <- 200
p<- 100
chi <- sim.factor.M1(n,p)
xi <- sim.idio(n,p)
sample.data <- chi$data + xi$data
fnets(sample.data, q=2, idio.method = "lasso")
