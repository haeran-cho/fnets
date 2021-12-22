require(fnets)
set.seed(222)
n <- 200
p<- 100
chi <- sim.factor.M2(n,p)
xi <- sim.idio(n,p)
sample.data <- chi$data + xi$data
model <- fnets(sample.data, q=2, idio.method = "lasso")
pr <- predict(model,sample.data, common.method = "static")
cpre <- common.predict(model,sample.data, common.method = "static")
ip <- idio.predict(model,sample.data, cpre)
