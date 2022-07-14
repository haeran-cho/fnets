library(fnets)

set.seed(123)
n <- 500
p <- 50
common <- sim.common2(n, p)
idio <- sim.var(n, p)
x <- common$data * apply(idio$data, 1, sd) / apply(common$data, 1, sd) + idio$data

bn <- bn.factor.number(x, q.max = NULL, center = FALSE, do.plot = TRUE)
bn$q.hat
